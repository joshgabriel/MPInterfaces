# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import


"""
Utility functions
"""

from six.moves import range, zip

import itertools as it
from functools import reduce
import linecache
import sys
import os
import math
import socket
import time
import subprocess as sp
import logging
from collections import OrderedDict, Counter
import yaml
from argparse import ArgumentParser
from glob import glob
import shutil as shu
import numpy as np
import pandas as pd

from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn

from pymatgen.core.sites import PeriodicSite
from pymatgen import Structure, Lattice, Element
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Incar,Poscar
from pymatgen.core.composition import Composition
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import _pt_data
from pymatgen.io.vasp.outputs import Vasprun, Oszicar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.elasticity.strain import Strain,Deformation

from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler

from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from ase.build import surface

from mpinterfaces.default_logger import get_default_logger
from mpinterfaces import *

__author__ = "Kiran Mathew, Joshua J. Gabriel, Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

logger = get_default_logger(__name__)

ELEMENT_RADII = {i: Element(i).atomic_radius for i in _pt_data}


def get_ase_slab(pmg_struct, hkl=(1, 1, 1), min_thick=10, min_vac=10):
    """
    takes in the intial structure as pymatgen Structure object
    uses ase to generate the slab
    returns pymatgen Slab object

    Args:
        pmg_struct: pymatgen structure object
        hkl: hkl index of surface of slab to be created
        min_thick: minimum thickness of slab in Angstroms
        min_vac: minimum vacuum spacing
    """
    ase_atoms = AseAtomsAdaptor().get_atoms(pmg_struct)
    pmg_slab_gen = SlabGenerator(pmg_struct, hkl, min_thick, min_vac)
    h = pmg_slab_gen._proj_height
    nlayers = int(math.ceil(pmg_slab_gen.min_slab_size / h))
    ase_slab = surface(ase_atoms, hkl, nlayers)
    ase_slab.center(vacuum=min_vac / 2, axis=2)
    pmg_slab_structure = AseAtomsAdaptor().get_structure(ase_slab)
    return Slab(lattice=pmg_slab_structure.lattice,
                species=pmg_slab_structure.species_and_occu,
                coords=pmg_slab_structure.frac_coords,
                site_properties=pmg_slab_structure.site_properties,
                miller_index=hkl, oriented_unit_cell=pmg_slab_structure,
                shift=0., scale_factor=None, energy=None)


def slab_from_file(hkl, filename):
    """
    reads in structure from the file and returns slab object.
    useful for reading in 2d/substrate structures from file.
    Args:
         hkl: miller index of the slab in the input file.
         filename: structure file in any format
                   supported by pymatgen
    Returns:
         Slab object
    """
    slab_input = Structure.from_file(filename)
    return Slab(slab_input.lattice,
                slab_input.species_and_occu,
                slab_input.frac_coords,
                hkl,
                Structure.from_sites(slab_input, to_unit_cell=True),
                shift=0,
                scale_factor=np.eye(3, dtype=np.int),
                site_properties=slab_input.site_properties)


def get_magmom_init(poscar,mag_init_file='mag_inits.yaml',is_spin_orbit=False):
    """
    Sets the mag init to ions found in the given poscar
    according to the values given in the yaml file 
    """
    magmom = []
    mag_init =0
    mag_inits=[]
    try:
       if type(mag_init_file)==str:
           print ('opening mag_init.yaml')
           mag_inits = yaml.load(open(mag_init_file))
       elif type(mag_init_file)==float:
           print ('mag init float')
           mag_init = mag_init_file
       elif type(mag_init_file)==list and not is_spin_orbit:
           magmom = mag_init_file
    except:
       logger.warn('mag_inits.yaml not found and manual initialization not specified,\
                   so 0.5 is being used as initialization \
                   for each element in your poscar file')
       mag_inits = []
       
    sites_dict = poscar.as_dict()['structure']['sites']

    print (poscar.comment,magmom)

    if len(magmom) == len(sites_dict):
       return magmom

    for n, s in enumerate(sites_dict):

        if s['label'] in mag_inits and not is_spin_orbit:
            magmom.append(mag_inits[s['label']])
        elif mag_init and not is_spin_orbit:
            magmom.append(mag_init)
        elif not is_spin_orbit:
            magmom.append(0.5)
        elif s['label'] in mag_inits and is_spin_orbit:
            magmom += [mag_inits[s['label']] for i in [0,1,2]]
        elif mag_init and is_spin_orbit:
            magmom += [mag_init for i in [0,1,2]]
        elif type(mag_init_file)==list and is_spin_orbit:
            magmom += [mag_init_file[n] for i in [0,1,2]]
        else:
            magmom += [0.5 for i in [0,1,2]]
    print (magmom)
    return magmom
    

def get_magmom_string(structure):
    """
    Based on a POSCAR, returns the string required for the MAGMOM
    setting in the INCAR. Initializes transition metals with 6.0
    bohr magneton and all others with 0.5.

    Args:
        structure (Structure): Pymatgen Structure object

    Returns:
        string with INCAR setting for MAGMOM according to mat2d
        database calculations
    """
    magmoms, considered = [], []
    for s in structure.sites:
        if s.specie not in considered:
            amount = int(structure.composition[s.specie])
            if s.specie.is_transition_metal:
                magmoms.append('{}*6.0'.format(amount))
            else:
                magmoms.append('{}*0.5'.format(amount))
            considered.append(s.specie)
    return ' '.join(magmoms)


def get_magmom_mae(poscar, mag_init):
    """
    mae
    """
    mae_magmom = []

    sites_dict = poscar.as_dict()['structure']['sites']

    # initialize a magnetic moment on the transition metal
    # in vector form on the x-direction
    for n, s in enumerate(sites_dict):

        if Element(s['label']).is_transition_metal:
            mae_magmom.append([0.0, 0.0, mag_init])
        else:
            mae_magmom.append([0.0, 0.0, 0.0])

    return sum(mae_magmom, [])


def get_magmom_afm(poscar, database=None):
    """
    returns the magmom string which is an N length list
    """

    afm_magmom = []
    orig_structure_name = poscar.comment

    if len(poscar.structure) % 2 != 0:

        if database == 'twod':
            # no need for more vacuum spacing
            poscar.structure.make_supercell([2, 2, 1])
        else:
            # for bulk structure
            poscar.structure.make_supercell([2, 2, 2])

    sites_dict = poscar.as_dict()['structure']['sites']

    for n, s in enumerate(sites_dict):

        if Element(s['label']).is_transition_metal:
            if n % 2 == 0:
                afm_magmom.append(6.0)
            else:
                afm_magmom.append(-6.0)

        else:
            if n % 2 == 0:
                afm_magmom.append(0.5)
            else:
                afm_magmom.append(-0.5)

    return afm_magmom, Poscar(structure=poscar.structure,
                              comment=orig_structure_name)

def get_magmom_hunds(poscar):
    """
    returns an Incar magmom string according to the 
    Hunds rule configuration of unpaired electrons
    """
    pass

def get_run_cmmnd(nnodes=1, ntasks=16, walltime='10:00:00', job_bin=None,
                  job_name=None, mem=None, partition='hpg1-compute'):
    """
    returns the fireworks CommonAdapter based on the queue
    system specified by mpint_config.yaml and the submit
    file template also specified in mpint_config.yaml
    NOTE: for the job_bin, please specify the mpi command as well:
          Eg: mpiexec /path/to/binary
    """
    d = {}
    job_cmd = None
    qtemp_file = open(QUEUE_TEMPLATE+os.sep+'qtemplate.yaml')
    qtemp = yaml.load(qtemp_file)
    qtemp_file.close()
    qtemp.update({'nnodes': nnodes, 'ntasks':ntasks, 'walltime': walltime, \
                  'rocket_launch': job_bin, 'job_name':job_name,'mem':mem})
    # SLURM queue
    if QUEUE_SYSTEM == 'slurm':
        qtemp['ntasks_per_node'] = int(qtemp['ntasks']/qtemp['nnodes'])
        # say nnodes = 1, gives ntasks_per_node = 32, divide into 16 per socket
        # say nnodes = 2, gives ntasks_per_node = 16, divide into 8 per socket
        # say nnodes = 4, gives ntasks_per_node = 8, divide into 4 per socket

        if qtemp['partition'] == 'hpg1-compute':
            qtemp['ntasks_per_socket'] = int(qtemp['ntasks_per_node']/4)

        elif qtemp['partition'] == 'hpg2-compute':
            qtemp['ntasks_per_socket'] = int(qtemp['ntasks_per_node']/2)

        else:
            del qtemp['partition']

        print ('Q DONE')

        if job_bin is None:
            job_bin = VASP_STD_BIN
        else:
            job_bin = job_bin
        d = {'type': 'SLURM',
             'params': qtemp}
    # PBS queue
    elif QUEUE_SYSTEM == 'pbs':
        if job_bin is None:
            job_bin = VASP_STD_BIN
        else:
            job_bin = job_bin
        d = {'type': 'PBS',
             'params': qtemp}
    else:
        job_cmd = ['ls', '-lt']
    if d:
        return (CommonAdapter(d['type'], **d['params']), job_cmd)
    else:
        return (None, job_cmd)


def get_job_state(job):
    """
    Args:
        job: job

    Returns:
           the job state and the job output file name
    """
    ofname = None

    # pbs
    if QUEUE_SYSTEM == 'pbs':# in hostname:
        try:
            output = sp.check_output(['qstat', '-i', job.job_id])
            state = output.rstrip('\n').split('\n')[-1].split()[-2]
        except:
            logger.info('Job {} not in the que'.format(job.job_id))
            state = "00"
        err_file = glob(job.job_dir+os.sep+'*.error')[-1]
        out_file = err_file.split('.')[0] + '.out'

    # slurm
    elif QUEUE_SYSTEM == 'slurm':
        try:
            output = str(sp.check_output(['squeue', '--job', job.job_id]))
            state = output.rstrip('\n').split('\n')[-1].split()[-4]
        except:
            logger.info('Job {} not in the que.'.format(job.job_id))
            logger.info(
                'This could mean either the batchsystem crashed(highly unlikely) or the job completed a long time ago')
            state = "00"
        err_file = glob(job.job_dir+os.sep+'*.error')[-1]
        out_file = err_file.split('.')[0] + '.out'


    # no batch system
    else:
        state = 'XX'
        ofname = 'MPInt_Job.out'
    return state, ofname


def update_checkpoint(job_ids=None, jfile=None, **kwargs):
    """
    rerun the jobs with job ids in the job_ids list. The jobs are
    read from the json checkpoint file, jfile.
    If no job_ids are given then the checkpoint file will
    be updated with corresponding final energy

    Args:
        job_ids: list of job ids to update or q resolve
        jfile: check point file
    """
    cal_log = loadfn(jfile, cls=MontyDecoder)
    cal_log_new = []
    all_jobs = []
    run_jobs = []
    handlers = []
    final_energy = None
    incar = None
    kpoints = None
    qadapter = None
    # if updating the specs of the job
    for k, v in kwargs.items():
        if k == 'incar':
            incar = v
        if k == 'kpoints':
            kpoints = v
        if k == 'que':
            qadapter = v
    for j in cal_log:
        job = j["job"]
        job.job_id = j['job_id']
        all_jobs.append(job)
        if job_ids and (j['job_id'] in job_ids or job.job_dir in job_ids):
            logger.info('setting job {0} in {1} to rerun'.format(j['job_id'],
                                                                 job.job_dir))
            contcar_file = job.job_dir + os.sep + 'CONTCAR'
            poscar_file = job.job_dir + os.sep + 'POSCAR'
            if os.path.isfile(contcar_file) and len(
                    open(contcar_file).readlines()) != 0:
                logger.info('setting poscar file from {}'
                            .format(contcar_file))
                job.vis.poscar = Poscar.from_file(contcar_file)
            else:
                logger.info('setting poscar file from {}'
                            .format(poscar_file))
                job.vis.poscar = Poscar.from_file(poscar_file)
            if incar:
                logger.info('incar overridden')
                job.vis.incar = incar
            if kpoints:
                logger.info('kpoints overridden')
                job.vis.kpoints = kpoints
            if qadapter:
                logger.info('qadapter overridden')
                job.vis.qadapter = qadapter
            run_jobs.append(job)
    if run_jobs:
        c = Custodian(handlers, run_jobs, max_errors=5)
        c.run()
    for j in all_jobs:
        final_energy = j.get_final_energy()
        cal_log_new.append({"job": j.as_dict(),
                            'job_id': j.job_id,
                            "corrections": [],
                            'final_energy': final_energy})
    dumpfn(cal_log_new, jfile, cls=MontyEncoder, indent=4)


def jobs_from_file(filename='calibrate.json'):
    """
    read in json file of format caibrate.json(the default logfile
    created when jobs are run through calibrate) and return the
    list of job objects.

    Args:
        filename: checkpoint file name

    Returns:
           list of all jobs
    """
    caljobs = loadfn(filename, cls=MontyDecoder)
    all_jobs = []
    for j in caljobs:
        job = j["job"]
        job.job_id = j['job_id']
        job.final_energy = j['final_energy']
        all_jobs.append(job)
    return all_jobs

def continue_job_inputs(chkpt_files=None, old_job_dir=None, user_filters=None,\
                     rerun_nonconverged_also=False):
   """
   utility function for any workflow that can be used to decide which
   materials pass to the next stage of a workflow
   """
   if chkpt_files:
      # Fuel source from checkpoint files
      rerun_paths = []
      chkpts = sum( [ jobs_from_file(j) for j in \
                    [chks for chks in chkpt_files]], []  )

      all_converged = [ (j.parent_job_dir+os.sep+j.job_dir,\
        j.job_dir.split('/')[-1].split('_')[0],\
        j.job_dir.split('/')[-1].split('_')[-1])\
        for j in chkpts if j.final_energy ]

      restarts = [ (j.parent_job_dir+os.sep+j.job_dir,\
        j.job_dir.split('/')[-1].split('_')[0],\
        j.job_dir.split('/')[-1].split('_')[-1])\
        for j in chkpts if not j.final_energy ]

      if user_filters:
           # filter based on compound and structure name tags
           # can be on more involved factors like % of jobs done in the step
           filtered_run = user_filters#yaml.load(open(user_filters))
           all_jobs = chkpts
           #print (list(filtered_run.keys()))
           converged = [p[0] for p in all_converged]
           if 'Completion' in list(filtered_run.keys()):
               if len(converged)/len(all_jobs) > filtered_run['Completion']:
                   rerun_paths = converged
               else:
                   print ('{0} of jobs in step {1} not complete yet'\
                           .format(filtered_run['Completion'], c))
           elif 'Condition' in list(filtered_run.keys()):
               script = filtered_run['Condition']
               p = subprocess.Popen(['python',script, '-i', converged], \
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
               stdout, stderr = p.communicate()
               if stdout: 
                  rerun_paths = stdout
               else:
                  print ('Condition in {0} not met by any job in {1}'.format(script, c))    
           elif 'Dir' in list(filtered_run.keys()):
               rerun_paths = filtered_run['Dir']
               print (rerun_paths)

                   
#           for (job_path, compound, structure) in all_converged:
#              for f in filtered_run.keys():
#                 for s in filtered_run[f]:
#                    if f == compound and s == structure:
#                       rerun_paths.append(job_path)
      else:
          rerun_paths = [p[0] for p in all_converged]

   elif old_job_dir:
       # Fuel source from list of directories
       rerun_paths = old_job_dirs

   if rerun_nonconverged_also:
       return rerun_paths, restarts
   else:
       return rerun_paths

def check_job_with_custodian(handlers,logfile,job=None,jdir=None):
    """
    function which checks a given checkpoint job
    or a job directory for errors and performs
    either just in time correction (writes STOPCAR to a job
    that is running with errors and re-submits to the queue
    with a corrected input set.
    or correction
    """
    corrected={'ErrorDir':[],'Error':[],'Correction':[]}
    if job:
      j = job
      with open(j.parent_job_dir+os.sep+j.job_dir+os.sep+j.output_file) as w:
        j_num= w.read()
      w.close()
      ofiles = glob(j.parent_job_dir+os.sep+j.job_dir+os.sep+'*-*.out')
      if len(ofiles)>1:
         j_nums = [int(o.replace('.out','').split('-')[-1]) for o in ofiles]
         j_num = max(j_nums)
         print (j_num)
      out_file = '-'.join([j.vis.qadapter['job_name'],str(j_num)]) + '.out'
      err_file = '-'.join([j.vis.qadapter['job_name'],str(j_num)]) + '.error'
      print (out_file)
      output_path = j.parent_job_dir+os.sep+j.job_dir
      if os.path.exists(output_path+os.sep+out_file):
        for h in handlers:
           h.output_filename = out_file
           os.chdir(j.job_dir)
           try:
               if h.check():
                   logfile.info('Making correction to {0} \n \
                               and detected errors are {1} \n \
                               Custodian correction is {2}'.\
                           format(output_path,h.errors, h.correct()))
                   #logfile.info('Detected errors are: {}'.format(h.errors))
                   #logfile.info('Correction from Custodian is {}'.format(h.correct()))
                   corrected['ErrorDir'].append(output_path)
                   corrected['Error'].append(h.errors)
                   corrected['Correction'].append(h.correct())
               else:
                   if 'Exceeded job memory' in sp.run(['tail', '-n2', err_file], stdout=sp.PIPE).stdout.decode('utf-8'):
                      logfile.info('Memory Error {}'.format(output_path))
                      corrected['ErrorDir'].append(output_path)
                      corrected['Error'].append('Memory Error')
                      corrected['Correction'].append('Change Memory settings for qparams')
                   elif 'Walltime' in sp.run(['tail', '-n2', err_file], stdout=sp.PIPE).stdout.decode('utf-8'):
                      logfile.info('Walltime exceeded {}'.format(output_path))
                      corrected['ErrorDir'].append(output_path)
                      corrected['Error'].append('Walltime Exceeded')
                      corrected['Correction'].append('Copy CONTCAR to POSCAR and Update the Walltime')
                   else:
                      logfile.info('Error with no Custodian Correction found {} printing last 5 lines of err and out file'.format(output_path))
                      print ('Here')
                      errf = sp.run(['tail', '-n5', err_file],stdout=sp.PIPE).stdout.decode('utf-8')
                      outf = 'Job {}\n'.format(str(j_num)) + sp.run(['tail', '-n10', out_file],stdout=sp.PIPE).stdout.decode('utf-8') 
                      corrected['ErrorDir'].append(output_path)
                      corrected['Error'].append({'out':outf,'err':errf})
                      corrected['Correction'].append('Check Directory')
                      print (corrected)
           except:
                 print (print_exception())
                 logfile.info('Custodian encountered an error')
                 corrected['ErrorDir'].append(output_path)
                 corrected['Error'].append(sp.run(['tail', '-n2', err_file], stdout=sp.PIPE).stdout.decode('utf-8'))
                 corrected['Correction'].append('XX')

           os.chdir(j.parent_job_dir)
      else:
        logfile.warn('File does not exist: {}'.format(output_path+os.sep+out_file))

      if corrected:
          return corrected
      else:
          return None

    elif jdir:
      pass

    else:
      logger.warn('No job or job directory specified for custodian to check...')
      return None

def launch_daemon(steps, interval, handlers=None, ld_logger=None):
    """
    run all the 'steps' in daemon mode
    checks job status every 'interval' seconds
    also runs all the error handlers
    """
    if ld_logger:
        global logger
        logger = ld_logger
    chkpt_files_prev = None
    for step in steps:
        chkpt_files = step(checkpoint_files=chkpt_files_prev)
        chkpt_files_prev = chkpt_files
        if not chkpt_files:
            return None
        while True:
            done = []
            reruns = []
            for cf in chkpt_files:
                time.sleep(3)
                update_checkpoint(job_ids=reruns, jfile=cf)
                all_jobs = jobs_from_file(cf)
                for j in all_jobs:
                    state, ofname = get_job_state(j)
                    if j.final_energy:
                        done = done + [True]
                    elif state == 'R':
                        logger.info('job {} running'.format(j.job_id))
                        done = done + [False]
                    elif state in ['C', 'CF', 'F', '00']:
                        logger.error(
                            'Job {0} in {1} cancelled or failed. State = {2}'.
                            format(j.job_id, j.job_dir, state))
                        done = done + [False]
                        if handlers:
                            logger.info('Investigating ... ')
                            os.chdir(j.job_dir)
                            if ofname:
                                if os.path.exists(ofname):
                                    for h in handlers:
                                        h.output_filename = ofname
                                        if h.check():
                                            logger.error(
                                                'Detected vasp errors {}'.format(
                                                    h.errors))
                                            # TODO: correct the error and mark the job for rerun
                                            # all error handling must done using proper errorhandlers
                                            # h.correct()
                                            # reruns.append(j.job_id)
                                else:
                                    logger.error(
                                        'stdout redirect file not generated, job {} will be rerun'.format(
                                            j.job_id))
                                    reruns.append(j.job_id)
                            os.chdir(j.parent_job_dir)
                    else:
                        logger.info(
                            'Job {0} pending. State = {1}'.format(j.job_id,
                                                                  state))
                        done = done + [False]
            if all(done):
                logger.info(
                    'all jobs in {} done. Proceeding to the next one'.format(
                        step.__name__))
                time.sleep(5)
                break
            logger.info(
                'all jobs in {0} NOT done. Next update in {1} seconds'.format(
                    step.__name__, interval))
            time.sleep(interval)


def get_convergence_data(jfile, params=('ENCUT', 'KPOINTS')):
    """
    returns data dict in the following format
    {'Al':
          {'ENCUT': [ [500,1.232], [600,0.8798] ],
            'KPOINTS':[ [], [] ]
          },
     'W': ...
    }

    Note: processes only INCAR parmaters and KPOINTS
    """
    cutoff_jobs = jobs_from_file(jfile)
    data = {}
    for j in cutoff_jobs:
        jdir = os.path.join(j.parent_job_dir, j.job_dir)
        poscar_file = os.path.join(jdir, 'POSCAR')
        struct_m = Structure.from_file(poscar_file)
        species = ''.join([tos.symbol for tos in struct_m.types_of_specie])
        if data.get(species):
            for p in params:
                if j.vis.incar.get(p):
                    data[species][p].append([j.vis.incar[p],
                                             j.final_energy / len(struct_m)])
                elif p == 'KPOINTS':
                    data[species]['KPOINTS'].append([j.vis.kpoints.kpts,
                                                     j.final_energy / len(
                                                         struct_m)])
                else:
                    logger.warn(
                        'dont know how to parse the parameter {}'.format(p))
        else:
            data[species] = {}
            for p in params:
                data[species][p] = []
                data[species][p] = []
    return data


def get_opt_params(data, species, param='ENCUT', ev_per_atom=0.001):
    """
    return optimum parameter
    default: 1 meV/atom
    """
    sorted_list = sorted(data[species][param], key=lambda x: x[1])
    sorted_array = np.array(sorted_list)
    consecutive_diff = np.abs(
        sorted_array[:-1, 1] - sorted_array[1:, 1] - ev_per_atom)
    min_index = np.argmin(consecutive_diff)
    return sorted_list[min_index][0]


# PLEASE DONT CHANGE THINGS WITHOUT UPDATING SCRIPTS/MODULES THAT DEPEND
# ON IT
# get_convergence_data and get_opt_params moved to *_custom
def get_convergence_data_custom(jfile, params=('ENCUT', 'KPOINTS')):
    """
    returns data dict in the following format
    {'Al':
          {'ENCUT': [ [500,1.232], [600,0.8798] ],
            'KPOINTS':[ [], [] ]
          },
     'W': ...
    }

    Note: processes only INCAR parmaters and KPOINTS
    Sufficient tagging of the data assumed from species,Poscar
    comment line and potcar functional
    """
    cutoff_jobs = jobs_from_file(jfile)
    data = {}
    for j in cutoff_jobs:
        jdir = os.path.join(j.parent_job_dir, j.job_dir)
        poscar_file = os.path.join(jdir, 'POSCAR')
        struct_m = Structure.from_file(poscar_file)

        species = ''.join([tos.symbol for tos in struct_m.types_of_specie])
        tag = '_'.join([species, Poscar.from_file(poscar_file).comment,
                        j.vis.potcar.functional])
        if data.get(tag):
            for p in params:
                if j.vis.incar.get(p):
                    data[tag][p].append([j.vis.incar[p],
                                         j.final_energy / len(struct_m),
                                         j.vis.potcar, j.vis.poscar])
                    #                print(j.vis.potcar.functional,j.vis.poscar)
                elif p == 'KPOINTS':
                    data[tag]['KPOINTS'].append([j.vis.kpoints.kpts,
                                                 j.final_energy / len(
                                                     struct_m), j.vis.potcar,
                                                 j.vis.poscar])
                else:
                    logger.warn(
                        'dont know how to parse the parameter {}'.format(p))
        else:
            data[tag] = {}
            for p in params:
                data[tag][p] = []
                data[tag][p] = []
    return data


def get_opt_params_custom(data, tag, param='ENCUT', ev_per_atom=1.0):
    """
    Args:
        data:  dictionary of convergence data
        tag:   key to dictionary of convergence dara
        param: parameter to be optimized
        ev_per_atom: minimizing criterion in eV per unit

    Returns
        [list] optimum parameter set consisting of tag, potcar object,
        poscar object, list of convergence data energies sorted according to
        param

    default criterion: 1 meV/atom
    """
    sorted_list = sorted(data[tag][param], key=lambda x: x[0])
    # sorted array data
    t = np.array(sorted_list)[:, 1]
    # print(sorted_array[:-1,1], sorted_array[1:,1], ev_per_atom)
    consecutive_diff = [float(j) - float(i) - ev_per_atom for i, j in
                        zip(t[:-1], t[1:])]
    # print("Consecutive_diff",consecutive_diff)
    min_index = np.argmin(consecutive_diff)
    # return the tag,potcar object, poscar object, incar setting and
    # convergence data for plotting that is optimum
    return [tag, data[tag][param][min_index][2],
            data[tag][param][min_index][3], sorted_list[min_index][0], t]


def partition_jobs(turn_knobs, max_jobs):
    """
    divide turn_knobs into smaller turn_knobs so that each one of
    them has smaller max_jobs jobs
    """
    params_len = [len(v) for k, v in turn_knobs.items()]
    n_total_jobs = reduce(lambda x, y: x * y, params_len)
    partition_size = int(n_total_jobs / max_jobs)
    max_index = np.argmax(params_len)
    max_len = max(params_len)
    max_key = list(turn_knobs.items())[max_index][0]
    partition = range(0, max_len, max(1, int(max_len / partition_size)))
    partition_1 = partition[1:] + [max_len]
    logger.info(
        '{0} list of length {1} will be partitioned into {2} chunks'.format(
            max_key, max_len, len(partition)))
    turn_knobs_list = []
    name_list = []
    for i, j in zip(partition, partition_1):
        ordered_list = []
        for k, v in turn_knobs.items():
            if k == max_key:
                tk_item = (k, v[i:j])
            else:
                tk_item = (k, v)
            ordered_list.append(tk_item)
        turn_knobs_list.append(OrderedDict(ordered_list))
        name_list.append('_'.join([str(i), str(j)]))
    return turn_knobs_list, name_list


def get_logger(log_file_name):
    """
    writes out logging file.
    Very useful project logging, recommended for use
    to monitor the start and completion of steps in the workflow
    Arg:
        log_file_name: name of the log file, log_file_name.log
    """
    loggr = logging.getLogger(log_file_name)
    loggr.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    fh = logging.FileHandler(log_file_name + '.log', mode='a')
    fh.setFormatter(formatter)
    loggr.addHandler(fh)
    return loggr


def set_sd_flags(poscar_input=None, n_layers=2, top=True, bottom=True,
                 poscar_output='POSCAR2'):
    """
    set the relaxation flags for top and bottom layers of interface.
    The upper and lower bounds of the z coordinate are determined
    based on the slab.
    Args:
         poscar_input: input poscar file name
         n_layers: number of layers to be relaxed
         top: whether n_layers from top are be relaxed
         bottom: whether n_layers from bottom are be relaxed
         poscar_output: output poscar file name
    Returns:
         None
         writes the modified poscar file
    """
    poscar1 = Poscar.from_file(poscar_input)
    sd_flags = np.zeros_like(poscar1.structure.frac_coords)
    z_coords = poscar1.structure.frac_coords[:, 2]
    z_lower_bound, z_upper_bound = None, None
    if bottom:
        z_lower_bound = np.unique(z_coords)[n_layers - 1]
        sd_flags[np.where(z_coords <= z_lower_bound)] = np.ones((1, 3))
    if top:
        z_upper_bound = np.unique(z_coords)[-n_layers]
        sd_flags[np.where(z_coords >= z_upper_bound)] = np.ones((1, 3))
    poscar2 = Poscar(poscar1.structure, selective_dynamics=sd_flags.tolist())
    poscar2.write_file(filename=poscar_output)


def print_exception():
    """
    Error exception catching function for debugging
    can be a very useful tool for a developer
    move to utils and activate when debug mode is on
    """
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno,
                                                       line.strip(), exc_obj))


def is_converged(directory):
    """
    Check if a relaxation has converged.

    Args:
        directory (str): path to directory to check.

    Returns:
        boolean. Whether or not the job is converged.
    """

    try:
        return Vasprun('{}/vasprun.xml'.format(directory)).converged
    except:
        return False


def get_spacing(structure):
    """
    Returns the interlayer spacing for a 2D material or slab.

    Args:
        structure (Structure): Structure to check spacing for.
        cut (float): a fractional z-coordinate that must be within
            the vacuum region.

    Returns:
        float. Spacing in Angstroms.
    """

    structure = align_axis(structure)
    structure = center_slab(structure)
    max_height = max([s.coords[2] for s in structure.sites])
    min_height = min([s.coords[2] for s in structure.sites])
    return structure.lattice.c - (max_height - min_height)


def center_slab(structure):
    """
    Centers the atoms in a slab structure around 0.5
    fractional height.

    Args:
        structure (Structure): Structure to center
    Returns:
        Centered Structure object.
    """

    center = np.average([s._fcoords[2] for s in structure.sites])
    translation = (0, 0, 0.5 - center)
    structure.translate_sites(range(len(structure.sites)), translation)
    return structure


def add_vacuum(structure, vacuum):
    """
    Adds padding to a slab or 2D material.

    Args:
        structure (Structure): Structure to add vacuum to
        vacuum (float): Vacuum thickness to add in Angstroms
    Returns:
        Structure object with vacuum added.
    """
    structure = align_axis(structure)
    coords = [s.coords for s in structure.sites]
    species = [s.specie for s in structure.sites]
    lattice = structure.lattice.matrix
    lattice[2][2] += vacuum
    structure = Structure(lattice, species, coords, coords_are_cartesian=True)
    return center_slab(structure)


def ensure_vacuum(structure, vacuum):
    """
    Adds padding to a slab or 2D material until the desired amount
    of vacuum is reached.

    Args:
        structure (Structure): Structure to add vacuum to
        vacuum (float): Final desired vacuum thickness in Angstroms
    Returns:
        Structure object with vacuum added.
    """

    structure = align_axis(structure)
    spacing = get_spacing(structure)
    structure = add_vacuum(structure, vacuum - spacing)
    return center_slab(structure)


def get_rotation_matrix(axis, theta):
    """
    Find the rotation matrix associated with counterclockwise rotation
    about the given axis by theta radians.
    Credit: http://stackoverflow.com/users/190597/unutbu

    Args:
        axis (list): rotation axis of the form [x, y, z]
        theta (float): rotational angle in radians

    Returns:
        array. Rotation matrix.
    """

    axis = np.array(list(axis))
    axis = axis / np.linalg.norm(axis)
    axis *= -np.sin(theta/2.0)
    a = np.cos(theta/2.0)
    b, c, d = tuple(axis.tolist())
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def align_axis(structure, axis='c', direction=(0, 0, 1)):
    """
    Rotates a structure so that the specified axis is along
    the [001] direction. This is useful for adding vacuum, and
    in general for using vasp compiled with no z-axis relaxation.

    Args:
        structure (Structure): Pymatgen Structure object to rotate.
        axis: Axis to be rotated. Can be 'a', 'b', 'c', or a 1x3 vector.
        direction (vector): Final axis to be rotated to.
    Returns:
        structure. Rotated to align axis along direction.
    """

    if axis == 'a':
        axis = structure.lattice._matrix[0]
    elif axis == 'b':
        axis = structure.lattice._matrix[1]
    elif axis == 'c':
        axis = structure.lattice._matrix[2]
    proj_axis = np.cross(axis, direction)
    if not(proj_axis[0] == 0 and proj_axis[1] == 0):
        theta = (
            np.arccos(np.dot(axis, direction)
            / (np.linalg.norm(axis) * np.linalg.norm(direction)))
        )
        R = get_rotation_matrix(proj_axis, theta)
        rotation = SymmOp.from_rotation_and_translation(rotation_matrix=R)
        structure.apply_operation(rotation)
    if axis == 'c' and direction == (0, 0, 1):
        structure.lattice._matrix[2][2] = abs(structure.lattice._matrix[2][2])

    return structure


def get_structure_type(structure, tol=0.1, seed_index=0,
                       write_poscar_from_cluster=False):

    """
    This is a topology-scaling algorithm used to describe the
    periodicity of bonded clusters in a bulk structure.
    Args:
        structure (structure): Pymatgen structure object to classify.
        tol (float): Additional percent of atomic radii to allow
            for overlap, thereby defining bonds
            (0.1 = +10%, -0.1 = -10%)
        seed_index (int): Atom number to start the cluster.
        write_poscar_from_cluster (bool): Set to True to write a
            POSCAR file from the sites in the cluster.
    Returns:
        string. "molecular" (0D), "chain" (1D), "layered" (2D), or
            "conventional" (3D). Also includes " heterogeneous"
            if the cluster's composition is not equal to that
            of the overal structure.
    """

    # Get conventional structure to orthogonalize the lattice as
    # much as possible. A tolerance of 0.1 Angst. was suggested by
    # pymatgen developers.
    s = SpacegroupAnalyzer(structure, 0.1).get_conventional_standard_structure()
    heterogeneous = False

    noble_gases = ["He", "Ne", "Ar", "Kr", "Xe", "Rn"]
    if len([e for e in structure.composition if e.symbol in noble_gases]) != 0:
        type = "noble gas"
    else:
        # make 2x2x2 supercell to ensure sufficient number of atoms
        # for cluster building.
        s.make_supercell(2)

        # Distance matrix (rowA, columnB) shows distance between
        # atoms A and B, taking PBCs into account.
        distance_matrix = s.distance_matrix

        # Fill diagonal with a large number, so the code knows that
        # each atom is not bonded to itself.
        np.fill_diagonal(distance_matrix, 100)

        # Rows (`radii`) and columns (`radiiT`) of radii.
        radii = [ELEMENT_RADII[site.species_string] for site in s.sites]
        radiiT = np.array(radii)[np.newaxis].T
        radii_matrix = radii + radiiT*(1+tol)

        # elements of temp that have value less than 0 are bonded.
        temp = distance_matrix - radii_matrix
        # True (1) is placed where temp < 0, and False (0) where
        # it is not.
        binary_matrix = (temp < 0).astype(int)

        # list of atoms bonded to the seed atom of a cluster
        seed = set((np.where(binary_matrix[seed_index]==1))[0])
        cluster = seed
        NEW = seed
        while True:
            temp_set = set()
            for n in NEW:
                # temp_set will have all atoms, without duplicates,
                # that are connected to all atoms in NEW.
                temp_set.update(set(np.where(binary_matrix[n]==1)[0]))

            if temp_set.issubset(cluster):
                # if temp_set has no new atoms, the search is done.
                break
            else:
                NEW = temp_set - cluster # List of newly discovered atoms
                cluster.update(temp_set) # cluster is updated with new atoms

        if len(cluster) == 0:  # i.e. the cluster is a single atom.
            cluster = [seed_index]  # Make sure it's not empty to write POSCAR.
            type = "molecular"

        elif len(cluster) == len(s.sites): # i.e. all atoms are bonded.
            type = "conventional"

        else:
            cmp = Composition.from_dict(Counter([s[l].specie.name for l in
                                        list(cluster)]))
            if cmp.reduced_formula != s.composition.reduced_formula:
                # i.e. the cluster does not have the same composition
                # as the overall crystal; therefore there are other
                # clusters of varying composition.
                heterogeneous = True

            old_cluster_size = len(cluster)
            # Increase structure to determine whether it is
            # layered or molecular, then perform the same kind
            # of cluster search as before.
            s.make_supercell(2)
            distance_matrix = s.distance_matrix
            np.fill_diagonal(distance_matrix,100)
            radii = [ELEMENT_RADII[site.species_string] for site in s.sites]
            radiiT = np.array(radii)[np.newaxis].T
            radii_matrix = radii + radiiT*(1+tol)
            temp = distance_matrix-radii_matrix
            binary_matrix = (temp < 0).astype(int)

            seed = set((np.where(binary_matrix[seed_index]==1))[0])
            cluster = seed
            NEW = seed
            check = True
            while check:
                temp_set = set()
                for n in NEW:
                    temp_set.update(set(np.where(binary_matrix[n]==1)[0]))

                if temp_set.issubset(cluster):
                    check = False
                else:
                    NEW = temp_set - cluster
                    cluster.update(temp_set)

            if len(cluster) != 4 * old_cluster_size:
                type = "molecular"
            else:
                type = "layered"
    if heterogeneous:
        type += " heterogeneous"

    cluster_sites = [s.sites[n] for n in cluster]
    if write_poscar_from_cluster:
        s.from_sites(cluster_sites).get_primitive_structure().to("POSCAR",
                                                                 "POSCAR")

    return type


def write_potcar(pot_path=VASP_PSP, types='None'):
    """
    Writes a POTCAR file based on a list of types.

    Args:
        pot_path (str): can be changed to override default location
            of POTCAR files.
        types (list): list of same length as number of elements
            containing specifications for the kind of potential
            desired for each element, e.g. ['Na_pv', 'O_s']. If
            left as 'None', uses the defaults in the
            'potcar_symbols.yaml' file in the package root.
    """

    if pot_path == None:
        # This probably means the mpint_config.yaml file has not
        # been set up.
        pass
    else:
        poscar = open('POSCAR', 'r')
        lines = poscar.readlines()
        elements = lines[5].split()
        poscar.close()

        potcar_symbols = loadfn(
            os.path.join(PACKAGE_PATH, 'mat2d', 'potcar_symbols.yaml')
        )

        if types == 'None':
            sorted_types = [potcar_symbols[elt] for elt in elements]
        else:
            sorted_types = []
            for elt in elements:
                for t in types:
                    if t.split('_')[0] == elt:
                        sorted_types.append(t)

        potentials = []

        # Create paths, open files, and write files to
        # POTCAR for each potential.
        for potential in sorted_types:
            potentials.append('{}/{}/POTCAR'.format(pot_path, potential))
        outfile = open('POTCAR', 'w')
        for potential in potentials:
            infile = open(potential)
            for line in infile:
                outfile.write(line)
            infile.close()
        outfile.close()


def write_circle_mesh_kpoints(center=[0, 0, 0], radius=0.1, resolution=20):
    """
    Create a circular mesh of k-points centered around a specific
    k-point and write it to the KPOINTS file. Non-circular meshes
    are not supported, but would be easy to code. All
    k-point weights are set to 1.

    Args:
        center (list): x, y, and z coordinates of mesh center.
            Defaults to Gamma.
        radius (float): Size of the mesh in inverse Angstroms.
        resolution (int): Number of mesh divisions along the
            radius in the 3 primary directions.
    """

    kpoints = []
    step = radius / resolution

    for i in range(-resolution, resolution):
        for j in range(-resolution, resolution):
            if i**2 + j**2 <= resolution**2:
                kpoints.append([str(center[0]+step*i),
                                str(center[1]+step*j), '0', '1'])
    with open('KPOINTS', 'w') as kpts:
        kpts.write('KPOINTS\n{}\ndirect\n'.format(len(kpoints)))
        for kpt in kpoints:
            kpts.write(' '.join(kpt))
            kpts.write('\n')


def get_markovian_path(points):
    """
    Calculates the shortest path connecting an array of 2D
    points. Useful for sorting linemode k-points.

    Args:
        points (list): list/array of points of the format
            [[x_1, y_1, z_1], [x_2, y_2, z_2], ...]

    Returns:
        list: A sorted list of the points in order on the markovian path.
    """

    def dist(x, y):
        return math.hypot(y[0] - x[0], y[1] - x[1])

    paths = [p for p in it.permutations(points)]
    path_distances = [
        sum(map(lambda x: dist(x[0], x[1]), zip(p[:-1], p[1:])))
        for p in paths]
    min_index = np.argmin(path_distances)

    return paths[min_index]


def remove_z_kpoints():
    """
    Strips all linemode k-points from the KPOINTS file that include a
    z-component, since these are not relevant for 2D materials and
    slabs.
    """
    kpoint_file = open('KPOINTS')
    kpoint_lines = kpoint_file.readlines()
    kpoint_file.close()

    twod_kpoints = []
    labels = {}
    i = 4

    while i < len(kpoint_lines):
        kpt_1 = kpoint_lines[i].split()
        kpt_2 = kpoint_lines[i+1].split()
        if float(kpt_1[2]) == 0.0 and [float(kpt_1[0]),
                                       float(kpt_1[1])] not in twod_kpoints:
            twod_kpoints.append([float(kpt_1[0]), float(kpt_1[1])])
            labels[kpt_1[4]] = [float(kpt_1[0]), float(kpt_1[1])]

        if float(kpt_2[2]) == 0.0 and [float(kpt_2[0]),
                                       float(kpt_2[1])] not in twod_kpoints:
            twod_kpoints.append([float(kpt_2[0]), float(kpt_2[1])])
            labels[kpt_2[4]] = [float(kpt_2[0]), float(kpt_2[1])]
        i += 3

    kpath = get_markovian_path(twod_kpoints)

    with open('KPOINTS', 'w') as kpts:
        for line in kpoint_lines[:4]:
            kpts.write(line)

        for i in range(len(kpath)):
            label_1 = [l for l in labels if labels[l] == kpath[i]][0]
            if i == len(kpath) - 1:
                kpt_2 = kpath[0]
                label_2 = [l for l in labels if labels[l] == kpath[0]][0]
            else:
                kpt_2 = kpath[i+1]
                label_2 = [l for l in labels if labels[l] == kpath[i+1]][0]

            kpts.write(' '.join([str(kpath[i][0]), str(kpath[i][1]), '0.0 !',
                                label_1]))
            kpts.write('\n')
            kpts.write(' '.join([str(kpt_2[0]), str(kpt_2[1]), '0.0 !',
                                label_2]))
            kpts.write('\n\n')
    kpts.close()

def update_submission_template(default_template, qtemplate):
    """
    helper function for writing a CommonAdapter template fireworks
    submission file based on a provided default_template which
    contains hpc resource allocation information and the qtemplate
    which is a yaml of commonly modified user arguments
    """
    pass

def check_errors(chkfile=None,logfile_name=None,jdirs=None,handlers=None):
    """
    Manual function to help check errors in job dirs
    of a checkpoint file or in a list of job dirs
    directly
    """
    logf = get_logger(logfile_name)
    if not handlers:
       handlers = [VaspErrorHandler()]

    if chkfile:
        js = update_checkpoint(jfile=chkfile)
        chk_jobs = [j for j in jobs_from_file(chkfile) if not j.final_energy]
        err_chks = [check_job_with_custodian(job=j,logfile=logf,handlers=handlers) for j in chk_jobs]
        return err_chks
    elif jdirs:
        pass

def rerun_jobs(job_cmd=None, chkfiles=None, job_dirs=None):
    """
    Reruns the jobs after performing corrections by Custodian
    reading jobs either from checkpoint file (recommended)
    or job directories.
    """
    logR = get_logger('RerunJobReport')
    if chkfiles:
        checks  = [check_errors(chkfile=f) for f in chkfiles]
        CustodianReport = {'ErrorDir':[c['ErrorDir'][0] for c in sum(checks,[])],
                           'Error':[c['Error'][0] for c in sum(checks,[])],
                           'Correction':[c['Correction'][0] for c in sum(checks,[])]}
        pd.DataFrame(CustodianReport).to_csv('CustodianReport.csv')
        logR.info('Checks result is: {}'.format(checks))
        jobs = sum([jobs_from_file(f) for f in glob('*.json')],[])
        not_done_jobs = [j.job_dir for j in jobs if not j.final_energy]
        logR.info('There are {0} Not done jobs and they are {1}'.format(len(not_done_jobs),not_done_jobs))
        if job_cmd:
            for j in np.unique(not_done_jobs):
                os.chdir(j)
                os.system(job_cmd)
                os.chdir('../')
    if job_dirs:
        pass

def process_to_dataframe(identifiers, metadata = ['energy','volume','kpoints'], chkfiles=glob('*.json'),job_dirs=None):
    """
    utility function that processes data from a list of checkpoints
    or directories into a pandas dataframe

    identifiers: Creates file names for the dataframes
                 1. dict type: job directory concatenation
                example: {0:(('_',0),('/',1)),1:('__',1),'JoinBy':'_','OtherNames':'PBE'}
                         will split the job_dirs elements by
                         '_' and take position 0 followed by a split by '__' take
                         position 1

    """
    if chkfiles:
       chk_jobs = sum([jobs_from_file(c) for c in chkfiles],[])
       energies = [j.final_energy for j in chk_jobs]

       if type(identifiers)==dict:

           def rec_split(string, split_rule):
               print (split_rule)
               for s in split_rule:
                   print (s)
                   string = string.split(s[0])[s[1]]
               return string

           dir_names = [j.job_dir for j in chk_jobs]
           print (list(identifiers.keys()))
           print (identifiers[0])
           dir_splits = [ [rec_split(d,identifiers[k]) for k in list(identifiers.keys()) \
                         if type(k)==int] for d in dir_names]
           dir_identifiers = [ identifiers['JoinBy'].join(d) for d in dir_splits]
           print (dir_identifiers)

def parse_script():
    """
    parser for output post processing
    """
    description = 'Input the checkpoint file or a list of directories enclosed in [..] or as a python glob'
    user_parser = ArgumentParser(description=description)

    user_parser.add_argument('-i', '--input', help='Input json file or a list of directories, see command man')
    user_parser.add_argument('-o', '--output', help='Output data file name in csv')
    if '.json' in sys.argv[2]:
        jdirs = [j.parent_job_dir + os.sep + j.job_dir \
                 for j in jobs_from_file(sys.argv[2]) if j.final_energy]
        print (len(jdirs),'out of',len([j.parent_job_dir + os.sep + j.job_dir \
                                       for j in jobs_from_file(sys.argv[2])]))
    elif 'glob(' in sys.argv[2]:
        jdirs = glob('{}*'.format(sys.argv[2].replace('glob(', '').replace(')','')))
    elif '[' in sys.argv[2]:
        jdirs = sys.argv[2].replace('[','').replace(']','').split(',')
    else:
        jdirs = sys.argv[2]

    output = sys.argv[4]

    return jdirs, output

def parse_material_input():
    """
    general parser for material inputs, takes a source API and 
    a MATERIALS_LIST.txt file. 
    """
    structures = []
    description = 'Input the MATERIALS_LIST.txt and a source. Default is Materials Project API'
    user_parser = ArgumentParser(description=description)

    user_parser.add_argument('-m', '--materials', help='Input MATERIALS_LIST.txt')
    user_parser.add_argument('-a', '--api', help='Source: MPAPI or MWAPI supported')
    if '.txt' in sys.argv[2]:
        materials_file = open(sys.argv[2])
        materials = [m for m in materials_file.readlines()]
        if not materials:
           logger.warn('Materials not read from {}'.format(sys.argv[2]))
    else:
        logger.warn('Please supply an input file ending with .txt')
    if sys.argv[4]:
       api = sys.argv[4]
    else:
       api = 'MPAPI'

    if api == 'MPAPI':
       for m in materials:
          if not 'mp' in m:
             print ('querying {}'.format(m))
             structures.append(get_struct_from_mp(m))
          else:
             pass
    elif api == 'MWAPI':
       for m in materials:
          if 'mw' in m:
             structures.append(get_struct_from_mw(m))
          else:
             pass

    return structures
          
def load_config_vars(config_dict):
    """
    Loads the given dict of config variables to the environment
    """
    if not os.path.exists(SETTINGS_FILE):
       user_configs = {key:None for key in ['username','bulk_binary','twod_binary',\
               'sol_binary','ncl_binary','custom_binary',\
               'vdw_kernel','potentials','MAPI_KEY', 'queue_system', 'queue_template']}

       user_configs['queue_system'] = 'slurm'
       user_configs['queue_template'] = PACKAGE_PATH

       with open(os.path.join(os.path.expanduser('~'),'.mpint_config.yaml'),'w') as config_file:
          yaml.dump(user_configs, config_file, default_flow_style=False)

    current_config = yaml.load(open(SETTINGS_FILE))
    current_config.update(config_dict)
    with open(SETTINGS_FILE,'w') as new_config:
       yaml.dump(current_config, new_config, default_flow_style=False)

    if os.path.exists(PACKAGE_PATH):
       print ('updating the Fireworks template file for the batch submission script based on the \n'\
               'the available qtemplate params in {}'.format(PACKAGE_PATH+'/qtemplate.yaml'))
       fireworks_path = PACKAGE_PATH.replace('/mpinterfaces','/fireworks/user_objects/')
       for f in glob(PACKAGE_PATH+'/*.txt'):
          os.system('cp {0} {1}'.format(f, fireworks_path))
 

def write_pbs_runjob(name, nnodes, nprocessors, pmem, walltime, binary):
    """
    writes a runjob based on a name, nnodes, nprocessors, walltime,
    and binary. Designed for runjobs on the Hennig group_list on
    HiperGator 1 (PBS).

    Args:
        name (str): job name.
        nnodes (int): number of requested nodes.
        nprocessors (int): number of requested processors.
        pmem (str): requested memory including units, e.g. '1600mb'.
        walltime (str): requested wall time, hh:mm:ss e.g. '2:00:00'.
        binary (str): absolute path to binary to run.
    """
    runjob = open('runjob', 'w')
    runjob.write('#!/bin/sh\n')
    runjob.write('#PBS -N {}\n'.format(name))
    runjob.write('#PBS -o test.out\n')
    runjob.write('#PBS -e test.err\n')
    runjob.write('#PBS -r n\n')
    runjob.write('#PBS -l walltime={}\n'.format(walltime))
    runjob.write('#PBS -l nodes={}:ppn={}\n'.format(nnodes, nprocessors))
    runjob.write('#PBS -l pmem={}\n'.format(pmem))
    runjob.write('#PBS -W group_list=hennig\n\n')
    runjob.write('cd $PBS_O_WORKDIR\n\n')
    runjob.write('mpirun {} > job.log\n\n'.format(binary))
    runjob.write('echo \'Done.\'\n')
    runjob.close()


def write_slurm_runjob(name, ntasks, pmem, walltime, binary):
    """
    writes a runjob based on a name, nnodes, nprocessors, walltime, and
    binary. Designed for runjobs on the Hennig group_list on HiperGator
    2 (SLURM).

    Args:
        name (str): job name.
        ntasks (int): total number of requested processors.
        pmem (str): requested memory including units, e.g. '1600mb'.
        walltime (str): requested wall time, hh:mm:ss e.g. '2:00:00'.
        binary (str): absolute path to binary to run.
    """

    nnodes = int(np.ceil(float(ntasks) / 32.0))

    runjob = open('runjob', 'w')
    runjob.write('#!/bin/bash\n')
    runjob.write('#SBATCH --job-name={}\n'.format(name))
    runjob.write('#SBATCH -o out_%j.log\n')
    runjob.write('#SBATCH -e err_%j.log\n')
    runjob.write('#SBATCH --qos=hennig-b\n')
    runjob.write('#SBATCH --nodes={}\n'.format(nnodes))
    runjob.write('#SBATCH --ntasks={}\n'.format(ntasks))
    runjob.write('#SBATCH --mem-per-cpu={}\n'.format(pmem))
    runjob.write('#SBATCH -t {}\n\n'.format(walltime))
    runjob.write('cd $SLURM_SUBMIT_DIR\n\n')
    runjob.write('module load intel/2016.0.109\n')
    runjob.write('module load openmpi/1.10.1\n')
    runjob.write('module load vasp/5.4.1\n\n')
    runjob.write('mpirun {} > job.log\n\n'.format(binary))
    runjob.write('echo \'Done.\'\n')
    runjob.close()

def add_mem_submit_file(orig_file,pmem,qtype='slurm'):
    #print ('adding memory')
    #print (orig_file)
    runjob = open(orig_file,'r')
    orig_lines = [l for l in runjob.readlines() if 'mem-per-cpu' not in l]
    non_queue_lines = [l for l in orig_lines if '#SBATCH' not in l]
    queue_lines = [l for l in orig_lines if '#SBATCH' in l]
    runjob.close()
    with open(orig_file,'w') as runjob:
       runjob.write(non_queue_lines[0])
       for l in queue_lines:
          runjob.write(l)
       if QUEUE_SYSTEM=='slurm':
          runjob.write('#SBATCH --mem-per-cpu={}'.format(pmem))
       elif QUEUE_SYSTEM=='pbs':
          runjob.write('#PBS --mem-per-cpu={}'.format(pmem))
       for l in non_queue_lines[1:]:
          runjob.write(l)
    runjob.close()

def add_walltime(orig_file,walltime):
    print ('adding walltime')
    runjob = open(orig_file,'r')
    job_path = orig_file.replace('/submit_script','')
    shu.copy(job_path+os.sep+'CONTCAR',job_path+os.sep+'POSCAR')
    orig_lines = [l for l in runjob.readlines() if 'walltime' or 'time' not in l]
    runjob.close()
    non_queue_lines = [l for l in orig_lines if '#SBATCH' or '#PBS' not in l]
    queue_lines = [l for l in orig_lines if '#SBATCH' or '#PBS' in l]
    runjob.close()
    with open(orig_file,'w') as runjob:
       runjob.write(non_queue_lines[0])
       for l in queue_lines:
          runjob.write(l)
       if QUEUE_SYSTEM=='slurm':
          runjob.write('#SBATCH --walltime={}'.format(walltime))
       elif QUEUE_SYSTEM=='pbs':
          runjob.write('#PBS --walltime={}'.format(walltime))
       for l in non_queue_lines[1:]:
          runjob.write(l)
    runjob.close()

def decode_log_file(logfile):                
    jdirs = []
    jids = []
    jname = []
    with open(logfile, 'r') as read_file:
       w_read = read_file.readlines()
       for l in w_read:
             if 'running job' in l:
                j_line = l.split(' ')
                w_pos_j  = [n for n,w in enumerate(j_line) if 'job' in w]
                jids.append(j_line[w_pos_j[0]+1])
                try:
                   jname.append(j_line[w_pos_j[0]+2])
                   jdirs.append(j_line[w_pos_j[0]+4].replace('\n',''))
                except:
                   jdirs.append(j_line[w_pos_j[0]+3].replace('\n',''))
    return len(jids),jids,jdirs,jname 

def decode_q(jid,jdir):
    if QUEUE_SYSTEM == 'slurm':
        try:
            output = str(sp.check_output(['squeue', '--job', jid]))
            state = output.rstrip('\n').split('\n')[-1].split()[-4]
        except:
            state = "00"
        if state == 'R' or state == 'CG':
            oszi = Oszicar(jdir+os.sep+'OSZICAR')
            nsw = Incar.from_file(jdir+os.sep+'INCAR')['NSW']
            return state, oszi, nsw
        elif state == 'PD':
            return 'PD', 'no OSZICAR yet','no nsw'
        else:
            nsw = Incar.from_file(jdir+os.sep+'INCAR')['NSW']
            try:
               oszi = Oszicar(jdir+os.sep+'OSZICAR')
               #nsw = Incar.from_file(jdir+os.sep+'INCAR')['NSW']
               return state, oszi, nsw
            except:
               return state, 'Failed OSZICAR', nsw

def get_defo_structure(struct,strain=0.001,strain_direction='N11'):
    idty = np.zeros((3,3),dtype=float)

    idty[0][0] = 1.0
    idty[1][1] = 1.0
    idty[2][2] = 1.0

    si = {'N11':(0,0),'N22':(1,1),'N33':(2,2)}
    si_c = si[strain_direction]
    idty[si_c[0]][si_c[1]] += strain

    defo_grad = idty
    d = Deformation(defo_grad)
    defo_struct = d.apply_to_structure(struct)

    return Poscar(defo_struct,comment=str(strain).replace('.','_'))


