"""
Objective:

Workflow steps overview:
1. Relax 2D and Bulk Competitors
2. Calculate Stability diagram
3. Prep a dummy IBZKPT
4. Run HSE band structure
5. Run +U calibrations
6. Run AFM calcualtions
7. Run Anisotropy calculations

"""

#### IMPORTS ####

##GENERAL python imports##
import os
import sys
from collections import OrderedDict
import yaml
from glob import glob

## MP imports ##
# Inputs
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.inputs import Potcar, Kpoints
# Error handling
from custodian.vasp.handlers import VaspErrorHandler

## MPInterfaces calibrate engine imports
from mpinterfaces.calibrate import Calibrate
from mpinterfaces.utils import *
from mpinterfaces import *

import numpy as np

# vasp binaries
vasp_config = {'twod_binary': VASP_TWOD_BIN,
               'bulk_binary': VASP_STD_BIN,
               'ncl_binary': VASP_NCL_BIN,
               'sol_binary': VASP_SOL_BIN,
               'custom_binary': VASP_CUSTOM_BIN}


print ('current vasp binaries are: ', vasp_config)
try:
   POTCAR_SPEC = yaml.load(open('potcar_symbols.yaml'))
except:
   try:
      POTCAR_SPEC = yaml.load(open(PACKAGE_PATH+os.sep+'potcar_symbols.yaml'))
   except:
      POTCAR_SPEC = None
      warnings.warn('specify a potcar_spec before starting calculations')

### STEP 1#####
def StepVASP0(my_project, struct_list,order_key=0):
    """
    1. reads a poscar list of a 2D motif(s)
    2. sets the default incar, kpoints and potcar as per mw
       database

    Returns:
        checkpoint Motif_like_2D_relax.json, Motif_like_bulk_relax.json
    """
    WORKFLOWS = my_project['Workflow']
    Workflow_Params = WORKFLOWS['Steps'][order_key]
    Workflow_name = Workflow_Params['NAME']
    job_dir = my_project['NAME'] + Workflow_Params['NAME']
    logger = get_logger(job_dir)
    chkpt = job_dir + '.json'

    incar_dict = my_project['Incar_General']
    incar_init = Incar.from_dict(incar_dict)
    kpoints_init = Kpoints.monkhorst_automatic([18,18,18]) # an initialization to anything sensible

    # Decide the Kpoints density per atom setting here
    turn_knobs= OrderedDict({'POSCAR': struct_list,'KPOINTS':Workflow_Params['KPOINTS'],\
                             'POTCAR_pseudopotential':Workflow_Params['PSEUDOPOTENTIAL']})

    if Workflow_Params['Other_Knobs']:
        #for k in Workflow_Params['Other_Knobs']: 
            #print (type(Workflow_Params['Other_Knobs']))
            turn_knobs.update(Workflow_Params['Other_Knobs'])
    job_bin = vasp_config[Workflow_Params['Queue']['Bin']]
    qdict = Workflow_Params['Queue']
    # Decide general queue settings for all runs in this step
    #print (turn_knobs)
    qadapter, job_cmd = get_run_cmmnd(partition=qdict['Partition'],ntasks=qdict['Ntasks'],\
                        nnodes = qdict['Nnodes'],walltime=qdict['Walltime'],job_bin=job_bin,\
                        mem=qdict['Memory'], job_name=job_dir)

    # run the jobs in this step
    run_cal(turn_knobs, qadapter, job_cmd, job_dir, logger,
            chkpt, incar=incar_init, kpoints=kpoints_init,
            poscar=struct_list[0], magnetism=Workflow_Params['Magnetism'],\
            is_matrix=Workflow_Params['Matrix'],\
            Grid_type=Workflow_Params['Kpt_Grid'])
    return [chkpt]

def StepVASP1(my_project,order_key):
   """
   performs a general VASP continuation step with the following actions:
   1. Uses relaxed CONTCARS of the Continue Source chkpt
   2. Updates the INCAR if so specified in the yaml file
   3. Updates the KPOINTS generation if so specified in the yaml
   4. Copies the necessary re-use files
   """
   WORKFLOWS = my_project['Workflow']
   Workflow_Params = WORKFLOWS['Steps'][order_key]
   Workflow_name = Workflow_Params['NAME']
   job_dir = my_project['NAME'] + Workflow_Params['NAME']
   logger= get_logger(job_dir)
   chkpt = job_dir + '.json'
   prev_filter = Workflow_Params['Continue']['Filter']
   prev_chkpt = Workflow_Params['Continue']['Source']

   rerun_paths = continue_job_inputs(chkpt_files= prev_chkpt,\
    user_filters=prev_filter)
   reuse = Workflow_Params['Reuse']

   #incar_init = Incar.from_file( rerun_paths[0]+os.sep+'INCAR' )
   kpoints_init = Kpoints.from_file( rerun_paths[0]+os.sep+'KPOINTS' )
   poscar_init = Poscar.from_file( rerun_paths[0]+os.sep+'POSCAR' )

   # Example update original INCAR (why? preserve settings like cutoff for this
   # step)
   #if Workflow_Params['Incar_Update']:
   incar_dict = Workflow_Params['Incar_Update']
   reuse_incar = Workflow_Params['Reuse_Incar']
   incar = Incar.from_dict(incar_dict)
   incar_remove = Workflow_Params['Incar_Remove']
   #print (incar)
   #else:
   #   incar = Incar.from_file( rerun_paths[0]+os.sep+'INCAR' )

   # Update for Kpoints to be taken care of in run_cal's Grid_type

   turn_knobs = {'POSCAR': rerun_paths}

   if Workflow_Params['Other_Knobs']:
            #print (Workflow_Params['Other_Knobs'])
            turn_knobs.update(Workflow_Params['Other_Knobs'])
            is_matrix = True
   else:
       is_matrix = False

   #print (turn_knobs)

   job_bin = vasp_config[Workflow_Params['Queue']['Bin']]
   qdict = Workflow_Params['Queue']
   # Decide general queue settings for all runs in this step
   qadapter, job_cmd = get_run_cmmnd(partition=qdict['Partition'],ntasks=qdict['Ntasks'],\
                        nnodes = qdict['Nnodes'],walltime=qdict['Walltime'],job_bin=job_bin,\
                        mem=qdict['Memory'], job_name=job_dir)
   if isinstance(Workflow_Params['Reuse'],list):
      reuse = Workflow_Params['Reuse']
   else:
      reuse = True
   #print ('makes to run_cal')
   run_cal(turn_knobs, qadapter, job_cmd, job_dir, logger,
            chkpt, incar=incar, kpoints=kpoints_init, poscar=poscar_init,
            reuse=reuse, reuse_incar=reuse_incar, incar_remove=incar_remove,\
            magnetism=Workflow_Params['Magnetism'],is_matrix=is_matrix, \
            Grid_type=Workflow_Params['Kpt_Grid'])

   return [chkpt]

def Non_VASP_Script(my_project):
   """
   performs an electronic structure prep calcaultion
   generating the optimized IBZKPT for HSE
   """

   WORKFLOWS = my_project['Workflow']
   Workflow_Params = WORKFLOWS['Steps'][2]
   Workflow_name = Workflow_Params['NAME']
   job_dir = my_project['NAME'] + Workflow_Params['NAME']
   chkpt = job_dir + '.json'
   prev_filter = Workflow_Params['Continue']['Filter']
   prev_chkpt = Workflow_Params['Continue']['Source']
   Script = Workflow_Params['Script']
   executable = Script['Executable']
   non_arg_inputs = Script['NonArgInput']
   arg_inputs = Script['ArgInput']

   rerun_paths = continue_job_inputs(chkpt_files= prev_chkpt,\
    user_filters=prev_filter)

   # Run the script now at the rerun_paths
   for r in rerun_paths:
       if inputs:
          shutil.copy(inputs, r)
       os.chdir(r)
       print ('Running {0} in {1}'.format(executable, r))
       script_output = sp.run([executable]+ arg_inputs, stdout=sp.PIPE).stdout.decode('utf-8')
       

   return None

def process_to_dataframe(chkpts, tags=['final_energy', 'job_name']):
    """
    general post-processing tool that will
    make use of MPInterfaces measurement module
    taking a list of checkpoint files
    and performing MeasurementAnalysis tasks on the
    concerned files as dictated by the tags
    Args:
        chkpts: (list) of checkpoint files
        tags : (list) of properties of interest
               Eg: final_energy (default), job_name,
                   band_gap
                   effective_mass, magnetization_orbs
                   magnetization_total
    """

    pass

def run_cal(turn_knobs, qadapter, job_cmd, job_dir, logger, checkpoint_file,incar_remove=None,
            incar=None, poscar=None, potcar=None, kpoints=None, reuse=None,
            Grid_type='G',pseudopotential='PBE',is_matrix=True, cal_type='Cal', reuse_incar = None,
            n_layers=None, magnetism=False, database=None, mappings_override=POTCAR_SPEC):
    """
    calibrate job launching function
    """
    cal = Calibrate(incar, poscar, potcar, kpoints, reuse=reuse,
                    reuse_incar = reuse_incar, is_matrix=is_matrix,
                    turn_knobs=turn_knobs, qadapter=qadapter,
                    job_cmd = job_cmd, job_dir=job_dir, reuse_override = False,
                    Grid_type=Grid_type,pseudopotential=pseudopotential,magnetism=magnetism,
                    checkpoint_file=checkpoint_file, cal_logger=logger, database=database,
                    mappings_override=mappings_override,incar_remove=incar_remove)
    cal.setup()
    cal.run()

if __name__=='__main__':
    """
    Workflow interface with input yaml and
    construction of sequence
    """

    ## Most General inputs for the project from
    # project.yaml

    my_project = yaml.load(open('my_project.yaml')) ## this will be the only CLI input

    NAME = my_project['NAME']

    INCAR_GENERAL = my_project['Incar_General']
    POTCAR_SPEC = yaml.load(open(my_project['Potcar_Spec']))

    MATERIALS_LIST = my_project['Insilico_Fab']['Material_List']
    struct_list = [Poscar.from_file(poscar) for poscar in glob('StructsDir/POSCAR*') \
                   if 'StructsDir' in MATERIALS_LIST] + \
                  [Poscar(get_struct_from_mp(p)) for p in MATERIALS_LIST \
                   if 'StructsDir' not in p]

    WORKFLOWS = my_project['Workflow']
    
    project_log=get_logger(NAME+"_InSilico_Materials")

    # general structure creation manipulation module , for example make slabs
    # make 2D material, or use a different source like GASP as well
    # insilico_fab(creator = SYNTH['creator'], script=SYNTH['script'])

    error_handler = [VaspErrorHandler()]
    steps= my_project['Workflow']['Steps'].keys() # testing [Relax]
    #print (steps)
    #print (struct_list)
    #print ('Reached past steps')
    Relax()

    #launch_daemon([Relax], interval=30,handlers=error_handler, ld_logger=project_log)
