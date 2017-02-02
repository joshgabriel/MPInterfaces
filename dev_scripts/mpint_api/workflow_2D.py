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
from twod_materials.utils import *
# Other helper TwodMaterials functions for calibrate inputs 
from twod_materials.stability.startup import get_competing_phases

## MPInterfaces measurement imports 
from mpinterfaces.measurement import MeasurementAnalysis 
# twod materials plotting commands 
from twod_materials.electronic_structure.analysis import * 

## General Data analysis tools 
import numpy as np
import pandas as pd

# logger

## GLOBAL INCAR and KPOINTS settings##
## Do we need these globals all the time ? 
# potcar , que , incar, kpoints 

# all these files must be copied into the project directory 
# will be taken care of by the mpint CLI 


incar_general = Incar.from_dict(yaml.load(open(INCAR_GENERAL)))
kpoints_general = Kpoints.from_dict(yaml.load(opne(KPOINTS_GENERAL)))
potcar_dict = yaml.load(open(POTCAR_SPEC))

# vasp binaries
config = yaml.load(open('config_mine.yaml')) 
twod_bin = config['twod_binary']
bulk_bin = config['bulk_binary']
ncl_bin = config['ncl_binary']
sol_bin = config['sol_binary']
custom_bin = config['custom_binary']

## this function can possibly be a utils function 
def insilico_fab(creator, inputs='script'):
    """
    General function that takes care of all structural 
    input creation. Should give support for MPInterfaces
    Interface module, the TSA algorithm and possibly GASP
    and path to any python script that may generate 
    poscar files in a customized way that can't be generalized

    Args:
 
       creator : (str) path to structure creating script OR
                  'Interface' to make use of 
                  mpinterfaces.interface once it is extended to 
                  command all generalizable interface creation 
                  OR directory of structure files ready to run
 
       inputs  : (str or dict) 
                  string option is only 'script' meaning creator is 
                  an independently customized script. This is 
                  the default and the inputs to the script 
                  must be taken care of independently by the user 

                  NOTE: This script must write out the structure files 
                  in the default poscar format into the 
                  directory STRUCTURE_INPUTS

                  STRUCTURE_INPUTS must be named one of 
                  BulkDir, TwodDir, LigandDir, SlabDir, InterfaceDir 

                  dict depends on the required inputs by the creation 
                  algorithm being invoked.

                  Eg: if creator == 'Interface':

                         inputs = {'Bulk': {}, 'Ligand': {}, 'Slab': {}, 
                        'Heterostructure': {}, 'TwoD': {} }

                       Example for Bulk dict: 
                        {'Bulk': {'BulkDir': '/path/to/Bulk_Source'} }
                         which will use the BulkDir path 
                         and append to bulk the POSCAR*.vasp files
                         OR 
                         {'Bulk': {'ChemMP': 'PbS', 'all_structs': True} }
                         pulls all structures from MP using REST API or 
                         alternatively pulls single structures lowest hull 

                       Example for Slab dict: 
                        make a slab of PbS in 100 direction for 
                        15 A thickness and vacuum each 
                        {'Slab': {'bulk': {'ChemMP': 'PbS'},
                                  'hkl': [1,0,0]
                                  'min_vacuum': 15,
                                  'min_thick': 15} }

                        NOTE: The above input for inputs can be extended 
                              for multiple interfaces of a mixture of types
                              either by multiple calls to the function 
                              OR .. inputs = {'Interface1':{..}, 'Interface2':{..} }
    """

    if inputs == 'script':
        
       project_log.info("assuming custom script {} generating structures \
              into  STRUCTURE_INPUTS".format(creator))
       project_log.info("executing {} ... ".format(creator))

       try:
           os.system('python {}'.format(creator))
           if os.isdir('BulkDir') or os.isdir('TwodDir') or os.isdir('LigandDir')\
              os.isdir('InterfaceDir') or os.isdir('SlabDir'):
              project_log.info("script ran successfully! ".format(creator))
              project_log.warn("Please check if the structures are as expected.. ")
           else:
              project_log.info("wrong directroy destination possibly \
              ... exiting ".format(creator))
              sys.exit()
       except:
           project_log.warn('Script failed ..  exiting ..')
           sys.exit()

    pass

# this function could also be a utils function 
def continue_job_inputs(chkpt_files=None, old_job_dir=None, user_filters=None,\
                     must_converge=True, run_restarts_also=False):
   """
   utility function for any workflow that can be used to decide which 
   materials pass to the next stage of a workflow 
   """
   if chkpt_files:
   # Fuel source from checkpoint files 
      rerun_paths = [] 
      chkpts = sum( [ jobs_from_file(j) for j in\
                    [chks for chks in chkpt_files], [] )

      all_converged = [ (j.parent_job_dir+os.sep+j.job_dir,\
        j.job_dir.split('/')[-1].split('_')[0],\
        j.job_dir.split('/')[-1].split('_')[-1])\
        for j in chkpts if j.final_energy ]

      restarts = [ (j.parent_job_dir+os.sep+j.job_dir,\
        j.job_dir.split('/')[-1].split('_')[0],\
        j.job_dir.split('/')[-1].split('_')[-1])\
        for j in chkpts if not j.final_energy ]

      if user_filters and must_converge:

           filtered_run = yaml.load(open('To_Run.yaml'))
           for (job_path, compound, structure) in all_converged:
              for f in filtered_run.keys():
                 for s in filtered_run[f]:
                    if f == compound and s == structure:
                       rerun_paths.append(job_path)
      else:
          rerun_paths = [p[0] for p in all_converged] 

   elif old_job_dir:
       # Fuel from old directory structure 
       rerun_paths = old_job_dirs

   if run_restarts_also:
       return rerun_paths, restarts
   else:
       return rerun_paths


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

### STEP 1#####   
def Relax(**kwargs):
    """
    1. reads a poscar list of a 2D motif(s)
    2. sets the default incar, kpoints and potcar as per mw
       database
    
    Returns:
        checkpoint Motif_like_2D_relax.json, Motif_like_bulk_relax.json
    """
    job_dir = NAME+'_TwoD_Relax'
    incar_dict = incar_general

    incar_dict['SYSTEM']= 'TwoD'
    incar_init = Incar.from_dict(incar_dict)
    kpoints_init = kpoints_general
   
    twod_list = [p for p in glob('TwodDir/POSCAR*.vasp') ] 
    poscar_init = twod_list[0]
    turn_knobs= {'POSCAR': twod_list}

    chkpt = NAME+'_TwoD.json'
    job_bin = twod_bin
    # feature request for naming of job_name according to the turn_knobs 
    # poscar also this is the only fireworks dependency in production 
    qadapter, job_cmd = get_run_cmmnd(ntasks=16, nnodes = 1,\
        walltime='00:02:00',job_bin=job_bin, mem=1800, job_name='TwoD')

    run_cal(turn_knobs, qadapter, job_cmd, job_dir,
            chkpt, incar=incar_init, kpoints=kpoints_init, 
            poscar=poscar_init, is_matrix=False, Grid_type='relax_2D', 
            database='twod')

    job_dir = NAME+'_Competitors_Relax'
    bulk_list = [p for p in glob('BulkDir/POSCAR*.vasp') ]
    poscar_init = bulk_list[0]
    turn_knobs_competitors = {'POSCAR': bulk_list}
    competitors_chkpt = NAME+'_Competitors.json'

    job_bin = bulk_bin
    run_cal(turn_knobs, qadapter, job_cmd, job_dir,
            chkpt, incar=incar_init, kpoints=kpoints_init, 
            poscar=poscar_init, is_matrix=False, Grid_type='relax_3D', 
            database='twod')

    return [twoD_chkpt]


def calculate_stability(**kwargs):
   """
   plots the publication quality formation energy diagram 
   using as inputs the two_relax.json and competitors.json 
   same checkpoints can be used to submit to the database 
   """
   pass


def HSE_prep(**kwargs):
   """
   performs an electronic structure prep calcaultion 
   generating the optimized IBZKPT for HSE 
   """

   rerun_paths = continue_job_inputs(chkpt_files=glob('chkpts/*.json'),\
    user_filters='To_Run.yaml') 

   incar_init = Incar.from_file( rerun_paths[0]+os.sep+'INCAR' ) 
   kpoints_init = Kpoints.from_file( rerun_paths[0]+os.sep+'KPOINTS' )
   poscar_init = Poscar.from_file( rerun_paths[0]+os.sep+'POSCAR' )

   chkpt = NAME+'_hse_prep.json'
   turn_knobs = {'POSCAR': rerun_paths}
   job_dir = NAME+"_HSE"
 
   job_bin = twod_bin
   qadapter, job_cmd = get_run_cmmnd(ntasks=16, walltime='00:02:00',\
                         job_bin=job_bin, mem=1800, job_name='pHSE_2')

   run_cal(turn_knobs, qadapter, job_cmd, job_dir,
            chkpt, incar=incar_init, kpoints=kpoints_init, poscar=poscar_init, 
            reuse=True, is_matrix=False, Grid_type='hse_bands_2D_prep', 
            database='twod', mappings_override = potcar_dict )
   
   return [prep_chkpt] 


def HSE(**kwargs):
   """
   runs the actual HSE electronic structure run 
   """

   # Reuse fuel for this step is always the hse_prep calculation 
   rerun_paths = continue_job_inputs(chkpt_files=[NAME+'_HSE_elec_prep_afm.json'])

   incar_init = Incar.from_file(reuse_list[0]+os.sep+'INCAR')
   kpoints_init = Kpoints.from_file(reuse_list[0]+os.sep+'KPOINTS')
   poscar_init = Poscar.from_file(reuse_list[0]+os.sep+'POSCAR') 

   job_dir = NAME+"_HSE_AFM"
   chkpt = NAME+'_hse_electronics.json' 
   # TODO can additionally measure out information on bandgap 
   turn_knobs_2D = {'POSCAR': rerun_paths, 'NSW':[0]} 

   job_bin = twod_bin
   qadapter, job_cmd = get_run_cmmnd(nnodes = 2, ntasks=64,\
   walltime='72:00:00', job_bin=job_bin, mem=1800, job_name='HSE_AFM')

   # is_matrix set to True here to ease the directory structure 
   run_cal(turn_knobs, qadapter, job_cmd, job_dir,\
           chkpt, incar=incar_init, kpoints=kpoints_init, magnetism = 'AFM',
           poscar=poscar_init, is_matrix=True, Grid_type='hse_bands_2D',\
           reuse=True, database='twod',mappings_override = potcar_dict)

   return [electronic_chkpt]


def AFM(**kwargs):
   """
   runs AFM solution on the structures 
   """
   rerun_paths = continue_job_inputs(old_job_dirs= [NAME+'_AFM/POS/VN2_599_VN2_2H'])

   incar_init = Incar.from_file( rerun_paths[0]+os.sep+'INCAR'  )
   kpoints_init = Kpoints.from_file( rerun_paths[0]+os.sep+'KPOINTS' )
   poscar_init = Poscar.from_file( rerun_paths[0]+os.sep+'CONTCAR' )

   job_dir = NAME+"_HSE_prep_AFM"
   chkpt = NAME+'_HSE_elec_afm.json'
   turn_knobs = {'POSCAR': rerun_paths, 'NSW':[0]}

   job_bin = twod_binary 
   qadapter, job_cmd = get_run_cmmnd(ntasks=16, walltime='10:00:00',\
                          job_bin=job_bin, mem=1800, job_name='hse_prep')


   run_cal(turn_knobs, qadapter, job_cmd, job_dir,
            chkpt, incar=incar, kpoints=kpoints, poscar=poscar_init,
            is_matrix=True, database='twod', magnetism = 'AFM',
            Grid_type='hse_bands_2D', reuse=True,
            mappings_override = potcar_dict)
   

def U_calib(**kwargs):
   """
   lda calibration which turn_knobs over +U 
   """

   rerun_paths = continue_job_inputs(chkpt_files=glob('chkpts/*.json'),\
                                     user_filters='To_Run.yaml')

   incar_init = Incar.from_file( rerun_paths[0]+os.sep+'INCAR'  )
   kpoints_init = Kpoints.from_file( rerun_paths[0]+os.sep+'KPOINTS' )
   poscar_init = Poscar.from_file( rerun_paths[0]+os.sep+'POSCAR' )

   # Example of updating an old incar 
   job_dir = NAME+"_vdWU_cal"
   chkpt = NAME+_'vdWU_calibration.json'

   LDA_Settings = {'LDAU': True,
                   'LDAUTYPE':2,
                   'LDAUL': [2,-1],
                   'LDAUPRINT':2,
                   'LDAUU': [2.5, 0.0],
                   'NSW': 50,
                   'LCHARG':True,
                   'LWAVE':False }
   incar_dict = incar.as_dict()
   ## Remove old settings -- check if this is still required 
   #vdW_tags = ('GGA', 'AGGAC', 'LUSE_VDW', 'PARAM1', 'PARAM2')
   #for key in vdW_tags:
   #   if key in incar_dict:
   #       del incar_dict[key]
   incar_dict.update(LDA_Settings) 
   incar_init = Incar.from_dict(incar_dict)
   turn_knobs_2D = {'POSCAR': rerun_paths,'POTCAR_functional':['PBE'], 
                    'LDAUU': [[0.5*x, 0.0] for x in range(1,20)]}

   job_bin = twod_bin
   qadapter, job_cmd = get_run_cmmnd(ntasks=16, walltime='15:00:00',
                                 job_bin=job_bin, mem=1800, job_name='vdWU')


   run_cal(turn_knobs_2D, qadapter, job_cmd, job_dir,
            chkpt, incar=incar, kpoints=kpoints, poscar=poscar_init,
            is_matrix=True, database = 'twod', reuse = True, reuse_incar='update',  
            mappings_override = potcar_dict)

   return [chkpt]

def MAE(**kwargs):
   """
   takes relaxed structure, CHGCAR
   updates NBANDS, LSORBIT, SAXIS, MAGMOM. ISIF
   runs MAE as a function of a user defined angle theta or
   the default 0 and 90 degs. 
   """

   rerun_paths = continue_job_inputs(old_job_dirs=\
                 ['VN2_117_VN2_2H__POT__LDA__LDAUU__4_5x0_0'])

   Saxis = [] 
   phis = [] 
   thetas = []

   for theta in range(0, 91, 5):
      rad_theta = theta * np.pi / 180
      for phi in range(0, 91, 5):
        rad_phi = phi * np.pi / 180

        saxis = \
         [np.sin(rad_theta) * np.cos(rad_phi), \
          np.sin(rad_theta) * np.sin(rad_phi),\
          np.cos(rad_theta)]

        Saxis.append(saxis)
        phis.append(phis)
        thetas.append(thetas)
        print "SAXIS: {0} PHI: {1} THETA: {2}".format(saxis, phi, theta)

   print (len(Saxis), len(phis), len(thetas))
   saxis_angle_translation={'Saxis':Saxis, 'phis': phis , 'thetas': thetas} 
   saxis_dats = pd.DataFrame(saxis_angle_translation)
   saxis_dats.to_csv('SAXIS.csv')   

   incar_dict = Incar.from_file(rerun_paths[0]+os.sep+'INCAR').as_dict()
   # remove the dict of vdW entries - already done in calibrate can be removed
   # mostly  
   vdW_tags = ('GGA', 'AGGAC', 'LUSE_VDW', 'PARAM1', 'PARAM2')
   for key in vdW_tags:
       if key in incar_dict:
            del incar_dict[key]
   incar_update = {'NSW':0, 'EDIFF': 1e-08, 'ICHARG':11}
   incar_dict.update(incar_update) 
   #############################################################
   incar_init = Incar.from_dict(incar_dict)
   kpoints_init = Kpoints.from_file(rerun_paths[0]+os.sep+'KPOINTS')
   poscar_init = Poscar.from_file( rerun_paths[0]+os.sep+'CONTCAR' )

   job_dir = NAME+"_MAE_Angles"
   chkpt = NAME+'_MAE_LDAU_Angles_Test_VN2.json'

   turn_knobs = {'POSCAR': rerun_paths, 'POTCAR_functional': ['LDA'], 'SAXIS': Saxis}

   job_bin = ncl_bin
   qadapter, job_cmd = get_run_cmmnd(ntasks=16, walltime='03:00:00',
                                 job_bin=job_bin, mem=1800, job_name='MAE_U')

   run_cal(turn_knobs, qadapter, job_cmd, job_dir,
            chkpt, incar=incar, kpoints=kpoints, poscar=poscar_init,
            is_matrix=True, database = 'twod', reuse = ['CHGCAR'], 
            reuse_incar='update',magnetism='MAE', mappings_override = potcar_dict)
   

def Restart_failed_jobs(**kwargs):
    """
    helper function to restart jobs that failed .. can do some 
    error diagnostics also perhaps or treat the case 
    according to a user's prescribed set of diagnostics 
    """
    rerun_paths = continue_job_inputs(run_restarts_also=True)
    pass 

def run_cal(turn_knobs, qadapter, job_cmd, job_dir, checkpoint_file,
            incar=None, poscar=None, potcar=None, kpoints=None, reuse=None,
            Grid_type='G',functional='PBE',is_matrix=True, cal_type='Cal', reuse_incar = None,
            n_layers=None, magnetism=False, database='twod', mappings_override=None):
    """
    calibrate job launching function 
    """
    cal = Calibrate(incar, poscar, potcar, kpoints, reuse=reuse,
                    reuse_incar = reuse_incar, is_matrix=is_matrix, 
                    turn_knobs=turn_knobs, qadapter=qadapter,
                    job_cmd = job_cmd, job_dir=job_dir, reuse_override = False,
                    Grid_type=Grid_type,functional=functional,magnetism=magnetism,
                    checkpoint_file=checkpoint_file, cal_logger=logger, database=database)
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

    NAME = my_project['Name']

    INCAR_GENERAL = my_project['Incar_General']
    KPOINTS_GENERAL = my_project['Kpoints_General']
    POTCAR_SPEC = my_project['Potcar_Spec']

    SYNTH = my_project['In_Silico_Fab']

    project_log=get_logger(NAME+"_InSilico_Materials")
    insilico_fab(creator = SYNTH['creator'], script=SYNTH['script'])

    error_handler = [VaspErrorHandler()]
    steps= my_project['workflow'].keys() # testing [Relax]

    launch_daemon(steps, interval=30,handlers=error_handler, ld_logger=project_log)
  
