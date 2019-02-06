# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from six.moves import range
from argparse import ArgumentParser
import argcomplete

def mpint_parse_arguments(args):
    m_description = """
Management tool for vasp projects, starting from
encut, kpoint or other parameter optimization of till the slab solvation.

it takes 3 arguments: input yaml file, type of calculation and the run mode
example:
    mpint -i naf.yaml -t bulk_calibrate run
this will read in the specifications for 'bulk_calibrate' job from the
input yaml file, naf.yaml, and runs the job i.e submits to the que.

run modes supported:
    1. run : submits job to the que

Everytime jobs are submitted or its sttaus queried, information such as job ids,
job folders etc are written to the log file 'mpint.log'. This makes it easier to
identify job ids and their corresponding to job folders.

Note: use your own materials project key to download the required
structure
"""
    parser = ArgumentParser(description=m_description)
    
    #parser.add_argument('-i', '--input', help="yaml input file")
    #parser.add_argument('-t', '--type', help="type of calculation")
    #parser.add_argument('-m', '--mem',help="new memory to add per cpu")

    subparsers = parser.add_subparsers(help='command', dest='command')
    argcomplete.autocomplete(parser)
    #print ('Here')

    project_parser1 = subparsers.add_parser('start_project', help='start project \
                                                                                  according to project.yaml')

    project_parser1.add_argument('-i', type=str,help='name of project file')

    project_parser2 = subparsers.add_parser('check_project', help='check project with custodian\
                                                                                  according to project.yaml')

    project_parser2.add_argument('-i', type=str,help='name of project file')

    project_parser2.add_argument('-c', type=str,help='checkpoint')

    project_parser3 = subparsers.add_parser('analyze_project', help='analyze/post_process \
                                                                                  according to project.yaml')

    project_parser3.add_argument('-i', type=str,help='name of project file')



    project_parser4 = subparsers.add_parser('rerun_project', help='rerun jobs from CustodianReport.yaml file')

    project_parser4.add_argument('-m', type=str, help="new memory to add per cpu for memory failed jobs")
    project_parser4.add_argument('-w',type=str,help="new walltime to apply for TIMEOUT jobs")
    project_parser4.add_argument('-s',type=str,help="new submit_file to use for all jobs")
       

    project_parser4.add_argument('-i', type=str,help='name of *CustodianReport.yaml file')

    project_parser4.add_argument('-inc', type=str,help='Incar update')
    project_parser4.add_argument('-dinc',type=str,help='Incar remove tags')
    project_parser4.add_argument('-kpt',type=str,help='new kpoint setting tuple of Monkhorst Pack or a pra input')

    project_parser4.add_argument('-c',type=str,help='choose checkpoint')
    

    project_parser5 = subparsers.add_parser('archive_project', help='archive all checkpoints, results, vaspruns \
                                                                                  according to project.yaml')


    project_parser5.add_argument('-i', type=str,help='name of project file')



    project_parser6 = subparsers.add_parser('continue_project', help='continue to the next workflow step in the project.yaml')

    project_parser6.add_argument('-i', type=str,help='name of project file')

    project_parser7 = subparsers.add_parser('load_settings', help='configure the mpint_config.yaml and re-load the new configuration variables for use')
    project_parser7.add_argument('-i', type=str, help='dict of configuration')

    project_parser8 = subparsers.add_parser('qcheck_project',help='check q and report status of project jobs, percent in q, n_electronics done, n_ionics done, any q err msg')

    project_parser8.add_argument('-i', type=str,help='name of project file')

    
    project_parser9 = subparsers.add_parser('cancel_project',help='check q and report status of project jobs, percent in q, n_electronics done, n_ionics done, any q err msg')

    project_parser9.add_argument('-i', type=str,help='name of project file')

#    cal_parser = subparsers.add_parser('run', help='run the specified job')
#    update_parser = subparsers.add_parser('update', help='update/rerun the checkpoint file calibrate.json ')
#    update_parser.add_argument('jids', type=str, nargs='*', help='list of job ids')

    return parser.parse_args(args)#;

