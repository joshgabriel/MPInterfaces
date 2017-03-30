# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from six.moves import range
from argparse import ArgumentParser


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

Note: this script submits jobs only to the PBS ques such as hipergator
"""
    parser = ArgumentParser(description=m_description)
    
    parser.add_argument('-i', '--input', help="yaml input file")
    parser.add_argument('-t', '--type', help="type of calculation")

    subparsers = parser.add_subparsers(help='command', dest='command')


    project_parser = subparsers.add_parser('make_project', help='make project \
                                                                                  according to project.yaml')

    project_parser.add_argument('i', type=str,help='name of project file')

    cal_parser = subparsers.add_parser('run', help='run the specified job')
    update_parser = subparsers.add_parser('update', help='update/rerun the checkpoint file calibrate.json ')
    update_parser.add_argument('jids', type=str, nargs='*', help='list of job ids')

    return parser.parse_args(args);

