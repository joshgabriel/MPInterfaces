# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from six.moves import range
from argparse import ArgumentParser



def create_parser():
    parser = ArgumentParser(description=m_description)
    
    parser.add_argument('-i', '--input', help="yaml input file")
    parser.add_argument('-t', '--type', help="type of calculation")

    subparsers = parser.add_subparsers(help='command', dest='command')


    project_parser = subparsers.add_parser('make_project', help='make project \
                                                                                  according to project.yaml')

    project_parser.add_argument('i', type=str,help='name of project file')

    cal_parser = subparsers.add_parser('run', help='run the specified job')
    cal_parser.set_defaults(func=self.run)
    update_parser = subparsers.add_parser('update', help='update/rerun the checkpoint file calibrate.json ')
    update_parser.add_argument('jids', type=str, nargs='*', help='list of job ids')
    update_parser.set_defaults(func=self.update)

    return parser;

