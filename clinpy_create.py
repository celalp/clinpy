# this script will create the necessary files folders theses are

# manage.py this is similar to django manage.py it takes configurations that are related to the project
# assays folder it will have 3 files with at 4th optional one per assay import/methods/tables
# need to find a way to use exisiting assays w/o copying all the stuffs

# utils for general table creation will NOT be here that will be used by manage.py though

#! python3

import argparse as arg
from datetime import datetime
import yaml
import sys
from config import *

from clinpy.utils.utils import dict_to_table, dict_to_engine
from clinpy import  project

from sqlalchemy import MetaData
from sqlalchemy.orm import Session

parser = arg.ArgumentParser(description='add to a project database with genome annotations junctions'
                                            'and expression data')
parser.add_argument('-c', '--config', help="config.py file see readme for details", type=str, action="store",
                    default="config.py")
args = parser.parse_args()

#TODO import config

db=dict_to_engine(database)

for assay in data.keys()