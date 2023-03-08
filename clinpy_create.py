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

import pandas as pd
import openpyxl

from clinpy.utils.utils import dict_to_table, dict_to_engine
from clinpy import  project

from sqlalchemy import MetaData
from sqlalchemy.orm import Session

if __name__ == "__main__":
    parser = arg.ArgumentParser(description='add to a project database with genome annotations junctions'
                                            'and expression data')
    parser.add_argument('-c', '--config', help="config.py file see readme for details", type=str, action="store",
                    default="config.py")
    args = parser.parse_args()
    
    #TODO import config
    exec(open("config.py").read())
    
    db=dict_to_engine(database['dbtype'],database['name'])
    project=Project(db_params=database)
    
    create=True#should be one of the arguments in the input
    
    #temp setting
    for assay in data.keys():
        if assay == 'base':
        print('starting base')
        data[assay]['assay'].import_data(data[assay]['filename'], project=project,
                                       read_fun=pd.read_excel,read_fun_params={},
                                      assay_params=data[assay]['assay_params'], create_assay=create)
        else:
        print(assay)
        data[assay]['assay'].import_data(data[assay]['filename'], project=project,
                                       meta_read_fun=pd.read_excel,
                                       read_fun=data[assay]['read_fun'],
                                      assay_params=data[assay]['assay_params'], create_assay=create)