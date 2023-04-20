import pandas as pd
from datetime import datetime
from sqlalchemy import Table, select
import gc
import yaml
from clinpy.utils.utils import dict_to_table
import numpy as np

def create_tables(params, project, create=True):
    tables = []
    for tablename in params.keys():
        tab=dict_to_table(params[tablename], tablename, project.metadata)
        id_cols = [str(col.name) for col in tab._columns if col.primary_key is True]
        tables.append([tablename, tab, id_cols])
        if create:
            tab.create()

    return tables

def import_interesting_expr(file, read_fun, samplename, engine, read_fun_params):
    dat = read_fun(file, **read_fun_params)
    #dat = dat[:100] test purpose
    dat["samplename"] = samplename

    #select only true and then drop the outlier column

    dat.columns = ["gene", "TPM_z_score","expression_outlier_z","direction","samplename"]
    dat = dat.loc[dat['expression_outlier_z']]
    dat = dat.drop('expression_outlier_z', axis=1)
    print(dat.shape)
    
    dat.drop_duplicates(subset=['gene'], inplace=True, keep='first')

    dat.to_sql('interesting_expr', engine, if_exists="append", index=False)


def import_data(file, project, meta_read_fun, read_fun, assay_params, create_assay=True):
    with open(assay_params["config"]) as y:
        config = yaml.safe_load(y)
        
    sample_col=config["sample_col"]
    expr_col=config["expr_col"]
    #mapping_file=meta_read_fun(file, **config["meta_read_fun_params"])
    mapping_file=meta_read_fun(file, sheet_name=config['sheetname'], comment="#")
    if create_assay:
        tables=create_tables(config["tables"], project, create_assay)    
        print("first time create table")
    else:
        sample_table = Table('interesting_expr', project.metadata, autoload=True, autoload_with=project.db)#may need to change
        query_existed_sample = select(sample_table.c.samplename).distinct()#make sure table has sample_id
        existed_sample = [id for id, in project.db.execute(query_existed_sample)]
        mapping_file = mapping_file[~mapping_file['sample_id'].isin(existed_sample)]
        if mapping_file.empty:
            print("no update on the {}".format(config['sheetname']))
            return None
        else:
            print("add new data to the existed table")

    for index, row in mapping_file.iterrows():
        print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "starting to add interesting expr for {}".format(row[sample_col]))

        import_interesting_expr(row[expr_col], read_fun, row[sample_col], project.db,
                          read_fun_params=config["read_fun_params"])