import pandas as pd
from clinpy.utils.utils import dict_to_table
import yaml

#TODO docstings

# here params is only table dict
def create_tables(params, project, create=True):
    tables = []
    for tablename in params.keys():
        tab=dict_to_table(params[tablename], tablename, project.metadata)
        id_cols = [str(col.name) for col in tab._columns if col.primary_key is True]
        tables.append([tablename, tab, id_cols])
        if create:
            tab.create()

    return tables

def import_expression(file, read_fun, samplename, engine, read_fun_params, gene=True,):
    dat = read_fun(file, **read_fun_params)
    dat["samplename"] = samplename
    if gene:
        dat = dat.drop(columns=["length", "effective_length", "transcript_id(s)"])
        dat.columns = ["gene", "expected_count", "tpm", "fpkm", "samplename"]
        dat.to_sql("gene_expression", engine, if_exists="append", index=False)
    else:
        dat = pd.read_csv(file, header=0, sep="\t")
        dat["samplename"] = samplename
        dat = dat.drop(columns=["length", "effective_length", "gene_id"])
        dat.columns = ["transcript", "expected_count", "tpm", "fpkm", "isopct", "samplename"]
        dat.to_sql("transcript_expression", engine, if_exists="append", index=False)


def import_data(file, project, meta_read_fun, read_fun, assay_params, create_assay=True):
    with open(assay_params["config"]) as y:
        config=yaml.safe_load(y)

    tables=create_tables(config["tables"], project, create_assay)

    mapping_file=meta_read_fun(file, **config["meta_read_fun_params"])
    
    if create_assay:
        #tables=create_tables(config["tables"], project, create_assay)
        print("first time create table")
    else:
        sample_table = Table('gene_expression', project.metadata, autoload=True, autoload_with=project.db)#may need to change
        query_existed_sample = select(sample_table.c.samplename).distinct()#make sure table has sample_id
        existed_sample = [id for id, in project.db.execute(query_existed_sample)]
        mapping_file = mapping_file[~mapping_file['sample_id'].isin(existed_sample)]
        print("add new data to the existed table")
    
    gene_col=config["gene_expression_column"]
    tx_col = config["transcript_expression_column"]
    sample_col=config["sample_col"]

    for index, row in mapping_file.iterrows():
        import_expression(row[gene_col], read_fun, row[sample_col], project.db, gene=True,
                          read_fun_params=config["read_fun_params"])
        import_expression(row[tx_col], read_fun, row[sample_col], project.db, gene=False,
                          read_fun_params=config["read_fun_params"])


