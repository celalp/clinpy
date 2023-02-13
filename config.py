# this module is used when a project is created along with manage py

# import your assays here
# assay do not have to come from clinpy as long as you can import them it does not matter
import pandas as pd
from clinpy.assays import base, junctions, expression, variants
from pytxdb import Genome
from clinpy.utils.utils import dict_to_engine

def get_annotated_introns(genome):
    introns=genome.introns(df=True)
    introns=introns[["chrom", "start", "end", "strand"]].copy()
    introns=introns.drop_duplicates()
    introns["annotated"]=True
    return introns

# database description
#TODO use an env file for this
database = {
    "name": "database.db",
    "dbtype": "sqlite",  # or postgresql, mariadb, or mysql, specify the name of the database (not file)
    "create": True, # otherwise append to existing tables this is important for assays with linker tables
    "user": None,
    "pwd" : None,
    "host": None,
    "port" : None
}

# create genome instance here
genome_params={
    "db":dict_to_engine(dbtype=None, name=None, username=None, pwd=None, host=None, port=None), #this is going to be a sqlalchemy connection
    "fasta":None,
    "mart":None
}

genome=Genome(**genome_params)

#TODO update this
data = {
    # this is mandatory you can only change the file name but the structure of the file needs to match the example
    "base": {
        "filename": "base_assay_tables.xlsx",
        "read_fun": pd.read_excel,
        "read_fun_params": {},
        "assay": base,
        "assay_params":{"config":"base.yaml"},  # you can pass arbitrary parameters with an assay
    },
    "junctions": {
        "filename": "junctions.csv",
        "read_fun": pd.read_csv,
        "read_fun_params": {"header": 0, "sep": ","},
        "assay": junctions,
        "sample_col": "sample_id",
        "assay_params": {"min_junc_reads":10,
                         "annotated_introns":get_annotated_introns(genome)},
    },
    "rna_expression": {
        "filename": "expression.csv",
        "read_fun": pd.read_csv,
        "read_fun_params": {"header": 0, "sep": "\t"},
        "assay": expression,
        "sample_col": "sample_id",
        "assay_params": {}
    },
    "rna_variants": {
        "filename": "rna_variants.csv",
        "read_fun": pd.read_csv,
        "read_fun_params": {"header": 0, "sep": "\t"},
        "assay": variants,
        "sample_col": "sample_id",
        "assay_params": {"prefix":"rna", "config": "vcf.yaml"}
    },
    # you can use the same assay instance multiple times if you want
    "dna_variants": {
        "filename": "dna_variants.csv",
        "read_fun": pd.read_csv,
        "read_fun_params": {"header": 0, "sep": "\t"},
        "assay": variants,
        "sample_col": "sample_id",
        "assay_params": {"config": "vcf.yaml"}
    }
}
