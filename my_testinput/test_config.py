# this module is used when a project is created along with manage py

# import your assays here
# assay do not have to come from clinpy as long as you can import them it does not matter
import pandas as pd
from clinpy.assays import base, junctions, expression, variants, interesting_junc
from pytxdb import Genome
from clinpy.utils.utils import dict_to_engine
import pybiomart as biomart


def get_annotated_introns(genome):
    introns=genome.introns(df=True)
    introns=introns[["chrom", "start", "end", "strand"]].copy()
    introns=introns.drop_duplicates()
    introns["annotated"]=True
    return introns

# database description
#TODO use an env file for this
database = {
    "name": "/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/test_clinpy2.db",
    "dbtype": "sqlite",  # or postgresql, mariadb, or mysql, specify the name of the database (not file)
    "create": True, # otherwise append to existing tables this is important for assays with linker tables
    "user": None,
    "pwd" : None,
    "host": None,
    "port" : None
}

# create genome instance here
genome_params={
    "db":dict_to_engine(dbtype='sqlite', name='/hpf/largeprojects/ccmbio/acelik_files/mital/genome.db', username=None, pwd=None, host=None, port=None), #this is going to be a sqlalchemy connection
    "fasta":"/hpf/largeprojects/ccmbio/acelik_files/clin_genomes/hg37/hs37d5_spikein.fa",
    "mart":biomart.Dataset(name='hsapiens_gene_ensembl',host='http://grch37.ensembl.org')
}

genome=Genome(**genome_params)

#TODO update this
data = {
    # this is mandatory you can only change the file name but the structure of the file needs to match the example
    "base": {
        "filename": "/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/input_example/full_meta.xlsx",
        "read_fun": {},
        "assay": base,
        "assay_params":{"config":"/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/input_example/base.yaml"},  # you can pass arbitrary parameters with an assay
    },
    "junctions": {
        "filename": "/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/input_example/full_meta.xlsx",
        "read_fun": pd.read_csv,
        "assay": junctions,
        "sample_col": "sample_id",
        "assay_params": {"config":"/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/input_example/junction.yaml",
                         "min_junc_reads":0,
                         "annotated_introns":get_annotated_introns(genome)},
    },
    "insteresting_junctions":{
        "filename": "/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/input_example/full_meta.xlsx",
        "read_fun": pd.read_csv,
        "assay": interesting_junc,
        "sample_col": "sample_id",
        "assay_params": {"config":"/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/input_example/filter_junc.yaml"}
    },
    "rna_expression": {
        "filename": "/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/input_example/full_meta.xlsx",
        "read_fun": pd.read_csv,
        "assay": expression,
        "sample_col": "sample_id",
        "assay_params": {"config":"/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/input_example/expression.yaml"}
    }  
}
