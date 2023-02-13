import pandas as pd
from sqlalchemy import Table, select
import gc
import yaml
from clinpy.utils.utils import dict_to_table

# TODO create table as

def modify_strand(df):
    if df["strand"]==0: #undefined
        return "."
    elif df["strand"]==1: #+
        return "+"
    elif df["strand"]==2: #-
        return "-"
    else:
        raise ValueError("unknown value in strand column")

def import_temp_junction(file, read_fun, read_fun_params, samplename, project, temp_table,
                         annotated_introns, min_junc_reads=10):
    j = read_fun(file, **read_fun_params)
    j["samplename"] = samplename
    # star annotated is useless because it is actually not the annotated but annotated+detected in
    # first pass
    j.columns = ["chrom", "start", "end", "strand", "motif", "annotated",
                  "uniq_map", "multi_map", "max_ohang", "samplename"]
    j = j.drop(columns=["max_ohang", "motif", "annotated"])
    j["strand"] = j.apply(modify_strand, axis=1)
    j = j[(j.uniq_map >= min_junc_reads) & (j.strand != ".")]

    j = j.merge(annotated_introns,
                how="left", on=["chrom", "start", "end", "strand"])

    mask = pd.isna(j.loc[:, "annotated"])
    j.loc[mask, "annotated"]=False
    j.to_sql(temp_table, project.db, if_exists="append", index=False)

def add_to_junction_tables(project, junc_table, mapping_table, temp_table):
    project.metadata.reflect()
    table = junc_table
    mapping = mapping_table
    temp = temp_table
    junc_table = Table(table, project.metadata, autoload=True, autoload_with=project.db)
    junc_temp = Table(temp, project.metadata, autoload=True, autoload_with=project.db)

    new_juncs = select(junc_temp.c.chrom, junc_temp.c.start, junc_temp.c.end,
                       junc_temp.c.strand, junc_temp.c.annotated).distinct().join(junc_table,
                                                           (junc_temp.c.chrom == junc_table.c.chrom) &
                                                           (junc_temp.c.start == junc_table.c.start) &
                                                           (junc_temp.c.end == junc_table.c.end) &
                                                           (junc_temp.c.strand == junc_table.c.strand),
                                                                                  isouter=True).\
        add_columns(junc_table.c.id).filter(junc_table.c.id == None)
    new_juncs = project.session.execute(new_juncs).fetchall()

    if len(new_juncs)>0:
        new_juncs = pd.DataFrame(new_juncs)
        new_juncs.columns = ["chrom", "start", "end", "strand", "annotated", "id"]
        new_juncs = new_juncs.drop(columns=["id"])
        new_juncs.to_sql(table, project.db, if_exists="append", index=False)

    query = select(junc_temp.c.samplename, junc_temp.c.uniq_map,
                   junc_temp.c.multi_map).join(junc_table, (junc_temp.c.chrom == junc_table.c.chrom) &
                                               (junc_temp.c.start == junc_table.c.start) &
                                               (junc_temp.c.end == junc_table.c.end) &
                                               (junc_temp.c.strand == junc_table.c.strand)). \
        add_columns(junc_table.c.id)
    junction_mapping = pd.DataFrame(project.db.execute(query).fetchall())
    junction_mapping.columns = ["samplename", "uniq_map", "multi_map", "junction"]
    junction_mapping.to_sql(mapping, project.db, index=False, if_exists="append")
    del (junction_mapping)
    gc.collect()

    project.db.execute("drop table if exists {}".format(temp))



def create_tables(params, project, create=True):
    tables = []
    for tablename in params.keys():
        tab=dict_to_table(params[tablename], tablename, project.metadata)
        id_cols = [str(col.name) for col in tab._columns if col.primary_key is True]
        tables.append([tablename, tab, id_cols])
        if create:
            tab.create()

    return tables
def import_data(file, project, meta_read_fun, read_fun, assay_params, create_assay=True):
    with open(assay_params["config"]) as y:
        config = yaml.safe_load(y)

    tables=create_tables(config["tables"], project, create_assay)
    mapping_file=meta_read_fun(file, **config["meta_read_fun_params"])

    sample_col=config["sample_col"]
    junction_col=config["junction_col"]

    for index, row in mapping_file.iterrows():
        # add each junction to the temp table
        import_temp_junction(row[junction_col], read_fun, config["read_fun_params"],
                             row[sample_col], project, config["temp_table"],
                             assay_params["annotated_introns"], config["min_junc_reads"])
        #this is here to reduce memory footrpint at the expense of runtime
        add_to_junction_tables(project, tables[0][0], tables[1][0], config["temp_table"])
