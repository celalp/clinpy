import pandas as pd
from sqlalchemy import Table, select
import gc

def modify_strand(df):
    if df["strand"]==0: #undefined
        return "."
    elif df["strand"]==1: #+
        return "+"
    elif df["strand"]==2: #-
        return "-"
    else:
        raise ValueError("unknown value in strand column")

def import_expression(file, samplename, engine, gene=True):
    dat = pd.read_csv(file, header=0, sep="\t")
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


def import_temp_junction(file, samplename, engine, min_junc_reads, filtered=True):
    j = pd.read_csv(file, header=None, sep="\t", low_memory=False)
    j["samplename"] = samplename
    # star annotated is useless because it is actually not the annotated but annotated+detected in
    # first pass
    j.columns = ["chrom", "start", "end", "strand", "motif", "annotated",
                  "uniq_map", "multi_map", "max_ohang", "samplename"]
    j = j.drop(columns=["max_ohang", "motif", "annotated"])
    j["strand"] = j.apply(modify_strand, axis=1)
    j = j[(j.uniq_map >= min_junc_reads) & (j.strand != ".")]
    if not filtered:
        j.to_sql("temp_all_junc", engine, if_exists="append", index=False)
    else:
        j.to_sql("temp_filt_junc", engine, if_exists="append", index=False)

def add_to_junction_tables(engine, meta, session, create=True, filtered=True):
    meta.reflect()
    if filtered:
        table = "junctions"
        mapping = "sample_to_junctions"
        temp = "temp_filt_junc"
        junc_table = Table("junctions", meta, autoload=True, autoload_with=engine)
        samp_to_junc = Table("sample_to_junction", meta, autoload=True, autoload_with=engine)
        junc_temp = Table("temp_filt_junc", meta, autoload=True, autoload_with=engine)
    else:
        table = "all_junctions"
        mapping = "sample_to_alljunctions"
        temp = "temp_all_junc"
        junc_temp = Table("temp_all_junc", meta, autoload=True, autoload_with=engine)
        junc_table = Table("all_junctions", meta, autoload=True, autoload_with=engine)
        samp_to_junc = Table("sample_to_alljunction", meta, autoload=True, autoload_with=engine)


    if create:
        # this means the all_junctions and junction tables are empty
        distinct_junc = select(junc_temp.c.chrom, junc_temp.c.start, junc_temp.c.end,
                                   junc_temp.c.strand).distinct()
        distinct_junc = session.execute(distinct_junc).fetchall()
        distinct_junc = pd.DataFrame(distinct_junc)
        distinct_junc.columns = ["chrom", "start", "end", "strand"]
        # this may take a while and run out of memory probably need a generator
        distinct_junc.to_sql(table, engine, index=False, if_exists="append")
        del (distinct_junc)
        gc.collect()
    else:
        # this means there are already files in junctions and alljunctions tables so we need to figure
        # out what the junction ids are if they are in the table and insert new ones
        new_juncs = select(junc_temp.c.chrom, junc_temp.c.start, junc_temp.c.end,
                           junc_temp.c.strand).distinct().join(junc_table,
                                                               (junc_temp.c.chrom == junc_table.c.chrom) &
                                                               (junc_temp.c.start == junc_table.c.start) &
                                                               (junc_temp.c.end == junc_table.c.end) &
                                                               (junc_temp.c.strand == junc_table.c.strand)). \
            add_columns(junc_table.c.id).filter(junc_table.c.id == None)
        new_juncs = session.execute(new_juncs).fetchall()

        if len(new_juncs)>0:
            new_juncs = pd.DataFrame(new_juncs)
            new_juncs.columns = ["chrom", "start", "end", "strand", "id"]
            new_juncs = new_juncs.drop(columns=["id"])
            new_juncs.to_sql(table, engine, if_exists="append", index=False)

    query = select(junc_temp.c.samplename, junc_temp.c.uniq_map,
                   junc_temp.c.multi_map).join(junc_table, (junc_temp.c.chrom == junc_table.c.chrom) &
                                               (junc_temp.c.start == junc_table.c.start) &
                                               (junc_temp.c.end == junc_table.c.end) &
                                               (junc_temp.c.strand == junc_table.c.strand)). \
        add_columns(junc_table.c.id)
    junction_mapping = pd.DataFrame(session.execute(query).fetchall())
    junction_mapping.columns = ["samplename", "uniq_map", "multi_map", "junction"]
    junction_mapping.to_sql(mapping, engine, index=False, if_exists="append")
    del (junction_mapping)
    gc.collect()

    engine.execute("drop table if exists {}".format(temp))

