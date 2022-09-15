from sqlalchemy import MetaData, Table, select, insert
from clinpy.database.base_tables import *
from clinpy.database.rna_tables import *
import argparse as arg
from sqlalchemy import create_engine
import pandas as pd
import os
import gc
from datetime import datetime

def modify_strand(df):
    if df["strand"]==0: #undefined
        return "."
    elif df["strand"]==1: #+
        return "+"
    elif df["strand"]==2: #-
        return "-"
    else:
        raise ValueError("unknown value in strand column")

if __name__=="__main__":
    parser = arg.ArgumentParser(description='add to a project database with genome annotations junctions'
                                            'and expression data')
    parser.add_argument('-s', '--samples', type=str, help='tab separated sample metadata must include columns sample_id and cohort',
                        action="store")
    parser.add_argument('-f', '--files', type=str, help='files file a tsv that contains the following columns and no header:'
                                                        'samplename, gene_expression_file_path, isoform_expression_file_path, junction_file_path, '
                                                        'filtered_junction_file_path cells either or both of the junction table can be empty')
    parser.add_argument('-o', '--output', help='project database, will be a sqlite database to be used later')
    parser.add_argument('-c', '--create', help="create database otherwise add to it", type=bool, action="store_true")
    parser.add_argument('--min_junc_reads', help="min uniq reads for a junction to be considered", default=10, action="store", type=int)
    args = parser.parse_args()

    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Creating database " + args.output)

    if not args.create and os.path.isfile(args.output):
        raise FileExistsError("database already exists")

    engine = create_engine("sqlite:///{}".format(args.output))
    project_meta = MetaData(bind=engine)
    session=Session(engine)

    sample_meta = pd.read_csv(arg.samples, header=0, sep="\t")
    sample_meta["user_annot"]=""
    cols=sample_meta.columns

    if "sample_id" not in cols:
        raise ValueError("sample_id column must be present in sample metadata")

    if "cohort" not in cols:
        raise ValueError("cohort column must be preset in sample metadata")
    if args.create:
        ProjectBase.metadata.create_all(engine)

    project_meta.reflect()

    # sample id and cohort are mandatory columns the rest will be converted to json and will be inserted that way
    man_columns=["sample_id", "cohort"]
    mandatory = sample_meta.loc[:, man_columns]
    mandatory = mandatory.to_dict(orient="records")
    others = sample_meta.drop(columns=man_columns)
    others = others.to_dict(orient="records")

    sample_meta = []
    for meta, man in zip(others, mandatory):
        man["sample_meta"] = meta
        sample_meta.append(man)


    sample_table=Table("samples", project_meta, autoload=True, autoload_with=engine)

    sample_meta = []
    for meta, man in zip(others, mandatory):
        man["sample_meta"] = meta
        sample_meta.append(man)

    # this will check the unique constraint
    sample_insert = insert(sample_table).values(sample_meta)
    engine.execute(sample_insert)

    files=pd.read_csv(args.files, sep="\t", header=None)
    samplenames=files[0].to_list()
    gene_expression=files[1].to_list()
    isoform_expression=files[2].to_list()
    all_juncs=files[3].to_list()
    filtered_juncs=files[4].to_list()

#insert sample data will also add additional user_annot field

# I'm not checking if your gene/tx names match the genome, there is no genome involved here
# that's up to the user
    for sample, gene, tx, all_junc, junc in zip(samplenames, gene_expression, isoform_expression,
                                                all_juncs, filtered_juncs):
        print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Adding gene expression for " + sample)
        g=pd.read_csv(gene, header=0, sep="\t")
        g["samplename"]=sample
        g=g.drop(columns=["length", "effective_length", "transcript_id(s)"])
        g.columns=["samplename", "gene", "expected_count", "tpm", "fpkm"]
        g.to_sql("gene_expression", engine, if_exists="append", index=False)

        print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Adding isoform expression for " + sample)
        t=pd.read_csv(tx, header=0, sep="\t")
        t["samplename"]=sample
        t=t.drop(columns=["length", "effective_length", "gene_id"])
        t.columns=["samplename", "transcript", "expected_count", "tpm", "fpkm", "isopct"]
        t.to_sql("transcript_expression", engine, if_exists="append", index=False)

        print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Adding junctions for " + sample)
        ja=pd.read_csv(all_junc, header=None, sep="\t", low_memory=False)
        ja["samplename"]=sample
        ja.columns=["chrom", "start", "end", "strand", "motif", "annotated",
                    "uniq_map", "multi_map", "max_ohang", "samplename"]
        # star annotated is useless because it is actually not the annotated but annotated+detected in
        # first pass
        if not pd.isna(junc):
            ja=ja.drop(columns=["max_ohang", "motif", "annotated"])
            ja["strand"]=ja.apply(modify_strand, axis=1)
            ja=ja[(ja.uniq_map>=args.min_junc_reads) & (junc.strand!=".")]
            ja.to_sql("temp_all_junc", engine, if_exists="append", index=False)

        if not pd.isna(junc):
        #maybe check the column names?
            print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Adding filtered junctions for " + sample)
            jf=pd.read_csv(junc, header=None, sep="\t")
            jf["samplename"]=sample
            jf.columns = ["samplename", "chrom", "start", "end", "strand", "motif", "annotated",
                      "uniq_map", "multi_map", "max_ohang"]
            jf=jf.drop("max_ohang", "motif", "annotated")
            jf["strand"]=jf.apply(modify_strand, axis=1)
            jf.to_sql("temp_filt_junc", engine, if_exists="append", index=False)
        else:
            print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "No filtered junctions for " + sample)

    #init tables
    project_meta.reflect()
    junc_temp = Table("temp_all_junc", project_meta, autoload=True, autoload_with=engine)
    all_junc_table = Table("all_junctions", project_meta, autoload=True, autoload_with=engine)
    filt_junc_table = Table("junctions", project_meta, autoload=True, autoload_with=engine)
    samp_to_alljunc = Table("sample_to_alljunction", project_meta, autoload=True, autoload_with=engine)
    samp_to_junc = Table("sample_to_junction", project_meta, autoload=True, autoload_with=engine)

    if "temp_filt_junc" in project_meta.tables.keys():
        filt_temp = Table("temp_filt_junc", project_meta, autoload=True, autoload_with=engine)

    if args.create:
        print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Creating junctions table")
        # this means the all_junctions and junction tables are empty
        distinct_all_junc=select(junc_temp.c.chrom, junc_temp.c.start, junc_temp.c.end,
                                  junc_temp.c.strand).distinct()
        distinct_all_junc=session.execute(distinct_all_junc).fetchall()
        distinct_all_junc=pd.DataFrame(distinct_all_junc)
        distinct_all_junc.columns=["chrom", "start", "end", "strand"]
        # this may take a while and run out of memory probably need a generator
        distinct_all_junc.to_sql("all_junctions", engine, index=False, if_exists="append")
        del(distinct_all_junc)
        gc.collect()

        if "temp_filt_junc" in project_meta.tables.keys():
            print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Creating filtered junctions table")
            distinct_filt_junc = select(filt_temp.c.chrom, filt_temp.c.start, filt_temp.c.end,
                                       filt_temp.c.strand).distinct()
            distinct_filt_junc = session.execute(distinct_filt_junc).fetchall()
            if len(distinct_filt_junc) > 0:
                distinct_filt_junc = pd.DataFrame(distinct_filt_junc)
                distinct_filt_junc.columns = ["chrom", "start", "end", "strand"]
                distinct_filt_junc.to_sql("junctions", engine, index=False, if_exists="append")
                del (distinct_filt_junc)
                gc.collect()

    else:
        # this means there are already files in junctions and alljunctions tables so we need to figure
        # out what the junction ids are if they are in the table and insert new ones
        print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Creating junctions table")
        new_juncs=select(junc_temp.c.chrom, junc_temp.c.start, junc_temp.c.end,
                                  junc_temp.c.strand).distinct().join(all_junc_table, (junc_temp.c.chrom == all_junc_table.c.chrom) &
                                                                      (junc_temp.c.start == all_junc_table.c.start) &
                                                                      (junc_temp.c.end == all_junc_table.c.end) &
                                                                      (junc_temp.c.strand == all_junc_table.c.strand)).\
            add_columns(all_junc_table.c.id).filter(all_junc_table.c.id==None)
        new_juncs=pd.DataFrame(session.execute(new_juncs).fetchall())
        new_juncs.columns=["chrom", "start", "end", "strand", "id"]
        new_juncs=new_juncs.drop(columns=["id"])
        new_juncs.to_sql("all_junctions", engine, if_exists="append", index=False)

        # repeat for filtered junctions
        if "temp_filt_junc" in project_meta.tables.keys():
            print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Creating filtered junctions table")
            distinct_filt_junc = select(filt_temp.c.chrom, filt_temp.c.start, filt_temp.c.end,
                                        filt_temp.c.strand).distinct()
            distinct_filt_junc = session.execute(distinct_filt_junc).fetchall()
            if len(distinct_filt_junc) > 0:
                new_juncs = select(filt_temp.c.chrom, filt_temp.c.start, filt_temp.c.end,
                                   filt_temp.c.strand).distinct().join(filt_junc_table,
                                                                       (filt_temp.c.chrom == filt_junc_table.c.chrom) &
                                                                       (filt_temp.c.start == filt_junc_table.c.start) &
                                                                       (filt_temp.c.end == filt_junc_table.c.end) &
                                                                       (filt_temp.c.strand == filt_junc_table.c.strand)). \
                    add_columns(filt_junc_table.c.id).filter(filt_junc_table.c.id == None)
                new_juncs = pd.DataFrame(session.execute(new_juncs).fetchall())
                new_juncs.columns = ["chrom", "start", "end", "strand", "id"]
                new_juncs = new_juncs.drop(columns=["id"])
                new_juncs.to_sql("junctions", engine, if_exists="append", index=False)


    # this is going to be very slow as numbers increase
    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Creating junction mapping for all junctions")
    query = select(junc_temp.c.samplename, junc_temp.c.uniq_map,
                   junc_temp.c.multi_map).join(all_junc_table, (junc_temp.c.chrom == all_junc_table.c.chrom) &
                                               (junc_temp.c.start == all_junc_table.c.start) &
                                               (junc_temp.c.end == all_junc_table.c.end) &
                                               (junc_temp.c.strand == all_junc_table.c.strand)). \
        add_columns(all_junc_table.c.id)
    junction_mapping = pd.DataFrame(session.execute(query).fetchall())
    junction_mapping.columns=["samplename", "uniq_map", "multi_map", "junction"]
    junction_mapping.to_sql("sample_to_alljunction", engine, index=False, if_exists="append")
    del (junction_mapping)
    gc.collect()

    engine.execute("drop table if exists temp_all_junc")

    if "temp_filt_junc" in project_meta.tables.keys():
        print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Creating junction mapping for filtered junctions")
        query = select(junc_temp.c.samplename, junc_temp.c.uniq_map,
                       junc_temp.c.multi_map).join(all_junc_table, junc_temp.c.chrom == all_junc_table.c.chrom,
                                                   junc_temp.c.start == all_junc_table.c.start,
                                                   junc_temp.c.end == all_junc_table.c.end,
                                                   junc_temp.c.strand == all_junc_table.c.strand). \
            add_columns(all_junc_table.c.id)
        samp_to_junc = Table("sample_to_junction", project_meta, autoload=True, autoload_with=engine)
        junction_mapping = pd.DataFrame(session.execute(query).fetchall())
        junction_mapping.columns = ["samplename", "uniq_map", "multi_map", "junction"]
        junction_mapping.to_sql("sample_to_junction", engine, index=False, if_exists="append")
        del (junction_mapping)
        gc.collect()

        engine.execute("drop table if exists temp_filt_junc")

    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Done!")

