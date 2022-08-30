from sqlalchemy import MetaData
from sqlalchemy import Column
from sqlalchemy import ForeignKey
from sqlalchemy import Integer
from sqlalchemy import String
from sqlalchemy import Float
from sqlalchemy import JSON
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship
import argparse as arg
from sqlalchemy import create_engine
import pandas as pd

# TODO all junction functions
# TODO dataframe generator to insert in chunks
# TODO subparser or some other form to specify to create database or add to database

if __name__=="__main__":
    parser = arg.ArgumentParser(description='add to a project database with genome annotations junctions'
                                            'and expression data')
    parser.add_argument('-s', '--samples', type=str, help='tab separated sample metadata', action="store")
    parser.add_argument('-j', '--junctions', type=str, help='filtered junctions file',
                        action="store")
    parser.add_argument('-a', '--all_junctions', type=str, help='all junctions file (optional)',
                        action="store")
    parser.add_argument('-g', '--gene_expression_file', type=str, help='merged gene expression file', action="store")
    parser.add_argument('-t', '--transcript_expression', type=str, help='merged isoform expression file',
                        action="store")
    parser.add_argument('-o', '--output', help='project database, will be a sqlite database to be used later')
    args = parser.parse_args()


    engine = create_engine("sqlite:///{}".format(args.output))
    project_meta = MetaData(bind=engine)
    ProjectBase = declarative_base(metadata=project_meta)

    class FilteredJunctions(ProjectBase):  # these are the junctions that pass the intense filtering described elsewhere
        __tablename__ = "junctions"
        id = Column(Integer, primary_key=True)
        chrom = Column(String())
        start = Column(Integer)
        end = Column(Integer)
        strand = Column(String(1))

    class AllJunctions(ProjectBase):  # these are junctions that pass some basic QC that's it They are not processed the way filtered junctions are processed this is just for record keeping
        __tablename__ = "all_junctions"
        id = Column(Integer, primary_key=True)
        chrom = Column(String)
        start = Column(Integer)
        end = Column(Integer)
        strand = Column(String(1))

    class Samples(ProjectBase):
        __tablename__ = "samples"
        id = Column(Integer, primary_key=True, index=True)
        study_id = Column(String(50))
        tm_id = Column(String(50))
        cohort = Column(String(3))
        user_annot=Column(JSON())
        junctions = relationship("sample_to_junction", back_populates="samples")
        gene_expression = relationship("gene_expression", back_populates="samples")
        tx_expression = relationship("transcript_expression", back_populates="samples")


    # many to many relationships
    class SampleToJunction(ProjectBase):
        __tablename__ = "sample_to_junction"
        samplename = Column(ForeignKey("samples.id"), primary_key=True, index=True)
        junction = Column(ForeignKey("junctions.id"), primary_key=True, index=True)
        uniq_map = Column(Integer)
        multi_map = Column(Integer)


    class SampleToAllJunction(ProjectBase):
        __tablename__ = "sample_to_alljunction"
        samplename = Column(ForeignKey("samples.id"), primary_key=True, index=True)
        junction = Column(ForeignKey("all_junctions.id"), primary_key=True, index=True)


    class GeneExpression(ProjectBase):
        __tablename__ = "gene_expression"
        samplename = Column(ForeignKey("samples.id"), primary_key=True, index=True)
        gene = Column(String, primary_key=True, index=True)
        expected_count = Column(Float)
        tpm = Column(Float)


    class TranscriptExpression(ProjectBase):
        __tablename__ = "transcript_expression"
        samplename = Column(ForeignKey("samples.id"), primary_key=True, index=True)
        transcript = Column(String, primary_key=True, index=True)
        expected_count = Column(Float)
        tpm = Column(Float)
        isopct = Column(Float)

    ProjectBase.metadata.create_all(engine)

    sample_meta = pd.read_csv(arg.samples, header=0, sep="\t")
    sample_meta.to_sql("samples", engine, if_exists="append", index=False)

    # takes aggregated data not ideal i know
    gene_expression = pd.read_csv(args.gene_expression, header=0, sep="\t")
    gene_expression.to_sql("gene_expression", engine, if_exists="append", index=False)

    isoform_expression = pd.read_csv(args.transcript_expression, header=0, sep="\t")
    isoform_expression.to_sql("transcript_expression", engine, if_exists="append", index=False)

    filtered_junctions = pd.read_csv(arg.junctions, header=0, sep="\t")
    uniq_junc = filtered_junctions.loc[:, ["chrom", "start", "end", "strand"]].drop_duplicates()
    uniq_junc.to_sql("junctions", engine, if_exists="append", index=False)

    # I need to return the autoincrement junc ids to match the junctions, need to a join, kinda redundant but
    # necessary if you want to add samples later on. Sqlite does not support returning so I need to query the db
    # this is still incomplete in the sense that if you have common junctions with the new samples it messes it up.
    # but it's a start. I think I can make this a lot more robust if needed

    inserted = pd.read_sql("select id,chrom,start,end,strand from junctions", engine)
    # now merge with the other junctions, being safe
    merged = pd.merge(filtered_junctions, inserted, on=["chrom", "start", "end", "strand"], how="inner")
    sample_to_junc = merged.loc[:, ["samplename", "id", "uniq_map", "multi_map"]]
    sample_to_junc = sample_to_junc.rename(columns={"id": "junction"})
    sample_to_junc.to_sql("sample_to_junction", engine, if_exists="append", index=False)

