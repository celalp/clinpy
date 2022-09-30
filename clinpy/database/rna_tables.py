from sqlalchemy import Column, ForeignKey, Integer, String, Float
from clinpy.database.base_tables import ProjectBase



class FilteredJunctions(ProjectBase):  # these are the junctions that pass the intense filtering described elsewhere
    __tablename__ = "junctions"
    id = Column(Integer, primary_key=True, autoincrement=True)
    chrom = Column(String())
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(String(1))


class AllJunctions(ProjectBase):  # these are junctions that pass some basic QC that's it They are not processed
    # the way filtered junctions are processed this is just for record keeping
    __tablename__ = "all_junctions"
    id = Column(Integer, primary_key=True, autoincrement=True)
    chrom = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(String(1))


# many to many relationships
class SampleToJunction(ProjectBase):
    __tablename__ = "sample_to_junction"
    samplename = Column(ForeignKey("samples.sample_id"), primary_key=True, index=True)
    junction = Column(ForeignKey("junctions.id"), primary_key=True, index=True)
    uniq_map = Column(Integer)
    multi_map = Column(Integer)


class SampleToAllJunction(ProjectBase):
    __tablename__ = "sample_to_alljunction"
    samplename = Column(ForeignKey("samples.sample_id"), primary_key=True, index=True)
    junction = Column(ForeignKey("all_junctions.id"), primary_key=True, index=True)
    uniq_map = Column(Integer)
    multi_map = Column(Integer)


class GeneExpression(ProjectBase):
    __tablename__ = "gene_expression"
    samplename = Column(ForeignKey("samples.sample_id"), primary_key=True, index=True)
    gene = Column(String, primary_key=True, index=True)
    expected_count = Column(Float)
    tpm = Column(Float)
    fpkm = Column(Float)


class TranscriptExpression(ProjectBase):
    __tablename__ = "transcript_expression"
    samplename = Column(ForeignKey("samples.sample_id"), primary_key=True, index=True)
    transcript = Column(String, primary_key=True, index=True)
    expected_count = Column(Float)
    tpm = Column(Float)
    fpkm = Column(Float)
    isopct = Column(Float)

class RNAVariants(ProjectBase):
    __tablename__="rna_variants"
    variant_id=Column(Integer, primary_key=True, index=True)
    chrom=Column(String, index=True)
    pos=Column(Integer)
    id=Column(String)
    ref=Column(String)
    alt=Column(String)


class FilteredRNAVariants(ProjectBase):
    __tablename__="filtered_rna_variants"
    variant_id=Column(Integer, primary_key=True, index=True)
    chrom=Column(String, index=True)
    pos=Column(Integer)
    id=Column(String)
    ref=Column(String)
    alt=Column(String)