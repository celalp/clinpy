from sqlalchemy import Column, ForeignKey, Integer, String, Float
from .base_tables import ProjectBase


class Variants(ProjectBase):
    __tablename__ = "variants"
    id = Column(Integer, primary_key=True, autoincrement=True, Index=True)
    chrom = Column(String(), Index=True)
    pos = Column(Integer, Index=True)
    id = Column(String)
    ref = Column(String)
    alt = Column(String)

class VariantImpacts(ProjectBase):
    __tablename__="impacts"
    id=Column(ForeignKey("variants.id"), Index=True)
    # TODO find a way to get columns and types

class SampleVariants(ProjectBase):
    __tablename__="sample_to_variant"
    variant_id=Column(ForeignKey("variants.id"), index=True, primary_key=True)
    sample_id=Column(ForeignKey("samples.sample_id"), index=True, primary_key=True)
    qual=Column()
    #TODO figure out a way to get the column types