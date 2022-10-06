from sqlalchemy import Column, ForeignKey, Integer, String, Float
from .base_tables import ProjectBase


class Variants(ProjectBase):
    __tablename__ = "variants"
    variant_id = Column(Integer, primary_key=True, index=True)
    chrom = Column(String, index=True)
    pos = Column(Integer, index=True)
    id = Column(String)
    ref = Column(String)
    alt = Column(String)

class FilteredVariants(ProjectBase):
    __tablename__ = "filtered_variants"
    variant_id = Column(Integer, primary_key=True, index=True)
    chrom = Column(String, index=True)
    pos = Column(Integer, index=True)
    id = Column(String)
    ref = Column(String)
    alt = Column(String)
