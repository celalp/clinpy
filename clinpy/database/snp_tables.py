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
