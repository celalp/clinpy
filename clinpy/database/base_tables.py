from sqlalchemy.orm import declarative_base
from sqlalchemy import String, JSON, Column

ProjectBase = declarative_base()


class Samples(ProjectBase):
    __tablename__ = "samples"
    sample_id = Column(String(), primary_key=True, index=True)
    cohort = Column(String(), index=True)
    sample_meta = Column(JSON())

