from sqlalchemy import MetaData, Table
from sqlalchemy.orm import Session


class Assay:
    """
    This is a generic assay class it will hold just an init method for database connections
    """

    def __init__(self, db, genome=None):
        self.db = db
        self.session = Session(self.project)
        self.metadata = MetaData(self.project)
        self.metadata.reflect(bind=self.project)

        self.sample_table = Table("samples", self.metadata, autoload=True, autoload_with=self.db)
        if genome is not None:
            self.genome = genome
