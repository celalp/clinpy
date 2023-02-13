from sqlalchemy import MetaData, Table
from sqlalchemy.orm import Session
from sqlalchemy.orm import declarative_base
from clinpy.utils.utils import dict_to_engine

ProjectBase = declarative_base()

#TODO need to make sure that the project loads the assays that are neeeded
#TODO need to find a way to check if a given assay is there and it's tables are thered
class Project:
    """
    This is a generic class it will hold just an init method for database connections and the assays that are in the project
    """

    def __init__(self, db_params, genome=None):
        """
        Simple init method to store the connection and get all the tables in the database
        :param db: this is a sqlalchemy engine
        :param genome: this is a pytxdb instance, it can reside in the same database or someplace else, if a
        server is used for project and a sqlite for genome should generate a warning
        """
        db=dict_to_engine(db_params["dbtype"], db_params["name"], db_params["user"],
                          db_params["pwd"], db_params["host"], db_params["port"])
        self.db = db
        self.session = Session(self.db)
        self.metadata = MetaData(self.db)
        self.metadata.reflect(bind=self.db)

        if genome is not None:
            self.genome=genome

    def check_assays(self, config):
        pass

