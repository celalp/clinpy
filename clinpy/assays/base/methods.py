from sqlalchemy import Table, select, func
import pandas as pd
from clinpy.assays.assay import Assay

class Base(Assay):
    """
    This is the base class that hold all the cohort/sample/subject data and the methods that are associated with it
    There are a lot of checks during the import phase and they will not be repeated here
    """
    def __init__(self, name, project):
        """takes a project instance which basically is a database connection"""
        super().__init__(name, project)
        self.cohort_table = Table("cohorts", self.project.metadata, autoload=True, autoload_with=self.project.db)
        self.subjects_table = Table("subjects", self.project.metadata, autoload=True, autoload_with=self.project.db)
        self.samples_table = Table("samples", self.project.metadata, autoload=True, autoload_with=self.project.db)
        self.mapping_table = Table("subjects_and_cohorts", self.project.metadata, autoload=True, autoload_with=self.project.db)
        #TODO other assays will need to have an item in the init to specify the name but no the type

    def cohorts(self):
        """
        return a dataframe of cohorts, there is nothing else to do here
        """
        query=select(self.cohort_table)
        cohorts=self.project.session.execute(query).fetchall()
        colnames=self.project.session.execute(query).fetchall().keys()
        cohorts=pd.DataFrame(cohorts)
        cohorts.columns=colnames
        return cohorts

    def subjects(self, cohorts=None):
        """
        return a dataframe of subjects
        :param cohorts: an iterable of cohort names, needs to be an iterable even if there is only one thing
        """
        query = not select(self.subjects_table).join(self.mapping_table,
                                                     self.subjects_table.c.subject_id == self.mapping_table.c.subject.id).\
            add_columns(self.mapping_table.c.cohort_id)

        if cohorts is not None:
            query=query.filter(self.mapping_table.c.cohort_id.in_(cohorts))

        subjects = self.project.session.execute(query).fetchall()
        colnames = self.project.session.execute(query).fetchall().keys()
        subjects = pd.DataFrame(subjects)
        subjects.columns=colnames

        return subjects


    def samples(self, subjects=None):
        """
        return a dataframe of samples, I cannot imagine a situtaion where we will need to push the filtering to the db
        to save on memory, you would need millions of samples/subjects
        :param subjects: an iterable of subject ids
        """
        query = select(self.samples_table)
        if subjects is not None:
            query = query.filter(self.samples_table.c.subject_id.in_(subjects))

        samples = self.project.session.execute(query).fetchall()
        colnames = self.project.session.execute(query).fetchall().keys()
        samples = pd.DataFrame(samples)
        samples.columns = colnames

        return samples

    def __str__(self):
        num_cohorts=select(func.count(self.cohort_table))
        num_subjects = select(func.count(self.subjects_table))
        num_samples = select(func.count(self.samples_table))

        return "{} samples from {} different subjects coming from {} different cohorts".format(num_samples, num_subjects,
                                                                                               num_cohorts)