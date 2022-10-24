import pandas as pd
from sqlalchemy import Table, select, update, func
from clinpy.assays.assay_base import Assay


class Project(Assay):
    def __init__(self, db, genome, assays):
        super().__init__(db, genome)
        #TODO add assays


    def view_meta_fields(self):
        table = Table("samples", self.metadata, autoload=True, autoload_with=self.db)
        query=select(table.c.sample_meta)
        fields=self.session.execute(query).fetchone()
        fields=list(fields[0].keys())
        fields=["sample_id", "cohort"] + fields
        return fields

    def samples(self, cohort=None):
        """
        search for samples
        :param cohort: name of the cohort if none all samples an interable
        :return:
        """
        table = Table("samples", self.metadata, autoload=True, autoload_with=self.db)
        query = select(table.c.study_id, table.c.tm_id, table.c.cohort, table.c.user_annot)
        if cohort is not None:
            query = query.where(table.c.cohort.in_(cohort))

        results = self.session.execute(query).fetchall()
        results = pd.DataFrame(results)
        results.columns = list(self.session.execute(query).keys())

        mandatory=results.loc[:,["sample_id", "cohort"]]
        others=pd.DataFrame.from_records(results.sample_meta.to_list())
        final=pd.concat(mandatory, others, axis=1)

        return final

    def add_annotation(self, to, annot):
        """
        add some arbitrary annotation to a sample and save it to the database.
        If run multiple times overwrites the existing annotation
        :param to: study_id of the sample
        :param annot: dict of annotations
        :return: nothing adds the dict to the json field of either gene or transcript table
        """
        table = Table("samples", self.metadata, autoload=True, autoload_with=self.db)

        command = update(table).values(sample_meta={"user_annot": annot}).where(table.c.sample_id == to)
        self.session.execute(command)
        self.session.commit()


    def __str__(self):
        """
        :return: summary stats, number of sample in each cohort
        """
        table = Table("samples", self.metadata, autoload=True, autoload_with=self.db)
        query = select(table.c.cohort, func.count(table.c.cohort)).group_by(table.c.cohort)
        results = pd.DataFrame(self.session.execute(query).fetchall())
        cohorts = ",".join([str(x) for x in results[0].to_list()])
        num_samples = ",".join([str(x) for x in results[1].to_list()])
        return "{} cohorts with names {} and {} samples respectively".format(len(results[0].to_list()),
                                                                             cohorts, num_samples)


