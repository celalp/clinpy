import pandas as pd
import sqlalchemy as sql
from sqlalchemy.orm import Session
from .junction import Junction

class Project:
    def __init__(self, db, genome):
        """
        initate the project object
        :param db: this is the project db it's a connection
        :param genome: this is the genome class instance from pytxdb
        """
        self.db = db
        self.session = Session(self.db)
        self.metadata = sql.MetaData(self.db)
        self.metadata.reflect(bind=self.db)
        self.genome = genome

    def samples(self, cohort=None):
        """
        search for samples
        :param cohort: name of the cohort if none all samples an interable
        :return:
        """
        table = sql.Table("samples", self.metadata, autoload=True, autoload_with=self.db)
        query = sql.select(table.c.study_id, table.c.tm_id, table.c.cohort, table.c.user_annot)
        if cohort is not None:
            query = query.where(table.c.cohort.in_(cohort))

        results = self.session.execute(query).fetchall()
        results = pd.DataFrame(results)
        results.columns = list(self.session.execute(query).keys())
        return results

    def add_annotation(self, to, annot):
        """
        add some arbitrary annotation to a sample and save it to the database.
        If run multiple times overwrites the existing annotation
        :param to: study_id of the sample
        :param annot: dict of annotations
        :return: nothing adds the dict to the json field of either gene or transcript table
        """
        table = sql.Table("samples", self.metadata, autoload=True, autoload_with=self.db)

        command = sql.insert(table).values(user_annot=annot).where(table.c.study_id == to)
        self.session.execute(command)
        self.session.commit()

    def expression(self, cohort=None, samples=None, gene=True, names=None, long=True, what=None):
        """
        return expression values
        :param cohort: samples from a cohort
        :param samples: specific samples
        :param gene: is this gene or isoform
        :param names: names of the genes/transcripts if none will return everything might take a while
        :param long: in long or wide format
        :param what: "TPM"/"counts" or "isopct" for transcripts optional if long=True
        :return:
        """
        if gene:
            table = sql.Table("gene_expression", self.metadata, autoload=True, autoload_with=self.db)
            query = sql.select(table.c.samplename, table.c.gene, table.c.expected_count, table.c.tpm)
        else:
            table = sql.Table("transcript_expression", self.metadata, autoload=True, autoload_with=self.db)
            query = sql.select(table.c.samplename, table.c.gene, table.c.expected_count, table.c.tpm,
                               table.c.isopct)

        if cohort is not None:
            sample_table = sql.Table("samples", self.metadata, autoload=True, autoload_with=self.db)
            samples = sql.select(sample_table.c.study_id).where(sample_table.c.cohort).in_(cohort). \
                scalar_subquery()
            query = query.where(table.samplename.in_(samples))

        if samples is not None:
            query = query.where(table.c.samplename.in_(samples))

        if names is not None:
            query = query.where(table.c.gene.in_(names))

        results = self.session.execute(query)
        results = pd.DataFrame(results)
        results.columns = list(self.session.execute(query).keys())

        if long:
            return results
        else:
            if what is not None:
                return results.pivot(index="gene", columns="samplename", values=what)
            else:
                raise ValueError(
                    "for wide format you need to specify either TPM, expectd count or isopct (transcripts only)")

    def junctions(self, cohort=None, uniq=False, samples=None, df=True):
        """
        returns a dataframe of junctions this is different than the junction class
        instances returned from Sample.junctions
        :param cohort: name of the cohort to get
        :param uniq: return unique junctions not grouped by sample
        :param samples: return junctions from those samples
        :param df: return a dataframe otherwise Junction instance
        :return: a dataframe of junctions and reads mapping to them, if uniq is
        selected mappings will not be there since different samples will have different
        numbers of reads mapping to each junction.
        """

        table = sql.Table("junctions", self.metadata, autoload=True, autoload_with=self.db)
        sample_to_junction = sql.Table("sample_to_junction", self.metadata,
                                       autoload=True, autoload_with=self.db)
        query = sql.select(table.c.chrom, table.c.start, table.c.end, table.c.strand)

        if not uniq:
            query = query.add_columns(sample_to_junction.c.samplename, sample_to_junction.c.uniq_map,
                                      sample_to_junction.c.multi_map)
            query = query.join(sample_to_junction, table.c.id == sample_to_junction.c.junction)

        if cohort is not None:
            sample_table = sql.Table("samples", self.metadata, autoload=True, autoload_with=self.db)
            samples = sql.select(sample_table.c.study_id).filter(sample_table.c.cohort.in_(cohort))

            junctions = sql.select(sample_to_junction.c.junction).distinct(). \
                filter(sample_to_junction.c.samplename.in_(samples))

            query = query.filter(table.c.id.in_(junctions))

        if samples is not None:
            junctions = sql.select(sample_to_junction.c.junction). \
                filter(sample_to_junction.c.samplename.in_(samples)).distinct()

            query = query.filter(table.c.id.in_(junctions))

        results = self.session.execute(query).fetchall()
        results = pd.DataFrame(results)
        results.columns = list(self.session.execute(query).keys())

        if df:
            return results
        elif not df and not uniq:
            juncs = []
            for chrom, start, end, strand, uniq_map, multi_map in zip(results.chrom, results.start, results.end,
                                                                      results.strand, results.uniq, results.multi):
                juncs.append(Junction(chrom, start, end, strand, uniq_map, multi_map))
                return juncs
        else:
            raise NotImplementedError("returning unique junctions as a junction class is not implemented")

    def __str__(self):
        """
        :return: summary stats, number of sample in each cohort
        """
        table = sql.Table("samples", self.metadata, autoload=True, autoload_with=self.db)
        query = sql.select(table.c.cohort, sql.func.count(table.c.cohort)).group_by(table.c.cohort)
        results = pd.DataFrame(self.session.execute(query).fetchall())
        cohorts = ",".join([str(x) for x in results[0].to_list()])
        num_samples = ",".join([str(x) for x in results[1].to_list()])
        return "{} cohorts with names {} and {} samples respectively".format(len(results[0].to_list()),
                                                                             cohorts, num_samples)

    def search_for_junctions(self, gr, samples=None, unique=False):
        """
        search a given region for junctions
        :param gr: a pyranges
        :param samples: sample id, if none search all samples otherwise limit to the list of samples provided an interable
        :param  unique: whether to return unique junctions, sample infomation will not be present, otherwise return
        other data as well
        :return: a dataframe of junctions that are in a given region
        """
        junctions = sql.Table("junctions", self.metadata, autoload=True, autoload_with=self.db)
        sample_to_junction = sql.Table("sample_to_junction", self.metadata,
                                       autoload=True, autoload_with=self.db)

        query = sql.select(junctions.c.chrom, junctions.c.start, junctions.c.end, junctions.c.strand).filter(
            sql.and_(junctions.c.chrom==str(gr.Chromosomes), junctions.c.strand==str(gr.Strand))).filter(
            sql.and_(junctions.c.end >= int(gr.Start), junctions.c.start<=int(gr.End)))

        if samples is not None:
            sample_junctions=sql.select(sample_to_junction.c.junction).filter(sample_to_junction.c.samplename.in_(samples))
            query=query.filter(junctions.c.id.in_(sample_junctions))

        if not unique:
            query = query.add_columns(sample_to_junction.c.samplename, sample_to_junction.c.uniq_map,
                                      sample_to_junction.c.multi_map)
            query = query.join(sample_to_junction, junctions.c.id == sample_to_junction.c.junction)

        results=self.session.execute(query).fetchall()
        columns=list(self.session.execute(query).keys())
        results=pd.DataFrame(results)
        results.columns=columns

        return results
