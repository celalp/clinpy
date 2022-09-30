import pandas as pd
from sqlalchemy import Table, MetaData
from sqlalchemy.orm import Session

class Variants:
    def __init__(self, db, genome, rna=False, filtered=False):
        """
         this is a simple class for interacting with the databse in terms of searching for variants, and
         filtering based on arbitrary criteria, this is not the variant class this is more about getting variants then
         looking at specifiic attributes of a variant
        :param db: Project class
        :param genome: Genome class from pytxdb
        :param rna: are you referring to RNA-Seq variants
        :param filtered: are referring to filtered variants
        """
        self.project=db
        self.session = Session(self.db)
        self.metadata = MetaData(self.db)
        self.metadata.reflect(bind=self.db)
        self.genome = genome
        table = "variants"
        mapping = "sample_variants"
        impacts = "variant_impacts"

        if rna:
            table = "rna_" + table
            mapping = "rna_" + mapping
            impacts = "rna_" + impacts

        if filtered:
            table = "filtered_" + table
            mapping = "filtered_" + mapping
            impacts = "filtered_" + impacts

        self.variants_table = Table(table, self.metadata, autoload=True, autoload_with=self.db)
        self.mapping_table = Table(mapping, self.metadata, autoload=True, autoload_with=self.db)
        self.impacts_table = Table(impacts, self.metadata, autoload=True, autoload_with=self.db)

    def list_impacts(self):
        """
        list the impact columns
        :return: list to table columns to be used for filtering
        """
        return(list(self.impacts_table.columns.keys()))

    def list_variant_quals(self):
        """
        list the columns like qual, filter and the genotype columns from the format column of the vcf
        :return: list of column names
        """
        return(list(self.mapping_table.columns.keys()))

    def filter(self, impacts, formats, save=False, action="append", distinct=False):
        """
        filter unfiltered variants, the filters need to be a dict where the key is the field name and value is what you are searching for
        for categorical need to provide a list for numerircal if you want exact values you can provide a list again if you want comparisons
        you would need to provide
        :param impacts: filters from impacts fields like gnome_af etc.
        :param formats: format filters
        :param to_table: save results to "filtered variants table(s).
        :param action: either append to the filtered variants tables or overwrite them, to_table must be true otherwise ignored
        :param distinct return the variants that match the description, if the impact is on multiple transcripts all will be returned
        but no sample information, if false sample information (samplename, GT, qual DP etc. will be returned as well.)
        :return: a dataframe of variants with all the applicable columns
        """
        pass

    def __str__(self):
        pass



class Variant:
    """
    same class will be used for both rna and dna variants
    this one is for learning about a variant and searching other samples for the same variant. You can have an iterable of these
    if you want, this class is mostly eager that is everything is stored in memory whereas variants class only stores connections and
    tables and nothing is retriveved until requested.
    """
    def __init__(self, project, variant_id, chrom, pos, ref, alt):
        pass

    @property
    def count(self, samples=None, cohort=None):
        pass

    @property
    def freq(self, samples=None, cohort=None):
        pass

    @property
    def samples(self, cohort=None, genotype="both"):
        pass

    @property
    def impact(self, highest=True):
        pass

    def __str__(self):
        pass