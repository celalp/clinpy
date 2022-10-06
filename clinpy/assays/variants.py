import pandas as pd
from sqlalchemy import Table, select, and_, func, distinct
from clinpy.assays.assay_base import Assay
from sqlalchemy_filters import apply_filters
from clinpy.utils.snp_functions import add_to_variant_tables

class Variants(Assay):
    def __init__(self, db, rna=False, filtered=False):
        """
         this is a simple class for interacting with the databse in terms of searching for variants, and
         filtering based on arbitrary criteria, this is not the variant class this is more about getting variants then
         looking at specifiic attributes of a variant
        :param db: Project class
        :param genome: Genome class from pytxdb
        :param rna: are you referring to RNA-Seq variants
        :param filtered: are referring to filtered variants
        """
        super.__init__(db)
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
        self.rna=rna

    def list_impacts(self):
        """
        list the impact columns
        :return: list to table columns to be used for filtering
        """
        return list(self.impacts_table.columns.keys())

    def list_variant_quals(self):
        """
        list the columns like qual, filter and the genotype columns from the format column of the vcf
        :return: list of column names
        """
        return list(self.mapping_table.columns.keys())

    #TODO need a more memory efficient way to get this done
    def filter(self, impacts=None, formats=None):
        """
        filter unfiltered variants, the filters need to be a dict where the key is the field name and value is what you are searching for
        for categorical need to provide a list for numerircal if you want exact values you can provide a list again if you want comparisons
        you would need to provide
        :param impacts: impact filters a SQLAlchem filters json
        :param formats: same for fomat fields
        :param save: save results to "filtered variants table(s).
        :return: a dataframe of variants with all the applicable columns
        """
        #TODO put it back in the database
        query=select(self.variants_table).join(self.impacts_table).join(self.impacts_table)

        if impacts is None and formats is None:
            raise ValueError("Both impacts and formats are none, there are no filters specified")

        if impacts is not None:
            query=apply_filters(query, impacts)

        if formats is not None:
            query=apply_filters(query, formats)

        # if too many variants are there this might fail but not sure how many is too many
        variants=self.session.execute(query).fetchall()
        variants=pd.DataFrame(variants)
        variants.columns=self.session.execute(query).keys()

        return variants

    def search_region(self, gr, samples=None):
        """
        search a given genomic region for variants
        :param gr: a pyranges object strand is ignored if there
        :param samples search specific samples for variants
        :return: a dataframe of variants and impact, a variant may have more than one impact
        """
        format_cols = self.list_variant_quals()[1:]

        query = select(self.variants_table).filter(and_(self.variants_table.c.pos <= int(gr.Start),
                                                        self.variants_table.c.pos >= int(gr.End)),
                                                   self.variants_table.c.chrom == gr.Chromosomes). \
            join(self.mapping_table, self.mapping_table.c.variant_id == self.variants_table.c.variant_id). \
            add_columns(*[self.mapping_table.c[col] for col in format_cols])

        if samples is not None:
            query.filter(self.mapping_table.c.variant_id.in_(samples))

        results = self.session.execute(query).fetchall()

        if len(results) > 0:
            results = pd.DataFrame(results)
            colnames = self.session.execute(query).keys()
            results.columns = colnames
        else:
            results = None

        return results

    def __str__(self):
        num_samples = self.session.query(func.count(distinct(self.mapping_table.c.samplename)))
        num_variants = self.session.query(func.count(self.variants_table.c.variant_id))
        num_impacts = self.session.query(func.count(self.impacts_table.c.variant_id))

        desc = "{} variants from {} samples with {} different impact annotations".format(num_samples,
                                                                                         num_variants,
                                                                                         num_impacts)


class Variant(Variants):
    """
    same class will be used for both rna and dna variants
    this one is for learning about a variant and searching other samples for the same variant. You can have an iterable of these
    if you want. Like all other assays this is lazy to prevent unnecessary queries
    """

    def __init__(self, project, rna, filtered, variant_id, chrom, pos, ref, alt):
        self.id=variant_id,
        self.chrom=chrom
        self.pos=pos
        self.ref=ref
        self.alt=alt
        super().__init__(project, rna, filtered)

    #TODO
    @property
    def counts(self, samples=None, cohort=None):
        ac_query=None
        af_query=None
        num_hets=None
        num_homs=None
        pass

    @property
    def samples(self, cohort=None, genotype="both"):
        query=select(self.mapping_table).filter(self.mapping_table.c.variant_id==self.id)
        results = self.session.execute(query).fetchall()
        results = pd.DataFrame(results)
        results.columns = self.session.execute(query).keys()

        if genotype=="both":
            return results
        elif genotype=="hom":
            results=results[results["GT"]=="(1, 1)"]
            return results
        elif genotype=="het":
            results = results[results["GT"] == "(0, 1)"]
            return results


    @property
    def impact(self):
        query=select(self.impacts_table).filter(self.impacts_table.c.variant_id==self.id)
        results=self.session.execute(query).fetchall()
        results=pd.DataFrame(results)
        results.columns=self.session.execute(query).keys()

        return results

    def __str__(self):
        return "A variant in {}:{} with from {} to {}".format(self.chrom, self.pos,
                                                              self.ref, self.alt)
