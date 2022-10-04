from clinpy.assays.assay_base import Assay
from sqlalchemy import Table, select
import pandas as pd


class Expression(Assay):
    """
    This is the expression class it is used to query expression tables and return either wide
    or long format expression table
    """
    def __init__(self, db, genome):
        super.__init__(db, genome)
        gene_table=Table("gene_expression", self.metadata, autoload=True, autoload_with=self.db)
        transcript_table=Table("transcript_expression", self.metadata, autoload=True, autoload_with=self.db)


    def get_expression(self, cohort=None, samples=None, gene=True, names=None, long=True, what=None):
        """
       return expression values
       :param cohort: samples from a cohort
       :param samples: specific samples
       :param gene: is this gene or isoform
       :param names: names of the genes/transcripts if none will return everything might take a while
       :param long: in long or wide format
       :param what: "TPM"/"FPKM/"counts" or "isopct" for transcripts optional if long=True
       :return: a dataframe of expression
       """
        if gene:
            table=self.gene_table
            query = select(table.c.samplename, table.c.gene, table.c.expected_count,
                               table.c.FPKM, table.c.tpm)
        else:
            table=self.transcript_table
            query = select(table.c.samplename, table.c.gene, table.c.expected_count,
                           table.c.fpkm, table.c.tpm, table.c.isopct)

        if cohort is not None:

            samples = select(self.sample_table.c.study_id).where(self.sample_table.c.cohort).in_(cohort). \
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
                    "for wide format you need to specify either TPM, FPKM, expected count or isopct (transcripts only)")


