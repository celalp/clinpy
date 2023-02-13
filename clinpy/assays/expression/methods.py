from sqlalchemy import Table, select
import pandas as pd


class Expression:
    """
    This is the expression class it is used to query expression tables and return either wide
    or long format expression table
    """
    def __init__(self, gene_expression_table, transcript_expression_table, project, name):
        self.project=project
        gene_table = Table(gene_expression_table, self.project.metadata, autoload=True, autoload_with=self.project.db)
        transcript_table = Table(transcript_expression_table, self.project.metadata, autoload=True, autoload_with=self.project.db)
        self.gene_table=gene_table
        self.transcript_table=transcript_table
        self.assay_type="expression"
        self.assay_name=name


    def get_expression(self, samples=None, gene=True, names=None, long=True, what=None):
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
                               table.c.fpkm, table.c.tpm)
        else:
            table=self.transcript_table
            query = select(table.c.samplename, table.c.transcript, table.c.expected_count,
                           table.c.fpkm, table.c.tpm, table.c.isopct)

        if samples is not None:
            query = query.where(table.c.samplename.in_(samples))

        if names is not None:
            if gene:
                query = query.where(table.c.gene.in_(names))
            else:
                query = query.where(table.c.transcript.in_(names))

        results = self.project.session.execute(query)
        results = pd.DataFrame(results)
        results.columns = list(self.project.session.execute(query).keys())

        if long:
            return results
        else:
            if what is not None:
                if gene:
                    return results.pivot(index="gene", columns="samplename", values=what)
                else:
                    return results.pivot(index="transcript", columns="samplename", values=what)
            else:
                raise ValueError(
                    "for wide format you need to specify either TPM, FPKM, expected count or isopct (transcripts only)")


