import pandas as pd
import sqlalchemy as sql
from sqlalchemy.orm import Session
import pyranges
from .utils import calc_overlap
import numpy as np


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

    def samples(self, cohort=None, ids=None):
        """
        search for samples
        :param cohort: name of the cohort if none all samples an interable
        :param ids: names of the samples if none all samples an iterable
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


class Junction:
    def __init__(self, chrom, start, end, strand, uniq_map, multi_map):
        """
        initiate a junction instance bare bones
        :param chrom chrom
        :param chrom chrom
        :param start start
        :param end end
        :param strand strand
        :param uniq_map: uniq map reads
        :param multi_map multi_map reads
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.uniq_map = uniq_map
        self.multi_map = multi_map

    def genes(self, genome, df=False):
        """
        return overapping genes
        :param df return a dataframe or pyranges
        :param genome: genome class instance
        :return: genes that match the start and end of the junction and if there are other genes in the middle
        """
        table = sql.Table("genes", genome.metadata, autoload=True, autoload_with=genome)
        # find the start and end genes, there might be others in between
        query = sql.select(table).filter(sql.and_(table.c.chrom == self.chrom,
                                                  table.c.strand == self.strand)). \
            filter(sql.or_(sql.and_(table.c.start < self.start, table.c.end > self.start),
                           sql.and_(table.c.start < self.end, table.c.end > self.end)))

        results = pd.DataFrame(genome.session.execute(query).fetchall())
        results.columns = list(genome.session.execute(query).keys())

        if df:
            return results
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=results.iloc[:, 1],
                                            starts=results.iloc[:, 2],
                                            ends=results.iloc[:, 3],
                                            strands=results.iloc[:, 4])

            gr.id = results["id"]
            gr.name = results["name"]
            gr.description = results["description"]
            gr.biotype = results["biotype"]
            return gr

    def transcripts(self, genome, df=False, biotype=None):
        """
        get the transcript the junction overlaps
        :param genome: genome class instance
        :param df: return a dataframe or pyranges
        :param biotype: limit to a certain biotype like "protein coding"
        :return: a pyranges or dataframe with all the transcript that match the description.
        """

        gene_table = sql.Table("genes", genome.metadata, autoload=True, autoload_with=genome)
        tx_table = sql.Table("transcripts", genome.metadata,
                             autoload=True, autoload_with=genome)

        subq = sql.select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()
        query = sql.select(tx_table, subq.c.chrom, subq.c.strand). \
            join(subq, tx_table.c.gene == subq.c.id). \
            filter(sql.or_(sql.and_(tx_table.c.start < self.start, tx_table.c.end > self.start),
                           sql.and_(tx_table.c.start < self.end, tx_table.c.end > self.end)),
                   sql.and_(subq.c.chrom == self.chrom, subq.c.strand == self.strand))

        if biotype is not None:
            query = query.where(tx_table.c.biotype.in_(biotype))

        results = genome.session.execute(query).fetchall()
        results = pd.DataFrame(results)
        results.columns = genome.session.execute(query).keys()

        if df:
            return results
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=results.iloc[:, 6],
                                            starts=results.iloc[:, 2],
                                            ends=results.iloc[:, 3],
                                            strands=results.iloc[:, 7])
            gr.tx_id = results.id
            gr.gene_id = results.gene
            gr.biotype = results.biotype
            return gr

    def features(self, genome, txs):
        """
        return which intron/exon the transcript matches to
        :param genome: genome instance
        :param txs: names of the transcripts because an exon/intron might belong to multiple transccript
        :return: a dict with keys "Start" and "End" values are tuples first one is either "intron" or "exon" the second
        one is a df row with all the details
        """
        start_interval = pd.Interval(self.start, self.start, closed="both")
        end_interval = pd.Interval(self.end, self.end, closed="both")

        introns = genome.introns(names=txs, df=True, level="transcript")
        intron_intervals = pd.arrays.IntervalArray.from_arrays(introns.start.to_numpy(), introns.end.to_numpy(),
                                                               closed="both")
        exons = genome.exons(names=txs, df=True)
        exon_intervals = pd.arrays.IntervalArray.from_arrays(exons.start.to_numpy(), exons.end.to_numpy(),
                                                             closed="both")
        intervals = {"intron": intron_intervals, "exon": exon_intervals}
        dfs = {"intron": introns, "exon": exons}

        for interval in intervals.keys():
            if any(intervals[interval].overlaps(start_interval)):
                start_loc = np.where(intervals[interval].overlaps(start_interval) == True)[0]
                start_type = interval
            else:
                continue

            if any(intervals[interval].overlaps(end_interval)):
                end_loc = np.where(intervals[interval].overlaps(end_interval) == True)[0]
                end_type = interval
            else:
                continue

        return {"Start": (start_type, dfs[start_type].iloc[start_loc, :]),
                "End": (end_type, dfs[end_type].iloc[end_loc, :])}

    def samples(self, project, tolerance=None, overlap=None, reciprocal=False):
        """
        find all the samples in the project that have the same junction
        :param project: project instance
        :param tolerance: a tuple how much from the 5' and 3' end is allowed in nucleo
        :param overlap: fraction overlap with self
        :param reciprocal: if True reciprocal overlap must be >= overlap value
        :return:
        """
        junctions_table = sql.Table("junctions", project.metadata, autoload=True, autoload_with=project)
        sample_to_junction = sql.Table("sample_to_junction", project.metadata, autoload=True,
                                       autoload_with=project)
        subq = sql.select(junctions_table.c.id).filter(sql.and_(junctions_table.c.chrom == self.chrom,
                                                                junctions_table.c.strand == self.strand))
        if overlap is None:
            if tolerance is None:
                subq = subq.filter(sql.and_(junctions_table.c.start == self.start,
                                            junctions_table.c.end == self.end)).scalar_subquery()
            else:
                subq = subq.filter(sql.and_(junctions_table.c.start >= self.start - tolerance[0],
                                            junctions_table.c.end <= self.end + tolerance[1])).scalar_subquery()

            query = sql.select(sample_to_junction.c.samplename).where(sample_to_junction.c.junction.in_(subq))
            sample_ids = project.session.execute(query).fetchall()
            return sample_ids

        else:
            if tolerance is not None:
                raise ValueError("you can specify tolerance OR overlap")
            else:
                # need to get all junctions and perform overlaps and not get self
                query = subq.add_columns(junctions_table.c.start, junctions_table.c.end)
                junctions = project.session.execute(query).fetchall()
                overlapping = []
                for junction in junctions:
                    if reciprocal:
                        olap1 = calc_overlap((self.start, self.end), (junction[1], junction[2]))
                        olap2 = calc_overlap((junction[1], junction[2]), (self.start, self.end))
                        if olap1 >= overlap and olap2 >= overlap:
                            overlapping.append(junction[0])  # just the id
                        else:
                            continue
                    else:
                        olap = calc_overlap((self.start, self.end), (junction[1], junction[2]))
                        if olap >= overlap:
                            overlapping.append(junction[0])
                        else:
                            continue
                sample_query = sql.select(sample_to_junction.c.samplename).where(
                    sample_to_junction.c.junction.in_(overlapping))
                sample_ids = project.session.execute(sample_query).fetchall()
                return sample_ids

    def new_transcript(self, features, genome, sequence=True, type="nuc"):
        """
        return the new transcript sequence assuming all junctions are the same
        :param features: a features dict from the features method
        :param genome: genome class instance
        :param sequence: if true return sequence otherwise return df
        :param type: if sequence return either nuc or aa
        :return: returns either a dataframe or sequence of the new transcript
        """
        start_exons = genome.exons(features["Start"][1].transcript.to_list(), df=True)
        end_exons = genome.exons(features["End"][1].transcript.to_list(), df=True)

        if self.strand == "+":
            start_exons = start_exons[start_exons["start"] < self.start]
            start_exons["end"][start_exons.shape[0]] = self.start

            end_exons = end_exons[end_exons["end"] > self.end]
            end_exons["start"][end_exons.shape[0]] = self.end
        else:
            start_exons = start_exons[start_exons["start"] > self.start]
            start_exons["start"][start_exons.shape[0]] = self.start

            end_exons = end_exons[end_exons["end"] < self.end]
            end_exons["end"][end_exons.shape[0]] = self.end

        new_exons = pd.concat([start_exons, end_exons]).sort_values(by=["start"])

        if sequence:
            new_exons_gr = pyranges.pyranges.PyRanges(chromosomes=new_exons.iloc[:, 6],
                                                      starts=new_exons.iloc[:, 2],
                                                      ends=new_exons.iloc[:, 3],
                                                      strands=new_exons.iloc[:, 7])
            new_exons = genome.get_sequence(new_exons_gr, type=type, concat=True)

        return new_exons

    def __str__(self):
        """
        print method
        :return: string
        """
        "A junction in {} starting at {} and ending at {} on strand {}".format(self.chrom, self.start,
                                                                               self.end, self.strand)
