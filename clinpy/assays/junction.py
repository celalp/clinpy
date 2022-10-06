import pandas as pd
from sqlalchemy import Table, select, and_, or_
import pyranges
from clinpy.utils.utils import calc_overlap
import numpy as np
from clinpy.assays.assay_base import Assay
from functools import partial


class Junctions(Assay):
    def __init__(self, db, genome):
        super.__init__(db, genome)

    def select(self, cohort=None, uniq=False, samples=None, df=True, filtered=True):
        """
        returns a dataframe of junctions this is different than the junction class
        instances returned from Sample.junctions
        :param cohort: name of the cohort to get
        :param uniq: return unique junctions not grouped by sample
        :param samples: return junctions from those samples
        :param df: return a dataframe otherwise Junction instance
        :param filtered: use the filtered junction table instead of all junctions
        :return: a dataframe of junctions and reads mapping to them, if uniq is
        selected mappings will not be there since different samples will have different
        numbers of reads mapping to each junction.
        """

        if filtered:
            table = Table("junctions", self.metadata, autoload=True, autoload_with=self.db)
            sample_to_junction = Table("sample_to_junction", self.metadata,
                                       autoload=True, autoload_with=self.db)
        else:
            table = Table("all_junctions", self.metadata, autoload=True, autoload_with=self.db)
            sample_to_junction = Table("sample_to_alljunction", self.metadata,
                                           autoload=True, autoload_with=self.db)

        query = select(table.c.chrom, table.c.start, table.c.end, table.c.strand)

        if not uniq:
            query = query.add_columns(sample_to_junction.c.samplename, sample_to_junction.c.uniq_map,
                                      sample_to_junction.c.multi_map)
            query = query.join(sample_to_junction, table.c.id == sample_to_junction.c.junction)

        if cohort is not None:
            sample_table = Table("samples", self.metadata, autoload=True, autoload_with=self.db)
            samples = select(sample_table.c.study_id).filter(sample_table.c.cohort.in_(cohort))

            junctions = select(sample_to_junction.c.junction).distinct(). \
                filter(sample_to_junction.c.samplename.in_(samples))

            query = query.filter(table.c.id.in_(junctions))

        if samples is not None:
            junctions = select(sample_to_junction.c.junction). \
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

    def search(self, gr, samples=None, unique=False, filtered=True):
        """
        search a given region for junctions
        :param gr: a pyranges
        :param samples: sample id, if none search all samples otherwise limit to the list of samples provided an interable
        :param  unique: whether to return unique junctions, sample infomation will not be present, otherwise return
        other data as well
        :return: a dataframe of junctions that are in a given region
        """
        if filtered:
            junctions = Table("junctions", self.metadata, autoload=True, autoload_with=self.db)
            sample_to_junction = Table("sample_to_junction", self.metadata,
                                           autoload=True, autoload_with=self.db)
        else:
            junctions = Table("all_junctions", self.metadata, autoload=True, autoload_with=self.db)
            sample_to_junction = Table("sample_to_alljunction", self.metadata,
                                           autoload=True, autoload_with=self.db)

        query = select(junctions.c.chrom, junctions.c.start, junctions.c.end, junctions.c.strand).filter(
            and_(junctions.c.chrom == str(gr.Chromosomes), junctions.c.strand == str(gr.Strand))).filter(
            and_(junctions.c.end >= int(gr.Start), junctions.c.start <= int(gr.End)))

        if samples is not None:
            sample_junctions = select(sample_to_junction.c.junction).filter(
                sample_to_junction.c.samplename.in_(samples))
            query = query.filter(junctions.c.id.in_(sample_junctions))

        if not unique:
            query = query.add_columns(sample_to_junction.c.samplename, sample_to_junction.c.uniq_map,
                                      sample_to_junction.c.multi_map)
            query = query.join(sample_to_junction, junctions.c.id == sample_to_junction.c.junction)

        results = self.session.execute(query).fetchall()
        columns = list(self.session.execute(query).keys())
        results = pd.DataFrame(results)
        results.columns = columns

        return results

    #TODO
    def filter(self, junc_func, skip_existing=True, **kwargs):
        """
        take a function, fill it's arguments and apply to each sample junctions, if the sample is already
        in the fitered junctions skip that sample if skip existing otherwise overwrite. The function needs to take
        a dataframe of known columns and return a dataframe (with fewer rows hopefully) with the same columns.

        The function can accept an arbitrary number of parameters but one of them must be named "df" this is the
        dataframe of unfiltered junctions. The parameters are passed normally and then handled by functools partial

        The columns are: chrom, start, end, strand, uniq_map, multi_map, samplename

        :param junc_func: a function takes a dataframe returns a dataframe
        :param skip_existing: if a sample is already in filtered junctions skip otherwise overwrite
        :param kwargs: arguments for the function except for the dataframe
        :return:
        """
        filled_func = partial(junc_func, kwargs)


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
        table = Table("genes", genome.metadata, autoload=True, autoload_with=genome)
        # find the start and end genes, there might be others in between
        query = select(table).filter(and_(table.c.chrom == self.chrom,
                                                  table.c.strand == self.strand)). \
            filter(or_(and_(table.c.start <= self.start, table.c.end >= self.start),
                           and_(table.c.start <= self.end, table.c.end >= self.end)))

        results=genome.session.execute(query).fetchall()
        if len(results)>0:
            results = pd.DataFrame(results)
            results.columns = list(genome.session.execute(query).keys())
        else:
            return None

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

        gene_table = Table("genes", genome.metadata, autoload=True, autoload_with=genome)
        tx_table = Table("transcripts", genome.metadata,
                             autoload=True, autoload_with=genome)

        subq = select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()
        query = select(tx_table, subq.c.chrom, subq.c.strand). \
            join(subq, tx_table.c.gene == subq.c.id). \
            filter(or_(and_(tx_table.c.start <= self.start, tx_table.c.end >= self.start),
                           and_(tx_table.c.start <= self.end, tx_table.c.end >= self.end)),
                   and_(subq.c.chrom == self.chrom, subq.c.strand == self.strand))

        if biotype is not None:
            query = query.where(tx_table.c.biotype.in_(biotype))

        results = genome.session.execute(query).fetchall()
        if len(results)>0:
            results = pd.DataFrame(results)
            results.columns = genome.session.execute(query).keys()
        else:
            return None

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

    def samples(self, project, tolerance=None, overlap=None, reciprocal=False, return_junctions=False):
        """
        find all the samples in the project that have the same junction with some overlap or tolerance in either end. if return_junction is
        True then a dataframe of sampleids, junction coords and mapping reads will be returned.
        :param project: project instance
        :param tolerance: a tuple how much from the 5' and 3' end is allowed in nucleo
        :param overlap: fraction overlap with self
        :param reciprocal: if True reciprocal overlap must be >= overlap value
        :param return_junctions: if True return a tuple with sample id and list of junctions that match
        :return:
        """
        junctions_table = Table("junctions", project.metadata, autoload=True, autoload_with=project)
        sample_to_junction = Table("sample_to_junction", project.metadata, autoload=True,
                                       autoload_with=project)
        subq = select(junctions_table.c.id).filter(and_(junctions_table.c.chrom == self.chrom,
                                                                junctions_table.c.strand == self.strand))
        if overlap is None:
            if tolerance is None:
                subq = subq.filter(and_(junctions_table.c.start == self.start,
                                            junctions_table.c.end == self.end)).scalar_subquery()
            else:
                subq = subq.filter(and_(junctions_table.c.start >= self.start - tolerance[0],
                                            junctions_table.c.end <= self.end + tolerance[1])).scalar_subquery()

            query = select(sample_to_junction.c.samplename, sample_to_junction.c.junction,
                               sample_to_junction.c.uniq_map, sample_to_junction.c.multi_map).where(sample_to_junction.c.junction.in_(subq))
            samples_junctions = project.session.execute(query).fetchall()
            columns = list(project.session.execute(query).keys())

        else:
            if tolerance is not None:
                raise NotImplementedError("you can specify tolerance OR overlap")

            else:
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
                sample_query = select(sample_to_junction.c.samplename, sample_to_junction.c.junction,
                                          sample_to_junction.c.uniq_map, sample_to_junction.c.multi_map).where(
                    sample_to_junction.c.junction.in_(overlapping))

                samples_junctions = project.session.execute(sample_query).fetchall()
                columns = list(project.session.execute(sample_query).keys())

        samples_junctions=pd.DataFrame(samples_junctions)
        samples_junctions.columns=columns

        if return_junctions:
            junc_ids=samples_junctions["junction"].drop_duplicates().to_list()
            junc_coords_q=select(junctions_table.c.id, junctions_table.c.chrom, junctions_table.c.start,
                                  junctions_table.c.end, junctions_table.c.strand).filter(junctions_table.c.id.in_(junc_ids))
            junc_coords=project.session.execute(junc_coords_q).fetchall()
            junc_coords=pd.DataFrame(junc_coords)
            junc_coords.columns = list(project.session.execute(junc_coords_q).keys())

            samples_junctions=samples_junctions.merge(junc_coords, how="left", left_on="junction", right_on="id")
            return samples_junctions

        else:
            return samples_junctions["samplename"].drop_duplicates().to_list()



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
        return "A junction in {} starting at {} and ending at {} on strand {}".format(self.chrom, self.start,
                                                                               self.end, self.strand)