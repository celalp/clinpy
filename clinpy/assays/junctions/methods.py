import pandas as pd
from sqlalchemy import Table, select, and_, or_
import pyranges
import numpy as np
from functools import partial
from clinpy.assays.junctions.utils import calc_overlap


class Junctions:
    def __init__(self, project, junction_table, mapping_table, name):
        self.project = project
        self.junction_table = Table(junction_table, self.project.metadata, autoload=True,
                                            autoload_with=self.project.db)
        self.mapping_table = Table(mapping_table, self.project.metadata,
                                   autoload=True, autoload_with=self.project.db)
        self.assay_name=name,
        self.assay_type="junctions"

    def select(self, uniq=False, samples=None, df=True):
        """
        returns a dataframe of junctions this is different than the junction class
        instances returned from Sample.junctions
        :param uniq: return unique junctions not grouped by sample
        :param samples: return junctions from those samples
        :param df: return a dataframe otherwise Junction instance
        :return: a dataframe of junctions and reads mapping to them, if uniq is
        selected mappings will not be there since different samples will have different
        numbers of reads mapping to each junction.
        """

        query = select(self.junction_table.c.chrom, self.junction_table.c.start, self.junction_table.c.end,
                       self.junction_table.c.strand)

        if not uniq:
            query = query.add_columns(self.mapping_table.c.samplename, self.mapping_table.c.uniq_map,
                                      self.mapping_table.c.multi_map)
            query = query.join(self.mapping_table, self.junction_table.c.id == self.mapping_table.c.junction)

        if samples is not None:
            junctions = select(self.mapping_table.c.junction). \
                filter(self.mapping_table.c.samplename.in_(samples)).distinct()

            query = query.filter(self.junction_table.c.id.in_(junctions))

        results = self.project.session.execute(query).fetchall()
        results = pd.DataFrame(results)
        results.columns = list(self.project.session.execute(query).keys())

        if df:
            return results
        elif not df and not uniq:
            juncs = []
            for chrom, start, end, strand, uniq_map, multi_map in zip(results.chrom, results.start, results.end,
                                                                      results.strand, results.uniq_map, results.multi_map):
                juncs.append(Junction(chrom, start, end, strand, uniq_map, multi_map, self))
            return juncs
        else:
            raise NotImplementedError("returning unique junctions as a junction class is not implemented")

    def search(self, gr, samples=None, unique=False):
        """
        search a given region for junctions
        :param gr: a pyranges
        :param samples: sample id, if none search all samples otherwise limit to the list of samples provided an interable
        :param  unique: whether to return unique junctions, sample infomation will not be present, otherwise return
        other data as well
        :return: a dataframe of junctions that are in a given region
        """

        query = select(self.junction_table.c.chrom, self.junction_table.c.start, self.junction_table.c.end,
                       self.junction_table.c.strand).filter(
            and_(self.junction_table.c.chrom == str(gr.Chromosome[0]),
                 self.junction_table.c.strand == str(gr.Strand[0]))).filter(
            or_(self.junction_table.c.end >= int(gr.Start[0]), self.junction_table.c.start <= int(gr.End[0])))

        if samples is not None:
            sample_junctions = select(self.mapping_table.c.junction).filter(
                self.mapping_table.c.samplename.in_(samples))
            query = query.filter(self.junction_table.c.id.in_(sample_junctions))

        if not unique:
            query = query.add_columns(self.mapping_table.c.samplename, self.mapping_table.c.uniq_map,
                                      self.mapping_table.c.multi_map)
            query = query.join(self.mapping_table, self.junction_table.c.id == self.mapping_table.c.junction)

        results = self.project.session.execute(query).fetchall()
        if len(results) > 0:
            columns = list(self.project.session.execute(query).keys())
            results = pd.DataFrame(results)
            results.columns = columns
            return results
        else:
            print("Did not find any junctions in that range")

    # TODO
    def filter(self, junc_func, add_to_db=True, **kwargs):
        """
        take a function, fill its arguments and apply to each sample junctions, if the sample is already
        in the fitered junctions skip that sample if skip existing otherwise overwrite. The function needs to take
        a dataframe of known columns and return a dataframe (with fewer rows hopefully) with the same columns.

        The function can accept an arbitrary number of parameters but one of them must be named "df" this is the
        dataframe of unfiltered junctions. The parameters are passed normally and then handled by functools partial

        The columns are: chrom, start, end, strand, uniq_map, multi_map, samplename

        :param junc_func: a function takes a dataframe returns a dataframe
        :param kwargs: arguments for the function except for the dataframe
        :return:
        """
        filled_func = partial(junc_func, kwargs)


    def show_filters(self):
        pass

    def to_junction(self, df):
        junction_list=[Junction(row[0], row[1], row[2], row[3], row[4], row[5], self) for row in
                        zip(df["chrom"], df["start"], df["end"], df["strand"], df["uniq_map"],
                            df["multi_map"])]
        return junction_list


class Junction:
    def __init__(self, chrom, start, end, strand, uniq_map, multi_map, junctions):
        """
        initiate a junction instance bare bones
        :param chrom: chrom
        :param start: start
        :param end: end
        :param strand: strand
        :param uniq_map: uniq map reads
        :param multi_map: multi_map reads
        :param junctions a junctions instance with all the methods and tables, this is not a superclass it is just an
        instance
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.uniq_map = uniq_map
        self.multi_map = multi_map
        self.junctions = junctions
        self.project=self.junctions.project

    def genes(self, df=False):
        """
        return overapping genes
        :param df return a dataframe or pyranges
        :return: genes that match the start and end of the junction and if there are other genes in the middle
        """
        genes_table = Table("genes", self.junctions.project.genome.metadata, autoload=True,
                            autoload_with=self.junctions.project.genome.metadata)
        # find the start and end genes, there might be others in between
        query = select(genes_table).filter(and_(genes_table.c.chrom == self.chrom,
                                                genes_table.c.strand == self.strand)). \
            filter(or_(and_(genes_table.c.start <= self.start, genes_table.c.end >= self.start),
                       and_(genes_table.c.start <= self.end, genes_table.c.end >= self.end)))

        results = self.junctions.project.genome.session.execute(query).fetchall()
        if len(results) > 0:
            results = pd.DataFrame(results)
            results.columns = list(self.junctions.project.genome.session.execute(query).keys())
        else:
            print("This junction does not correspond to any genes")

        if df:
            return results
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=results.iloc[:, 1],
                                            starts=results.iloc[:, 2],
                                            ends=results.iloc[:, 3],
                                            strands=results.iloc[:, 4])

            gr.id = results["id"]
            return gr

    def transcripts(self, df=False):
        """
        get the transcript the junction overlaps
        :param df: return a dataframe or pyranges
        :return: a pyranges or dataframe with all the transcript that match the description.
        """

        gene_table = Table("genes", self.junctions.project.genome.metadata, autoload=True,
                           autoload_with=self.junctions.project.genome)
        tx_table = Table("transcripts", self.junctions.project.genome.metadata,
                         autoload=True, autoload_with=self.junctions.project.genome)

        subq = select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()
        query = select(tx_table, subq.c.chrom, subq.c.strand). \
            join(subq, tx_table.c.gene == subq.c.id). \
            filter(or_(and_(tx_table.c.start <= self.start, tx_table.c.end >= self.start),
                       and_(tx_table.c.start <= self.end, tx_table.c.end >= self.end)),
                   and_(subq.c.chrom == self.chrom, subq.c.strand == self.strand))

        results = self.junctions.project.genome.session.execute(query).fetchall()
        if len(results) > 0:
            results = pd.DataFrame(results)
            results.columns = self.junctions.project.genome.session.execute(query).keys()
        else:
            print("This junction does not correspond to any transcripts")

        if df:
            return results
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=results.iloc[:, 5],
                                            starts=results.iloc[:, 2],
                                            ends=results.iloc[:, 3],
                                            strands=results.iloc[:, 6])
            gr.tx_id = results.id
            return gr

    def features(self, tx):
        """
        return which intron/exon the transcript matches to
        :param genome: genome instance
        :param txs: names of the transcripts because an exon/intron might belong to multiple transccript
        :return: a dict with keys "Start" and "End" values are tuples first one is either "intron" or "exon" the second
        one is a df row with all the details
        """
        start_interval = pd.Interval(self.start, self.start, closed="both")
        end_interval = pd.Interval(self.end, self.end, closed="both")

        introns = self.junctions.project.genome.introns(names=tx, df=True, level="transcript")
        intron_intervals = pd.arrays.IntervalArray.from_arrays(introns.start.to_numpy(), introns.end.to_numpy(),
                                                               closed="both")
        exons = self.junctions.project.genome.exons(names=tx, df=True)
        exon_intervals = pd.arrays.IntervalArray.from_arrays(exons.start.to_numpy(), exons.end.to_numpy(),
                                                             closed="both")
        intervals = {"intron": intron_intervals, "exon": exon_intervals}
        dfs = {"intron": introns, "exon": exons}

        start_type = None
        end_type = None

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

        features_dict = {}
        if start_type is not None:
            features_dict["Start"] = (start_type, dfs[start_type].iloc[start_loc, :])
        else:
            features_dict["Start"] = (None, None)

        if end_type is not None:
            features_dict["End"] = (end_type, dfs[end_type].iloc[end_loc, :])
        else:
            features_dict["End"] = (None, None)

        return features_dict

    def samples(self, tolerance=None, overlap=None, reciprocal=False, return_junctions=False):
        """
        find all the samples in the project that have the same junction with some overlap or tolerance in either end.
        if return_junction is
        True then a dataframe of sampleids, junction coords and mapping reads will be returned.
        :param tolerance: a tuple how much from the 5' and 3' end is allowed in nucleo
        :param overlap: fraction overlap with self
        :param reciprocal: if True reciprocal overlap must be >= overlap value
        :param return_junctions: if True return a tuple with sample id and list of junctions that match
        :return:
        """
        junctions_table = Table(self.junctions.junction_table, self.junctions.project.metadata, autoload=True,
                                autoload_with=self.project.db)
        sample_to_junction = Table(self.junctions.mapping_table, self.junctions.project.metadata, autoload=True,
                                   autoload_with=self.junctions.project.db)
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
                           sample_to_junction.c.uniq_map, sample_to_junction.c.multi_map).where(
                sample_to_junction.c.junction.in_(subq))
            samples_junctions = self.junctions.project.session.execute(query).fetchall()
            columns = list(self.junctions.project.session.execute(query).keys())

        else:
            if tolerance is not None:
                raise NotImplementedError("you can specify tolerance OR overlap")

            else:
                query = subq.add_columns(junctions_table.c.start, junctions_table.c.end)
                junctions = self.junctions.project.session.execute(query).fetchall()
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

                samples_junctions = self.junctions.project.session.execute(sample_query).fetchall()
                columns = list(self.junctions.project.session.execute(sample_query).keys())

        samples_junctions = pd.DataFrame(samples_junctions)
        samples_junctions.columns = columns

        if return_junctions:
            junc_ids = samples_junctions["junction"].drop_duplicates().to_list()
            junc_coords_q = select(junctions_table.c.id, junctions_table.c.chrom, junctions_table.c.start,
                                   junctions_table.c.end, junctions_table.c.strand).filter(
                junctions_table.c.id.in_(junc_ids))
            junc_coords = self.junctions.project.session.execute(junc_coords_q).fetchall()
            junc_coords = pd.DataFrame(junc_coords)
            junc_coords.columns = list(self.junctions.project.session.execute(junc_coords_q).keys())

            samples_junctions = samples_junctions.merge(junc_coords, how="left", left_on="junction", right_on="id")
            return samples_junctions

        else:
            return samples_junctions["samplename"].drop_duplicates().to_list()

    def new_transcript(self, features, sequence=True, type="nuc"):
        """
        return the new transcript sequence assuming all junctions are the same
        :param features: a features dict from the features method
        :param sequence: if true return sequence otherwise return df
        :param type: if sequence return either nuc or aa
        :return: returns either a dataframe or sequence of the new transcript
        """
        start_exons = self.junctions.genome.exons(features["Start"][1].transcript.to_list(), df=True)
        end_exons = self.junctions.genome.exons(features["End"][1].transcript.to_list(), df=True)

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
            new_exons = self.junctions.genome.get_sequence(new_exons_gr, type=type, concat=True)

        return new_exons

    def __str__(self):
        """
        print method
        :return: string
        """
        return "A junction in {} starting at {} and ending at {} on strand {}".format(self.chrom, self.start,
                                                                                      self.end, self.strand)