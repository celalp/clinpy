import warnings

import pandas as pd
from sqlalchemy import Table, select, and_, func, distinct, insert

from clinpy.assays.variants.utils import compound_filter, apply_filter


class Variants:
    def __init__(self, project, variants_table, mapping_table, impacts_table,
                 filters_table, filtered_variants_table, name):
        """
         this is a simple class for interacting with the databse in terms of searching for variants, and
         filtering based on arbitrary criteria, this is not the variant class this is more about getting variants then
         looking at specifiic attributes of a variant
        :param project: Project class
        """
        self.project = project
        self.variants_table = Table(variants_table, self.project.metadata, autoload=True, autoload_with=self.project.db)
        self.mapping_table = Table(mapping_table, self.project.metadata, autoload=True, autoload_with=self.project.db)
        self.impacts_table = Table(impacts_table, self.project.metadata, autoload=True, autoload_with=self.project.db)
        self.filters_table = Table(filters_table, self.project.metadata, autoload=True, autoload_with=self.project.db)
        self.filtered_variants_table = Table(filtered_variants_table, self.project.metadata, autoload=True,
                                             autoload_with=self.project.db)
        self.assay_name = name,
        self.assay_type = "variants"

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

    # TODO need a more memory efficient way to get this done may be a temp table
    def filter_variants(self, filters=None, filter_name=None, save=True, df=False, return_results=True, rerun=False,
                        return_ops=False):
        """
        filter unfiltered variants, the filters need to be a dict where the key is the field name and value is what you are searching for
        for categorical need to provide a list for numerircal if you want exact values you can provide a list again if you want comparisons
        you would need to provide
        :param filters: a dictionary of filters they can be nested
        :param filter_name: a name for the filter as a reminder
        :param save: save results to "filtered variants table along with the filters
        :param df: return a dataframe of variants with fitered impacts othewise a list of variant class instances
        :param return_results: return the results of the filter
        :param return_ops: return available operators and quit
        :return: a dataframe of variants with all the applicable columns
        """

        if return_ops:
            print("available operators are '==', '>', '>=', '<', '<=', '!=', "
                  "'like', 'ilike', 'notlike', 'notilike', 'in', 'not_in' you can specify"
                  "null with the string null for in and not in you need to have a list of elements"
                  "the type check is hardcoded")
            return None

        else:
            if filter_name is None:
                raise ValueError("please provide a filter name")
            else:
                # check if filter is in the table
                query = select(self.filters_table).where(self.filters_table.c.name == filter_name)
                filt_res = self.project.session.execute(query).fetchall()
                if rerun:
                    if len(filt_res) == 0:
                        raise ValueError("No filters found with that name")

                    if filters is not None:
                        warnings.warn("during a re-run filters will be retrieved from the "
                                      "database, ignoring set filters")
                        filters = filt_res[0][2]
                else:
                    if len(filt_res) > 0:
                        raise ValueError("There is already a filter with that name, if you want to "
                                         "re-run the filtering please set rerun=True")

            # I dont want to have the variant_id column twice
            impact_cols = [col for col in self.list_impacts() if col != "variant_id"]

            if df:
                # this might be too much for a large study returning all the variant impacts
                query = select(self.variants_table).join(self.impacts_table,
                                                         self.variants_table.c.variant_id == self.impacts_table.c.variant_id). \
                    add_columns(*[self.impacts_table.c[col] for col in impact_cols])
            else:
                query = select(self.variants_table).join(self.impacts_table,
                                                         self.variants_table.c.variant_id == self.impacts_table.c.variant_id)

            parsed_filters = []
            for key in filters.keys():
                parsed = compound_filter(key, filters[key], self.project.metadata)
                parsed_filters.append(parsed)

            query = apply_filter(query, parsed_filters)
            results = self.project.session.execute(query).fetchall()

            if save:
                # sqlite and possibly other engines do not support returning clause (other than postgres)
                # this means I have to do another query to get the max (presumably what we
                # inserted). This would break if someone wants to uuids
                if rerun:
                    filter_id = filt_res[0][0]
                    # i have to go through this song and dance because while some newer versions of sqlite is
                    # are supporting upserts there is no way of knowing the user is using that version
                    # also the same for other databases. The only one I know that has been supporting upserts for a
                    # long time is postgres
                    existing = select(self.filtered_variants_table.c.variant_id). \
                        filter(self.filtered_variants_table.c.filter_id == filter_id)
                    existing = self.project.session.execute(existing).fetchall()
                    existing = [var[0] for var in existing]
                    results = [res for res in results if res[0] not in existing]
                else:
                    command = insert(self.filters_table).values(filter=filters, name=filter_name)
                    self.project.session.execute(command)
                    self.project.session.commit()
                    filter_id = select(func.max(self.filters_table.c.id))
                    filter_id = self.project.session.execute(filter_id).fetchone()[0]

                # TODO this assumes the first column is the variant id and there is no check for it
                filtered_ids = [(var[0], filter_id) for var in results]
                command = insert(self.filtered_variants_table).values(filtered_ids)
                self.project.session.execute(command)
                self.project.session.commit()
                print("Filtering complete and all the variants have been saved in the database")

            if return_results:
                if df:
                    filtered_variants = pd.DataFrame(results)
                    filtered_variants.columns = [desc["name"] for desc in query.column_descriptions]
                else:
                    filtered_variants = []
                    for variant in results:
                        var = Variant(variant_id=variant[0], chrom=variant[1], pos=variant[2],
                                      ref=variant[3], alt=variant[4], variants=self)
                        filtered_variants.append(var)

                return filtered_variants

    # TODO search only among the filtered variants
    def search_region(self, gr, samples=None, df=True, filtered=False, filter_id=None, include_impacts=False):
        """
        search a given genomic region for variants
        :param gr: a pyranges object strand is ignored if there
        :param samples search specific samples for variants
        :param df: return a DataFrame? Otherwise will return variant class instance
        :param filtered: limit results to filtered variants, if filter id is none will pick the latest one
        :param include_impacts: also add the variant impacts, this will significantly increase the # of rows
        returned since a variant can have several impacts
        :return: a dataframe of variants and impact, a variant may have more than one impact
        """
        format_cols = self.list_variant_quals()[1:]
        impact_cols = self.list_impacts()

        query = select(self.variants_table.c.variant_id, self.variants_table.c.chrom,
                       self.variants_table.c.pos, self.variants_table.c.ref,
                       self.variants_table.c.alt).filter(and_(self.variants_table.c.pos >= int(gr.Start),
                                                              self.variants_table.c.pos <= int(gr.End)),
                                                         self.variants_table.c.chrom == str(gr.chromosomes[0])). \
            join(self.mapping_table, self.mapping_table.c.variant_id == self.variants_table.c.variant_id). \
            add_columns(*[self.mapping_table.c[col] for col in format_cols])

        if samples is not None:
            query.filter(self.mapping_table.c.variant_id.in_(samples))

        if include_impacts:
            query = query.join(self.impacts_table, self.variants_table.c.variant_id == self.impacts_table.c.variant_id). \
                add_columns(*[self.impacts_table.c[col] for col in impact_cols])

        if filtered:
            if filter_id is None:
                print("No filter id is provided selecting the latest one")
                latest = select(func.max(self.filters_table.c.id))
                filter = select(self.filters_table.c.filter).where(self.filters_table.c.filter_id == latest)
                filter = self.project.session.execute(filter).fetchone()[0]
                filtered_variants = select(self.filtered_variants_table.c.variant_id). \
                    filter(self.filtered_variants_table.c.filter_id == filter)
                query = query.filter(self.variants_table.c.variant_id.in_(filtered_variants))
            else:
                filter = select(self.filters_table.c.filter).where(self.filters_table.c.filter_id == filter_id)
                filter = self.project.session.execute(filter).fetchall()
                if len(filter) == 0:
                    raise KeyError("No filters found with that filter please check your filter id")
                else:
                    filtered_variants = select(self.filtered_variants_table.c.variant_id). \
                        filter(self.filtered_variants_table.c.filter_id == filter)
                    query = query.filter(self.variants_table.c.variant_id.in_(filtered_variants))

        results = self.project.session.execute(query).fetchall()

        if len(results) > 0:
            if df:
                results = pd.DataFrame(results)
                colnames = self.project.session.execute(query).keys()
                results.columns = colnames
            else:
                variants = []
                for variant in results:
                    var = Variant(variant_id=variant[0], chrom=variant[1], pos=variant[2],
                                  ref=variant[3], alt=variant[4], variants=self)
                    variants.append(var)
                results = variants
        else:
            results = None

        return results

    def __str__(self):
        num_samples = select(func.count(distinct(self.mapping_table.c.samplename)))
        num_samples = self.project.session.execute(num_samples).fetchone()[0]

        num_variants = select(func.count(self.variants_table.c.variant_id))
        num_variants = self.project.session.execute(num_variants).fetchone()[0]

        num_impacts = select(func.count(self.impacts_table.c.variant_id))
        num_impacts = self.project.session.execute(num_impacts).fetchone()[0]

        desc = "{} variants from {} samples with {} different impact annotations".format(num_variants,
                                                                                         num_samples, num_impacts)
        return desc

    def to_variant(self, df):
        """
        convert a dataframe with variant information to a list of Variant classes
        """
        variant_list = [Variant(row[0], row[1], row[2], row[3], row[4], self) for row in
                        zip(df["variant_id"], df["chrom"], df["pos"], df["ref"], df["alt"])]
        return variant_list

    def show_filters(self):
        """
        show what kind of filters have been applied
        """
        query = select(self.filters_table)
        filters = self.project.session.execute(query).fetchall()
        if len(filters) == 0:
            print("No filters have been applied yet")
            filters = None
        else:
            filters = pd.DataFrame(filters)
            filters.columns = list(self.filters_table.columns.keys())

        return filters


class Variant:
    """
    same class will be used for both rna and dna variants
    this one is for learning about a variant and searching other samples for the same variant. You can have an iterable of these
    if you want. Like all other assays this is lazy to prevent unnecessary queries
    """

    def __init__(self, variants, variant_id, chrom, pos, ref, alt):
        self.id = variant_id
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.variants = variants

    def counts(self, samples=None):
        """
        return allele counts, number of homs and hets
        :param samples only search witihn these samples list (even if only one sample)
        : return a tuple of # hets, # homs, ac and af based on either the entire project or samples provided
        """

        num_hets = select(func.count(self.variants.mapping_table.c.variant_id)).filter(
            and_(self.variants.mapping_table.c.variant_id == self.id,
                 self.variants.mapping_table.c.gt == "(0, 1)"))

        num_homs = select(func.count(self.variants.mapping_table.c.variant_id)).filter(
            and_(self.variants.mapping_table.c.variant_id == self.id,
                 self.variants.mapping_table.c.gt == "(1, 1)"))

        if samples is not None:
            num_hets = num_hets.filter(self.variants.mapping_table.c.samplename.in_(samples))
            num_homs = num_homs.filter(self.variants.mapping_table.c.samplename.in_(samples))
            total_samples = len(samples)
        else:
            total_samples = select(func.count(self.variants.project.samples_table.c.sample_id))
            total_samples = int(self.variants.project.session.execute(total_samples).fetchone()[0])

        num_hets = self.variants.project.session.execute(num_hets).fetchone()[0]
        num_homs = self.variants.project.session.execute(num_homs).fetchone()[0]

        ac = num_hets + num_homs * 2
        af = ac / total_samples

        return (num_hets, num_homs, ac, af)

    def samples(self, in_samples=None, genotype="both"):
        """
        return samples containing the exact variant
        :param in_samples: an iterable with the samples to restrict the search to
        :param genotype: return everything, het or hom
        returns a list of sample ids
        """
        query = select(self.variants.mapping_table.c.samplename).filter(
            self.variants.mapping_table.c.variant_id == self.id)

        if genotype == "het":
            query = query.filter(self.variants.mapping_table.c.gt == "(0, 1)")
        elif genotype == "hom":
            query = query.filter(self.variants.mapping_table.c.gt == "(1, 1)")
        elif genotype not in ["hom", "het", "both"]:
            raise ValueError("genotype can only be hom, het or both")

        if in_samples is not None:
            query = query.filter(self.variants.mapping_table.c.samplename.in_(in_samples))

        results = self.variants.project.session.execute(query).fetchall()

        samples = [result[0] for result in results]

        return samples

    def impact(self):
        query = select(self.variants.impacts_table).filter(self.variants.impacts_table.c.variant_id == self.id)
        results = self.variants.project.session.execute(query).fetchall()
        results = pd.DataFrame(results)
        results.columns = self.variants.project.session.execute(query).keys()

        return results

    def __str__(self):
        return "A variant in {}:{} with from {} to {}".format(self.chrom, self.pos,
                                                              self.ref, self.alt)
