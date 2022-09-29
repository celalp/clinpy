import gc
import os
from functools import reduce

import pandas as pd
import pysam
from sqlalchemy import Table, select

from clinpy.utils.utils import dict_to_table


def compare_fields(files, info_name="CSQ", not_same="error", info_sep="|"):
    """
    gets a list of vcf file paths and compares the descriptions of consequence info, if they are not all the same
    either throws an error, gets their union, or their intersection
    :param files: list of vcf paths
    :param field_name: info fileld name to investigate
    :param not_same: what to do if they are not identical, "error", "union", "intersection"
    :return: list of fields to get
    """
    descriptions = []
    formats = []
    for file in files:
        if not os.path.isfile(file):
            raise FileNotFoundError("{} file does not exist".format(file))
        else:
            vcf = pysam.VariantFile(file)
            header = vcf.header
            info = header.info[info_name]
            descriptions.append(info.description)
            formats.append(list(header.formats.keys()))

    if not all(descriptions) or not all(formats):
        if not_same == "error":
            raise ValueError("the descriptions of fields or formats are not the same in all files")
        elif not_same == "union":
            split_fields = [desc.split(info_sep)[1:] for desc in
                            descriptions]  # with vep the first one is allele we already have that
            fields = list(reduce(set.union, [set(item) for item in split_fields]))
            formats = list(reduce(set.union, [set(item) for item in formats]))
        elif not_same == "intersection":
            split_fields = [desc.split(info_sep)[1:] for desc in descriptions][1:]
            fields = list(reduce(set.intersection, [set(item) for item in split_fields]))
            formats = list(reduce(set.intersection, [set(item) for item in formats]))
        else:
            raise NotImplementedError("not_same field can only be 'error', 'union' or 'intersection'")
    else:
        fields = descriptions[0].split(info_sep)[1:]
        formats = formats[0]

    return fields, formats


def coerce(effects, type_dict):
    """
    take a description from config dict and convert them to appropriatie types
    :param effects: effects for a variant
    :param type_dict: type dict from config
    :return: return the same dict but values converted
    """
    new_var = {}
    for key in effects.keys():
        if key not in type_dict.keys():
            continue
        else:
            if effects[key] == '':
                new_var[key] = None  # type independent
            else:
                if type_dict[key]["type"] == "str":
                    new_var[key] = effects[key]
                elif type_dict[key]["type"] == "int":
                    new_var[key] = int(effects[key])
                elif type_dict[key]["type"] == "bool":
                    new_var[key] = bool(effects[key])
                elif type_dict[key]["type"] == "float":
                    new_var[key] = float(effects[key])
    return new_var


def parse_vcf(file, fields, formats, type_dict, field_name="CSQ", field_split="|", ignore=True):
    """
    get the consequences of variants based on the field name
    :param file: vcf file connection via pysam
    :param fields: output of compare fields
    :param formats: output of compare fileds
    :param type_dict the dict of data types and whether to index from vcf.yaml
    :param field_name: name of the variant consequence field default is CSQ for VEP
    :param field_split: split char of the specifiic info
    :return: a dataframe this is to be inserted as a temp table and then further processed like junctions
    """
    header = file.header
    info = header.info[field_name]
    split_fields = info.description.split(field_split)[1:]
    split_fields = [field.lower() for field in split_fields]
    fields = [field.lower() for field in fields]
    variants = []
    for var in file:  # go over each variant
        var_details = {"chrom": var.chrom, "pos": var.pos, "id": var.id, "ref": var.ref, "alt": var.alts[0],
                       "qual": var.qual, "filter": var.filter.keys()[0]}  # these are mandatory vcf fields

        for item in var.samples[0].keys():  # assuming one sample per vcf, otherwise will get the first one not ideal
            if item in formats:
                var_details[item.lower()] = str(var.samples[0][item])

        consequences = []  # there will be a separate consequence for each transcript
        consqs = var.info[field_name]
        consqs = [cons.split(field_split)[1:] for cons in
                  consqs]  # because we are removing the first one above this is vep specific not ideal
        for impact in consqs:  # go over each impact
            consequence = {}
            for i in range(len(split_fields)):
                if split_fields[i] in fields:
                    consequence[split_fields[i]] = impact[
                        i]  # having a hard time with variable names TODO refactor for better namin
                else:
                    if ignore:
                        continue
                    else:
                        raise ValueError("{} is not defined in the vcf config".format(impact))
            consequences.append(coerce(consequence, type_dict))
        var_details["CSQ"] = consequences
        variants.append(var_details)

    variants = pd.DataFrame(variants)
    variants = variants.explode("CSQ", ignore_index=True)
    consequence_df = pd.json_normalize(variants["CSQ"].to_list())
    variants = pd.concat([variants, consequence_df], axis=1)
    variants = variants.drop(columns="CSQ")
    return variants


# TODO foreign key constraint on sample id and variant id
def generate_variant_tables(vcf_params, fields, formats, meta, name_type=str, rna=False, filtered=False):
    """

    :param vcf_params: vcf yaml containing all the information about the vcf file
    :param fields: output from compare fields
    :param formats: ouptut from compare fiields
    :param meta: sqlachemy metadata for the project db
    :param name_type: name type of the sample id (usually str or int)
    :param rna: are these variants from RNA-Seq
    :param filtered: are these filtered variants
    :return: nothing creates the necessary tables in the project database and quits
    """
    impacts = {}
    for field in fields:
        field = field.lower()
        if field in vcf_params["variant_impacts"].keys():
            impacts[field] = vcf_params["variant_impacts"][field]
        else:
            continue
    impacts["variant_id"] = {"type": "int", "index": True}
    impact_name = "variant_impacts"
    samples_name = "sample_variants"
    if rna:
        impact_name = "rna_" + impact_name
        samples_name = "rna_" + samples_name
    if filtered:
        impact_name = "filtered_" + impact_name
        samples_name = "filtered_" + samples_name

    impacts_table = dict_to_table(impacts, impact_name, meta)
    impacts_table.create()

    sample_variants = {}
    sample_variants["variant_id"] = {"type": "int", "index": True, "pk": True}
    sample_variants["samplename"] = {"type": name_type,
                                     "index": True, "pk": True}
    sample_variants["qual"] = {"type": "float", "index": False}
    sample_variants["filter"] = {"type": "str", "index": True}

    for format in formats:
        if format == "GT":
            sample_variants[format] = {"type": "str", "index": True}
        else:
            sample_variants[format] = {"type": "str", "index": False}

    sample_variants_table = dict_to_table(sample_variants, samples_name, meta)
    sample_variants_table.create()


def import_temp_variants(variants, samplename, engine, filtered=False):
    """
    insert into a temp table to further process the variants table
    :param variants: dataframe of variants, output of parse_vcf
    :param samplename: name of the sample
    :param engine sqlalchemy connection
    :param filtered are these filtered variants
    :return:
    """
    variants["samplename"] = samplename
    if filtered:
        variants.to_sql("temp_filt_variants", engine, if_exists="append", index=False)
    else:
        variants.to_sql("temp_all_variants", engine, if_exists="append", index=False)


def add_to_variant_tables(engine, meta, session, fields, formats, create=True, filtered=False, rna=False):
    """
    process the temp variants table and then remove it
    :param engine: slqalchemy connection
    :param meta: project db metadade
    :param session: session object
    :param create: is this a new project or added to an existing one
    :param filtered: are these filtered variants?
    :param rna is this from rnaseq data
    :return:
    """

    fields = [field.lower() for field in fields]
    formats = [format.lower() for format in formats]

    meta.reflect()
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
        temp = "temp_filt_variants"
    else:
        temp = "temp_all_variants"

    temp_table = Table(temp, meta, autoload=True, autoload_with=engine)
    variants_table = Table(table, meta, autoload=True, autoload_with=engine)

    if create:
        distinct_vars = select(temp_table.c.chrom, temp_table.c.pos, temp_table.c.id, temp_table.c.ref,
                               temp_table.c.alt).distinct()
        distinct_vars = pd.DataFrame(session.execute(distinct_vars).fetchall())
        distinct_vars.columns = ["chrom", "pos", "id", "ref", "alt"]
        distinct_vars.to_sql(table, engine, index=False, if_exists="append")
        del distinct_vars
        gc.collect()
    else:
        new_vars = select(temp_table.c.chrom, temp_table.c.pos, temp_table.c.id, temp_table.c.ref,
                          temp_table.c.alt).distinct().join(
            variants_table,
            (temp_table.c.chrom == variants_table.c.chrom) &
            (temp_table.c.pos == variants_table.c.pos) &
            (temp_table.c.ref == variants_table.c.ref) &
            (temp_table.c.alt == variants_table.c.alt)
        ).filter(variants_table.c.variant_id is None)

        new_vars = session.execute(new_vars)

        if len(new_vars) > 0:
            new_vars = pd.DataFrame(new_vars)
            new_vars.columns = ["chrom", "pos", "id", "ref", "alt", "variant_id"]
            new_vars = new_vars.drop(columns=["id"])
            new_vars.to_sql(table, engine, index=False, if_exists="append")

    mapping_cols = ["samplename", "qual", "filter"] + formats
    # TODO this is very slow because of the join, would need to find a different way to get the data
    mapping_query = select(*[temp_table.c[col] for col in mapping_cols]).distinct().join(
        variants_table, (temp_table.c.chrom == variants_table.c.chrom) &
                        (temp_table.c.pos == variants_table.c.pos) &
                        (temp_table.c.ref == variants_table.c.ref) &
                        (temp_table.c.alt == variants_table.c.alt)
    ).add_columns(variants_table.c.variant_id)
    variant_mapping = pd.DataFrame(session.execute(mapping_query).fetchall())
    mapping_cols.append("variant_id")
    variant_mapping.columns = mapping_cols
    variant_mapping.to_sql(mapping, engine, index=False, if_exists="append")
    del variant_mapping
    gc.collect()

    # this will be the variant impacts

    impacts_query = select(*[temp_table.c[col] for col in fields]).join(
        variants_table, (temp_table.c.chrom == variants_table.c.chrom) &
                        (temp_table.c.pos == variants_table.c.pos) &
                        (temp_table.c.ref == variants_table.c.ref) &
                        (temp_table.c.alt == variants_table.c.alt)
    ).add_columns(variants_table.c.variant_id)
    variant_impacts = pd.DataFrame(session.execute(impacts_query).fetchall())
    fields.append("variant_id")
    variant_impacts.columns = fields
    variant_impacts.to_sql(impacts, engine, index=False, if_exists="append")

    del variant_impacts
    gc.collect()

    engine.execute("drop table if exists {}".format(temp))
