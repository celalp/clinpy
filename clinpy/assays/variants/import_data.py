import gc
import os
import pandas as pd
import pysam
from functools import reduce
from sqlalchemy import Table, select
import yaml
from pysam import VariantFile

from clinpy.utils.utils import dict_to_table

# genotypes
genotypes={"(1, 1)": ("ALT", "ALT"),
           "(0, 1)": ("REF", "ALT"),
           "(1, 0)": ("ALT", "REF"),
           "(None, 1)": ("UNK", "ALT"),
           "(1, None)": ("ALT", "UNK"),
           "(None, None)": ("UNK", "UNK"),
           "(1)": ("ALT", "DEL")}


def parse_genotype(var, genotypes):
    """
    change the gt field so that there is a better way of displaying the data, I will include the dict here
    since I don't really see any other possible combinations of genotypes for snps
    """
    gt=genotypes[str(var["GT"])]
    if var.phased:
        sep="|"
    else:
        sep="/"
    gt = "{}{}{}".format(gt[0], sep, gt[1])
    return gt


def compare_fields(files, info_name="CSQ",  info_sep="|", not_same="error"):
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

    fields = [field.lower() for field in fields]
    formats = [format.lower() for format in formats]
    return fields, formats

# TODO there is an issue with some of the indel where VEP outputs 2 or more annotations in the frequency
# columns which then messes up type coersion i'm not sure if this is an issuel with mital vcfs or a general
# thing with VEP need to as nour/madeline this
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

# TODO properly identify genotype and phasing assume gt fields are reserved and throw and error if gt is not there
# TODO multisample vcf

def parse_vcf(file, fields, formats, type_dict, info_name="CSQ", info_sep="|", ignore=True, genotypes=genotypes):
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
    if genotypes is None:
        genotypes = genotypes
    header = file.header
    info = header.info[info_name]
    split_fields = info.description.split(info_sep)[1:]
    split_fields = [field.lower() for field in
                    split_fields]  # this is re-done in the vcf to compare with the union/intersection stuff

    variants = []
    for var in file:  # go over each variant
        var_details = {"chrom": var.chrom, "pos": var.pos, "id": var.id, "ref": var.ref, "alt": var.alts[0],
                       "qual": var.qual, "filter": var.filter.keys()[0]}  # these are mandatory vcf fields

        for item in var.samples[0].keys():  # assuming one sample per vcf, otherwise will get the first one not ideal
            if item=="GT":
                gt=parse_genotype(var.samples[0], genotypes)
                var_details["gt"]=gt
                #var_details[item.lower()] = str(var.samples[0][item])
            elif item.lower() in formats:
                var_details[item.lower()] = str(var.samples[0][item])
            else:
                continue

        consequences = []  # there will be a separate consequence for each transcript
        consqs = var.info[info_name]
        consqs = [cons.split(info_sep)[1:] for cons in
                  consqs]  # because we are removing the first one above this is vep specific not ideal
        for impact in consqs:  # go over each impact
            consequence = {}
            for i in range(len(split_fields)):
                if split_fields[i] in fields:
                    consequence[split_fields[i]] = impact[i]  # TODO refactor for better naming
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


def create_tables(params, project, create=True):
    tables = []
    for tablename in params.keys():
        tab=dict_to_table(params[tablename], tablename, project.metadata)
        id_cols = [str(col.name) for col in tab._columns if col.primary_key is True]
        tables.append([tablename, tab, id_cols])
        if create:
            tab.create()

    return tables


def add_to_variant_tables(project, fields, formats, impacts_table, mapping_table, variants_table, temp_table):
    """
    process the temp variants table and then remove it
    :param engine: slqalchemy connection
    :param meta: project db metadade
    :param session: session object
    :param create: is this a new project or added to an existing one
    :return:
    """
    project.metadata.reflect()
    table = variants_table
    mapping = mapping_table
    impacts = impacts_table

    temp_table = Table(temp_table, project.metadata, autoload=True, autoload_with=project.db)
    variants_table = Table(table, project.metadata, autoload=True, autoload_with=project.db)

    new_vars = select(temp_table.c.chrom, temp_table.c.pos, temp_table.c.id, temp_table.c.ref,
                      temp_table.c.alt).distinct().join(
        variants_table,
        (temp_table.c.chrom == variants_table.c.chrom) &
        (temp_table.c.pos == variants_table.c.pos) &
        (temp_table.c.ref == variants_table.c.ref) &
        (temp_table.c.alt == variants_table.c.alt), isouter=True
    ).add_columns(variants_table.c.variant_id).filter(variants_table.c.variant_id==None)

    new_vars = project.session.execute(new_vars).fetchall()

    if len(new_vars) > 0:
        new_vars = pd.DataFrame(new_vars)
        new_vars.columns = ["chrom", "pos", "id", "ref", "alt", "variant_id"]
        new_vars = new_vars.drop(columns=["variant_id"])
        new_vars.to_sql(table, project.db, index=False, if_exists="append")

    mapping_cols = ["samplename", "qual", "filter"] + formats
    # TODO this is very slow because of the join, would need to find a different way to get the data
    mapping_query = select(*[temp_table.c[col] for col in mapping_cols]).distinct().join(
        variants_table, (temp_table.c.chrom == variants_table.c.chrom) &
                        (temp_table.c.pos == variants_table.c.pos) &
                        (temp_table.c.ref == variants_table.c.ref) &
                        (temp_table.c.alt == variants_table.c.alt)
    ).add_columns(variants_table.c.variant_id)
    variant_mapping = pd.DataFrame(project.session.execute(mapping_query).fetchall())
    mapping_cols.append("variant_id")
    variant_mapping.columns = mapping_cols
    variant_mapping.to_sql(mapping, project.db, index=False, if_exists="append")
    del variant_mapping
    gc.collect()

    # this will be the variant impacts

    impacts_query = select(*[temp_table.c[col] for col in fields]).join(
        variants_table, (temp_table.c.chrom == variants_table.c.chrom) &
                        (temp_table.c.pos == variants_table.c.pos) &
                        (temp_table.c.ref == variants_table.c.ref) &
                        (temp_table.c.alt == variants_table.c.alt)
    ).add_columns(variants_table.c.variant_id)
    variant_impacts = pd.DataFrame(project.session.execute(impacts_query).fetchall())
    fields.append("variant_id")
    variant_impacts.columns = fields
    variant_impacts.to_sql(impacts, project.db, index=False, if_exists="append")

    del variant_impacts
    gc.collect()

    project.db.execute("drop table if exists temp_variants")


# this is specific to what i'm doing it's not flexible and will need some refactoring to be made more useful
def import_data(file, project, meta_read_fun, assay_params, create_assay=True):

    with open(assay_params["config"]) as y:
        config = yaml.safe_load(y)

    tables = create_tables(config["tables"], project, create_assay)
    mapping_file = meta_read_fun(file, **config["meta_read_fun_params"])

    info_name=config["parsing_params"]["info_name"]
    info_sep=config["parsing_params"]["info_sep"]

    all_vcfs=mapping_file[config["vcf_col"]].tolist()

    fields, formats=compare_fields(all_vcfs, info_name, info_sep ,
                                   config["parsing_params"]["not_same"])

    for index, row in mapping_file.iterrows():
        print(row["samplename"])
        vcf=VariantFile(row[config["vcf_col"]])
        type_dict= config["tables"][config["parsing_params"]["impact_table"]]
        variants=parse_vcf(vcf, fields, formats,type_dict,
                           info_name, info_sep, config["parsing_params"]["ignore"])
        variants["samplename"]=row[config["sample_col"]]
        variants.to_sql(config["temp_table"], project.db, index=False, if_exists="append")
        add_to_variant_tables(project, fields, formats, tables[2][0], tables[1][0], tables[0][0], config["temp_table"])





