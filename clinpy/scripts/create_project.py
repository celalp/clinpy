#TODO need to come up with a way to get the default python and see if compatible

#! python3.9


import argparse as arg
from datetime import datetime
from pysam import VariantFile
import yaml

from sqlalchemy import MetaData
from sqlalchemy.orm import Session

from clinpy.database.rna_tables import *
from clinpy.database.snp_tables import *
from clinpy.utils.rna_functions import *
from clinpy.utils.snp_functions import *
from clinpy.utils.utils import dict_to_engine


if __name__ == "__main__":
    parser = arg.ArgumentParser(description='add to a project database with genome annotations junctions'
                                            'and expression data')
    parser.add_argument('-y', '--yaml', help="config.yaml file see readme for details", type=str, action="store",
                        default="config.yaml")
    args = parser.parse_args()

    if args.yaml is None:
        raise ValueError("no yaml file proveded")

    if not os.path.isfile(args.yaml):
        raise FileNotFoundError("count not find {}".format(args.yaml))

    with open(args.yaml) as y:
        params = yaml.safe_load(y)

    engine = dict_to_engine(params["output"])
    project_meta = MetaData(bind=engine)
    session = Session(engine)

    if params["output"]["create"]:
        if os.path.isfile(params["output"]["name"]):
            raise FileExistsError("database already exists")
        else:  # create the database
            project_meta.reflect()
            sample_table = dict_to_table(params["sample_meta"]["columns"], "samples", project_meta)
            sample_table.create()
            print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Creating database " + params["output"]["name"])
    else:
        # just get what kind of tables there
        project_meta.reflect()

    sample_meta = pd.read_csv(params["sample_meta"]["file"], header=0, sep="\t")
    file_cols = list(sample_meta.columns)
    file_cols = [col.lower() for col in file_cols]
    sample_meta.columns = file_cols


    # these are mandatory columns if they are not there there will be an error
    if "sample_id" not in file_cols:
        raise ValueError("sample_id column must be present in sample metadata")

    if "cohort" not in file_cols:
        raise ValueError("cohort column must be preset in sample metadata")

    # pandas will check for columns, if there are extra ones they will be omitted if missing ones there will be an error
    sample_meta = sample_meta[list(params["sample_meta"]["columns"].keys())]
    # if there are extra samples or duplicates this will fail with a unique constraint error
    sample_meta.to_sql("samples", engine, index=False, if_exists="append")

    mods = params["data"]["modalities"].keys()
    for modality in mods:
        if modality == "rna":
            file = params["data"]["modalities"][modality]["file"]
            columns = params["data"]["modalities"][modality]["columns"]
            files = pd.read_csv(file, header=0, sep="\t")

            if "samplename" not in files.columns:
                raise ValueError(
                    "samplename column is missing from {}".format(params["data"]["modalities"][modality]["file"]))

            if "rna_variants" in columns or "filtered_rna_variants" in columns:
                if "vcf_config" not in params["data"]["modalities"][modality].keys():
                    raise ValueError("did not provide a vcf config file")
                else:
                    with open(params["data"]["modalities"][modality]["vcf_config"]) as y:
                        vcf_params = yaml.safe_load(y)

            # go over the columns and create tables
            ProjectBase.metadata.reflect(engine)  # so we know the samples table is there
            for column in columns:
                if column == "samplename":
                    continue
                elif column == "gene_expression":
                    GeneExpression.__table__.create(engine)
                elif column == "isoform_expression":
                    TranscriptExpression.__table__.create(engine)
                elif column == "unfiltered_junctions":
                    AllJunctions.__table__.create(engine)
                    SampleToAllJunction.__table__.create(engine)
                elif column == "filtered_junctions":
                    FilteredJunctions.__table__.create(engine)
                    SampleToJunction.__table__.create(engine)
                elif column == "rna_variants":
                    RNAVariants.__table__.create(engine)
                    vcf_files = [file for file in files["rna_variants"].to_list() if not pd.isna(file)]
                    fields, formats = compare_fields(vcf_files, vcf_params["info"]["name"], vcf_params["not_same"],
                                                     vcf_params["info"]["sep"])
                    project_meta.reflect()

                    generate_variant_tables(vcf_params, fields, formats, project_meta,
                                            params["sample_meta"]["columns"]["sample_id"]["type"], rna=True,
                                            filtered=False)
                elif column == "filtered_rna_variants":
                    FilteredRNAVariants.__table__.create(engine)
                    vcf_files = [file for file in files["filtered_rna_variants"].to_list() if not pd.isna(file)]
                    filt_fields, filt_formats = compare_fields(vcf_files, vcf_params["info"]["name"],
                                                               vcf_params["not_same"],
                                                               vcf_params["info"]["sep"])
                    project_meta.reflect()

                    generate_variant_tables(vcf_params, filt_fields, filt_formats, project_meta,
                                            params["sample_meta"]["columns"]["sample_id"]["type"], rna=True,
                                            filtered=True)

            files = files.to_dict(orient="records")  # create a dict
            for data in files:

                for dat_type in data.keys():
                    if dat_type == "samplename":
                        sample = data["samplename"]
                    elif dat_type == "gene_expression" and not pd.isna(data[dat_type]):
                        print("[" + datetime.now().strftime(
                            "%Y/%m/%d %H:%M:%S") + "] " + "Adding gene expression for " + str(sample))
                        import_expression(data[dat_type], sample, engine, gene=True)

                    elif dat_type == "isoform_expression" and not pd.isna(data[dat_type]):
                        print("[" + datetime.now().strftime(
                            "%Y/%m/%d %H:%M:%S") + "] " + "Adding isoform expression for " + str(sample))
                        import_expression(data[dat_type], sample, engine, gene=False)

                    elif dat_type == "unfiltered_junctions" and not pd.isna(data[dat_type]):
                        print("[" + datetime.now().strftime(
                            "%Y/%m/%d %H:%M:%S") + "] " + "Adding junctions for " + str(sample))

                        import_temp_junction(data[dat_type], sample, engine,
                                             min_junc_reads=params["data"]["modalities"]["rna"]["min_junction_reads"],
                                             filtered=False)


                    elif dat_type == "filtered_junction" and not pd.isna(data[dat_type]):
                        print("[" + datetime.now().strftime(
                            "%Y/%m/%d %H:%M:%S") + "] " + "Adding filtered junctions for " + str(sample))

                        import_temp_junction(data[dat_type], sample,
                                             params["data"]["modalities"]["rna"]["min_junction_reads"],
                                             filtered=True)


                    elif dat_type == "rna_variants":
                        if not pd.isna(data[dat_type]):
                            print("[" + datetime.now().strftime(
                                "%Y/%m/%d %H:%M:%S") + "] " + "Adding variants for " + str(sample))
                            vcf = VariantFile(data[dat_type])
                            variants = parse_vcf(vcf, fields, formats,
                                                 type_dict=vcf_params["variant_impacts"],
                                                 field_name=vcf_params["info"]["name"],
                                                 field_split=vcf_params["info"]["sep"],
                                                 ignore=vcf_params["missing_impact"])
                            import_temp_variants(variants, sample, engine, filtered=False)


                    elif dat_type == "filtered_rna_variants":
                        if not pd.isna(data[dat_type]):
                            print("[" + datetime.now().strftime(
                                "%Y/%m/%d %H:%M:%S") + "] " + "Adding filtered variants for " + str(sample))
                            variants = parse_vcf(data[dat_type], fields, formats,
                                                 type_dict=vcf_params["variant_impacts"],
                                                 field_name=vcf_params["info"]["name"],
                                                 field_split=vcf_params["info"]["sep"],
                                                 ignore=vcf_params["missing_impact"])
                            import_temp_variants(variants, sample, engine, filtered=True)

                    else:
                        raise NotImplementedError(
                            "{} has not been implemented in clinpy you can create a feature request at "
                            "https://github.com/celalp/clinpy".format(
                                dat_type))

            # this is the second iteration, now that all the files are in the temp tables we can split
            # them and put them where they belong
            for column in columns:
                if column=="unfiltered_junctions":
                    print("[" + datetime.now().strftime(
                        "%Y/%m/%d %H:%M:%S") + "] " + "Normalizing unfiltered junction tables")
                    add_to_junction_tables(engine, project_meta, session=session, create=params["output"]["create"],
                                           filtered=False)
                elif column == "filtered_junctions":
                    print("[" + datetime.now().strftime(
                        "%Y/%m/%d %H:%M:%S") + "] " + "Normalizing filtered junction tables")
                    add_to_junction_tables(engine, project_meta, session=session, create=params["output"]["create"],
                                           filtered=True)
                elif column == "rna_variants":
                    print("[" + datetime.now().strftime(
                        "%Y/%m/%d %H:%M:%S") + "] " + "Normalizing rna variant tables")
                    add_to_variant_tables(engine, project_meta, session, fields, formats,
                                          create=params["output"]["create"], filtered=False,
                                          rna=True)
                elif column == "filtered_rna_variants":
                    print("[" + datetime.now().strftime(
                        "%Y/%m/%d %H:%M:%S") + "] " + "Normalizing filtered rna variant tables")
                    add_to_variant_tables(engine, project_meta, session, fields, formats,
                                          create=params["output"]["create"], filtered=True,
                                          rna=True)
                else:
                    continue

        elif modality == "snps":

            file = params["data"]["modalities"][modality]["file"]
            columns = params["data"]["modalities"][modality]["columns"]
            files = pd.read_csv(file, header=0, sep="\t")

            if "samplename" not in files.columns:
                raise ValueError(
                    "samplename columns is missing from {}".format(params["data"]["modalities"][modality]["file"]))

            if "vcf_config" not in params["data"]["modalities"][modality].keys():
                raise ValueError("did not provide a vcf config file")
            else:
                with open(params["data"]["modalities"][modality]["vcf_config"]) as y:
                    vcf_params = yaml.safe_load(y)

            for column in columns:
                ProjectBase.metadata.reflect(engine) # to update things as we go along
                if column == "samplename":
                    continue
                elif column == "variants":
                    Variants.__table__.create(engine)
                    vcf_files = [file for file in files["variants"].to_list() if not pd.isna(file)]
                    fields, formats = compare_fields(vcf_files, vcf_params["info"]["name"], vcf_params["not_same"],
                                                 vcf_params["info"]["sep"])
                    generate_variant_tables(vcf_params, fields, formats, project_meta,
                                        params["sample_meta"]["columns"]["sample_id"]["type"], rna=False,
                                        filtered=False)
                elif column == "filtered_variants":
                    FilteredVariants.__table__.create(engine)
                    vcf_files = [file for file in files["filtered_variants"].to_list() if not pd.isna(file)]
                    filt_fields, filt_formats = compare_fields(vcf_files, vcf_params["info"]["name"],
                                                           vcf_params["not_same"],
                                                           vcf_params["info"]["sep"])
                    generate_variant_tables(vcf_params, filt_fields, filt_formats, project_meta,
                                        params["sample_meta"]["columns"]["sample_id"]["type"], rna=False,
                                        filtered=True)

            files = files.to_dict(orient="records")  # create a dict
            for data in files:
                for dat_type in data.keys():
                    if dat_type == "samplename":
                        sample = data["samplename"]

                    elif dat_type == "variants":
                        if not pd.isna(data[dat_type]):
                            print("[" + datetime.now().strftime(
                            "%Y/%m/%d %H:%M:%S") + "] " + "Adding variants for " + str(sample))
                            vcf=VariantFile(data[dat_type])
                            variants = parse_vcf(vcf, fields, formats,
                                             type_dict=vcf_params["variant_impacts"],
                                             field_name=vcf_params["info"]["name"],
                                             field_split=vcf_params["info"]["sep"],
                                             ignore=vcf_params["missing_impact"])
                            import_temp_variants(variants, sample, engine, filtered=False)


                    elif dat_type == "filtered_variants":
                        if not pd.isna(data[dat_type]):
                            print("[" + datetime.now().strftime(
                            "%Y/%m/%d %H:%M:%S") + "] " + "Adding filtered variants for " + str(sample))
                            variants = parse_vcf(data[dat_type], fields, formats,
                                             type_dict=vcf_params["variant_impacts"],
                                             field_name=vcf_params["info"]["name"],
                                             field_split=vcf_params["info"]["sep"],
                                             ignore=vcf_params["missing_impact"])
                            import_temp_variants(variants, sample, engine, filtered=True)

                    else:
                        raise NotImplementedError(
                            "{} has not been implemented in clinpy you can create a feature request at "
                            "https://github.com/celalp/clinpy".format(
                                type))

            for column in columns:
                if column == "variants":
                    add_to_variant_tables(engine, project_meta, session, fields, formats,
                                      create=params["output"]["create"], filtered=False,
                                      rna=False)
                elif column == "filtered_variants":
                    add_to_variant_tables(engine, project_meta, session, fields, formats,
                                  create=params["output"]["create"], filtered=True,
                                  rna=False)
                else:
                    continue

        else:  # the list will increase as time goes on
            raise NotImplementedError(
                "{} has not been implemented in clinpy you can create a feature request at "
                "https://github.com/celalp/clinpy".format(
                    modality))

    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Done!")
