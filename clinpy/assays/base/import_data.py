from sqlalchemy import select, Table
from clinpy.utils.utils import dict_to_table
from datetime import datetime
import yaml



def isin(list1, list2, tablename):
    """
    This function is there to check columns between the yaml file and the actual file
    """
    if len(list2) > 0:
        for item in list1:
            if item in list2:
                raise ValueError("{} is already in {}".format(item, tablename))
            else:
                continue


def import_meta(read_fun, excel, project, table, sheet, id_cols, read_fun_params):
    data=read_fun(excel, sheet_name=sheet, comment="#", **read_fun_params)
    if data[id_cols].duplicated().sum() > 0: #this also serves as column check
        raise ValueError("There are duplicates in {} id please check your file".format(sheet))

    if table in project.metadata.tables.keys(): #there is a table check for duplicates
        data_table=Table(table, project.metadata, autoload=True, autoload_with=project.db)
        query = select(*[data_table.c[col] for col in id_cols])

        current = project.session.execute(query).fetchall()
        if data.shape[0] == 0:
            print("There are no {} specified but there are existing {} in the database".format(sheet, sheet))
        else:
            isin(data[id_cols].values.tolist(), [dat[0] for dat in current], table)
    else: #there is no table
        if data.shape[0] > 0:
            raise ValueError("There are no {} specified in the file and there are no {} in "
                         "database".format(sheet, sheet))

    data.to_sql(table, project.db, index=False, if_exists="append")
    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "{} have been imported".format(sheet))

def create_tables(params, project, create=True):
    """
    This creates an issue while in newer pythons the dict remember order so this will have to be
    python >3.6 which is fine since 3.7 is being deprecated.

    But also raises the issue on the end user where you need to define tables in the order you want to
    create them insertion will happen in the same order

    :param params: dict of descriptions from the config yaml file
    :param project: project class instance
    :return a list of dicts there the key is the table name and value is the table object and pk_cols

    """
    tables=[]
    for tablename in params.keys():
        tab=dict_to_table(params[tablename], tablename, project.metadata)
        id_cols=[str(col.name) for col in tab._columns if col.primary_key is True]
        tables.append([tablename, tab, id_cols])
        if create:
            tab.create()

    return tables


def import_data(file, project, read_fun, assay_params, read_fun_params, create_assay=True):
    """
    :param file: this is the excel file form base["filename"]
    :param project: project class instance
    :param read_fun: fuction to read file not string but the function itself
    :param read_fun_params: parameters for the read_function
    :param assay_params: table descriptions for the base assay
    """
    with open(assay_params["config"]) as y:
        config=yaml.safe_load(y)



    tables=create_tables(config, project, create_assay)

    for table in tables:
        import_meta(read_fun, file, project, table[0], table[0], table[2], read_fun_params)
