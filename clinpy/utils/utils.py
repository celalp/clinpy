from sqlalchemy import Table, Column, Integer, String, \
    Float, Date, Boolean, JSON, ForeignKey, BLOB, select
from sqlalchemy import create_engine
from functools import partial

def str_to_type(st, *args):
    """
    return sqlalchemy types based on a string, this is used to create tables dynamically
    :param st: a string
    :param args tablename and column name if you are creating a foreing key
    :return: a sqlachemy type supported types are str, int, float, date, bool and json there is
    limited support for json if not in the list return an error. you can create a foreign key if you provide the
    tablename and column name in *args in that order
    """
    if st == "int":
        coltype = Integer
    elif st == "str":
        coltype = String
    elif st == "float":
        coltype = Float
    elif st == "date":
        coltype = Date
    elif st == "bool":
        coltype = Boolean
    elif st == "json":
        coltype = JSON
    elif st == "blob":
        coltype = BLOB
    elif st == "fk":
        coltype = ForeignKey("{}.{}".format(args[0], args[1]))
    else:
        raise NotImplementedError("{} is not an implemented column type".format(st))
    return coltype


def dict_to_table(dict, tablename, meta):
    """
    create a sqlalchemy table from a dictionary, currently only a handful of things are supported
    The dictionanry structure is a as follows:
    {colname:
        type: str, int, bool, json, date (only one)
        pk: true
        fk:
            table:table for the foreign key
            column: column name for the foreignkey
    }

    all error checks (like if the table exists etc) are done by sqlalchemy

    :param dict:  a dictionary (see above)
    :param tablename: name of the table
    :param meta: table metadata so they can be added but will not be created until called table.create()
    :return: a sqlalchemy table with metadata (engine data) associated with it
    """
    colnames = list(dict.keys())
    coltypes = []
    for col in colnames:
        if dict[col]["type"] != "fk":
            coltypes.append(str_to_type(dict[col]["type"]))
        else:
            coltypes.append(str_to_type(dict[col]["type"], dict[col]["fk"]["table"], dict[col]["fk"]["column"]))
    idxs = [dict[colname]["index"] for colname in colnames]
    pks = [True if "pk" in dict[colname].keys() else False for colname in colnames]
    table = Table(tablename, meta,
                  *(Column(colname, coltype, primary_key=pk, index=idx)
                    for colname, coltype, pk, idx in zip(colnames, coltypes, pks, idxs))
                  , extend_existing=True)

    return table

def dict_to_engine(dbtype, name, username=None, pwd=None, host=None, port=None):
    """
    create a database connection based on the yaml description
    :param params: the output section of the config.py
    :param username of the database (not used in sqlite)
    :param pwd: the password for the user (not used in sqlite)
    :param host: the ip or url of the host (not used in sqlite)
    :param port: the port of the database in the host (not used in sqlite)
    :return: a sqlalchemy engine
    """
    if dbtype == "sqlite":
        dbstring = "sqlite:///{}".format(name)
    elif dbtype == "postgresql":
        dbstring = "postgresql://{}:{}@{}:{}/{}".format(username, pwd, host, port, name)
    elif dbtype == "mariadb":
        dbstring = "mariadb+mariadbconnector:://{}:{}@{}:{}/{}".format(username,pwd, host, port, name)
    elif dbtype == "mysql":
        dbstring = "mysql+mysqlconnector://{}:{}@{}:{}/{}".format(username, pwd, host, port, name)
    else:
        raise NotImplementedError("There is no current support for {}".format(dbtype))
    engine = create_engine(dbstring)
    return engine


#TODO test
def prep_read_fun(fun, **kwargs):
    """
    """
    read_fun=partial(fun, kwargs)
    return read_fun

