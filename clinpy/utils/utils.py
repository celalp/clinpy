from sqlalchemy import Table, Column, Integer, String, Float, Date, Boolean, JSON, ForeignKey
from sqlalchemy import or_, and_, not_, select
from sqlalchemy import create_engine, MetaData
import sqlalchemy.sql.operators as ops


def calc_overlap(int1, int2):
    """
    takes 2 tuples 0 is start 1 is end, assumes that they are in the same chromosome
    :param int1: interval 1
    :param int2: interval 2
    :return:
    """
    len1 = int1[1] - int1[0]
    len2 = int2[1] - int2[0]
    if int2[0] > int1[1] or int1[0] > int2[1]:
        return 0
    elif int2[0] <= int1[0] and int2[1] >= int1[1]:  # complete coverage
        return 1
    elif int2[0] >= int1[0] and int2[1] <= int1[1]:  # 1 covers 2
        return len2 / len1
    elif int2[0] >= int1[0] and int2[1] > int1[1]:  # 3' overlap
        return (int1[1] - int2[0]) / len1
    elif int2[0] <= int1[0] and int2[1] < int1[1]:  # 5' overlap
        return (int2[1] - int1[0]) / len1


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
            colname: column name for the foreignkey
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


def dict_to_engine(params, **kwargs):
    #TODO support postgres, mariadb and mysql
    # means adding an .env file to the create_project.py
    """
    create a database connection based on the yaml description
    :param dict: the output section of the config yaml
    :param kwargs: if type is not sqlite in this order username, password, host, port
    :return: a sqlalchemy engine
    """
    if params["type"] == "sqlite":
        dbstring = "sqlite:///{}".format(params["name"])
    else:
        raise NotImplementedError("currently only sqlite is supported")
        # dbstring = "{}://{}:{}@{}:{}/{}".format(kwargs, params["name"])
    engine = create_engine(dbstring)
    return engine

def str_to_op(st):
    if st=="==":
        op=ops.eq
    elif st == ">":
        op=ops.gt
    elif st == ">=":
        op=ops.ge
    elif st == "<":
        op=ops.lt
    elif st == "<=":
        op=ops.le
    elif st == "!=":
        op=ops.ne
    elif st == "like":
        op=ops.like_op
    elif st == "ilike":
        op=ops.ilike_op
    elif st == "in":
        op=ops.in_op
    elif st == "not_in":
        op=ops.not_in_op
    else:
        raise NotImplementedError("{} is not an available comparison".format(st))
    return op

def single_filter(tup, meta):
    """

    """
    table=meta.tables[tup[0]]
    col=table.c[tup[1]]
    op=str_to_op(tup[2])
    filter=op(col, tup[3])
    return filter

def compound_filter(dct, meta):
    pass