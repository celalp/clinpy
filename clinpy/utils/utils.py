from sqlalchemy import Table, Column, Integer, String, Float, Date, Boolean, JSON
from sqlalchemy import create_engine

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


def str_to_type(st):
    """
    return sqlalchemy types based on a string, this is used to create tables dynamically
    :param st: a string
    :return: a sqlachemy type supported types are str, int, float, date, bool and json there is
    limited support for json if not in the list return an error
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
    else:
        raise NotImplementedError("{} is not an implemented column type".format(st))
    return coltype


def dict_to_table(dict, tablename, meta):
    colnames = list(dict.keys())
    coltypes = [str_to_type(dict[colname]["type"]) for colname in colnames]
    idxs = [dict[colname]["index"] for colname in colnames]
    pks = [True if "pk" in dict[colname].keys() else False for colname in colnames]

    table = Table(tablename, meta,
                  *(Column(colname, coltype, primary_key=pk, index=idx)
                    for colname, coltype, pk, idx in zip(colnames, coltypes, pks, idxs))
                  , extend_existing=True)

    return table

def dict_to_engine(params, **kwargs):
    """
    create a database connection based on the yaml description
    :param dict: the output section of the config yaml
    :param kwargs: if type is not sqlite in this order username, password, host, port
    :return: a sqlalchemy engine
    """
    if params["type"]=="sqlite":
        dbstring="sqlite:///{}".format(params["name"])
    else:
        raise NotImplementedError("currently only sqlite is supported")
        #dbstring = "{}://{}:{}@{}:{}/{}".format(kwargs, params["name"])

    engine=create_engine(dbstring)
    return engine