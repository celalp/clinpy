import sqlalchemy.sql.operators as ops
from sqlalchemy import or_, and_
from sqlalchemy.sql import null


def str_to_op(st):
    if st == "==":
        op = ops.eq
    elif st == ">":
        op = ops.gt
    elif st == ">=":
        op = ops.ge
    elif st == "<":
        op = ops.lt
    elif st == "<=":
        op = ops.le
    elif st == "!=":
        op = ops.ne
    elif st == "like":
        op = ops.like_op
    elif st == "ilike":
        op = ops.ilike_op
    elif st == "notlike":
        op = ops.notlike_op
    elif st == "notilike":
        op = ops.not_ilike_op
    elif st == "in":
        op = ops.in_op
    elif st == "not_in":
        op = ops.not_in_op
    else:
        raise NotImplementedError("{} is not an available comparison".format(st))
    return op


def str_to_null(st):
    if st == "null":
        value = null
    else:
        raise ValueError("this is not a null")
    return value

def single_filter(tup, meta):
    """
    :param tup a tuple in this order (tablename, column name, comparison see str_to_op, and value)
    :param meta: database connection metadata
    """
    table = meta.tables[tup[0]]
    col = table.c[tup[1]]
    op = str_to_op(tup[2])
    if tup[2]=="in" or tup[2]=="not_in":
        if type(tup[3]) != list:
            raise TypeError("for in and not_in operations a list needs to be"
                            "provided")
    if tup[3]=="null":
        value=str_to_null(tup[3])
    else:
        value=tup[3]
    filter = op(col, value)
    return filter

def compound_filter(logic, filter_list, meta):
    """
    #TODO
    """
    filters = []
    for item in filter_list:
        if type(item) == tuple:
            filt = single_filter(item, meta)
            filters.append(filt)
        elif type(item) == dict:
            for key in item.keys():
                filt = compound_filter(key, item[key], meta)
                filters.append(filt)
        else:
            raise ValueError("Please check your filters, the types should only be tuple and dict")

    if logic == "and":
        compounded = and_(*filters).self_group()
    elif logic == "or":
        compounded = or_(*filters).self_group()
    else:
        raise ValueError("only 'and' or 'or' is supported")

    return compounded


def apply_filter(query, filters):
    """
    take the output of compound filters and apply them to the query one by one
    :param query: this is a sqlalchemy select, it can be anything, any number of joins, pre applied filters etc are fine
    you can take the output of this function and run it again with some different filters as well
    :param filters: the output of compound filters
    :return another select with the filters applied using sqlalchemy.filter()
    """
    for filt in filters:
        query = query.filter(filt)

    return query
