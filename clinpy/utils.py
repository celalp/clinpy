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
