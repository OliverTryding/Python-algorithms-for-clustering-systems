import networkx as nx
from matplotlib import pyplot as plt
from itertools import combinations

def overlaps(A,B):
    """returns True if sets A and B overlap.
    
    Paramters
    -----------------
    A,B : numeric sets.

    Returns
    -----------------
    Boolean
        True if A and B overlap.
        False otherwise.
    """
    if A.issubset(B):
        return(False)
    elif B.issubset(A):
        return(False)
    elif A.isdisjoint(B):
        return(False)
    else:
        return(True)

def closed_help(clusterSys,intersect1,intersect2,intersect3):
    if not intersect1 or (intersect1 in clusterSys):
        pass
    else:
        return(False)

    if not intersect2 or (intersect2 in clusterSys):
        pass
    else:
        return(False)

    if not intersect3 or (intersect3 in clusterSys):
        pass
    else:
        return(False)
    return(True)

def property_L_help(c,intersect1,intersect2,intersect3):
    if overlaps(c[0],c[1]) and overlaps(c[0],c[2]):
        if intersect1 == intersect2:
            pass
        else:
            return(False)

    if overlaps(c[1],c[0]) and overlaps(c[1],c[2]):
        if intersect1 == intersect3:
            pass
        else:
            return(False)

    if overlaps(c[2],c[0]) and overlaps(c[2],c[1]):
        if intersect2 == intersect3:
            pass
        else:
            return(False)
    return(True)

def is_level_1(clusterSys):
    """Returns True if clusterSys fullfills the 'closed' and the 'L' properties 
    and therefore are compatible with some level-1 network.

    Paramters
    -----------------
    clusterSys : a set of numeric sets representing a clustering system.

    Returns
    -----------------
    Boolean
        True if clusterSys fulfills 'closed' and 'property L'.
        False otherwise.
    """


    for c in combinations(clusterSys,3):
        intersect1 = c[0].intersection(c[1])
        intersect2 = c[0].intersection(c[2])
        intersect3 = c[1].intersection(c[2])
        #---------------closed-----------------------
        if closed_help(clusterSys,intersect1,intersect2,intersect3):
            pass
        else:
            return(False)
        #---------------closed-----------------------

        #---------------property L-------------------
        if property_L_help(c,intersect1,intersect2,intersect3):
            pass
        else:
            return(False)
        #---------------property L-------------------
    return(True)