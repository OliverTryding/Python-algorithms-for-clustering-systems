import networkx as nx
from matplotlib import pyplot as plt
from networkx.classes.function import neighbors
#from itertools import izip

def get_leaves(clusterSys):
    X = []
    i = len(clusterSys) - 1
    while i >= 0:
        if len(clusterSys[i]) == 1:
            X.append(clusterSys[i])
        else:
            break
        i -= 1
    return(X)

def set_comp(A,B):
    """Compares two sets.

    First priority in the comparison is the number of elements in the set; 
    A set with more elements should come first in the sorted clustering set. 
    If two sets are of equal length then 
    the elements of each array are pairwise compared first to last 
    and smaller elements are prioritised 
    with the corresponding set placed first.

    Paramters
    -----------------
    A,B : sets

    Returns
    -----------------
    Boolean
        True if A would be prioritized over B.
        False otherwise.
    """

    if len(A) == len(B): #Longer array is prioritised
        for (a,b) in zip(A,B):
            if a < b: #Smaller numbers come first.
                return(True)
            elif a > b:
                return(False)
    elif len(A) > len(B):
        return(True)

    return(False)

def clusterSysSort(clusterSys):
    """A sorting algorithm for clustering systems.
    Works similarly to mergesort.
    
    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.

    Returns
    -----------------
    Null
    """

    if len(clusterSys) > 1:

        m = len(clusterSys) // 2
        L = clusterSys[:m]
        R = clusterSys[m:]
        clusterSysSort(L)
        clusterSysSort(R)

        i = j = k = 0
        while i < len(L) and j < len(R):
            if set_comp(L[i],R[j]):
                clusterSys[k] = L[i]
                i += 1
            else:
                clusterSys[k] = R[j]
                j += 1
            k += 1

        while i < len(L):
            clusterSys[k] = L[i]
            i += 1
            k += 1
  
        while j < len(R):
            clusterSys[k] = R[j]
            j += 1
            k += 1

def includes(clusterSys,A):
    """includes returns true if and only if set A is in the clustering system.

    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.
    A : a numeric set.

    Returns
    -----------------
    Boolean
        True if A is in clusterSys.
        False otherwise.
    """

    m = len(clusterSys) // 2
    M = clusterSys[m]
    L = clusterSys[:m]
    R = clusterSys[m:]

    if A == M:
        return(True)
    elif len(clusterSys) > 1:
        if set_comp(A,M):
            return includes(L,A)
        else:
            return includes(R,A)
    return(False)

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

def weak_hierarchy_help(clusterSys,intersections,i,j,k):
    if frozenset((i,j)) in intersections.keys():
        intersection_1 = intersections.get(frozenset((i,j)))
    else:
        intersection_1 = clusterSys[i].intersection(clusterSys[j])
        intersections[frozenset((i,j))] = intersection_1

    if frozenset((i,k)) in intersections.keys():
        intersection_2 = intersections.get(frozenset((i,k)))
    else:
        intersection_2 = clusterSys[i].intersection(clusterSys[k])
        intersections[frozenset((i,k))] = intersection_2

    if frozenset((j,k)) in intersections.keys():
        intersection_3 = intersections.get(frozenset((j,k)))
    else:
        intersection_3 = clusterSys[j].intersection(clusterSys[k])
        intersections[frozenset((j,k))] = intersection_3

    intersection = intersection_1.intersection(clusterSys[k])
    if not intersection or intersection in [intersection_1,intersection_2,intersection_3]:
        return(True)
    else:
        return(False)

def test_weak_hierarchy(clusterSys):
    """Returns True if the network corresponding to 
    clusterSys fulfills the property 'weak hierarchy'.
    
    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.

    Returns
    -----------------
    Boolean
        True if clusterSys is a weak hierarchy.
        False otherwise.
    """

    intersections = {}

    i = 0
    while i < len(clusterSys) - 2:
        j = i + 1
        while j < len(clusterSys) - 1:
            k = j + 1
            while k < len(clusterSys):
                if weak_hierarchy_help(clusterSys,intersections,i,j,k):
                    pass
                else:
                    return(False)
                k += 1
            j += 1
        i += 1
    return(True)

def closed_help(clusterSys,intersections,i,j):
    if frozenset((i,j)) in intersections.keys():
        intersection = intersections.get(frozenset((i,j)))
    else:
        intersection = clusterSys[i].intersection(clusterSys[j])
        intersections[frozenset((i,j))] = intersection

    if not intersection or includes(clusterSys,intersection):
        return(True)
    else:
        return(False)

def test_closed(clusterSys):
    """Returns True if clusterSys fulfills the 'closed' property.

    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.

    Returns
    -----------------
    Boolean
        True if clusterSys is closed.
        False otherwise.
    """

    i = 0
    while i < len(clusterSys)-1:
        j = i + 1
        while j < len(clusterSys):
            #---------------closed-----------------------
            intersection = clusterSys[i].intersection(clusterSys[j])

            if not intersection or includes(clusterSys,intersection):
                pass
            else:
                return(False)
            #---------------closed-----------------------
            j += 1
        i += 1
    return(True)

def property_L_help(clusterSys,intersections,i,j,k):
    if frozenset((i,j)) in intersections.keys():
        intersection_1 = intersections.get(frozenset((i,j)))
    else:
        intersection_1 = clusterSys[i].intersection(clusterSys[j])
        intersections[frozenset((i,j))] = intersection_1

    if frozenset((i,k)) in intersections.keys():
        intersection_2 = intersections.get(frozenset((i,k)))
    else:
        intersection_2 = clusterSys[i].intersection(clusterSys[k])
        intersections[frozenset((i,k))] = intersection_2

    if intersection_1 == intersection_2:
        return(True)
    elif overlaps(clusterSys[i],clusterSys[j]) and overlaps(clusterSys[i],clusterSys[k]):
        return(False)
    else:
        return(True)

def test_property_L(clusterSys):
    """Returns True if clusterSys fulfills the 'L' property.

    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.

    Returns
    -----------------
    Boolean
        True if clusterSys fulfills 'property L'.
        False otherwise.
    """

    intersections = {}

    i = 0
    while i < len(clusterSys):
        j = 0
        while j < len(clusterSys) - 1:
            k = j + 1
            while k < len(clusterSys):
                if property_L_help(clusterSys,intersections,i,j,k):
                    pass
                else:
                    return(False)
                k += 1
            j += 1
        i += 1
    return(True)

def is_level_1(clusterSys):
    """Returns True if clusterSys fulfills the 'closed' and the 'L' properties 
    and therefore are compatible with some level-1 network.

    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.

    Returns
    -----------------
    Boolean
        True if clusterSys fulfills 'closed' and 'property L'.
        False otherwise.
    """

    intersections = {}

    i = 0
    while i < len(clusterSys):
        j = 0
        while j < len(clusterSys)-1:

            if i < j:
                if closed_help(clusterSys,intersections,i,j):
                    pass
                else:
                    return(False)
            k = j
            while k < len(clusterSys):
                if property_L_help(clusterSys,intersections,i,j,k):
                    pass
                else:
                    return(False)
                k += 1
            j += 1
        i += 1
    return(True)

def clusterSys_to_overlap(clusterSys):
    """Returns a overlap graph from a clustering system.
    
    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.

    Returns
    -----------------
    NetworkX graph
    """

    size = len(clusterSys)

    g = nx.Graph()

    i = 0
    while i < size-1:
        j = i + 1
        while j < size:
            if overlaps(clusterSys[i],clusterSys[j]):
                g.add_edge(i,j)
            j += 1
        i += 1

    return(g)

def test_N3O(G):
    """Returns true if the overlap graph does not contain any triangles, 
    making the corresponding clustering system fulfill the 'N3O' property.
    
    Paramters
    -----------------
    G : a networkx graph where edges between vertices represent an overlap of the vertices.

    Returns
    -----------------
    Boolean
        True if G contains no triangles, 
        meaning that the corresponding clustering system 
        has three pairwise overlaping clusters.

        False otherwise.
    """
    E = G.edges()
    for e in E:
        if not set(G[e[0]]).isdisjoint(set(G[e[1]])):
            return(False)
    return(True)

def test_hierarchy(G):
    """Returns true if the overlap graph does not contain any edges, 
    making the corresponding clustering system fulfill the 'hierarchy' property.
    
    Paramters
    -----------------
    G : a networkx graph where edges between vertices represent an overlap of the vertices.

    Returns
    -----------------
    Boolean
        True if G contains no edges, 
        meaning that the corresponding clustering system has no overlaps.

        False otherwise.
    """
    if not G.edges:
        return(True)
    else:
        return(False)

def test_paired_hierarchy(G):
    """Returns true if no vertices in the overlap graph 
    has a degree of larger than 1. 
    Meaning no cluster in the corresponding clustering system 
    overlaps with more than one other cluster.
    
    Paramters
    -----------------
    G : a networkx graph where edges between vertices represent an overlap of the vertices.

    Returns
    -----------------
    Boolean
        True if every vertex in G as a degree of at most 1.
        False otherwise.
    """

    V = G.nodes()
    for v in V:
        if G.degree[v] > 1:
            return(False)
    return(True)

def inc_Min(clusterSys, pair):
    """Returns the unique inclusion minimal set 
    in clusterSys with respect to a pair or 
    returns an empty set if there is no such unique inclusion minimal set.
    
    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.

    Returns
    -----------------
    set : if possible, a unique inclusion minimal set 
    in clusterSys for the pair. Otherwise an empty set.
    """

    k = len(clusterSys) - 1
    out = set()
    while k >= 0:
        C = clusterSys[k]

        if pair[0] in C and pair[1] in C:
            if not out:
                out = C
            elif out < C:
                pass
            else:
                return({})

        k -= 1
    return(out)

def test_prebinary(clusterSys, X):
    """Returns True if clusterSys fulfills the 'prebinary' property.

    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.
    X          : a list of numeric sets representing the leaves.

    Returns
    -----------------
    Boolean
        True if clusterSys fulfills 'prebinary'.
        False otherwise.
    """

    i = 0
    while i < len(X) - 1:

        x = next(iter(X[i]))

        j = i + 1
        while j < len(X):

            y = next(iter(X[j]))
            
            C = inc_Min(clusterSys, (x,y))
            if not C:
                return(False)

            j += 1
        i += 1
    return(True)

def test_binary(clusterSys, X):
    """Returns True if clusterSys fulfills the 'binary' property.

    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.
    X          : a list of numeric sets representing the leaves.

    Returns
    -----------------
    Boolean
        True if clusterSys fulfills 'binary'.
        False otherwise.
    """

    CS = set()

    i = 0
    while i < len(X) - 1:

        x = next(iter(X[i]))

        j = i + 1
        while j < len(X):

            y = next(iter(X[j]))
            
            C = inc_Min(clusterSys, (x,y))
            if not C:
                return(False)
            else:
                CS.add(frozenset(C))

            j += 1
        i += 1
    if len(CS) == len(clusterSys) - len(X):
        return(True)
    else:
        return(False)

def clusterSys_to_Hasse(clusterSys):
    """Returns a networkX graph representing the corresponding 
    Hasse diagram to the input clustering system.
    
    Paramters
    -----------------
    clusterSys : a list of numeric sets representing a clustering system.

    Returns
    -----------------
    NetworkX graph
    """

    size = len(clusterSys)

    g = nx.DiGraph()

    i = 0
    while i < size-1:
        j = i + 1
        while j < size:
            if clusterSys[i] > clusterSys[j]:
                g.add_edge(i,j)
            j += 1
        i += 1
    tg = nx.transitive_reduction(g)

    return(tg)

def test_2Inc(N):
    """Returns True if the corresponding clustering system 
    fulfills the '2-Inc' property.

    Paramters
    -----------------
    N : A Hasse diagram corresponding to the clustering system.

    Returns
    -----------------
    Boolean
        True if the corresponding clustering system fulfills '2-Inc'.
        False otherwise.
    """

    V = N.nodes()

    i = 0
    while i < len(V):
        if N.out_degree[i] <= 2 and N.in_degree[i] <= 2:
            pass
        else:
            return(False)
        i += 1

    return(True)

def draw_network(N):
    V = N.nodes()
    size = len(V)

    if N == None:
        print("No such graph!")
    #pos = nx.planar_layout(N)
    labels = {0: "r"}
    i = j = 1
    while i < size:
        if N.out_degree[i] != 0:
            labels[i] = ""
        else:
            labels[i] = j
            j += 1
        i += 1
    nx.draw_planar(N, labels = labels)
    plt.show()
    #layout = g.layout("kk")
    #plot(g,vertex_label_size = 25, vertex_label_dist=2, bbox = (1000, 1000), margin = 100, layout = layout)
