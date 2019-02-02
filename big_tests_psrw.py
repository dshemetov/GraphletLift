#import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.isomorphism as iso
import sympy
import numpy as np
import random
import time
import itertools
import math
#from IPython.display import clear_output

def graphlet_list(N):
    assert N > 0
    foo = 1
    loc_graphlet_list = {n: [] for n in range(1,N+1)}
    while True:
        G = nx.graph_atlas(foo)
        n = G.number_of_nodes()
        if n>N:
            break
        if nx.is_connected(G):
            loc_graphlet_list[n].append(G)
        foo += 1
    return loc_graphlet_list


def find_type_match(T):
    n = T.number_of_nodes()
    if n==1:
        return((0, {u: 0 for u in T.nodes()}))
    if n==2:
        return((0, {u: i for i,u in enumerate(T.nodes())}))
    if n==3:
        if T.number_of_edges()==2:
            u0 = next((node for node in T.nodes() if T.degree(node)==2))
            (u1,u2) = (node for node in T.neighbors(u0))
            return((0, {u0: 0, u1: 1, u2: 2}))
        if T.number_of_edges()==3:
            return((1,{u:i for i,u in enumerate(T.nodes())}))
    if n==4:
        e_num = T.number_of_edges()
        max_degree = max((T.degree(node) for node in T.nodes()))
        if e_num==3 and max_degree==3:
            u3 = next((node for node in T.nodes() if T.degree(node)==3))
            (u0,u1,u2) = (node for node in T.neighbors(u3))
            return((0, {u0:0, u1:1, u2:2, u3:3}))
        if e_num==3 and max_degree==2:
            (u0,u1) = (node for node in T.nodes() if T.degree(node)==2)
            u2 = next((node for node in T.neighbors(u1) if node!=u0))
            u3 = next((node for node in T.neighbors(u0) if node!=u1))
            return((1, {u0:0, u1:1, u2:2, u3:3}))
        if e_num==4 and max_degree==3:
            u3 = next((node for node in T.nodes() if T.degree(node)==3))
            (u1,u2) = (node for node in T.nodes() if T.degree(node)==2)
            u0 = next((node for node in T.nodes() if T.degree(node)==1))
            return((2, {u0:0, u1:1, u2:2, u3:3}))
        if e_num==4 and max_degree==2:
            u0 = next((node for node in T.nodes()))
            (u1,u3) = (node for node in T.neighbors(u0))
            u2 = next((node for node in T.neighbors(u1) if node!=u0))
            return((3, {u0:0, u1:1, u2:2, u3:3}))
        if e_num==5:
            (u0,u2) = (node for node in T.nodes() if T.degree(node)==3)
            (u1,u3) = (node for node in T.nodes() if T.degree(node)==2)
            return((4, {u0:0, u1:1, u2:2, u3:3}))
        if e_num==6:
            (u0,u1,u2,u3) = (node for node in T.nodes())
            return((5, {u0:0, u1:1, u2:2, u3:3}))
    # Improve matching procedure here for n>4.
    GM = next((i, iso.GraphMatcher(T,T_))
              for (i,T_) in enumerate(CACHED_GRAPHLET_LIST[n])
              if iso.GraphMatcher(T,T_).is_isomorphic())
    assert GM[1].is_isomorphic()
    return((GM[0],GM[1].mapping))

def find_type(T):
    n = T.number_of_nodes()
    if n==1:
        return 0
    if n==2:
        return 0
    if n==3:
        if T.number_of_edges()==2:
            return 0
        if T.number_of_edges()==3:
            return 1
    if n==4:
        e_num = T.number_of_edges()
        max_degree = max((T.degree(node) for node in T.nodes()))
        if e_num==3 and max_degree==3:
            return 0
        if e_num==3 and max_degree==2:
            return 1
        if e_num==4 and max_degree==3:
            return 2
        if e_num==4 and max_degree==2:
            return 3
        if e_num==5:
            return 4
        if e_num==6:
            return 5
    # Improve matching procedure here at least for n=4.
    GM = next((i
              for (i,T_) in enumerate(CACHED_GRAPHLET_LIST[n])
              if iso.GraphMatcher(T,T_).is_isomorphic()))
    return GM

def subgraph(G, nodes):
    list_nodes = list(nodes)
    T = nx.Graph()
    T.add_nodes_from(nodes)
    for i in range(len(nodes)):
        for j in range(i):
            if list_nodes[i] in G.neighbors(list_nodes[j]):
                T.add_edge(list_nodes[i],list_nodes[j])
    return T

def sub_edge_num(k, T_type):
    if k < 3:
        return k
    T = CACHED_GRAPHLET_LIST[k][T_type]
    count = 0
    for u in T.nodes():
        S = subgraph(T, T.nodes()-{u})
        if not nx.is_connected(S):
            continue
        for v in S.nodes():
            if not nx.is_connected(subgraph(S, S.nodes()-{u})):
                continue
            if not nx.is_connected(subgraph(T, T.nodes()-{v})):
                continue
            count+=1
    return count

def lift(G, vert, k):
    graphlet = set([vert])
    if k==1:
        return graphlet
    u = vert
    neig_list = []
    for n in range(2, k+1):
        neig_list = ([v for v in neig_list if v!=u]
                     + [v for v in G.neighbors(u) if v not in graphlet])
        u = random.choice(neig_list)
        graphlet.add(u)
    return graphlet

def SRW_step(G, graphlet):
    del_ins_list = []
    for u in graphlet:
        if not nx.is_connected(subgraph(G, graphlet-{u})):
            continue
        neigh = set()
        for v in graphlet-{u}:
            neigh.update(set(G.neighbors(v)))
        neigh = neigh - graphlet - {u}
        for w in neigh:
            del_ins_list.append((u,w))
    pair = random.choice(del_ins_list)
    new_graphlet = graphlet-{pair[0]}; new_graphlet.add(pair[1])
    return (new_graphlet, len(del_ins_list))

def load_graph(name, N=3):
    ground_truth = None
    G = None

    if name=='bio-celegansneural':
        G = nx.read_edgelist(
            'Graphs/bio-celegansneural.edgelist',
            create_using = nx.Graph())

    if name=='ia-email-univ':
        G = nx.read_edgelist(
            'Graphs/ia-email-univ.edgelist',
            create_using = nx.Graph())

    if name=='misc-polblogs':
        G = nx.read_edgelist(
            'Graphs/misc-polblogs.edgelist',
            create_using = nx.Graph())

    if name=='misc-as-caida':
        G = nx.read_edgelist(
            'Graphs/misc-as-caida.mtx',
            create_using = nx.Graph())

    if name=='misc-fullb':
        G = nx.read_edgelist(
            'Graphs/misc-fullb.mtx',
            create_using = nx.Graph())

    if graph_name == 'socfb-B-anon':
        graph = nx.read_edgelist(
            'Graphs/socfb-B-anon.edgelist',
            create_using=nx.Graph())

    if G is None:
        raise KeyError

    return {'graph': G, 'ground_truth': ground_truth}

def psrw_mixing_variance(G, k, steps_num=1000, burn_in_limit=20):
    v = random.choice(list(G.nodes()))
    init_graphlet = lift(G, v, k-1)
    old_graphlet = init_graphlet
    graphlet_num = len(CACHED_GRAPHLET_LIST[k])
    exp_counter = {i:0 for i in range(graphlet_num)}
    var_counter = {i:0 for i in range(graphlet_num)}
    type_counter = {i:0 for i in range(graphlet_num)}
    pair_counter = {i:
                    {burn_in: 0
                     for burn_in in range(0,burn_in_limit)}
                    for i in range(graphlet_num)}
    corr_counter = {i:
                    {burn_in: 0
                     for burn_in in range(0,burn_in_limit)}
                    for i in range(graphlet_num)}
    expectation = {i:0 for i in range(graphlet_num)}
    variance = {i:0 for i in range(graphlet_num)}
    cached_sub_edge_num = {T_type: sub_edge_num(k, T_type)
                           for T_type in range(len(CACHED_GRAPHLET_LIST[k]))
                          }

    memory = [None for _ in range(burn_in_limit)]
    for _ in range(steps_num):
        new_graphlet = SRW_step(G, old_graphlet)[0]
        T = old_graphlet.union(new_graphlet)
        old_graphlet = new_graphlet
        assert len(T)==k
        T_type = find_type(subgraph(G, T))
        T_prob = cached_sub_edge_num[T_type]
        type_counter[T_type] += 1
        exp_counter[T_type] += (T_prob)**(-1)
        var_counter[T_type] += (T_prob)**(-2)
        ind = 0
        while ind < burn_in_limit and memory[ind] is not None:
            S_type, S_prob = memory[ind]
            if T_type==S_type:
                pair_counter[T_type][ind] += 1
                corr_counter[T_type][ind] += (T_prob*S_prob)**(-1)
            ind+=1
        memory = [(T_type, T_prob)] + memory[:-1]

    for i in range(graphlet_num):
        expectation[i] = exp_counter[i]*steps_num**(-1)
        variance[i] = (var_counter[i]*steps_num**(-1)
                       - expectation[i]**2)

    correlation = {i: [(corr_counter[i][burn_in]*(steps_num-burn_in)**(-1)- expectation[i]**2)
                       *(variance[i])**(-1)
                       for burn_in in range(burn_in_limit)]
                   for i in range(graphlet_num)
                   if variance[i]!=0}

#     print("Expectation")
#     for i in range(graphlet_num):
#         print(expectation[i])

#     print("Normalized Variance")
#     for i in range(graphlet_num):
#         if expectation[i]!=0:
#             print(variance[i]*expectation[i]**(-2))
#         else:
#             print("No graphlets found")

#     for i in range(graphlet_num):
#         print ("Correlation for Graphlet ID{}".format(i+1))
#         if expectation[i]!=0:
#             for burn_in, val in enumerate(correlation[i]):
#                 print("({0}, {1:.5f})".format(burn_in+1, val))
#         else:
#             print("No graphlets found")
    return (expectation, variance, correlation)

def psrw_count(G, k, steps_num=1000, burn_in=10):
    v = random.choice(list(G.nodes()))
    init_graphlet = lift(G, v, k-1)
    old_graphlet = init_graphlet
    graphlet_num = len(CACHED_GRAPHLET_LIST[k])
    exp_counter = {i:0 for i in range(graphlet_num)}
    cached_sub_edge_num = {T_type: sub_edge_num(k, T_type)
                           for T_type in range(len(CACHED_GRAPHLET_LIST[k]))
                          }

    for _ in range(steps_num):
        new_graphlet = SRW_step(G, old_graphlet)[0]
        T = old_graphlet.union(new_graphlet)
        old_graphlet = new_graphlet
        assert len(T)==k
        T_type = find_type(subgraph(G, T))
        T_prob = cached_sub_edge_num[T_type]
        exp_counter[T_type] += (T_prob)**(-1)

    exp_counter = {i: exp_counter[i]*(steps_num)**(-1)
                   for i in range(graphlet_num)}

    return exp_counter

def run_psrw(graph_name, k, steps_num):
    G = load_graph(graph_name, k)['graph']
    expectation, variance, correlation = psrw_mixing_variance(G, k, steps_num)
    import pickle
    with open("experiments/psrw/" + graph_name
              + "_" + str(k) + '_expectation.pickle', 'wb') as f:
        pickle.dump(expectation, f)
    with open("experiments/psrw/" + graph_name
              + "_" + str(k) + '_variance.pickle', 'wb') as f:
        pickle.dump(variance, f)
    with open("experiments/psrw/" + graph_name
              + "_" + str(k) + '_correlation.pickle', 'wb') as f:
        pickle.dump(correlation, f)

CACHED_GRAPHLET_LIST = graphlet_list(6)
NUM_STEPS = 10**2

graph_names = [
    "bio-celegansneural",
    # "ia-email-univ",
    # "misc-fullb",
    # "as-caida",
]
for graph_name in graph_names:
    run_psrw(graph_name, 5, steps_num = NUM_STEPS)
    run_psrw(graph_name, 6, steps_num = NUM_STEPS)

# run_psrw("socfb-B-anon", 3, steps_num = 10**7)
# run_psrw("socfb-B-anon", 4, steps_num = 10**7)
# run_psrw("socfb-B-anon", 5, steps_num = 10**7)
# run_psrw("socfb-B-anon", 6, steps_num = 10**7)
