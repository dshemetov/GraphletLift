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

def graphlet_list(k):
    assert k > 0
    foo = 1
    loc_graphlet_list = {n: [] for n in range(1,k+1)}
    while True:
        G = nx.graph_atlas(foo)
        n = G.number_of_nodes()
        if n>k:
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
              for (i,T_) in enumerate(cached_graphlet_list[n])
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
              for (i,T_) in enumerate(cached_graphlet_list[n])
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

def random_walk_nodes(G, v0, steps_num):
    curr_vert = v0
    for _ in range(steps_num):
        curr_vert = random.choice(list(G.neighbors(curr_vert)))
    return curr_vert

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
            'Graphs/polblogs.edgelist',
            create_using = nx.Graph())

    if name=='misc-as-caida':
        G = nx.read_edgelist(
            'Graphs/as-caida.edgelist',
            create_using = nx.Graph())

    if name=='misc-fullb':
        G = nx.read_edgelist(
            'Graphs/fullb.edgelist',
            create_using = nx.Graph())

    if name == 'socfb-B-anon':
        G = nx.read_edgelist(
            'Graphs/socfb-B-anon.edgelist',
            create_using=nx.Graph())

    if G is None:
        raise KeyError

    return G

def waddling_mixing_variance(G, k, steps_num=1000, burn_in_limit=20):
    if k==4:
        return waddling_mixing_variance_4(G, steps_num=steps_num, burn_in_limit=burn_in_limit)
    if k==3:
        return waddling_mixing_variance_3(G, steps_num=steps_num, burn_in_limit=burn_in_limit)

def waddling_mixing_variance_3(G, steps_num=1000, burn_in_limit=20):
    k=3
    longest_paths = {0:2, 1:6}
    v0 = random.choice(list(G.nodes()))
    graphlet_num = len(cached_graphlet_list[k])
    exp_counter = {i:0 for i in range(graphlet_num)}
    type_counter = {i:0 for i in range(graphlet_num)}
    var_counter = {i:0 for i in range(graphlet_num)}
    corr_counter = {i:
                    {burn_in: 0
                     for burn_in in range(0,burn_in_limit)}
                    for i in range(graphlet_num)}
    variance = {i:0 for i in range(graphlet_num)}

    memory = [None for _ in range(burn_in_limit)]
    for _ in range(steps_num):
        v1 = random_walk_nodes(G, v0, 1)
        v2 = random_walk_nodes(G, v1, 1)
        T = {v0, v1, v2}
        if len(T)==k:
            T_type = find_type(subgraph(G,T))
            T_prob = (longest_paths[T_type] *
                      G.degree(v1)**(-1) *
                      (2*cached_edge_number)**(-1))
            type_counter[T_type] += 1
            exp_counter[T_type] += (T_prob)**(-1)
            var_counter[T_type] += (T_prob)**(-2)
            ind = 0
            while ind < burn_in_limit and memory[ind] is not None:
                S_type, S_prob = memory[ind]
                if T_type==S_type:
                    #pair_counter[T_type][ind] += 1
                    corr_counter[T_type][ind] += 2*(T_prob*S_prob)**(-1)
                ind+=1
            memory = [(T_type, T_prob)] + memory[:-1]
        else:
            memory = [(-1, 0)] + memory[:-1]
        v0 = random_walk_nodes(G, v2, 1)

    beta_coeff = {i: [abs(corr_counter[i][burn_in]
                          *(steps_num - burn_in)
                          *exp_counter[i]**(-2) - 1)
                      for burn_in in range(burn_in_limit)]
                  for i in range(graphlet_num)
                  if type_counter[i]!=0}

    for i in range(graphlet_num):
        if exp_counter[i]!=0:
            variance[i] = (var_counter[i]*(steps_num)
                           *(exp_counter[i])**(-2)-1)

    print("Expectation")
    for i in range(graphlet_num):
        print(exp_counter[i]*(steps_num)**(-1))

    print("Variance")
    for i in range(graphlet_num):
        print(variance[i])

    for i in range(graphlet_num):
        print ("Graphlet ID{}".format(i))
        if exp_counter[i]!=0:
            for burn_in, val in enumerate(beta_coeff[i]):
                print("({0}, {1:.5f})".format(burn_in+1, val))
    return (beta_coeff, variance)

def waddling_mixing_variance_4(G, steps_num=1000, burn_in_limit=20):
    k=4
    longest_paths = {0:6, 1:2, 2:4, 3:8, 4:12, 5:24}
    v0 = random.choice(list(G.nodes()))
    graphlet_num = len(cached_graphlet_list[k])
    exp_counter = {i:0 for i in range(graphlet_num)}
    type_counter = {i:0 for i in range(graphlet_num)}
    var_counter = {i:0 for i in range(graphlet_num)}
    corr_counter = {i:
                    {burn_in: 0
                     for burn_in in range(0,burn_in_limit)}
                    for i in range(graphlet_num)}
    variance = {i:0 for i in range(graphlet_num)}
    expectation = {i:0 for i in range(graphlet_num)}

    memory = [None for _ in range(burn_in_limit)]
    for _ in range(steps_num/2):
        v1 = random_walk_nodes(G, v0, 1)
        v2 = random_walk_nodes(G, v1, 1)
        v3 = random_walk_nodes(G, v2, 1)
        T = {v0, v1, v2, v3}
        if len(T)==k:
            T_type = find_type(subgraph(G,T))
            T_prob = (longest_paths[T_type] *
                      (G.degree(v1)*G.degree(v2))**(-1) *
                      (2*cached_edge_number)**(-1))
            type_counter[T_type] += 1
            exp_counter[T_type] += (T_prob)**(-1)
            var_counter[T_type] += (T_prob)**(-2)
            ind = 0
            while ind < burn_in_limit and memory[ind] is not None:
                S_type, S_prob = memory[ind]
                if T_type==S_type:
                    #pair_counter[T_type][ind] += 1
                    corr_counter[T_type][ind] += (T_prob*S_prob)**(-1)
                ind+=1
            memory = [(T_type, T_prob)] + memory[:-1]
        else:
            memory = [(-1, 0)] + memory[:-1]
        v0 = random_walk_nodes(G, v3, 1)

    memory = [None for _ in range(burn_in_limit)]
    for _ in range(steps_num/2):
        v1 = random_walk_nodes(G, v0, 1)
        v2 = random_walk_nodes(G, v1, 1)
        v3 = random.choice(list(G.neighbors(v1)))
        T = {v0,v1,v2,v3}
        if len(T)==4 and find_type(subgraph(G,T))==0:
            T_prob = (longest_paths[0] *
                      (G.degree(v1))**(-2) *
                      (2*cached_edge_number)**(-1))
            type_counter[0] += 1
            exp_counter[0] += (T_prob)**(-1)
            var_counter[0] += (T_prob)**(-2)
            ind = 0
            while ind < burn_in_limit and memory[ind] is not None:
                S_type, S_prob = memory[ind]
                if S_type==0:
                    #pair_counter[0][ind] += 1
                    corr_counter[0][ind] += (T_prob*S_prob)**(-1)
                ind+=1
            memory = [(0, T_prob)] + memory[:-1]
        else:
            memory = [(-1, 0)] + memory[:-1]
        v0 = random_walk_nodes(G, v2, 1)

    for i in range(graphlet_num):
        expectation[i] = exp_counter[i]*(steps_num/2)**(-1)
        variance[i] = (var_counter[i]*(steps_num/2)**(-1)
                       - expectation[i]**2)

    correlation = {i: [(corr_counter[i][burn_in]*(steps_num/2-burn_in)**(-1)- expectation[i]**2)
                       *(variance[i])**(-1)
                       for burn_in in range(burn_in_limit)]
                   for i in range(graphlet_num)
                   if variance[i]!=0}

    print("Expectation")
    for i in range(graphlet_num):
        print(expectation[i])

    print("Normalized Variance")
    for i in range(graphlet_num):
        if expectation[i]!=0:
            print(variance[i]*expectation[i]**(-2))
        else:
            print("No graphlets found")

    for i in range(graphlet_num):
        print ("Correlation for Graphlet ID{}".format(i+1))
        if expectation[i]!=0:
            for burn_in, val in enumerate(correlation[i]):
                print("({0}, {1:.5f})".format(burn_in+1, val))
        else:
            print("No graphlets found")
    return (correlation, variance)
#
# def waddling_count(G, k, steps_num, burn_in, ground_truth):
#     if k==4:
#         return waddling_count_4(G, steps_num=steps_num, burn_in=burn_in, ground_truth=ground_truth)
#     if k==3:
#         return waddling_count_3(G, steps_num=steps_num, burn_in=burn_in, ground_truth=ground_truth)

def waddling_count_3(G, steps_num, burn_in):
    k=3
    longest_paths = {0:2, 1:6}
    v0 = random.choice(list(G.nodes()))
    graphlet_num = len(cached_graphlet_list[k])
    exp_counter = {i:0 for i in range(graphlet_num)}
    samples = []
    #errors = {i:[] for i in range(graphlet_num)}
    for step in range(1, steps_num+1):
        v1 = random_walk_nodes(G, v0, 1)
        v2 = random_walk_nodes(G, v1, 1)
        T = {v0, v1, v2}
        if len(T)==k:
            if v0 in G.neighbors(v2):
                T_type = 1
            else:
                T_type = 0
            T_prob = (longest_paths[T_type] *
                      G.degree(v1)**(-1) *
                      (2*cached_edge_number)**(-1))
            exp_counter[T_type] += (T_prob)**(-1)
            samples.append([T,T_type,T_prob])
        v0 = random_walk_nodes(G, v2, burn_in)

    return exp_counter, samples

def waddling_count_4(G, steps_num, burn_in):
    # assert k==4
    v0 = random.choice(list(G.nodes()))
    graphlet_num = len(cached_graphlet_list[4])
    exp_counter = {i:0 for i in range(graphlet_num)}
    longest_paths = {0:6, 1:2, 2:4, 3:8, 4:12, 5:24}
    samples = []
    for _ in range(int(steps_num/2)):
        v1 = random_walk_nodes(G, v0, 1)
        v2 = random_walk_nodes(G, v1, 1)
        v3 = random_walk_nodes(G, v2, 1)
        T = {v0, v1, v2, v3}
        if len(T)==4:
            T_type = find_type(subgraph(G,T))
            assert T_type != 0
            T_prob = (longest_paths[T_type] *
                      (G.degree(v1)*G.degree(v2))**(-1) *
                      (2*cached_edge_number)**(-1))
            exp_counter[T_type] += (T_prob)**(-1)
            samples.append([T, T_type, T_prob])
        v0 = random_walk_nodes(G, v3, burn_in)

    for _ in range(int(steps_num/2)):
        v1 = random_walk_nodes(G, v0, 1)
        v2 = random_walk_nodes(G, v1, 1)
        v3 = random.choice(list(G.neighbors(v1)))
        T = {v0,v1,v2,v3}
        if len(T)==4 and find_type(subgraph(G,T))==0:
            T_prob = (longest_paths[0] *
                      (G.degree(v1))**(-2) *
                      (2*cached_edge_number)**(-1))
            exp_counter[0] += (T_prob)**(-1)
            samples.append([T, T_type, T_prob])
        v0 = random_walk_nodes(G, v2, burn_in)

    expectation = {i: exp_counter[i] * (steps_num/2)**(-1)
                   for i in range(graphlet_num)}
    return expectation, samples
#
# def waddling_count_5(graph, steps_num, burn_in):
#     v0 = random.choice(list(graph.nodes()))
#     graphlet_num = len(cached_graphlet_list[k])
#     exp_counter = {i:0 for i in range(graphlet_num)}
#     longest_paths = {i:get_coefficient(g) for i, g in enumerate(gl[k])}
#
#     for _ in range(steps_num):
#         path = [v0]
#         for i in range(5):
#             path.append(random_walk_nodes(graph, v0, 1))
#         unique_nodes = len(set(path))
#         if unique_nodes == 3:
#             neighbors = get_neighbors_path(path, graph)
#             choices = [random.choice(neighbors) for i in range(2)]
#             if len(set(choices)) == 2:
#                 T_type = 0
#                 T_prob = ((G.degree(path[1])**3) *
#                           (2*cached_edge_number) /
#                           longest_paths[T_type])
#                 exp_counter[T_type] += T_prob
#
#         if unique_nodes == 4:
#             neighbors = get_neighbors_path(path, graph)
#             choices = [random.choice(neighbors) for i in range()]
#             if len(set(choices)) == 2:
#                 T_type = 0
#                 T_prob = ((G.degree(path[1])**3) *
#                           (2*cached_edge_number) /
#                           longest_paths[T_type])
#                 exp_counter[T_type] += T_prob
#
#         v0 = random_walk_nodes(G, path[-1], burn_in)

import pickle

graph_names = [
    "bio-celegansneural",
    "ia-email-univ",
    "misc-fullb",
    "misc-polblogs",
    "misc-as-caida",
    "socfb-B-anon"
]
cached_graphlet_list = graphlet_list(6)

for graph_name in graph_names:
    G = load_graph(graph_name)
    cached_edge_number = G.number_of_edges()
    NUM_STEPS = 10**7
    BURN_IN = 5
    expectation, samples = waddling_count_3(G, NUM_STEPS, BURN_IN)
    with open("experiments/waddle/" + graph_name
              + "_" + str(3) + '_expectation.pickle', 'wb') as f:
        pickle.dump(expectation, f)
    with open("experiments/waddle/" + graph_name
              + "_" + str(3) + '_samples.pickle', 'wb') as f:
        pickle.dump(samples, f)

    expectation, samples = waddling_count_4(G, NUM_STEPS, BURN_IN)
    with open("experiments/waddle/" + graph_name
              + "_" + str(4) + '_expectation.pickle', 'wb') as f:
        pickle.dump(expectation, f)
    with open("experiments/waddle/" + graph_name
              + "_" + str(4) + '_samples.pickle', 'wb') as f:
        pickle.dump(samples, f)
