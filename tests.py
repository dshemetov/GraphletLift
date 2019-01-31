"""
Unit tests for the LiftSRW module.
"""
import networkx as nx
import numpy as np
import lift as lt

NUM_STEPS = 10000
VERBOSE = True
THRESHOLD = 0.1

def get_percent_error(v1, v2):
    """
    Obtain the percent error between two values.
    """
    return abs(v1 - v2) / float(v2)

def run_test(graph, k):
    """
    Simply runs the new code for a given graph and a given graphlet size.
    """
    lift_unordered = lt.Lift(graph, k, lift_type="unordered")
    graphlet_counts = lift_unordered.get_graphlet_count(
        num_steps=NUM_STEPS)

    return graphlet_counts

def test_graph(graph):
    """
    Tests that the sampling method obtains correct counts for the path graph,
    wheel graph, and ladder graph (see networkx documentation for graph
    details). Large scale test that gives a basic sanity check of the whole sampling class.
    """
    if graph == "path":
        path_graph = nx.path_graph(5)
        lift_unordered = lt.Lift(path_graph, 3, lift_type="unordered")
        graphlet_counts = lift_unordered.get_graphlet_count(
            num_steps=NUM_STEPS)
        assert graphlet_counts["wedge"] == 3
        assert graphlet_counts["triangle"] == 0
        lift_unordered = lt.Lift(path_graph, 2, lift_type="unordered")
        graphlet_counts = lift_unordered.get_graphlet_count(
            num_steps=NUM_STEPS)
        assert graphlet_counts["2-path"] == 4

        print("Path graph passed.")
    elif graph == "wheel":
        wheel_graph = nx.wheel_graph(6) # this is a 6-cycle with a star center node
        lift_unordered = lt.Lift(wheel_graph, 3, lift_type="unordered")
        graphlet_counts = lift_unordered.get_graphlet_count(
            num_steps=NUM_STEPS)
        assert graphlet_counts["wedge"] == 10
        assert graphlet_counts["triangle"] == 5
        lift_unordered = lt.Lift(wheel_graph, 2, lift_type="unordered")
        graphlet_counts = lift_unordered.get_graphlet_count(
            num_steps=NUM_STEPS)
        assert graphlet_counts["2-path"] == 10
        print("Wheel graph passed.")
    elif graph == "ladder":
        ladder_graph = nx.ladder_graph(4) # this is two 6-paths joined one to one
        lift_unordered = lt.Lift(ladder_graph, 3, lift_type="unordered")
        graphlet_counts = lift_unordered.get_graphlet_count(
            num_steps=NUM_STEPS)
        assert graphlet_counts["wedge"] == 16
        assert graphlet_counts["triangle"] == 0
        lift_unordered = lt.Lift(ladder_graph, 2, lift_type="unordered")
        graphlet_counts = lift_unordered.get_graphlet_count(
            num_steps=NUM_STEPS)
        assert graphlet_counts["2-path"] == 10
        print("Ladder graph passed.")
    elif graph == "bio-celegansneural":
        graphlet_counts = run_test("bio-celegansneural", 3)
        actual_triangle_count = 12.6 * 10**3 / 3
        assert ((graphlet_counts["triangle"] - actual_triangle_count)
                / actual_triangle_count < THRESHOLD)
        print(graph + " passed.")
        print(graphlet_counts, "\n")
    elif graph == "ia-email-univ":
        graphlet_counts = run_test("ia-email-univ", 3)
        actual_triangle_count = 16000 / 3
        assert ((graphlet_counts["triangle"] - actual_triangle_count)
                / actual_triangle_count < THRESHOLD)
        print(graph + " passed.")
        print(graphlet_counts, "\n")
    elif graph == "misc-fullb":
        graphlet_counts = run_test("misc-fullb", 3)
        actual_triangle_count = 180.6 * 10**6 / 3
        assert ((graphlet_counts["triangle"] - actual_triangle_count)
                / actual_triangle_count < THRESHOLD)
        print(graph + " passed.")
        print(graphlet_counts, "\n")
    elif graph == "misc-polblogs":
        graphlet_counts = run_test("misc-polblogs", 3)
        actual_triangle_count = 459.4 * 10**3 / 3
        assert ((graphlet_counts["triangle"] - actual_triangle_count)
                / actual_triangle_count < THRESHOLD)
        print(graph + " passed.")
        print(graphlet_counts, "\n")
    else:
        print("Graph unknown.")
#
# test_graph("path")
# test_graph("wheel")
# test_graph("ladder")
# test_graph("bio-celegansneural")
# test_graph("ia-email-univ")
# # test_graph("misc-fullb")
# # test_graph("misc-polblogs")
#
# # ICYMI: nx.star_graph(4) has 5 nodes.
# graphlet_counts = run_test(nx.star_graph(4), 5)
# assert graphlet_counts[0] == 1
# print("Star graph passed.")
#
# graphlet_counts = run_test(nx.complete_graph(5), 5)
# assert graphlet_counts[20] == 1
#
# graphlet_counts = run_test(nx.complete_graph(10), 5)
# assert graphlet_counts[20] == 252
# print("Complete graph passed.")

import time
lift = lt.Lift("bio-celegansneural", 4)
times = []
for i in range(100):
    start = time.time()
    lift.get_graphlet_count(num_steps=1)
    times.append(time.time() - start)
print(
    "Average time taken for a single iteration: ",
    sum(times)/100
    )

# # Pynauty tests.
# import networkx as nx
# import pynauty as na
#
# g = na.Graph(number_of_vertices=8, directed=False,
#           adjacency_dict = { 0: [1,2],
#                              1: [0,4],
#                              2: [0,3,4,5],
#                              3: [2,4,5,6],
#                              4: [1,2,3,5],
#                              5: [4,2,3,7],
#                              6: [3,7],
#                              7: [5,6] })
#
# h = na.Graph(number_of_vertices=8, directed=False,
#           adjacency_dict = { 0: [1,2],
#                              1: [0,4],
#                              2: [0,3,4,5],
#                              3: [2,4,5,7],
#                              4: [1,2,3,5],
#                              5: [4,2,3,6],
#                              6: [5,7],
#                              7: [3,6] })
#
# i = na.Graph(number_of_vertices=8, directed=False,
#           adjacency_dict = { 0: [2,4,3,6],
#                              1: [5,4],
#                              2: [0,3,4,5],
#                              3: [0,2,4,7],
#                              4: [0,1,2,3],
#                              5: [1,2],
#                              6: [0,7],
#                              7: [3,6] })
#
# j = na.Graph(number_of_vertices=8, directed=False,
#           adjacency_dict = { 0: [1,2],
#                              1: [0,3],
#                              2: [0,3,4,5],
#                              3: [1,2,4,5],
#                              4: [2,3,5,6],
#                              5: [2,3,4,5,6,7],
#                              6: [4,5,7],
#                              7: [5,6] })
#
# assert all([na.certificate(g) == na.certificate(h),
#            na.certificate(g) == na.certificate(i),
#            na.certificate(g) != na.certificate(j)])
#
# num_graphs = 4
# random_nx_graphs = [ nx.gnp_random_graph(8,0.4) for i in range(num_graphs) ]
# random_na_graphs = [ na.Graph(number_of_vertices = 8, directed = False,
#                               adjacency_dict = { n: list(nbrdict.keys()) for n, nbrdict in graph.adjacency() }
#                               ) for graph in random_nx_graphs
#                    ]
#
# nx_iso = [ nx.is_isomorphic(random_nx_graphs[i],random_nx_graphs[j]) for i in range(num_graphs) for j in range(i,num_graphs) ]
# na_iso = [ na.certificate(random_na_graphs[i]) == na.certificate(random_na_graphs[j]) for i in range(num_graphs) for j in range(i,num_graphs) ]
#
# assert nx_iso == na_iso
# print("Pynauty tests passed.")
# #print("What do the isomorphism matrices look like?\n", nx_iso, na_iso)
# #print("These are the na graphs:\n",random_na_graphs)
# #print("These are the nx graphs:\n",random_nx_graphs)
