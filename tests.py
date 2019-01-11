"""
Unit tests for the LiftSRW module.
"""
import networkx as nx
import lift as lt

STEPS_NUM = 20000
VERBOSE = True

def get_percent_error(vertex1, vertex2):
    """
    Obtain the percent error between two values.
    """
    return abs(vertex1 - vertex2) / float(vertex2)

def run_test(graph, k):
    """
    Simply runs the new code for a given graph and a given graphlet size.
    """
    lift_unordered = lt.Lift(graph, k, lift_type="unordered")
    graphlet_counts = lift_unordered.graphlet_count(steps_num=STEPS_NUM)
    #graphlet_names = graphlet_counts.keys()
    #for name in graphlet_names:
    #    print("Counts for " + str(name), graphlet_counts[name])
    return graphlet_counts

def test_path_graph():
    path_graph = nx.path_graph(5)
    lift_unordered = lt.Lift(path_graph, 3, lift_type="unordered")
    graphlet_counts = lift_unordered.graphlet_count(steps_num=STEPS_NUM)
    assert graphlet_counts["wedge"] == 3
    assert graphlet_counts["triangle"] == 0
    lift_unordered = lt.Lift(path_graph, 2, lift_type="unordered")
    graphlet_counts = lift_unordered.graphlet_count(steps_num=STEPS_NUM)
    assert graphlet_counts["2-path"] == 4

    print("Path graph passed.")

def test_wheel_graph():
    wheel_graph = nx.wheel_graph(6) # this is a 6-cycle with a star center node
    lift_unordered = lt.Lift(wheel_graph, 3, lift_type="unordered")
    graphlet_counts = lift_unordered.graphlet_count(steps_num=STEPS_NUM)
    assert graphlet_counts["wedge"] == 10
    assert graphlet_counts["triangle"] == 5
    lift_unordered = lt.Lift(wheel_graph, 2, lift_type="unordered")
    graphlet_counts = lift_unordered.graphlet_count(steps_num=STEPS_NUM)
    assert graphlet_counts["2-path"] == 10
    print("Wheel graph passed.")

def test_ladder_graph():
    ladder_graph = nx.ladder_graph(4) # this is two 6-paths joined one to one
    lift_unordered = lt.Lift(ladder_graph, 3, lift_type="unordered")
    graphlet_counts = lift_unordered.graphlet_count(steps_num=STEPS_NUM)
    assert graphlet_counts["wedge"] == 16
    assert graphlet_counts["triangle"] == 0
    lift_unordered = lt.Lift(ladder_graph, 2, lift_type="unordered")
    graphlet_counts = lift_unordered.graphlet_count(steps_num=STEPS_NUM)
    assert graphlet_counts["2-path"] == 10
    print("Ladder graph passed.")

def test_big_graph(name, triangle_count, verbose=False):
    graphlet_counts = run_test(name, 3)
    triangle_count = triangle_count / 3.0
    threshold = .1 * triangle_count
    if verbose:
        print("Given: ", triangle_count,
              "; empirical: ", graphlet_counts["triangle"])
    else:
        assert abs(graphlet_counts["triangle"] - triangle_count) < threshold
    print(name + " triangles passed.")


# TEST THE TESTS
# print("Path graph:")
# run_test(nx.path_graph(5),3)
# run_test(nx.path_graph(5),2)
# print("Wheel graph:")
# run_test(nx.wheel_graph(6),3)
# run_test(nx.wheel_graph(6),2)
# print("Ladder graph:")
# run_test(nx.ladder_graph(4),3)
# run_test(nx.ladder_graph(4),2)scrip


# SMALL GRAPH TESTS
test_path_graph()
test_wheel_graph()
test_ladder_graph()


# BIG GRAPH TESTS
# # http://networkrepository.com/bio-celegansneural.php
test_big_graph("bio-celegansneural", 9900, verbose=VERBOSE)
# # http://networkrepository.com/ia-email-univ.php
test_big_graph("ia-email-univ", 16000, verbose=VERBOSE)
# http://networkrepository.com/fullb.php
test_big_graph("misc-fullb", 180.6 * 10**6, verbose=VERBOSE)
#
# test_big_graph("misc-polblogs", 180.6 * 10**6, verbose=VERBOSE)
#
# test_big_graph("misc-as-caida", 180.6 * 10**6, verbose=VERBOSE)
