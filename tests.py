"""
Unit tests for the LiftSRW module.
"""
import networkx as nx
import lift as lt

NUM_STEPS = 13000
VERBOSE = True

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
    graphlet_counts = lift_unordered.graphlet_count(num_steps=NUM_STEPS)

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
        graphlet_counts = lift_unordered.graphlet_count(num_steps=NUM_STEPS)
        assert graphlet_counts["wedge"] == 3
        assert graphlet_counts["triangle"] == 0
        lift_unordered = lt.Lift(path_graph, 2, lift_type="unordered")
        graphlet_counts = lift_unordered.graphlet_count(num_steps=NUM_STEPS)
        assert graphlet_counts["2-path"] == 4

        print("Path graph passed.")
    elif graph == "wheel":
        wheel_graph = nx.wheel_graph(6) # this is a 6-cycle with a star center node
        lift_unordered = lt.Lift(wheel_graph, 3, lift_type="unordered")
        graphlet_counts = lift_unordered.graphlet_count(num_steps=NUM_STEPS)
        assert graphlet_counts["wedge"] == 10
        assert graphlet_counts["triangle"] == 5
        lift_unordered = lt.Lift(wheel_graph, 2, lift_type="unordered")
        graphlet_counts = lift_unordered.graphlet_count(num_steps=NUM_STEPS)
        assert graphlet_counts["2-path"] == 10
        print("Wheel graph passed.")
    elif graph == "ladder":
        ladder_graph = nx.ladder_graph(4) # this is two 6-paths joined one to one
        lift_unordered = lt.Lift(ladder_graph, 3, lift_type="unordered")
        graphlet_counts = lift_unordered.graphlet_count(num_steps=NUM_STEPS)
        assert graphlet_counts["wedge"] == 16
        assert graphlet_counts["triangle"] == 0
        lift_unordered = lt.Lift(ladder_graph, 2, lift_type="unordered")
        graphlet_counts = lift_unordered.graphlet_count(num_steps=NUM_STEPS)
        assert graphlet_counts["2-path"] == 10
        print("Ladder graph passed.")
    else:
        print "Graph unknown."

test_graph("path")
test_graph("wheel")
test_graph("ladder")


def test_graph_file(name, triangle_count, verbose=False):
    """
    Similar to the above, except it counts triangles on large graphs from files.
    """
    graphlet_counts = run_test(name, 3)
    triangle_count = triangle_count / 3.0
    threshold = .1 * triangle_count
    if verbose:
        print("Given: ", triangle_count,
              "; empirical: ", graphlet_counts["triangle"])
    else:
        assert abs(graphlet_counts["triangle"] - triangle_count) < threshold
    print(name + " triangles passed.")

# http://networkrepository.com/bio-celegansneural.php
test_graph_file("bio-celegansneural", 9900, verbose=VERBOSE)
# http://networkrepository.com/ia-email-univ.php
test_graph_file("ia-email-univ", 16000, verbose=VERBOSE)
# http://networkrepository.com/fullb.php
test_graph_file("misc-fullb", 180.6 * 10**6, verbose=VERBOSE)
# http://networkrepository.com/polblogs.php
test_graph_file("misc-polblogs", 459.4 * 10**3, verbose=VERBOSE)
