"""
Unit tests for the LiftSRW module.
"""
import networkx as nx
import lift as lt

NUM_STEPS = 20000
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
    elif graph == "bio-celegansneural":
        graphlet_counts = run_test("bio-celegansneural", 3)
        actual_triangle_count = 12.6 * 10**3 / 3
        assert ((graphlet_counts["triangle"] - actual_triangle_count)
                / actual_triangle_count < THRESHOLD)
        print(graph + " passed.")
        print(graphlet_counts, "\n")
    elif graph == "ia-email-univ":
        graphlet_counts = run_test("bio-celegansneural", 3)
        actual_triangle_count = 16000 / 3
        assert ((graphlet_counts["triangle"] - actual_triangle_count)
                / actual_triangle_count < THRESHOLD)
        print(graph + " passed.")
        print(graphlet_counts, "\n")
    elif graph == "misc-fullb":
        graphlet_counts = run_test("bio-celegansneural", 3)
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

test_graph("path")
test_graph("wheel")
test_graph("ladder")
test_graph("bio-celegansneural")
test_graph("ia-email-univ")
test_graph("misc-fullb")
test_graph("misc-polblogs")

# # ICYMI: nx.star_graph(4) has 5 nodes.
# graphlet_counts = run_test(nx.star_graph(4), 5)
# assert graphlet_counts[0] == 1

graphlet_counts = run_test(nx.complete_graph(5), 5)
assert graphlet_counts[20] == 1

graphlet_counts = run_test(nx.complete_graph(10), 5)
assert graphlet_counts[20] == 252
