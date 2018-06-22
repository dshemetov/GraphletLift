"""
Unit tests for the LiftSRW module.
"""
import lift as lt


def get_percent_error(vertex1, vertex2):
    """
    Obtain the percent error between two values.
    """
    return abs(vertex1 - vertex2) / float(vertex2)

def run_test(graph_name, k):
    """
    Simply runs the new code for a given graph and a given graphlet size.
    """
    lift_unordered = lt.Lift(graph_name, k, lift_type="unordered")
    graphlet_counts = lift_unordered.graphlet_count(steps_num=20000)
    graphlet_names = graphlet_counts.keys()
    for name in graphlet_names:
        print("Counts for " + str(name), graphlet_counts[name])

run_test("bio-celegansneural", 5)
run_test("bio-celegansneural", 1)
