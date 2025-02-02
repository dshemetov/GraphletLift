import pickle
import networkx as nx
import numpy as np
import lift_np as lt

NUM_STEPS = 200000

def run_longtest(graph_name, k):
    lift = lt.Lift(graph_name, k)
    lift.get_graphlet_count(num_steps=NUM_STEPS)
    with open("experiments/" + graph_name
              + "_" + str(k) + '_samples.pickle', 'wb') as f:
        pickle.dump(lift.graphlet_samples, f)
    with open("experiments/" + graph_name
              + "_" + str(k) + '_graphletcounts.pickle', 'wb') as f:
        pickle.dump(lift.graphlet_counts, f)

graph_names = [
    "bio-celegansneural",
    "ia-email-univ",
    "misc-fullb",
    "as-caida",
    "socfb-B-anon"
]
for graph_name in graph_names:
    run_longtest(graph_name, 3)
    run_longtest(graph_name, 4)
    run_longtest(graph_name, 5)
    run_longtest(graph_name, 6)
