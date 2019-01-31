import pickle
import networkx as nx
import numpy as np
import lift as lt

def run_longtest(graph_name, k):
    lift = lt.Lift(graph_name, k)
    lift.get_graphlet_count(num_steps=200000)
    with open("experiments/" + graph_name
              + "_" + str(k) + '_samples.pickle', 'wb') as f:
        pickle.dump(lift.graphlet_samples, f)
    with open("experiments/" + graph_name
              + "_" + str(k) + '_graphletcounts.pickle', 'wb') as f:
        pickle.dump(lift.graphlet_counts, f)

graph_names = [
    "bio-celegansneural",
    "fullb",
    "as-caida",
    "ia-email-univ",
    "socfb-B-anon"
]
for graph_name in graph_names:
    run_longtest(graph_name, 3)
    run_longtest(graph_name, 4)
