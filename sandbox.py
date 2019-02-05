# A Hydrogen space for debugging.
import networkx as nx
import lift as lt
import numpy as np
import sys
import inspect

graph = lt.load_graph_fromfile('misc-polblogs')

print(list(nx.isolates(graph)))
