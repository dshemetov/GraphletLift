import networkx as nx
import pynauty as na

g = na.Graph(number_of_vertices=8, directed=False,
          adjacency_dict = { 0: [1,2],
                             1: [0,4],
                             2: [0,3,4,5],
                             3: [2,4,5,6],
                             4: [1,2,3,5],
                             5: [4,2,3,7],
                             6: [3,7],
                             7: [5,6] })

h = na.Graph(number_of_vertices=8, directed=False,
          adjacency_dict = { 0: [1,2],
                             1: [0,4],
                             2: [0,3,4,5],
                             3: [2,4,5,7],
                             4: [1,2,3,5],
                             5: [4,2,3,6],
                             6: [5,7],
                             7: [3,6] })

i = na.Graph(number_of_vertices=8, directed=False,
          adjacency_dict = { 0: [2,4,3,6],
                             1: [5,4],
                             2: [0,3,4,5],
                             3: [0,2,4,7],
                             4: [0,1,2,3],
                             5: [1,2],
                             6: [0,7],
                             7: [3,6] })

j = na.Graph(number_of_vertices=8, directed=False,
          adjacency_dict = { 0: [1,2],
                             1: [0,3],
                             2: [0,3,4,5],
                             3: [1,2,4,5],
                             4: [2,3,5,6],
                             5: [2,3,4,5,6,7],
                             6: [4,5,7],
                             7: [5,6] })

print("g = h? ", na.certificate(g) == na.certificate(h), "g = i? ", na.certificate(g) == na.certificate(i), "g = h? ", na.certificate(g) == na.certificate(j))

num_graphs = 4
random_nx_graphs = [ nx.gnp_random_graph(8,0.4) for i in range(num_graphs) ]
random_na_graphs = [ na.Graph(number_of_vertices = 8, directed = False,
                              adjacency_dict = { n: list(nbrdict.keys()) for n, nbrdict in graph.adjacency() }
                              ) for graph in random_nx_graphs
                   ]

nx_iso = [ nx.is_isomorphic(random_nx_graphs[i],random_nx_graphs[j]) for i in range(num_graphs) for j in range(i,num_graphs) ]
na_iso = [ na.certificate(random_na_graphs[i]) == na.certificate(random_na_graphs[j]) for i in range(num_graphs) for j in range(i,num_graphs) ]

print("Are the isomorphisms the same? ", nx_iso == na_iso)
print("What do the isomorphism matrices look like?\n", nx_iso, na_iso)
print("These are the na graphs:\n",random_na_graphs)
print("These are the nx graphs:\n",random_nx_graphs)
