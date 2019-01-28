import networkx as nx


def remove_self_loops(graph):
    """
    Removes self loops from a graph.

    Trusted - DS.
    """
    for node in graph.nodes():
        if graph.has_edge(node, node):
            graph.remove_edge(node, node)

def convert_mtx(file_in, file_out1, file_out2):
    with open(file_out1, 'w') as fout1:
        with open(file_out2, 'w') as fout2:
            with open(file_in, 'r') as fin:
                lines = fin.readlines()
                comment = True
                i = 0
                while comment:
                    if lines[i][0] == "%":
                        i += 1
                        continue
                    else:
                        comment = False

                graph = nx.parse_edgelist(
                    lines[i+1:],
                    create_using=nx.Graph(),
                    data=False
                    )
                remove_self_loops(graph)
                fout1.write(str(graph.number_of_nodes()))
                fout1.write(" ")
                fout1.write(str(graph.number_of_edges()))
                fout1.write("\n")
                for edge in graph.edges():
                    fout1.write(str(edge[0]))
                    fout1.write(" ")
                    fout1.write(str(edge[1]))
                    fout1.write("\n")
                    fout2.write(str(edge[0]))
                    fout2.write(" ")
                    fout2.write(str(edge[1]))
                    fout2.write("\n")

directory = "/Users/dmitron/Documents/GradSchool/GraphletLift/code/LiftSRW/Graphs/"
network_files = [
    "bio-celegansneural",
    "fullb",
    "ia-email-univ",
    "polblogs"
]
for network_file in network_files:
    convert_mtx(
        directory + network_file + ".mtx",
        directory + network_file + ".edges",
        directory + network_file + ".edgelist"
        )
