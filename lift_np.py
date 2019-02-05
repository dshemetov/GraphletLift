"""
Graphlet Lifting Module.
"""
import random
import networkx as nx

# CONSTANTS
# Names for small graphs
SMALL_ISOMORPHIC_GRAPHS_DICT = {(1, 0): 0,
                                (2, 1): 0,
                                (3, 2): 0,
                                (3, 3): 1,
                                (4, 3, 2): 1,
                                (4, 3, 3): 0,
                                (4, 4, 2): 3,
                                (4, 4, 3): 2,
                                (4, 5): 4,
                                (4, 6): 5}
BURN_IN = 20

def get_graphlet_dict(k):
    """
    Generate a dict of lists of all graphlets of size up to 'k', keyed by
    size of graphlet. Uses networkx atlas.
    """
    from networkx.generators.atlas import graph_atlas_g
    assert k > 0
    atlas = graph_atlas_g()[1:]
    graphlet_dict = {i:[] for i in range(1, k+1)}
    for graph in atlas:
        n = graph.number_of_nodes()
        if n > k:
            break
        if nx.is_connected(graph):
            graphlet_dict[n].append(graph)
    return graphlet_dict

def get_graphlet_names(k):
    """
    Return a list of names corresponding to graphlets from 'graphlet_dict'.
    """
    if k == 1:
        return ['singleton']
    elif k == 2:
        return ['2-path']
    elif k == 3:
        return ['2-star', '3-cycle']
    elif k == 4:
        return ['3-star', '4-path', '4-tailedtriangle', '4-cycle',
                '4-chordcycle', '4-clique']
    else:
        graphlet_dict = get_graphlet_dict(k)
        return list(range(len(graphlet_dict[k])))

def get_subgraph(graph, nodes):
    """
    DEPRECATED: I tested the claim that this is faster than the networkx
    subgraph method and found it to be empirically false.
    ===
    Manually constructs the induced subgraph given a list of nodes from the full graph.
    Returns a new networkx graph object.
    Helper function for shotgun method and probability functions in the unordered method.


    NOTE:
        We use this because the networkx subgraph method is very slow.
    """
    list_nodes = list(nodes)
    subgraph = nx.Graph()
    subgraph.add_nodes_from(nodes)
    for i, node in enumerate(list_nodes):
        neighbors = list(graph.neighbors(node))
        for j in range(i+1, len(list_nodes)):
            if list_nodes[j] in neighbors:
                subgraph.add_edge(node, list_nodes[j])
    return subgraph

def find_type_match(nx_graphlet_dict, na_graphlet_cert_dict, graph):
    """
    Given a graph, find an isomorphism with one of the canonical graphs from
    'graphlet_list'.
    Return index of the corresponding graph from 'graphlet_list' and a
    match dictionary.
    The match dictionary has format {u_i: v_i}, 'u_i' are nodes from 'graph'
    and 'v_i' are nodes from canonical graph.
    Helper function for 'prob_functions' for unordered method.
    """
    import networkx.algorithms.isomorphism as iso
    import pynauty as na

    nodes = graph.nodes()
    node_num = len(nodes)
    nx_graphlet_list = nx_graphlet_dict[node_num]

    if node_num == 1:
        # trivial graph: relabel the node to zero.
        return (0, {u: 0 for u in nodes})
    if node_num == 2:
        # 2-path graph: graph is symmetric, choose one of two isomorphisms.
        return (0, {u: i for i, u in enumerate(nodes)})
    if node_num == 3:
        if graph.number_of_edges() == 2:
            # 3-path (or wedge): map the root to zero, the other two are
            # interchangeable.
            u0 = next((node for node in nodes if graph.degree(node) == 2))
            (u1, u2) = (node for node in graph.neighbors(u0))
            return (0, {u0: 0, u1: 1, u2: 2})
        if graph.number_of_edges() == 3:
            # 3-clique (or triangle): all three are interchangeable.
            return (1, {u: i for i, u in enumerate(nodes)})
    if node_num == 4:
        e_num = graph.number_of_edges()
        max_degree = max((graph.degree(node) for node in nodes))
        # 3-star
        if e_num == 3 and max_degree == 3:
            u3 = next((node for node in nodes if graph.degree(node) == 3))
            (u0, u1, u2) = tuple(graph.neighbors(u3))
            return (0, {u0: 0, u1: 1, u2: 2, u3: 3})
        # 4-path
        if e_num == 3 and max_degree == 2:
            (u0, u1) = (node for node in nodes if graph.degree(node) == 2)
            u2 = next((node for node in graph.neighbors(u1) if node != u0))
            u3 = next((node for node in graph.neighbors(u0) if node != u1))
            return (1, {u0: 0, u1: 1, u2: 2, u3: 3})
        # 4-tailedtriangle
        if e_num == 4 and max_degree == 3:
            u3 = next((node for node in nodes if graph.degree(node) == 3))
            (u1, u2) = (node for node in nodes if graph.degree(node) == 2)
            u0 = next((node for node in nodes if graph.degree(node) == 1))
            return (2, {u0: 0, u1: 1, u2: 2, u3: 3})
        # 4-cycle
        if e_num == 4 and max_degree == 2:
            u0 = next((node for node in nodes))
            (u1, u3) = tuple(graph.neighbors(u0))
            u2 = next((node for node in graph.neighbors(u1) if node != u0))
            return (3, {u0: 0, u1: 1, u2: 2, u3: 3})
        # 4-chordcycle
        if e_num == 5:
            (u0, u2) = (node for node in nodes if graph.degree(node) == 3)
            (u1, u3) = (node for node in nodes if graph.degree(node) == 2)
            return (4, {u0: 0, u1: 1, u2: 2, u3: 3})
        # 4-clique
        if e_num == 6:
            (u0, u1, u2, u3) = tuple(nodes)
            return (5, {u0: 0, u1: 1, u2: 2, u3: 3})
        raise ValueError("wrong graphlet format")
    else:
        # Use pynauty for n > 4.
        na_graph = nxgraph_to_relabeled_nagraph(graph)
        (_, ind) = na_graphlet_cert_dict[na.certificate(na_graph)]

        #import pdb; pdb.set_trace()
        matcher = iso.GraphMatcher(graph, nx_graphlet_list[ind])
        mapping = next(matcher.match())
        return (ind, mapping)

def find_type(graph, na_graphlet_cert_dict):
    """
    Given graph T, find an isomorphic graph from 'graphlet_list'.
    Returns the index of the isomorphic graph in 'graphlet_list'.
    """
    import networkx.algorithms.isomorphism as iso
    import pynauty as na

    edge_num = graph.number_of_edges()
    node_num = graph.number_of_nodes()
    #na_graphlet_list = na_graphlet_dict[node_num]

    if node_num == 1 or node_num == 2 or node_num == 3:
        return SMALL_ISOMORPHIC_GRAPHS_DICT[(node_num, edge_num)]
    elif node_num == 4:
        max_degree = max((graph.degree(node) for node in graph.nodes()))
        if edge_num == 3 or edge_num == 4:
            return SMALL_ISOMORPHIC_GRAPHS_DICT[(node_num, edge_num,
                                                 max_degree)]
        elif edge_num == 5 or edge_num == 6:
            return SMALL_ISOMORPHIC_GRAPHS_DICT[(node_num, edge_num)]
    else:
        graph_cert = na.certificate(nxgraph_to_relabeled_nagraph(graph))
        (_, graph_index) = na_graphlet_cert_dict[graph_cert]
        return graph_index
        # Improve matching procedure here for n=5
        # for (i, graph_) in enumerate(graphlet_list):
        #     if iso.GraphMatcher(graph, graph_).is_isomorphic():
        #         graph_name = i
        #         break

def adjacency_to_nagraph(adjacency):
    import pynauty as na
    return na.Graph(number_of_vertices=len(adjacency.keys()),
                    adjacency_dict=adjacency)

def nxgraph_to_relabeled_nagraph(graph):
    node_mapping = {node: i for (i, node) in enumerate(graph.nodes())}
    graph_dict_mapped = {node_mapping[node] : [node_mapping[node2]
                                               for node2 in graph.neighbors(node)]
                         for node in graph.nodes()}
    return adjacency_to_nagraph(graph_dict_mapped)

def load_graph_fromfile(graph_name):
    """
    Load graph using the networkx 'read_edgelist' method.
    All files are organized as a set of edges, possibly with weight values.
    Supported graphs: bio-celegansneural, ia-email-univ, misc-polblogs,
                      misc-as-caida, misc-fullb.
    All graphs were downloaded from http://networkrepository.com/networks.php.

    TODO:
        --get rid of networkx
        --load from other graph repositories

    Trusted - DS.
    """
    graph = None
    if graph_name == 'bio-celegansneural':
        graph = nx.read_edgelist(
            'Graphs/bio-celegansneural.edgelist',
            create_using=nx.Graph())

    if graph_name == 'ia-email-univ':
        graph = nx.read_edgelist(
            'Graphs/ia-email-univ.edgelist',
            create_using=nx.Graph())

    if graph_name == 'misc-fullb':
        graph = nx.read_edgelist(
            'Graphs/fullb.edgelist',
            create_using=nx.Graph())

    if graph_name == 'misc-polblogs':
        graph = nx.read_edgelist(
            'Graphs/polblogs.edgelist',
            create_using=nx.Graph())

    if graph_name == 'as-caida':
        graph = nx.read_edgelist(
            'Graphs/as-caida.edgelist',
            create_using=nx.Graph())

    if graph_name == 'socfb-B-anon':
        graph = nx.read_edgelist(
            'Graphs/socfb-B-anon.edgelist',
            create_using=nx.Graph())

    if graph_name == 'ia-wiki-Talk-dir':
        graph = nx.read_edgelist(
            'Graphs/ia-wiki-Talk-dir.edgelist',
            create_using=nx.Graph())

    if graph is None:
        raise KeyError("Graph name not found")

    return graph

def remove_self_loops(graph):
    """
    Removes self loops from a graph.

    Trusted - DS.
    """
    for node in graph.nodes():
        if graph.has_edge(node, node):
            graph.remove_edge(node, node)

def get_degree_list(graph, graphlet_match):
    invert_match = {j: i for i, j in graphlet_match.items()}
    degree_list = [graph.degree(invert_match[i])
                   for i in range(len(invert_match.keys()))]
    return degree_list

def get_vertex_prob_sympy(graph, vertex_distribution):
    """
    Returns a sympy expression for the probability distribution of nodes.
    Variable 'x_0' represents the degree of the node.
    Possible distributions:
        -- edge_uniform: stationary probability of the standard ranom walk,
        -- node_uniform: uniform distributions among edges.

    Helper function for 'prob_functions' method
    """
    import sympy

    if vertex_distribution == "edge_uniform":
        return sympy.var('x_0') / (2 * graph.size())
    if vertex_distribution == "node_uniform":
        return sympy.Integer(1) / graph.number_of_nodes()
    else:
        raise NotImplementedError

def get_vertex_prob(graph, vertex_degree, vertex_distribution):
    """
    Numerical value of the probability of a vertex given its degree.

    Helper function for 'lift_shotgun' method
    """
    if vertex_distribution == "edge_uniform":
        return (vertex_degree *
                (2 * graph.size()) ** (-1))
    if vertex_distribution == "node_uniform":
        return graph.number_of_nodes() ** (-1)
    else:
        raise NotImplementedError

def sample_vertex(graph, extra=None, burn_in=BURN_IN):
    """
    Samples a vertex from a graph by a random walk of a fixed length.

    Trusted - DS.
    """
    vertex0 = random.choice(list(graph.nodes()))
    vertex = randomwalk_from_vertex(graph, vertex0, burn_in)
    return vertex

def randomwalk_from_vertex(graph, vertex0, burn_in=BURN_IN):
    """
    Random walk from a vertex.

    Trusted - DS.
    """
    current_vertex = vertex0
    for _ in range(burn_in):
        current_vertex = random.choice(list(graph.neighbors(current_vertex)))
    return current_vertex

def get_graphlet_cert_dict(na_graphlet_dict):
    """
    Builds the certificate dictionary for the pynauty graphlets.
    """
    import pynauty as na
    na_graphlet_cert_dict = {}
    for num_nodes in na_graphlet_dict.keys():
        for (ind, na_graphlet) in enumerate(na_graphlet_dict[num_nodes]):
            na_graphlet_cert_dict[na.certificate(na_graphlet)] = (num_nodes,
                                                                  ind)

    return na_graphlet_cert_dict

def get_graphlet_prob(
        graph, nx_graphlet_dict, na_graphlet_cert_dict,
        prob_functions, graphlet_nodes):
    """
    Helper function for the unordered method. Uses probability functions.
    """
    subgraph = nx.Graph(graph.subgraph(graphlet_nodes))
    graphlet_type, graphlet_match = find_type_match(
        nx_graphlet_dict, na_graphlet_cert_dict, subgraph
        )
    # invert_match = {j: i for i, j in graphlet_match.items()}
    # degree_list = [self.graph.degree(invert_match[i])
    #                for i in range(self.k)]
    degree_list = get_degree_list(graph, graphlet_match)
    # import pdb; pdb.set_trace()
    prob_function = prob_functions[graphlet_type]
    graphlet_prob = prob_function(*degree_list)
    return (graphlet_type, graphlet_prob)

def sample_unordered_lift(
        graph, k, vertex, burn_in=BURN_IN):
    """
    Attempts a lift at the vertex v. If stuck in a disconnected component of
    size less than k, the vertex is resampled.
    Returns a list of graphlet nodes.

    Trusted - DS.
    """
    graphlet_nodes = sample_unordered_lift_once(graph, k, vertex)
    while len(graphlet_nodes) != k:
        vertex = sample_vertex(graph)
        graphlet_nodes = sample_unordered_lift_once(graph, k, vertex)
    return graphlet_nodes

def sample_unordered_lift_once(graph, k, init_vertex):
    """
    Lift procedure for ordered and unordered method.
    To get a 'k'-graphlet from a 'k-1'-graphlet, the next vertex is chosen
    uniformly from the neighbors of the 'k-1'-graphlet.

    Trusted - DS.
    """
    graphlet_nodes = set([init_vertex])
    if k == 1:
        return graphlet_nodes
    u = init_vertex
    neighbor_list = []
    for _ in range(1, k):
        neighbor_list = ([v for v in neighbor_list if v != u]
                         + [v for v in graph.neighbors(u)
                            if v not in graphlet_nodes])
        u = random.choice(neighbor_list)
        graphlet_nodes.add(u)
    return graphlet_nodes

class Lift():
    """
    A class for running a 'graphlet_count' method on the graph.
    Arguments for initialization:
        -- 'graph_name'- string of one of the graph names available in
        'load_graph',
        -- 'k'- the size of graphlets to be counted,
        -- 'lift_type'- different lifting methods: ordered, unordered(default)
        or shotgun
        -- 'vertex_distribution'- option for different vertex distributions.
    EXAMPLE:
        python: graph_name = "bio-celegansneural"
        ...: lift_unordered = Lift(graph_name, k=4, lift_type="unordered")
        ...: lift_unordered.graphlet_count(num_steps=15000)
        {'4-chordcycle': 23878, '3-star': 634781, '4-tailedtriangle': 192482,
         '4-cycle': 15824, '4-clique': 2254, '4-path': 514435}
    """
    def __init__(
            self, graph, k, lift_type="unordered",
            vertex_choice = "uniform"):
        """
        Prepares a graph for lifting. Loads the graph, computes the probability
        functions, and prepares a list of k-graphlets for isomorphism.

        Trusted - DS.
        """
        if isinstance(graph, str):
            self.graph = load_graph_fromfile(graph)
        else:
            self.graph = graph
        if k > self.graph.number_of_nodes():
            raise ValueError("Graphlet size bigger than graph size.")

        self.k = k
        self.type = lift_type
        self.vertex_choice = vertex_choice
        # Build graphlet library
        self.nx_graphlet_dict = get_graphlet_dict(k)
        self.na_graphlet_dict = {i:[nxgraph_to_relabeled_nagraph(graph)
                                    for graph in self.nx_graphlet_dict[i]]
                                 for i in range(1, k+1)}
        self.na_graphlet_cert_dict = get_graphlet_cert_dict(self
                                                            .na_graphlet_dict)
        self.graphlet_counts = None
        self.total_samples = 0
        self.graphlet_samples = []

        if self.type == "unordered":
            self.set_prob_functions()
        else:
            self.prob_functions = None

    def get_graphlet_count(
            self, num_steps=1000, burn_in=BURN_IN,
            num_epoch=1):
        """
        Wrapper for different 'get_graphlet_count' methods.

        Trusted - DS.
        """
        if self.type == "unordered":
            samples, graphlet_samples = self._get_graphlet_count_unordered(
                num_steps, burn_in)
            self.graphlet_samples = zip(samples, graphlet_samples)
            #import pdb; pdb.set_trace()
            self.update_counts(samples)
            if self.k <= 4:
                renamed_counts = {
                    get_graphlet_names(self.k)[key]: self.graphlet_counts[key]
                    for key in self.graphlet_counts.keys()
                    }
                return renamed_counts
            else:
                return self.graphlet_counts
        if self.type == "shotgun":
            return self._get_graphlet_count_shotgun(
                num_steps, burn_in)
        if self.type == "ordered":
            return self._get_graphlet_count_ordered(
                num_steps, burn_in)
        raise ValueError("wrong lift type")

    def update_counts(self, samples):
        if self.graphlet_counts is None:
            graphlet_num = len(self.nx_graphlet_dict[self.k])
            graphlet_counts = {
                i: 0 for i in range(graphlet_num)
                }
            for sample in samples:
                graphlet_type, graphlet_prob = sample[0], sample[1]
                graphlet_counts[graphlet_type] += 1 / graphlet_prob
            for graphlet_type in graphlet_counts.keys():
                graphlet_counts[graphlet_type] = int(round(
                    graphlet_counts[graphlet_type] / len(samples)
                ))
            self.graphlet_counts = graphlet_counts
            self.total_samples = len(samples)
        else:
            graphlet_num = len(self.nx_graphlet_dict[self.k])
            graphlet_counts = {
                i: 0 for i in range(graphlet_num)
                }
            for sample in samples:
                graphlet_type, graphlet_prob = sample[0], sample[1]
                graphlet_counts[graphlet_type] += 1 / graphlet_prob
            for graphlet_type in graphlet_counts.keys():
                self.graphlet_counts[graphlet_type] = int(round(
                    (self.graphlet_counts[graphlet_type]*self.total_samples
                     + graphlet_counts[graphlet_type])
                    / (self.total_samples + len(samples))
                    ))
            self.total_samples += len(samples)

    def _get_graphlet_count_unordered(
            self, num_steps, burn_in=BURN_IN):
        """
        Unordered lifting method.
        Lifts a graphlet from initial vertex and calculates its frequency.

        Trusted - DS.
        """
        # Sample vertices.
        if self.vertex_choice == "random walk":
            vertices = [sample_vertex(self.graph, burn_in=BURN_IN)]
            for _ in range(num_steps-1):
                vertices.append(randomwalk_from_vertex(
                    self.graph, vertices[-1], burn_in))
        else:
            vertices = [sample_vertex(self.graph, burn_in=BURN_IN)
                        for i in range(num_steps)]

        # Non-parallel code.
        # Sample graphlets at those vertices.
        graphlet_samples = [
            sample_unordered_lift(self.graph, self.k, vertex, burn_in)
            for vertex in vertices
        ]

        # Convert graphlets to probabilities and types.
        # Non-parallel code.
        samples = [
            get_graphlet_prob(
                self.graph, self.nx_graphlet_dict,
                self.na_graphlet_cert_dict,
                self.prob_functions, graphlet_sample)
            for graphlet_sample in graphlet_samples
            ]

        return (samples, graphlet_samples)

    def _get_graphlet_count_ordered(self, num_steps,
                                    burn_in=BURN_IN):
        raise NotImplementedError

    def _get_graphlet_count_shotgun(
            self, num_steps, burn_in=BURN_IN):
        """
        Shotgun lifting method.
        Lifts a vertex to 'k-1' graphlet, and inludes all its neighbors into
        the estimation.
        """
        graphlet_num = len(self.nx_graphlet_dict[self.k])

        if self.k == 4:
            co = {0: 12, 1: 8, 2: 12, 3: 16, 4: 20, 5: 24}
        elif self.k == 3:
            co = {0: 4, 1: 6}
        elif self.k == 2:
            co = {0: 2}

        v = sample_vertex(self.graph, burn_in=BURN_IN)
        graphlet_counts = {i: 0 for i in range(graphlet_num)}
        for __ in range(num_steps):
            (subgraph, subgraph_prob,
             neighbor_list) = (self.sample_lift_shotgun(v))
            neighbor_set = set(neighbor_list)
            for u in neighbor_set:
                graph = nx.Graph(self.graph.subgraph(subgraph.union({u})))
                graph_type = find_type(graph, self.na_graphlet_cert_dict)
                graphlet_counts[graph_type] += (subgraph_prob)**(-1)
            v = randomwalk_from_vertex(self.graph, v, burn_in)

        expectation = {}
        graphlet_names_list = get_graphlet_names(self.k)
        for i in range(graphlet_num):
            expectation[graphlet_names_list[i]] = (int(
                graphlet_counts[i] * (num_steps * co[i]) ** (-1)))
        return expectation

    def sample_lift_shotgun(self, init_vertex):
        """
        Lift procedure for the shotgun method.
        Lifts the initial vertex to a graphlet 'S' of size 'k-1'.
        Returns the vertices of graphlet 'S', the probability of 'S', and all
        nodes in the neighborhood of 'S'.
        """
        vertex_neighbors = list(self.graph.neighbors(init_vertex))
        prob = get_vertex_prob(
            self.graph, len(vertex_neighbors), "edge_uniform")
        graphlet = set([init_vertex])
        e_num = 0
        if self.k == 1:
            return (graphlet, e_num, prob, vertex_neighbors)
        u = init_vertex
        neighbor_list = vertex_neighbors
        for _ in range(1, self.k-1):
            u = random.choice(neighbor_list)
            vertex_neighbors = list(self.graph.neighbors(u))
            subgraph_degree = len([v for v in vertex_neighbors
                                   if v in graphlet])
            prob = prob * subgraph_degree * (len(neighbor_list))**(-1)
            neighbor_list = ([v for v in neighbor_list if v != u]
                             + [v for v in vertex_neighbors if v not in
                                graphlet])
            graphlet.add(u)
            e_num += subgraph_degree
        return (graphlet, prob, neighbor_list)

    def set_prob_functions(self):
        """
        Calculate the symbolic expressions for vertex probability.
        """
        import sympy

        variables = [sympy.var('x_{}'.format(i)) for i in range(self.k)]
        self.prob_functions = {ind: sympy.lambdify(variables, func)
                               for ind, func in
                               self.calculate_prob_functions().items()}

    def calculate_prob_functions(self):
        """
        Construct sympy formulas for graphlet distributions in unordered
        method.
        'k' is the number of nodes in graphlets,
        'vertex distribution' is a sympy formula for the distribution of
        initial vertex.
        Variable 'x_i' corresponds to the degree of i-th node in the graphlet.
        The ordering of nodes in graphlets is the same as in 'graphlet_dict'.
        Returns a dictionary {ind: probability}, with 'ind' being the index of
        the graphlet in 'graphlet_dict'.
        EXAMPLE:
            > graph = nx.path_graph(5)
            > prob_functions(graph, 3)
            {0: 1/(4*(x_0 + x_2 - 2)) + 1/(4*(x_0 + x_1 - 2)),
             1: 1/(2*(x_1 + x_2 - 2)) + 1/(2*(x_0 + x_2 - 2))
                + 1/(2*(x_0 + x_1 - 2))

        Helper function for the unordered method.
        """
        import sympy

        k = self.k
        x = {n: sympy.var('x_{}'.format(n)) for n in range(k+1)}
        y = {n: sympy.var('y_{}'.format(n)) for n in range(k+1)}
        graphlet_probs = {
            0: get_vertex_prob_sympy(self.graph, "edge_uniform")
            }
        for n in range(2, k+1):
            subgraph_probs = graphlet_probs # Contains P(S_(n-1))
            graphlet_probs = {} # Builds P(S_n)
            for graph_index, graph in enumerate(self.nx_graphlet_dict[n]):
                graphlet_probs[graph_index] = 0
                for u in graph.nodes():
                    # We sum the conditional probabilities
                    # P(S_(n-1)) * P(S_n | S_(n-1)) for each connected subgraph
                    # S_(n-1) of S_n and for n in {2, ..., k}.
                    subgraph = nx.Graph(graph.subgraph(graph.nodes()-{u}))
                    if not nx.is_connected(subgraph):
                        continue
                    subgraph_index, subgraph_match = find_type_match(
                        self.nx_graphlet_dict,
                        self.na_graphlet_cert_dict,
                        subgraph
                        )
                    subgraph_prob = (subgraph_probs[subgraph_index]
                                     .subs({x[i]: y[i] for i in range(n-1)})
                                     .subs({y[j]: x[i] for i, j in
                                            subgraph_match.items()}))
                    subgraph_deg = (sum(x[i] for i in subgraph.nodes()) -
                                    2 * subgraph.number_of_edges())
                    graphlet_probs[graph_index] += (subgraph_prob *
                                                    graph.degree(u) /
                                                    subgraph_deg)
        return graphlet_probs
