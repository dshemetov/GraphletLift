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
TRANSIENT_LENGTH = 20

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
        return ['wedge', 'triangle']
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

def find_type_match(graph, nx_graphlet_dict, na_graphlet_cert_dict):
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
            # wedge-graph: map the root to zero, the other two are
            # interchangeable.
            u0 = next((node for node in nodes if graph.degree(node) == 2))
            (u1, u2) = (node for node in graph.neighbors(u0))
            return (0, {u0: 0, u1: 1, u2: 2})
        if graph.number_of_edges() == 3:
            # triangle: all three are interchangeable.
            return (1, {u: i for i, u in enumerate(nodes)})
    if node_num == 4:
        e_num = graph.number_of_edges()
        max_degree = max((graph.degree(node) for node in nodes))
        if e_num == 3 and max_degree == 3:
            u3 = next((node for node in nodes if graph.degree(node) == 3))
            (u0, u1, u2) = tuple(graph.neighbors(u3))
            return (0, {u0: 0, u1: 1, u2: 2, u3: 3})
        if e_num == 3 and max_degree == 2:
            (u0, u1) = (node for node in nodes if graph.degree(node) == 2)
            u2 = next((node for node in graph.neighbors(u1) if node != u0))
            u3 = next((node for node in graph.neighbors(u0) if node != u1))
            return (1, {u0: 0, u1: 1, u2: 2, u3: 3})
        if e_num == 4 and max_degree == 3:
            u3 = next((node for node in nodes if graph.degree(node) == 3))
            (u1, u2) = (node for node in nodes if graph.degree(node) == 2)
            u0 = next((node for node in nodes if graph.degree(node) == 1))
            return (2, {u0: 0, u1: 1, u2: 2, u3: 3})
        if e_num == 4 and max_degree == 2:
            u0 = next((node for node in nodes))
            (u1, u3) = tuple(graph.neighbors(u0))
            u2 = next((node for node in graph.neighbors(u1) if node != u0))
            return (3, {u0: 0, u1: 1, u2: 2, u3: 3})
        if e_num == 5:
            (u0, u2) = (node for node in nodes if graph.degree(node) == 3)
            (u1, u3) = (node for node in nodes if graph.degree(node) == 2)
            return (4, {u0: 0, u1: 1, u2: 2, u3: 3})
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
            'Graphs/misc-fullb.edgelist',
            create_using=nx.Graph())

    if graph_name == 'misc-polblogs':
        graph = nx.read_edgelist(
            'Graphs/misc-polblogs.edgelist',
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

def get_vertex_prob(graph, vertex_distribution, vertex_degree):
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

def sample_vertex(graph, vertex_distribution="edge_uniform",
                         transient_length=TRANSIENT_LENGTH):
    """
    Samples a vertex from a graph by a random walk of a fixed length.

    Trusted - DS.
    """
    vertex0 = random.choice(list(graph.nodes()))
    if vertex_distribution == "edge_uniform":
        vertex = random_walk_nodes(graph, vertex0, transient_length)
    else:
        raise NotImplementedError
    return vertex

def random_walk_nodes(graph, vertex0, transient_length=TRANSIENT_LENGTH):
    """
    Random walk used to pick an initial vertex.

    Trusted - DS.
    """
    current_vertex = vertex0
    for _ in range(transient_length):
        current_vertex = random.choice(list(graph.neighbors(current_vertex)))
    return current_vertex

def get_graphlet_cert_dict(na_graphlet_dict):
    import pynauty as na
    na_graphlet_cert_dict = {}
    for num_nodes in na_graphlet_dict.keys():
        for (ind, na_graphlet) in enumerate(na_graphlet_dict[num_nodes]):
            na_graphlet_cert_dict[na.certificate(na_graphlet)] = (num_nodes,
                                                                  ind)

    return na_graphlet_cert_dict


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
    def __init__(self, graph, k, lift_type="unordered",
                 vertex_distribution="edge_uniform"):
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
            raise ValueError("Graphlets bigger than size of graph.")

        self.k = k
        self.type = lift_type
        self.vertex_distribution = vertex_distribution
        # Build graphlet library
        self.nx_graphlet_dict = get_graphlet_dict(k)
        self.na_graphlet_dict = {i:[nxgraph_to_relabeled_nagraph(graph)
                                    for graph in self.nx_graphlet_dict[i]]
                                 for i in range(1, k+1)}
        self.na_graphlet_cert_dict = get_graphlet_cert_dict(self
                                                            .na_graphlet_dict)

        self.graphlet_num = len(self.nx_graphlet_dict[self.k])
        self.exp_counter = None

        if self.type == "unordered":
            self.set_prob_functions()
        else:
            self.prob_functions = None

    def change_type(self, lift_type):
        """
        Option to change types without creating a new 'Lift' instance.

        Untrusted - DS.
        """
        self.type = lift_type
        if (lift_type == "unordered") and (self.prob_functions is None):
            self.set_prob_functions()

    def graphlet_count(self, num_steps=1000,
                       transient_length=TRANSIENT_LENGTH, num_epoch=1):
        """
        Wrapper for different 'graphlet_count' methods.

        Trusted - DS.
        """
        if self.type == "unordered":
            return self._graphlet_count_unordered(num_steps, transient_length,
                                                  num_epoch)
        if self.type == "shotgun":
            return self._graphlet_count_shotgun(num_steps, transient_length,
                                                num_epoch)
        if self.type == "ordered":
            return self._graphlet_count_ordered(num_steps, transient_length,
                                                num_epoch)
        raise ValueError("wrong lift type")

    def _graphlet_count_unordered(self, num_steps,
                                  transient_length=TRANSIENT_LENGTH,
                                  num_epoch=1):
        """
        Unordered lifting method.
        Lifts a graphlet from initial vertex and calculates its frequency.

        Trusted - DS.
        """
        if self.vertex_distribution != "edge_uniform":
            raise NotImplementedError

        # Perform sampling.
        for _ in range(num_epoch):
            vertex = sample_vertex(self.graph, self.vertex_distribution,
                                   transient_length)
            exp_counter = {i: 0 for i in range(self.graphlet_num)}
            for __ in range(num_steps):
                graphlet_nodes = self.sample_unordered_lift(vertex,
                                                            transient_length)
                graphlet_type, graphlet_probs = self.get_graphlet_prob(
                                                     graphlet_nodes)
                exp_counter[graphlet_type] += (graphlet_probs)**(-1)
                vertex = random_walk_nodes(self.graph, vertex, transient_length)

        # Calculate the expected frequency.
        expectation = {}
        graphlet_names_list = get_graphlet_names(self.k)
        for i in range(self.graphlet_num):
            expectation[graphlet_names_list[i]] = (int(round(
                exp_counter[i] * (num_steps * num_epoch) ** (-1))))

        self.exp_counter = exp_counter
        return expectation

    def _graphlet_count_ordered(self, num_steps,
                                transient_length=TRANSIENT_LENGTH,
                                num_epoch=1):
        raise NotImplementedError

    def _graphlet_count_shotgun(self, num_steps,
                                transient_length=TRANSIENT_LENGTH,
                                num_epoch=1):
        """
        Shotgun lifting method.
        Lifts a vertex to 'k-1' graphlet, and inludes all its neighbors into
        the estimation.
        """
        if self.k == 4:
            co = {0: 12, 1: 8, 2: 12, 3: 16, 4: 20, 5: 24}
        elif self.k == 3:
            co = {0: 4, 1: 6}
        elif self.k == 2:
            co = {0: 2}
        for _ in range(num_epoch):
            v = sample_vertex(self.graph, self.vertex_distribution,
                              transient_length)
            exp_counter = {i: 0 for i in range(self.graphlet_num)}
            for __ in range(num_steps):
                (subgraph, subgraph_prob,
                 neighbor_list) = (self.sample_lift_shotgun(v))
                neighbor_set = set(neighbor_list)
                for u in neighbor_set:
                    graph = nx.Graph(self.graph.subgraph(subgraph.union({u})))
                    graph_type = find_type(graph, self.na_graphlet_cert_dict)
                    exp_counter[graph_type] += (subgraph_prob)**(-1)
                v = random_walk_nodes(self.graph, v, transient_length)
        expectation = {}
        graphlet_names_list = get_graphlet_names(self.k)
        for i in range(self.graphlet_num):
            expectation[graphlet_names_list[i]] = (int(
                exp_counter[i] * (num_steps * co[i]) ** (-1)))
        return expectation

    def sample_unordered_lift(self, vertex, transient_length=TRANSIENT_LENGTH):
        """
        Attempts a lift at the vertex v. If a graphlet of smaller size is
        obtained, the node is resampled until we get a graphlet of size k.
        Returns a list of graphlet nodes.

        Trusted - DS.
        """
        graphlet_nodes = self.sample_unordered_lift_once(vertex)
        while len(graphlet_nodes) != self.k:
            vertex = sample_vertex(self.graph, self.vertex_distribution,
                                   transient_length)
            graphlet_nodes = self.sample_unordered_lift_once(vertex)
        return graphlet_nodes

    def sample_unordered_lift_once(self, init_vertex):
        """
        Lift procedure for ordered and unordered method.
        To get a 'k'-graphlet from a 'k-1'-graphlet, the next vertex is chosen
        uniformly from the neighbors of the 'k-1'-graphlet.

        Trusted - DS.
        """
        graphlet_nodes = set([init_vertex])
        if self.k == 1:
            return graphlet_nodes
        u = init_vertex
        neighbor_list = []
        for _ in range(1, self.k):
            neighbor_list = ([v for v in neighbor_list if v != u]
                             + [v for v in self.graph.neighbors(u)
                                if v not in graphlet_nodes])
            u = random.choice(neighbor_list)
            graphlet_nodes.add(u)
        return graphlet_nodes

    def sample_lift_shotgun(self, init_vertex):
        """
        Lift procedure for the shotgun method.
        Lifts the initial vertex to a graphlet 'S' of size 'k-1'.
        Returns the vertices of graphlet 'S', the probability of 'S', and all
        nodes in the neighborhood of 'S'.
        """
        vertex_neighbors = list(self.graph.neighbors(init_vertex))
        prob = get_vertex_prob(self.graph, self.vertex_distribution,
                               len(vertex_neighbors))
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
        graphlet_probs = {0: get_vertex_prob_sympy(self.graph,
                                                   self.vertex_distribution)}
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
                    subgraph_index, subgraph_match = find_type_match(subgraph,
                                                    self.nx_graphlet_dict,
                                                    self.na_graphlet_cert_dict)
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

    def get_graphlet_prob(self, graphlet_nodes):
        """
        Helper function for an unordered method, uses probability functions.
        """
        assert self.type == "unordered"
        subgraph = nx.Graph(self.graph.subgraph(graphlet_nodes))
        graphlet_type, graphlet_match = find_type_match(subgraph,
                                                    self.nx_graphlet_dict,
                                                    self.na_graphlet_cert_dict)
        # invert_match = {j: i for i, j in graphlet_match.items()}
        # degree_list = [self.graph.degree(invert_match[i])
        #                for i in range(self.k)]
        degree_list = get_degree_list(self.graph, graphlet_match)
        # import pdb; pdb.set_trace()
        prob_function = self.prob_functions[graphlet_type]
        graphlet_prob = prob_function(*degree_list)
        return (graphlet_type, graphlet_prob)
