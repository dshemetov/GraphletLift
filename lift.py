"""
Graphlet Lifting Module.
"""

import random
import sympy
import networkx as nx
import networkx.algorithms.isomorphism as iso


# CONSTANTS

# Graph names for small graphs
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


# Helper function for initialization
def load_graph(graph_name):
    """
    Load graph using networkx 'read_edgelist' method.
    All files are organized as a set of edges, possibly with weight values.
    Supported graphs: bio-celegansneural, ia-email-univ, misc-polblogs,
                      misc-as-caida, misc-fullb.
    All graphs were downloaded from http://networkrepository.com/networks.php

    TODO:
        --get rid of networkx
        --load from other graph repositories
    """
    graph = None
    if graph_name == 'bio-celegansneural':
        graph = nx.read_edgelist(
            'Graphs/bio-celegansneural.mtx',
            create_using=nx.Graph(), data=(('weight', float),))

    if graph_name == 'ia-email-univ':
        graph = nx.read_edgelist(
            'Graphs/ia-email-univ.mtx',
            create_using=nx.Graph())

    if graph_name == 'misc-polblogs':
        graph = nx.read_edgelist(
            'Graphs/misc-polblogs.mtx',
            create_using=nx.Graph(), data=(('weight', float),))

    if graph_name == 'misc-as-caida':
        graph = nx.read_edgelist(
            'Graphs/misc-as-caida.mtx',
            create_using=nx.Graph(), data=(('weight', float),))

    if graph_name == 'misc-fullb':
        graph = nx.read_edgelist(
            'Graphs/misc-fullb.mtx',
            create_using=nx.Graph())
        for node in graph.nodes():
            graph.remove_edge(node, node)

    if graph is None:
        raise KeyError("Graph name not found")

    return graph


def get_graphlet_list(k):
    """
    Generate list of all graphlets of size 'k'.
    List is taken from graph_atlas of networkx.
    """
    from networkx.generators.atlas import graph_atlas_g
    assert k > 0
    atlas = graph_atlas_g()[1:]
    graphlet_list = []
    for graph in atlas:
        n = graph.number_of_nodes()
        if n < k:
            continue
        if n > k:
            break
        if nx.is_connected(graph):
            graphlet_list.append(graph)
    return graphlet_list


def get_graphlet_names(k):
    """
    Return a list of names corresponding to graphlets from 'graphlet_list'.
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
    return list(range(len(get_graphlet_list(k))))


# Helper function for shotgun method and probability functions
# in unordered method.
def get_subgraph(graph, nodes):
    """
    Manually construct a subgraph of 'T' given the list of nodes from 'T'.
    Returns a new networkx graph object.

    NOTE:
        Don't use networkx subgraph method because it's very slow.
    """
    list_nodes = list(nodes)
    subgraph = nx.Graph()
    subgraph.add_nodes_from(nodes)
    for i, node in enumerate(list_nodes):
        neighbors = list(graph.neighbors(node))
        for j in range(i, len(list_nodes)):
            if list_nodes[j] in neighbors:
                subgraph.add_edge(node, list_nodes[j])
    return subgraph


def find_type(graph, graphlet_list):
    """
    Given graph T, find an isomorphic graph from 'graphlet_list'.
    Returns the index of the isomorphic graph in 'graphlet_list'.
    """
    edge_num = graph.number_of_edges()
    node_num = graph.number_of_nodes()
    if node_num == 1 or node_num == 2 or node_num == 3:
        return SMALL_ISOMORPHIC_GRAPHS_DICT[(node_num, edge_num)]
    if node_num == 4:
        max_degree = max((graph.degree(node) for node in graph.nodes()))
        if edge_num == 3 or edge_num == 4:
            return SMALL_ISOMORPHIC_GRAPHS_DICT[(node_num, edge_num,
                                                 max_degree)]
        elif edge_num == 5 or edge_num == 6:
            return SMALL_ISOMORPHIC_GRAPHS_DICT[(node_num, edge_num)]
    # Improve matching procedure here for n=5
    for (i, graph_) in enumerate(graphlet_list):
        if iso.GraphMatcher(graph, graph_).is_isomorphic():
            graph_name = i
            break
    return graph_name


# Helper function for 'prob_functions' for unordered method.
def find_type_match(graph, graphlet_list):
    """
    Given a graph, find an isomorphism with one of the canonical graphs from
    'graphlet_list'.
    Return index of the corresponding graph from 'graphlet_list' and a
    match dictionary.
    The match dictionary has format {u_i: v_i}, 'u_i' are nodes from 'graph'
    and 'v_i' are nodes from canonical graph.
    """
    nodes = graph.nodes()
    n = len(nodes)
    if n == 1:
        # trivial graph: just send it to zero!
        return (0, {u: 0 for u in nodes})
    if n == 2:
        # 2-path graph: both nodes are equal, pick a random isomorphism
        return (0, {u: i for i, u in enumerate(nodes)})
    if n == 3:
        if graph.number_of_edges() == 2:
            # wedge-graph: find root, other two are arbitrary
            u0 = next((node for node in nodes if graph.degree(node) == 2))
            (u1, u2) = (node for node in graph.neighbors(u0))
            return (0, {u0: 0, u1: 1, u2: 2})
        if graph.number_of_edges() == 3:
            # triangle: all three are arbitrary
            return (1, {u: i for i, u in enumerate(nodes)})
    if n == 4:
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
    # Improve matching procedure here for n>4.
    for (i, graph_) in enumerate(graphlet_list):
        graph_matcher = iso.GraphMatcher(graph, graph_)
        if graph_matcher.is_isomorphic():
            break
    #assert graph_id[1].is_isomorphic()
    return (i, graph_matcher.mapping)


class LiftGraph(object):
    """
    Helper class that includes various methods needed for lifting.
    The loading of graph is in the init.
    Probability functions are calculated here.
    """
    def __init__(self, graph_name, k,
                 vertex_distribution="edge_uniform"):
        self.graph = load_graph(graph_name)
        self.edge_num = self.graph.size()
        self.node_num = self.graph.number_of_nodes()
        self.k = k
        self.vertex_distribution = vertex_distribution
        self.graphlet_list = get_graphlet_list(k)

    # Helper function for 'prob_functions' method
    def vertex_prob_sympy(self):
        """
        Returns a sympy expression for the probability distribution of nodes.
        Variable 'x_0' represents the degree of the node.
        Possible distributions:
            -- edge_uniform: stationary probability of the standard ranom walk,
            -- node_uniform: uniform distributions among edges.
        """
        if self.vertex_distribution == "edge_uniform":
            return sympy.var('x_0') / (2 * self.edge_num)
        if self.vertex_distribution == "node_uniform":
            return sympy.Integer(1) / self.node_num
        else:
            raise NotImplementedError

    # Helper function for 'lift_shotgun' method
    def vertex_prob(self, vertex_degree):
        """
        Numerical value of the probability of a vertex given its degree.
        """
        if self.vertex_distribution == "edge_uniform":
            return (vertex_degree *
                    (2 * self.edge_num)**(-1))
        if self.vertex_distribution == "node_uniform":
            return self.node_num ** (-1)
        else:
            raise NotImplementedError

    # Helper function for the unordered method.
    def prob_functions(self):
        """
        Construct sympy formulas for graphlet distributions in unordered
        method.
        'k' is the number of nodes in graphlets,
        'vertex distribution' is a sympy formula for the distribution of
        initial vertex.
        Variable 'x_i' corresponds to the degree of i-th node in the graphlet.
        The ordering of nodes in graphlets is the same as in 'graphlet_list'.
        Returns a dictionary {ind: probability}, with 'ind' being the index of
        the graphlet in 'graphlet_list'.
        EXAMPLE:
            python: prob_functions(graph, 3)
            {0: 0.001/(x_0 + x_2 - 2) + 0.001/(x_0 + x_1 - 2),
             1: 0.002/(x_1 + x_2 - 2) + 0.002/(x_0 + x_2 - 2) +
                0.002/(x_0 + x_1 - 2)}
        """
        x = {0: sympy.var('x_0')}
        y = {0: sympy.var('y_0')}
        k = self.k
        graphlet_prob = {0: self.vertex_prob_sympy()}
        for n in range(2, k+1):
            x[n-1] = sympy.var('x_{}'.format(n-1))
            y[n-1] = sympy.var('y_{}'.format(n-1))
            subgraph_prob_weight = graphlet_prob
            graphlet_prob = {}
            for graph_ind, graph in enumerate(get_graphlet_list(n)):
                graphlet_prob[graph_ind] = 0
                for u in graph.nodes():
                    # We sum the conditional probabilities P(S)*P(T|S) for each
                    # connected subgraph S of T.
                    subgraph = get_subgraph(graph, graph.nodes()-{u})
                    if not nx.is_connected(subgraph):
                        continue
                    subgraph_ind, subgraph_match = find_type_match(
                        subgraph, self.graphlet_list)
                    subgraph_prob = (subgraph_prob_weight[subgraph_ind]
                                     .subs({x[i]: y[i] for i in range(n-1)})
                                     .subs({y[j]: x[i] for i, j in
                                            subgraph_match.items()}))
                    subgraph_deg = (sum(x[i] for i in subgraph.nodes()) -
                                    2 * subgraph.number_of_edges())
                    graphlet_prob[graph_ind] += (subgraph_prob *
                                                 graph.degree(u) / subgraph_deg)
        return graphlet_prob

    def random_walk_nodes(self, v0, burn_in):
        """
        Random walk used to pick an initial vertex.
        """
        current_vertex = v0
        for _ in range(burn_in):
            current_vertex = random.choice(list(self.graph
                                                .neighbors(current_vertex)))
        return current_vertex

    def lift_unordered(self, init_vertex):
        """
        Lift procedure for ordered and unordered method.
        To get a 'k'-graphlet from a 'k-1'-graphlet, the next vertex is chosen uniformly from the
        neighbors of the 'k-1'-graphlet.
        """
        graphlet_nodes = set([init_vertex])
        if self.k == 1:
            return graphlet_nodes
        u = init_vertex
        neighbor_list = []
        for _ in range(1, self.k):
            neighbor_list = ([v for v in neighbor_list if v != u]
                             + [v for v in self.graph.neighbors(u) if v not in
                                graphlet_nodes])
            u = random.choice(neighbor_list)
            graphlet_nodes.add(u)
        return graphlet_nodes

    def lift_shotgun(self, init_vertex):
        """
        Lift procedure for the shotgun method.
        Lifts the initial vertex to a graphlet 'S' of size 'k-1'.
        Returns the vertices of graphlet 'S', the probability of 'S', and all
        nodes in the neighborhood of 'S'.
        """
        vertex_neighbors = list(self.graph.neighbors(init_vertex))
        prob = self.vertex_prob(len(vertex_neighbors))
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


class Lift(object):
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
        ...: lift_unordered.graphlet_count(steps_num=15000)
        {'4-chordcycle': 23878, '3-star': 634781, '4-tailedtriangle': 192482,
         '4-cycle': 15824, '4-clique': 2254, '4-path': 514435}
    """

    def __init__(self, graph_name, k, lift_type="unordered",
                 vertex_distribution="edge_uniform"):
        """
        Initialization of Lift object.
        Loading graph, construction of probability functions and different
        lifting procedures are relegated to the 'LiftGraph' class
        """
        self.lift_graph = LiftGraph(graph_name, k, vertex_distribution)
        self.graph = self.lift_graph.graph
        self.k = k
        self.type = lift_type
        self.graphlet_list = self.lift_graph.graphlet_list
        self.graphlet_num = len(self.graphlet_list)
        self.vertex_distribution = vertex_distribution
        if self.type == "unordered":
            self.set_prob_functions()
        else:
            self.prob_functions = None

    def set_prob_functions(self):
        """
        Calculate the symbolic expressions for vertex probability.
        """
        x = [sympy.var('x_{}'.format(i)) for i in range(self.k)]
        self.prob_functions = {ind: sympy.lambdify(x, func)
                               for ind, func in
                               self.lift_graph.prob_functions().items()}

    def change_type(self, lift_type):
        """
        Option to change types without creating a new 'Lift' instance.
        """
        self.type = lift_type
        if (lift_type == "unordered") and (self.prob_functions is None):
            self.set_prob_functions()

    def find_type_prob(self, graphlet_nodes):
        """
        Helper function for an unordered method, uses probability functions
        calculated in 'LiftGraph'.
        """
        assert self.type == "unordered"
        graphlet_type, graphlet_match = find_type_match(
            get_subgraph(self.graph, graphlet_nodes), self.graphlet_list)
        inv_match = {i: j for j, i in graphlet_match.items()}
        degree_list = [self.graph.degree(inv_match[i]) for i in range(self.k)]
        graphlet_prob = self.prob_functions[graphlet_type](*degree_list)
        return (graphlet_type, graphlet_prob)

    def graphlet_count(self, steps_num=1000, burn_in=5, epoch_num=1):
        """
        Wrapper for different 'graphlet_count' methods.
        """
        if self.type == "unordered":
            return self._graphlet_count_unordered(steps_num, burn_in, epoch_num)
        if self.type == "shotgun":
            return self._graphlet_count_shotgun(steps_num, burn_in, epoch_num)
        if self.type == "ordered":
            return self._graphlet_count_ordered(steps_num, burn_in, epoch_num)
        raise ValueError("wrong lift type")

    def _graphlet_count_unordered(self, steps_num, burn_in, epoch_num):
        """
        Unordered lifting method.
        Lifts a graphlet from initial vertex, and calculates its probability.
        """
        assert self.type == "unordered"
        if self.vertex_distribution != "edge_uniform":
            raise NotImplementedError
        for _ in range(epoch_num):
            v0 = random.choice(list(self.graph.nodes()))
            if self.vertex_distribution == "edge_uniform":
                v = self.lift_graph.random_walk_nodes(v0, 100)
            else:
                raise NotImplementedError
            exp_counter = {i: 0 for i in range(self.graphlet_num)}
            for __ in range(steps_num):
                graphlet_nodes = self.lift_graph.lift_unordered(v)
                v = self.lift_graph.random_walk_nodes(v, burn_in)
                graphlet_type, graphlet_prob = self.find_type_prob(
                    graphlet_nodes)
                exp_counter[graphlet_type] += (graphlet_prob)**(-1)
        expectation = {}
        graphlet_names_list = get_graphlet_names(self.k)
        for i in range(self.graphlet_num):
            expectation[graphlet_names_list[i]] = (int(
                exp_counter[i] * (steps_num * epoch_num) ** (-1)))
        return expectation

    def _graphlet_count_ordered(self, steps_num, burn_in, epoch_num):
        raise NotImplementedError

    def _graphlet_count_shotgun(self, steps_num, burn_in, epoch_num):
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
        for _ in range(epoch_num):
            v0 = random.choice(list(self.graph.nodes()))
            if self.vertex_distribution == "edge_uniform":
                v = self.lift_graph.random_walk_nodes(v0, 100)
            else:
                raise NotImplementedError
            exp_counter = {i: 0 for i in range(self.graphlet_num)}
            for __ in range(1, steps_num+1):
                (subgraph, subgraph_prob, neighbor_list) = (self.lift_graph
                                                            .lift_shotgun(v))
                neighbor_set = set(neighbor_list)
                for u in neighbor_set:
                    graph = get_subgraph(self.graph, subgraph.union({u}))
                    graph_type = find_type(graph, self.graphlet_list)
                    exp_counter[graph_type] += (subgraph_prob)**(-1)
                v = self.lift_graph.random_walk_nodes(v, burn_in)
        expectation = {}
        graphlet_names_list = get_graphlet_names(self.k)
        for i in range(self.graphlet_num):
            expectation[graphlet_names_list[i]] = (int(
                exp_counter[i] * (steps_num * co[i]) ** (-1)))
        return expectation
