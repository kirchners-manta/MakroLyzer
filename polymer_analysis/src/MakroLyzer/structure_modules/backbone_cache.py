class BackboneCache:
    """
    Shared backbone cache for static-topology analyses.
    Computes longest paths on a reduced graph once and reuses them.
    """

    def __init__(self):
        self._paths = None
        self._subgraph_paths = None

    def get_paths(self, graph):
        if self._paths is None:
            self._compute(graph)
        return self._paths

    def get_subgraph_paths(self, graph):
        if self._subgraph_paths is None:
            self._compute(graph)
        return self._subgraph_paths

    def _compute(self, graph):
        reduced = graph.remove_1order()
        reduced.surrounding()
        reduced.update_degree()
        subgraphs = reduced.get_subgraphs()

        paths = []
        subgraph_paths = []
        for subgraph in subgraphs:
            longest = subgraph.find_longest_path()
            if len(longest) < 2:
                continue
            paths.append(longest)
            subgraph_paths.append((subgraph, longest))

        self._paths = paths
        self._subgraph_paths = subgraph_paths

    def get_endpoints(self, graph):
        paths = self.get_paths(graph)
        return [(path[0], path[-1]) for path in paths if len(path) >= 2]
