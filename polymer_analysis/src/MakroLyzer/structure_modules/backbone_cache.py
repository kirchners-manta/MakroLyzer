class BackboneCache:
    """
    Shared backbone cache for static-topology analyses.
    Computes longest paths on a reduced graph once and reuses them.
    """

    def __init__(self):
        self._paths = None

    def get_paths(self, graph):
        if self._paths is not None:
            return self._paths

        reduced = graph.remove_1order()
        reduced.surrounding()
        reduced.update_degree()
        subgraphs = reduced.get_subgraphs()

        paths = []
        for subgraph in subgraphs:
            longest = subgraph.find_longest_path()
            if len(longest) < 2:
                continue
            paths.append(longest)

        self._paths = paths
        return paths

    def get_endpoints(self, graph):
        paths = self.get_paths(graph)
        return [(path[0], path[-1]) for path in paths if len(path) >= 2]
