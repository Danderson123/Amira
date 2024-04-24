import statistics

class Unitig:

    def __init__(self,
                path,
                path_node_coverages,
                unitig_ID):
        self._path = path
        self._min_coverage = min(path_node_coverages)
        self._max_coverage = max(path_node_coverages)
        self._mean_coverage = statistics.mean(path_node_coverages)
        self._unitig_ID = unitig_ID
        self._terminal_nodes = {path[0], path[1]}

    def get_path(self):
        return self._path

    def get_min_coverage(self):
        return self._min_coverage

    def get_max_coverage(self):
        return self._max_coverage

    def get_mean_coverage(self):
        return self._mean_coverage

    def get_unitig_ID(self):
        return self._unitig_ID

    def get_terminal_nodes(self):
        return self._terminal_nodes

    def shares_terminals(self, other_unitig):
        if len(self.get_terminal_nodes().union(other_unitig.get_terminal_nodes())) > 0:
            return 1
        else:
            return 0