import unittest

from amira.construct_graph import GeneMerGraph
from amira.path_finding_operations import (
    cluster_adjacent_paths,
    construct_suffix_tree
)
from suffix_tree import Tree


class TestPathFindingConstructor(unittest.TestCase):

    def test___find_full_paths_linear_simple(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene6",
            "+gene7",
        ]
        genes2 = [
            "-gene2",
            "+gene3",
            "-gene4",
        ]
        genes3 = [
            "-gene4",
            "-gene6",
            "+gene7",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes3}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 1)
        # assertion
        self.assertEqual(len(f_paths), 1)
        for k in f_paths:
            self.assertEqual(len(f_paths[k]), 4)

    def test___find_full_paths_no_adjacent_paths(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene4",
            "-gene4",
            "-gene4",
            "+gene7",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes1, "read4": genes1}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 1)
        # assertion
        self.assertEqual(len(f_paths), 1)
        for k in f_paths:
            self.assertEqual(len(f_paths[k]), 4)

    def test___find_full_paths_one_adjacent_path(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene4", "-gene4", "-gene4", "+gene7", "-gene8"]
        genes2 = [
            "+gene-1",
            "-gene0",
            "+gene1",
            "-gene2",
            "+gene4",
            "-gene4",
            "-gene4",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10" "+gene11",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 1)
        # assertion
        self.assertEqual(len(f_paths), 1)
        for k in f_paths:
            self.assertEqual(len(f_paths[k]), 4)

    def test___find_full_paths_linear_path_duplicates_simple(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene4",
            "-gene6",
            "+gene7",
        ]
        genes2 = [
            "-gene2",
            "+gene3",
            "-gene4",
        ]
        genes3 = [
            "+gene4",
            "-gene6",
            "+gene7",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes3}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 1)
        # assertion
        self.assertEqual(len(f_paths), 1)
        for k in f_paths:
            self.assertEqual(len(k), 7)
            self.assertEqual(len(f_paths[k]), 2)

    def test___find_full_paths_linear_path_contained(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "+gene7", "-gene8"]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 1)
        # assertion
        self.assertEqual(len(f_paths), 2)
        for k in f_paths:
            self.assertTrue(len(k) == 11 or len(k) == 8)
            self.assertEqual(len(f_paths[k]), 2)

    def test___find_full_paths_linear_path_contained_two(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene3",
            "-gene4",
            "+gene5",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 1)
        # assertion
        self.assertEqual(len(f_paths), 2)
        for k in f_paths:
            self.assertTrue(len(k) == 11 or len(k) == 14)
            self.assertEqual(len(f_paths[k]), 2)

    def test___find_full_paths_terminate_at_junction(self):
        # setup
        genes1 = [
            "-trbC",
            "-trbB",
            "-group_1081",
            "-group_6156",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "-group_5175",
            "+group_5625",
        ]
        genes2 = [
            "-alkB",
            "-ada",
            "-apbE",
            "-ompC",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+rcsD",
            "+rcsB",
            "-rcsC",
            "+atoS",
            "+atoC",
            "+atoD",
            "+atoA",
            "+atoE",
            "+atoB",
            "-yfaP",
            "-yfaQ",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "blaCMY54NG_0488491", 1)
        # assertion
        self.assertEqual(len(f_paths), 2)
        for k in f_paths:
            self.assertTrue(len(k) == 18 or len(k) == 24)
            self.assertEqual(len(f_paths[k]), 2)

    def test___find_full_paths_terminate_and_start_at_junction(self):
        # setup
        genes1 = [
            "-trbC",
            "-trbB",
            "-group_1081",
            "-group_6156",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+sugE",
            "-blc",
            "-group_5175",
            "+group_5625",
        ]
        genes2 = [
            "-alkB",
            "-ada",
            "-apbE",
            "-ompC",
            "+sugE",
            "-blc",
            "-blaCMY54NG_0488491",
            "+rcsD",
            "+rcsB",
            "-rcsC",
            "+atoS",
            "+atoC",
            "+atoD",
            "+atoA",
            "+atoE",
            "+atoB",
            "-yfaP",
            "-yfaQ",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "blaCMY54NG_0488491", 1)
        # assertion
        self.assertEqual(len(f_paths), 2)
        for k in f_paths:
            self.assertTrue(len(k) == 18 or len(k) == 26)
            self.assertEqual(len(f_paths[k]), 2)

    def test___find_full_paths_branching_path(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene5",
            "-gene6",
            "+gene7",
            "-gene4",
            "-gene6",
            "+gene7",
            "-gene10",
            "-gene11",
        ]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene5",
            "-gene6",
            "+gene7",
            "-gene4",
            "-gene13",
            "+gene14",
            "-gene15",
            "-gene16",
        ]
        genes3 = [
            "+gene7",
            "-gene4",
            "-gene13",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes2, "read5": genes3}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 1)
        # assertion
        self.assertEqual(len(f_paths), 2)
        for k in f_paths:
            self.assertTrue(len(k) == 3 or len(k) == 11)
            self.assertTrue(len(f_paths[k]) == 2 or len(f_paths[k]) == 3)

    def test___find_full_paths_triangle(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene7",
            "+gene5",
            "+gene7",
            "+gene5",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        genes2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        genes3 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "+gene7",
            "+gene5",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
        ]
        graph = GeneMerGraph(
            {
                "read1": genes1,
                "read2": genes1,
                "read3": genes2,
                "read4": genes2,
                "read5": genes3,
                "read6": genes3,
            },
            3,
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 1)
        # assertion
        self.assertEqual(len(f_paths), 3)
        for k in f_paths:
            self.assertTrue(len(k) == 12 or len(k) == 14 or len(k) == 10)
            self.assertEqual(len(f_paths[k]), 2)

    def test___find_full_paths_linear_path_duplicates_long_reads(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "+gene4",
            "-gene10",
            "+gene11",
            "-gene12",
        ]
        genes2 = [
            "-gene2",
            "+gene3",
            "-gene4",
        ]
        genes3 = [
            "+gene4",
            "-gene10",
            "+gene11",
        ]
        graph = GeneMerGraph(
            {"read1": genes1, "read2": genes1, "read3": genes2, "read4": genes3}, 3
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 1)
        # assertion
        self.assertEqual(len(f_paths), 1)
        for k in f_paths:
            self.assertEqual(len(k), 13)
            self.assertEqual(len(f_paths[k]), 2)

    def test___find_full_paths_complex_one(self):
        # setup
        import json

        call_file = "tests/complex_gene_calls_one.json"
        position_file = "tests/complex_gene_positions_one.json"
        with open(call_file) as i:
            calls = json.load(i)
        with open(position_file) as i:
            positions = json.load(i)
        filtered_calls = {}
        for r in calls:
            if any(g[1:] == "mphANG_0479861" for g in calls[r]):
                filtered_calls[r] = calls[r]
        graph = GeneMerGraph(filtered_calls, 3, positions)
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "mphANG_0479861", 1)
        # assertion
        self.assertEqual(len(f_paths), 2)
        for k in f_paths:
            self.assertTrue(len(k) == 16 or len(k) == 43)
            self.assertTrue(len(f_paths[k]) == 2 or len(f_paths[k]) == 7)

    def test___find_full_paths_diverging_paths_at_terminals(self):
        # setup
        genes1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
            "-gene12",
            "+gene13",
            "-gene14",
            "+gene15",
        ]
        genes2 = [
            "+gene16",
            "-gene17",
            "+gene18",
            "-gene19",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
            "-gene10",
            "+gene11",
            "-gene20",
            "+gene21",
            "-gene22",
            "+gene23",
        ]
        graph = GeneMerGraph(
            {
                "read1": genes1,
                "read2": genes1,
                "read3": genes1,
                "read4": genes2,
                "read5": genes2,
                "read6": genes2,
            },
            3,
        )
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene7", 1)
        # assertion
        self.assertEqual(len(f_paths), 2)
        for k in f_paths:
            self.assertEqual(len(k), 15)
            self.assertEqual(len(f_paths[k]), 3)

    def test___find_full_paths_edge_case_one(self):
        # setup
        import json

        with open("tests/complex_gene_calls_five.json") as i:
            calls = json.load(i)
        with open("tests/complex_gene_positions_five.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "blaCTXM110NG_0489052", 1)
        self.assertEqual(len(f_paths), 1)
        for f in f_paths:
            self.assertEqual(len(f), 1)
            self.assertEqual(len(f_paths[f]), 44)

    def test___find_full_paths_variant(self):
        # setup
        import json

        with open("tests/complex_gene_calls_six.json") as i:
            calls = json.load(i)
        with open("tests/complex_gene_positions_six.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "blaTEM239NG_0766451", 1)
        # assertion
        self.assertEqual(len(f_paths), 2)
        for f in f_paths:
            self.assertEqual(len(f), 2)
            self.assertTrue(len(f_paths[f]) == 30 or len(f_paths[f]) == 9)

    def test___find_full_paths_multi_variants(self):
        # setup
        import json

        with open("/home/daniel/Documents/GitHub/Amira/test_data/WM/SRR30896056/corrected_gene_calls_after_filtering.json") as i:
            calls = json.load(i)
        with open("/home/daniel/Documents/GitHub/Amira/test_data/WM/SRR30896056/corrected_gene_positions_after_filtering.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "blaTEM239NG_0766451", 10)
        # assertion
        self.assertEqual(len(f_paths), 2)
        for f in f_paths:
            self.assertTrue(len(f) == 8 or len(f) == 19)
            self.assertTrue(len(f_paths[f]) == 28 or len(f_paths[f]) == 49)

    def test___find_full_paths_multi_tandem(self):
        # setup
        read1 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
        ]
        read2 = [
            "+gene1",
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "+gene5",
            "-gene6",
            "+gene7",
            "-gene8",
            "+gene9",
        ]
        read3 = [
            "-gene2",
            "+gene3",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "+gene5",
            "-gene6",
        ]
        read4 = [
            "+gene3",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "-gene4",
            "+gene5",
        ]
        read5 = ["+gene1", "-gene2", "+gene3", "-gene4", "-gene4", "-gene4"]
        read6 = ["-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "+gene5", "-gene6"]
        read7 = ["-gene4", "-gene4", "-gene4"]
        read8 = ["+gene3", "-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "+gene5"]
        read9 = ["-gene10", "+gene4", "-gene11"]
        read10 = ["-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "-gene4", "-gene4"]
        calls = {
            "read1": read1,
            "read2": read2,
            "read3": GeneMerGraph({}, 3).reverse_list_of_genes(read3),
            "read4": GeneMerGraph({}, 3).reverse_list_of_genes(read4),
            "read5": read5,
            "read6": read6,
            "read7": read7,
            "read8": read8,
            "read9": read9,
            "read10": read10,
            "read11": read10,
            "read12": read10,
            "read13": read10,
        }
        graph = GeneMerGraph(calls, 3)
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "gene4", 2)
        self.assertEqual(len(f_paths), 2)
        for f in f_paths:
            self.assertTrue(len(f) == 8 or len(f) == 9)
            self.assertTrue(len(f_paths[f]) == 3 or len(f_paths[f]) == 6)

    def test___find_full_paths_end_with_self_loop(self):
        # setup
        import json

        with open("tests/complex_gene_calls_seven.json") as i:
            calls = json.load(i)
        with open("tests/complex_gene_positions_seven.json") as i:
            positions = json.load(i)
        graph = GeneMerGraph(calls, 3, positions)
        calls = graph.get_reads().copy()
        rc_reads = {}
        for r in calls:
            rc_reads[r + "_reverse"] = graph.reverse_list_of_genes(calls[r])
        calls.update(rc_reads)
        tree = Tree(calls)
        # execution
        f_paths, f_path_coverages = graph.get_paths_for_gene(tree, "blaIMI9NG_0491711", 5)
        # assertion
        self.assertEqual(len(f_paths), 1)
        for f in f_paths:
            self.assertEqual(len(f), 5)
            self.assertEqual(len(f_paths[f]), 9)

    def test___cluster_adjacent_paths(self):
        # setup
        adjacent_paths = {
            (0, 1, 2, 3, 4): {"read1"},
            (1, 2, 3, 4): {"read2"},
            (2, 3, 4): {"read3"},
            (5, 6, 3, 4): {"read4", "read5"},
            (6, 3, 4): {"read6"},
            (5, 3, 2, 4): {"read7"},
            (3, 4): {"read8"},
        }
        # execution
        actual_clusters = cluster_adjacent_paths(adjacent_paths)
        # assertion
        self.assertEqual(len(actual_clusters), 3)
        self.assertTrue((2, 3, 4) in actual_clusters)
        self.assertTrue((6, 3, 4) in actual_clusters)
        self.assertTrue((5, 3, 2, 4) in actual_clusters)

    def test___cluster_adjacent_paths_overlapping(self):
        # setup
        adjacent_paths = {
            (0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12): {"read1", "read2"},
            (5, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12): {"read3", "read4", "read5"},
            (5, 6, 2, 3, 4, 7, 8, 9, 10, 11, 12): {"read6", "read7"},
            (2, 3, 4, 7, 8, 9, 10): {"read8"},
        }
        # execution
        actual_clusters = cluster_adjacent_paths(adjacent_paths)
        # assertion
        self.assertEqual(len(actual_clusters), 3)
        self.assertTrue((0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12) in actual_clusters)
        self.assertTrue((5, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12) in actual_clusters)
        self.assertTrue((5, 6, 2, 3, 4, 7, 8, 9, 10, 11, 12) in actual_clusters)
