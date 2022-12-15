import unittest
import sys
sys.path.insert(0, "amira_prototype")

from construct_node import Node
from construct_read import Read

class TestNodeConstructor(unittest.TestCase):

    def test___init_Node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                genes)
        geneMer = [x for x in read1.get_geneMers(3)][0]
        # execution
        actual_node = Node(geneMer)
        actual_canonicalGeneMer = actual_node.get_canonical_geneMer()
        actual_geneMerHash = actual_node.__hash__()
        actual_nodeCoverage = actual_node.get_node_coverage()
        actual_listOfReads = actual_node.get_list_of_reads()
        actual_forwardEdgeHashes = actual_node.get_forward_edge_hashes()
        actual_backwardEdgeHashes = actual_node.get_backward_edge_hashes()
        # assertion
        expected_canonicalGeneMer = geneMer.get_canonical_geneMer()
        expected_geneMerHash = geneMer.__hash__()
        expected_nodeCoverage = 0
        expected_listOfReads = []
        expected_forwardEdgeHashes = []
        expected_backwardEdgeHashes = []
        self.assertEqual(actual_canonicalGeneMer, expected_canonicalGeneMer)
        self.assertEqual(actual_geneMerHash, expected_geneMerHash)
        self.assertEqual(actual_nodeCoverage, expected_nodeCoverage)
        self.assertEqual(actual_listOfReads, expected_listOfReads)
        self.assertEqual(actual_forwardEdgeHashes, expected_forwardEdgeHashes)
        self.assertEqual(actual_backwardEdgeHashes, expected_backwardEdgeHashes)

    def test___node_increment_node_coverage(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        read1 = Read("read1",
                genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        # execution
        for g in geneMers:
            node = Node(g)
            # assertion
            self.assertEqual(node.get_node_coverage(), 0)
            for i in range(5):
                actual_coverage = node.increment_node_coverage()
                # assertion
                expected_coverage = i + 1
                self.assertEqual(actual_coverage, expected_coverage)

    def test___add_read_to_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        read1 = Read("read1",
                    genes)
        read2 = Read("read2",
                    genes)
        read3 = Read("read3",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        for g in geneMers:
            node = Node(g)
            readList = [read1, read1, read2, read3, read3]
            for read in readList:
                # execution
                node.add_read(read)
            # assertion
            self.assertTrue(all(r in node.get_reads() for r in readList))
            self.assertTrue(all([x for x in node.get_reads()].count(r) == 1 for r in readList))

    def test___remove_read_in_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        read1 = Read("read1",
                    genes)
        read2 = Read("read2",
                    genes)
        read3 = Read("read3",
                    genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        for g in geneMers:
            node = Node(g)
            readList = [read1, read2, read3]
            for read in readList:
                node.add_read(read)
            # execution
            for read in readList:
                node.remove_read(read)
                # assertion
                self.assertFalse(read in node.get_reads())

    def test___remove_read_not_in_node(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
        read1 = Read("read1",
                    genes)
        read2 = Read("read2",
                    genes)
        read3 = Read("read3",
                    genes)
        read4 = Read("read4",
                genes)
        geneMers = [x for x in read1.get_geneMers(3)]
        readList = [read1, read1, read2, read3, read3]
        for g in geneMers:
            node = Node(g)
            for read in readList:
                node.add_read(read)
        # assertion
        self.assertRaises(AssertionError, node.remove_read, read4)
        self.assertTrue(all(r in node.get_reads() for r in readList))
        self.assertTrue(all([x for x in node.get_reads()].count(r) == 1 for r in readList))

    def test_add_forward_edge_hash_to_empty(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        # execution
        actual_updated_node = node.add_forward_edge_hash(12345)
        actual_list_of_forward_hashes = actual_updated_node.get_forward_edge_hashes()
        # assertion
        expected_list_of_forward_hashes = [12345]
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_list_of_forward_hashes, expected_list_of_forward_hashes)

    def test_add_forward_edge_hash_to_non_empty(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        node.add_forward_edge_hash(12345)
        # execution
        actual_updated_node = node.add_forward_edge_hash(56789)
        actual_list_of_forward_hashes = actual_updated_node.get_forward_edge_hashes()
        # assertion
        expected_list_of_forward_hashes = [12345, 56789]
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_list_of_forward_hashes, expected_list_of_forward_hashes)

    def test_add_forward_duplicate_forward_edge_hash(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        node.add_forward_edge_hash(12345)
        # execution
        actual_updated_node = node.add_forward_edge_hash(12345)
        actual_list_of_forward_hashes = actual_updated_node.get_forward_edge_hashes()
        # assertion
        expected_list_of_forward_hashes = [12345]
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_list_of_forward_hashes, expected_list_of_forward_hashes)

    def test_add_backward_edge_hash_to_empty(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        # execution
        actual_updated_node = node.add_backward_edge_hash(12345)
        actual_list_of_backward_hashes = actual_updated_node.get_backward_edge_hashes()
        # assertion
        expected_list_of_backward_hashes = [12345]
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_list_of_backward_hashes, expected_list_of_backward_hashes)

    def test_add_backward_edge_hash_to_non_empty(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        node.add_backward_edge_hash(12345)
        # execution
        actual_updated_node = node.add_backward_edge_hash(56789)
        actual_list_of_backward_hashes = actual_updated_node.get_backward_edge_hashes()
        # assertion
        expected_list_of_backward_hashes = [12345, 56789]
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_list_of_backward_hashes, expected_list_of_backward_hashes)

    def test_add_backward_duplicate_forward_edge_hash(self):
        # setup
        genes = ["+gene1", "-gene2", "+gene3"]
        read1 = Read("read1",
                    genes)
        node = [Node(x) for x in read1.get_geneMers(3)][0]
        node.add_backward_edge_hash(12345)
        # execution
        actual_updated_node = node.add_backward_edge_hash(12345)
        actual_list_of_backward_hashes = actual_updated_node.get_backward_edge_hashes()
        # assertion
        expected_list_of_backward_hashes = [12345]
        self.assertEqual(actual_updated_node, node)
        self.assertEqual(actual_list_of_backward_hashes, expected_list_of_backward_hashes)
