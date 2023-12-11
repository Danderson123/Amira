import unittest

from amira_prototype.construct_edge import Edge
from amira_prototype.construct_node import Node
from amira_prototype.construct_read import Read


class TestEdgeConstructor(unittest.TestCase):
    def test___init_Read(self):
        # setup
        sourceNode = 0
        targetNode = 1
        sourceNodeDirection = 1
        targetNodeDirection = -1
        # execution
        actual_edge = Edge(sourceNode, targetNode, sourceNodeDirection, targetNodeDirection)
        actual_edge_sourceNode = actual_edge.get_sourceNode()
        actual_edge_targetNode = actual_edge.get_targetNode()
        actual_sourceNodeDirection = actual_edge.get_sourceNodeDirection()
        actual_targetNodeDirection = actual_edge.get_targetNodeDirection()
        actual_edge_coverage = actual_edge.get_edge_coverage()
        # assertion
        expected_edge_sourceNode = 0
        expected_edge_targetNode = 1
        expected_sourceNodeDirection = 1
        expected_targetNodeDirection = -1
        expected_edge_coverage = 0
        self.assertEqual(actual_edge_sourceNode, expected_edge_sourceNode)
        self.assertEqual(actual_edge_targetNode, expected_edge_targetNode)
        self.assertEqual(actual_sourceNodeDirection, expected_sourceNodeDirection)
        self.assertEqual(actual_targetNodeDirection, expected_targetNodeDirection)
        self.assertEqual(actual_edge_coverage, expected_edge_coverage)

    def test___increment_edge_coverage(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1", genes1)
        nodes = [Node(g) for g in read1.get_geneMers(3)]
        directions = [g.get_geneMerDirection() for g in read1.get_geneMers(3)]
        for n in range(len(nodes) - 1):
            edge = Edge(nodes[n], nodes[n + 1], directions[n], directions[n + 1])
            for i in range(5):
                # execution
                actual_coverage = edge.increment_edge_coverage()
                # assertion
                expected_coverage = i + 1
                self.assertEqual(actual_coverage, expected_coverage)

    def test___eq_edge_reads_identical(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        nodes1 = [Node(g) for g in read1.get_geneMers(3)]
        directions1 = [g.get_geneMerDirection() for g in read1.get_geneMers(3)]
        nodes2 = [Node(g) for g in read2.get_geneMers(3)]
        directions2 = [g.get_geneMerDirection() for g in read2.get_geneMers(3)]
        edges1 = []
        for n in range(len(nodes1) - 1):
            edge = Edge(nodes1[n], nodes1[n + 1], directions1[n], directions1[n + 1])
            edges1.append(edge)
        edges2 = []
        for n in range(len(nodes2) - 1):
            edge = Edge(nodes2[n], nodes2[n + 1], directions2[n], directions2[n + 1])
            edges2.append(edge)
        # assertion
        self.assertTrue(
            all(
                edges2[i].__eq__(edges1[i]) and edges1[i].__eq__(edges2[i])
                for i in range(len(edges1))
            )
        )

    def test___eq_edge_different(self):
        # setup
        genes1 = ["+gene1", "-gene2"]
        genes2 = ["-gene2", "+gene3"]
        genes3 = ["-gene2", "+gene4"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        read3 = Read("read3", genes3)
        node1 = [Node(g) for g in read1.get_geneMers(2)][0]
        direction1 = [g.get_geneMerDirection() for g in read1.get_geneMers(2)][0]
        node2 = [Node(g) for g in read2.get_geneMers(2)][0]
        direction2 = [g.get_geneMerDirection() for g in read2.get_geneMers(2)][0]
        node3 = [Node(g) for g in read3.get_geneMers(2)][0]
        direction3 = [g.get_geneMerDirection() for g in read3.get_geneMers(2)][0]
        edge1 = Edge(node1, node2, direction1, direction2)
        edge2 = Edge(node2, node3, direction2, direction3)
        # assertion
        self.assertFalse(edge1.__eq__(edge2))
        self.assertFalse(edge2.__eq__(edge1))

    def test___eq_edge_reversed(self):
        # setup
        genes1 = ["+gene1", "-gene2"]
        genes2 = ["-gene2", "+gene3"]
        genes3 = ["+gene2", "-gene1"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        read3 = Read("read3", genes3)
        node1 = [Node(g) for g in read1.get_geneMers(2)][0]
        direction1 = [g.get_geneMerDirection() for g in read1.get_geneMers(2)][0]
        node2 = [Node(g) for g in read2.get_geneMers(2)][0]
        direction2 = [g.get_geneMerDirection() for g in read2.get_geneMers(2)][0]
        node3 = [Node(g) for g in read3.get_geneMers(2)][0]
        direction3 = [g.get_geneMerDirection() for g in read3.get_geneMers(2)][0]
        edge1 = Edge(node1, node2, direction1, direction2)
        edge2 = Edge(node2, node3, direction2, direction3)
        # assertion
        self.assertTrue(edge1.__eq__(edge2))
        self.assertTrue(edge2.__eq__(edge1))

    def test___hash_edge_same_edges(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        nodes1 = [Node(g) for g in read1.get_geneMers(3)]
        directions1 = [g.get_geneMerDirection() for g in read1.get_geneMers(3)]
        nodes2 = [Node(g) for g in read2.get_geneMers(3)]
        directions2 = [g.get_geneMerDirection() for g in read2.get_geneMers(3)]
        edgeHashes1 = []
        for n in range(len(nodes1) - 1):
            edge = Edge(nodes1[n], nodes1[n + 1], directions1[n], directions1[n + 1])
            edgeHashes1.append(edge.__hash__())
        edgesHashes2 = []
        for n in range(len(nodes2) - 1):
            edge = Edge(nodes2[n], nodes2[n + 1], directions2[n], directions2[n + 1])
            edgesHashes2.append(edge.__hash__())
        # assertion
        self.assertTrue(all(edgesHashes2[i] == edgeHashes1[i] for i in range(len(edgeHashes1))))
        self.assertTrue(all(edgeHashes1[i] == edgesHashes2[i] for i in range(len(edgesHashes2))))

    def test___hash_edge_same_complement_edges(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["+gene4", "-gene3", "+gene2", "-gene1"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        nodes1 = [Node(g) for g in read1.get_geneMers(3)]
        directions1 = [g.get_geneMerDirection() for g in read1.get_geneMers(3)]
        nodes2 = [Node(g) for g in read2.get_geneMers(3)]
        directions2 = [g.get_geneMerDirection() for g in read2.get_geneMers(3)]
        edgeHashes1 = []
        for n in range(len(nodes1) - 1):
            edge = Edge(nodes1[n], nodes1[n + 1], directions1[n], directions1[n + 1])
            edgeHashes1.append(edge.__hash__())
        edgesHashes2 = []
        for n in range(len(nodes2) - 1):
            edge = Edge(nodes2[n + 1], nodes2[n], directions2[n + 1], directions2[n])
            edgesHashes2.append(edge.__hash__())
        # assertion
        self.assertTrue(
            all(
                nodes1[i].__hash__() == list(reversed(nodes2))[i].__hash__()
                for i in range(len(nodes1))
            )
        )
        self.assertTrue(
            all(
                nodes2[i].__hash__() == list(reversed(nodes1))[i].__hash__()
                for i in range(len(nodes2))
            )
        )
        self.assertTrue(all(edgesHashes2[i] == edgeHashes1[i] for i in range(len(edgeHashes1))))
        self.assertTrue(all(edgeHashes1[i] == edgesHashes2[i] for i in range(len(edgesHashes2))))

    def test___hash_edge_different_edges(self):
        # setup
        genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
        genes2 = ["+gene4", "-gene5", "+gene6", "-gene7"]
        read1 = Read("read1", genes1)
        read2 = Read("read2", genes2)
        nodes1 = [Node(g) for g in read1.get_geneMers(3)]
        directions1 = [g.get_geneMerDirection() for g in read1.get_geneMers(3)]
        nodes2 = [Node(g) for g in read2.get_geneMers(3)]
        directions2 = [g.get_geneMerDirection() for g in read2.get_geneMers(3)]
        edgeHashes1 = []
        for n in range(len(nodes1) - 1):
            edge = Edge(nodes1[n], nodes1[n + 1], directions1[n], directions1[n + 1])
            edgeHashes1.append(edge.__hash__())
        edgesHashes2 = []
        for n in range(len(nodes2) - 1):
            edge = Edge(nodes2[n], nodes2[n + 1], directions2[n], directions2[n + 1])
            edgesHashes2.append(edge.__hash__())
        # assertion
        self.assertTrue(all(edgesHashes2[i] != edgeHashes1[i] for i in range(len(edgeHashes1))))
        self.assertTrue(all(edgeHashes1[i] != edgesHashes2[i] for i in range(len(edgesHashes2))))
