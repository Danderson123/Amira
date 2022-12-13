import sys
sys.path.insert(0, "../amira_prototype")

from construct_node import Node
from construct_read import Read
from construct_edge import Edge

def test_node_increment_coverage():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
    read = Read("read1",
                genes)
    geneMers = [x for x in read.get_geneMers(3)]
    for g in geneMers:
        node = Node(g)
        for i in range(5):
            coverage = node.increment_coverage()
            assert coverage == i + 2

def test_node_add_read():
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
            node.add_read(read)
        assert all(r in node.get_reads() for r in readList), "A read is missing from the node"
        assert all([x for x in node.get_reads()].count(r) == 1 for r in readList), "A read occurs more than once in the read list"

def test_node_remove_read():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
    read1 = Read("read1",
                genes)
    read2 = Read("read2",
                genes)
    read3 = Read("read3",
                genes)
    geneMers = [x for x in read1.get_geneMers(3)]
    # test removing a read that is in the list
    for g in geneMers:
        node = Node(g)
        readList = [read1, read1, read2, read3, read3]
        for read in readList:
            node.add_read(read)
        for read in readList:
            node.remove_read(read)
            assert not read in node.get_reads()
    # test removing a read that is not in the list
    read4 = Read("read4",
                genes)
    for g in geneMers:
        node = Node(g)
        readList = [read1, read1, read2, read3, read3]
        for read in readList:
            node.add_read(read)
    node.remove_read(read4)
    assert all(r in node.get_reads() for r in readList), "A read was incorrectly removed from the Node"
    assert all([x for x in node.get_reads()].count(r) == 1 for r in readList), "A read occurs more than once in the read list"

def test_node_add_forward_edge():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
    read1 = Read("read1",
                genes)
    geneMers = [x for x in read1.get_geneMers(3)]
    # test the edge has not been added before
    nodes = []
    edges = []
    sourceNodes = []
    targetNodes = []
    targetNodeDirections = []
    for n in range(len(geneMers) - 1):
        sourceNode = Node(geneMers[n])
        targetNode = Node(geneMers[n+1])
        edge = Edge(sourceNode,
                    targetNode,
                    1,
                    geneMers[n+1].get_geneMerDirection())
        sourceNode.add_forward_edge(edge)
        edges.append(edge)
        sourceNodes.append(sourceNode)
        targetNodes.append(targetNode)
        targetNodeDirections.append(geneMers[n+1].get_geneMerDirection())
    for n in range(len(sourceNodes)):
        forwardEdge = sourceNodes[n].get_forward_edges()[0]
        assert forwardEdge.get_sourceNode().__eq__(sourceNodes[n]) and sourceNodes[n].__eq__(forwardEdge.get_sourceNode())
        assert forwardEdge.get_targetNode().__eq__(targetNodes[n]) and targetNodes[n].__eq__(forwardEdge.get_targetNode())
        assert forwardEdge.get_sourceNodeDirection() == 1
        assert forwardEdge.get_targetNodeDirection() == targetNodeDirections[n]
    # test the edge has been added before
    for e in range(len(edges)):
        returnedEdge = sourceNodes[e].add_forward_edge(edges[e])
        assert returnedEdge.get_sourceNode().__eq__(sourceNodes[e]) and sourceNodes[e].__eq__(returnedEdge.get_sourceNode())
        assert returnedEdge.get_targetNode().__eq__(targetNodes[e]) and targetNodes[e].__eq__(returnedEdge.get_targetNode())
        assert returnedEdge.get_sourceNodeDirection() == 1
        assert returnedEdge.get_targetNodeDirection() == targetNodeDirections[e]
        # ensure the coverage has been increased by 1
        assert returnedEdge.edgeCoverage == 2

def test_node_add_backward_edge():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
    read1 = Read("read1",
                genes)
    geneMers = [x for x in read1.get_geneMers(3)]
    # test the edge has not been added before
    edges = []
    sourceNodes = []
    targetNodes = []
    targetNodeDirections = []
    for n in range(len(geneMers) - 1):
        sourceNode = Node(geneMers[n])
        targetNode = Node(geneMers[n+1])
        edge = Edge(sourceNode,
                    targetNode,
                    -1,
                    geneMers[n+1].get_geneMerDirection())
        sourceNode.add_backward_edge(edge)
        edges.append(edge)
        sourceNodes.append(sourceNode)
        targetNodes.append(targetNode)
        targetNodeDirections.append(geneMers[n+1].get_geneMerDirection())
    for n in range(len(sourceNodes)):
        backwardEdge = sourceNodes[n].get_backward_edges()[0]
        assert backwardEdge.get_sourceNode().__eq__(sourceNodes[n]) and sourceNodes[n].__eq__(backwardEdge.get_sourceNode())
        assert backwardEdge.get_targetNode().__eq__(targetNodes[n]) and targetNodes[n].__eq__(backwardEdge.get_targetNode())
        assert backwardEdge.get_sourceNodeDirection() == -1
        assert backwardEdge.get_targetNodeDirection() == targetNodeDirections[n]
    # test the edge has been added before
    for e in range(len(edges)):
        returnedEdge = sourceNodes[e].add_backward_edge(edges[e])
        assert returnedEdge.get_sourceNode().__eq__(sourceNodes[e]) and sourceNodes[e].__eq__(returnedEdge.get_sourceNode())
        assert returnedEdge.get_targetNode().__eq__(targetNodes[e]) and targetNodes[e].__eq__(returnedEdge.get_targetNode())
        assert returnedEdge.get_sourceNodeDirection() == -1
        assert returnedEdge.get_targetNodeDirection() == targetNodeDirections[e]
        # ensure the coverage has been increased by 1
        assert returnedEdge.edgeCoverage == 2

def test_node_add_edge():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
    read1 = Read("read1",
                genes)
    geneMers = [x for x in read1.get_geneMers(3)]
    # test adding a forward edge that has not been added before
    edges = []
    sourceNodes = []
    targetNodes = []
    targetNodeDirections = []
    for n in range(len(geneMers) - 1):
        sourceNode = Node(geneMers[n])
        targetNode = Node(geneMers[n+1])
        edge = Edge(sourceNode,
                    targetNode,
                    1,
                    geneMers[n+1].get_geneMerDirection())
        sourceNode.add_edge(edge)
        edges.append(edge)
        sourceNodes.append(sourceNode)
        targetNodes.append(targetNode)
        targetNodeDirections.append(geneMers[n+1].get_geneMerDirection())
    for n in range(len(sourceNodes)):
        forwardEdge = sourceNodes[n].get_forward_edges()[0]
        assert forwardEdge.get_sourceNode().__eq__(sourceNodes[n]) and sourceNodes[n].__eq__(forwardEdge.get_sourceNode())
        assert forwardEdge.get_targetNode().__eq__(targetNodes[n]) and targetNodes[n].__eq__(forwardEdge.get_targetNode())
        assert forwardEdge.get_sourceNodeDirection() == 1
        assert forwardEdge.get_targetNodeDirection() == targetNodeDirections[n]
    # test adding a forward edge that has been added before
    for e in range(len(edges)):
        returnedEdge = sourceNodes[e].add_edge(edges[e])
        assert returnedEdge.get_sourceNode().__eq__(sourceNodes[e]) and sourceNodes[e].__eq__(returnedEdge.get_sourceNode())
        assert returnedEdge.get_targetNode().__eq__(targetNodes[e]) and targetNodes[e].__eq__(returnedEdge.get_targetNode())
        assert returnedEdge.get_sourceNodeDirection() == 1
        assert returnedEdge.get_targetNodeDirection() == targetNodeDirections[e]
        # ensure the coverage has been increased by 1
        assert returnedEdge.edgeCoverage == 2
    # test adding a backward edge that has not been added before
    edges = []
    sourceNodes = []
    targetNodes = []
    targetNodeDirections = []
    for n in range(len(geneMers) - 1):
        sourceNode = Node(geneMers[n])
        targetNode = Node(geneMers[n+1])
        edge = Edge(sourceNode,
                    targetNode,
                    -1,
                    geneMers[n+1].get_geneMerDirection())
        sourceNode.add_edge(edge)
        edges.append(edge)
        sourceNodes.append(sourceNode)
        targetNodes.append(targetNode)
        targetNodeDirections.append(geneMers[n+1].get_geneMerDirection())
    for n in range(len(sourceNodes)):
        backwardEdge = sourceNodes[n].get_backward_edges()[0]
        assert backwardEdge.get_sourceNode().__eq__(sourceNodes[n]) and sourceNodes[n].__eq__(backwardEdge.get_sourceNode())
        assert backwardEdge.get_targetNode().__eq__(targetNodes[n]) and targetNodes[n].__eq__(backwardEdge.get_targetNode())
        assert backwardEdge.get_sourceNodeDirection() == -1
        assert backwardEdge.get_targetNodeDirection() == targetNodeDirections[n]
    # test the edge has been added before
    for e in range(len(edges)):
        returnedEdge = sourceNodes[e].add_edge(edges[e])
        assert returnedEdge.get_sourceNode().__eq__(sourceNodes[e]) and sourceNodes[e].__eq__(returnedEdge.get_sourceNode())
        assert returnedEdge.get_targetNode().__eq__(targetNodes[e]) and targetNodes[e].__eq__(returnedEdge.get_targetNode())
        assert returnedEdge.get_sourceNodeDirection() == -1
        assert returnedEdge.get_targetNodeDirection() == targetNodeDirections[e]
        # ensure the coverage has been increased by 1
        assert returnedEdge.edgeCoverage == 2

sys.stderr.write("Testing construct_node: Node.increment_coverage\n")
test_node_increment_coverage()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_node: Node.add_read\n")
test_node_add_read()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_node: Node.remove_read\n")
test_node_remove_read()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_node: Node.add_forward_edge\n")
test_node_add_forward_edge()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_node: Node.add_backward_edge\n")
test_node_add_backward_edge()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_node: Node.add_edge\n")
test_node_add_edge()
sys.stderr.write("Test passed\n")