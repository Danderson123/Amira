import sys
sys.path.insert(0, "../amira_prototype")

from construct_edge import Edge
from construct_node import Node
from construct_read import Read

def test_node_increment_coverage():
    genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
    read1 = Read("read1",
                genes1)
    nodes = []
    for g in [x for x in read1.get_geneMers(3)]:
        node = Node(g, 1)
        nodes.append(node)
    for n in range(len(nodes) - 1):
        edge = Edge(nodes[n],
                    nodes[n+1])
        for i in range(5):
            coverage = edge.increment_coverage()
            assert coverage == i + 1

def get_two_read_edges(read1,
                    read2):
    nodes1 = []
    for g in [x for x in read1.get_geneMers(3)]:
        node = Node(g, 1)
        nodes1.append(node)
    nodes2 = []
    for g in [x for x in read2.get_geneMers(3)]:
        node = Node(g, 1)
        nodes2.append(node)
    edges1 = []
    for n in range(len(nodes1) - 1):
        edge = Edge(nodes1[n],
                    nodes1[n+1])
        edges1.append(edge.__hash__())
    edges2 = []
    for n in range(len(nodes2) - 1):
        edge = Edge(nodes2[n],
                    nodes2[n+1])
        edges2.append(edge.__hash__())
    return edges1, edges2

def test_hash_edge():
    # check gene mers on both reads the same
    genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
    genes2 = ["+gene1", "-gene2", "+gene3", "-gene4"]
    read1 = Read("read1",
                genes1)
    read2 = Read("read2",
                genes2)
    edges1, edges2 = get_two_read_edges(read1,
                                        read2)
    assert all(edges2[i] == edges1[i] for i in range(len(edges1)))
    # check gene mer and rc the same
    genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
    genes2 = ["+gene4", "-gene3", "+gene2", "-gene1"]
    read1 = Read("read1",
                genes1)
    read2 = Read("read2",
                genes2)
    edges1, edges2 = get_two_read_edges(read1,
                                        read2)
    assert all(edges2[i] == edges1[i] for i in range(len(edges1)))
    # check different gene mers
    genes1 = ["+gene1", "-gene2", "+gene3", "-gene4"]
    genes2 = ["+gene4", "-gene5", "+gene6", "-gene7"]
    read1 = Read("read1",
                genes1)
    read2 = Read("read2",
                genes2)
    edges1, edges2 = get_two_read_edges(read1,
                                        read2)
    assert all(edges2[i] != edges1[i] for i in range(len(edges1)))

sys.stderr.write("Testing construct_edge: Edge.increment_coverage\n")
test_node_increment_coverage()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_edge: Edge.__hash__\n")
test_hash_edge()
sys.stderr.write("Test passed\n")