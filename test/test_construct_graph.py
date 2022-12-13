import sys
sys.path.insert(0, "../amira_prototype")

from construct_graph import GeneMerGraph
from construct_read import Read

def test_graph_increment_nodeId():
    graph = GeneMerGraph([],
                        3,
                        1,
                        1)
    assert graph.currentNodeId == 0
    assert graph.increment_nodeID() == 1 and graph.currentNodeId == 1

def test_graph_determine_node_index():
    # test gene mer in the graph
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
    read1 = Read("read1",
                genes)
    geneMers = [x for x in read1.get_geneMers(2)]
    graph = GeneMerGraph([],
                        3,
                        1,
                        1)
    graph.nodeHashes = [g.__hash__() for g in geneMers]
    assert all(graph.determine_node_index(geneMers[i]) == i for i in range(len(geneMers)))
    # test gene-mer not in graph
    try:
        graph.determine_node_index(1)
    except Exception as e:
        assert isinstance(e, AssertionError)

def test_graph_add_node():
    # add a node to an empty graph
    genes = ["+gene1", "-gene2", "+gene3"]
    read1 = Read("read1",
                genes)
    geneMer = [x for x in read1.get_geneMers(3)][0]
    graph = GeneMerGraph([],
                        3,
                        1,
                        1)
    returnedNode = graph.add_node(geneMer,
                                "read1")
    # ensure node and gene mer hashes are the same
    assert returnedNode.__hash__() == geneMer.__hash__()
    # ensure the node id assignment is correct
    assert returnedNode.get_nodeId() == 0 and graph.currentNodeId == 1
    # ensure the returned node is identical to the node added to the graph
    assert returnedNode.__hash__() == graph.graph[0].__hash__()
    # ensure the read id is added to the Node
    assert [x for x in returnedNode.get_reads()][0] == "read1"
    # ensure the hash is in the list of nodeHashes exactly once
    assert graph.nodeHashes.count(returnedNode.__hash__()) == 1
    # add a different node to a non-empty graph
    genes2 = ["+gene4", "-gene3", "+gene2"]
    read2 = Read("read2",
                genes2)
    geneMer2 = [x for x in read2.get_geneMers(3)][0]
    returnedNode2 = graph.add_node(geneMer2,
                                "read2")
    # ensure node and gene mer hashes are the same
    assert returnedNode2.__hash__() == geneMer2.__hash__()
    # ensure the node id assignment is correct
    assert returnedNode2.get_nodeId() == 1 and graph.currentNodeId == 2
    # ensure the returned node is identical to the node added to the graph
    assert returnedNode2.__hash__() == graph.graph[1].__hash__()
    # ensure the read id is added to the Node
    assert [x for x in returnedNode2.get_reads()][0] == "read2"
    # ensure the hash is in the list of nodeHashes exactly once
    assert graph.nodeHashes.count(returnedNode2.__hash__()) == 1
    # add the same node to a non-empty graph
    genes3 = ["-gene3", "+gene2", "-gene1"]
    read3 = Read("read3",
                genes3)
    geneMer3 = [x for x in read3.get_geneMers(3)][0]
    returnedNode3 = graph.add_node(geneMer3,
                                "read3")
    # ensure node and gene mer hashes are the same
    assert returnedNode3.__hash__() == geneMer3.__hash__() and returnedNode3.__hash__() == geneMer.__hash__()
    # ensure the node id assignment is correct
    assert returnedNode3.get_nodeId() == 0 and graph.currentNodeId == 2
    # ensure the returned node is identical to the node added to the graph
    assert returnedNode3.__hash__() == graph.graph[0].__hash__()
    # ensure the read id is added to the Node
    assert [x for x in returnedNode3.get_reads()] == ["read1", "read3"] or [x for x in returnedNode3.get_reads()] == ["read3", "read1"]
    # ensure the hash is in the list of nodeHashes exactly once
    assert graph.nodeHashes.count(returnedNode3.__hash__()) == 1

def test_graph_get_node():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
    read1 = Read("read1",
                genes)
    geneMers = [x for x in read1.get_geneMers(3)]
    graph = GeneMerGraph([],
                        3,
                        1,
                        1)
    for gmer in geneMers:
        graph.add_node(gmer,
                    "read1")
    #test gene-mer is in the graph
    ids = [0, 1, 2, 3, 4, 5, 0]
    for gmer in range(len(geneMers)):
        assert graph.get_node(geneMers[gmer]).get_nodeId() == ids[gmer]
    # test gene-mer is not in the graph
    genes2 = ["+gene7", "-gene8", "+gene9"]
    read2 = Read("read2",
                genes2)
    geneMer2 = [x for x in read2.get_geneMers(3)][0]
    try:
        graph.get_node(geneMer2)
    except Exception as e:
        assert isinstance(e, AssertionError)

def test_graph_get_nodes_containing():
    graph = GeneMerGraph({"read1": ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]},
                        3,
                        1,
                        1)
    allNodes = graph.graph
    # select two random genes and check all nodes contain the genes
    selectedGenes = ["gene2", "gene6"]
    selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
    nodeGenes = [[g.get_name() for g in n.get_geneMer().get_canonical_geneMer()] for n in selectedNodes]
    assert all(any(g in ng for g in selectedGenes) for ng in nodeGenes)
    # confirm all nodes are returned when all genes are specified
    selectedGenes = ["gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene3", "gene2", "gene1"]
    selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
    nodeGenes = [[g.get_name() for g in n.get_geneMer().get_canonical_geneMer()] for n in selectedNodes]
    assert all(any(g in ng for g in selectedGenes) for ng in nodeGenes)
    # confirm no nodes are returned when no genes in the graph are specified
    selectedGenes = ["gene10"]
    selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
    nodeGenes = [[g.get_name() for g in n.get_geneMer().get_canonical_geneMer()] for n in selectedNodes]
    assert nodeGenes == []
    selectedGenes = []
    selectedNodes = [n for n in graph.get_nodes_containing(selectedGenes)]
    nodeGenes = [[g.get_name() for g in n.get_geneMer().get_canonical_geneMer()] for n in selectedNodes]
    assert nodeGenes == []
    # confirm there is an assertion error when the specified genes include strand information
    selectedGenes = ["gene2", "+gene6"]
    try:
        [n for n in graph.get_nodes_containing(selectedGenes)]
    except Exception as e:
        assert isinstance(e, AssertionError)

def test_graph_add_edge():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6", "-gene3", "+gene2", "-gene1"]
    read1 = Read("read1",
                genes)
    geneMers = [x for x in read1.get_geneMers(3)]
    graph = GeneMerGraph([],
                        3,
                        1,
                        1)
    for gmer in geneMers:
        graph.add_node(gmer,
                    "read1")
    allGraphNodes = [n for n in graph.all_nodes()]
    for n in range(len(allGraphNodes) - 1):
        startSourceCoverage = allGraphNodes[n].get_coverage()
        startTargetCoverage = allGraphNodes[n+1].get_coverage()
        graph.add_edge(allGraphNodes[n],
                    allGraphNodes[n + 1])
        # ensure the node coverage is not modified as the nodes already exist in the graph
        assert allGraphNodes[n].get_coverage() == startSourceCoverage, "Source node coverage changed when an edge was added"
        assert allGraphNodes[n+1].get_coverage() == startTargetCoverage, "Target node coverage changed when an edge was added"

sys.stderr.write("Testing construct_graph: Graph.increment_nodeId\n")
test_graph_increment_nodeId()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_graph: Graph.determine_node_index\n")
test_graph_determine_node_index()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_graph: Graph.add_node\n")
test_graph_add_node()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_graph: Graph.get_node\n")
test_graph_get_node()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_graph: Graph.get_nodes_containing\n")
test_graph_get_nodes_containing()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_graph: Graph.add_edge\n")
test_graph_add_edge()
sys.stderr.write("Test passed\n")