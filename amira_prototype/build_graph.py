class NodeData:
    def __init__(self,
                nodeId,
                nodeName,
                rc_nodeName,
                nodeCoverage,
                interestingNode,
                forwardEdgeList = [],
                reverseEdgeList = [],
                forwardEdgeIds = [],
                reverseEdgeIds = [],
                readIds = set()):
        self.nodeId = nodeId # the integer ID of this node (also the index in the graphData list)
        self.nodeName = nodeName # the geneMer we have chosen as the canonical geneMer
        self.rc_nodeName = rc_nodeName
        self.nodeCoverage = nodeCoverage # the number of times we have seen this geneMer in the data
        self.interestingNode = interestingNode # whether this node contains genes we are interested in
        self.forwardEdgeList = forwardEdgeList # edges in the forward direction from this node
        self.reverseEdgeList = reverseEdgeList # edges in the reverse direction from this node
        self.forwardEdgeIds = forwardEdgeIds # the list of target node IDs we see in the data to check if we have seen an edge before
        self.reverseEdgeIds = reverseEdgeIds # the list of target node IDs we see in the data to check if we have seen an edge before
        self.readIds = readIds

    def get_nodeId(self):
        return self.nodeId
    def get_nodeName(self):
        return self.nodeName
    def get_rc_nodeName(self):
        return self.rc_nodeName
    def get_nodeCoverage(self):
        return self.nodeCoverage
    def get_forwardEdgeList(self):
        return self.forwardEdgeList
    def get_reverseEdgeList(self):
        return self.reverseEdgeList
    def get_forwardEdgeIds(self):
        return self.forwardEdgeIds
    def get_reverseEdgeIds(self):
        return self.reverseEdgeIds
    def get_interestingNode(self):
        return self.interestingNode
    def get_readIds(self):
        return list(self.readIds)

class EdgeData:

    def __init__(self,
                targetNodeId,
                sourceDirection,
                targetDirection,
                edgeWeight):
        self.targetNodeId = targetNodeId # the target node ID (also the index in the graphData list)
        self.sourceDirection = sourceDirection # whether the target is the canonical kmer of the reverse complement of it
        self.targetDirection = targetDirection # direction that the edge leaves the source node
        self.edgeWeight = edgeWeight # the number of times we see the edge in the reads (in either orientation)

    def get_targetNodeId(self):
        return self.targetNodeId
    def get_targetDirection(self):
        return self.targetDirection
    def get_sourceDirection(self):
        return self.sourceDirection
    def get_edgeWeight(self):
        return self.edgeWeight
    def increment_weight(self):
        self.edgeWeight = self.edgeWeight + 1

def rc_geneMer(geneMer):
    """ returns the reverse complement of a list of integers gene identifiers. """
    # get the reverse complement of the geneMer
    reverse_geneMer = reversed(geneMer)
    # the strand of the gene is indicated by whether the integer gene ID is positive or negative
    reverse_geneMer = [str(int(g)*-1) for g in reverse_geneMer]
    return reverse_geneMer

def kmerise(geneMer,
            geneMerStrands):
    """ returns a ',' delimited string of gene integer identifiers and the strand of the geneMer. This is the canonical geneMer for the node. """
    geneMer = [str(geneMer[g] * geneMerStrands[g]) for g in range(len(geneMer))]
    reverse_geneMer = rc_geneMer(geneMer)
    # this selects the canonical geneMer by choosing the smallest of the seen geneMer and its reverse complement
    sorted_complements = sorted([geneMer, reverse_geneMer])
    canonical_geneMer = sorted_complements[0]
    rc_canonical_geneMer = sorted_complements[1]
    # if the canonical geneMer is not the seen geneMer then we say the strand is -1
    if not canonical_geneMer == geneMer:
        strand = -1
    else:
        strand = +1
    return ",".join(canonical_geneMer), ",".join(rc_canonical_geneMer), strand

def define_node(geneMerList,
                geneMerStrands,
                geneMerNameIdDict,
                geneMerIdNameDict,
                current_geneMerId,
                refined_genes,
                readId):
    """ returns a NodeData object for this node """
    # this decides what the canonical geneMer is
    canonical_geneMer, rc_canonical_geneMer, strand = kmerise(geneMerList,
                                                            geneMerStrands) # strand here means whether we see the canonical geneMer or reverse complement of it in this case
    # this assigns node IDs to each canonical geneMer
    if not canonical_geneMer in geneMerNameIdDict:
        # if we have not seen the canonical geneMer before we get a new nodeID
        geneMerNameIdDict[canonical_geneMer] = current_geneMerId
        geneMerIdNameDict[current_geneMerId] = canonical_geneMer
        nodeId = current_geneMerId
        current_geneMerId += 1
        # if we have seen the canonical geneMer before the nodeID is consistent
    else:
        nodeId = geneMerNameIdDict[canonical_geneMer]
    # this checks to see if any of the genes we are interested in are found in this node
    if any(g in geneMerList or -1 * g in geneMerList for g in refined_genes):
        interesting = True
    else:
        interesting = False
    this_node = NodeData(nodeId,
                        canonical_geneMer,
                        rc_canonical_geneMer,
                        1, # if this is a new node then the coverage will be 1
                        interesting, # whether or not this node contains a gene we are interested in
                        [], # this is empty as it will be populated with edge information later
                        [], # this is empty as it will be populated with edge information later ZAM added this
                        [], # this is empty as it will be populated with edge information later Zam added this
                        [], # this is empty as it will be populated with edge information later
                        set([readId]))
    return this_node, strand, geneMerNameIdDict, geneMerIdNameDict, current_geneMerId

def add_edge_to_graph(sourceNodeData,
                    sourceNodeDirection,
                    sourceToTargetEdge,
                    targetNodeData,
                    revCompEdgeFromTargetToSource,
                    graphData):
    """ this modifies the graphData list by index to add an edge between the source and target and then the target and the source """
    sourceNodeId = sourceNodeData.get_nodeId()
    targetNodeId = targetNodeData.get_nodeId()
    # ensure the source and target IDs are less than the max allowed length
    assert sourceNodeId < 50001
    assert targetNodeId < 50001
    # add an edge to the target node to the list of edges for the source node
    # we initialised the graphData list earlier so here we are changing the entry from None to the node information if this is the first time we see the geneMer
    if graphData[sourceNodeId] == None:
        graphData[sourceNodeId] = sourceNodeData
    else:
        # we see every geneMer as a sourceNode so add 1 to the coverage if we see it again
        graphData[sourceNodeId].nodeCoverage += 1
        # we add the readId to the list of reads for this node
        graphData[sourceNodeId].readIds.add(sourceNodeData.get_readIds()[0])
    # if we have not added an edge to the target node from the source node before then we add one here
    if sourceNodeDirection == 1:
        if not targetNodeId in graphData[sourceNodeId].get_forwardEdgeIds():
            graphData[sourceNodeId].forwardEdgeList.append(sourceToTargetEdge)
            graphData[sourceNodeId].forwardEdgeIds.append(targetNodeId)
        else:
            # increase the edge weight by 1 because we have seen this edge before
            graphData[sourceNodeId].forwardEdgeList[graphData[sourceNodeId].forwardEdgeIds.index(targetNodeId)].increment_weight()
    else:
        if not targetNodeId in graphData[sourceNodeId].get_reverseEdgeIds():
            graphData[sourceNodeId].reverseEdgeList.append(sourceToTargetEdge)
            graphData[sourceNodeId].reverseEdgeIds.append(targetNodeId)
        else:
            #get sourceToTargetEdge from graphData[sourceNodeId].reverseEdgeList and use increment_weight() on it
            graphData[sourceNodeId].reverseEdgeList[graphData[sourceNodeId].reverseEdgeIds.index(targetNodeId)].increment_weight()
    # add an edge to the source node to the list of edges for the target node
    # we initialised the graphData list earlier so here we are changing the entry from None to the node information if this is the first time we see the geneMer
    if graphData[targetNodeId] == None:
        graphData[targetNodeId] = targetNodeData
    return graphData

def define_edge(sourceNodeData,
                sourceNodeDirection,
                targetNodeData,
                targetNodeDirection,
                graphData):
    # this defines the edge leaving the source (may be in + or - direction) forward edge which we add to the list of edges for the source node
    sourceToTargetEdge = EdgeData(targetNodeData.get_nodeId(),  # the id of the target node
                                sourceNodeDirection,
                                targetNodeDirection,
                                1)
    revCompEdgeFromTargetToSource = EdgeData(targetNodeData.get_nodeId(), # the id of the target node
                                        -1 * sourceNodeDirection,
                                        -1 * targetNodeDirection,
                                        1)
    # this modifies the node data in the graphData list by index
    graphData = add_edge_to_graph(sourceNodeData,
                                sourceNodeDirection,
                                sourceToTargetEdge,
                                targetNodeData,
                                revCompEdgeFromTargetToSource,
                                graphData)
    return graphData

def build_graph(readData,
                kmer_size,
                refined_genes):
    """ this function returns a list of NodeData objects where each element represents a unique geneMer and its reverse complement """
    # initialise the graphData list so we can refer to nodes by index
    graphData = [None for i in range(50000)]
    # keep track of geneMers that we have seen before and ensure they are given consistent node IDs
    geneMerNameIdDict = {}
    geneMerIdNameDict = {}
    current_geneMerId = 0
    # for each read in the readData dictionary
    for read in readData:
        # get the id of the read
        readId = read.get_readId()
        # get the list of gene integers for each read
        annotatedGenes = read.get_geneIdList()
        # get the list of strands corresponding to each gene identifier
        annotatedGeneStrands = read.get_geneStrandList()
        # filter out reads that do not have at least k genes annotated
        if len(annotatedGenes) >= kmer_size:
            for i in range(len(annotatedGenes) - (kmer_size - 1)):
                # iterate through the list of genes by index and add a source node for all possible geneMers
                sourceNodeData, sourceNodeDirection, geneMerNameIdDict, geneMerIdNameDict, current_geneMerId = define_node(annotatedGenes[i: i + kmer_size],
                                                                                                                    annotatedGeneStrands[i: i + kmer_size],
                                                                                                                    geneMerNameIdDict,
                                                                                                                    geneMerIdNameDict,
                                                                                                                    current_geneMerId,
                                                                                                                    refined_genes,
                                                                                                                    readId)
                # we can only add a target node and therefore and edge up to the gene k genes away from the end
                if not i == (len(annotatedGenes) - (kmer_size)):
                    # define the data for the target node
                    targetNodeData, targetNodeDirection, geneMerNameIdDict, geneMerIdNameDict, current_geneMerId = define_node(annotatedGenes[i + 1: i + kmer_size + 1],
                                                                                                                        annotatedGeneStrands[i + 1: i + kmer_size + 1],
                                                                                                                        geneMerNameIdDict,
                                                                                                                        geneMerIdNameDict,
                                                                                                                        current_geneMerId,
                                                                                                                        refined_genes,
                                                                                                                        readId)
                    # add a forward edge from the source node to the target node and a backward edge from the target node to the source node
                    graphData = define_edge(sourceNodeData,
                                            sourceNodeDirection,
                                            targetNodeData,
                                            targetNodeDirection,
                                            graphData)
    return [node for node in graphData if node], geneMerIdNameDict