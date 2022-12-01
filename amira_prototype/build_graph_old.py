class GraphData:

    def __init__(self,
                nodeId,
                interesting_node,
                nodeName,
                nodeCoverage,
                nodeStrandList=[],
                edgeList=[],
                edgeIds = [],
                readIdSet=set()):
        self.nodeId = nodeId
        self.interesting_node = interesting_node
        self.nodeName = nodeName
        self.nodeCoverage = nodeCoverage
        self.nodeStrandList = nodeStrandList
        self.edgeList = edgeList
        self.edgeIds = edgeIds
        self.readIdSet = readIdSet

    def get_nodeId(self):
        return self.nodeId
    def get_interestingNode(self):
        return self.interesting_node
    def get_nodeName(self):
        return self.nodeName
    def get_rc_nodeName(self):
        nodeList = list(reversed(self.nodeName.split(",")))
        return ",".join(str(int(g)*-1) for g in nodeList)
    def get_nodeCoverage(self):
        return self.nodeCoverage
    def get_nodeStrandList(self):
        return self.nodeStrandList
    def get_edgeList(self):
        return self.edgeList
    def get_edgeIds(self):
        return self.edgeIds
    def get_readIdSet(self):
        return self.readIdSet

class EdgeData:

    def __init__(self,
                targetNodeId,
                targetNodeName,
                strand,
                edgeDirection,
                edgeWeight):
        self.targetNodeId = targetNodeId
        self.targetNodeName = targetNodeName
        self.strand = strand # whether the target is the complement or reverse complement
        self.edgeDirection = edgeDirection # direction of edge in the graph
        self.edgeWeight = edgeWeight

    def get_targetNodeId(self):
        return self.targetNodeId
    def get_targetNodeName(self):
        return self.targetNodeName
    def get_strand(self):
        return self.strand
    def get_edgeDirection(self):
        return self.edgeDirection
    def get_edgeWeight(self):
        return self.edgeWeight

def initialise_building_data_structures():
    return 0, {}, {}, []

def reverse_complement(gene):
    reverse_gene = reversed(gene)
    reverse_gene = [str(int(g)*-1) for g in reverse_gene]
    sorted_strands =  sorted([",".join(gene), ",".join(reverse_gene)])
    if not sorted_strands[0] == ",".join(gene):
        strand = -1
    else:
        strand = 1
    return sorted_strands[0], strand

def kmerise(gene_mer,
            gene_mer_strands,
            rc):
    if len(gene_mer) > 0:
        gene_mer = [str(gene_mer[g] * gene_mer_strands[g]) for g in range(len(gene_mer))]
        return ",".join(gene_mer), +1
    else:
        return None

def rc_geneMer(gene):
    reverse_gene = reversed(gene)
    reverse_gene = [str(int(g)*-1) for g in reverse_gene]
    return ",".join(reverse_gene)

def define_node(gene_mer,
                gene_mer_strands,
                geneMerId,
                geneMerNameIdDict,
                geneMerIdNameDict,
                readId,
                rc,
                refined_genes):
    if any(g in gene_mer or -1 * g in gene_mer for g in refined_genes):
        interesting = True
    else:
        interesting = False
    gene_mer_name, strand = kmerise(gene_mer,
                                    gene_mer_strands,
                                    rc)
    if rc:
        rc_genes = rc_geneMer(gene_mer_name.split(","))
        if not (gene_mer_name in geneMerNameIdDict or rc_genes in geneMerNameIdDict):
            this_geneMerId = geneMerId
            geneMerNameIdDict[gene_mer_name] = this_geneMerId
            geneMerIdNameDict[this_geneMerId] = gene_mer_name
            geneMerId += 1
            strand = 1
        else:
            if gene_mer_name in geneMerNameIdDict:
                this_geneMerId = geneMerNameIdDict[gene_mer_name]
                strand = 1
            else:
                this_geneMerId = geneMerNameIdDict[rc_genes]
                strand = -1
    else:
        if not gene_mer_name in geneMerNameIdDict:
            this_geneMerId = geneMerId
            geneMerNameIdDict[gene_mer_name] = this_geneMerId
            geneMerIdNameDict[this_geneMerId] = gene_mer_name
            geneMerId += 1
            strand = 1
        else:
            this_geneMerId = geneMerNameIdDict[gene_mer_name]
            strand = 1
    this_node = GraphData(this_geneMerId,
                        interesting,
                        gene_mer_name,
                        0,
                        [strand],
                        [],
                        [],
                        set([readId]))
    return this_node, geneMerId, geneMerNameIdDict, geneMerIdNameDict

def define_sourceNode(gene_mer,
                    gene_mer_strands,
                    geneMerId,
                    geneMerNameIdDict,
                    geneMerIdNameDict,
                    graphData,
                    readId,
                    rc,
                    refined_genes):
    sourceNodeData, geneMerId, geneMerNameIdDict, geneMerIdNameDict = define_node(gene_mer,
                                                                                gene_mer_strands,
                                                                                geneMerId,
                                                                                geneMerNameIdDict,
                                                                                geneMerIdNameDict,
                                                                                readId,
                                                                                rc,
                                                                                refined_genes)
    if len(graphData) == sourceNodeData.get_nodeId():
        sourceNodeData.nodeCoverage = 1
        graphData.append(sourceNodeData)
    else:
        graphData[sourceNodeData.get_nodeId()].nodeCoverage += 1
        graphData[sourceNodeData.get_nodeId()].nodeStrandList += sourceNodeData.get_nodeStrandList()
        graphData[sourceNodeData.get_nodeId()].readIdSet.update(sourceNodeData.get_readIdSet())
    return sourceNodeData, graphData, geneMerId, geneMerNameIdDict, geneMerIdNameDict

def define_targetNode(gene_mer,
                    gene_mer_strands,
                    geneMerId,
                    geneMerNameIdDict,
                    geneMerIdNameDict,
                    graphData,
                    readId,
                    rc,
                    refined_genes):
    targetNodeData, geneMerId, geneMerNameIdDict, geneMerIdNameDict = define_node(gene_mer,
                                                                                gene_mer_strands,
                                                                                geneMerId,
                                                                                geneMerNameIdDict,
                                                                                geneMerIdNameDict,
                                                                                readId,
                                                                                rc,
                                                                                refined_genes)
    if len(graphData) == targetNodeData.get_nodeId():
        targetNodeData.nodeCoverage = 0
        graphData.append(targetNodeData)
    else:
        graphData[targetNodeData.get_nodeId()].nodeStrandList += targetNodeData.get_nodeStrandList()
    return targetNodeData, graphData, geneMerId, geneMerNameIdDict, geneMerIdNameDict

def define_nodeEdges(sourceNodeData,
                    targetNodeData,
                    graphData):
    if targetNodeData.get_nodeId() > sourceNodeData.get_nodeId():
        edgeDirection = +1
    else:
        edgeDirection = -1
    if not targetNodeData.get_nodeStrandList()[0] == sourceNodeData.get_nodeStrandList()[0]:
        relative_node_strand = -1
    else:
        relative_node_strand = 1
    forward_edge = EdgeData(targetNodeData.get_nodeId(),
                        targetNodeData.get_nodeName(),
                        targetNodeData.get_nodeStrandList()[0],
                        edgeDirection * 1,
                        1)
    if not targetNodeData.get_nodeId() in graphData[sourceNodeData.get_nodeId()].get_edgeIds():
        graphData[sourceNodeData.get_nodeId()].edgeList.append(forward_edge)
        graphData[sourceNodeData.get_nodeId()].edgeIds.append(forward_edge.get_targetNodeId())
    else:
        graphData[sourceNodeData.get_nodeId()].edgeList[graphData[sourceNodeData.get_nodeId()].edgeIds.index(targetNodeData.get_nodeId())].edgeWeight += 1
    backward_edge = EdgeData(sourceNodeData.get_nodeId(),
                            sourceNodeData.get_nodeName(),
                            sourceNodeData.get_nodeStrandList()[0],
                            edgeDirection * -1,
                            1)
    if not sourceNodeData.get_nodeId() in graphData[targetNodeData.get_nodeId()].get_edgeIds():
        graphData[targetNodeData.get_nodeId()].edgeList.append(backward_edge)
        graphData[targetNodeData.get_nodeId()].edgeIds.append(backward_edge.get_targetNodeId())
    else:
        graphData[targetNodeData.get_nodeId()].edgeList[graphData[targetNodeData.get_nodeId()].edgeIds.index(sourceNodeData.get_nodeId())].edgeWeight += 1
    assert targetNodeData.get_nodeId() in graphData[sourceNodeData.get_nodeId()].edgeIds
    assert sourceNodeData.get_nodeId() in graphData[targetNodeData.get_nodeId()].edgeIds
    return graphData

def build_graph(readData,
                kmer_size,
                rc,
                refined_genes):
    geneMerId, geneMerNameIdDict, geneMerIdNameDict, graphData = initialise_building_data_structures()
    for read in readData:
        readId = read.get_readId()
        annotatedGenes = read.get_geneIdList()
        annotatedGeneStrands = read.get_geneStrandList()
        if len(annotatedGenes) >= kmer_size:
            for i in range(len(annotatedGenes) - (kmer_size - 1)):
                sourceNodeData, graphData, geneMerId, geneMerNameIdDict, geneMerIdNameDict = define_sourceNode(annotatedGenes[i: i + kmer_size],
                                                                                                            annotatedGeneStrands[i: i + kmer_size],
                                                                                                            geneMerId,
                                                                                                            geneMerNameIdDict,
                                                                                                            geneMerIdNameDict,
                                                                                                            graphData,
                                                                                                            readId,
                                                                                                            rc,
                                                                                                            refined_genes)
                if not i == (len(annotatedGenes) - (kmer_size)):
                    targetNodeData, graphData, geneMerId, geneMerNameIdDict, geneMerIdNameDict = define_targetNode(annotatedGenes[i + 1: i + kmer_size + 1],
                                                                                                                annotatedGeneStrands[i + 1: i + kmer_size + 1],
                                                                                                                geneMerId,
                                                                                                                geneMerNameIdDict,
                                                                                                                geneMerIdNameDict,
                                                                                                                graphData,
                                                                                                                readId,
                                                                                                                rc,
                                                                                                                refined_genes)
                    graphData = define_nodeEdges(sourceNodeData,
                                                targetNodeData,
                                                graphData)
    return graphData, geneMerIdNameDict
