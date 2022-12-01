class NodeData:
    def __init__(self,
                nodeId,
                nodeName,
                rc_nodeName,
                forwardEdgeList=[],
                reverseEdgeList=[],
                forwardEdgeIds = [],
                reverseEdgeIds = []):
        self.nodeId = nodeId # the integer ID of this node (also the index in the graphData list)
        self.nodeName = nodeName # the geneMer we have chosen as the canonical geneMer
        self.rc_nodeName = rc_nodeName
        self.forwardEdgeList = forwardEdgeList # edges in the forward direction from this node
        self.reverseEdgeList = reverseEdgeList # edges in the reverse direction from this node
        self.forwardEdgeIds = forwardEdgeIds # the list of target node IDs we see in the data to check if we have seen an edge before
        self.reverseEdgeIds = reverseEdgeIds # the list of target node IDs we see in the data to check if we have seen an edge before

    def get_nodeId(self):
        return self.nodeId
    def get_nodeName(self):
        return self.nodeName

    #Zam: the rc nodename is immutable, it never changes after construction. Could just store it?
    #also you have another function below (not a member of this class) which also does this rev comp?        
    def get_rc_nodeName(self):
        nodeList = list(reversed(self.nodeName.split(",")))
        return ",".join(str(int(g)*-1) for g in nodeList)


    def get_forwardEdgeList(self):
        return self.forwardEdgeList
    def get_reverseEdgeList(self):
        return self.reverseEdgeList
    def get_forwardEdgeIds(self):
        return self.forwardEdgeIds
    def get_reverseEdgeIds(self):
        return self.reverseEdgeIds


class EdgeData:

    def __init__(self,
                sourceDirection,
                targetNodeId,
                targetDirection,
                edgeWeight): 

        self.targetNodeId = targetNodeId # the target node ID (also the index in the graphData list)
        self.targetDirection = targetDirection # whether the target is the canonical kmer of the reverse complement of it
        self.sourceDirection = sourceDirection # direction that the edge leaves the source node
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
        self.edgeWeight = self.edgeWeight+1

#unchanged by zam
def rc_geneMer(geneMer):
    """ returns the reverse complement of a list of integers gene identifiers. """
    # get the reverse complement of the geneMer
    reverse_geneMer = reversed(geneMer)
    # the strand of the gene is indicated by whether the integer gene ID is positive or negative
    reverse_geneMer = [str(int(g)*-1) for g in reverse_geneMer]
    return reverse_geneMer



#Zam: what type of object is passed in? i've not modified this code.
def kmerise(geneMer):
 """ returns a ',' delimited string of gene integer kmers (each of which is guaranteed canonical), alternating with directions (+/-) """
     #Zam: i think it returns an alternating list of canonical geneMers and strands.
    geneMer = [str(geneMer[g]) for g in range(len(geneMer))]
    reverse_geneMer = rc_geneMer(geneMer)
    # this selects the canonical geneMer by choosing the smallest of the seen geneMer and its reverse complement
    canonical_geneMer = sorted([geneMer, reverse_geneMer])[0]
    # if the canonical geneMer is not the seen geneMer then we say the strand is -1
    if not canonical_geneMer == geneMer:
        strand = -1
    else:
        strand = +1
    return ",".join(canonical_geneMer), strand


#Zam: i'm confused by this. A node just has one geneMer, right? And this is basically a create-new-node-or-return-existing-one-if-it-already-exists function?
#     why are we passing in a list, and returning just one? Or are you passing in the whole list plus the index/id into that list?
def define_node(geneMerList,
                geneMerIds,
                current_geneMerId):
    """ returns a NodeData object for this node """
    # this decides what the canonical geneMer is
    canonical_geneMer, strand = kmerise(geneMerList) # strand here means whether we see the canonical geneMer or reverse complement of it in this case
    # this assigns node IDs to each canonical geneMer
    if canonical_geneMer in geneMerIds:
        # if we have seen the canonical geneMer before the nodeID is consistent
        nodeId = geneMerIds[canonical_geneMer]
    else:
        # if we have not seen the canonical geneMer before we get a new nodeID
        nodeId = current_geneMerId
        current_geneMerId += 1
    this_node = NodeData(nodeId,
                        canonical_geneMer,
                        [], # this is empty as it will be populated with edge information later
                        [], # this is empty as it will be populated with edge information later ZAM added this
                        [], # this is empty as it will be populated with edge information later Zam added this
                        []) # this is empty as it will be populated with edge information later
    return this_node, strand, geneMerIds, current_geneMerId



#Zam comment. I'd like to avoid double meanings of words. In your original version of this function,
#you have sourceForwardEdge and targetBackward edge, which makes sense, but we also have forward/reverse
#for orientation with respect to canonical, at each node. 

def add_edge_to_graph(sourceNodeData,
                    sourceDirection,
                    sourceToTargetEdge,
                    targetNodeData,
                    revCompEdgeFromTargetToSource,
                    graphData):

    """ this modifies the graphData list by index to add an edge between the source and target and then the target and the source """
    # add an edge to the target node to the list of edges for the source node
    sourceNodeId = sourceNodeData.get_nodeId()
    # we initialised the graphData list earlier so here we are changing the entry from None to the node information if this is the first time we see the geneMer

    #Zam added
    if sourceNodeId > 50000
        #error out, you've exceeded the pre-initialised length of the graphData structure--> or do something better so this not needed

    if graphData[sourceNodeId] == None:
        graphData[sourceNodeId] = sourceNodeData

    targetNodeId = targetNodeData.get_nodeId()
    if graphData[targetNodeId] == None:
        graphData[targetNodeId] = targetNodeData

    if (sourceDirection==1)
        if not targetNodeId in graphData[sourceNodeId].get_forwardEdgeIds():
            graphData[sourceNodeId].forwardEdgeList.append(sourceToTargetEdge)
            graphData[sourceNodeId].forwardEdgeIds.append(targetNodeId)
        else 
            #get sourceToTargetEdge from graphData[sourceNodeId].forwardEdgeList and use increment_weight() on it

    else #we are leaving the source node in the reverse direction
       if not targetNodeId in graphData[sourceNodeId].get_reverseEdgeIds():
            graphData[sourceNodeId].reverseEdgeList.append(sourceToTargetEdge)
            graphData[sourceNodeId].reverseEdgeIds.append(targetNodeId)           
        else 
            #get sourceToTargetEdge from graphData[sourceNodeId].reverseEdgeList and use increment_weight() on it
    return graphData
    ##Zam - given we pass graphData in, why do we need to return it?



#Zam - comment; i've systematically removed "strand" and  put direction in instead. The reason is
# strand is a bit confusing as you are led to think about strand of the whole molecule/sequence.
#whereas this is just canonical-direction-or-not. Not sure if it really helps, but that was my reasoning

##i've also got rid of all the relative strand stuff which is just confusing, and i think wrong.
## all you need to know is, when you arrive at the target, do you arrive in the forward or reverse direction
## of that node (with respect to its canonical kmer)
def define_edge(sourceNodeData,
                sourceNodeDirection,
                targetNodeData,
                targetNodeDirection,
                graphData):


    # this defines the edge leaving the source (may be in + or - direction) forward edge which we add to the list of edges for the source node


    sourceToTargetEdge = EdgeData(sourceNodeDirection,
                                 targetNodeData.get_nodeId(), # the id of the target node
                                targetNodeDirection)
    revCompEdgeFromTargetToSource = EdgeData(-targetNodeDirection, #Zam: what type is this direction. can i just minus it?
                                targetNodeData.get_nodeId(), # the id of the source node
                                sourceNodeData.get_nodeId(),
                                -sourceNodeDirection) #Zam: what type is this direction. can i just minus it?


def add_edge_to_graph(sourceNodeData,
                    sourceDirection,
                    sourceToTargetEdge,
                    targetNodeData,
                    revCompEdgeFromTargetToSource,
                    graphData):


    # this modifies the node data in the graphData list by index
    graphData = add_edge_to_graph(sourceNodeData,
                                sourceToTargetEdge,
                                targetNodeData,
                                revCompEdgeFromTargetToSource,
                                graphData)
    return graphData

def build_graph(readData,
                kmer_size):
    """ this function returns a list of NodeData objects where each element represents a unique geneMer and its reverse complement """
    # initialise the graphData list so we can refer to nodes by index
    graphData = [None for i in range(50000)]
    # keep track of geneMers that we have seen before and ensure they are given consistent node IDs
    geneMerIds = {}
    current_geneMerId = 0
    # for each read in the readData dictionary
    for read in readData:
        # get the list of gene integers for each read
        annotatedGenes = read.get_geneIdList()
        # filter out reads that do not have at least k+1 genes annotated
        if len(annotatedGenes) > kmer_size:  # Zam - changed to k+1 here

            for i in range(len(annotatedGenes) - (kmer_size - 1)):
                # iterate through the list of genes by index and add a source node for all possible geneMers
                sourceNodeData, sourceNodeDirection, geneMerIds, current_geneMerId = define_node(annotatedGenes[i: i + kmer_size],
                                                                                            geneMerIds,
                                                                                            current_geneMerId)
                # define the data for the target node
                targetNodeData, targetNodeDirection, geneMerIds, current_geneMerId = define_node(annotatedGenes[i + 1: i + kmer_size + 1],
                                                                                                geneMerIds,
                                                                                                current_geneMerId)
                    # add an edge from the source node to the target node and a reverse complement edge from the target node to the source node
                    graphData = define_edge(sourceNodeData,
                                            sourceNodeDirection,
                                            targetNodeData,
                                            targetNodeDirection,
                                            graphData)
    return [node for node in graphData if node]