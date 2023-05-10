import gzip
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import os
import subprocess
import sys
from tqdm import tqdm

from construct_graph import GeneMerGraph
from construct_node import Node
from construct_gene_mer import define_rc_geneMer
from construct_gene import convert_int_strand_to_string

class Unitigs:
    def __init__(self,
                graph: GeneMerGraph,
                listOfGenes: list):
            self._graph = graph
            self._listOfGenes = listOfGenes
    def get_graph(self):
        """ returns the geneMerGraph """
        return self._graph
    def get_selected_genes(self):
        """ returns the list of selected genes """
        return self._listOfGenes
    def get_nodes_of_interest(self,
                            geneOfInterest):
        """ extracts the graph nodes containing the genes of interest and returns them as a list """
        return self.get_graph().get_nodes_containing(geneOfInterest)
    # def get_forward_node_from_node(self,
    #                             sourceNode: Node,
    #                             allAMRHashes) -> list:
    #     """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
    #     # get the list of forward edge hashes for this node
    #     nodeForwardEdges = sourceNode.get_forward_edge_hashes()
    #     if len(nodeForwardEdges) > 0:
    #         targetNodes = []
    #         targetNodeDirections = []
    #         # iterate through the edge hashes
    #         for edge_hash in nodeForwardEdges:
    #             # get the edge object corresponding to this edge hash
    #             edge = self.get_graph().get_edges()[edge_hash]
    #             # get the target node for this edge
    #             targetNode = edge.get_targetNode()
    #             # check if the target node contains an AMR gene
    #             if targetNode.__hash__() in allAMRHashes:
    #                 targetNodes.append(edge.get_targetNode())
    #                 # get the direction we are going into the target node
    #                 targetNodeDirections.append(edge.get_targetNodeDirection())
    #         return targetNodes, targetNodeDirections
    #     else:
    #         return None, None
    # def get_backward_node_from_node(self,
    #                             sourceNode: Node,
    #                             allAMRHashes) -> list:
    #     """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
    #     # get the list of forward edge hashes for this node
    #     nodeBackwardEdges = sourceNode.get_backward_edge_hashes()
    #     if len(nodeBackwardEdges) > 0:
    #         targetNodes = []
    #         targetNodeDirections = []
    #         # iterate through the edge hashes
    #         for edge_hash in nodeBackwardEdges:
    #             # get the edge object corresponding to this edge hash
    #             edge = self.get_graph().get_edges()[edge_hash]
    #             # get the target node for this edge
    #             targetNode = edge.get_targetNode()
    #             # check if the target node contains an AMR gene
    #             if targetNode.__hash__() in allAMRHashes:
    #                 targetNodes.append(edge.get_targetNode())
    #                 # get the direction we are going into the target node
    #                 targetNodeDirections.append(edge.get_targetNodeDirection())
    #         return targetNodes, targetNodeDirections
    #     else:
    #         return None, None
    # def get_node_geneMer(self,
    #                     node,
    #                     nodeDirection):
    #     if nodeDirection == 1:
    #         return node.get_canonical_geneMer()
    #     else:
    #         return node.get_geneMer().get_rc_geneMer()
    def get_all_nodes_containing_AMR_genes(self):
        """ return a dictionary of nodes containing AMR genes """
        AMRNodes = {}
        # iterate through the list of specified genes
        for geneOfInterest in self.get_selected_genes():
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_of_interest(geneOfInterest)
            # add nodes of interest to the AMR node set
            for n in nodesOfInterest:
                AMRNodes[n.__hash__()] = n
        return AMRNodes
    def get_AMR_anchors(self,
                        AMRNodes):
        nodeAnchors = set()
        nodeJunctions = set()
        # get nodes that are anchors for traversing in the forward direction of the graph
        for nodeHash in AMRNodes:
            # get the forward edges for this node
            forward_edges = [self.get_graph().get_edge_by_hash(h) for h in AMRNodes[nodeHash].get_forward_edge_hashes()]
            # get the backward edges for this node
            backward_edges = [self.get_graph().get_edge_by_hash(h) for h in AMRNodes[nodeHash].get_backward_edge_hashes()]
            # get the forward edges that contain an AMR gene
            forwardAMRNodes = [e.get_targetNode() for e in forward_edges] #if e.get_targetNode().__hash__() in AMRNodes]
            # get the backward edges that contain an AMR gene
            backwardAMRNodes = [e.get_targetNode() for e in backward_edges] #if e.get_targetNode().__hash__() in AMRNodes]
            # if the number of backward neighbors is not 1 then this is a forward anchor
            if len(forward_edges) + len(backward_edges) > 2:
                nodeJunctions.add(nodeHash)
            else:
                if len(backwardAMRNodes) == 0 or len(forwardAMRNodes) == 0:
                    nodeAnchors.add(nodeHash)
        return nodeAnchors, nodeJunctions
    def write_clusters(self,
                    output_dir,
                    clusters):
        fastqContent = parse_fastq("/home/daniel/Documents/GitHub/amira_data/E.coli/pandora.paper.data/data/samples/Escherichia_coli_MINF_9A/Escherichia_coli_MINF_9A.nanopore.fastq.gz")
        pathId = 1
        for nodeHash in tqdm(clusters):
            if len(clusters[nodeHash]) > 1:
                # make the output directory
                if not os.path.exists(os.path.join(output_dir, str(pathId))):
                    os.mkdir(os.path.join(output_dir, str(pathId)))
                # subset the reads for this unitig
                unitigReads = list(clusters[nodeHash].keys())
                subsettedReadData = {}
                for r in unitigReads:
                    subsettedReadData[r] =  fastqContent[r]
                # write the per unitig fastq data
                readFileName = os.path.join(output_dir, str(pathId), str(pathId) + ".fastq.gz")
                write_fastq(readFileName,
                            subsettedReadData)
                # increment the path ID
                pathId += 1
    def make_histogram(self,
                    outputDir,
                    numbersToPlot,
                    y_label,
                    x_label,
                    histogramFile):
        import numpy as np
        # make the output dir if it doesn't exist
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)
        # plot a bar chart of the number of genes per read
        q25, q75 = np.percentile(numbersToPlot, [25, 75])
        bin_width = 2 * (q75 - q25) * len(numbersToPlot) ** (-1/3)
        numBins = round((max(numbersToPlot) - min(numbersToPlot)) / bin_width)
        fig, ax = plt.subplots()
        ax.hist(numbersToPlot, bins=numBins, density = False)
        plt.margins(x=0)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        fig.savefig(os.path.join(outputDir, histogramFile))
        plt.close(fig)
    def get_informative_reads(self,
                            AMRNodes,
                            allReadNodes,
                            graph,
                            allAnchors):
        informativeReads = {}
        anchorCounts = []
        seenReads = set()
        for nodeHash in tqdm(AMRNodes):
            node = graph.get_node_by_hash(nodeHash)
            # get the reads for this node
            nodeReads = [r for r in node.get_reads()]
            # get the nodes covered by each read
            nodesOnReads = [allReadNodes[r] for r in nodeReads]
            # get the anchors and junctions on each read
            anchorsOnReads = [[n for n in r if n in allAnchors] for r in nodesOnReads]
            # get the read lengths for this node
            readLengths = [len(r) for r in anchorsOnReads]
            for r in range(len(readLengths)):
                if not nodeReads[r] in seenReads:
                    anchorCounts.append(readLengths[r])
                    seenReads.add(nodeReads[r])
            # the longest reads for this node will be become the representatives
            for r in range(len(readLengths)):
                if readLengths[r] > 1 or (readLengths[r] == 1 and len(nodesOnReads[r]) == 1):
                    startMask = nodesOnReads[r].index(anchorsOnReads[r][0])
                    endMask = len(nodesOnReads[r]) - 1 - nodesOnReads[r][::-1].index(anchorsOnReads[r][-1])
                    informativeReads[nodeReads[r]] = nodesOnReads[r][startMask: endMask + 1]
        return informativeReads
    def assign_to_regions(self,
                        informativeReads,
                        allAnchors):
        regions = set()
        for readId in tqdm(informativeReads):
            this_cluster = set([readId])
            allRegionNodeAnchors = set([n for n in informativeReads[readId] if n in allAnchors])
            for otherReadId in informativeReads:
                if not readId == otherReadId:
                    otherNodeAnchors = [n for n in informativeReads[otherReadId] if n in allAnchors]
                    if any(n in allRegionNodeAnchors for n in otherNodeAnchors):
                        for nodeHash in otherNodeAnchors:
                            allRegionNodeAnchors.add(nodeHash)
                        this_cluster.add(otherReadId)
            this_cluster = tuple(sorted(list(this_cluster)))
            regions.add(this_cluster)
        return list(regions)
    def contains_sublist(self,
                        allNodes,
                        subNodes):
        n = len(subNodes)
        if any((subNodes == allNodes[i:i+n]) for i in range(len(allNodes)-n+1)):
            return True
        else:
            return False
    def cluster_region_reads(self,
                            regions,
                            informativeReads,
                            allReadNodes):
        clusters = {}
        seenReads = set()
        for regionReads in tqdm(regions):
            # convert the read tuple to a list
            regionReadList = list(regionReads)
            # get the actual nodes for all the reads in this read list
            #regionReadNodes = [informativeReads[r] if r in informativeReads else allReadNodes[r] for r in regionReadList]
            regionReadNodes = [allReadNodes[r] for r in regionReadList]
            for r1 in range(len(regionReadNodes)):
                readId = regionReads[r1]
                if not readId in clusters:
                    clusters[readId] = set([readId])
                for r2 in range(len(regionReadNodes)):
                    otherReadId = regionReads[r2]
                    if not readId == otherReadId:
                        if self.contains_sublist(regionReadNodes[r1], regionReadNodes[r2]):
                            clusters[readId].add(otherReadId)#[otherReadId] = regionReadNodes[r2]
                            seenReads.add(otherReadId)
                        else:
                            if self.contains_sublist(regionReadNodes[r1], list(reversed(regionReadNodes[r2]))):
                                clusters[readId].add(otherReadId)#[otherReadId] = list(reversed(regionReadNodes[r2]))
                                seenReads.add(otherReadId)
        # for r1 in clusters:
        #     for r2 in clusters[r1]:
        #         #if r2 in informativeReads:
        #           #  pathNodes = str(r2) + ":\t" + ", ".join([str(n) for n in informativeReads[r2]])
        #         #else:
        #         pathNodes = str(r2) + ":\t" + ", ".join([str(n) for n in allReadNodes[r2]])
        #         with open("pathClusters_single_included.txt", "a+") as o:
        #             o.write(pathNodes + "\n")
        #     with open("pathClusters_single_included.txt", "a+") as o:
        #             o.write("\n")
        for r in seenReads:
            if r in clusters:
                del clusters[r]
        # for r1 in clusters:
        #     for r2 in clusters[r1]:
        #         #if r2 in informativeReads:
        #             #pathNodes = str(r2) + ":\t" + ", ".join([str(n) for n in informativeReads[r2]])
        #         #else:
        #         pathNodes = str(r2) + ":\t" + ", ".join([str(n) for n in allReadNodes[r2]])
        #         with open("pathClusters_trimmed_single_included.txt", "a+") as o:
        #             o.write(pathNodes + "\n")
        #     with open("pathClusters_trimmed_single_included.txt", "a+") as o:
        #             o.write("\n")
        return clusters
    def supplement_regions(self,
                        regions,
                        informativeReads,
                        allReadNodes):
        sys.stderr.write("\nSupplementing regions\n")
        supplementedReads = set()
        for readId in tqdm(allReadNodes):
            if not readId in informativeReads:
                readNodes = set(allReadNodes[readId])
                if len(readNodes) > 0:
                    for reg in range(len(regions)):
                        regionReads = set(list(regions[reg]))
                        regionNodes = [informativeReads[r] for r in regionReads if not r in supplementedReads]
                        for regionReadNodes in regionNodes:
                            if any(n in readNodes for n in regionReadNodes):
                                regionReads.add(readId)
                                supplementedReads.add(readId)
                                continue
                        regions[reg] = tuple(list(regionReads))
        return regions
    def recluster_regions(self,
                        regions):
        newRegions = {}
        seenRegions = set()
        for regionReads in tqdm(regions):
            newRegions[regionReads] = list(regionReads)
            regionReadSet = set(regionReads)
            for otherRegionReads in regions:
                if not regionReads == otherRegionReads:
                    if any(r in regionReadSet for r in list(otherRegionReads)):
                        newRegions[regionReads] += list(otherRegionReads)
                        seenRegions.add(otherRegionReads)
                        continue
        for r in seenRegions:
            del newRegions[r]
        regions = set()
        for c in newRegions:
            regions.add(tuple(set(sorted(list(newRegions[c])))))
        return list(regions)
    def old_get_unitigs_of_interest(self):
        # isolate nodes containing AMR genes
        AMRNodes = self.get_all_nodes_containing_AMR_genes()
        # get anchors and junctions
        nodeAnchors, nodeJunctions = self.get_AMR_anchors(AMRNodes)
        allAnchors = nodeAnchors.union(nodeJunctions)
        # get the graph
        graph = self.get_graph()
        # get the reads as nodes
        allReadNodes = graph.get_readNodes()
        # get the longest read for each AMR node
        informativeReads = self.get_informative_reads(AMRNodes,
                                                    allReadNodes,
                                                    graph,
                                                    allAnchors)
        # cluster nodes into regions of the graph
        regions = self.assign_to_regions(informativeReads,
                                        allAnchors)
        regions = self.recluster_regions(regions)
        # for c in regions:
        #     clusterList = list(c)
        #     clusterNodes = [str(r) + ":\t" + ", ".join([str(h) for h in informativeReads[r]]) for r in clusterList]
        #     with open("readClusters_informative.txt", "a+") as o:
        #         o.write("\n".join(clusterNodes) + "\n\n")
        # supplement regions with non informative reads
        regions = self.supplement_regions(regions,
                                        informativeReads,
                                        allReadNodes)
        # for c in regions:
        #     clusterList = list(c)
        #     clusterNodes = []
        #     for r in clusterList:
        #         if r in informativeReads:
        #             nodePath = str(r) + ":\t" + ", ".join([str(h) for h in informativeReads[r]])
        #         else:
        #             nodePath = str(r) + ":\t" + ", ".join([str(h) for h in allReadNodes[r]])
        #         clusterNodes.append(nodePath)
        #     with open("readClusters_single_included.txt", "a+") as o:
        #         o.write("\n".join(clusterNodes) + "\n\n")
        # cluster the reads in each region based on the path they follow
        clusteredRegions = self.cluster_region_reads(regions,
                                                    informativeReads,
                                                    allReadNodes)
        nodes = {}
        for r in clusteredRegions:
            nodes[r] = [0]
        return nodes, clusteredRegions
        end
    def get_unitigs_of_interest(self):
        # isolate nodes containing AMR genes
        AMRNodes = self.get_all_nodes_containing_AMR_genes()
        # get anchors and junctions
        nodeAnchors, nodeJunctions = self.get_AMR_anchors(AMRNodes)
        # get the graph
        graph = self.get_graph()
        # keep track of what clusters reads are assigned to
        clusterReads = {}
        clusterReferences = {}
        clusterId = 1
        # parse the fastq
        fastqContent = parse_fastq("/home/daniel/Documents/GitHub/amira_data/E.coli/pandora.paper.data/data/samples/063_STEC/063_STEC.nanopore.fastq.gz")
        # iterate through nodes containing AMR genes
        for nodeHash in tqdm(AMRNodes):
            # get the reads for this node
            node = graph.get_node_by_hash(nodeHash)
            nodeReads = [r for r in node.get_reads()]
            this_cluster = clusterId
            clusterReferences[this_cluster] = set()
            clusterReads[this_cluster] = set()
            for readId in nodeReads:
                # get the nodes that are anchors for this readID
                readAnchors = [n.__hash__() for n in graph.get_nodes_containing_read(readId) if n.__hash__() in nodeAnchors or n.__hash__() in nodeJunctions]
                if len(readAnchors) > 2:
                    clusterReferences[this_cluster].add(readId)
                else:
                    clusterReads[this_cluster].add(readId)
            clusterId += 1
        empty = 0
        not_empty = 0
        for cluster in tqdm(clusterReferences):
            if not (len(clusterReferences[cluster]) == 0 and len(clusterReads[cluster]) == 0):
                subprocess.run("mkdir -p " + os.path.join("/home/daniel/Documents/GitHub/amira_prototype/amira.output/063_STEC.nanopore/read_clusters", str(cluster)),
                            shell=True,
                            check=True)
                not_empty += 1
                referenceSequences = []
                for readId in clusterReferences[cluster]:
                    referenceSequences.append(">" + readId + "\n" + fastqContent[readId]["sequence"])
                with open(os.path.join("/home/daniel/Documents/GitHub/amira_prototype/amira.output/063_STEC.nanopore/read_clusters", str(cluster), "reference_reads.fasta"), "w") as o:
                    o.write("\n".join(referenceSequences))
                querySequences = {}
                for readId in clusterReads[cluster]:
                    querySequences[readId] = fastqContent[readId]
                write_fastq(os.path.join("/home/daniel/Documents/GitHub/amira_prototype/amira.output/063_STEC.nanopore/read_clusters", str(cluster), "query_reads.fastq.gz"),
                            querySequences)
                subprocess.run("minimap2 -a --MD -t 8 -x asm20 " + os.path.join("/home/daniel/Documents/GitHub/amira_prototype/amira.output/063_STEC.nanopore/read_clusters", str(cluster), "reference_reads.fasta") + " " + os.path.join("/home/daniel/Documents/GitHub/amira_prototype/amira.output/063_STEC.nanopore/read_clusters", str(cluster), "query_reads.fastq.gz") + " > " + os.path.join("/home/daniel/Documents/GitHub/amira_prototype/amira.output/063_STEC.nanopore/read_clusters", str(cluster), "mapped_to_references.sam"),
                            shell=True,
                            check=True)
            else:
                empty += 1
        print(not_empty, empty)
        end

    # def get_unitigs_of_interest(self):
    #     # isolate nodes containing AMR genes
    #     AMRNodes = self.get_all_nodes_containing_AMR_genes()
    #     # get anchors and junctions
    #     nodeAnchors, nodeJunctions = self.get_AMR_anchors(AMRNodes)
    #     # count how many red nodes have another red node next to them
    #     graph = self.get_graph()
    #     sharedReads = {}
    #     for nodeHash in tqdm(nodeJunctions):
    #         node = graph.get_node_by_hash(nodeHash)
    #         fw_neighbors = [graph.get_edge_by_hash(e).get_targetNode() for e in node.get_forward_edge_hashes()]
    #         bw_neighbors = [graph.get_edge_by_hash(e).get_targetNode() for e in node.get_backward_edge_hashes()]
    #         for fw_node in fw_neighbors:
    #             fw_reads = set([r for r in fw_node.get_reads()])
    #             for bw_node in bw_neighbors:
    #                 bw_reads = set([r for r in bw_node.get_reads()])
    #                 readIntersection = fw_reads.intersection(bw_reads)
    #                 if len(readIntersection) > 4:
    #                     sharedReads[tuple([fw_node.__hash__(), bw_node.__hash__()])] = readIntersection
    #     # # go through the shared read sets and see if any share a read
    #     # newPaths = {}
    #     # for p1 in tqdm(sharedReads):
    #     #     intersects = False
    #     #     for p2 in sharedReads:
    #     #         if not p1 == p2:
    #     #             readIntersection = sharedReads[p1].intersection(sharedReads[p2])
    #     #             if len(readIntersection) > 0:
    #     #                 pathCombined = list(p1) + list(p2)
    #     #                 newPaths[tuple(sorted(pathCombined))] = readIntersection
    #     #                 intersects= True
    #     #     if not intersects:
    #     #         newPaths[p1] = sharedReads[p1]
    #     # import statistics
    #     # print(len(newPaths), statistics.mean(list(newPaths.values())), statistics.mode(list(newPaths.values())))
    #     readClusterId = {}
    #     clusterId = 0
    #     readsAsClusters = {}
    #     for p1 in tqdm(sharedReads):
    #         readSet = sharedReads[p1]
    #         if any(r in readClusterId for r in readSet):
    #             for read in readSet:
    #                 if read in readClusterId:
    #                     this_id = readClusterId[read]
    #                     continue
    #         else:
    #             this_id = clusterId
    #             clusterId += 1
    #         for read in readSet:
    #             readClusterId[read] = this_id
    #     for read in readClusterId:
    #         if not readClusterId[read] in readsAsClusters:
    #             readsAsClusters[readClusterId[read]] = set()
    #         readsAsClusters[readClusterId[read]].add(read)
    #     import statistics
    #     print(len(readsAsClusters), statistics.mean([len(l) for l in list(readsAsClusters.values())]), statistics.mode([len(l) for l in list(readsAsClusters.values())]))
    #     fastqContent = parse_fastq("/home/daniel/Documents/GitHub/amira_data/E.coli/pandora.paper.data/data/samples/063_STEC/063_STEC.nanopore.fastq.gz")
    #     for cluster in readsAsClusters:
    #         subsettedReadData = {}
    #         for r in readsAsClusters[cluster]:
    #             subsettedReadData[r] =  fastqContent[r]
    #         # write the per unitig fastq data
    #         if not os.path.exists(os.path.join("new_clustering_amira_output/unitigs", str(cluster))):
    #             os.mkdir(os.path.join("new_clustering_amira_output/unitigs", str(cluster)))
    #         readFileName = os.path.join("new_clustering_amira_output/unitigs", str(cluster), str(cluster) + ".fastq.gz")
    #         write_fastq(readFileName,
    #                     subsettedReadData)
        end


class UnitigTools:
    def __init__(self,
                graph,
                genesOfInterest):
        self._unitigs = Unitigs(graph,
                            genesOfInterest)
        self._unitigsOfInterest, self._unitigReadsOfInterest = self.get_unitigs().get_unitigs_of_interest()
    def get_unitigs(self) -> Unitigs:
        """ returns the Unitigs object containing the specified AMR genes """
        return self._unitigs
    def get_unitigsOfInterest(self) -> dict:
        """ returns all unitigs containing the specified AMR genes """
        return self._unitigsOfInterest
    def get_unitigReadsOfInterest(self) -> dict:
        """ returns all unitigs containing the specified AMR genes """
        return self._unitigReadsOfInterest
    def get_uniqueReads(self) -> dict:
        """ returns a dictionary of reads that are unique to a single path """
        return self._uniqueReads
    def convert_paths_to_genes(self,
                            pathGeneMers):
        geneMerGenes = [convert_int_strand_to_string(gene.get_strand()) + gene.get_name() for gene in pathGeneMers[0]]
        for geneMer in pathGeneMers[1:]:
            gene_strand = convert_int_strand_to_string(geneMer[-1].get_strand())
            gene_name = geneMer[-1].get_name()
            geneMerGenes.append(gene_strand + gene_name)
        return geneMerGenes
    def separate_reads(self,
                output_dir,
                readFile):
        fastqContent = parse_fastq(readFile)
        readFiles = []
        # get the untig genes and reads
        unitigsOfInterest = self.get_unitigsOfInterest()
        unitigsReadsOfInterest = self.get_unitigReadsOfInterest()
        # store the unitig hash to integer mappings
        unitig_mapping = {}
        # make the output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        # assign integers to paths for readability
        pathId = 1
        # iterate through the unitigs
        for path in tqdm(unitigsOfInterest):
            if not len(unitigsOfInterest[path]) == 0:
                # make the output directory
                if not os.path.exists(os.path.join(output_dir, str(pathId))):
                    os.mkdir(os.path.join(output_dir, str(pathId)))
                # get a readable list of genes in this path
                #readableGenes = self.convert_paths_to_genes(unitigsOfInterest[path])
                # write out the list of genes
                #with open(os.path.join(output_dir, str(pathId), "annotated_genes.txt"), "w") as outGenes:
                #    outGenes.write("\n".join(readableGenes))
                # subset the reads for this unitig
                unitigReads = unitigsReadsOfInterest[path]
                subsettedReadData = {}
                for r in unitigReads:
                    subsettedReadData[r] =  fastqContent[r]
                # write the per unitig fastq data
                readFileName = os.path.join(output_dir, str(pathId), str(pathId) + ".fastq.gz")
                write_fastq(readFileName,
                            subsettedReadData)
                readFiles.append(readFileName)
                # update the unitig IDs
                unitig_mapping[path] = pathId
                # increment the path ID
                pathId += 1
        return readFiles, unitig_mapping
    def multithread_flye(self,
                        readFiles,
                        flye_path,
                        threads):

        def run_flye(inputFastq,
                flye_path,
                flye_threads):
            flye_command = " ".join([flye_path,
                            "--nano-raw",
                            inputFastq,
                            "-t",
                            str(flye_threads),
                            "--out-dir",
                            os.path.join(os.path.dirname(inputFastq), "flye_output")])
            try:
                subprocess.run(flye_command, shell=True, check=True)
            except:
                pass

        job_list = [
            readFiles[i:i + threads] for i in range(0, len(readFiles), threads)
        ]
        for subset in tqdm(job_list):
            Parallel(n_jobs=threads)(delayed(run_flye)(r,
                                                flye_path,
                                                1) for r in subset)
    def multithread_raven(self,
                        readFiles,
                        raven_path,
                        threads):

        def run_raven(inputFastq,
                    raven_path,
                    raven_threads):
            outputConsensus = os.path.join(os.path.dirname(inputFastq), "raven_assembly.fa")
            raven_command = " ".join([raven_path,
                                    "-t",
                                    str(raven_threads),
                                    inputFastq,
                                    ">",
                                    outputConsensus])
            try:
                subprocess.run(raven_command, shell=True, check=True)
            except:
                pass

        job_list = [
            readFiles[i:i + threads] for i in range(0, len(readFiles), threads)
        ]
        for subset in tqdm(job_list):
            Parallel(n_jobs=threads)(delayed(run_raven)(r,
                                                    raven_path,
                                                    1) for r in subset)
    def polish_pandora_consensus(self,
                            readFiles,
                            racon_path,
                            consensusFastq,
                            genesOfInterest,
                            threads):

        def run_racon(file,
                    genesOfInterest,
                    fastqContent):
            # get the AMR genes in this unitig
            with open(os.path.join(os.path.dirname(file), "annotated_genes.txt"), "r") as inGenes:
                unitigGenes = [g[1:] for g in inGenes.read().split("\n") if g[1:] in genesOfInterest]
            # define the output file
            outputDir = os.path.join(os.path.dirname(file), "pandora.polished.consensus")
            if not os.path.exists(outputDir):
                os.mkdir(outputDir)
            # make a temp file for each AMR gene consensus sequence
            geneCounts = {}
            for gene in unitigGenes:
                if unitigGenes.count(gene) > 1:
                    if not gene in geneCounts:
                        geneCounts[gene] = 0
                    else:
                        geneCounts[gene] += 1
                    gene = gene + "_" + str(geneCounts[gene])
                # make the gene output dir
                if not os.path.exists(os.path.join(outputDir, gene)):
                    os.mkdir(os.path.join(outputDir, gene))
                # write out the pandora consensus
                try:
                    with open(os.path.join(outputDir, gene, "01.pandora.consensus.fasta"), "w") as outFasta:
                        outFasta.write(">" + gene + "\n" + "N"*500 + fastqContent[gene]["sequence"] + "N"*500 + "\n")
                except KeyError:
                    continue
                # map the reads to the consensus file
                map_command = "minimap2 -a --MD -t 1 -x asm20 "
                map_command += os.path.join(outputDir, gene, "01.pandora.consensus.fasta") + " " + file
                map_command += " > " + os.path.join(outputDir, gene, "02.read.mapped.sam")
                subprocess.run(map_command, shell=True, check=True)
                # polish the pandora consensus
                racon_command = racon_path + " -t 1 -w " + str(len(fastqContent[gene]["sequence"])) + " "
                racon_command += file + " " + os.path.join(outputDir, gene, "02.read.mapped.sam") + " "
                racon_command += os.path.join(outputDir, gene, "01.pandora.consensus.fasta") + " "
                racon_command += "> " + os.path.join(outputDir, gene, "03.polished.consensus.fasta")
                try:
                    subprocess.run(racon_command, shell=True, check=True)
                    # trim the buffer
                    with open(os.path.join(outputDir, gene, "03.polished.consensus.fasta"), "r") as i:
                        racon_seq = i.read()
                    header = racon_seq.split("\n")[0]
                    sequence = "\n".join(racon_seq.split("\n")[1:]).replace("N", "")
                    with open(os.path.join(outputDir, gene, "04.polished.consensus.trimmed.fasta"), "w") as outTrimmed:
                        outTrimmed.write(">" + header + "\n" + sequence + "\n")
                except:
                    pass

        # load the pandora consensus fastq
        fastqContent = parse_fastq(consensusFastq)
        # iterate through the reads for each unitig
        job_list = [
            readFiles[i:i + threads] for i in range(0, len(readFiles), threads)
        ]
        for subset in tqdm(job_list):
            Parallel(n_jobs=threads)(delayed(run_racon)(r,
                                                    genesOfInterest,
                                                    fastqContent) for r in subset)
    def initialise_plots(self,
                        unitigCount):
        #plt.rc('font', size=2)
        fig, ax = plt.subplots()
        # set the x-spine
        ax.spines['left'].set_position('zero')
        # turn off the right spine/ticks
        ax.spines['right'].set_color('none')
        # set the y-spine
        ax.spines['bottom'].set_position('zero')
        # turn off the top spine/ticks
        ax.spines['top'].set_color('none')
        ax.set_ylim([0, unitigCount + 1])
        # set the acis labels
        ax.set_ylabel("Unitig ID")
        ax.set_xlabel("Unitig length (bp)")
        return fig, ax
    def get_gene_colour(self,
                        currentGene,
                        allAMRGenes):
        if currentGene in allAMRGenes:
            color = "red"
        else:
            color = "lightgrey"
        return color
    def get_gene_coordinates(self,
                            geneStrand,
                            previousCoordinates,
                            geneLength,
                            unitigCount):
        if geneStrand == "+":
            x = sum(previousCoordinates)
            y = unitigCount
            dx = geneLength
            dy = 0
        else:
            x = sum(previousCoordinates) + geneLength
            y = unitigCount
            dx = geneLength * -1
            dy = 0
        return x, y, dx, dy
    def get_minimum_gene_length(self,
                                geneLengths,
                                listOfGenes):
        unitiglengths = [max(geneLengths[g[1:]]) for g in listOfGenes]
        return unitiglengths
    def visualise_unitigs(self,
                        geneLengths: dict,
                        unitig_mapping: dict,
                        output_dir: str):
        """ generate a figure to visualise the genes, order of genes, direction and lengths of genes on unitigs containing AMR genes """
        allAMRGenes = self.get_unitigs().get_selected_genes()
        unitigGenesOfInterest = self.get_unitigsOfInterest()
        # initialise the unitig figure
        fig, ax = self.initialise_plots(len(unitigGenesOfInterest))
        # get minimum length of all genes
        minLength = min(min(list(geneLengths.values())))
        # iterate through the unitigs
        for unitig in tqdm(unitigGenesOfInterest):
            if not len(unitigGenesOfInterest[unitig]) == 0:
                # get the ID for this unitig hash
                unitigId = unitig_mapping[unitig]
                # get a readable list of genes in this path
                listOfGenes = self.convert_paths_to_genes(unitigGenesOfInterest[unitig])
                # get the lengths of genes on this unitig
                unitiglengths = self.get_minimum_gene_length(geneLengths,
                                                        listOfGenes)
                # we need to keep track of the y values in the plot
                height = []
                # iterate through the genes on this unitig
                for gene in range(len(listOfGenes)):
                    # decide what colour we will make this gene in the figure
                    geneColor = self.get_gene_colour(listOfGenes[gene][1:],
                                                    allAMRGenes)
                    # get the coordinate for this gene
                    x, y, dx, dy = self.get_gene_coordinates(listOfGenes[gene][0],
                                                            height,
                                                            unitiglengths[gene],
                                                            unitigId)
                    # add an arrow for each gene
                    ax.arrow(x,
                            y,
                            dx,
                            dy,
                            width = 0.015,
                            head_width = 0.13,
                            head_length=50,
                            color = geneColor,
                            length_includes_head=True)
                    height.append(unitiglengths[gene])
        plt.yticks([i for i in range(1, len(unitigGenesOfInterest) + 1)])
        plotFilename = os.path.join(output_dir, "context_plots") + ".pdf"
        fig.set_size_inches(10, (len(unitigGenesOfInterest)/3))
        fig.savefig(plotFilename, dpi=400)

def parse_fastq_lines(fh):
    # Initialize a counter to keep track of the current line number
    line_number = 0
    # Iterate over the lines in the file
    for line in fh:
        # Increment the line number
        line_number += 1
        # If the line number is divisible by 4, it's a sequence identifier line
        if line_number % 4 == 1:
            # Extract the identifier from the line
            identifier = line.split(" ")[0][1:]
        # If the line number is divisible by 4, it's a sequence line
        elif line_number % 4 == 2:
            sequence = line.strip()
        elif line_number % 4 == 0:
            # Yield the identifier, sequence and quality
            yield identifier, sequence, line.strip()

def parse_fastq(fastq_file):
    # Initialize an empty dictionary to store the results
    results = {}
    # Open the fastq file
    with gzip.open(fastq_file, 'rt') as fh:
        # Iterate over the lines in the file
        for identifier, sequence, quality in parse_fastq_lines(fh):
            # Add the identifier and sequence to the results dictionary
            results[identifier] = {"sequence": sequence,
                                "quality": quality}
    # Return the dictionary of results
    return results

def write_fastq(fastq_file,
                data):
    # Open the fastq file
    with gzip.open(fastq_file, 'wt') as fh:
        # Iterate over the data
        for identifier, value in data.items():
            # Write the identifier line
            fh.write(f'@{identifier}\n')
            # Write the sequence line
            fh.write(f'{value["sequence"]}\n')
            # Write the placeholder quality lines
            fh.write('+\n')
            fh.write(f'{value["quality"]}\n')