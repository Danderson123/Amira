import gzip
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import os
import subprocess
from tqdm import tqdm

from construct_graph import GeneMerGraph
from construct_gene import convert_int_strand_to_string

class UnitigTools:
    def __init__(self,
                graph: GeneMerGraph,
                listOfGenes: list,
                readFile: str,
                output_dir: str):
        self._graph = graph
        self._listOfGenes = listOfGenes
        self._fastqContent = parse_fastq(readFile)
        self._output_dir = output_dir
    def get_graph(self):
        """ returns the geneMerGraph """
        return self._graph
    def get_selected_genes(self):
        """ returns the list of selected genes """
        return self._listOfGenes
    def get_fastqContent(self):
        """ returns a dictionary of all read data in the input fastq """
        return self._fastqContent
    def get_output_dir(self):
        """ returns a string of the output directory """
        return self._output_dir
    def get_read_data(self,
                    readId):
        """ return a dictionary of the data for this read """
        return self.get_fastqContent()[readId]
    def get_nodes_of_interest(self,
                            geneOfInterest):
        """ extracts the graph nodes containing the genes of interest and returns them as a list """
        return self.get_graph().get_nodes_containing(geneOfInterest)
    def get_all_nodes_containing_AMR_genes(self):
        """ return a dictionary of nodes containing AMR genes """
        AMRNodes = {}
        # iterate through the list of specified genes
        for geneOfInterest in tqdm(self.get_selected_genes()):
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_of_interest(geneOfInterest)
            # add nodes of interest to the AMR node set
            for n in nodesOfInterest:
                AMRNodes[n.__hash__()] = n
        return AMRNodes
    def get_junctions(self):
        nodeJunctions = set()
        # get nodes that are anchors for traversing in the forward direction of the graph
        for node in self.get_graph().all_nodes():
            # if the degree of this node is more then this is a junction
            if self.get_graph().get_degree(node) > 2:
                nodeJunctions.add(node.__hash__())
        return nodeJunctions
    def get_AMR_anchors_and_junctions(self,
                                    AMRNodes):
        # get the graph
        graph = self.get_graph()
        # initialise the node anchor and junction sets
        nodeAnchors = set()
        nodeJunctions = set()
        # get nodes that are anchors for traversing in the forward direction of the graph
        for nodeHash in AMRNodes:
            node = graph.get_node_by_hash(nodeHash)
            if not graph.get_degree(node) > 2:
                # get the forward neighbors for this node
                forward_neighbors = graph.get_forward_neighbors(node)
                # get the backward neighbors for this node
                backward_neighbors = graph.get_backward_neighbors(node)
                # get the forward neighbors that contain an AMR gene
                forwardAMRNodes = [n for n in forward_neighbors if n.__hash__() in AMRNodes]
                # get the backward neighbors that contain an AMR gene
                backwardAMRNodes = [n for n in backward_neighbors if n.__hash__() in AMRNodes]
                if len(backwardAMRNodes) == 0 or len(forwardAMRNodes) == 0:
                    nodeAnchors.add(nodeHash)
            # if the number of backward neighbors is not 1 or 0 then this is a forward anchor
            else:
                nodeJunctions.add(nodeHash)
        return nodeAnchors, nodeJunctions
    def resolve_one_junction_clusters(self,
                                    reads_to_separate,
                                    clusterJunctions,
                                    graph):
        """ separate reads in intermediate clusters based on the path they follow through junctions """
        paths = {}
        for r in range(len(reads_to_separate)):
            nodeHashes = graph.get_readNodes()[reads_to_separate[r]]
            # if the only node on a read is a junction node it is not useful for us
            if not len(nodeHashes) == 1:
                indicesOfJunctions = [i for i, x in enumerate(nodeHashes) if x in clusterJunctions[r]]
                adjancentNodes = []
                for i in indicesOfJunctions:
                    if not i == 0:
                        adjancentNodes.append(nodeHashes[i-1])
                    if not i == len(nodeHashes) - 1:
                        adjancentNodes.append(nodeHashes[i+1])
                paths[reads_to_separate[r]] = adjancentNodes
        # a cluster with one junction will have 2 possible paths
        subclusters = {}
        for readId in paths:
            paths[readId].sort()
            pathTuple = tuple(paths[readId])
            if not pathTuple in subclusters:
                subclusters[pathTuple] = []
            subclusters[pathTuple].append(readId)
        # there is only 1 junction so two potential paths through it
        return [v for v in subclusters.values()]
    def assess_resolvability(self,
                            clusters,
                            amrJunctions):
        # get the graph
        graph = self.get_graph()
        easy = {}
        intermediate = {}
        difficult = {}
        clusterId = 1
        for c in clusters:
            # express all the reads in this cluster as a list of nodes
            clusterJunctions = []
            for readId in clusters[c]:
                clusterJunctions.append([n for n in graph.get_readNodes()[readId] if n in amrJunctions])
            # get the number of junction nodes in this cluster
            junctionCount = max([len(j) for j in clusterJunctions])
            clusters[c] = list(clusters[c])
            # if there are 0 junctions, this is easy to resolve
            if junctionCount == 0:
                easy[clusterId] = clusters[c]
            # if there is 1 junction, this is intermediate to resolve
            elif junctionCount == 1 or junctionCount == 2:
                intermediate[clusterId] = clusters[c]
            # if there is more than 1 junction, this is difficult to resolve
            else:
                intermediate[clusterId] = clusters[c]
            clusterId += 1
        return easy, intermediate, difficult, clusterId
    def write_subset_of_reads(self,
                        output_dir,
                        readIds):
        # make the output directory
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        # subset the reads for this unitig
        subsettedReadData = {}
        for r in readIds:
            subsettedReadData[r] = self.get_read_data(r)
        # write the per unitig fastq data
        readFileName = os.path.join(output_dir, "cluster_reads.fastq.gz")
        write_fastq(readFileName,
                    subsettedReadData)
        return readFileName
    def resolve_easy_clusters(self,
                        easyClusters):
        # initialise a list of read file paths
        readFiles = []
        # get the graph
        graph = self.get_graph()
        # iterate through the easy to resolve clusters
        for cluster in tqdm(easyClusters):
            # get all the node hashes for this cluster
            clusterReads = set()
            for readId in easyClusters[cluster]:
                readNodes = [graph.get_node_by_hash(h) for h in graph.get_readNodes()[readId]]
                for node in readNodes:
                    for read in node.get_reads():
                        clusterReads.add(read)
            # write a fastq of the reads in this path
            outputDir = os.path.join(self.get_output_dir(), str(cluster))
            readFileName = self.write_subset_of_reads(outputDir,
                                                clusterReads)
            # keep track of the read file
            readFiles.append(readFileName)
        return readFiles
    def bin_reads_by_path(self,
                        clusterReadNodes,
                        pathsOfInterest):
        # bin the reference reads based on which paths they follow
        subclusteredReads = {}
        for readId in clusterReadNodes:
            for path in pathsOfInterest:
                #if self.contains_sublist(clusterReadNodes[readId],
                #                        list(path)):
                if all(p in clusterReadNodes[readId] for p in list(path)):
                    if not path in subclusteredReads:
                        subclusteredReads[path] = []
                    subclusteredReads[path].append(readId)
        return subclusteredReads
    def get_all_reads_in_path(self,
                            uniqueNodes,
                            graph):
        # get the reads for each of these useful nodes
        clusterReads = set()
        for nodeHash in uniqueNodes:
            node = graph.get_node_by_hash(nodeHash)
            for readId in node.get_reads():
                clusterReads.add(readId)
        return clusterReads
    def resolve_one_junction(self,
                            junction,
                            clusterReadNodes,
                            graph):
        pathsOfInterest = set()
        # get the node for the hash
        node = graph.get_node_by_hash(junction)
        # get the neighbors for the junction
        fw_neighbors = graph.get_forward_neighbors(node)
        bw_neighbors = graph.get_backward_neighbors(node)
        # we want to separate reads based on the node immediately next to the junction
        if len(fw_neighbors) == 2:
            next_nodes = fw_neighbors
        if len(bw_neighbors) == 2:
            next_nodes = bw_neighbors
        for next_node in next_nodes:
            pathsOfInterest.add((junction, next_node.__hash__()))
        # split the reference reads into subclusters based on the path they follow
        subclusteredReferenceReads = self.bin_reads_by_path(clusterReadNodes,
                                                        pathsOfInterest)
        # get the nodes that are unique to each path
        subclusterAllReads = {}
        for subcluster in subclusteredReferenceReads:
            uniqueNodes = set()
            for readId in subclusteredReferenceReads[subcluster]:
                readNodes = clusterReadNodes[readId]
                # get the indices of the nodes in the path
                junctionIndex = readNodes.index(subcluster[0])
                adjacentNodeIndex = readNodes.index(subcluster[1])
                assert not (junctionIndex == -1 or adjacentNodeIndex == -1)
                # trim the read nodes so we are only adding nodes that occur after the junction
                if junctionIndex > adjacentNodeIndex:
                    trimmedReadNodes = readNodes[:junctionIndex]
                else:
                    trimmedReadNodes = readNodes[adjacentNodeIndex:]
                uniqueNodes.update(trimmedReadNodes)
            # get the reads for each of these useful nodes
            subclusterAllReads[subcluster] = self.get_all_reads_in_path(uniqueNodes,
                                                                    graph)
        return subclusterAllReads
    def resolve_two_junctions(self,
                            junctions,
                            clusterReadNodes,
                            graph):
        pathsOfInterest = set()
        # if the junctions are adjacent, we need to use 4 nodes to resolve them
        nodes = [graph.get_node_by_hash(h) for h in junctions]
        # get a list of the neighbors that occur immediately after each junction
        neighbors = [sorted([graph.get_forward_neighbors(junctionNode),
                        graph.get_backward_neighbors(junctionNode)], key=len, reverse=True)[0] for junctionNode in nodes]
        # now we get all the possible paths we can see through the junctions
        for j1_neighbor in neighbors[0]:
            for j2_neighbor in neighbors[1]:
                pathsOfInterest.add((j1_neighbor.__hash__(), j2_neighbor.__hash__()))
        # split the reference reads into subclusters based on the path they follow
        subclusteredReferenceReads = self.bin_reads_by_path(clusterReadNodes,
                                                    pathsOfInterest)
        # if the junction nodes are not immediately adjacent
        if not graph.check_if_nodes_are_adjacent(graph.get_node_by_hash(junctions[0]),
                                                                    graph.get_node_by_hash(junctions[1])):
            # the junctions are not immediately adjacent
            junctionsAdjacent = False
        else:
            # the junctions are immediately adjacent
            junctionsAdjacent = True
        # get the nodes that are unique to each path
        junctionSet = set(junctions)
        subclusterAllReads = {}
        for subcluster in subclusteredReferenceReads:
            uniqueNodes = set()
            for readId in subclusteredReferenceReads[subcluster]:
                readNodes = clusterReadNodes[readId]
                # get the indices of the junctions in the path
                firstJunctionIndex = readNodes.index(junctions[0])
                secondJunctionIndex = readNodes.index(junctions[1])
                assert not (firstJunctionIndex == -1 or secondJunctionIndex == -1)
                if not junctionsAdjacent:
                    # trim the read nodes so we are only adding nodes that occur after the junction
                    if firstJunctionIndex > secondJunctionIndex:
                        trimmedReadNodes = readNodes[secondJunctionIndex + 1: firstJunctionIndex]
                    else:
                        trimmedReadNodes = readNodes[firstJunctionIndex + 1: secondJunctionIndex]
                else:
                    # trim the read nodes so that we are only adding nodes outside of the junctions
                    trimmedReadNodes = [n for n in readNodes if not n in junctionSet]
                uniqueNodes.update(trimmedReadNodes)
            # get the reads for each of these useful nodes
            subclusterAllReads[subcluster] = self.get_all_reads_in_path(uniqueNodes,
                                                                    graph)
        return subclusterAllReads
    def get_unitigs_on_read(self,
                        nodesOnRead,
                        anchorsAndJnctions):
        output = []
        chunk = []
        for element in nodesOnRead:
            if element in anchorsAndJnctions:
                if len(nodesOnRead) == 1:
                    chunk.append(element)  # add current delimiter to chunk
                    output.append(chunk)
                    continue
                if chunk:  # if the chunk already started, add it to the output
                    chunk.append(element)  # add current delimiter to chunk
                    output.append(chunk)
                chunk = [element]  # start a new chunk
            elif chunk:  # if we are inside a chunk, add the element to it
                chunk.append(element)
        return output
    def resolve_intermediate_clusters(self,
                                    intermediateClusters,
                                    amrJunctions,
                                    amrAnchors,
                                    clusterId):
        # initialise a list of read file paths
        readFiles = []
        # get the graph
        graph = self.get_graph()
        # join the junction and anchor sets
        junctionsAndAnchors = amrJunctions.union(amrAnchors)
        # iterate through the intermediate clusters
        for cluster in tqdm(intermediateClusters):
            readsAsUnitigs = {}
            unitigMapping = {}
            unitigId = 0
            for readId in intermediateClusters[cluster]:
                # get the read node hashes
                readNodeHashes = graph.get_readNodes()[readId]
                # split the readNodes into chunks
                chunked = self.get_unitigs_on_read(readNodeHashes,
                                                junctionsAndAnchors)
                readsAsUnitigs[readId] = []
                for c in chunked:
                    tup = tuple(c)
                    rev_tup = tuple(list(reversed(c)))
                    if tup in unitigMapping:
                        readsAsUnitigs[readId].append(unitigMapping[tup])
                    elif rev_tup in unitigMapping:
                        readsAsUnitigs[readId].append(unitigMapping[rev_tup])
                    else:
                        unitigMapping[tup] = unitigId
                        readsAsUnitigs[readId].append(unitigMapping[tup])
                        unitigId += 1
            print(sorted([v for v in readsAsUnitigs.values()], key=len, reverse=True))
    def old_resolve_intermediate_clusters(self,
                                    intermediateClusters,
                                    amrJunctions,
                                    clusterId):
        # initialise a list of read file paths
        readFiles = []
        # get the graph
        graph = self.get_graph()
        # iterate through the intermediate clusters
        for cluster in tqdm(intermediateClusters):
            clusterReadNodes = {}
            junctions = set()
            for readId in intermediateClusters[cluster]:
                # get the read node hashes
                readNodeHashes = graph.get_readNodes()[readId]
                # store the nodes on these reads
                clusterReadNodes[readId] = readNodeHashes
                print(readNodeHashes)
                # we also want to keep track of all the junctions that we see in the reads
                junctions.update([h for h in readNodeHashes if h in amrJunctions])
            end
            junctions = list(junctions)
            assert len(junctions) < 3
            # if there is only 1 junction we can resolve using 1 node
            if len(junctions) == 1:
                subclusterAllReads = self.resolve_one_junction(junctions[0],
                                                        clusterReadNodes,
                                                        graph)
            else:
                subclusterAllReads = self.resolve_two_junctions(junctions,
                                                            clusterReadNodes,
                                                            graph)
            # iterate through the subclusters
            for subcluster in subclusterAllReads:
                # write a fastq of the reads in this path
                outputDir = os.path.join(self.get_output_dir(), str(clusterId))
                readFileName = self.write_subset_of_reads(outputDir,
                                                    subclusterAllReads[subcluster])
                # keep track of the read file
                readFiles.append(readFileName)
                clusterId += 1
        return readFiles
    def find(self,
            x,
            parent):
        if x not in parent:
            parent[x] = x
        elif parent[x] != x:
            parent[x] = self.find(parent[x], parent)
        return parent[x]
    def union(self,
            x,
            y,
            parent):
        parent[self.find(x, parent)] = self.find(y, parent)
    def cluster_anchor_reads(self,
                            anchorReads):
        parent = {}
        nodeReads = {}
        for readId, nodeHashes in anchorReads.items():
            for nodeHash in nodeHashes:
                parent.setdefault(nodeHash, nodeHash)
                nodeReads.setdefault(nodeHash, set()).add(readId)
        # Union genes that share a read
        for nodeHashes in anchorReads.values():
            for i in range(len(nodeHashes) - 1):
                self.union(nodeHashes[i], nodeHashes[i+1], parent)
        # Create clusters
        clusters = {}
        for nodeHash, p in parent.items():
            root = self.find(p, parent)
            if root not in clusters:
                clusters[root] = set()
            clusters[root].update(nodeReads[nodeHash])
        return {i+1: cluster for i, cluster in enumerate(clusters.values())}
    def get_anchoring_reads(self,
                        nodeAnchors):
        """ get reads that contain AMR genes and that contain at least 2 AMR anchors """
        # get the graph
        graph = self.get_graph()
        anchorReads = {}
        for nodeHash in nodeAnchors:
            node = graph.get_node_by_hash(nodeHash)
            for readId in node.get_reads():
                if not readId in anchorReads:
                    anchorReads[readId] = []
                anchorReads[readId].append(nodeHash)
        # get reads that have at least 2 anchors
        to_delete = []
        for readId in anchorReads:
            if len(anchorReads[readId]) < 2:
                to_delete.append(readId)
        # remove reads with fewer than 2 anchors
        for r in to_delete:
            del anchorReads[r]
        return anchorReads
    def contains_sublist(self,
                        allNodes,
                        subNodes):
        n = len(subNodes)
        if any((subNodes == allNodes[i:i+n] or list(reversed(subNodes)) == allNodes[i:i+n]) for i in range(len(allNodes)-n+1)):
            return True
        else:
            return False
    def separate_paralogs(self):
        # isolate nodes containing AMR genes
        AMRNodes = self.get_all_nodes_containing_AMR_genes()
        # get the AMR nodes that are anchors and nodes
        amrAnchors, amrJunctions = self.get_AMR_anchors_and_junctions(AMRNodes)
        # make a dictionary of AMR anchor nodes per read
        anchorReads = self.get_anchoring_reads(amrAnchors)
        # cluster the anchor reads if they have at least 2 anchors share at least one anchor
        anchor_clusters = self.cluster_anchor_reads(anchorReads)
        # separate the reads into easy, intermediate and difficult to resolve clusters
        easy, intermediate, difficult, clusterId = self.assess_resolvability(anchor_clusters,
                                                                            amrJunctions)
        # resolve the easy clusters
        readFiles = self.resolve_easy_clusters(easy)
        # resolve the intermediate clusters
        readFiles += self.resolve_intermediate_clusters(intermediate,
                                                        amrJunctions,
                                                        amrAnchors,
                                                        clusterId)
        return readFiles
    def convert_paths_to_genes(self,
                            pathGeneMers):
        geneMerGenes = [convert_int_strand_to_string(gene.get_strand()) + gene.get_name() for gene in pathGeneMers[0]]
        for geneMer in pathGeneMers[1:]:
            gene_strand = convert_int_strand_to_string(geneMer[-1].get_strand())
            gene_name = geneMer[-1].get_name()
            geneMerGenes.append(gene_strand + gene_name)
        return geneMerGenes
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
                if os.path.exists(os.path.join(os.path.dirname(inputFastq), "flye_output", "assembly.fasta")):
                    # map the reads to the consensus file
                    map_command = "minimap2 -a --MD -t 1 "
                    map_command += os.path.join(os.path.dirname(inputFastq), "flye_output", "assembly.fasta") + " " + inputFastq
                    map_command += " > " + os.path.join(os.path.dirname(inputFastq), "reads_mapped_to_consensus.sam")
                    subprocess.run(map_command, shell=True, check=True)
                    # polish the pandora consensus
                    racon_command = "panRG_building_tools/racon/build/bin/racon" + " -t 1 "
                    racon_command += inputFastq + " " + os.path.join(os.path.dirname(inputFastq), "reads_mapped_to_consensus.sam") + " "
                    racon_command += os.path.join(os.path.dirname(inputFastq), "flye_output", "assembly.fasta") + " "
                    racon_command += "> " + os.path.join(os.path.dirname(inputFastq), "racon_polished_assembly.fasta")
                    subprocess.run(racon_command, shell=True, check=True)
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