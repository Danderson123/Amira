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
    def get_forward_node_from_node(self,
                                sourceNode: Node,
                                allAMRHashes) -> list:
        """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
        # get the list of forward edge hashes for this node
        nodeForwardEdges = sourceNode.get_forward_edge_hashes()
        if len(nodeForwardEdges) > 0:
            targetNodes = []
            targetNodeDirections = []
            # iterate through the edge hashes
            for edge_hash in nodeForwardEdges:
                # get the edge object corresponding to this edge hash
                edge = self.get_graph().get_edges()[edge_hash]
                # get the target node for this edge
                targetNode = edge.get_targetNode()
                # check if the target node contains an AMR gene
                if targetNode.__hash__() in allAMRHashes:
                    targetNodes.append(edge.get_targetNode())
                    # get the direction we are going into the target node
                    targetNodeDirections.append(edge.get_targetNodeDirection())
            return targetNodes, targetNodeDirections
        else:
            return None, None
    def get_backward_node_from_node(self,
                                sourceNode: Node,
                                allAMRHashes) -> list:
        """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
        # get the list of forward edge hashes for this node
        nodeBackwardEdges = sourceNode.get_backward_edge_hashes()
        if len(nodeBackwardEdges) > 0:
            targetNodes = []
            targetNodeDirections = []
            # iterate through the edge hashes
            for edge_hash in nodeBackwardEdges:
                # get the edge object corresponding to this edge hash
                edge = self.get_graph().get_edges()[edge_hash]
                # get the target node for this edge
                targetNode = edge.get_targetNode()
                # check if the target node contains an AMR gene
                if targetNode.__hash__() in allAMRHashes:
                    targetNodes.append(edge.get_targetNode())
                    # get the direction we are going into the target node
                    targetNodeDirections.append(edge.get_targetNodeDirection())
            return targetNodes, targetNodeDirections
        else:
            return None, None
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
        # get nodes that are anchors for traversing in the forward direction of the graph
        for nodeHash in AMRNodes:
            # get the forward edges for this node
            forward_edges = [self.get_graph().get_edge_by_hash(h) for h in AMRNodes[nodeHash].get_forward_edge_hashes()]
            # get the backward edges for this node
            backward_edges = [self.get_graph().get_edge_by_hash(h) for h in AMRNodes[nodeHash].get_backward_edge_hashes()]
            # get the forward edges that contain an AMR gene
            forwardAMRNodes = [e.get_targetNode() for e in forward_edges if e.get_targetNode().__hash__() in AMRNodes]
            # get the backward edges that contain an AMR gene
            backwardAMRNodes = [e.get_targetNode() for e in backward_edges if e.get_targetNode().__hash__() in AMRNodes]
            # if the number of backward neighbors is not 1 then this is a forward anchor
            if not len(forwardAMRNodes) == 1:
                nodeAnchors.add(nodeHash)
            if not len(backwardAMRNodes) == 1:
                nodeAnchors.add(nodeHash)
        return nodeAnchors
    def get_node_geneMer(self,
                        node,
                        nodeDirection):
        if nodeDirection == 1:
            return node.get_canonical_geneMer()
        else:
            return node.get_geneMer().get_rc_geneMer()
    def check_if_node_in_anchors(self,
                                nodeHash,
                                nodeAnchors):
        if nodeHash in nodeAnchors:
            return False
        else:
            return True
    def get_path_start_and_ends(self,
                            complex_paths):
            starts = {}
            ends = {}
            for pathHash in complex_paths:
                path = complex_paths[pathHash]
                # let's get the start and end node of the path
                start_node = tuple(path[0])
                end_node = tuple(path[-1])
                # add the starts and to a dictionary
                if not start_node in starts:
                    starts[start_node] = []
                starts[start_node].append(pathHash)
                if not end_node in ends:
                    ends[end_node] = []
                ends[end_node].append(pathHash)
            return starts, ends
    def get_path_pairs_to_merge(self,
                            starts,
                            ends):
        pathsToMerge = []
        for nodeHash in ends:
            if nodeHash in starts:
                for firstPath in ends[nodeHash]:
                    for secondPath in starts[nodeHash]:
                        if not firstPath == secondPath:
                            pathsToMerge.append((firstPath, secondPath))
        return pathsToMerge
    def merge_paths(self,
                    pathsToMerge,
                    complex_paths,
                    allPathReads):
        allMergedPaths = {}
        allMergedPathReads = {}
        seenPathHashes = set()
        pathId = 1
        for p in pathsToMerge:
            # get the first path to merge
            firstPath = complex_paths[p[0]]
            # get the first path reads
            firstPathReads = allPathReads[p[0]]
            # get the second path to merge
            secondPath = complex_paths[p[1]]
            # get the second path reads
            secondPathReads = allPathReads[p[1]]
            # check to see that there is at least 1 read supporting merging these paths
            if not set(firstPathReads).isdisjoint(secondPathReads):
                # merge the paths
                mergedPath = firstPath + secondPath[1:]
                # merge the path reads
                mergedReads = list(set(firstPathReads + secondPathReads))
                # add the merged path information to a dictionary
                allMergedPaths[pathId] = mergedPath
                # add the merged read information to a dictionary
                allMergedPathReads[pathId] = mergedReads
                # keep track of which paths we have merged
                seenPathHashes.add(p[0])
                seenPathHashes.add(p[1])
                # increment the path ID
                pathId += 1
        for pathHash in complex_paths:
            # retrieve information for paths that have not been merged
            if not pathHash in seenPathHashes:
                # store unmerged paths in the path dict
                allMergedPaths[pathId] = complex_paths[pathHash]
                # store unmerged path reads in the read dict
                allMergedPathReads[pathId] = allPathReads[pathHash]
                seenPathHashes.add(pathHash)
                # increment path ID
                pathId += 1
        return allMergedPaths, allMergedPathReads
    def simplify_paths(self,
                    complex_paths,
                    all_path_reads):
        """ returns a dictionary of merged & enumerated paths """
        starts, ends = self.get_path_start_and_ends(complex_paths)
        pathsToMerge = self.get_path_pairs_to_merge(starts,
                                                    ends)
        mergedPaths, mergedReads = self.merge_paths(pathsToMerge,
                                        complex_paths,
                                        all_path_reads)
        return mergedPaths, mergedReads
    def get_unitigs_of_interest(self):
        # isolate nodes containing AMR genes
        AMRNodes = self.get_all_nodes_containing_AMR_genes()
        # get AMR nodes that can act as anchors
        nodeAnchors = self.get_AMR_anchors(AMRNodes)
        # get the forward path for all forward anchors until we hit a backward anchor
        complex_paths = {}
        simple_paths = {}
        all_path_reads = {}
        # get a set of all node hashes we are interested in
        allAMRHashes = set(AMRNodes.keys())
        for nodeHash in tqdm(list(nodeAnchors)):
            # get forward paths
            fw_node_list, fw_node_direction_list = self.get_forward_node_from_node(AMRNodes[nodeHash],
                                                                                allAMRHashes)
            bw_node_list, bw_node_direction_list = self.get_backward_node_from_node(AMRNodes[nodeHash],
                                                                                allAMRHashes)
            if fw_node_list and bw_node_list:
                fw_node_list += bw_node_list
                fw_node_direction_list += bw_node_direction_list
            if not fw_node_list and bw_node_list:
                fw_node_list = bw_node_list
                fw_node_direction_list = bw_node_direction_list
            if fw_node_list:
                for n in range(len(fw_node_list)):
                    next_node = fw_node_list[n]
                    next_node_direction = fw_node_direction_list[n]
                    next_node_geneMer = self.get_node_geneMer(next_node,
                                                        next_node_direction)
                    path = [AMRNodes[nodeHash].get_canonical_geneMer(), next_node_geneMer]
                    pathNodeHashes = [AMRNodes[nodeHash].__hash__(), next_node.__hash__()]
                    pathReads = set([r for r in next_node.get_reads()])#set([r for r in AMRNodes[nodeHash].get_reads()] + [r for r in next_node.get_reads()])
                    pathNodes = [AMRNodes[nodeHash], next_node]
                    # check if a node is an anchor, if so get the path
                    extend = self.check_if_node_in_anchors(next_node.__hash__(),
                                                        nodeAnchors)
                    while extend:
                        if next_node_direction == 1:
                            next_node_list, next_node_direction_list = self.get_forward_node_from_node(next_node,
                                                                                                    allAMRHashes)
                        else:
                            next_node_list, next_node_direction_list = self.get_backward_node_from_node(next_node,
                                                                                                    allAMRHashes)
                        next_node = next_node_list[0]
                        next_node_direction = next_node_direction_list[0]
                        next_node_geneMer = self.get_node_geneMer(next_node,
                                                            next_node_direction)
                        path.append(next_node_geneMer)
                        pathNodes.append(next_node)
                        pathNodeHashes.append(next_node.__hash__())
                        # if the next node is an anchor we end the extension
                        extend = self.check_if_node_in_anchors(next_node.__hash__(),
                                                        nodeAnchors)
                        if extend:
                            # add the reads to the read set
                            for read in next_node.get_reads():
                                pathReads.add(read)
                    # get a hash for the path. This is the same for the forward and reverse path
                    pathNodeHashes = sorted([pathNodeHashes, list(reversed(pathNodeHashes))])[0]
                    pathHash = hash(tuple(pathNodeHashes))
                    # we can say a path is simple or complicated depending on whether it includes a branching node at either end
                    if not (self.get_graph().get_degree(pathNodes[0]) > 2 or self.get_graph().get_degree(pathNodes[-1]) > 2):
                        path_dict = simple_paths
                    else:
                        path_dict = complex_paths
                    # add a list of genes and their strands to the relevant path dictionary
                    path_dict[pathHash] = path
                    # initialise the dictionary of reads per path
                    if not pathHash in all_path_reads:
                        all_path_reads[pathHash] = list(pathReads)
            else:
                if self.get_graph().get_degree(AMRNodes[nodeHash]) == 0:
                    # this is a singleton
                    pathHash = hash(AMRNodes[nodeHash].__hash__())
                    simple_paths[pathHash] = [AMRNodes[nodeHash].get_canonical_geneMer()]#[convert_int_strand_to_string(gene.get_strand()) + gene.get_name() for gene in AMRNodes[nodeHash].get_canonical_geneMer()]
                    all_path_reads[pathHash] = [r for r in AMRNodes[nodeHash].get_reads()]
        complex_paths.update(simple_paths)
        # if a path is complex due to a junction, we can make a unitig longer by enumerating all possible paths through the junction
        mergedPaths, mergedPathReads = self.simplify_paths(complex_paths,
                                                        all_path_reads)
        return simple_paths, all_path_reads

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
                readableGenes = self.convert_paths_to_genes(unitigsOfInterest[path])
                # write out the list of genes
                with open(os.path.join(output_dir, str(pathId), "annotated_genes.txt"), "w") as outGenes:
                    outGenes.write("\n".join(readableGenes))
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