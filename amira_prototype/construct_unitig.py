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
                                sourceNode: Node) -> list:
        """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
        # get the list of forward edge hashes for this node
        nodeForwardEdges = sourceNode.get_forward_edge_hashes()
        if len(nodeForwardEdges) > 0:
            # iterate through the edge hashes
            for edge_hash in nodeForwardEdges:
                # get the edge object corresponding to this edge hash
                edge = self.get_graph().get_edges()[edge_hash]
                # get the target node for this edge
                targetNode = edge.get_targetNode()
                # get the degree of the target node
                targetNodeDegree = self.get_graph().get_degree(targetNode)
                # if the degree of the target node is 1 or 2 we can extend the linear path to the next node
                if targetNodeDegree == 2 or targetNodeDegree == 1:
                    # get the direction we are going into the target node
                    targetNodeDirection = edge.get_targetNodeDirection()
                    return True, targetNode, targetNodeDirection
                # else we cannot extend the linear path
                else:
                    return False, None, None
        else:
            return False, None, None
    def get_backward_node_from_node(self,
                                sourceNode: Node) -> list:
        """ returns a list of nodes in the forward direction from this node until a branch or end of unitig is reached """
        # get the list of forward edge hashes for this node
        nodeBackwardEdges = sourceNode.get_backward_edge_hashes()
        if len(nodeBackwardEdges) > 0:
            # iterate through the edge hashes
            for edge_hash in nodeBackwardEdges:
                # get the edge object corresponding to this edge hash
                edge = self.get_graph().get_edges()[edge_hash]
                # get the target node for this edge
                targetNode = edge.get_targetNode()
                # get the degree of the target node
                targetNodeDegree = self.get_graph().get_degree(targetNode)
                # if the degree of the target node is 1 or 2 we can extend the linear path to the next node
                if targetNodeDegree == 2 or targetNodeDegree == 1:
                    # get the direction we are going into the target node
                    targetNodeDirection = edge.get_targetNodeDirection()
                    return True, targetNode, targetNodeDirection
                else:
                    # else we cannot extend the linear path
                    return False, None, None
        else:
            return False, None, None
    def get_forward_path_from_node(self,
                                node):
        forward_nodes_from_node = []
        forward_reads = []
        # get the next node in the forward direction
        forwardExtend, forwardNode, forwardNodeDirection = self.get_forward_node_from_node(node)
        # if we are extending further in the forward direction, get the next canonical gene mer
        while forwardExtend:
            forward_reads.append([r for r in forwardNode.get_reads()])
            # if we enter the next node in the forward direction, we get the next forward node
            if forwardNodeDirection == 1:
                forward_nodes_from_node.append(forwardNode.get_canonical_geneMer())
                forwardExtend, forwardNode, forwardNodeDirection = self.get_forward_node_from_node(forwardNode)
            # if we enter the next node in the backward direction, we get the next backward node
            else:
                forward_nodes_from_node.append(forwardNode.get_geneMer().get_rc_geneMer())
                forwardExtend, forwardNode, forwardNodeDirection = self.get_backward_node_from_node(forwardNode)
        return forward_nodes_from_node, forward_reads
    def get_backward_path_from_node(self,
                                    node):
        backward_nodes_from_node = []
        backward_reads = []
        # get the next node in the backward direction
        backwardExtend, backwardNode, backwardNodeDirection = self.get_backward_node_from_node(node)
        # if we are extending further in the backward direction, get the next canonical gene mer
        while backwardExtend:
            backward_reads.insert(0, [r for r in backwardNode.get_reads()])
            if backwardNodeDirection == -1:
                backward_nodes_from_node.insert(0, backwardNode.get_geneMer().get_canonical_geneMer())
                backwardExtend, backwardNode, backwardNodeDirection = self.get_backward_node_from_node(backwardNode)
            # if we enter the next node in the forward direction, we get the next forward node
            else:
                backward_nodes_from_node.insert(0, backwardNode.get_geneMer().get_rc_geneMer())
                backwardExtend, backwardNode, backwardNodeDirection = self.get_forward_node_from_node(backwardNode)
        return backward_nodes_from_node, backward_reads
    def get_unitig_for_node(self,
                            node):
        """ builds a unitig starting from the node of interest and expanding in both directions """
        if self.get_graph().get_degree(node) == 2 or self.get_graph().get_degree(node) == 1 or self.get_graph().get_degree(node) == 0:
            # get the forward nodes from this node
            forward_nodes_from_node, forward_reads = self.get_forward_path_from_node(node)
            # get the backward nodes from this node
            backward_nodes_from_node, backward_reads = self.get_backward_path_from_node(node)
            # join the backward nodes, this node, and the forward nodes to get the path of nodes
            unitig = backward_nodes_from_node + [node.get_canonical_geneMer()] + forward_nodes_from_node
            unitig_reads = backward_reads + [[r for r in node.get_reads()]] + forward_reads
            return unitig, unitig_reads
        else:
            return None, None
    def hash_unitig(self,
                    unitigGeneMers: list,
                    reverseUnitigGeneMers: list) -> int:
        """ return the hash of a unitig to check if we have seen it before """
        geneMerHashes = [hash(tuple([gene.__hash__() for gene in geneMer])) for geneMer in unitigGeneMers]
        reverseGeneMerHashes = [hash(tuple([gene.__hash__() for gene in geneMer])) for geneMer in reverseUnitigGeneMers]
        untigHash = sorted([hash(tuple(geneMerHashes)), hash(tuple(reverseGeneMerHashes))])[0]
        return untigHash
    def get_unitigs_of_interest(self):
        """ returns a dictionary of unitigs containing AMR genes of interest and a dictionary of reads per unitig"""
        unitigGenesOfInterest = {}
        unitigsReadsOfInterest = {}
        # iterate through the list of specified genes
        for geneOfInterest in self.get_selected_genes():
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_of_interest(geneOfInterest)
            # iterate through the nodes containing this gene
            for node in nodesOfInterest:
                # get the linear path for this node
                node_unitig, unitig_reads = self.get_unitig_for_node(node)
                if node_unitig:
                    reversed_node_unitig = list(reversed([define_rc_geneMer(n) for n in node_unitig]))
                    untigHash = self.hash_unitig(node_unitig,
                                                reversed_node_unitig)
                    if not untigHash in unitigGenesOfInterest:
                        unitigGenesOfInterest[untigHash] = []
                        unitigsReadsOfInterest[untigHash] = []
                    # add the read to the reads for this unitig if it doesn't exist already
                    for geneMerReads in unitig_reads:
                        for read in geneMerReads:
                            if not read in unitigsReadsOfInterest[untigHash]:
                                unitigsReadsOfInterest[untigHash].append(read)
                    # add the gene and it's strand to the genes for this unitig
                    unitigGenesOfInterest[untigHash] = [convert_int_strand_to_string(gene.get_strand()) + gene.get_name() for gene in node_unitig[0]]
                    for geneMer in node_unitig[1:]:
                        gene_strand = convert_int_strand_to_string(geneMer[-1].get_strand())
                        gene_name = geneMer[-1].get_name()
                        unitigGenesOfInterest[untigHash].append(gene_strand + gene_name)
        return unitigGenesOfInterest, unitigsReadsOfInterest

class UnitigTools:
    def __init__(self,
                graph,
                genesOfInterest):
        self._unitigs = Unitigs(graph,
                            genesOfInterest)
        self._unitigGenesOfInterest, self._unitigsReadsOfInterest = self.get_unitigs().get_unitigs_of_interest()
    def get_unitigs(self) -> Unitigs:
        """ returns the Unitigs object containing the specified AMR genes """
        return self._unitigs
    def get_unitigGenesOfInterest(self) -> dict:
        """ returns all unitigs containing the specified AMR genes """
        return self._unitigGenesOfInterest
    def get_unitigReadsOfInterest(self) -> dict:
        """ returns all unitigs containing the specified AMR genes """
        return self._unitigsReadsOfInterest
    def separate_reads(self,
                output_dir,
                readFile):
        fastqContent = parse_fastq(readFile)
        readFiles = []
        # get the untig genes and reads
        unitigGenesOfInterest = self.get_unitigGenesOfInterest()
        unitigsReadsOfInterest = self.get_unitigReadsOfInterest()
        # assign an integer id to each unitig
        unitig_count = 1
        # store the unitig hash to integer mappings
        unitig_mapping = {}
        # iterate through the unitigs
        for unitig in tqdm(unitigGenesOfInterest):
            if not len(unitigGenesOfInterest[unitig]) == 0:
                # make the output directory
                if not os.path.exists(os.path.join(output_dir, str(unitig_count))):
                    os.mkdir(os.path.join(output_dir, str(unitig_count)))
                # write out the list of genes
                with open(os.path.join(output_dir, str(unitig_count), "annotated_genes.txt"), "w") as outGenes:
                    outGenes.write("\n".join(unitigGenesOfInterest[unitig]))
                # subset the reads for this unitig
                unitigReads = unitigsReadsOfInterest[unitig]
                subsettedReadData = {}
                for r in unitigReads:
                    subsettedReadData[r] =  fastqContent[r]
                # write the per unitig fastq data
                readFileName = os.path.join(output_dir, str(unitig_count), str(unitig_count) + ".fastq.gz")
                write_fastq(readFileName,
                            subsettedReadData)
                readFiles.append(readFileName)
                # update the unitig IDs
                unitig_mapping[unitig] = unitig_count
                unitig_count += 1
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
        unitigGenesOfInterest = self.get_unitigGenesOfInterest()
        # initialise the unitig figure
        fig, ax = self.initialise_plots(len(unitigGenesOfInterest))
        # get minimum length of all genes
        minLength = min(min(list(geneLengths.values())))
        # iterate through the unitigs
        for unitig in tqdm(unitigGenesOfInterest):
            if not len(unitigGenesOfInterest[unitig]) == 0:
                # get the ID for this unitig hash
                unitigId = unitig_mapping[unitig]
                # get the list of genes on this unitig
                listOfGenes = unitigGenesOfInterest[unitig]
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