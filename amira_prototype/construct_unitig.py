import gzip
import matplotlib.pyplot as plt
import os
import subprocess

from construct_graph import GeneMerGraph
from construct_node import Node
from construct_gene_mer import define_rc_geneMer

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
        if self.get_graph().get_degree(node) == 2 or self.get_graph().get_degree(node) == 1:
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
    def get_unitigs_of_interest(self):
        """ returns a dictionary of genes of interest and their linear paths in the graph """
        unitigsOfInterest = {}
        # iterate through the list of specified genes
        for geneOfInterest in self.get_selected_genes():
            # get the graph nodes containing this gene
            nodesOfInterest = self.get_nodes_of_interest(geneOfInterest)
            all_unitigs = []
            all_unitig_reads = []
            # iterate through the nodes containing this gene
            for node in nodesOfInterest:
                # get the linear path for this node
                node_unitig, unitig_reads = self.get_unitig_for_node(node)
                if node_unitig:
                    reversed_node_unitig = list(reversed([define_rc_geneMer(n) for n in node_unitig]))
                    if not (node_unitig in all_unitigs or reversed_node_unitig in all_unitigs):
                        all_unitigs.append(node_unitig)
                        all_unitig_reads.append(unitig_reads)
            # populate a dictionary with the gene name and the corresponding unitigs
            unitigsOfInterest[geneOfInterest] = {"unitigs": all_unitigs,
                                                "reads": all_unitig_reads}
        return unitigsOfInterest

class UnitigBuilder:

    def __init__(self,
                unitigsOfInterest):
        self._unitigsOfInterest = unitigsOfInterest
    def get_unitigs(self) -> dict:
        """ returns all unitigs containing the specified AMR genes """
        return self._unitigsOfInterest
    def separate_reads(self,
                output_dir,
                readFile):
        fastqContent = parse_fastq(readFile)
        readFiles = []
        for geneOfInterest in self.get_unitigs():
            if not os.path.exists(os.path.join(output_dir, geneOfInterest)):
                os.mkdir(os.path.join(output_dir, geneOfInterest))
            for i in range(len(self.get_unitigs()[geneOfInterest]["unitigs"])):
                unitigReads = self.get_unitigs()[geneOfInterest]["reads"][i]
                subsettedReadData = {}
                for u in unitigReads:
                    for r in u:
                        subsettedReadData[r] =  fastqContent[r]
                # write the per unitig fastq data
                filename = os.path.join(output_dir, geneOfInterest, geneOfInterest + "_" + str(i) + ".fastq.gz")
                write_fastq(filename,
                            subsettedReadData)
                readFiles.append(filename)
        return readFiles
    def run_raven(self,
                inputFastq,
                raven_path):
        outputConsensus = inputFastq.replace(".fastq.gz", ".fasta")
        raven_command = " ".join([raven_path,
                                "-t",
                                "8",
                                inputFastq,
                                ">",
                                outputConsensus])
        subprocess.run(raven_command, shell=True, check=True)
    def initialise_plots(self,
                        unitigCount,
                        paralog):
        plt.rc('font', size=6)
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
        ax.set_ylabel(paralog)
        ax.set_xlabel("Unitig length (bp)")
        return fig, ax
    def get_geneMer_unionAll(self,
                            unitigGenes):
        unionAll = unitigGenes[0]
        for geneMer in unitigGenes[1:]:
            unionAll.append(geneMer[-1])
        return unionAll
    def get_gene_colour(self,
                        currentGene,
                        geneOfInterest):
        if currentGene == geneOfInterest:
            color = "r"
        else:
            color = "g"
        return color
    def get_gene_coordinates(self,
                            geneStrand,
                            previousCoordinates,
                            geneLength,
                            paralogCount):
        if geneStrand == 1:
            x = sum(previousCoordinates)
            y = paralogCount
            dx = geneLength
            dy = 0
        else:
            x = sum(previousCoordinates) + geneLength
            y = paralogCount
            dx = geneLength * -1
            dy = 0
        return x, y, dx, dy
    def get_minimum_gene_length(self,
                                geneLengths,
                                unionAll):
        unitiglengths = [max(geneLengths[g.get_name()]) for g in unionAll]
        return unitiglengths, min(unitiglengths)
    def visualise_unitigs(self,
                        geneLengths: dict,
                        output_dir: str):
        """ generate a figure to visualise the genes, order of genes, direction and lengths of genes on unitigs containing AMR genes """
        for paralog in self.get_unitigs():
            fig, ax = self.initialise_plots(len(self.get_unitigs()[paralog]["unitigs"]),
                                            paralog)
            count = 1
            for unitig in self.get_unitigs()[paralog]["unitigs"]:
                height = []
                unionAll = self.get_geneMer_unionAll(unitig)
                unitiglengths, minLength = self.get_minimum_gene_length(geneLengths,
                                                                unionAll)
                for gene in range(len(unionAll)):
                    geneLength = unitiglengths[gene]
                    geneColor = self.get_gene_colour(unionAll[gene].get_name(),
                                                    paralog)
                    x, y, dx, dy = self.get_gene_coordinates(unionAll[gene].get_strand(),
                                                            height,
                                                            geneLength,
                                                            count)
                    ax.arrow(x,
                            y,
                            dx,
                            dy,
                            width = ((len(self.get_unitigs()[paralog]["unitigs"]) + 1) / 1000),
                            head_width = ((len(self.get_unitigs()[paralog]["unitigs"]) + 1) / 200),
                            head_length=minLength,
                            color = geneColor,
                            length_includes_head=True)
                    height.append(geneLength)
                count +=1
            plt.yticks([i for i in range(1, count)])
            plotFilename = os.path.join(output_dir, paralog, paralog) + ".pdf"
            plt.savefig(plotFilename)

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