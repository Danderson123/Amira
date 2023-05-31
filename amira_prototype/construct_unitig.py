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
    def get_reads_per_amr_node(self,
                            AMRNodes):
        # make a dictionary of AMR nodes per read
        amr_gene_reads = {}
        for nodeHash in tqdm(AMRNodes):
            node = AMRNodes[nodeHash]
            amr_gene_reads[nodeHash] = []
            for readId in node.get_reads():
                amr_gene_reads[nodeHash].append(readId)
        return amr_gene_reads
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
    def cluster_reads(self,
                    read_to_genes):
        parent = {}
        for readId in read_to_genes:
            for nodeHash in read_to_genes[readId]:
                parent[nodeHash] = nodeHash
        # Union genes that share a read
        for genes in read_to_genes.values():
            for i in range(len(genes) - 1):
                self.union(genes[i], genes[i+1], parent)
        # Create clusters
        clusters = {}
        for gene, p in parent.items():
            cluster = self.find(p, parent)
            if cluster not in clusters:
                clusters[cluster] = []
            clusters[cluster].append(gene)
        return {i: cluster for i, cluster in enumerate(clusters.values())}
    def get_junctions(self):
        nodeJunctions = set()
        # get nodes that are anchors for traversing in the forward direction of the graph
        for node in self.get_graph().all_nodes():
            # if the degree of this node is more then this is a junction
            if self.get_graph().get_degree(node) > 2:
                nodeJunctions.add(node.__hash__())
        return nodeJunctions
    def assess_resolvability(self,
                            clusters):
        # get nodes that are at junction
        nodeJunctions = self.get_junctions()
        # get the graph
        graph = self.get_graph()
        easy = {}
        intermediate = {}
        difficult = {}
        for c in tqdm(clusters):
            # express all the reads in this cluster as a list of nodes
            clusterNodes = set()
            for readId in clusters[c]:
                clusterNodes.update(graph.get_readNodes()[readId])
            # get the number of junction nodes in this cluster
            junctionCount = len(nodeJunctions.intersection(clusterNodes))
            # if there are 0 junctions, this is easy to resolve
            if junctionCount == 0:
                easy[c] = clusters[c]
            # if there is 1 junction, this is intermediate to resolve
            elif junctionCount == 1:
                intermediate[c] = intermediate[c]
            # if there is more than 1 junction, this is difficult to resolve
            else:
                difficult[c] = difficult[c]
        return easy, intermediate, difficult
    def separate_paralogs(self):
        # isolate nodes containing AMR genes
        AMRNodes = self.get_all_nodes_containing_AMR_genes()
        # make a dictionary of AMR nodes per read
        amr_gene_reads = self.get_reads_per_amr_node(AMRNodes)
        # cluster the reads if they share at least one gene
        clusters = self.cluster_reads(amr_gene_reads)
        # separate the reads into easy, intermediate and difficult to resolve clusters
        easy, intermediate, difficult = self.assess_resolvability(clusters)
        print(len(easy), len(intermediate), len(difficult))
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
        # get the untig genes and reads
        unitigsOfInterest = self.get_unitigsOfInterest()
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