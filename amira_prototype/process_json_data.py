import json

class ReadData:

    def __init__(self,
                readName,
                readId,
                geneIdList=None,
                geneStrandList=None):
        self.readName = readName
        self.readId = readId
        self.geneIdList = geneIdList
        self.geneStrandList = geneStrandList

    def get_readName(self):
        return self.readName
    def get_readId(self):
        return self.readId
    def get_geneIdList(self):
        return self.geneIdList
    def get_geneStrandList(self):
        return self.geneStrandList

class GeneData:

    def __init__(self,
            geneName,
            geneId,
            allReadIds=[],
            allReadStrands=[]):
        self.geneName = geneName
        self.geneId = geneId
        self.allReadIds = allReadIds
        self.allReadStrands = allReadStrands

    def get_geneName(self):
        return self.geneName
    def get_geneId(self):
        return self.geneId
    def get_readIdList(self):
        return self.allReadIds
    def get_readStrands(self):
        return self.allReadStrands

def load_json_data(json_path):
    """import the annotated read data as a dictionary"""
    with open(json_path, "r") as i:
        annotated_reads = json.loads(i.read())
    return annotated_reads

def initialise_processing_data_structures():
    """initialise the read and gene data structures"""
    return {}, {}, 1, []

def get_read_strand(gene):
    """determine the strand of the gene on the read and remove strand from gene string"""
    if gene[0] == "+":
        strand = +1
    if gene[0] == "-":
        strand = -1
    gene = gene[1:]
    return gene, strand

def determine_geneId(gene,
                    geneId,
                    geneNameIdDict,
                    geneIdNameDict):
    if not gene in geneNameIdDict:
        this_geneId = geneId
        geneNameIdDict[gene] = this_geneId
        geneIdNameDict[this_geneId] = gene
        geneId += 1
    else:
        this_geneId = geneNameIdDict[gene]
    return geneId, this_geneId, geneNameIdDict, geneIdNameDict

def append_geneData(gene,
                    this_geneId,
                    readId,
                    geneData,
                    strand):
    if len(geneData) == this_geneId - 1:
        this_gene = GeneData(gene,
                            this_geneId,
                            [readId],
                            [strand])
        geneData.append(this_gene)
    else:
        this_gene = geneData[this_geneId - 1]
        this_gene.allReadIds.append(readId)
        this_gene.allReadStrands.append(strand)
        geneData[this_geneId - 1] = this_gene
    return geneData

def extract_geneData(geneId,
                    geneNameIdDict,
                    geneIdNameDict,
                    geneData,
                    readDict,
                    read,
                    readId):
    geneIdList = []
    geneStrandList = []
    for gene in readDict[read]:
        gene, strand = get_read_strand(gene)
        geneId, this_geneId, geneNameIdDict, geneIdNameDict = determine_geneId(gene,
                                                                            geneId,
                                                                            geneNameIdDict,
                                                                            geneIdNameDict)
        geneData = append_geneData(gene,
                                this_geneId,
                                readId,
                                geneData,
                                strand)
        geneIdList.append(this_geneId)
        geneStrandList.append(strand)
    return geneId, geneIdList, geneStrandList, geneNameIdDict, geneIdNameDict, geneData

def extract_readData(read,
                    readId,
                    geneIdList,
                    geneStrandList,
                    readData):
    this_read = ReadData(read,
                        readId,
                        geneIdList,
                        geneStrandList)
    readData.append(this_read)
    readId += 1
    return readId, readData

def isolate_reads(readDict,
                refine):
    #refinedReadDict = {}
    if refine:
        with open(refine, "r") as i:
            required_genes = i.read().splitlines()
        refined_genes = set()
        for g in required_genes:
            refined_genes.add("+" + g)
            refined_genes.add("-" + g)
    else:
        refined_genes = set()
        for read in readDict:
            for a in readDict[read]:
                refined_genes.add(a)
    #for read in readDict:
     #   if any(a in refined_genes for a in readDict[read]):
      #      refinedReadDict[read] = readDict[read]
    return refined_genes

def classify_read_json(readDict,
                    refine):
    readNameIdDict, readIdNameDict, readId, readData = initialise_processing_data_structures()
    geneNameIdDict, geneIdNameDict, geneId, geneData = initialise_processing_data_structures()
    refined_genes = isolate_reads(readDict,
                                refine)
    for read in readDict:
        readNameIdDict[read] = readId
        readIdNameDict[readId] = read
        geneId, geneIdList, geneStrandList, geneNameIdDict, geneIdNameDict, geneData = extract_geneData(geneId,
                                                                                                        geneNameIdDict,
                                                                                                        geneIdNameDict,
                                                                                                        geneData,
                                                                                                        readDict,
                                                                                                        read,
                                                                                                        readId)
        readId, readData = extract_readData(read,
                                            readId,
                                            geneIdList,
                                            geneStrandList,
                                            readData)
    interesting_genes = set()
    for gene in list(refined_genes):
        if gene[1:] in geneNameIdDict:
            interesting_genes.add(geneNameIdDict[gene[1:]])
    return readData, readIdNameDict, geneData, geneIdNameDict, list(interesting_genes)