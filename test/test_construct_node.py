import sys
sys.path.insert(0, "../amira_prototype")

from construct_node import Node
from construct_read import Read

def test_node_increment_coverage():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
    read = Read("read1",
                genes)
    geneMers = [x for x in read.get_geneMers(3)]
    for g in geneMers:
        node = Node(g)
        for i in range(5):
            coverage = node.increment_coverage()
            assert coverage == i + 1

def test_add_read():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
    read1 = Read("read1",
                genes)
    read2 = Read("read2",
                genes)
    read3 = Read("read3",
                genes)
    geneMers = [x for x in read1.get_geneMers(3)]
    for g in geneMers:
        node = Node(g)
        readList = [read1, read1, read2, read3, read3]
        for read in readList:
            node.add_read(read)
        assert all(r in node.get_reads() for r in readList), "A read is missing from the node"
        assert all([x for x in node.get_reads()].count(r) == 1 for r in readList), "A read occurs more than once in the read list"

def test_remove_read():
    genes = ["+gene1", "-gene2", "+gene3", "-gene4", "+gene5", "-gene6"]
    read1 = Read("read1",
                genes)
    read2 = Read("read2",
                genes)
    read3 = Read("read3",
                genes)
    geneMers = [x for x in read1.get_geneMers(3)]
    # test removing a read that is in the list
    for g in geneMers:
        node = Node(g)
        readList = [read1, read1, read2, read3, read3]
        for read in readList:
            node.add_read(read)
        for read in readList:
            node.remove_read(read)
            assert not read in node.get_reads()
    # test removing a read that is not in the list
    read4 = Read("read4",
                genes)
    for g in geneMers:
        node = Node(g)
        readList = [read1, read1, read2, read3, read3]
        for read in readList:
            node.add_read(read)
    node.remove_read(read4)
    assert all(r in node.get_reads() for r in readList), "A read was incorrectly removed from the Node"
    assert all([x for x in node.get_reads()].count(r) == 1 for r in readList), "A read occurs more than once in the read list"


sys.stderr.write("Testing construct_node: Node.increment_coverage\n")
test_node_increment_coverage()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_node: Node.add_read\n")
test_add_read()
sys.stderr.write("Test passed\n")
sys.stderr.write("Testing construct_node: Node.remove_read\n")
test_remove_read()
sys.stderr.write("Test passed\n")
