import gzip
import statistics


class Unitig:

    def __init__(self, path, path_node_coverages, unitig_ID):
        self._path = path
        self._min_coverage = min(path_node_coverages)
        self._max_coverage = max(path_node_coverages)
        self._mean_coverage = statistics.mean(path_node_coverages)
        self._unitig_ID = unitig_ID
        self._terminal_nodes = {path[0], path[1]}

    def get_path(self):
        return self._path

    def get_min_coverage(self):
        return self._min_coverage

    def get_max_coverage(self):
        return self._max_coverage

    def get_mean_coverage(self):
        return self._mean_coverage

    def get_unitig_ID(self):
        return self._unitig_ID

    def get_terminal_nodes(self):
        return self._terminal_nodes

    def shares_terminals(self, other_unitig):
        if len(self.get_terminal_nodes().union(other_unitig.get_terminal_nodes())) > 0:
            return 1
        else:
            return 0


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
    if ".gz" in fastq_file:
        try:
            with gzip.open(fastq_file, "rt") as fh:
                # Iterate over the lines in the file
                for identifier, sequence, quality in parse_fastq_lines(fh):
                    # Add the identifier and sequence to the results dictionary
                    results[identifier] = {"sequence": sequence, "quality": quality}
            return results
        except:
            pass
    with open(fastq_file, "r") as fh:
        # Iterate over the lines in the file
        for identifier, sequence, quality in parse_fastq_lines(fh):
            # Add the identifier and sequence to the results dictionary
            results[identifier] = {"sequence": sequence, "quality": quality}
    # Return the dictionary of results
    return results


def write_fastq(fastq_file, data):
    # Open the fastq file
    with gzip.open(fastq_file, "wt") as fh:
        # Iterate over the data
        for identifier, value in data.items():
            # Write the identifier line
            fh.write(f"@{identifier}\n")
            # Write the sequence line
            fh.write(f'{value["sequence"]}\n')
            # Write the placeholder quality lines
            fh.write("+\n")
            fh.write(f'{value["quality"]}\n')
