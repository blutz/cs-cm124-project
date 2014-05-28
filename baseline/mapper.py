# baseline read mapper with insertions
# USAGE: ....py REFERENCE_FILE READS_FILE
# Output (to stdout): Chromosome #, insertion sequence, insertion location

# Expecting reads to be length 50

import sys
import resource
from copy import deepcopy

MAX_INSERTION_LENGTH = 5

def find_insertions(read, reference):
    # First find where it needs to be aligned
    align_index = None
    align_mismatches = 0
    align_inserts = [{"seq":"","pos":None},{"seq":"","pos":None}]
    for i in range(len(reference)-len(read)+(2*MAX_INSERTION_LENGTH)):
        mismatches = [0, 0, 0]
        insertions = 0
        insertion_seq = [{"seq":"","pos":None},{"seq":"","pos":None}]
        bad_exit = False # Gets set when the loop breaks without finding a match
        for j in range(len(read)-1):
            # We've reached the end of the read
            if j+sum(mismatches) >= len(read):
                break
            # We found a matching character
            if read[j] == reference[j+i-sum(mismatches)]:
                # Advance the insertion number if needed
                if mismatches[insertions] > 0:
                    insertions = insertions + 1
            # We found a mismatch
            else:
                if insertions > 1:
                    bad_exit = True
                    break
                mismatches[insertions] = mismatches[insertions] + 1
                # Save this mismatch
                insertion_seq[insertions]["seq"]+=read[j]
                if not insertion_seq[insertions]["pos"]:
                    insertion_seq[insertions]["pos"] = i+j-sum(mismatches)
                if mismatches[insertions] > MAX_INSERTION_LENGTH:
                    bad_exit = True
                    break
        if insertions < 2 or (insertions == 2 and mismatches[2] == 0):
            # We've found a match with 2 or fewer insertions
            if not align_index and not bad_exit:
                align_index = i
                align_mismatches = sum(mismatches)
                align_inserts = deepcopy(insertion_seq)
            elif sum(mismatches) < align_mismatches and not bad_exit:
                align_index = i
                align_mismatches = sum(mismatches)
                align_inserts = deepcopy(insertion_seq)
    if align_inserts[1]['pos'] == None:
        align_inserts.pop(1)
    if align_inserts[0]['pos'] == None:
        align_inserts.pop(0)

    return align_inserts
                

def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Wrong number of arguments\n")
        sys.exit()

    chromosome_number = 0
    reference_genome = ""

    print str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576) + "MB"

    # These will throw errors if the file doesn't exist. Good enough for the baseline
    with open(sys.argv[1], "r") as ref_file:
        for line in ref_file:
            # This line is a base sequence
            if line[0] in {"A", "G", "T", "C"}:
                reference_genome = reference_genome + line.strip()
                continue
            # This line is info about the file
            elif line[0] == ">":
                if line[1:3] == "chr":
                    chromosome_number = line[4]
                continue
            # This file is not structured properly
            else:
                sys.stderr.write("Improper reference file format\n")
                sys.exit()

    print str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576) + "MB"

    # Now reference_genome is the genome and reads is an array of all reads

    print str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576) + "MB"
    reads = []
    with open(sys.argv[2], "r") as reads_file:
        for line in reads_file:
            if line[0] not in {"A", "G", "T", "C"}:
                continue
            # If we reach this point, line is a read sequence
            reads.extend(line.strip().split(","))

    print str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576) + "MB"
    insertions = []
    for read in reads:
        print str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576) + "MB"
        sys.exit()
        read_insertions = find_insertions(read, reference_genome)
        if read_insertions:
            insertions.extend(read_insertions)
        print insertions

    for ins in sorted(insertions, key=lambda ins: ins["pos"]):
        print str(chromosome_number)+","+str(ins['seq'])+","+str(ins['pos'])


if __name__ == "__main__":
    main()
