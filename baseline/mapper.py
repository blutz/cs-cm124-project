# baseline read mapper with insertions
# USAGE: ....py REFERENCE_FILE READS_FILE
# Output (to stdout): Chromosome #, insertion sequence, insertion location

# Expecting reads to be length 50

import sys

MAX_INSERTION_LENGTH = 5

def find_insertions(read, reference):
    print "Testing read " + read
    # First find where it needs to be aligned
    align_index = None
    align_mismatches = 0
    for i in range(len(reference)-len(read)+(2*MAX_INSERTION_LENGTH)):
        print "Checking index " + str(i)
        mismatches = [0, 0, 0]
        insertions = 0
        bad_exit = False # Gets set when the loop breaks without finding a match
        if i > 8004:
            sys.exit()
        for j in range(len(read)-1):
            # We've reached the end of the read
            if j+sum(mismatches) >= len(read):
                print "1"
                break
            # We found a matching character
            print " " + read[j] + " " + str(insertions) + " " + str(mismatches)
            if read[j] == reference[j+i-sum(mismatches)]:
                # Advance the insertion number if needed
                if mismatches[insertions] > 0:
                    insertions = insertions + 1
            # We found a mismatch
            else:
                if insertions > 1:
                    bad_exit = True
                    print "2"
                    break
                mismatches[insertions] = mismatches[insertions] + 1
                if mismatches[insertions] > MAX_INSERTION_LENGTH:
                    bad_exit = True
                    print "3"
                    break
        if insertions < 2 or (insertions == 2 and mismatches[2] == 0):
            # We've found a match with 2 or fewer insertions
            if not align_index and not bad_exit:
                align_index = i
                align_mismatches = sum(mismatches)
                print "NEW ALIGNMENT AT " + str(i) + " WITH " + str(sum(mismatches)) + " MISMATCHES"
            elif sum(mismatches) < align_mismatches and not bad_exit:
                align_index = i
                align_mismatches = sum(mismatches)
    print read
    print "FINAL ALIGNMENT AT " + str(align_index) + " WITH " + str(align_mismatches) + " MISMATCHES"
    sys.exit()
                

def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Wrong number of arguments\n")
        sys.exit()

    chromosome_number = 0
    reference_genome = ""

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


    # Now reference_genome is the genome and reads is an array of all reads

    reads = []
    with open(sys.argv[2], "r") as reads_file:
        for line in reads_file:
            if line[0] not in {"A", "G", "T", "C"}:
                continue
            # If we reach this point, line is a read sequence
            reads.extend(line.strip().split(","))

    insertions = []
    for read in reads:
        read_insertions = find_insertions(read, reference_genome)
        if read_insertions:
            insertions.extend(read_insertions)

    print sorted(insertions, key=lambda ins: ins[2])


if __name__ == "__main__":
    main()
