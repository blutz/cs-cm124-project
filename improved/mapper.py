# baseline read mapper with insertions
# USAGE: ....py REFERENCE_FILE READS_FILE
# Output (to stdout): Chromosome #, insertion sequence, insertion location

# Expecting reads to be length 50

import sys
import resource
import datetime

MAX_INSERTION_LENGTH = 5
MAX_NUM_INSERTIONS = 2
SUFFIX_TREE_SIZE = 10
READ_LENGTH = 50
READ_DIVISIONS = READ_LENGTH/SUFFIX_TREE_SIZE

DEBUG=True

# Returns str/pos inserts needed to make fragment equal to the
# start of text
def min_insertions(fragment, text, end=False):
    if fragment == "" or text == "":
        return ([], True)
    if end:
        fragment = fragment[::-1]
        text = text[::-1]

    i = 0
    inserts = []
    inserts_cur = 0
    inserts_len = 0
    while i < len(fragment):
        if (i-inserts_len) >= len(text):
            return([], False)
        # The boring case -- the text matches
        if text[i-inserts_len] == fragment[i]:
            if len(inserts) > inserts_cur:
                inserts_cur = len(inserts)
        else:
            # We need to insert a new placeholder into inserts
            if len(inserts) <= inserts_cur:
                inserts.append({"pos":i, "str":""})
            inserts[inserts_cur]['str'] += fragment[i]
            inserts_len += 1
        i += 1

    if end:
        for ins in inserts:
            ins["str"] = ins["str"][::-1]
            ins["pos"] = len(fragment) - ins["pos"] - len(ins["str"])

    valid = True
    for ins in inserts:
        if len(ins["str"]) > 5:
            valid = False

    return inserts, valid

# Returns a list of dictionaries with entries "pos" (where the insertion is) and "str"
# (what the insertion is)
def find_insertions(read, reference_genome, stree):
    print read
    ret_size = 10
    ret_val = []
    # We're going to go through each of the five sections of the read
    # and look them up in the suffix tree until we get a match
    for i in range(READ_DIVISIONS):
        print read[i*SUFFIX_TREE_SIZE:(i+1)*SUFFIX_TREE_SIZE]
        suffixes = stree[read[i*SUFFIX_TREE_SIZE:(i+1)*SUFFIX_TREE_SIZE]]
        for suffix in suffixes:
            print suffix
            inserts_left, valid_left = min_insertions(read[:i*SUFFIX_TREE_SIZE], reference_genome[:suffix], True)
            print inserts_left
            print valid_left
            inserts_right, valid_right = min_insertions(read[(i+1)*SUFFIX_TREE_SIZE:], reference_genome[suffix+SUFFIX_TREE_SIZE:], False)
            print inserts_right
            print valid_right
            print "-----"
            if (((len(inserts_left) + len(inserts_right)) <= 2) and valid_left and
                valid_right and ((len(inserts_left) + len(inserts_right) < ret_size) or
                ret_val == None)):
                for ins in inserts_left:
                    ins["pos"] += i
                for ins in inserts_right:
                    ins["pos"] += i
        print "****"

    return ret_val

def all_gene_combos(togo):
    # Base case
    if togo == 1:
        return ['A', 'G', 'C', 'T']

    a = []
    g = []
    c = []
    t = []
    for b in all_gene_combos(togo-1):
        a.append("A"+b)
        g.append("G"+b)
        c.append("C"+b)
        t.append("T"+b)

    return a + g + c + t

def build_suffix_tree(genome):
    stree = {}
    # Build tree framework
    for b in all_gene_combos(10):
        stree[b] = []
    # Now actually put data in the tree
    for i in range(len(genome)-1-SUFFIX_TREE_SIZE):
        stree[genome[i:i+SUFFIX_TREE_SIZE]].append(i)
    if DEBUG and False:
        for a, b in stree.iteritems():
            print str(a) + ": "
            print " " + str(b)
    return stree

# Main setup to read in the files
def main():
    d1 = datetime.datetime.now()
    if len(sys.argv) < 3:
        sys.stderr.write("Wrong number of arguments\n")
        sys.exit()

    chromosome_number = 0
    reference_genome = ""
    
    if DEBUG:
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

    d2 = datetime.datetime.now()
    stree = build_suffix_tree(reference_genome)
    if DEBUG:
        print str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576) + "MB"
    d3 = datetime.datetime.now()

    # Now reference_genome is the genome and reads is an array of all reads

    if DEBUG:
        print str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576) + "MB"
    reads = []
    with open(sys.argv[2], "r") as reads_file:
        for line in reads_file:
            if line[0] not in {"A", "G", "T", "C"}:
                continue
            # If we reach this point, line is a read sequence
            reads.extend(line.strip().split(","))

    d4 = datetime.datetime.now()
    di = 0
    if DEBUG:
        print str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576) + "MB"
    insertions = []
    for read in reads:
        read_insertions = find_insertions(read, reference_genome, stree)
        if DEBUG:
            print str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576) + "MB"
        d5 = datetime.datetime.now()
        if di >= 500 and DEBUG:
            print d4-d1
            print d5-d4
            sys.exit()
        if read_insertions:
            insertions.extend(read_insertions)
        if DEBUG:
            print insertions
        di += 1

    for ins in sorted(insertions, key=lambda ins: ins["pos"]):
        print str(chromosome_number)+","+str(ins['seq'])+","+str(ins['pos'])


if __name__ == "__main__":
    main()
