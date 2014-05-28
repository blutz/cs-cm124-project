# baseline read mapper with insertions
# USAGE: ....py REFERENCE_FILE READS_FILE
# Output (to stdout): Chromosome #, insertion sequence, insertion location

import sys

print "Hello"

if len(sys.argv < 3):
    sys.stderr.write("Wrong number of arguments\n")
    sys.exit()
