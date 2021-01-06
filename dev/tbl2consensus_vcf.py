import pysam
import sys
import pdb

if(len(sys.argv) == 1):
    sys.stderr.write("Needs {gt table} {ct table} {normal ID} {tumor ID}\n")
    exit(1)
