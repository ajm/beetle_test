import sys
import os
import glob
import collections
import re
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment


def remove_gaps(fname) :

    msa = AlignIO.read(open(fname), "fasta")

    for record in msa :
	record.seq = record.seq.ungap("-")
	#print record.seq

    return msa


def write_gapless_file(msa, fname) :
    f = open(fname, "w")
    AlignIO.write(msa, f, "fasta")
    f.close()

def main() :

    if len(sys.argv) != 3 :
        print >> sys.stderr, "Usage: %s <input.fasta> <output.fasta>" % sys.argv[0]
        sys.exit(1)

    inputfname = sys.argv[1]
    outputfname = sys.argv[2]

    msa = remove_gaps(inputfname)
    write_gapless_file(msa, outputfname)
    


if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by User..."
        sys.exit(1)

