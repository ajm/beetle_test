import sys
import os
import glob
import collections
import re
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment


def get_msa_species(fname) :
    beetle_species = set()
    tc_genes = []

    #print >> sys.stderr, "\nreading %s ... " % fname

    for s in SeqIO.parse(fname, 'fasta') :
        header = s.description.split()
        species = header[1].split('=')[1]
        geneid = header[2].split('=')[1]

        if species == 'tribolium_castaneum' :
            tc_genes.append(geneid)
            continue

        beetle_species.add(species)

    #print >> sys.stderr, "  contains %d beetle species, %d tribolium castaneum genes" % (len(beetle_species), len(tc_genes))

    return list(beetle_species),tc_genes


def remove_tribolium(fname) :

    msa = AlignIO.read(open(fname), "fasta")
    newmsa = MultipleSeqAlignment([])
    for record in msa :

	header = record.description.split()
	species = header[1].split('=')[1]
	geneid = header[2].split('=')[1]
	if species != 'tribolium_castaneum' :
	    newmsa.append(record)

    return newmsa


def write_triboliumless_file(msa, fname) :
    f = open(fname, "w")
    AlignIO.write(msa, f, "fasta")
    f.close()


def main() :

    if len(sys.argv) != 3 :
        print >> sys.stderr, "Usage: %s <input.fasta> <output.fasta>" % sys.argv[0]
        sys.exit(1)

    inputfname = sys.argv[1]
    outputfname = sys.argv[2]

    msa = remove_tribolium(inputfname)
    write_triboliumless_file(msa, outputfname)



if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by User..."
        sys.exit(1)


