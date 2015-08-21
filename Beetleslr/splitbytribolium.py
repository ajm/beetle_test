import sys
import os
import glob
import collections
import re
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment


def split_msa(fname) :
    msa = AlignIO.read(open(fname), "fasta")
    msalist = []

    specieslist = []
    newmsa = MultipleSeqAlignment([])

    for record in msa :

        header = record.description.split()
        species = header[1].split('=')[1]

        if not 'tribolium_castaneum' in specieslist :
            newmsa.append(record)
	    specieslist.append(species)
	elif species != 'tribolium_castaneum' :
	    newmsa.append(record)
	    specieslist.append(species)
	else :
	    msalist.append(newmsa)
	    newmsa = MultipleSeqAlignment([])
	    newmsa.append(record)
	    specieslist.append(species)

    if newmsa :
	msalist.append(newmsa)

    return msalist


def write_split_files(msalist, fname) :
    counter = 0;

    for msa in msalist :
	counter += 1
	thisfname = '%sv%s.fasta' % (fname, counter)
	f = open(thisfname, "w")
    	AlignIO.write(msa, f, "fasta")
	f.close()


def main() :

    if len(sys.argv) != 3 :
        print >> sys.stderr, "Usage: %s <input.fasta> <output>" % sys.argv[0]
        sys.exit(1)

    inputfname = sys.argv[1]
    output = sys.argv[2]

    msalist = split_msa(inputfname)
    write_split_files(msalist, output)

    #msa = remove_tribolium(inputfname)
    #write_triboliumless_file(msa, outputfname)




if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by User..."
        sys.exit(1)

