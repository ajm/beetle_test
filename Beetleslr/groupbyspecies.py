import sys
import os
import glob
import collections
import re
import shutil
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment


def get_msa_species(fname) :
    beetle_species = set()

    for s in SeqIO.parse(fname, 'fasta') :
        header = s.description.split()
        species = header[1].split('=')[1]

        if species == 'tribolium_castaneum' :
            continue

        beetle_species.add(species)

    return len(beetle_species)


def main() :

    if len(sys.argv) != 3 :
        print >> sys.stderr, "Usage: %s <input.fasta> <outputfolder>" % sys.argv[0]
        sys.exit(1)

    fname = sys.argv[1]
    outputfolder = sys.argv[2]
    species = get_msa_species(fname)
    newdestination = '%s/%s' % (outputfolder, str(species))
    shutil.move(fname, newdestination)


if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by User..."
        sys.exit(1)
