from Bio import SeqIO
from sys import stderr, exit, argv
from glob import glob
from os.path import join, isfile, splitext, basename, dirname
from os import system, getcwd
import os

from Bio import AlignIO


raxml_seed = 37

slr_config = """seqfile: %s
treefile: %s
outfile: %s
"""

paml_config = """seqfile = %s
treefile = %s
outfile = %s

noisy = 9
verbose = 0
runmode = 0 * user defined tree

seqtype = 1 * codons
CodonFreq = 0 * equal

model = 0 * one omega for all branches

NSsites = 0 * one omega ratio

icode = 0 

fix_kappa = 0
kappa = 2

fix_omega = 0
omega = 0.1
"""

def print_out_details(fname, species) :
    if species :
        print fname, len(species)

def count_species_fasta(fname, ref_species) :
    aligned_species = set()

    for s in SeqIO.parse(fname, 'fasta') :
        species = s.description.split()[1].split('=')[1]
        if species == ref_species :
            continue
            
        aligned_species.add(species)

    return len(aligned_species)

def rm_f(files) :
    for f in files :
        try :
            os.remove(f)
        except :
            pass

def build_seqfile(fasta_fname, out_type='phylip-sequential') :
    phylip_fname = _new_ext(fasta_fname, 'phy')

    fin = open(fasta_fname)
    fout = open(phylip_fname, 'w')

    alignments = AlignIO.parse(fin, 'fasta')
    AlignIO.write(alignments, fout, out_type)

    fin.close()
    fout.close()

    print "\twritten %s" % phylip_fname

    return phylip_fname

def _new_ext(fname, ext) :
    return splitext(fname)[0] + '.' + ext

def build_tree(phylip_in) :
    global raxml_seed
    
    out = splitext(phylip_in)[0]
    command = "raxmlHPC-SSE3 -m GTRGAMMA -# 10 -p %d -s %s -n %s" % (raxml_seed, phylip_in, out)
    command += " > /dev/null 2> /dev/null"

    print "\trunning raxml..."

    if system(command) != 0 :
        print >> stderr, "FAIL!\n" + " ".join(command.split())
        exit(1)

    #resultfile = 'RAxML_result.' + out
    resultfile = 'RAxML_bestTree.' + out

    if not isfile(resultfile) :
        print >> stderr, "%s not found" % resultfile
        exit(1)

    newresultfile = _new_ext(phylip_in, 'tree')
    os.rename(resultfile, newresultfile)

    rm_f(glob("RAxML_*"))

    print "\tgot tree file %s" % newresultfile

    return newresultfile

def build_config(program, seqfile, treefile) :
    global paml_config
    global slr_config

    configfile = _new_ext(seqfile, 'config')
    resultfile = _new_ext(seqfile, program)

    config = slr_config if program == 'slr' else paml_config

    with open(configfile, 'w') as f :
        f.write(config % (seqfile, treefile, resultfile))

    print "\twritten %s, expect results in %s" % (configfile, resultfile)

    return configfile, resultfile

def paml(configfile, resultfile) :
    command = "codeml %s" % configfile
    command +=  " > /dev/null 2> /dev/null"

    print "\trunning paml..."

    if system(command) != 0 :
        print >> stderr, "FAIL!\n" + " ".join(command.split())
        exit(1)

    if not isfile(resultfile) :
        print >> stderr, "%s not found" % resultfile
        exit(1)

    print "\tgot paml results %s" % resultfile

    return resultfile

def slr(configfile, resultfile) :
    command = "Slr %s" % configfile
    command += " > /dev/null 2> /dev/null"

    print "\trunning slr..."

    if system(command) != 0 :
        print >> stderr, "FAIL!\n" + " ".join(command.split())
        exit(1)

    if not isfile(resultfile) :
        print >> stderr, "%s not found" % resultfile
        exit(1)

    print "\tgot slr results %s" % resultfile

    return resultfile

def run_programs(fname) :
    print "running " + fname
    
    seqfile = build_seqfile(fname)
    treefile = build_tree(seqfile)

    configfile,resultfile = build_config('slr', seqfile, treefile)
    slr(configfile, resultfile) 
    
    configfile,resultfile = build_config('paml', seqfile, treefile)
    paml(configfile, resultfile)

def main() :
    
    if len(argv) != 2 :
        print >> stderr, "Usage: %s <msa>" % argv[0]
        exit(1)

    msa = argv[1]

    msa_dir = dirname(msa)
    fname = basename(msa)

    # fucking statisticians can't code...
    pwd = getcwd()
    os.chdir(msa_dir)

    run_programs(fname)

    os.chdir(pwd)

    return 0

if __name__ == '__main__' :
    try :
        exit(main())
    except KeyboardInterrupt :
        print >> stderr, "Killed by user..."
        exit(1)

