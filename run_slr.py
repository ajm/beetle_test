from Bio import SeqIO
from sys import stderr, exit, argv
from glob import glob
from os.path import join, isfile, splitext, basename, dirname
from os import system, getcwd
import os
from Bio import Phylo
import tempfile
import pexpect

from Bio import AlignIO

RUN_RAXML = False

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

    alignments = AlignIO.read(fin, 'fasta')

    AlignIO.write(alignments, fout, out_type)

    fin.close()
    fout.close()

    #print "\twritten %s" % phylip_fname

    return phylip_fname

def _new_ext(fname, ext) :
    return splitext(fname)[0] + '.' + ext

def build_tree(phylip_in) :
    global raxml_seed
    
    out = splitext(phylip_in)[0]
    command = "raxmlHPC-SSE3 -m GTRGAMMA -# 10 -p %d -s %s -n %s" % (raxml_seed, phylip_in, out)
    #command += " > /dev/null 2> /dev/null"

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

    print "\texpect results in %s" % (resultfile.replace('_', ''))

    return configfile, resultfile

def paml(configfile, resultfile) :
    command = "codeml %s" % configfile
    #command +=  " > /dev/null 2> /dev/null"

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

    print "\trunning Slr..."

    command = "Slr %s" % configfile
    #command += " > /dev/null 2> /dev/null"
    search = ''

    op = pexpect.spawn(command, timeout=None)
    op.maxsize = 1
    originalname = configfile.replace("_.config", "")

    for line in op :
       	lst = line.split()
    	print line,

    	if len(lst) >= 3 :
            search = lst[0] + lst[1]
            if search == "Initialf:" :
                if lst[2] ==  "inf" :
                    file_cleanup(originalname, "Fail: Initial f was inf.")
                    print "\nInitial f was inf, exiting Slr of %s.\n" % originalname
                    exit()
                elif len(str(lst[2])) > 50 :
                    file_cleanup(originalname, "Fail: Initial f was too big.")
                    print "\nInitial f was way too big, exiting Slr of %s.\n" % originalname
                    exit()
        if len(lst) >= 6 :
            search = lst[0] + lst[1] + lst[2] + lst[3] + lst[4] + lst[5]
            if search == "Alignmentcontainsstopcodons.Cannotcontinue." :
                print "\nAlignment contained stop codons, exiting Slr of %s.\n" % originalname
                file_cleanup(originalname, "Fail: Alignment contained stop codons.")
                exit()

        if len(lst) == 4 :
            #search = lst[1] + lst[2] + lst[3] + lst[4] + lst[5]
            search = lst[1] + lst[2] + lst[3]
            if search == "lnL=inf" :
                print "\nlnL was inf, exiting Slr of %s.\n" % originalname
                file_cleanup(originalname, "Fail: lnL was inf.")
                exit()


    op.close()


    if not isfile(resultfile) :
        print >> stderr, "Got no Slr results for %s" % resultfile.replace("_", "")
        exit(1)

    print "\tgot slr results %s" % resultfile.replace("_", "")

    return resultfile

def generate_new_files(fname) :
    # to get gene names that slr can handle (short enough)
    newfname = fname.replace(".", "_.")

    # generate a fasta file with new ids
    d = {}
    sequences = []

    runningids = 1
    for record in SeqIO.parse(fname, 'fasta') :
        d[record.id] = "flyg%s" % runningids
        record.id = d[record.id]
        record.name = ""
        record.description = ""
        sequences.append(record)
        runningids += 1

    SeqIO.write(sequences, newfname, "fasta")
    

    if not RUN_RAXML :
    # generate a treefile with new ids
        treefile = fname.replace("fasta", "tree")
        newtreefile = newfname.replace("fasta", "tree")

        tree = Phylo.read(treefile, 'newick')

        for node in tree.get_terminals():
            node.name = d[node.name]

        Phylo.write(tree, newtreefile, 'newick')

    return newfname

def run_programs(fname) :
    print "\nrunning " + fname.replace('_', '')
    
    seqfile = build_seqfile(fname)

    if RUN_RAXML :
        treefile = build_tree(seqfile)
    else :
        treefile = fname.replace("fasta", "tree")

    configfile,resultfile = build_config('slr', seqfile, treefile)
    slr(configfile, resultfile) 
    
    #configfile,resultfile = build_config('paml', seqfile, treefile)
    #paml(configfile, resultfile)

def file_cleanup(originalname, message) :
    os.rename('%s_.config' % originalname, '%s.config' % originalname)
    os.remove('%s_.fasta' % originalname)
    os.remove('%s_.tree' % originalname)
    os.rename('%s_.phy' % originalname, '%s.phy' % originalname)
    if os.path.isfile('%s_.slr' % originalname) :
        os.rename('%s_.slr' % originalname, '%s.slr' % originalname)

    statusfile = '%s.stats' % originalname
    f = open(statusfile, 'w')
    message = originalname + "\t" + message + "\n"
    f.write(message)
    f.close()

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

    fname = generate_new_files(fname)

    run_programs(fname)

    originalname = fname.replace("_.fasta", "")
    file_cleanup(originalname, "Success: Slr went fine.")

    os.chdir(pwd)

    return 0

if __name__ == '__main__' :
    try :
        exit(main())
    except KeyboardInterrupt :
        print >> stderr, "Killed by user..."
        exit(1)

