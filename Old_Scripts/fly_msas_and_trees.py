import sys
import urllib2
import urllib
import httplib
import xml.sax.saxutils
import os
import tempfile
import glob
import collections
import json
import re
import time

from Bio import Phylo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

DEBUG = False
the_name = ""

biomart_url = "http://metazoa.ensembl.org/biomart/martservice"

biomart_query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "metazoa_mart_26" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.7" >
            
    <Dataset name = "tcastaneum_eg_gene" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "%s" />
    </Dataset>
</Query>
"""

biomart_query_attributes = (
        "dananassae_eg_gene",
        "derecta_eg_gene",
        "dgrimshawi_eg_gene",
        "dmelanogaster_eg_gene",
        "dmojavensis_eg_gene",
        "dpersimilis_eg_gene",
        "dpseudoobscura_eg_gene",
        "dsechellia_eg_gene",
        "dsimulans_eg_gene",
        "dvirilis_eg_gene",
        "dwillistoni_eg_gene",
        "dyakuba_eg_gene"
    )

restapi_genetree_url = "http://rest.ensemblgenomes.org/genetree/member/symbol/%s/%s?content-type=text/x-phyloxml;sequence=cdna;aligned=1"
restapi_homology_url = "http://rest.ensemblgenomes.org/homology/symbol/%s/%s?content-type=application/json;compara=metazoa;aligned=1;sequence=cdna"

fly_species = (
        "drosophila_ananassae",
        "drosophila_erecta",
        "drosophila_grimshawi",
        "drosophila_melanogaster",
        "drosophila_mojavensis",
        "drosophila_persimilis",
        "drosophila_pseudoobscura",
        "drosophila_sechellia",
        "drosophila_simulans",
        "drosophila_virilis",
        "drosophila_willistoni",
        "drosophila_yakuba"
    )

extra_html_unescapes = {
    "&quot;" : '"',
    "&apos;" : "'"
}

def unescape_html(s) :
    return xml.sax.saxutils.unescape(s, extra_html_unescapes)

def remove_gap_columns(msa) :
    global DEBUG
    global the_name
    tmp = None

    for i in range(0, msa.get_alignment_length(), 3) :
        col1 = msa[:,i]
	col2 = msa[:,i+1]
	col3 = msa[:,i+2]

        if col1.count('-') == len(col1) and col2.count('-') == len(col2) and col3.count('-') == len(col3) :
            continue

        if tmp :
            tmp += msa[:, i:i+3]
        else :
            tmp = msa[:, i:i+3]

    #print >> sys.stderr, "  remove gap columns %dx%d --> %dx%d" % (
    #            len(msa), msa.get_alignment_length(), 
    #            len(tmp), tmp.get_alignment_length())

    if tmp.get_alignment_length() % 3 != 0 :
	print the_name, "The alignment was not dividable by 3, exiting..."
	sys.exit(1)
 
    if DEBUG :
        tmp = tmp[:, :24]
    return tmp

def filter_alignment_from_phyloxml(tree, species_list) :
    msa = MultipleSeqAlignment([])
    flies = set()

    #print >> sys.stderr, "  filtering alignment to fly species..."

    for node in tree.get_terminals() :
        include = False

        for prop in node.properties :
            if (prop.ref == 'Compara:genome_db_name') and (prop.value in species_list) :
                flies.add(prop.value)
                include = True
                break
        if include :
            assert len(node.sequences) == 1
            msa.append(node.sequences[0].to_seqrecord())


    for a in msa :
        a.id = a.id.replace('Ensembl:', '')

    return tuple(flies), remove_gap_columns(msa)

def filter_phyloxml(tree, species_list) :
    msa = MultipleSeqAlignment([])
    flies = set()
    
    #print >> sys.stderr, "  filtering phyloxml to fly species..."

    for node in tree.get_terminals() :
        include = False

        for prop in node.properties :
            if (prop.ref == 'Compara:genome_db_name') and (prop.value in species_list) :
                flies.add(prop.value)
                include = True
                break

        if include :
            assert len(node.sequences) == 1
            sqrcd = node.sequences[0].to_seqrecord()
            sqrcd.id = node.name
            sqrcd.description = ""
            msa.append(sqrcd)

        if not include :
            tree.prune(node)

    return tuple(flies), remove_gap_columns(msa), tree

def get_genetree(species, symbol) :
    #print >> sys.stderr, "\ngetting genetree (species=%s symbol=%s)" % (species, symbol)

    catcher = 0

    while catcher < 10 :

        try :
            f = urllib2.urlopen(restapi_genetree_url % (species, symbol))
            break

        except (IOError, httplib.HTTPException) :
            print "\n\033[93mrestapi error with %s, %s, trying again in 5 seconds ...\033[0m" % (species, symbol)
            catcher += 1
            time.sleep(5)

    if catcher == 10 :
        sys.exit(1)

    tree = Phylo.read(f, 'phyloxml')
    f.close()

    return filter_phyloxml(tree, fly_species)

def get_protein_to_gene_table(species, symbol) :
    #print >> sys.stderr, "\ngetting homology (species=%s symbol=%s)" % (species, symbol)

    try :
        f = urllib2.urlopen(restapi_homology_url % (species, symbol))

    except urllib2.URLError, ue :
        print >> sys.stderr, "REST api error " + str(ue)
        exit(1)

    data = json.load(f)

    f.close()

    prot2gene = {}

    for d in data["data"][0]["homologies"] :
        if d["target"]["species"] not in fly_species :
            continue
        
        prot2gene[d["source"]["protein_id"]] = d["source"]["id"]
        prot2gene[d["target"]["protein_id"]] = d["target"]["id"]

    return prot2gene

def get_homology_information() :
    species = "tribolium castaneum"
    print >> sys.stderr, "\ngetting %s orthologues from biomart..." % species

    tmp = collections.defaultdict(lambda : collections.defaultdict(set))

    for index,attribute in enumerate(biomart_query_attributes) :
        print >> sys.stderr, "  from %s" % fly_species[index]

        try :
            f = urllib2.urlopen(biomart_url, urllib.urlencode({ 'query' : biomart_query % attribute }))

        except urllib2.URLError, ue :
            print >> sys.stderr, "Error getting %s orthologues from biomart (%s)" % (fly_species[index], str(ue))
            sys.exit(1)

        line_num = 0
        for line in f :
            line_num += 1
            data = line.strip().split(',')
            tc_gene,fly_gene = data

            if not fly_gene :
                continue

            #for index,geneid in enumerate(data[1:]) :
            #    if geneid == '' :
            #        continue
            #
            #    tmp[tc_gene][fly_species[index]].add(geneid)

            tmp[tc_gene][fly_species[index]].add(fly_gene)

            sys.stderr.write("\r    read %d lines for %d %s genes" % (line_num, len(tmp), species))
            sys.stderr.flush()    

        print >> sys.stderr, "\r    read %d lines for %d %s genes" % (line_num, len(tmp), species)

    for tc_gene in tmp :
        for d_gene in tmp[tc_gene] :
            tmp[tc_gene][d_gene] = list(tmp[tc_gene][d_gene])

    return tmp

def save_homology_information_tofile(fname, data) :
    data = json.dumps(data)
    o = open(fname, "w")
    o.write(data)
    #print ("Writing data to file %s...") % fname
    o.close()

def get_homology_information_fromfile(fname) :
    o = open(fname, "r")
    data = json.loads(o.read())
    o.close()
    return data

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

    #print >> sys.stderr, "  contains %d beetle species, %d tribolium casteneum genes" % (len(beetle_species), len(tc_genes))

    return list(beetle_species),tc_genes

def count_seq(fname) :
    count = 0

    for s in SeqIO.parse(fname, 'fasta') :
        header = s.description.split()
        species = header[1].split('=')[1]

        if species == 'tribolium_castaneum' :
            continue

        count += 1

    return count

def get_rid_of_bootstrap(treefile) :
    
    tree = Phylo.read(treefile, 'newick')

    Phylo.write(tree, treefile, 'newick', branch_length_only=True)

def main() :

    if len(sys.argv) != 4 :
        print >> sys.stderr, "Usage: %s <orthologues_input.json> <msa_fname.fasta> <fly output folder>" % sys.argv[0]
        sys.exit(1)

    json_input = sys.argv[1]
    msa_fname = sys.argv[2]
    msa_number = filter(str.isdigit, msa_fname)
    fly_directory = sys.argv[3]
    fly_fasta_path = ("%s/fly%s.fasta") % (fly_directory, msa_number)
    fly_tree_path = ("%s/fly%s.tree") % (fly_directory, msa_number)
    minimum_species = 10 # arbitrary

    global the_name
    the_name = msa_fname

    if not os.path.exists(msa_fname) :
        print >> sys.stderr, "Error: %s does not exist!" % msa_fname
        sys.exit(1)

    orthologues = get_homology_information_fromfile(json_input)
    #orthologues = get_homology_information()

    msa_species,tc_genes = get_msa_species(msa_fname)
    tc_gene = tc_genes[0]

    if tc_gene not in orthologues :
        #print "\nSkipped %s, %s: missing from orthologues" % (msa_fname, tc_gene)
        sys.exit(1)

    """if len(msa_species) < minimum_species :
        #print >> sys.stderr, "\033[93m" + "\nskipping: not enough beetle species..." + "\033[0m"
        print "\nSkipped %s, %s: not enough beetle species" % (msa_fname, tc_gene)
        sys.exit(1)"""

    if len(orthologues[tc_gene]) < minimum_species :
        #print >> sys.stderr, "\033[93m" + "\nskipping: not enough fly species..." + "\033[0m"
        #print "\nSkipped %s, %s: not enough fly species" % (msa_fname, tc_gene)
        sys.exit(1)

    tmp = orthologues[tc_gene]
    tmp_species = tmp.keys()[0]
    tmp_gene = tmp[tmp_species][0]
    tmp_flies, tmp_alignment, tree = get_genetree(tmp_species, tmp_gene)
    
    AlignIO.write(tmp_alignment, fly_fasta_path, 'fasta')

    Phylo.write(tree, fly_tree_path, 'newick')

    get_rid_of_bootstrap(fly_tree_path)

    #print >> sys.stderr, "\n\033[92mfly = %d genes, beetle = %d genes\033[0m" % (len(tmp_alignment), count_seq(msa_fname))
    #print "\n Wrote %s, %s (homologues of %s) \n" % (fly_fasta_path, fly_tree_path, tc_gene)

    return 0

if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by User..."
        sys.exit(1)

