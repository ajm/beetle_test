import sys
import urllib2
import urllib
import xml.sax.saxutils
import os
import tempfile
import glob
import collections
import json

from Bio import Phylo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment



biomart_url = "http://metazoa.ensembl.org/biomart/martservice"

biomart_query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "metazoa_mart_26" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.7" >
            
    <Dataset name = "tcastaneum_eg_gene" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "dananassae_eg_gene" />
        <Attribute name = "derecta_eg_gene" />
        <Attribute name = "dgrimshawi_eg_gene" />
        <Attribute name = "dmelanogaster_eg_gene" />
        <Attribute name = "dmojavensis_eg_gene" />
        <Attribute name = "dpersimilis_eg_gene" />
        <Attribute name = "dpseudoobscura_eg_gene" />
        <Attribute name = "dsechellia_eg_gene" />
        <Attribute name = "dsimulans_eg_gene" />
        <Attribute name = "dvirilis_eg_gene" />
        <Attribute name = "dwillistoni_eg_gene" />
        <Attribute name = "dyakuba_eg_gene" />
    </Dataset>
</Query>
"""

biomart_query2 = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "metazoa_mart_26" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.7" >
            
    <Dataset name = "tcastaneum_eg_gene" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "%s" />
    </Dataset>
</Query>
"""

biomart_query2_attributes = (
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
    tmp = None
    last = 0

    for i in range(msa.get_alignment_length()) :
        col = msa[:,i]

        if col.count('-') == len(col) :
            continue

        if tmp :
            tmp += msa[:, i:i+1]
        else :
            tmp = msa[:, i:i+1]

    print >> sys.stderr, "  remove gap columns %dx%d --> %dx%d" % (
                len(msa), msa.get_alignment_length(), 
                len(tmp), tmp.get_alignment_length())

    return tmp

def filter_alignment_from_phyloxml(tree, species_list) :
    msa = MultipleSeqAlignment([])
    flies = set()

    print >> sys.stderr, "  filtering alignment to fly species..."

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

    return tuple(flies), remove_gap_columns(msa)

def get_genetree(species, symbol) :
    print >> sys.stderr, "\ngetting genetree (species=%s symbol=%s)" % (species, symbol)

    try :
        f = urllib2.urlopen(restapi_genetree_url % (species, symbol))

    except urllib2.URLError, ue :
        print >> sys.stderr, "REST api error " + str(ue)
        exit(1)

    tree = Phylo.read(f, 'phyloxml')

    f.close()

    return filter_alignment_from_phyloxml(tree, fly_species)

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

    for index,attribute in enumerate(biomart_query2_attributes) :
        print >> sys.stderr, "  from %s" % fly_species[index]

        try :
            f = urllib2.urlopen(biomart_url, urllib.urlencode({ 'query' : biomart_query2 % attribute }))

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

def get_msa_species(fname) :
    beetle_species = set()
    tc_genes = []

    print >> sys.stderr, "\nreading %s ... " % fname

    for s in SeqIO.parse(fname, 'fasta') :
        header = s.description.split()
        species = header[1].split('=')[1]
        geneid = header[2].split('=')[1]

        if species == 'tribolium_castaneum' :
            tc_genes.append(geneid)
            continue

        beetle_species.add(species)

    print >> sys.stderr, "  contains %d beetle species, %d tribolium casteneum genes" % (len(beetle_species), len(tc_genes))

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

def main() :

    if len(sys.argv) != 2 :
        print >> sys.stderr, "Usage: %s <scaffold directory>" % sys.argv[0]
        sys.exit(1)

    scaffold_directory = sys.argv[1]
    minimum_species = 10 # arbitrary

    if not os.path.exists(scaffold_directory) :
        print >> sys.stderr, "Error: %s does not exist!" % scaffold_directory
        sys.exit(1)

    # step 1.   the ensembl rest api does not work for tribolium castaneum (tc genes do not appear to 
    #           have ensembl ids), so we need to use biomart to get homology information for all tribolium
    #           genes vs all fly species 
    orthologues = get_homology_information()

    # step 2.   iterate over all multiple sequence alignments output by glutton, these all have filenames
    #           like msaXXX.fasta (where XXX is an integer) 
    for msa_fname in glob.glob(os.path.join(scaffold_directory, "msa*.fasta")) :
        
        # step 3.  read the MSA file and get two pieces of information
        #               1. the names of beetle species apart from tribolium castaneum (msa_species)
        #               2. the names of the tribolium genes in the alignment (tc_genes)
        msa_species,tc_genes = get_msa_species(msa_fname)
        tc_gene = tc_genes[0]

        # step 4a.  reject multiple sequence alignments with less than a minimum number of beetle species
        if len(msa_species) < minimum_species :
            print >> sys.stderr, "\033[93m" + "\nskipping: not enough beetle species..." + "\033[0m"
            continue

        # step 4b.  reject gene trees with less than a minimum number of fly species
        if len(orthologues[tc_gene]) < minimum_species :
            print >> sys.stderr, "\033[93m" + "\nskipping: not enough fly species..." + "\033[0m"
            continue

        # step 5.   we need to get the alignment from ensembl genetrees for the fly species,
        #           we only need a single gene from a single fly species 
        tmp = orthologues[tc_gene]
        tmp_species = tmp.keys()[0]
        tmp_gene = tmp[tmp_species][0]

        tmp_flies, tmp_alignment = get_genetree(tmp_species, tmp_gene)

        # step 6.   now we have two alignments, one for beetles (msa_fname) and one for flies (tmp_alignment)
        #           so we can do some analysis (but for now just print out some metadata)
        #
        #AlignIO.write(tmp_alignment, 'filename.fasta', 'fasta')
        print >> sys.stderr, "\n\033[92mfly = %d genes, beetle = %d genes\033[0m" % (len(tmp_alignment), count_seq(msa_fname))

    return 0

if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by User..."
        sys.exit(1)

