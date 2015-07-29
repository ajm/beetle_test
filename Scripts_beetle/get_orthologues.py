import sys
import urllib2
import urllib
import xml.sax.saxutils
import os
import tempfile
import glob
import collections
import json
import re

from Bio import Phylo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment


biomart_url = "http://metazoa.ensembl.org/biomart/martservice"

biomart_query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "metazoa_mart_27" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.7" >
            
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
    print ("Writing data to file %s...") % fname
    o.close()


def get_homology_information_fromfile(fname) :
    o = open(fname, "r")
    data = json.loads(o.read())[0]
    o.close()
    return data


def main() :

    if len(sys.argv) != 2 :
        print >> sys.stderr, "Usage: %s <output_fname>" % sys.argv[0]
        sys.exit(1)

    output_fname = sys.argv[1]
    orthologues = get_homology_information()
    #orthologues = get_homology_information_fromfile(input_fname)
    save_homology_information_tofile(output_fname, orthologues)

    return 0


if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by User..."
        sys.exit(1)
