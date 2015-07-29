import sys
import os
import glob
from os import system
from sys import stderr
from Bio import AlignIO

READGAPS = False
all_omegas = []
all_means = []
filtered_all_means = []

def read_all_files(slr_directory, output) :
	flist = glob.glob("%s/*slr" % slr_directory)

	for f in flist :
		omegas = read_omegas(f)
		gaps = read_gaps(f)
		parse_omegas_with_gaps(omegas, gaps, output)

def read_omegas(fname) :
	global all_omegas
	global all_means
	global filtered_all_means
	file_omegas = []
	filtered_file_omegas = []

	with open(fname) as f :
		next(f)
		for line in iter(f) :

			# basic
			if line.endswith("Single char\n") :
				continue
			if line.endswith("All gaps\n") :
				continue
			w = line.split()
			omg = float(w[3])

			# all omegas in file omegas for counting the unfiltered mean
			file_omegas.append(omg)

			# filtering
			if omg <= 2 :
				all_omegas.append(omg)
				filtered_file_omegas.append(omg)

	# get the means of the file and add them to the total
	all_means.append(sum(file_omegas)/len(file_omegas))
	filtered_all_means.append(sum(filtered_file_omegas)/len(filtered_file_omegas))

	return file_omegas

def read_gaps(fname) :
	gaps = []
	if READGAPS :
		fname = fname.replace("slr", "fasta")
		msa = AlignIO.read(open(fname), "fasta")
	

		for i in range (0, msa.get_alignment_length(), 3) :
			c1 = msa[:,i]
			c2 = msa[:,i+1]
			c3 = msa[:,i+2]

			countgaps = 0
			for j in range(0, len(c1)) :
				if (c1[j] == '-') or (c2[j] == '-') or (c3[j] == '-') :
					countgaps += 1
			if countgaps < 11 :
				gaps.append(countgaps)

	return gaps

def write_files(output) :
	global all_omegas
	global all_means
	global filtered_all_means

	omegas_output = output + ".results"
	means_output = output + "_means.results"
	fltrd_means_output = output + "_fltrd_means.results"

	o = open(omegas_output, "w")
	for omg in all_omegas :
		o.write(str(omg))
		o.write("\n")
	o.close()

	o = open(means_output, "w")
	for mean in all_means :
		o.write(str(mean))
		o.write("\n")
	o.close()

	o = open(fltrd_means_output, "w")
	for fltrd_mean in filtered_all_means :
		o.write(str(fltrd_mean))
		o.write("\n")
	o.close()

def parse_omegas_with_gaps(omegas, gaps, output) :
	if READGAPS :
		print "parsing gaps not yet implemented"

def plotting(output) :
    scriptfolder = os.path.dirname(os.path.realpath(__file__))
    command = "Rscript %s/plot.R %s" % (scriptfolder, output)
    #command +=  " > /dev/null 2> /dev/null"

    print "\tplotting..."

    if system(command) != 0 :
        print >> stderr, "FAIL!\n" + " ".join(command.split())
        exit(1)

def main() :
	if len(sys.argv) != 3 :
		print >> sys.stderr, "%s <slr directory> <output>" % sys.argv[0]
		sys.exit(1)

	slr_directory = sys.argv[1]
	output = sys.argv[2]

	if not os.path.exists(slr_directory) :
		print >> sys.stderr, ""
	
	read_all_files(slr_directory, output)
	write_files(output)
	plotting(output)

	return 0

if __name__ == '__main__' :
	try :
		exit(main())

	except KeyboardInterrupt :
		print >> sys.stderr, "Killed by user..."
		sys.exit(1)
