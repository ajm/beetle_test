import sys
import os
import glob
from os import system
from sys import stderr

def read_all_files(paml_directory, output) :
    flist = glob.glob("%s/*paml" % paml_directory)

    omegas = []

    for f in flist :
	omega = read_omega(f)
        omegas.append(omega)

    return omegas


def read_omega(fname) :

    with open(fname) as f :
	for line in f :
	    if line.startswith("omega") :
		w = line.split()
		omega = w[3]
    return omega


def write_results(omegas, output) :
    
    o = open(output, "w")
    for omg in omegas :
	o.write(str(omg))
	o.write("\n")
    o.close()

def main() :
        if len(sys.argv) != 3 :
                print >> sys.stderr, "%s <paml_directory> <output>" % sys.argv[0]
                sys.exit(1)

        paml_directory = sys.argv[1]
        output = sys.argv[2]

        if not os.path.exists(paml_directory) :
                print >> sys.stderr, ""

        omegas = read_all_files(paml_directory, output)
        write_results(omegas, output)
        # plotting(output)

        return 0



if __name__ == '__main__' :
        try :
                exit(main())

        except KeyboardInterrupt :
                print >> sys.stderr, "Killed by user..."
                sys.exit(1)

