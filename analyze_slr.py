import sys
import os
import glob

global countingskipped

def read_all_files(in_directory, out_fname) :
    o = open(out_fname, "w")
    filelist = glob.glob("%s/*.slr" % in_directory)
    allomegas = []

    for f in filelist :
        omegas = read_file(f)
        for om in omegas:
            o.write(om)
            o.write("\n")

    o.close()


def read_file(in_fname) :

    # slr file should be like this
    # w[0]=Site, w[1]=Neutral, w[2]=Optimal, w[3]=Omega, w[4]=Lower, w[5]=Upper,
    # w[6]=LRT_Stat, w[7]=Pval, w[8]=Adj.Pval, w[9]=Q-value, w[10]=Result, w[11]=Note

    with open(in_fname) as f :
	# parsing the file, skipping all "Single char" and "All gaps" rows, also too large omegas
        next(f)
	global countingskipped
        omegas = []

        for line in iter(f) :
	    if line.endswith("Single char\n") :
		continue
	    if line.endswith("All gaps\n") :
		continue
            w = line.split()
	    if float(w[3]) <= 5 :
                omegas.append(w[3])
	    else :
		countingskipped += 1

        return omegas


def main() :
    global countingskipped
    global the_name
    countingskipped = 0
    if len(sys.argv) != 3 :
        # slr output files need to have a suffix ".slr" and in the same dir
        print >> sys.stderr, "Usage: %s <slr_directory> <out_fname>" % sys.argv[0]
        sys.exit(1)

    slr_directory = sys.argv[1]
    out_fname = sys.argv[2]

    if not os.path.exists(slr_directory) :
        print >> sys.stderr, "Error: %s does not exist!" % in_fname
        sys.exit(1)

    read_all_files(slr_directory, out_fname)
    print countingskipped, "skipped because they were bigger than 5"
    return 0


if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by User..."
        sys.exit(1)
