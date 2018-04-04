#
# Paul Scherrer Institute
# by Filip Leonarski (filip-karol.leonarski@psi.ch)                            
# with significant contribution by Kay Diederichs (Uni Konstanz) and Meitian Wang (PSI)
#

import sys

if len(sys.argv) != 3:
    print ""
    print "Usage python CorrDataset.py <INTEGRATE.HKL 1> <INTEGRATE.HKL 2> "
    print "Extracts intensities and highest pixel count (MAXC) of common reflections (no symmetry applied) between 2 datasets."
    print "INTEGRATE.HKL created by XDS are necessary."
    print ""
    print "INTEGRATE.HKL 1 - path to INTEGRATE.HKL file of the first dataset generated by XDS"
    print "INTEGRATE.HKL 2 - path to INTEGRATE.HKL file of the second dataset generated by XDS"
    print ""
    quit()

import numpy, sys

intensity1 = {}
intensity2 = {}
resolution1 = {}
resolution2 = {}

def analyze(intensity, fname, resolution):
    with open(fname,'r') as f:
        for line in f:
            if len(line) > 15 and line[:13] == '!UNIT_CELL_CO':
                splitl = line.split()
                A = float(splitl[1])
                B = float(splitl[2])
                C = float(splitl[3])
                alpha = float(splitl[4])
                beta = float(splitl[5])
                gamma = float(splitl[6])
	    if line[0] != '!': 
	        linespl = line.split()
	        h = int(linespl[0])
	        k = int(linespl[1])
	        l = int(linespl[2])
	        recip = (h,k,l)
	        Iobs = float(linespl[3])
                peak = float(linespl[9])

                if "INTEGRATE.HKL" in fname:
	            maxc = float(linespl[11])
                else:
                    maxc = 0.0
                if peak == 100.0:
		    intensity[recip] = (Iobs,maxc)
                    resolution[recip] = 1/numpy.sqrt((h/A)**2+(k/B)**2+(l/C)**2)

analyze(intensity1, sys.argv[1],resolution1)
analyze(intensity2, sys.argv[2],resolution2)

print "#h k l Iobs(1) MAXC(1) Iobs(2) MAXC(2) d"

for recip in intensity1:
    if recip in intensity2:
	print "%d %d %d %f %f %f %f %f"%(recip[0],recip[1],recip[2],intensity1[recip][0],intensity1[recip][1],intensity2[recip][0],intensity2[recip][1],resolution1[recip])


