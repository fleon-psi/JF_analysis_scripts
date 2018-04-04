#
# Paul Scherrer Institute
# by Filip Leonarski (filip-karol.leonarski@psi.ch)                            
# with significant contribution by Kay Diederichs (Uni Konstanz) and Meitian Wang (PSI)
#

import sys
import numpy,math

if len(sys.argv) != 5:
    print ""
    print "Usage python Rpixel.py <XDS_ASCII.HKL> <dmax> <dmin> <nshells>"
    print ""
    print "Out of the box only P1,P2, P21, P422, P41212, P43212, I23 and I213 space groups are implemented." 
    print "For all others P1 symmetry is asummed. Modify source code accordingly for other space groups."
    print ""
    print "XDS_ASCII.HKL - path to XDS_ASCII.HKL file"
    print "dmax - low resolution limit"
    print "dmin - high resolution limit"
    print "nshells - number of resolution shells"
    print ""
    print "Output is provided on standard output. Run: "
    print "python Rpixel.py <XDS_ASCII.HKL> <dmax> <dmin> <nshells> > out_file"
    print "to save data to a file. First row will provide a description of each reported column."
    print ""
    quit()

def res(recip,A,B,C):
    return 1/numpy.sqrt((recip[0]/A)**2+(recip[1]/B)**2+(recip[2]/C)**2)

def load_xds(in_filename, d1, d2):
    intensity = {}
    intensity_det = {}
    with open(in_filename,'r') as f:
        for line in f:
	    linespl = line.split()
	    if '!SPACE_GROUP_NUMBER=' == linespl[0]:
	        space_group = int(linespl[1])

	    elif '!UNIT_CELL_CONSTANTS=' == linespl[0]:
	        A = float(linespl[1])
	        B = float(linespl[2])
	        C = float(linespl[3])
	        a = float(linespl[4])
	        b = float(linespl[5])
	        g = float(linespl[6])

	    elif line[0] != '!': 
	        h = int(linespl[0])
	        k = int(linespl[1])
	        l = int(linespl[2])
	        Iobs = float(linespl[3])
	        s = numpy.abs(float(linespl[4])) # remove numpy.abs to exlcude outliers
	        x = float(linespl[5])
	        y = float(linespl[6])
	        z = float(linespl[7])
	        xpos = int(x*10+5)%10 # In XDS center of reflection is (0,0)
	        ypos = int(y*10+5)%10 # In XDS center of reflection is (0,0)
	        zpos = int(y*10)%10
                d = res((h,k,l),A,B,C)

# To add a space group use a template below and provide list of equivalent reflections
                if space_group == 96 or space_group == 92 or space_group == 89:
                    equiv = ((h,k,l),(-h,-k,l),(-k,h,l),(k,-h,l), (-h,k,-l), (h,-k,-l), (k,h,-l),(-k,-h,-l))
                elif space_group == 199 or space_group == 197:
                    equiv = ((h,k,l),(-h,-k,l),(-h,k,-l), (h,-k,-l), (l,h,k),(l,-h,-k), (-l,-h,k),(-l,h,-k),(k,l,h),(-k,l,-h),(k,-l,-h),(-k,-l,h))
                elif space_group == 3 or space_group == 4:
                    equiv = ((h,k,l),(-h,k,-l))
#               elif space_group == xxx:
#                   equiv = ((h,k,l), ... )
                else:
                    equiv = ((h,k,l),)

	        if s > 0.0 and Iobs/s > -3*s and d1 <= d < d2:
	            found = False
	            for hkl in equiv:
		        if hkl in intensity:
		            intensity[hkl].append(Iobs)
		            intensity_det[hkl].append((Iobs,s,xpos,ypos,zpos,y))
		            found = True
	            if not found:
		        intensity[hkl] = [Iobs]
		        intensity_det[hkl] = [(Iobs,s,xpos,ypos,zpos,y)]

    Delta_xy = numpy.zeros((10,10))
    R_meas = 0.0

    Iobs_sum = 0.0

    for hkl in intensity:
        avg = numpy.mean(intensity[hkl])
        n = float(len(intensity[hkl]))
        if n > 1:
            for I in intensity_det[hkl]:
                    Delta_xy[I[3],I[2]] += numpy.sqrt(n/(n-1))*((I[0]-avg))*100
                    R_meas += numpy.sqrt(n/(n-1))*numpy.abs((I[0]-avg))*100
                    Iobs_sum += avg
    
    R_pixel = 0.0
    R_meas = R_meas / Iobs_sum
    for i in range(10):
         for j in range(10):
             R_pixel += numpy.abs(Delta_xy[i,j])
    R_pixel = R_pixel / Iobs_sum
    print("%f %8.3f %8.3f"%(1/(d1**2),R_meas,R_pixel))


nshells = int(sys.argv[4])
dmin = float(sys.argv[3])
dmax = float(sys.argv[2])

print "#1/d^2      R_meas   R_pixel"

shells = [dmax]

one_over_dmin2 = 1/(dmin**2)
one_over_dmax2 = 1/(dmax**2)

for i in range(1,nshells+1):
    one_over_d2 = one_over_dmax2 + float(i) * (one_over_dmin2 - one_over_dmax2)/float(nshells)
    shells.append(1/numpy.sqrt(one_over_d2))

for i in range(0,nshells):
    load_xds(sys.argv[1],shells[i+1],shells[i])
