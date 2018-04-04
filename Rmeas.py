#
# Paul Scherrer Institute
# by Filip Leonarski (filip-karol.leonarski@psi.ch)                            
# with significant contribution by Kay Diederichs (Uni Konstanz) and Meitian Wang (PSI)
#

import sys,numpy

if len(sys.argv) != 5:
    print ""
    print "Usage python Rmeas.py <XDS_ASCII.HKL> <dmax> <dmin> <nshells>"
    print "Calculates X-ray data statistics for arbitrary number of resolution shells."
    print ""
    print "Out of the box only P1,P2, P21, P422, P41212, P43212, I23 and I213 space groups are implemented." 
    print "For all others P1 symmetry is asummed. Modify source code accordingly for other space groups."
    print "Systematic absences are only excluded from data statistics for P41212 and P43212 space groups."
    print "Resolution shells without reflections are not printed."
    print ""
    print "XDS_ASCII.HKL - path to XDS_ASCII.HKL file"
    print "dmax - low resolution limit"
    print "dmin - high resolution limit (if dmin is < 0, than dmin will be automatically extended to a value < |dmin|, if there reflections below the given limit, otherwise |dmin| will be used as cutoff)"
    print "nshells - number of resolution shells"
    print ""
    print "Output is provided on standard output. Run: "
    print "python Rmeas.py <XDS_ASCII.HKL> <dmax> <dmin> <nshells> > out_file"
    print "to save data to a file. First row will provide a description of each reported column."
    print ""
    quit()



intensity = {}
intensity_raw = {}
sigma = {}
Isigma = {}

dmax = float(sys.argv[2])
dmin = float(sys.argv[3])

if dmax <= numpy.abs(dmin):
    print ""
    print "Usage python Rmeas.py <XDS_ASCII.HKL> <dmax> <dmin> <nshells>"
    print ""
    print "dmax needs to be larger than dmin"
    print ""
    quit()

dmin_auto = 0

if dmin < 0:
    dmin = -dmin
    dmin_auto = 1

nshells = int(sys.argv[4])

def res(recip):
    return 1/numpy.sqrt((recip[0]/A)**2+(recip[1]/B)**2+(recip[2]/C)**2)

cnt = 0
cnt1 = 0
cnt2 = 0
sys_cnt = 0

with open(sys.argv[1],'r') as f:
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
	    s = float(linespl[4])
            rlp = float(linespl[8])
            peak = float(linespl[9])
            Iraw = Iobs*(peak/100.0)/rlp

            d = res((h,k,l))
            if dmin_auto == 1 and d < dmin:
                dmin = d
#
# To add a space group use a template below and provide list of equivalent reflections
#
            if space_group == 96 or space_group == 92 or space_group == 89:
                equiv = ((h,k,l),(-h,-k,l),(-k,h,l),(k,-h,l), (-h,k,-l), (h,-k,-l), (k,h,-l),(-k,-h,-l))
            elif space_group == 199 or space_group == 197:
                equiv = ((h,k,l),(-h,-k,l),(-h,k,-l), (h,-k,-l), (l,h,k),(l,-h,-k), (-l,-h,k),(-l,h,-k),(k,l,h),(-k,l,-h),(k,-l,-h),(-k,-l,h))
            elif space_group == 3 or space_group == 4:
                equiv = ((h,k,l),(-h,k,-l))
#           elif space_group == xxx:
#               equiv = ((h,k,l), ... )
            else:
                equiv = ((h,k,l),)


	    if Iobs >= -3*s and s > 0.0 and d >= dmin and d <= dmax:
                if (space_group != 92 and space_group != 96) or not ((h%2 != 0 and k == 0 and l == 0) or (k%2 != 0 and h == 0 and l == 0) or (l%4 != 0 and h == 0 and k == 0)):
	            found = False
	            for hkl in equiv:
		        if hkl in intensity:
                            if not found:
		                intensity[hkl].append(Iobs)
		                intensity_raw[hkl].append(Iraw)
		                sigma[hkl].append((Iobs,s))
		                found = True
	            if not found:
		        intensity[(h,k,l)] = [Iobs]
		        intensity_raw[(h,k,l)] = [Iraw]
		        sigma[(h,k,l)] = [(Iobs,s)]
                else:
                    sys_cnt +=1

for hkl in sigma:
    tmp1 = 0.0
    tmp2 = 0.0
    for x in sigma[hkl]:
        tmp1 += x[0]/(x[1]**2)
        tmp2 += 1/(x[1]**2)
    I = tmp1/tmp2
    s = numpy.sqrt(1/tmp2)
    Isigma[hkl] = (I,s)

shells = []

shell_len = numpy.zeros(nshells)
shell_Iobs = numpy.zeros(nshells)
shell_Iobs2 = numpy.zeros(nshells)
shell_Iraw = numpy.zeros(nshells)
shell_Iraw_min = numpy.zeros(nshells)
shell_Iraw_max = numpy.zeros(nshells)
shell_IS = numpy.zeros(nshells)
shell_IS_1 = numpy.zeros(nshells)
shell_IS_2 = numpy.zeros(nshells)
shell_ISunmrg_1 = numpy.zeros(nshells)
shell_ISunmrg_2 = numpy.zeros(nshells)
shell_ISunmrg_mult_1 = numpy.zeros(nshells)
shell_ISunmrg_3 = numpy.zeros(nshells)
shell_ISunmrg_4 = numpy.zeros(nshells)
shell_Rmerge_1 = numpy.zeros(nshells)
shell_Rmerge_2 = numpy.zeros(nshells)
shell_Rmerge = numpy.zeros(nshells)
shell_unique = numpy.zeros(nshells)
shell_unique_merged = numpy.zeros(nshells)

one_over_dmax2 = 1/(dmax**2)
one_over_dmin2 = 1/(dmin**2)

for i in range(1,nshells+1):
    one_over_d2 = one_over_dmax2 + float(i) * (one_over_dmin2 - one_over_dmax2)/float(nshells)
    shells.append(1/numpy.sqrt(one_over_d2))
    
for i in range(0,nshells):
    shell_Iraw_min[i] = 10e8

for hkl in intensity:
    d = res(hkl)
    sh = 0
    while (sh < nshells-1) and (d < shells[sh]):
       sh += 1

    shell_IS[sh] += Isigma[hkl][0]/Isigma[hkl][1]
    shell_IS_1[sh] += Isigma[hkl][0]
    shell_IS_2[sh] += Isigma[hkl][1]
    shell_unique_merged[sh] += 1.0

    Iobs = numpy.array(intensity[hkl])
    shell_len[sh] += len(Iobs)
    shell_unique[sh] += 1.0
    shell_Iobs[sh] += numpy.sum(Iobs)
    for x in intensity[hkl]:
        shell_Iobs2[sh] += x**2
    shell_Iraw[sh] += numpy.sum(numpy.array(intensity_raw[hkl]))
    Iraw_min = numpy.min(numpy.array(intensity_raw[hkl]))
    if Iraw_min < shell_Iraw_min[sh]:
        shell_Iraw_min[sh] = Iraw_min

    Iraw_max = numpy.max(numpy.array(intensity_raw[hkl]))
    if Iraw_max > shell_Iraw_max[sh]:
        shell_Iraw_max[sh] = Iraw_max

    for x in sigma[hkl]:
        shell_ISunmrg_1[sh] += x[0]/x[1]
        shell_ISunmrg_mult_1[sh] += numpy.sqrt(float(len(Iobs)))*x[0]/x[1]
        shell_ISunmrg_2[sh] += 1.0
        shell_ISunmrg_3[sh] += x[0]
        shell_ISunmrg_4[sh] += x[1]
    avg = numpy.mean(Iobs)
    N = float(len(Iobs))
    if N > 1:
        for x in Iobs:
            shell_Rmerge_1[sh] += numpy.sqrt(N/(N-1))*numpy.abs(avg-x)
            shell_Rmerge_2[sh] += x

for i in range(len(shells)):
    if shell_IS_2[i] > 0:
        shell_IS_1[i] = shell_IS_1[i]/shell_IS_2[i]
    if shell_ISunmrg_2[i] > 0:
        shell_ISunmrg_1[i] = shell_ISunmrg_1[i] / shell_ISunmrg_2[i]
        shell_ISunmrg_3[i] = shell_ISunmrg_3[i] / shell_ISunmrg_4[i]
    if shell_Rmerge_2[i] > 0.0:
        shell_Rmerge[i] = shell_Rmerge_1[i] / shell_Rmerge_2[i]
    else:
        shell_Rmerge[i] = -1.0

print "# 1/d^2      d       N           Nunq     <Iobs>   Rmeas <I/sigma>mrg  <I>/<s>mrg   <Iraw>   Iraw(max)   Iraw(min) <I/sigma>unmrg <I>/<sigma>unmrg "


for i in range(len(shells)):
    if shell_len[i] > 0 and shell_Rmerge[i] >= 0.0 and shell_unique_merged[i] > 0:
        print "%8.5f %6.2f %10d %10d %10.2f %6.2f%% %10.3f %10.3f %11.3f %11.3f %11.3f %11.4f %11.4f"%(1/(shells[i]**2), shells[i], int(shell_len[i]), int(shell_unique[i]), shell_Iobs[i]/shell_len[i], shell_Rmerge[i]*100, shell_IS[i]/ shell_unique_merged[i], shell_IS_1[i], shell_Iraw[i]/shell_len[i], shell_Iraw_max[i], shell_Iraw_min[i], shell_ISunmrg_1[i],shell_ISunmrg_3[i]) 

