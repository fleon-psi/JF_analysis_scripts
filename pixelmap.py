#
# Paul Scherrer Institute
# by Filip Leonarski (filip-karol.leonarski@psi.ch)
# with significant contribution by Kay Diederichs (Uni Konstanz) and Meitian Wang (PSI)
#

import sys
import numpy,math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) != 4:
    print ""
    print "Usage python pixelmap.py <XDS_ASCII.HKL> <dmax> <dmin> <nshells>"
    print ""
    print "Out of the box only P1,P2, P21, P422, P41212, P43212, I23 and I213 space groups are implemented."
    print "For all others P1 symmetry is asummed. Modify source code accordingly for other space groups."
    print ""
    print "XDS_ASCII.HKL - path to XDS_ASCII.HKL file"
    print "dmax - low resolution limit"
    print "dmin - high resolution limit"
    print ""
    print "Resulting pixel map will be saved as Pixelmap.pdf"
    print "To adjust image scale please modify the source code accordingly."
    print ""
    quit()

dmax = float(sys.argv[2])
dmin = float(sys.argv[3])

if dmax <= dmin:
    print ""
    print "Usage python pixelmap.py <XDS_ASCII.HKL> <dmax> <dmin> <nshells>"
    print ""
    print "dmax needs to be larger than dmin"
    print ""
    quit()
	
# Change to adjust scale range for the image
scale = 15

def res(recip,A,B,C):
    return 1/numpy.sqrt((recip[0]/A)**2+(recip[1]/B)**2+(recip[2]/C)**2)

def load_xds(in_filename, Delta_xy):
    ct = 0
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
	        s = numpy.abs(float(linespl[4]))
	        x = float(linespl[5])
	        y = float(linespl[6])
                z = float(linespl[7])
	        xpos = int(x*10+5)%10
	        ypos = int(y*10+5)%10
	        zpos = int(z*10)%10
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

	        if Iobs >= -3*s and s > 0.0 and d >= dmin and d <= dmax:
                    ct += 1
	            found = False
	            for hkl in equiv:
		        if hkl in intensity:
		            intensity[hkl].append(Iobs)
		            intensity_det[hkl].append((Iobs,s,xpos,ypos,x,y,z))
		            found = True
	            if not found:
		        intensity[hkl] = [Iobs]
		        intensity_det[hkl] = [(Iobs,s,xpos,ypos,x,y,z)]


    Delta_xy_cnt = numpy.zeros((10,10))

    for hkl in intensity:
        avg = numpy.mean(intensity[hkl])
        n = len(intensity[hkl])
        if n > 1:            
            for I in intensity_det[hkl]:
                    Delta_xy[I[3],I[2]] += numpy.sqrt(n/(n-1))*((I[0]-avg))*100
                    Delta_xy_cnt[I[3],I[2]] += numpy.abs(avg)

    for i in range(10):
        for j in range(10):
            if Delta_xy_cnt[i,j] == 0.0:
                Delta_xy[i,j] = 0.0
            else:
                Delta_xy[i,j] = Delta_xy[i,j]/Delta_xy_cnt[i,j]

JF1M_Delta_xy = numpy.zeros((10,10))

load_xds(sys.argv[1], JF1M_Delta_xy)

print "Min Delta_xy %f"%numpy.min(JF1M_Delta_xy)
print "Max Delta_xy %f"%numpy.max(JF1M_Delta_xy)

fig = plt.figure(figsize=(4,4))

font = {'family' : 'sans-serif','sans-serif':['Helvetica'],
        'size'   : 10}

mpl.rc('font', **font)
mpl.rc('text', usetex=True)

ax1 = fig.add_subplot(1,1,1)
imgplot = plt.imshow(JF1M_Delta_xy, interpolation='nearest', cmap='bwr',vmin=-scale, vmax=scale)
ax1.get_xaxis().set_ticks([])
ax1.get_yaxis().set_ticks([])
[i.set_linewidth(2) for i in ax1.spines.itervalues()]

mpl.rcParams['axes.linewidth'] = 1

plt.tight_layout()
cbar = plt.colorbar(ax=[ax1], shrink=0.545, ticks=[-scale,-scale/2.0, 0, scale/2.0, scale])
cbar.ax.set_yticklabels(['-%.1f\%%'%scale, '0\%','%.1f\%%'%scale])

if numpy.max(JF1M_Delta_xy) == numpy.min(JF1M_Delta_xy) and numpy.min(JF1M_Delta_xy) == 0.0:
    print "No reflections observed, please check input file"
else:
    plt.savefig("Pixelmap.pdf")
print "Saved to Pixelmap.pdf"
