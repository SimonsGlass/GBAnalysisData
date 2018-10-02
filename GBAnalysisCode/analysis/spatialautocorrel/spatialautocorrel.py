

# Tristan Sharp 20180402 21:04


doplot = 1 # for running on cluster
import numpy as np
if doplot:
    import matplotlib.pyplot as plt
import time

def smooth(y, box_pts):
            box = np.ones(box_pts)/box_pts
            y_smooth = np.convolve(y, box, mode='same')
            return y_smooth


# Settings.  If change them, use a new version

ver='06' # '00' #  for saving data
regenerate=1

# These defaults may get changed in the if statements below
consideronlytype = -1 # -1 is consider all, 0 is CNA unknown, 1 is FCC, 2-7? are other CNA types, 
# 10 (very disordered) 11 (3.4A near CNA=FCC/HCP) # Super strict filter # Use the aggresive near-FCC filtering
# ver 06 is same as 04 but consider only type 0

# bins of the distance between particles
nbins=6400.0 # 64.0 for 05 # 6400.0 for 04 # 40.0 for 03 #400.0 for 02 #70.0 for 01 #120.0
binwidth=0.025 # 2.86 for 05 # 0.025 for 04 # 2.0 for 03 #0.2 for 02 #0.05 for 01 #0.075#0.025

smoothamount = 12
print([1,2.5,3,4], smooth([1,2.5,3,4],1))

if ver == '00': 
    #filename="150cubesmall43000.cfg.phop.CNA0sof" # used ver none
    filename="50cube_small43000.cfg.phop.CNA0sof" # used ver 00 and 01
    nbins=120.0
    binwidth=0.025 # .075
if ver == '01': 
    filename="50cube_small43000.cfg.phop.CNA0sof" # used ver 00 and 01
    nbins=70.
    binwidth=0.05
if ver == '02':
    filename="50stripsmall43000.cfg.phop.CNA0sof" # used ver 02 03 04
    nbins=400.0
    binwidth=0.2
if ver == '03':
    filename="50stripsmall43000.cfg.phop.CNA0sof" # used ver 02 03 04
    nbins=40.0
    binwidth=2.0
if ver == '04':
    filename="50stripsmall43000.cfg.phop.CNA0sof" # used ver 02 03 04
if ver == '05':
    filename="30stripsmall43000.cfg.phop.CNA0sof" # used ver 04 05. Presented Mach conference 04 (maybe 05)
    nbins=64.0 
    binwidth=2.86 
if ver == '06':
    filename="50rod__small45000.cfg.phop.CNA0sof" # used ver 07 for CNA 0 and ver 08 for CNA 10
    consideronlytype = -1
if ver == '07':
    filename="50rod__small45000.cfg.phop.CNA0sof" # used ver 07 for CNA 0 and ver 08 for CNA 10
    consideronlytype = 0
if ver == '08':
    filename="50rod__small45000.cfg.phop.CNA0sof" # used ver 07 for CNA 0 and ver 08 for CNA 10
    consideronlytype = 10
if ver == '09':
    filename="150cubesmall45000.cfg.phop.CNA0sof" 
    consideronlytype = -1
if ver == '10':
    filename="150cubesmall45000.cfg.phop.CNA0sof" 
    consideronlytype = 0
if ver == '11':
    filename="150cubesmall45000.cfg.phop.CNA0sof"
    consideronlytype = 10




# Dont change below here
minbincenter=2.0
bincenters=minbincenter+np.arange(nbins)*binwidth
binedges=minbincenter+np.arange(nbins+1)*binwidth-binwidth/2.0
minbinedge=binedges[0]
maxbinedge=binedges[nbins]
maxbinedgesq=maxbinedge**2.0

#L=50.0 # system length, equal in x,y,z directions
#center=150.0 # coordinate of system middle, equal in x,y,z directions
nndist = 2.86 # Angstroms. From FCC of E minimized state

filearray = np.genfromtxt(filename, skip_header=9)
myid = filearray[:,0]
x = filearray[:,1]
y = filearray[:,2]
z = filearray[:,3]
phop = filearray[:,4]
soft = filearray[:,5]

CNA = np.array([])
tstep = '1'+filename[12:17]
print("ts", tstep)
if consideronlytype >= 0 and consideronlytype < 10:
    CNA = np.genfromtxt("CNAtrainingAve.0"+tstep+".cfg", skip_header=9)
if consideronlytype >= 10:
    CNA = np.genfromtxt("NearCNAtrainingAve.0"+tstep+".cfg", skip_header=9)
numatomsincluded = 0

numinbin = np.zeros(nbins)
phop0phoprSUM = np.zeros(nbins)
soft0softrSUM = np.zeros(nbins)

phopSUM=(0.0)
phopsquaredSUM=(0.0)
softSUM=(0.0)
softsquaredSUM=(0.0)

if regenerate == 1:
  #print(bincenters)
  #print(binedges)
  dolength = len(filearray)*1+0*1000 # reduce size for fast testing
  for iatom in range(0,dolength):

    if consideronlytype < 0 or (consideronlytype == CNA[-1 + myid[iatom]]):
      numatomsincluded += 1

      xi=(x[iatom])
      yi=(y[iatom])
      zi=(z[iatom])

      phopi=(phop[iatom])
      phopSUM += phopi
      phopsquaredSUM += phopi**2.0

      softi=(soft[iatom])
      softSUM += softi
      softsquaredSUM += softi**2.0
      
      for jatom in range(iatom+1,dolength):
              
          if consideronlytype < 0 or (consideronlytype == CNA[-1 + myid[jatom]]):

                  #if (x[jatom] - xi)**2.0 < maxbinedgesq: # with the strip, x and y are always nearby
                  #if (y[jatom] - yi)**2.0 < maxbinedgesq:
                  if (z[jatom] - zi)**2.0 < maxbinedgesq:
                      distancesq = (x[jatom] - xi)**2.0 + (y[jatom] - yi)**2.0 + (z[jatom] - zi)**2.0
                      if (distancesq < maxbinedgesq):
                          #distancesq = (x[jatom]-xi)**2.0 + (y[jatom]-yi)**2.0 + (z[jatom] - zi)**2.0
                          binindex = int((np.sqrt(distancesq) - minbinedge)/binwidth)
                          numinbin[binindex] += 1.0
                          phop0phoprSUM[binindex] += phop[iatom]*phop[jatom]
                          soft0softrSUM[binindex] += soft[iatom]*soft[jatom]

  meanphop = phopSUM/numatomsincluded
  mean_of_phopsq = phopsquaredSUM/numatomsincluded
  meanphop0phopr = phop0phoprSUM/numinbin

  meansoft = softSUM/numatomsincluded
  mean_of_softsq = softsquaredSUM/numatomsincluded
  meansoft0softr = soft0softrSUM/numinbin

  outputfile = open('outfiles/bincenters_v'+ver+'.txt','w')
  for item in bincenters:
          outputfile.write("%s\n" % str(item))
  outputfile = open('outfiles/numinbin_v'+ver+'.txt','w')
  for item in numinbin:
          outputfile.write("%s\n" % str(item))
  outputfile = open('outfiles/meanphop_v'+ver+'.txt','w')
  for item in [meanphop]:
          outputfile.write("%s\n" % str(item))
  outputfile = open('outfiles/mean_of_phopsq_v'+ver+'.txt','w')
  for item in [mean_of_phopsq]:
          outputfile.write("%s\n" % str(item))
  outputfile = open('outfiles/meanphop0phopr_v'+ver+'.txt','w')
  for item in meanphop0phopr:
          outputfile.write("%s\n" % str(item))
  outputfile = open('outfiles/meansoft_v'+ver+'.txt','w')
  for item in [meansoft]:
          outputfile.write("%s\n" % str(item))
  outputfile = open('outfiles/mean_of_softsq_v'+ver+'.txt','w')
  for item in [mean_of_softsq]:
          outputfile.write("%s\n" % str(item))
  outputfile = open('outfiles/meansoft0softr_v'+ver+'.txt','w')
  for item in meansoft0softr:
          outputfile.write("%s\n" % str(item))

else: # dont regenerate
  bincenters = np.genfromtxt('outfiles/bincenters_v'+ver+'.txt')
  numinbin = np.genfromtxt('outfiles/numinbin_v'+ver+'.txt')
  meanphop = np.genfromtxt('outfiles/meanphop_v'+ver+'.txt')
  mean_of_phopsq = np.genfromtxt('outfiles/mean_of_phopsq_v'+ver+'.txt')
  meanphop0phopr = np.genfromtxt('outfiles/meanphop0phopr_v'+ver+'.txt')
  meansoft = np.genfromtxt('outfiles/meansoft_v'+ver+'.txt')
  mean_of_softsq = np.genfromtxt('outfiles/mean_of_softsq_v'+ver+'.txt')
  meansoft0softr = np.genfromtxt('outfiles/meansoft0softr_v'+ver+'.txt')


if doplot:
    volumecontributingtobin =  4.0 * 3.141592 * bincenters**2.0 * binwidth
    numberofatomcenters = dolength / 2.0
    plt.plot(bincenters/nndist, numinbin / volumecontributingtobin * nndist**3.0 / numberofatomcenters)
    plt.plot(bincenters/nndist,np.zeros(len(bincenters)),marker='',lw=1, color='black')
    plt.ylabel('g(r)')
    plt.xlabel('r/d')
    plt.grid(b=True, which='major')
    plt.axis([0,4,0,10000])
    #plt.axis([0,50,0,3500])
    plt.savefig('outfiles/gr_v'+ver+'.pdf')
    plt.show()
    #plt.plot(bincenters/nndist, meanphop0phopr)
    #plt.ylabel('<phop(0) phop(r)>')
    
    normalizedcorrelphop = (meanphop0phopr - meanphop**2.0)/(mean_of_phopsq - meanphop**2.0)
    plt.plot(bincenters/nndist, np.log(smooth(normalizedcorrelphop,smoothamount)),marker='')
    #plt.plot(bincenters/nndist,np.zeros(len(bincenters)),marker='',lw=1, color='black')
    #plt.plot(bincenters/nndist, np.log10(normalizedcorrelphop))
    plt.ylabel('Log of Normalized Correl (Sci 2017 Eq. 2 p_hop for D^2_min)')
    plt.ylabel('Normalized Correl (Sci 2017 Eq. 2 p_hop for D^2_min)')
    plt.xlabel('r/d')
    plt.grid(b=True, which='major')
    #plt.axis([0,10,-0.25,5.0])
    #plt.axis([0.8,4,-0.25,5.0])
    plt.savefig('outfiles/zoomsmoothedAutophop_v'+ver+'.pdf')
    plt.show()
    
    normalizedcorrelsoft = (meansoft0softr - meansoft**2.0)/(mean_of_softsq - meansoft**2.0)
    #plt.plot(bincenters/nndist, smooth(normalizedcorrelsoft,smoothamount),marker='')
    #plt.plot(bincenters/nndist,np.zeros(len(bincenters)),marker='',lw=1, color='black')
    plt.plot(bincenters/nndist, np.log10(normalizedcorrelsoft))
    
    plt.ylabel('Log of Normalized Correl (Sci 2017 Eq. 3)')
    plt.ylabel('Normalized Correl (Sci 2017 Eq. 3)')
    plt.xlabel('r/d')
    plt.grid(b=True, which='major')
    #plt.axis([0,10,-2.5,1.1])
    #plt.axis([0.8,4,-2.5,1.1])
    plt.savefig('outfiles/zoomsmoothedAutoS_v'+ver+'.pdf')
    plt.show()



