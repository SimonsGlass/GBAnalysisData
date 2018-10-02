
'''
Tristan Sharp
12/1/16
Try to add CNA based on the input positions

To run: ovitos scriptToGetGBs filetoanalysze
/home/tsharp/bin/ovito/ovito-2.8.0-x86_64/bin/ovitos

12/7/2016 I verified that the output from ovitos matches the adaptive CNA from my Windows ovito

2018-07-04 I am going to analyze individual grain boundaries?
'''

import sys
from ovito.io import *
from ovito.modifiers import *
import numpy as np
import matplotlib.pyplot as plt

#if len(sys.argv) != 2:
#        print("Usage: ovitos scriptToGetGBs filetoanalysze")
#        sys.exit()

#infile = sys.argv[1]
nndist = 2.95 # in the length units of the file (Angstrom)
nparticles = 29422

# This is made from "insertGBIDsintodump.cpp"
#infile = "trainingAve.0123000.cfg.phop.CNA0soft.GBID"
infile = "620K/small.0130600.cfg.gaveids.phop.CNA0soft.GBID"
node =  import_file(infile)
#postsim_frame="trainingAve.0129000.cfg.phop.CNA0soft"
postsim_frame="620K/small.0138900.cfg.gaveids.phop.CNA0soft"
dmod = CalculateDisplacementsModifier()
dmod.reference.load(postsim_frame)
node.modifiers.append(dmod)
node.modifiers.append(CommonNeighborAnalysisModifier())

print("Ovito compute")
node.compute()
print("Done Ovito compute")

print(node.output.particle_properties)
dmag = node.output.particle_properties['Displacement Magnitude'].array
disp = node.output.particle_properties['Displacement'].array
print(len(disp[2,:]))

print(dmag[0:30])

# subtract mean displacement
for idir in range(3):
    #meandisp = 0.0
    #for iparticle in range(nparticles):
    #    meandisp += disp[iparticle,idir]
    meandisp = np.sum(disp[:,idir]) / nparticles
    #meandisp /= nparticles
    #for iparticle in range(nparticles):
    disp[:,idir] -= meandisp
dmag = np.sqrt( np.sum( disp, 1))
print(len(dmag))
print(dmag[0:30])

cna = node.output.particle_properties['Structure Type'].array
soft = node.output.particle_properties['soft'].array

GBID_vs_atomindex = node.output.particle_properties['gbid'].array
nGBs = int(np.max(GBID_vs_atomindex))

for isize in [200,600,1200,1800,2400]:
  meansoft = np.zeros(nGBs)
  msd = np.zeros(nGBs)

  for iGB in range(nGBs):
        #if GBsize_vs_GBID[iGB] > 200:
        if np.sum(GBID_vs_atomindex == iGB) > isize:
                print(iGB, "iGB")
                atomids_of_this_gb = np.where(GBID_vs_atomindex == iGB)
                meansoft[iGB] = np.mean(soft[atomids_of_this_gb])
                msd[iGB] = np.sum(dmag[atomids_of_this_gb] ** 2) / len(atomids_of_this_gb)                
  plt.plot(meansoft, msd, 'o',label=">"+str(isize))

xmodel = np.arange(30)/10. - 1.5
ymodel = 2e7*np.exp(1.0 * xmodel - 10.0)
plt.plot(xmodel,ymodel)
#plt.colorbar()
plt.legend()
plt.show()
#plt.savefig()


plt.show()
#plt.savefig()





plt.show()
#plt.savefig()
