
'''
Tristan Sharp
To run: ovitos scriptToGetGBs filetoanalysze
/home/tsharp/bin/ovito/ovito-2.8.0-x86_64/bin/ovitos

If do not need to regenerate MSD and mean softness, but just plot then ovitos not needed

12/1/16
Try to add CNA based on the input positions

12/7/2016 
I verified that the output from ovitos matches the adaptive CNA from my Windows ovito


2018-07-04 I am going to analyze individual grain boundaries?
'''

import sys
import numpy as np
import matplotlib.pyplot as plt

#if len(sys.argv) != 2:
#        print("Usage: ovitos scriptToGetGBs filetoanalysze")
#        sys.exit()

#infile = sys.argv[1]
nndist = 2.95 # in the length units of the file (Angstrom)
regenerate = 0
averageallGBs = 0 # 1 Hijacking with 1 to not look at individual GBs

fig1 = plt.figure(1,figsize=(5,5))
axs = fig1.add_subplot(111)
hotclr=['purple','b','teal','green','gold','r','k']
hotclr = hotclr+hotclr+hotclr

#case = sys.argv[1] # 5
for case in [5]: # [1,2,3,4,5,7]:#[7,5,4,3,2,1]:#[5,7]:#[1, 2, 3, 4, 5, 7]:
  if case == 0:
      # This is made from "insertGBIDsintodump.cpp"
      infile = "620K/small.0130600.cfg.gaveids.phop.CNA0soft.GBID"
      postsim_frame="620K/small.0138900.cfg.gaveids.phop.CNA0soft"
  elif case == 1:
      infile = "310K/trainingAve.0143000.cfg.phop.CNA0soft.GBID"
      postsim_frame="310K/trainingAve.0149000.cfg.phop.CNA0soft"
  elif case == 2:
      infile = "463K/trainingAve.0143000.cfg.phop.CNA0soft.GBID"
      postsim_frame="463K/trainingAve.0149000.cfg.phop.CNA0soft"
  elif case == 3:
      infile = "556K/trainingAve.0133000.cfg.phop.CNA0soft.GBID"
      postsim_frame="556K/trainingAve.0139000.cfg.phop.CNA0soft"
  elif case == 4:
      infile = "620K/trainingAve.0133000.cfg.phop.CNA0soft.GBID"
      postsim_frame="620K/trainingAve.0139000.cfg.phop.CNA0soft"
  elif case == 5:
      infile = "694K/trainingAve.0123000.cfg.phop.CNA0soft.GBID"
      postsim_frame="694K/trainingAve.0129000.cfg.phop.CNA0soft"
  elif case == 6:
      infile = "756K/trainingAve.0123000.cfg.phop.CNA0soft.GBID"
      postsim_frame="756K/trainingAve.0129000.cfg.phop.CNA0soft"
  elif case == 7:
      infile = "834K/trainingAve.0143000.cfg.phop.CNA0soft.GBID"
      postsim_frame="834K/trainingAve.0149000.cfg.phop.CNA0soft"
  
  if regenerate == 1:
    from ovito.io import *
    from ovito.modifiers import *
    node =  import_file(infile)
    
    dmod = CalculateDisplacementsModifier()
    dmod.reference.load(postsim_frame)
    node.modifiers.append(dmod)
    node.modifiers.append(CommonNeighborAnalysisModifier())
    
    print("Ovito compute")
    node.compute()
    print("Done Ovito compute")
    print(node.output.particle_properties)
    dmag = np.array(node.output.particle_properties['Displacement Magnitude'].array)
    disp = np.array(node.output.particle_properties['Displacement'].array)
    pos = node.output.particle_properties['Position'].array
    cna = node.output.particle_properties['Structure Type'].array
    soft = node.output.particle_properties['soft'].array
    GBID_vs_atomindex = node.output.particle_properties['gbid'].array
    
    import pickle
    f = open(infile+"temp_picklevars", 'wb')
    pickle.dump((dmag,disp,pos,cna,soft,GBID_vs_atomindex), f)
    f.close()
    
  else:
    import pickle
    f = open(infile+"temp_picklevars", 'rb')
    dmag,disp,pos,cna,soft,GBID_vs_atomindex = pickle.load(f)
    f.close()

  print(len(cna))
  print(cna[0])

  nGBs = int(np.max(GBID_vs_atomindex))
  print(nGBs, " is number of GBs of all sizes")
  
  if averageallGBs:
          nGBs=1

  for minsize in [40]:#[40,160,640,2560]:
      meansoft = np.zeros(nGBs)-1.0
      msd = np.zeros(nGBs)-1.0
      
      for iGB in range(nGBs):
          #if GBsize_vs_GBID[iGB] > 200:
          if np.sum(GBID_vs_atomindex == iGB) > minsize:
                  #print(iGB, "iGB has ", np.sum(GBID_vs_atomindex == iGB)," atoms")
                  gbcriterion = (GBID_vs_atomindex == iGB)
                  if averageallGBs:
                          gbcriterion = (GBID_vs_atomindex >= 1)
                  #gbcriterion = np.all([GBID_vs_atomindex == iGB, cna != 1],0)
                  atomids_of_this_gb = (np.where(gbcriterion))[0]
                  meansoft[iGB] = np.mean(soft[atomids_of_this_gb])
                  # ALTERNATIVE
                  #S_of_crystal = -2.25
                  # mygbsize = np.sum(soft[atomids_of_this_gb] - S_of_crystal)
                  # meansoft[iGB] = np.sum((soft[atomids_of_this_gb] - S_of_crystal)**2.0)/mygbsize
                  
                  # calculate MSD of the GB. First subtract mean motion of this GB
                  for idir in range(3):
                      meandisp = np.sum(disp[atomids_of_this_gb,idir]) / len(atomids_of_this_gb)
                      disp[atomids_of_this_gb,idir] -= meandisp
                      #print("Mean displacement was ", meandisp)
                  dmag_of_gbatoms = np.sqrt( np.sum( (disp[atomids_of_this_gb,:])** 2.0, 1))
                  msd[iGB] = np.sum(dmag_of_gbatoms ** 2.0) / len(atomids_of_this_gb)                
                  # ALTERNATIVE
                  # msd[iGB] = np.sum(dmag_of_gbatoms ** 2.0) / mygbsize                
      
      here=np.where(msd != -1.0)
      print("Num data points ", len((here)[0]))
      #label= ">"+str(minsize)+"_"+infile[0:4]+infile[17:23]
      label= infile[0:4]
      x = meansoft[here]
      y = msd[here] # np.log10(msd[here])
      simtime=6.0 #psec
      plt.plot(x, np.log10(y/simtime), 'o', color=hotclr[case], ms=4.0)
      plt.savefig("tempdel.pdf")
      plt.show()
    
      Sbinmidstart = -4
      Sbinsize = 0.25
      NSbins = 40
      binleftedges = np.arange(NSbins)*Sbinsize + Sbinmidstart - Sbinsize/2.0
      sumyvaluesinbin = np.zeros(NSbins)
      sumy2valuesinbin = np.zeros(NSbins)
      numinbin = np.zeros(NSbins)
      for idatapoint in range(len(x)):
        if np.abs(x[idatapoint]) > NSbins *Sbinsize + Sbinmidstart:
            print(x[idatapoint])
        foundbin = ((np.where(binleftedges > x[idatapoint]))[0])[0]
        sumyvaluesinbin[foundbin-1] += y[idatapoint]
        sumy2valuesinbin[foundbin-1] += y[idatapoint]*y[idatapoint]
        numinbin[foundbin-1] += 1.0
      
      meanyinbin = sumyvaluesinbin/numinbin
      meany2inbin = sumy2valuesinbin/numinbin
      #axs.errorbar(binleftedges+Sbinsize/2.0, 
      #             meanyinbin,
      #             #yerr = np.sqrt(meany2inbin - meanyinbin*meanyinbin),
      #             label=label, color=hotclr[case], marker = 'o')



      coefs = np.polyfit(x, y, 1)
      print("yoffset ", coefs[1])
      print("slope = ", coefs[0])
      plt.plot(x, 1.0*x*coefs[0] + coefs[1], color=hotclr[case])

      xmodel = np.arange(45)/10. - 1.
      msdmodel = (0.4*np.exp(1.0 * xmodel))
      #plt.plot(xmodel, np.log10(1.0/msdmodel), color='k')
      plt.plot(xmodel, msdmodel, color='k')
  
#plt.colorbar()

# From tags 50-56
Sall=np.array([-2.5,-2.25,-2.,-1.75,-1.5,-1.25,-1.,-0.75,-0.5,-0.25,0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25])
PR834=np.array([0.000e+00,0.000e+00,1.625e-05,3.725e-05,5.075e-05,1.045e-04,2.021e-04,4.152e-04,8.700e-04,1.719e-03,3.322e-03,5.791e-03,9.621e-03,1.458e-02,2.070e-02,2.809e-02,3.580e-02,4.439e-02,5.297e-02,6.339e-02,7.251e-02,8.394e-02,9.440e-02,1.040e-01])
PR756=np.array([0.000e+00,2.292e-06,2.860e-06,5.540e-06,2.254e-05,4.791e-05,9.905e-05,2.456e-04,5.814e-04,1.239e-03,2.155e-03,4.017e-03,6.491e-03,9.774e-03,1.304e-02,1.746e-02,2.292e-02,2.768e-02,3.488e-02,4.217e-02,4.990e-02,5.858e-02,6.720e-02,7.686e-02])
PR694=np.array([0.000e+00,1.208e-06,1.809e-06,5.162e-06,7.811e-06,2.294e-05,5.889e-05,1.848e-04,4.579e-04,1.014e-03,1.805e-03,2.775e-03,4.502e-03,6.898e-03,9.020e-03,1.224e-02,1.591e-02,2.042e-02,2.500e-02,3.104e-02,3.632e-02,4.351e-02,5.129e-02,5.861e-02])
PR620=np.array([0.000e+00,0.000e+00,8.153e-07,4.092e-07,5.847e-06,1.179e-05,3.800e-05,1.211e-04,3.312e-04,6.216e-04,1.004e-03,1.776e-03,2.695e-03,4.036e-03,5.845e-03,7.975e-03,1.031e-02,1.306e-02,1.658e-02,2.028e-02,2.559e-02,3.065e-02,3.606e-02,4.044e-02])
PR556=np.array([0.000e+00,0.000e+00,3.421e-07,8.726e-07,5.497e-06,9.100e-06,3.303e-05,8.978e-05,1.981e-04,3.459e-04,6.001e-04,1.091e-03,1.747e-03,2.385e-03,3.559e-03,5.076e-03,6.457e-03,9.015e-03,1.154e-02,1.480e-02,1.676e-02,2.245e-02,2.628e-02,3.076e-02])
PR463=np.array([0.000e+00,0.000e+00,0.000e+00,0.000e+00,2.431e-06,0.000e+00,5.422e-06,3.954e-05,8.755e-05,1.591e-04,2.260e-04,4.346e-04,7.087e-04,1.042e-03,1.497e-03,2.226e-03,3.134e-03,4.098e-03,5.795e-03,7.445e-03,8.582e-03,1.236e-02,1.435e-02,1.740e-02])
yscaleit = 1.5e2
plt.plot(Sall,np.log10(yscaleit*PR463),color=hotclr[2])
plt.plot(Sall,np.log10(yscaleit*PR556),color=hotclr[3])
plt.plot(Sall,np.log10(yscaleit*PR620),color=hotclr[4])
plt.plot(Sall,np.log10(yscaleit*PR694),color=hotclr[5])
#plt.plot(Sall,np.log10(yscaleit*PR756),color=hotclr[6])
plt.plot(Sall,np.log10(yscaleit*PR834),color=hotclr[7])
plt.legend(bbox_to_anchor=(0.25,0.85))
plt.show()
#plt.savefig()

##plt.show()
#plt.savefig()

plt.show()
#plt.savefig()
