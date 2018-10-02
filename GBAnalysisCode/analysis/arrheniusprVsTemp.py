
# TSharp
# 2016.12.20

# Script 

import numpy as np
import matplotlib.pyplot as plt

alummelttemp = 926.

temperatureLUT = 1.0*np.zeros(2000)
temperatureLUT[40] = 463.
temperatureLUT[41] = 556.
temperatureLUT[42] = 620.
temperatureLUT[43] = 694.
temperatureLUT[44] = 756.
temperatureLUT[45] = 834.
temperatureLUT[46] = 310.
temperatureLUT[50] = 463.
temperatureLUT[51] = 556.
temperatureLUT[52] = 620.
temperatureLUT[53] = 694.
temperatureLUT[54] = 756.
temperatureLUT[55] = 834.
temperatureLUT[56] = 310.
temperatureLUT[60] = 463.
temperatureLUT[61] = 556.
temperatureLUT[62] = 620.
temperatureLUT[63] = 694.
temperatureLUT[64] = 756.
temperatureLUT[65] = 834.
temperatureLUT[66] = 310.
temperatureLUT[90] = 463.
temperatureLUT[91] = 556.
temperatureLUT[92] = 620.
temperatureLUT[93] = 694.
temperatureLUT[94] = 756.
temperatureLUT[95] = 834.
temperatureLUT[96] = 310.
# P_R(CentroSymm) datasets generated from Qpredictive.py
temperatureLUT[20] = 463.
temperatureLUT[21] = 556.
temperatureLUT[22] = 620.
temperatureLUT[23] = 694.
temperatureLUT[24] = 756.
temperatureLUT[25] = 834.
temperatureLUT[26] = 310.
# P_R(Energy) datasets generated from Qpredictive.py
temperatureLUT[30] = 463.
temperatureLUT[31] = 556.
temperatureLUT[32] = 620.
temperatureLUT[33] = 694.
temperatureLUT[34] = 756.
temperatureLUT[35] = 834.
temperatureLUT[36] = 310.
# P_R(CentroSymm_InstantOf100fsavgdpositions) datasets generated from Qpredictive.py
temperatureLUT[80] = 463.
temperatureLUT[81] = 556.
temperatureLUT[82] = 620.
temperatureLUT[83] = 694.
temperatureLUT[84] = 756.
temperatureLUT[85] = 834.
temperatureLUT[86] = 310.
# P_R(VoroVolumeDeviationOf100fsavgdpositions)
temperatureLUT[100] = 463.
temperatureLUT[101] = 556.
temperatureLUT[102] = 620.
temperatureLUT[103] = 694.
temperatureLUT[104] = 756.
temperatureLUT[105] = 834.
temperatureLUT[106] = 310.
# P_R(PEOf100fsavgdpositions)
temperatureLUT[110] = 463.
temperatureLUT[111] = 556.
temperatureLUT[112] = 620.
temperatureLUT[113] = 694.
temperatureLUT[114] = 756.
temperatureLUT[115] = 834.
################################
# Main
#################################


# Used 20161225
datasets=[40,41,42,43,44,45]
softnesses = range(19,25) 
numLowestTempsToFit=4
consideronlytype=-1
folder2016=1
# Should only consider subset of particles?  -1 (consider all)  0 (disordered)  1 (FCC)

# Used 20170407
datasets=[40,41,42,43,44,45]
softnesses = range(15,29)
numLowestTempsToFit=4
consideronlytype=-1
folder2016=1

# Used 20170420
datasets=[50,51,52,53,55]
softnesses = range(18,29)
numLowestTempsToFit=4
consideronlytype=-1
folder2016=1

# Used 20170424
datasets=[50,54]#,51,52,53,54,55]
#softnesses = range(100)
softnesses = [44,48]+range(52,70,2)+[72,76]
numLowestTempsToFit=4
consideronlytype=-1
folder2016=1

# Used 20170502
datasets=[50,51,52,53,54,55]
softnesses = range(48,68,1)#2#1#2
numLowestTempsToFit=4
consideronlytype=-1
folder2016=1

"""
# Used 20170518
datasets=[20,21,22,23,24,25]
# actually CSymmetry param
softnesses = [9]+range(12,38,4)+[39]
numLowestTempsToFit=4
consideronlytype=-1
folder2016=0

# Used 20170518
datasets=[20,21,22,23,24,25]
# actually CSymmetry param
softnesses = range(12,28,2)+[28,32]
softnesses = range(12,33,1)
numLowestTempsToFit=4
consideronlytype=-1
folder2016=0

# Used 20170531
datasets=[30,31,32,33,34,35]
# actually local atomic energy
softnesses = range(12,36,2)
numLowestTempsToFit=4
consideronlytype=-1
folder2016=0



# Used 20170601 - 20170629 used for supplementary material dE vs S/CS plot
datasets=[80,81,82,83,84,85]
# actually instantaneous CSymmetry (single CS value calculated from the 100-fs avg'd positions
softnesses = range(8,36,3)
#softnesses = range(12,28,2)+[28,32]
#softnesses = range(12,33,1)
numLowestTempsToFit=4
consideronlytype=-1
folder2016=0
mean_soft_per_X_file='mean_soft_per_CS_T694K_CNA-1.txt'
mean_soft_per_X_file='mean_soft_per_CS_CNA-1.txt'


# Used 20170605
# Experimenting with CS-averaged-over-100fs and instant, and averaged Energy, etc
datasets=[90,91,92,93,94,95]
softnesses = range(0,41,3)#2#1#2
numLowestTempsToFit=4
consideronlytype=10 # 10 and above is near-crystal filtering (super strict filtering.)
folder2016=0


# Used 20170614
datasets=[56,50,51,52,53,54,55]
#datasets=[96,90,91,92,93,94,95]
softnesses = range(0,14,3)#2#1#2
numLowestTempsToFit=4
consideronlytype=-1 # 10 # 10 and above is near-crystal filtering (super strict filtering.)
folder2016=0
#datasets=[80,81,82,83,84,85]
#datasets=[50,51,52,53,54,55]
#datasets=[40,41,42,43,44,45]
#datasets=[80,81,82,83,84,85]
#softnesses = range(48,68,1)#2#1#2
#softnesses = range(48,68,2)#2#1#2


# Used 20170627. Also 6/29/17 for supplementary material dE vs E/CS plot
datasets=[50,51,52,53,54,55]
#datasets=[96,90,91,92,93,94,95]
softnesses = range(10,31,1)#2#1#2
numLowestTempsToFit=4
consideronlytype=-1 # 10 # 10 and above is near-crystal filtering (super strict filtering.)
folder2016=0
datasets=[100,101,102,103,104,105]#106 first



# Used 20171219 Voro Volume
datasets=[100,101,102,103,104,105]#106 first
softnesses = range(19,35,1)#2#1#2
numLowestTempsToFit=4
consideronlytype=-1 # -1 # 10 # 10 and above is near-crystal filtering (super strict filtering.)
folder2016=0


# Used 20171220 PE of 100-fs avg'd
datasets=[110,111,112,113,114,115]#116 first
softnesses = range(12,30,2)
softnesses = range(10,30,2)
numLowestTempsToFit=4
consideronlytype=-1 # -1 # 10 # 10 and above is near-crystal filtering (super strict filtering.)
folder2016=0
mean_soft_per_X_file='mean_soft_per_PEavgpos_CNA-1.txt'
"""

# Used 20180109.
datasets=[50,51,52,53,54,55]
#datasets=[96,90,91,92,93,94,95]
softnesses = range(10,30,2)#2#1#2
numLowestTempsToFit=4
consideronlytype=-1 # 10 # 10 and above is near-crystal filtering (super strict filtering.)
folder2016=0

# Used 20180419.
datasets=[56,50,51,52,53,54,55]
#datasets=[96,90,91,92,93,94,95]
softnesses = range(10,30,2)#2#1#2
numLowestTempsToFit=5
consideronlytype=-1 # 10 # 10 and above is near-crystal filtering (super strict filtering.)
folder2016=0

plotstyle=1
if plotstyle == 1:
    #myfig = plt.figure(1,figsize=(6,9))
    #myfig = plt.figure(1,figsize=(6,7))
    myfig = plt.figure(1,figsize=(6,6.33))
    axs = myfig.add_subplot(111)#plt.subplot(1,1,1)#, sharex=True)
    axs.set_yscale('log')

if plotstyle == 2:
    myfig = plt.figure(1,figsize=(4,2))
    axs = myfig.add_subplot(111)#plt.subplot(1,1,1)#, sharex=True)

if plotstyle == 3:
    myfig = plt.figure(1,figsize=(6,7))
    axs = myfig.add_subplot(111)#plt.subplot(1,1,1)#, sharex=True)

if plotstyle == 4:
    myfig = plt.figure(1,figsize=(6,6))
    myfig = plt.figure(1,figsize=(6,4))
    axs = myfig.add_subplot(111)#plt.subplot(1,1,1)#, sharex=True)

softnesses.reverse()
Sigma_vs_S = np.zeros(len(softnesses))
DEnergyOverkB_vs_S = np.zeros(len(softnesses))



print("We now assume all files are indexed by ID, and we do no additional sorting")

for iterofsoft, softnessdataindex in enumerate(softnesses):
  temperatures = np.array(0)
  frachopping = np.array(0)
  numhopping = np.array(0)

  for idataset in datasets:
    print("dataset =",idataset)

    ####################################################
    # Input parameters
    ####################################################
    
    
    # Plotting setup and colors
    savepdffile=1
    #hotclr=['purple','b','green','orange','r','k','teal','gold'
    #hotclr=['purple','violet','navy','b','teal','darkgreen','green','orange','salmon','gold','darkgoldenrod','r','gray','chocolate','k']
    #hotclr=['#0000ff','#0022dd','#0044bb','#006699','#008877','#00aa55','#00cc33','#00ee11','#00ff00','#22dd00','#44bb00','#669900','#887700','#aa5500','#cc3300','#ee1100']
    #hotclr=['#0000ff','#0022dd','#0044bb','#006699','#008877','#00aa55','#00cc33','#00ee11','#00ff00','#22dd00','#44bb00','#669900','#887700','#aa5500','#cc3300','#ee1100']
    #hotclr=['#0000ff','#0044ff','#0077ee','#0099cc','#00aabb','#00bbaa','#00cc99','#00ee77','#00ff44','#00ff00',
    #        '#00ff00','#44ff00','#77ee00','#99cc00','#aabb00','#bbaa00','#cc9900','#ee7700','#ff4400','#ff0000']
    hotclr=['#0000dd','#0044dd','#0077cc','#0088bb','#0099aa','#00aa99','#00bb88','#00cc77','#00dd44','#00dd00',
            '#00dd00','#44dd00','#77cc00','#88bb00','#99aa00','#aa9900','#bb8800','#cc7700','#dd4400','#dd0000']
    hotclr.extend(hotclr)
    hotclr.extend(hotclr)
    hotclr.extend(hotclr)
    
    # Restore data
    if folder2016:
      # normal style (softness)
      outfilenamebase = 'outfiles/pR_dset'+str(idataset)
    else:
      # X style (softness, centro symm, averaged-energy)
      outfilenamebase = 'outfiles/X'+str(idataset)+'_pR_dset'+str(idataset)
    #print("Opening: "+ outfilenamebase+'_bincenters_CNA'+str(consideronlytype)+'.txt')
    bin_centers = np.genfromtxt(outfilenamebase+'_bincenters_CNA'+str(consideronlytype)+'.txt')
    numberineach_bin = np.genfromtxt(outfilenamebase+'_numinbin_CNA'+str(consideronlytype)+'.txt')
    fractionhopping_vs_bin = np.genfromtxt(outfilenamebase+'_frachopping_CNA'+str(consideronlytype)+'.txt')
    numberhopping_vs_bin = fractionhopping_vs_bin * numberineach_bin
    numberhopping_vs_bin += 1e-12
    print(bin_centers, 'are bin centers')

    temperatures = np.append(temperatures,temperatureLUT[idataset])
    frachopping = np.append(frachopping,fractionhopping_vs_bin[softnessdataindex])
    numhopping = np.append(numhopping,numberhopping_vs_bin[softnessdataindex])
        

  temperatures = temperatures[1::]
  frachopping = frachopping[1::]
  numhopping = numhopping[1::]

  #print(temperatures)
  #print(frachopping+ 1e-12)
  #print(bin_centers[softnessdataindex])
  #print("S")
  #plt.semilogy(...)
  if plotstyle == 1:
      axs.errorbar(alummelttemp/temperatures, frachopping+ 1e-12, yerr=frachopping/np.sqrt(numhopping),
                    marker='o',
                    label=str(bin_centers[softnessdataindex]), color=hotclr[2*iterofsoft])

  # Extract the implied energy barriers and distribution from the Arrhenius-like curves
  # Adding this part 20170407
  # P_R = exp( Sigma - DEnergy / kT)
  # log_e P_R = Sigma - 1/T DEnergy/k
  # Calculate best-fit to log
  x=1.0/temperatures
  y=np.log(frachopping+ 1e-12)
  x=x[0:numLowestTempsToFit]
  y=y[0:numLowestTempsToFit]
  coefs = np.polyfit(x, y, 1)
  Sigma_vs_S[iterofsoft] = coefs[1]
  DEnergyOverkB_vs_S[iterofsoft] = -coefs[0]

  # Plot the fit line
  exTmOverT = np.array([-20,-10,-1,-0.1,0.1,0.5,1.0,1.5,2.0,2.5,5])
  ex1OverT = exTmOverT / alummelttemp
  # Plot on the Tm/T axes
  plt.plot(exTmOverT, np.exp(coefs[1] +coefs[0]*ex1OverT), '--', color=hotclr[2*softnesses.index(softnessdataindex)], linewidth=1)



if plotstyle == 1:
    if 1:
     plt.axis([1,3.05,1*7e-5,1*.2]) # 1+2.05
     plt.legend(title='S',loc=4,prop={'size':7},bbox_to_anchor=(0.9,0.03))
    elif 1:
     plt.axis([1,2.05,1*7e-5,1*.2]) # 1+2.05
     plt.legend(title='S',loc=4,prop={'size':7},bbox_to_anchor=(0.9,0.03))
    else:
     plt.axis([1,0+2.05,1e-3,.2]) # 2.05
     plt.legend(title='S',loc=4,prop={'size':7},bbox_to_anchor=(0.2,0.03))
    plt.xlabel(r'$T_m$'+'/'+'$T$')
    plt.grid(b=True, which='major')
    plt.ylabel(r'Pr(S,T)') # (r'$p_{hop}$')
    if savepdffile:
        plt.savefig("del_PE_fig20180420arrheniusprVsTempCNA"+str(consideronlytype)+".eps")
    plt.show()
    plt.clf()

if plotstyle == 2:
    #print("Caution, Using bin_centers of the last-loaded file to know what softness.  This assumes that theyre all the same.")
    kB = 8.6173e-2 # meV/K
    #print(bin_centers[softnesses],"softnesses")
    #print(DEnergyOverkB_vs_S * kB, "DE/kB")
    #print(Sigma_vs_S,"e^Sigma")
    axs.set_yscale('linear')
    #plt.plot(bin_centers[softnesses], np.exp(Sigma_vs_S), 'o-', color='black', linewidth=2, markersize=15)
    plt.legend(title='S',loc=4,prop={'size':7},bbox_to_anchor=(0.9,0.03))
    plt.xlabel(r'S')
    plt.grid(b=True, which='major')
    plt.ylabel('Sigma*kB(black)*alummelttemp (divided by 1) [meV]  and DEnergy(blue) [meV]')
    plt.plot(bin_centers[softnesses], Sigma_vs_S*kB*alummelttemp / 1.0, 'o-', color='black', linewidth=2, markersize=9)
    plt.plot(bin_centers[softnesses], DEnergyOverkB_vs_S * kB, 'o-', color='blue', linewidth=2,markersize=9)
    plt.axis([bin_centers[np.min(softnesses)],bin_centers[np.max(softnesses)],0,300])
    print([np.min(softnesses),np.max(softnesses),0,300])
    if savepdffile:
        plt.savefig("del_fig20180420DEnergyCNA"+str(consideronlytype)+".eps")
    plt.show()


if plotstyle == 3:
    #myfig = plt.figure(1,figsize=(4,4))
    #axs = myfig.add_subplot(111)#plt.subplot(1,1,1)#, sharex=True)
    #plt.xlabel('Voro volume [Ang^3]')
    #plt.xlabel('Soft')
    plt.xlabel('CS')
    #plt.xlabel('Pot. Energy [eV]')
    plt.ylabel('Sigma and DEnergy/kB/Tmelt(blue) ')
    plt.plot(bin_centers[softnesses], Sigma_vs_S, 'o-', color='black', linewidth=2, markersize=9)
    plt.plot(bin_centers[softnesses], DEnergyOverkB_vs_S / alummelttemp, 'o-', color='blue', linewidth=2,markersize=9)
    #plt.axis([bin_centers[np.min(softnesses)],bin_centers[np.max(softnesses)],-4.5,3.5])
    plt.grid(b=True, which='major')
    #plt.axis([-3,6,-10,5])
    #plt.axis([-.5,4.5,-5,5])
    plt.axis([bin_centers[min(softnesses)],bin_centers[max(softnesses)],-4,3])
    if savepdffile:
        plt.savefig("20180420SigmaEoverkBTCNA"+str(consideronlytype)+"_vs_CS.eps")
        #plt.savefig("20171220SigmaEoverkBTCNA"+str(consideronlytype)+"_vs_Soft.pdf")
        #plt.savefig("20171220SigmaEoverkBTCNA"+str(consideronlytype)+"_vs_VoroVol.pdf")
        #plt.savefig("20171220SigmaEoverkBTCNA"+str(consideronlytype)+"_vs_PE.pdf")
    plt.show()



# Centro symmetry on Softness axes (using mean softness of each CS value)
if plotstyle == 4:
    # Take the Energy barrier vs CS and plot Energy barrier against MeanSoftnessOf(CS)
    data = np.genfromtxt(mean_soft_per_X_file)
    #xaxisofCSwas to xaxis is Softness, is
    inCSbins = data[:,0]
    if len(inCSbins) != len(bin_centers):
        print("CS bins in E vs CS are not the same as in Softness vs CS")
        print(inCSbins)
        print(bin_centers)
        import sys
        sys.quit
    for iel in np.arange(len(inCSbins)):
        if np.abs(inCSbins[iel] - bin_centers[iel]) > 1e-3:
            print("CS bins in E vs CS are not the same as in Softness vs CS")
            print(inCSbins)
            print(bin_centers)
            import sys
            sys.quit

    inSoftvalues = data[:,1]

    axs.set_xlabel("Softness")
    #ax2 = axs.twiny()
    #new_tick_locations = inSoftvalues
    #ax2.set_xlim(axs.get_xlim())
    axs.set_xticks(inSoftvalues)
    axs.set_xticklabels(inCSbins)
    axs.set_xlabel("CS values of this mean Softness")
    plt.grid(b=True, which='major')
    plt.plot(inSoftvalues[softnesses], Sigma_vs_S, 's-', color='blue', linewidth=2, markersize=9)
    plt.plot(inSoftvalues[softnesses], DEnergyOverkB_vs_S / alummelttemp, 's-', color='blue', linewidth=2,markersize=9)
    plt.axis([-3,6,-10,5])
    plt.axis([-.5,4.5,-5,5])
    if savepdffile:
        plt.savefig("delete_fig20180420SigmaEoverkBTCNA"+str(consideronlytype)+".eps")
    plt.show()


