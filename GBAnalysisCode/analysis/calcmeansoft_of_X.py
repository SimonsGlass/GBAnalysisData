

# TSharp
# 20170605

# Script 
# Find the effective softness at each CS value
# This could be the mean softness, <S>_CS = Integral dS S P(S|CS)
# or S' such that <S>_CS = Integral dS S P(S|CS) PR(S)

import numpy as np
import matplotlib.pyplot as plt

def getjth_dump_col(filename, Nrecords, j, numheaderlines):

        print("Reading "+str(Nrecords)+" values in col "+str(j)+" from file ..." + filename[-50::])

	with open(filename,'r') as file:
                # Skip the header
                for iline in range(numheaderlines-1):
                        file.readline()
                if numheaderlines > 0:
                        print("Last header line (if any): ",file.readline())
                        
                # Read in the whole file, then crop
                if 1:
                    # Split each remaining line by spaces, extract jth col, cast as float
                    myfloatcol = [float((line.split())[j]) for line in file]
                    if Nrecords == -1:
                            Nrecords = len(myfloatcol)
                    myfloatcol = myfloatcol[0:Nrecords]

        print("col ",j," is ", myfloatcol[0], myfloatcol[1],myfloatcol[2], myfloatcol[3])
        print("")
        if (len(myfloatcol) != Nrecords):
                print("Length read in from file, "+str(len(myfloatcol))+
                      " doesnt match Input parameter "+str(Nrecords))
                import sys
                sys.exit(0)
        
        return myfloatcol



################################
# Main
#################################

numheaderlines=1

if 0:
    # Define the softness bins [-4.0 +/- 0.25], [-3.5 +/- 0.25], ... [4.0 +/- 0.25]
    binrange=10
    binmid=2
if 0:
    # Define the CS bins [-4.0 +/- 0.25], [-3.5 +/- 0.25], ... [4.0 +/- 0.25]
    binrange=16
    binmid=8
    files = ["datahistsoft/dset53/Centrosymm_vs_soft_Al694K_t123000_seeCNAall.txt"]
    outfilenamebase = 'mean_soft_per_CS_T694K'
    files = ["datahistsoft/dset50/Centrosymm_vs_soft_Al463K_t143000_seeCNAall.txt"]
    outfilenamebase = 'mean_soft_per_CS'
if 1:
    # Define the PE bins
    binrange=.6
    binmid=-3.2
    files = ["datahistsoft/dset50/PEavgpos_vs_soft_Al463K_t143000_seeCNAall.txt"]
    outfilenamebase = 'mean_soft_per_PEavgpos'

if 1:
    nbins=41#121#41
    binwidth = 1.0*binrange/(nbins-1) # bin range is from lo-bin-center to hi-bin-center
    nbinedges = nbins + 1
    binedgemin = binmid - binrange/2.0 - binwidth/2.0 
    bin_edges = binedgemin + binwidth*np.arange(nbinedges)
    bin_centers = (bin_edges[0:nbins] + bin_edges[1:1+nbins]) / 2.0


    ####################################################
    # Input parameters
    ####################################################
    # Should only consider subset of particles?  -1 (consider all)  0 (disordered, CNA!=FCC)  1 (CNA=FCC)
    consideronlytype = -1 # considering all atoms would give higher Q's but P_R(mu-sigma) = undefined
    # 10 (very disordered) 11 (3.4A near CNA=FCC/HCP) # Super strict filter # Use the aggresive near-FCC filtering  
    # Time (Number of frames * file spacign) between softness and last time we look for a hop
    #doshow =1
    #regenerate=1 # 0 or 1
    
    # Plotting setup and colors
    #axs = plt.subplot(1,1,1)#, sharex=True)
    fig1 = plt.figure(1,figsize=(6,6))
    axs = fig1.add_subplot(111)
    hotclr=['purple','b','green','orange','r','k']#'teal','gold'
    hotclr+=hotclr
    print(bin_edges,"bin edges")
    print(bin_centers,"bin centers")
    #import sys
    #sys.exit(0)

    numparticlestoanalyze = -1 # 50000 # -1 is all
    softness = getjth_dump_col(files[0], numparticlestoanalyze, 1, numheaderlines)
    centrosym = getjth_dump_col(files[0], numparticlestoanalyze, 0, numheaderlines)
    print("Opening softness-centrosymm file")
            

    print("==>  ", np.mean(softness), np.sqrt(np.var(softness)), "are mean and std dev of X of all input particles")
    ####################################################
    # Assign each particle a bin based on its softness
    ####################################################
    print("Binning by centrosym")
    bin_number_vs_id = (np.array(centrosym)-binedgemin)/binwidth
    bin_number_vs_id = np.clip(bin_number_vs_id.astype(int), 0, nbins-1)
    
    sumofSs_vs_bin = np.zeros(nbins)
    numofSs_vs_bin = np.zeros(nbins)
    meanofSs_vs_bin = np.zeros(nbins)
    sumofS2s_vs_bin = np.zeros(nbins)
    meanofS2s_vs_bin = np.zeros(nbins)
    stddevS_vs_bin = np.zeros(nbins)

    for i in np.arange(len(softness)):
        sumofS2s_vs_bin[bin_number_vs_id[i]] += (softness[i])**(2.0)
        sumofSs_vs_bin[bin_number_vs_id[i]] += softness[i]
        numofSs_vs_bin[bin_number_vs_id[i]] += 1

    for i in np.arange(nbins):
        if numofSs_vs_bin[i] != 0:
            meanofSs_vs_bin[i] = sumofSs_vs_bin[i] / numofSs_vs_bin[i]
            meanofS2s_vs_bin[i] =sumofS2s_vs_bin[i] / numofSs_vs_bin[i]
            stddevS_vs_bin[i] = np.sqrt(meanofS2s_vs_bin[i] - (meanofSs_vs_bin[i])**2.0) 

    axs.errorbar(bin_centers, meanofSs_vs_bin, yerr=stddevS_vs_bin)
    print(stddevS_vs_bin)

    
    outputfile = open(outfilenamebase+'_CNA'+str(consideronlytype)+'.txt','w')
    for i in np.arange(nbins):
        outputfile.write(str(bin_centers[i])+' '+str(meanofSs_vs_bin[i])+'\n')

    plt.savefig(outfilenamebase+'_CNA'+str(consideronlytype)+'.pdf')
    plt.show()




"""
        

    for 
    
        
    #for iline in range(nparticles*len(frames)):
    #    print("I think", softness[iline], " belongs in bin ", bin_number_vs_fr_id[iline], " which is " , bin_edges[bin_number_vs_fr_id[iline]], "-to-", bin_edges[1+bin_number_vs_fr_id[iline]])

    if 1:
        ####################################################
        # Learn which of the atoms hopped in the next frames
        ####################################################
        print("Open files to see CNA and if rearrange")
        particle_did_hop = np.zeros(numparticlestoanalyze*len(frames))
        pc_hopcutoff = 1.00
        
        cna_vs_fr_id = []
        for inum,frame in enumerate(frames):  
            # Get lattice type from first frame
            CNAfile = path+'CNAtrainingAve.0'+str(frame)+'.cfg'
            cnacol=0
            cnaheader=9
            if consideronlytype >= 10:
                temptemp = path[-5:-2]
                print(temptemp)
                CNAfile = path+'NearCrystalLabel_Al'+str(temptemp)+'K_t'+str(frame)+'.txt'
                cnacol=1
                cnaheader=1
            print("Opening ",CNAfile[-50::])
            cna_vs_id = getjth_dump_col(CNAfile, numparticlestoanalyze, cnacol, cnaheader)
            cna_vs_fr_id.extend(cna_vs_id)
            
            # HERE IS WHERE I LEFT OFF

            print("Checking mean softness")
            for delayframe in frame + framespacing*np.arange(1 + delaywindow/framespacing):
                phopfile = path+phopfileprefix+str(delayframe)+phopfilesuffix
                print("Opening ",phopfile[-50::])
                phop_vs_id = getjth_dump_col(phopfile, numparticlestoanalyze, 5, 9)
                #print(phopfile, ' opened and phop=',phop_vs_id)
                
                for iparticle in range(numparticlestoanalyze):
                    if (phop_vs_id[iparticle] > pc_hopcutoff):
                            particle_did_hop[iparticle+inum*numparticlestoanalyze]=1

        
        ####################################################
        # Get the percentage that hopped vs bin
        ####################################################
        #print('Find number and number hopping in each bin')
        numberineach_bin = 0.0*np.zeros(nbins)
        numberhopping_vs_bin = 0.0*np.zeros(nbins)
        fractionhopping_vs_bin = 0.0*np.zeros(nbins)
        for iframeparticle in range(numparticlestoanalyze*len(frames)):
            if (consideronlytype == -1 or (cna_vs_fr_id[iframeparticle] == consideronlytype)):
                numberineach_bin[bin_number_vs_fr_id[iframeparticle]] += 1.0
                if particle_did_hop[iframeparticle]:
                    numberhopping_vs_bin[bin_number_vs_fr_id[iframeparticle]] += 1.0
        
        #print(numberineach_bin,'number of particles in each bin')
        #print(numberhopping_vs_bin,'number of particles that hopped vs bin')
        
        for ibin in range(nbins):
            if numberhopping_vs_bin[ibin] > 0: 
                fractionhopping_vs_bin[ibin] = numberhopping_vs_bin[ibin] / numberineach_bin[ibin]



    
    ####################################################
    # Save
    ####################################################
    
    outfilenamebase = 'outfiles/X'+str(Xselection)+'_pR_dset'+dset
    if regenerate == 1:
        outputfile = open(outfilenamebase+'_bincenters_CNA'+str(consideronlytype)+'.txt','w')
        for item in bin_centers:
            outputfile.write("%s\n" % str(item))
        outputfile = open(outfilenamebase+'_numinbin_CNA'+str(consideronlytype)+'.txt','w')
        for item in numberineach_bin:
            outputfile.write("%s\n" % str(item))
        outputfile = open(outfilenamebase+'_frachopping_CNA'+str(consideronlytype)+'.txt','w')
        for item in fractionhopping_vs_bin:
            outputfile.write("%s\n" % str(item))
    else:
        bin_centers = np.genfromtxt(outfilenamebase+'_bincenters_CNA'+str(consideronlytype)+'.txt')
        numberineach_bin = np.genfromtxt(outfilenamebase+'_numinbin_CNA'+str(consideronlytype)+'.txt')
        fractionhopping_vs_bin = np.genfromtxt(outfilenamebase+'_frachopping_CNA'+str(consideronlytype)+'.txt')
        numberhopping_vs_bin = fractionhopping_vs_bin * numberineach_bin
    
    
    ####################################################
    # Plot
    ####################################################
    numberhopping_vs_bin += 1e-12
    print(np.sum(numberhopping_vs_bin)," is num hopping", np.sum(numberineach_bin)," is num in bins")
    #axs.errorbar(bin_centers, fractionhopping_vs_bin, yerr=fractionhopping_vs_bin/np.sqrt(numberhopping_vs_bin), 
    axs.errorbar(bin_centers, fractionhopping_vs_bin, yerr=1.0/numberineach_bin, 
                    label=path[-13::], color=hotclr[np.mod(Xselection,6)])
    print(bin_centers)
    print(fractionhopping_vs_bin)
    print(numberineach_bin)
    #plt.semilogy(...)






savelabel='_20170601E' # '_20170503'
axs.set_yscale('log')
#plt.axis([-10,10,1e-7,1e-0])
#plt.axis([-2,7,1e-6,1e-0])
plt.axis([min(bin_centers), max(bin_centers), 1e-6,1e-0])
plt.legend(title='Temperature ',loc=4,prop={'size':12},bbox_to_anchor=(0.85,0.02))
plt.xlabel(r'S')
plt.grid(b=True, which='major')
plt.ylabel(r'Pr(S)') # (r'$p_{hop}$')
plt.yticks([1*1e-6,2*1e-6,3*1e-6,4*1e-6,5*1e-6,6*1e-6,7*1e-6,8*1e-6,9*1e-6,
            1*1e-5,2*1e-5,3*1e-5,4*1e-5,5*1e-5,6*1e-5,7*1e-5,8*1e-5,9*1e-5,
            1*1e-4,2*1e-4,3*1e-4,4*1e-4,5*1e-4,6*1e-4,7*1e-4,8*1e-4,9*1e-4,
            1*1e-3,2*1e-3,3*1e-3,4*1e-3,5*1e-3,6*1e-3,7*1e-3,8*1e-3,9*1e-3,
            1*1e-2,2*1e-2,3*1e-2,4*1e-2,5*1e-2,6*1e-2,7*1e-2,8*1e-2,9*1e-2,
            1*1e-1,2*1e-1,3*1e-1,4*1e-1,5*1e-1,6*1e-1,7*1e-1,8*1e-1,9*1e-1])
if 1:
    plt.savefig(outfilenamebase+'_Q_of_X'+str(Xselection)+'_'+str(consideronlytype)+savelabel+'.pdf')
if doshow:
    plt.show()


# Histogram
avg= np.sum(bin_centers*numberineach_bin)/np.sum(numberineach_bin)
print("mean of CNA-filtered:", avg)
avg_of_sq = np.sum(bin_centers*bin_centers*numberineach_bin)/np.sum(numberineach_bin)
print("stdv of CNA-filtered:", np.sqrt(avg_of_sq -avg*avg))
plt.plot(bin_centers, numberineach_bin)
plt.grid(b=True, which='minor')
if doshow:
    plt.show()
"""
