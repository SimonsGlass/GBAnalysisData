
# TSharp
# 2016.09.18
# 2017.04.21

# Script 
# to read the Softness of each particle at a particular set of frames.
# Then we assign each particle to a bin based on the Softness.
# We then find the fraction of particles in each bin that rearrange within DELAY number of frames.
# That is, we open files starting from the softness frame onward to DELAY frames later, and see if phop exceeds a threshold.

# if enable FCC version, only considers atoms that are labeled with CNA = FCC (labelled 1)

import numpy as np
import matplotlib.pyplot as plt


########################################################################
# Function to read jth col from dump file. 
# Indexing stats at 0, but 1st column white space counts.
# Here used to read phop, in 5th col. (id type x y z phop)
########################################################################
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
                    myfloatcol = myfloatcol[0:Nrecords]
                else:
                    # This is untested
                    myfloatcol = []
                    for iline, line in enumerate(file):
                        if iline > Nrecords:
                                break
                        myfloatcol.extend(float((line.split())[j]))

        print("col ",j," is ", myfloatcol[0], myfloatcol[1],myfloatcol[2], myfloatcol[3])
        print("")
        if (len(myfloatcol) != Nrecords):
                print("Length read in from file, "+str(len(myfloatcol))+
                      " doesnt match Input parameter "+str(Nrecords))
                import sys
                sys.exit(0)
        
        return myfloatcol

        # Old method was 11x slower:
	#with open(filename,'r') as file:
        #        tempStructure = [y for y in file]
	#tempStructure = tempStructure[numheaderlines:Nrecords+numheaderlines]
        #print("Splitting structure...")		
	#tempStructure = [[float(y) for y in x.split()] for x in tempStructure]
        # Could sort: tempStructure.sort(key = lambda x: x[0]) # sort by ID
        #return [x[j] for x in tempStructure]




################################
# Main
#################################
datasets = range(40,46)+range(50,66)

#datasets=[55,56] # View now-fixed PBC problem
#datasets=[40,41,42,43,44,45] # Train on CNA=0 (GB) atoms
#datasets =[50,51,53,54,55] # Train on all atoms
#datasets+=[60,61,63,64,65] # Retrain at each temperature
#datasets=[55,56]
#datasets=[54,57,58,59]

datasets=range(50,56)
#print("RESET THE ABOVE LINE")
#datasets=[40, 41, 42, 43, 44, 45]

print("We assume all files data are sorted by atom ID.  We don't sort it.")

for idataset in datasets:
    print("dataset =",idataset)

    #############################
    # Data selection data base
    #############################
    
    # Frame of softness prediction
    # Path to softness file
    # Softness file name
    # PHop file name
    # Number of particles in the soft files.  Assumed that their IDs are 1-N consecutive
    if idataset == 1:
        frames = [1541000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/poly816ggAnneal/100fsframes/smallchunk/'
        softfileprefix = path+'KA_N40337_tag22exc1.00win4sam100slow7000sym00_Softness.'
        #phopfile = path+'smallCustom.'+str(frame)+'.phop.dump.gaveids'
        phopfileprefix = 'smallCustom.'
        phopfilesuffix = '.phop.dump.gaveids'
        nparticles = 40337
    elif idataset == 2:
        frames = [1540500]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/poly816ggAnneal/100fsframes/'
        softfileprefix = path+'KA_N5044374_tag21exc1.00win3sam100000slow7000_Softness.'
        phopfileprefix = 'trainingAnnealCustom.'
        phopfilesuffix = '.dump.phop'
        nparticles = 5044374
    elif idataset == 3:
        frames = [0120600]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'
        softfileprefix = path+'KA_N3873924_tag32exc1.00win2sam3000slow3600_Softness.'
        phopfileprefix = 'trainingAve.'
        phopfilesuffix = '.cfg.phop'
        nparticles = 3873924
    elif idataset == 4:
        frames = [120600]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/smallchunk/'
        softfileprefix = path+'KA_N29538_tag32exc1.00win2sam3000slow3600_Softness.'
        phopfileprefix = 'small.0'
        phopfilesuffix = '.cfg.gaveids.phop'
        nparticles = 29508
    elif idataset == 5:
        frames = [120600]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/smallchunk/'
        softfileprefix = path+'KA_N29043_tag34exc1.00win2sam3000slow3600_Softness.'
        phopfileprefix = 'small.0'
        phopfilesuffix = '.cfg.gaveids.phop'
        nparticles = 29043
    elif idataset == 6:
        frames = [120600]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/smallchunk/'
        softfileprefix = path+'KA_N29043_tag32exc1.00win2sam3000slow3600_Softness.'
        phopfileprefix = 'small.0'
        phopfilesuffix = '.cfg.gaveids.phop'
        nparticles = 29043
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 7:
        frames = [120600]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/questionabletemperatures/poly746ggMishin620K/smallchunk/'
        softfileprefix = path+'KA_N193008_tag35exc1.00win2sam3000slow3600_Softness.'
        phopfileprefix = 'small.0'
        phopfilesuffix = '.cfg.gaveids.phop'
        nparticles = 193008
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 8:
        frames = [120600]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/smallchunk/'
        softfileprefix = path+'KA_N191903_tag36exc1.00win2sam3000slow3600_Softness.'
        phopfileprefix = 'small.0'
        phopfilesuffix = '.cfg.gaveids.phop'
        nparticles = 191903

    # 6 Temperatures, Train on GB only #
    elif idataset == 40:
        frames = [143000]
        #frames = [141000]
        print("RESET THE ABOVE LINE")
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag40exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000 
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 41:
        frames = [133000]
        #frames = [131000]
        print("RESET THE ABOVE LINE")
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag41exc1.00win2sam10000slow1800_Softness.'
        softfileprefix = path+'KA_N3873924_tag40exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 42:
        frames = [133000]
        #frames = [131000]
        print("RESET THE ABOVE LINE")
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag42exc1.00win2sam10000slow1800_Softness.'
        softfileprefix = path+'KA_N3873924_tag40exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 43:
        frames = [123000]
        #frames = [121000]
        print("RESET THE ABOVE LINE")
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag43exc1.00win2sam10000slow1800_Softness.'
        softfileprefix = path+'KA_N3873924_tag40exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 44:
        frames = [123000]
        #frames = [121000]
        print("RESET THE ABOVE LINE")
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag44exc1.00win2sam10000slow1800_Softness.'
        softfileprefix = path+'KA_N3873924_tag40exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 45:
        frames = [143000]
        #frames = [141000]
        print("RESET THE ABOVE LINE")
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag45exc1.00win2sam10000slow1800_Softness.'
        softfileprefix = path+'KA_N3873924_tag40exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
        


    # 6 Temperatures, Train on all #
    elif idataset == 50:
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 51:
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 52:
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 53:
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 54:
        frames = [123000,125000,127000,129000] # 123000 was VERY different than the others, so now rerunning
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 55:
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'

    # Data set before PBC bug fix #
    elif idataset == 56: # data from before PBC bug fix. Compare 56 with 55
        frames = [140900]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'
        softfileprefix = path+'b4PBCbugfix/KA_N3873924_tag55exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 # 200000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'


    # Temporary, compare each subsection of dataset 54
    elif idataset == 57:
        frames = [125000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 58:
        frames = [127000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 59:
        frames = [129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'


    # Retrain at each temperature.  Tag is still 50-55. #
    elif idataset == 60:
        frames = [143000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 61:
        frames = [133000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag51exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 62:
        frames = [133000] # Frame 133000 is missing?
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag52exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 63:
        frames = [123000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag53exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 64:
        frames = [123000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag54exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'
    elif idataset == 65:
        frames = [143000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag55exc1.00win2sam10000slow1800_Softness.'
        phopfileprefix = 'trainingAve.0'
        phopfilesuffix = '.cfg.phop'
        nparticlestotal = 3873924 # 200000
        numparticlestoanalyze = 3873924 #  10000
        CNAfileprefix = 'CNAtrainingAve.0'
        CNAfilesuffix = '.cfg'




    else : 
        print("No idataset =", idataset)
        import sys
        sys.exit(0)
    
    

    ####################################################
    # Input parameters
    ####################################################
    
    # Should only consider subset of particles?  -1 (consider all)  0 (disordered, CNA!=FCC)  1 (CNA=FCC)
    consideronlytype = -1 # considering all atoms would give higher Q's but P_R(mu-sigma) = undefined
    # If use the aggresive near-FCC filtering  # 10 (very disordered) 11 (3.4A near CNA=FCC/HCP)
    # Time (Number of frames * file spacign) between softness and last time we look for a hop
    delaywindow = 0 # fsec (0 is usual.  200 was usual) units is file ids, in this case femtoseconds
    delaywindowoffset = 0 # 200 # fsec. Introduced 2018-05-31. Since phop written to first J window file. Should be reset to 0 to recover previous results
    print("RESET THE ABOVE LINE")
    framespacing = 100
    
    regenerate = 0 # 0 or 1

    # Define the bins [-4.0 +/- 0.25], [-3.5 +/- 0.25], ... [4.0 +/- 0.25]
    if numparticlestoanalyze > 1000000:
            nbins = 101
            binwidth=0.5/((nbins-1)/40)
            
    else:
            nbins = 41
            binwidth = 0.5
    nbinedges = nbins + 1
    binedgemin = -0.5*binwidth*(nbins) # centering the bins at zero
    bin_edges = binedgemin + binwidth*np.arange(nbinedges)
    bin_centers = (bin_edges[0:nbins] + bin_edges[1:1+nbins]) / 2.0
    
    # Plotting setup and colors
    #axs = plt.subplot(1,1,1)#, sharex=True)
    fig1 = plt.figure(1,figsize=(6,6))
    axs = fig1.add_subplot(111)
    hotclr=['purple','b','green','orange','r','k']#'teal','gold'
    hotclr+=hotclr
    #print(bin_edges,"bin edges")
    #print(bin_centers,"bin centers")
    #import sys
    #sys.exit(0)
    
    if regenerate == 1:
        ####################################################
        # Read softness files into softness_by_frame_and_id list
        ####################################################
        soft_vs_fr_id = []
        for frame in frames:
            print("Opening softness file", softfileprefix+str(frame))#+softfilesuffix)
            #softness = getjth_dump_col(softfileprefix+str(frame)+softfilesuffix, numparticlestoanalyze, column, numheaderlines)
            with open(softfileprefix+str(frame)) as f:
                softness = [map(float,num.split()) for num in f.readlines()]
        
            if (len(softness) != numparticlestoanalyze):
                print("Length of Soft file"+str(len(softness))+
                      " doesnt match Input parameter"+str(numparticlestoanalyze))
                import sys
                sys.exit(0)
            soft_vs_fr_id.extend(softness)

        print("==>  ", np.mean(softness), np.sqrt(np.var(softness)), "are mean and std dev of X of all input particles")
        ####################################################
        # Assign each particle a bin based on its softness
        ####################################################
        print("Binning by softness")
        bin_number_vs_fr_id = (np.array(soft_vs_fr_id)-binedgemin)/binwidth
        bin_number_vs_fr_id = np.clip(bin_number_vs_fr_id.astype(int), 0, nbins-1)
        
        #for iline in range(nparticles*len(frames)):
        #    print("I think", softness[iline], " belongs in bin ", bin_number_vs_fr_id[iline], " which is " , bin_edges[bin_number_vs_fr_id[iline]], "-to-", bin_edges[1+bin_number_vs_fr_id[iline]])
        
        ####################################################
        # Learn which of the atoms hopped in the next frames
        ####################################################
        print("Open files to see CNA and if rearrange")
        particle_did_hop = np.zeros(numparticlestoanalyze*len(frames))
        pc_hopcutoff = 1.00
        
        cna_vs_fr_id = []
        for inum,frame in enumerate(frames):  
            # Get lattice type from first frame
            CNAfile = path+CNAfileprefix+str(frame)+CNAfilesuffix
            print("Opening ",CNAfile[-50::])
            cna_vs_id = getjth_dump_col(CNAfile, numparticlestoanalyze, 0, 9)
            cna_vs_fr_id.extend(cna_vs_id)
            
            print("Checking if particles hop")
            for delayframe in frame + delaywindowoffset + framespacing*np.arange(1 + delaywindow/framespacing):
                phopfile = path+phopfileprefix+str(delayframe)+phopfilesuffix
                print("Opening ",phopfile[-50::])
                phop_vs_id = getjth_dump_col(phopfile, numparticlestoanalyze, 5, 9)
                #print(phopfile, ' opened and phop=',phop_vs_id)
                
                for iparticle in range(numparticlestoanalyze):
                    if (phop_vs_id[iparticle] > pc_hopcutoff):
                            particle_did_hop[iparticle+inum*numparticlestoanalyze]=1

        #print('particle_did_hop',particle_did_hop)
        
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

    outfilenamebase = 'outfiles/pR_dset'+str(idataset)
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
                    label=path[-13::], color=hotclr[datasets.index(idataset)])
    #plt.semilogy(...)

    print(path[-13::])
    np.set_printoptions(precision=3)
    print(bin_centers[40:64], fractionhopping_vs_bin[40:64])



savelabel='_2018delete'
axs.set_yscale('log')
plt.axis([-10,10,1e-7,1e-0])
plt.axis([-2,7,1e-6,1e-0])
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
    plt.savefig(outfilenamebase+'_PR_of_S_'+str(consideronlytype)+savelabel+'.pdf')
plt.show()

