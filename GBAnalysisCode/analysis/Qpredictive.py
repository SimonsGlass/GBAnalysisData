
# TSharp
# 2016.09.18
# 2017.04.21

# Script 
# to read the Quantity X of each particle at a particular set of frames.
# Then we assign each particle to a bin based on the X value.
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
#datasets = range(40,46)+range(50,66)

#datasets=[55,56] # View now-fixed PBC problem
#datasets=[40,41,42,43,44,45] # Train on CNA=0 (GB) atoms
#datasets =[50,51,53,54,55] # Train on all atoms
#datasets+=[60,61,63,64,65] # Retrain at each temperature
#datasets=[55,56]
#datasets=[54,57,58,59]

#datasets=range(50,56)

print("We assume all files data are sorted by atom ID.  We don't sort it.")

#for Xselection in [0,1,2,3,4,5,20,21,22,23,24,25,30,31,32,33,34,35,80,81,82,83,84,85]:#80,81,82,83,84,85]:
#for Xselection in [100,101,102,103,104,105]:#80,31,32,38,34,35]:
for Xselection in [56]:#[56,50,51,52,53,54,55]:#[56]#[91,92,93,94,95]:#80,31,32,38,34,35]:
#for Xselection in [100,101,102,103,104,105,106]:#[50]:#[86]:#,80,81,82,83,84,85]:
#for Xselection in [106]:#[110,111,112,113,114,115,116]:
    
    #idataset == 50:
    frames = [143000,145000,147000,149000]
    path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'#top200k/'
    #idataset == 54:
    #frames = [123000] #,125000,127000,129000] # 123000 was VERY different than the others, so now rerunning
    #path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
    phopfileprefix = 'trainingAve.0'
    phopfilesuffix = '.cfg.phop'
    nparticlestotal = 3873924 # 200000
    numparticlestoanalyze = 3873924 #  10000
    dset=''
    
    
    numheaderlines=9
    softfileprefix = path+'trainingAve.0'
    softfilesuffix = '.cfg'
    if Xselection == 0:
        print("Soft")
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
        softfilesuffix = ''
        numheaderlines=0
        binmid = 2
        binrange = 10
        column = 0
    elif Xselection == 1:
        print("100fs-averaged centrosymmetry values")
        dset='20'
        column = 5
        binmid = 8
        binrange = 16
    elif Xselection == 2:
        print("ave Energy")
        column = 6
        binmid = -3.2
        binrange = .6
    elif Xselection == 3:
        print("abs ( voro vol - volavg)")
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Vol_vs_soft_Al463K_t'
        softfilesuffix = '.txt'
        numheaderlines=1
        column = 0 # 0 # 1 (col 1 is softness, with which can verify ordering of file entries)
        binmid = -4 # 18 # 1
        binrange = 8 # 6 # 12
    elif Xselection == 4:
        print("centrosymmetry of 100-fs avg positions")
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Centrosymm_vs_soft_Al463K_t'
        softfilesuffix = '.txt'
        numheaderlines=1
        column = 0
        binmid = 8
        binrange = 16
    elif Xselection == 5:
        print("Soft, train GB only")
        softfileprefix = path+'KA_N3873924_tag40exc1.00win2sam10000slow1800_Softness.'
        softfilesuffix = ''
        numheaderlines=0
        binmid = 2
        binrange = 10
        column = 0
    
    
    # The following are like Xselection==1, but with different temperatures
    elif Xselection == 20:
        print("ave centrosymmetry")
        dset='20'
        print("ATTENTION THAT ONE OF THESE FRAMES, 143100 doesnt have aggresive filtering data")
        frames = [143000,145000,147000,149000]
        frames = [145000,147000,149000]
        column = 5
        binmid = 8
        binrange = 16
    elif Xselection == 21:
        print("ave centrosymmetry, 556K")
        dset='21'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 5
        binmid = 8
        binrange = 16
    elif Xselection == 22:
        print("ave centrosymmetry, 620K")
        dset='22'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 5
        binmid = 8
        binrange = 16
    elif Xselection == 23:
        print("ave centrosymmetry, 694K")
        dset='23'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 5
        binmid = 8
        binrange = 16
    elif Xselection == 24:
        print("ave centrosymmetry, 756K")
        dset='24'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 5
        binmid = 8
        binrange = 16
    elif Xselection == 25:
        print("ave centrosymmetry, 834K")
        dset='25'
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 5
        binmid = 8
        binrange = 16
    
    
    # The following are like Xselection 2, but with diff. temperatures and more frames 
    elif Xselection == 30:
        print("ave Energy, T0")
        dset='30'
        print("ATTENTION THAT ONE OF THESE FRAMES, 143100 doesnt have aggresive filtering data")
        frames = [143100,145000,147000,149000]
        frames = [143000,145000,147000,149000]
        column = 6
        binmid = -3.2
        binrange = .6
    elif Xselection == 31:
        print("ave Energy, T1")
        dset='31'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 6
        binmid = -3.2
        binrange = .6
    elif Xselection == 32:
        print("ave Energy, T2")
        dset='32'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 6
        binmid = -3.2
        binrange = .6
    elif Xselection == 33:
        print("ave Energy, T3")
        dset='33'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 6
        binmid = -3.2
        binrange = .6
    elif Xselection == 34:
        print("ave Energy, T4")
        dset='34'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 6
        binmid = -3.2
        binrange = .6
    elif Xselection == 35:
        print("ave Energy, T5")
        dset='35'
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        column = 6
        binmid = -3.2
        binrange = .6


    # These 50-55 are from plotSoftVsProbHop, but now are copied to here and given X labels, to unify and make running my scripts simpler
    # 6 Temperatures, Train on all #
    elif Xselection == 56:
        print("Softness T310K, train on all")
	dset='56'
        frames = [143000,145000,147000,149000]
        frames = [x for x in range(141000,150000+1,1000)]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin310K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 50:
        print("Softness T463K, train on all")
        dset='50'
        print("ATTENTION THAT ONE OF THESE FRAMES, 143100 doesnt have aggresive filtering data")
        frames = [143000,145000,147000,149000]
        #frames = [145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 51:
        print("Softness T543K, train on all")
        dset='51'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 52:
        print("Softness T620K, train on all")
	dset='52'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 53:
        print("Softness T694K, train on all")
	dset='53'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 54:
        print("Softness T756K, train on all")
	dset='54'
        frames = [123000,125000,127000,129000] # 123000 was VERY different than the others, so now reran
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 55:
        print("Softness T834K, train on all")
	dset='55'
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0



    # 6 Temperatures, Train on T620K, interior of GB #
    elif Xselection == 96:
        print("Softness T310K, train on interior GB")
	dset='96'
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin310K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag92exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 90:
        print("Softness T463K, train on interior GB")
        dset='90'
        print("ATTENTION THAT ONE OF THESE FRAMES, 143100 doesnt have aggresive filtering data")
        frames = [143000,145000,147000,149000]
        #frames = [141000]#,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag92exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 91:
        print("Softness T543K, train on interior GB")
        dset='91'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag92exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 92:
        print("Softness T620K, train on interior GB")
	dset='92'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag92exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 93:
        print("Softness T694K, train on interior GB")
	dset='93'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag92exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 94:
        print("Softness T756K, train on interior GB")
	dset='94'
        frames = [123000,125000,127000,129000] # 123000 was VERY different than the others, so now reran
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag92exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0
    elif Xselection == 95:
        print("Softness T834K, train on interior GB")
	dset='95'
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
        softfileprefix = path+'KA_N3873924_tag92exc1.00win2sam10000slow1800_Softness.'
	softfilesuffix=''
	numheaderlines=0
	binmid = 2
	binrange = 10
	column = 0




    # The following are like Xselection 4, but with diff. temperatures and more frames 
    elif Xselection == 86:
        print("centrosymmetry of 100-fs avg positions, T6")
        dset='86'
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin310K/'#top200k/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset56/Centrosymm_vs_soft_Al310K_t'
        softfilesuffix = '_seeCNAall.txt'
        numheaderlines=1
        column = 0
        binmid = 8
        binrange = 16
    elif Xselection == 80:
        print("centrosymmetry of 100-fs avg positions, T0")
        dset='80'
        print("ATTENTION THAT ONE OF THESE FRAMES, 143100 doesnt have aggresive filtering data")
        frames = [143000,145000,147000,149000]
        frames = [145000,147000,149000]
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Centrosymm_vs_soft_Al463K_t'
        softfilesuffix = '_seeCNAall.txt'
        numheaderlines=1
        column = 0
        binmid = 8
        binrange = 16
    elif Xselection == 81:
        print("centrosymmetry of 100-fs avg positions, T1")
        dset='81'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset51/Centrosymm_vs_soft_Al556K_t'
        softfilesuffix = '_seeCNAall.txt'
        numheaderlines=1
        column = 0
        binmid = 8
        binrange = 16
    elif Xselection == 82:
        print("centrosymmetry of 100-fs avg positions, T2")
        dset='82'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset52/Centrosymm_vs_soft_Al620K_t'
        softfilesuffix = '_seeCNAall.txt'
        numheaderlines=1
        column = 0
        binmid = 8
        binrange = 16
    elif Xselection == 83:
        print("centrosymmetry of 100-fs avg positions, T3")
        dset='83'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset53/Centrosymm_vs_soft_Al694K_t'
        softfilesuffix = '_seeCNAall.txt'
        numheaderlines=1
        column = 0
        binmid = 8
        binrange = 16
    elif Xselection == 84:
        print("centrosymmetry of 100-fs avg positions, T4")
        dset='84'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset54/Centrosymm_vs_soft_Al756K_t'
        softfilesuffix = '_seeCNAall.txt'
        numheaderlines=1
        column = 0
        binmid = 8
        binrange = 16
    elif Xselection == 85:
        print("centrosymmetry of 100-fs avg positions, T5")
        dset='85'
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset55/Centrosymm_vs_soft_Al834K_t'
        softfilesuffix = '_seeCNAall.txt'
        numheaderlines=1
        column = 0
        binmid = 8
        binrange = 16

    # Like Xselection == 3, Voro vol, but with multiple temperatures
    elif Xselection == 106:
        print("Voro vol of 100-fs avg positions, T6, abs ( voro vol - volavg)")
        dset='106'
        frames = [149000]#143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin310K/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Vol_vs_soft_Al310K_t'
        softfilesuffix = '.txt'
        numheaderlines=1
        column = 0 
        binmid = 18
        binrange = 12
    elif Xselection == 100:
        print("Voro vol of 100-fs avg positions, T0, abs ( voro vol - volavg)")
        dset='100'
        frames = [149000]#143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Vol_vs_soft_Al463K_t'
        softfilesuffix = '.txt'
        numheaderlines=1
        column = 0 
        binmid = 18
        binrange = 12
    elif Xselection == 101:
        print("Voro vol of 100-fs avg positions, T1, abs ( voro vol - volavg)")
        dset='101'
        frames = [139000]#143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Vol_vs_soft_Al556K_t'
        softfilesuffix = '.txt'
        numheaderlines=1
        column = 0 
        binrange = 12
    elif Xselection == 102:
        print("Voro vol of 100-fs avg positions, T2, abs ( voro vol - volavg)")
        dset='102'
        frames = [139000]#143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Vol_vs_soft_Al620K_t'
        softfilesuffix = '.txt'
        numheaderlines=1
        column = 0 
        binrange = 12
    elif Xselection == 103:
        print("Voro vol of 100-fs avg positions, T3, abs ( voro vol - volavg)")
        dset='103'
        frames = [129000]#143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Vol_vs_soft_Al694K_t'
        softfilesuffix = '.txt'
        numheaderlines=1
        column = 0 
        binrange = 12
    elif Xselection == 104:
        print("Voro vol of 100-fs avg positions, T4, abs ( voro vol - volavg)")
        dset='104'
        frames = [129000]#143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Vol_vs_soft_Al756K_t'
        softfilesuffix = '.txt'
        numheaderlines=1
        column = 0 
        binrange = 12
    elif Xselection == 105:
        print("Voro vol of 100-fs avg positions, T5, abs ( voro vol - volavg)")
        dset='105'
        frames = [149000]#143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'
        softfileprefix = '/home/tsharp/projsoft/softfind/tasanalysis/datahistsoft/dset50/Vol_vs_soft_Al834K_t'
        softfilesuffix = '.txt'
        numheaderlines=1
        column = 0 
        binrange = 12
        print("Worrysome that dset folder does not match dset??")


    # The following are like Xselection 2, but with diff. temperatures and more frames 
    elif Xselection == 110:
        print("Energy of avg, T0")
        dset='110'
        print("ATTENTION THAT ONE OF THESE FRAMES, 143100 doesnt have aggresive filtering data")
        frames = [143100,145000,147000,149000]
        frames = [143000,145000,147000,149000]
        frames = [143000,145000,147000,149000]
        softfilesuffix = '.cfg.idxyzPE'
        column = 4
        binmid = -3.2
        binrange = .6
    elif Xselection == 111:
        print("Energy of avg, T1")
        dset='111'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        softfilesuffix = '.cfg.idxyzPE'
        column = 4
        binmid = -3.2
        binrange = .6
    elif Xselection == 112:
        print("Energy of avg, T2")
        dset='112'
        frames = [133000,135000,137000,139000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        softfilesuffix = '.cfg.idxyzPE'
        column = 4
        binmid = -3.2
        binrange = .6
    elif Xselection == 113:
        print("Energy of avg, T3")
        dset='113'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        softfilesuffix = '.cfg.idxyzPE'
        column = 4
        binmid = -3.2
        binrange = .6
    elif Xselection == 114:
        print("Energy of avg, T4")
        dset='114'
        frames = [123000,125000,127000,129000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        softfilesuffix = '.cfg.idxyzPE'
        column = 4
        binmid = -3.2
        binrange = .6
    elif Xselection == 115:
        print("Energy of avg, T5")
        dset='115'
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        softfilesuffix = '.cfg.idxyzPE'
        column = 4
        binmid = -3.2
        binrange = .6
    elif Xselection == 116:
        print("Energy of avg, T6")
        dset='116'
        frames = [143000,145000,147000,149000]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin310K/'#top200k/'
        softfileprefix = path+'trainingAve.0'
        softfilesuffix = '.cfg.idxyzPE'
        column = 4
        binmid = -3.2
        binrange = .6

    else : 
        print("No Xselection =", Xselection)
        import sys
        sys.exit(0)
    
    

    ####################################################
    # Input parameters
    ####################################################

    doshow=1 # 0 or 1
    regenerate=0 # 0 or 1
    
    # Should only consider subset of particles?  -1 (consider all)  0 (disordered, CNA!=FCC)  1 (CNA=FCC)
    consideronlytype = -1 # 10 # -1 # considering all atoms would give higher Q's but P_R(mu-sigma) = undefined
    # 10 (very disordered) 11 (3.4A near CNA=FCC/HCP) # Super strict filter # Use the aggresive near-FCC filtering  

    # Define the bins [-4.0 +/- 0.25], [-3.5 +/- 0.25], ... [4.0 +/- 0.25]
    nbins = 41
    binwidth = 1.0*binrange/(nbins-1) # bin range is from lo-bin-center to hi-bin-center
    nbinedges = nbins + 1
    binedgemin = binmid - binrange/2.0 - binwidth/2.0 
    bin_edges = binedgemin + binwidth*np.arange(nbinedges)
    bin_centers = (bin_edges[0:nbins] + bin_edges[1:1+nbins]) / 2.0
    
    # Time (Number of frames * file spacign) between softness and last time we look for a hop
    delaywindow = 0 # fsec (0 is usual.  200 was usual) units is file ids, in this case femtoseconds
    framespacing = 100
    
    # Plotting setup and colors
    #axs = plt.subplot(1,1,1)#, sharex=True)
    fig1 = plt.figure(1,figsize=(6,6))
    axs = fig1.add_subplot(111)
    hotclr=['violet','purple','b','green','orange','r','k']#'teal','gold'
    hotclr+=hotclr
    print(bin_edges,"bin edges")
    print(bin_centers,"bin centers")
    #import sys
    #sys.exit(0)
    
    if regenerate == 1:
        ####################################################
        # Read softness files into softness_by_frame_and_id list
        ####################################################
        soft_vs_fr_id = []
        for frame in frames:
            print("Opening softness file", softfileprefix+str(frame)+softfilesuffix)
            softness = getjth_dump_col(softfileprefix+str(frame)+softfilesuffix, numparticlestoanalyze, column, numheaderlines)
            
            if Xselection == 3:
                print("Taking the abs of the difference of the input X field from its average")
                softness =  np.abs( softness - np.mean(softness))
        
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
            
            print("Checking if particles hop")
            for delayframe in frame + framespacing*np.arange(1 + delaywindow/framespacing):
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
                    label=path[-13::], color=hotclr[np.mod(Xselection,7)])
    print(bin_centers)
    print(fractionhopping_vs_bin)
    print(numberineach_bin)
    #plt.semilogy(...)






savelabel='_20180420' # '_20170601E' # '_20170503'
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
