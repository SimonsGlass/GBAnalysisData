
# TSharp
# 2016.12.19

# Script 
# to read the Softness of each particle at a particular frame.

import numpy as np
import matplotlib.pyplot as plt

regenerate = 0
extralabel = 'delete' # ''

# Plotting setup
fig1 = plt.figure(1,figsize=(5,5))
axs = fig1.add_subplot(111)
#axs = plt.subplot(1,1,1)#, sharex=True)
hotclr=['purple','b','green','orange','r','k']#'teal','gold'

# Define the bins [-4.0 +/- 0.25], [-3.5 +/- 0.25], ... [4.0 +/- 0.25]
nbins = 41 #401 # 41
nbinedges = nbins + 1
binwidth = 0.5 #.05
binedgemin = -0.5*binwidth*(nbins) # centering the bins at zero
bin_edges = binedgemin + binwidth*np.arange(nbinedges)

# Select datasets 
datasets=[40,41,42,43,44,45]
for idataset in datasets:
    if idataset == 0:
            softfile = "delete.txt"
            path = 'test'
    elif idataset == 40:
            frames = [140600]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'#top200k/'
            softfile = path+'KA_N3873924_tag40exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 41:
            frames = [130600]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
            softfile = path+'KA_N3873924_tag41exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 42:
            framess = [131000]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
            softfile = path+'KA_N3873924_tag42exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 43:
            frames = [121000]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
            softfile = path+'KA_N3873924_tag43exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 44:
            frames = [121000]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
            softfile = path+'KA_N3873924_tag44exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 45:
            frames = [141000]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
            softfile = path+'KA_N3873924_tag45exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 50:
            frames = [140600]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'#top200k/'
            softfile = path+'KA_N3873924_tag50exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 51:
            frames = [130600]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin556K/'#top200k/'
            softfile = path+'KA_N3873924_tag51exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 52:
            framess = [131000]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin620K/'#top200k/'
            softfile = path+'KA_N3873924_tag52exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 53:
            frames = [121000]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin694K/'#top200k/'
            softfile = path+'KA_N3873924_tag53exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 54:
            frames = [121000]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin756K/'#top200k/'
            softfile = path+'KA_N3873924_tag54exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    elif idataset == 55:
            frames = [141000]
            path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin834K/'#top200k/'
            softfile = path+'KA_N3873924_tag55exc1.00win2sam10000slow1800_Softness.'+str(frames[0])
    else: 
            print("No idataset =", idataset)
            import sys
            sys.exit(0)


    ####################################################
    # Read the text file
    ####################################################
    if regenerate == 1:
        print("Opening softness file")
        with open(softfile) as f:
            softness = [map(float,num.split()) for num in f.readlines()]
    
        # Histogram the data
        hist, bin_edges = np.histogram(softness, normed=True, bins=bin_edges)
        bin_centers = (bin_edges[0:nbins] + bin_edges[1:1+nbins]) / 2.0
        print(np.sum(hist*np.diff(bin_edges)))


    ####################################################
    # Plot
    ####################################################
    

    outfilenamebase = 'outfiles/hist_dset'+str(idataset)
    if regenerate == 1:
        outputfile = open(outfilenamebase+'_bincenters.txt'+extralabel,'w')
        for item in bin_centers:
                outputfile.write("%s\n" % str(item))
    else:
        bin_centers = np.genfromtxt(outfilenamebase+'_bincenters.txt'+extralabel)
    if regenerate == 1:
        outputfile = open(outfilenamebase+'_pdf.txt'+extralabel,'w')
        for item in hist:
            outputfile.write("%s\n" % str(item))
    else:
        hist = np.genfromtxt(outfilenamebase+'_pdf.txt'+extralabel)

    print(bin_centers,"bin centers")
    axs.errorbar(bin_centers, hist,
                    label=path[-13::], color=hotclr[np.mod(idataset,6)])
    #plt.semilogy(...)


#plt.axis([0,len(p_hop[0]),0,1.2])
plt.axis([-5,5,1e-5,6e-1])
plt.legend(title='Temperature ',loc=4,prop={'size':12},bbox_to_anchor=(0.95,0.3))
plt.xlabel('S, Softness')
plt.ylabel('P(S), Distribution of softness')
#plt.xlabel('S, Softness')
#plt.ylabel('P_R(S), Rearrange probability, within window')
plt.grid(b=True, which='major')
plt.savefig(outfilenamebase+'_pdf.'+extralabel+'.pdf')
plt.show()


