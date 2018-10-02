#paste into hophistograms
# make histograms of atoms that hops N1 or N2 times
import numpy as np
import matplotlib.pyplot as plt
natoms=3873924



distances = np.genfromtxt("movedistances.txt")
nhops_iatom = np.genfromtxt("num_of_0.5hops.txt")
#nhops_iatom = np.genfromtxt("num_of_1.0hops.txt")


for inumhop,inumhopmax in [(0,0),(1,1),(2,2),(4,4),(8,9),(14,40)]:
    print(inumhop, inumhopmax)
    truesandfalses = np.bitwise_and((nhops_iatom >= inumhop), (nhops_iatom <= inumhopmax))
    here = np.where(np.bitwise_and((truesandfalses), (distances < 200)) )

    print(inumhop)
    print(here[0])
    print(len(here))
    print("")
    
    # Define the bins [-4.0 +/- 0.25], [-3.5 +/- 0.25], ... [4.0 +/- 0.25]
    nbins = 41 # 101 #401 # 41
    nbinedges = nbins + 1
    binwidth = 0.025 # 0.05 # 0.2 #.05 #.5
    binedgemin = 0.0 
    #bin_edges = binedgemin + binwidth*np.arange(nbinedges)
    bin_edges = binedgemin + binwidth*np.arange(nbinedges)**(1.5)
    print(bin_edges)
    nbins = len(bin_edges)-1

    hist, bin_edges = np.histogram(distances[here], normed=True, bins=bin_edges)
    bin_centers = (bin_edges[0:nbins] + bin_edges[1:1+nbins]) / 2.0
    fig1 = plt.figure(1,figsize=(7,7))
    axs = fig1.add_subplot(111)
    #axs.errorbar(bin_centers, np.log(hist), label=str(inumhop))
    axs.errorbar(bin_centers, (hist), label=str(inumhop))

    print(inumhop, " has avg ", (distances[here]).mean())
    
    heredistances = distances[here]
    here2,  = np.where(heredistances < .4)
    print(len(here2), len(heredistances))
    print(inumhop, " has the fraction of atoms with distances less than .4 is  ", len(here2) / 1.0 / len(heredistances))

#plt.axis([0,4,-8,1.5])
plt.legend(title='Nhops ',loc=4,prop={'size':12},bbox_to_anchor=(0.95,0.3))
plt.xlabel('D, Distance moved')
plt.ylabel('log P(D), Distribution of distances')
plt.grid(b=True, which='major')
plt.savefig('out20180427_del.pdf' )
plt.show()

