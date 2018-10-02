
# TSharp
# 2017.05.09

# Python Script 
# Plot mean softness vs z-coordinate.
# After reading in a dump file of positions and softness

import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit



########################################################################
# Function to read jth col from dump file. 
# Indexing stats at 0, but 1st column white space counts.
# Here used to read phop, in 5th col. (id type x y z phop)
########################################################################
def getjth_dump_col(filename, Nrecords, j, numheaderlines):

        print("Reading "+str(Nrecords)+" values in col "+str(j)+" from file ..." + filename[-50::])

	with open(filename,'r') as file:
                # Skip the header
                for iline in range(numheaderlines):
                        file.readline()
                # Read in the whole file, then crop
                # Split each remaining line by spaces, extract jth col, cast as float
                myfloatcol = [float((line.split())[j]) for line in file]
                myfloatcol = myfloatcol[0:Nrecords]

        if (len(myfloatcol) != Nrecords):
                print("Length read in from file, "+str(len(filename))+
                      " doesnt match Input parameter "+str(Nrecords))
                import sys
                sys.exit(0)
        #print("First elements read: ", myfloatcol[0:10])
        return myfloatcol


#http://stackoverflow.com/questions/11507028/fit-a-gaussian-function
def gauss(x, *p):
        A, mu, sigma, offset = p
        return offset + A*numpy.exp(-(x-mu)**2/(2.*sigma**2))
def find_max(bin_centers, meanSthisframe):
        #return bin_centers[np.argmax(meanSthisframe)]
        p0 = [3.5,0.0,1.0,-2.5]
        coeff, var_matrix = curve_fit(gauss, bin_centres, meanSthisframe, p0=p0)

        # Get the fitted curve
        hist_fit = gauss(bin_centres, *coeff)
        plt.plot(bin_centres, hist, label='Test data')
        plt.plot(bin_centres, hist_fit, label='Fitted data')
        # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
        print 'Fitted mean = ', coeff[1]
        print 'Fitted standard deviation = ', coeff[2]
        plt.show()


regenerate = 1
extralabel = 'delete' # ''

# Plotting setup
fig1 = plt.figure(1,figsize=(5,5))
axs = fig1.add_subplot(111)
hotclr=['purple','b','green','orange','r','k']
hotclr=['purple','b','green','orange','r','k','gray','teal','skyblue','gold','violet']

# Define the z-coordinate bins [-4.0 +/- 0.25], [-3.5 +/- 0.25], ... [4.0 +/- 0.25]
nbins = 101 # 100 #401 # 41
nbinedges = nbins + 1
binwidth = 0.25 # 0.5 #.05
binedgemin = -0.5*binwidth*(nbins) # centering the bins at zero
bin_edges = binedgemin + binwidth*np.arange(nbinedges)
bin_centers = (bin_edges[0:nbinedges-1] + bin_edges[1:nbinedges])/2.0
print(bin_centers)

datasets = [0,1]
#datasets = range(20,38,2)
for idataset in datasets:
    if idataset == 0:
        Natoms = 43988
        frames = range(80000,80100,300)#100)
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/biX300Kp01n1n34.n3n14avg/'
    if idataset == 1:
        Natoms = 43988
        frames = range(90000,90100,300)
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/biX300Kp01n1n34.n3n14avg/'
    if idataset == 10:
        Natoms = 43722
        frames = range(80000,81000,100)
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/biX300Kp0152n7.25n7avg/'

    # Check time dependence, especially during migration (time 80k+)
    if idataset >= 20 and idataset < 40:
        frames = [80000+10*100*(idataset-20)]
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/biX300Kp01n1n34.n3n14avg/'
        path = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/biX300Kp0152n7.25n7avg/'
        Natoms = 43988
        Natoms = 43722

    if regenerate == 1:
        meanS = np.zeros(nbins)
        meanSstddev = np.zeros(nbins)
        numAtZ = np.zeros(nbins)
        for iframe in frames:
            filename = 'dump.00'+str(iframe)+'.cfglammps.CNA0soft'

            meanSthisframe = np.zeros(nbins)
            numAtZthisframe = np.zeros(nbins)
            ####################################################
            # Read the text file
            ####################################################
            print("Opening softness file")
            zcol = 3
            softcol = 13
            numheaderlines = 9
            zpos = getjth_dump_col(path+filename, Natoms, zcol, numheaderlines)
            softness = getjth_dump_col(path+filename, Natoms, softcol, numheaderlines)
            
            
            # Find the mean softness at each z, to find the bin with the maximum
            for ielement, izpos in enumerate(zpos):
                zbin = round((izpos - np.min(bin_edges))/binwidth)
                if ((zbin > -1) and (zbin < nbins)):
                    numAtZthisframe[zbin] += 1
                    meanSthisframe[zbin] += softness[ielement]
                    



            for i in np.arange(nbins):
                if numAtZthisframe[i] != 0:
                    meanSthisframe[i] /= numAtZthisframe[i]

            # Repeat, but now offset by zofGB
            #zofGB = find_max(bin_centers, meanSthisfram) 
            #zofGB = bin_centers[np.argmax(meanSthisframe)]
            #zofGB = 0.
            bulkSoftness = np.sum(meanSthisframe[0:10])/10
            Sdistribution = (meanSthisframe - bulkSoftness)**2.0/np.sum((meanSthisframe - bulkSoftness)**2.0)/binwidth
            zofGB = np.sum(Sdistribution*bin_centers*binwidth)
            print("zofGB", zofGB, "since bulk softness is", bulkSoftness)
            zpos -= zofGB
            for ielement, izpos in enumerate(zpos):
                zbin = round((izpos - np.min(bin_edges))/binwidth)
                if ((zbin > -1) and (zbin < nbins)):
                    numAtZ[zbin] += 1
                    meanS[zbin] += softness[ielement]
                    meanSstddev[zbin] += softness[ielement] * softness[ielement]
        
        for i in np.arange(nbins):
            if numAtZ[i] != 0:
                meanS[i] /= numAtZ[i]
                meanSstddev[i] = np.sqrt(meanSstddev[i]/numAtZ[i]-meanS[i]*meanS[i])


    ####################################################
    # Plot
    ####################################################
    # May later do each x-stripe separately, since GB has ~0.5 Angstrom roughness    

    #axs.errorbar(bin_centers, numAtZ/10, marker='o',
    #        label=path[-13::], color=hotclr[np.mod(idataset,len(hotclr))])
    axs.errorbar(bin_centers, meanS, yerr=meanSstddev , marker='',#'o'
            label=path[-13::], color=hotclr[np.mod(idataset,len(hotclr))])
    #print(bin_centers, "plotted", meanS)


axs.errorbar(bin_centers, -2.5+3.5*np.exp(-bin_centers*bin_centers/8), 
                marker='', label='Gaussian', color='gray', linewidth=1)
                
outfilenamebase = 'spatialdistribsoft/del_zdep_'

#plt.axis([0,len(p_hop[0]),0,1.2])
#plt.axis([-10,10,-3,3])
plt.legend(title='Temperature ',loc=4,prop={'size':12},bbox_to_anchor=(0.95,0.7))
plt.xlabel('Z position, relative to most soft plane')
plt.ylabel('S_bar(z), Mean softness')
#plt.xlabel('S, Softness')
#plt.ylabel('P_R(S), Rearrange probability, within window')
plt.grid(b=True, which='major')
plt.savefig(outfilenamebase+'_pdf.'+extralabel+'.pdf')
plt.show()




