import numpy as np
import matplotlib.pyplot as plt

if 0:
    datasets=range(50,56+1-1)
    labels = [463,556,620,694,756,834,310]
    hotclr=['purple','b','green','orange','r','k','violet']#'teal','gold'
elif 0:
    datasets=[50,52,54,55]
    labels = [463,620,756,834]
    hotclr=['purple','green','r','k']#'teal','gold'
elif 1:
    datasets=range(50,56+1-1)
    labels = [463,556,620,694,756,834]
    hotclr=['purple','b','green','orange','r','k']#'teal','gold'

#datasets.reverse()
#hotclr.reverse()

linesty=['-','-','--','--']#'--']    
fig1 = plt.figure(1,figsize=(5,5))
axs = fig1.add_subplot(111)
#axs.set_yscale('log')

for i,dataset in enumerate(datasets):
        datafolder = "datahistsoft/dset"+str(dataset)

        # Open combined 
        subset = "CNAall"
        filename = datafolder+"/PDFsoft_trainallistag50_see"+subset+".txt"
        print("Opening", filename)
        datatable = np.genfromtxt(filename,skip_header=1)
        bins = datatable[:,0]
        hist = datatable[:,1]
        # Assume equally spaced bins
        total = np.sum(hist)*(bins[1]-bins[0])
        print (total, " is area of distrib and the num atoms is ", np.sum(hist))

        #plotsubsets = ["CNA0","CNAordered","VoroTop0","VoroTop1"]
        plotsubsets = ["VoroTop0","VoroTop1"]
        #plotsubsets = ["CNA0","CNAordered"]
        for j,subset in enumerate(plotsubsets):
                filename = datafolder+"/PDFsoft_trainallistag50_see"+subset+".txt"
                print("Opening", filename)
                datatable = np.genfromtxt(filename,skip_header=1)
                bins = datatable[:,0]
                hist = datatable[:,1]
                if 0:
                    hist /= total
                else:
                    # Assume equally spaced bins
                    print("Now I normalize each curve to area 1.0")
                    totalsubset = np.sum(hist)*(bins[1]-bins[0])
                    print("this one sums to ", totalsubset, " and num atoms is ", np.sum(hist)) 
                    hist /= totalsubset
                    print("")

                if j == 0:
                    plt.plot(bins,hist,color=hotclr[i],linestyle=linesty[j],marker='',lw=2,label=str(labels[i]))
                else:
                    plt.plot(bins,hist,color=hotclr[i],linestyle=linesty[j],marker='',lw=2)




#axs.set_yscale('log')
if 1:
    plt.axis([-3,5,0,1.3])
else:
    plt.axis([-3,6,0,1.3])
plt.legend(title='Temperature (K)',loc=4,prop={'size':12},bbox_to_anchor=(0.85,0.7))
plt.xlabel(r'S')
plt.ylabel(r'P(S)')
plt.grid(b=True, which='major')

outlabel='20170627'
if 1:
        plt.savefig(datafolder+'/_P_of_S_'+outlabel+'.pdf')
plt.show()
