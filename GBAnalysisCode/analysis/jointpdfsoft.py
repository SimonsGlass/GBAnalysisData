
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


# Figure settings
##fig1 = plt.figure(1,figsize=(15, 15))
fig1 = plt.figure(figsize=(10, 10))
fig1.subplots_adjust(wspace=1)
#ax = fig1.add_subplot(111)
#figjointpdf = plt.figure(3,figsize=(15,25)) # (9, 15)) #+datasetnum
colors=['k','r','g','b','y','pink','purple']
colors=colors+colors
savepdf=1

#nsoftbins=24
nsoftbins=36
minsoft = -4
maxsoft = 8
nareabins=24

datasets = [1,2,3]
numdatasets = len(datasets)
for datasetnum, idataset in enumerate(datasets):


  pathsoft = 'jointPDF_Soft_and_X/'
  patharea = 'jointPDF_Soft_and_X/'
  pathcna='/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'
  filenamecna='CNAtrainingAve.0143000.cfg'
  numcnaheaderlines=9
  pathout = 'jointPDF_Soft_and_X/'
  numheaderlines=1
  minarea = 0
  maxarea = 18
  #maxarea = 14
  areacol=0
  softcol=1
  if idataset == 0:
    Natoms = 200
    filenamesoft = 'Centrosymm_vs_soft_Al463K_t143000slice.txt'
    filenamearea = 'Centrosymm_vs_soft_Al463K_t143000slice.txt'
  if idataset == 1:
    Natoms = 2000000
    filenamesoft = 'Centrosymm_vs_soft_Al463K_t143000.txt'#slice.txt
    filenamearea = 'Centrosymm_vs_soft_Al463K_t143000.txt'#slice.txt
  if idataset == 2:
    Natoms = 2000000
    filenamesoft = 'Vol_vs_soft_Al463K_t143000.txt'
    filenamearea = 'Vol_vs_soft_Al463K_t143000.txt'
    #minarea = 15
    minarea = 15
    #maxarea = 22
    maxarea = 22
  if idataset == 3:
    Natoms = 2000000
    numheaderlines=9
    pathsoft='/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'
    patharea='/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/'
    filenamesoft='trainingAve.0143000.cfg.phop.CNA0soft'
    filenamearea='trainingAve.0143000.cfg'
    areacol=6
    softcol=5
    #minarea = -3.45
    minarea = -3.45
    #maxarea = -2.95
    maxarea = -2.95
  
  softedges=1.0*np.arange(nsoftbins+1)/nsoftbins*(maxsoft - minsoft) + minsoft
  areaedges=1.0*np.arange(nareabins+1)/nareabins*(maxarea - minarea) + minarea
  X, Y = np.meshgrid(softedges, areaedges)
  
  area = getjth_dump_col(patharea+filenamearea, Natoms, areacol, numheaderlines)
  soft = getjth_dump_col(pathsoft+filenamesoft, Natoms, softcol, numheaderlines)
  cna = getjth_dump_col(pathcna+filenamecna, Natoms, 0, numcnaheaderlines)
  
  if 0:
      filteredarea = []
      filteredsoft = []
      for i in np.arange(Natoms):
        #if soft[i] > -0.25 and soft[i] < 0.25:
        #if soft[i] > 0.5:
        #if soft[i] > -1.5:
        if cna[i] == 0:
          filteredarea += [area[i]]
          filteredsoft += [soft[i]]
      area = np.array(filteredarea)
      soft = np.array(filteredsoft)
      print(np.mean(area))


  hist2d, xbin_edges, ybin_edges = np.histogram2d(area, soft, normed=True, bins=(areaedges, softedges))
  # fig = plt.figure(figsize=(8, 8))    
  # save local variables
  #np.savez(path+'jointhistsoft'+thick+anotparam+'_minarea'+str(minareacutoff)+'.npz', nsoftbins=nsoftbins, nareabins=nareabins, hist2d=hist2d, ae=ae, se=se )
  
  
  #plt.subplot(1,4,1)
  ax = fig1.add_subplot(2,2,datasetnum+1)#1,1,1)#4,2,datasetnum+1)
  #ax.tick_params(direction='in', length=2, width=2, colors='k')
  ax.set_title('Joint PDF p03D=')#+str(myp03Dlabel))
  #ax.set_title(thick+'Min area cutoff='+str(myminarealabel))
  if datasetnum == 0:
      ax.set_ylabel('centrosymm ')            
  if datasetnum == numdatasets-1:
      ax.set_xlabel('soft')
  noisefloor=np.min(hist2d[np.where(hist2d > 0)])
  print(noisefloor, "is noise floor")
  here=(hist2d == 0)
  hist2d[here] = noisefloor/1000.
  plt.grid(b=True,which='major')
  mappable = ax.pcolormesh(X, Y, np.log10(hist2d), cmap=plt.get_cmap('Blues'), vmin=-5.0, vmax=0.0, grid=True)#, vmax=10**(-.5)) # hist2d)    
  #ax.set_aspect(.5)
  fig1.colorbar(mappable)
  #plt.axis([minsoft,maxsoft,minarea,maxarea])
  #areabinwidth=1.0*(maxarea-minarea)/nareabins
  plt.axis([-4,5,minarea,maxarea])

  plt.grid(b=True,which='major')
if (savepdf == 1):
        plt.savefig(pathout+"jointPDFs_soft_CSPEVol.eps")
plt.show()
#plt.savefig(pdfpage, format='pdf')

