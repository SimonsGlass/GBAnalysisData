# paste into calcdistmoved.py
# Construct a file of the distance moved of each atom over the simulation time

import numpy as np

natoms=3873924
frstfile = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/trainingAve.0140400.cfg'
lastfile = '/archive/p_ajliu/tsharp/projsoft/fromspencer/aluminum/poly746ggMishin463K/trainingAve.0150300.cfg'

x=np.zeros(natoms)
y=np.zeros(natoms)
z=np.zeros(natoms)
distance=np.zeros(natoms)

file = open(frstfile,'r')
for iline, line in enumerate(file):
    if iline < 9:
        print(iline, line)
    else:
        col = line.split(" ")
        iatom = int(col[0])-1
        x[iatom] = col[2]
        y[iatom] = col[3]
        z[iatom] = col[4]
file.close()

file = open(lastfile,'r')
for iline, line in enumerate(file):
    if iline < 9:
        print(iline, line)
    else:
        col = line.split(" ")
        iatom = int(col[0])-1
        distance[iatom] = np.sqrt(
                        (np.float64(col[2]) - x[iatom])**(2.0) + 
                        (np.float64(col[3]) - y[iatom])**(2.0) + 
                        (np.float64(col[4]) - z[iatom])**(2.0))
file.close()

file = open("movedistances.txt","w")
for iatom in range(natoms):
    file.write(str(distance[iatom])+"\n")
file.close()

