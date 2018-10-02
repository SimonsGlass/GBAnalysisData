# assume pid starts from 0
import numpy as np

numberofatoms = 3873924

outfile = open("num_of_0.5hops.txt", 'w')
hopping_pids = np.genfromtxt("sorted_filtered0.5_ids_of_hops.txt")
nlinesout = 0
thecount = 0
for pid in hopping_pids:
    while nlinesout < pid:
        outfile.write(str(thecount)+"\n")
        nlinesout += 1
        thecount = 0
    thecount += 1

outfile.write(str(thecount)+"\n")
nlinesout += 1

while nlinesout < numberofatoms:
    outfile.write("0\n")
    nlinesout += 1

outfile.close()

