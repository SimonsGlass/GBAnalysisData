
'''
Tristan Sharp
12/1/16
Try to add CNA based on the input positions

To run: ovitos scriptToGetGBs filetoanalysze
/home/tsharp/bin/ovito/ovito-2.8.0-x86_64/bin/ovitos

12/7/2016 I verified that the output from ovitos matches the adaptive CNA from my Windows ovito

2018-07-04 I am going to analyze individual grain boundaries?
'''

import sys
from ovito.io import *
from ovito.modifiers import *
import numpy as np

#if len(sys.argv) != 2:
#        print("Usage: ovitos scriptToGetGBs filetoanalysze")
#        sys.exit()

#infile = sys.argv[1]
#infile = '694K/trainingAve.0123000.cfg.phop.CNA0soft'

case = int(sys.argv[1]) # 5
if case == 0:
    # This is made from "insertGBIDsintodump.cpp"
    infile = "620K/small.0130600.cfg.gaveids.phop.CNA0soft"
elif case == 1:
    infile = "310K/trainingAve.0143000.cfg.phop.CNA0soft"
elif case == 2:
    infile = "463K/trainingAve.0143000.cfg.phop.CNA0soft"
elif case == 3:
    infile = "556K/trainingAve.0133000.cfg.phop.CNA0soft"
elif case == 4:
    infile = "620K/trainingAve.0133000.cfg.phop.CNA0soft"
elif case == 5:
    infile = "694K/trainingAve.0123000.cfg.phop.CNA0soft"
elif case == 6:
    infile = "756K/trainingAve.0123000.cfg.phop.CNA0soft"
elif case == 7:
    infile = "834K/trainingAve.0143000.cfg.phop.CNA0soft"
elif case == 50:
    infile = "694K/trainingAve.0129000.cfg.phop.CNA0soft"
print("case = ", case)

basenameprefix = infile[infile.find("/")+1:infile.find("/")+6]
print(basenameprefix, "is prefix")
outfilegrainids = infile.replace(basenameprefix, "grainIDs_"+basenameprefix)
nndist = 2.95 # in the length units of the file (Angstrom)
print("reading ", infile)
node = import_file(infile)

####################
# Compute Grain IDS
print("Compute Grain IDS")
node.modifiers.append(WrapPeriodicImagesModifier())
node.modifiers.append(CommonNeighborAnalysisModifier())

modifier = SelectParticleTypeModifier(property = "Structure Type")
modifier.types = {0,2,3,4} # All but FCC=1 (0-disorder, 2 HCP dont include, 3...)
node.modifiers.append(modifier)
node.modifiers.append(ExpandSelectionModifier(cutoff=nndist*1.05))
node.modifiers.append(InvertSelectionModifier())
# Assigns numbers starting from 1, 0 is for atoms that get no cluster
node.modifiers.append(ClusterAnalysisModifier(cutoff=nndist*1.1, sort_by_size=True, only_selected=True))
node.compute()
#grainID_vs_atomindex = node.output.particle_properties['Cluster'].array

##################
# Save intermediate output (Grain IDs) for viewing 
print("Write Grain IDS into LAMMPS file")
export_file(node, outfilegrainids, "lammps_dump",
            columns = ["Particle Identifier","Position.X", "Position.Y", "Position.Z", "phop", "soft", "Cluster", "Structure Type"])


#################
# Find atoms that neighbor two (and only two) grains
node = import_file(outfilegrainids)
node.compute()
#node.modifiers.append(WrapPeriodicImagesModifier())
grainID_vs_atomindex = node.output.particle_properties['cluster'].array
ngrains = int(np.max(grainID_vs_atomindex))
print("ngrains", ngrains)
natoms = node.output.number_of_particles
print("natoms",natoms)
GBsize_vs_GBID = np.zeros(ngrains*ngrains)
GBID_vs_atomindex = [-1 for iten in range(natoms)]
closeby = nndist*3.5 # 4.5 # 3.5
print("Find atoms within ", closeby, " of two and only two grains")

iGB = 0 
# Compute Grain IDS
for igrain in range(min(ngrains,50)): # skip grain 0 which was all non-clustered atoms

    print("igrain = ", igrain)
    for jgrain in range(igrain+1,min(ngrains,50)):
        print("jgrain = ", jgrain)

        # Get StillGood array - from 
        #  node2 = nodeempty
        #  node2 = copy.deepcopy(node)
        #node = import_file(outfilegrainids)
        
        node.modifiers.append(SelectExpressionModifier(expression = 'cluster == '+str(igrain+1)))
        node.modifiers.append(ExpandSelectionModifier(cutoff=closeby))
        node.modifiers.append(ComputePropertyModifier(output_property="NearGrainI",expressions=["Selection"]))
        
        node.modifiers.append(SelectExpressionModifier(expression = 'cluster == '+str(jgrain+1)))
        node.modifiers.append(ExpandSelectionModifier(cutoff=closeby))
        node.modifiers.append(ComputePropertyModifier(output_property="NearGrainJ",expressions=["Selection"]))

        node.modifiers.append(ComputePropertyModifier(output_property="StillGood",expressions=["(NearGrainI && NearGrainJ)"]))
        
        node.compute()
        StillGood = node.output.particle_properties['StillGood'].array
        for imod in range(len(node.modifiers)):
            node.modifiers.pop()

        if np.sum(StillGood) > 10:
            #  node2 = nodeempty
            #  node2 = import_file(infile)
            #node = import_file(outfilegrainids)
            #  node2 = copy.deepcopy(node)
            for kgrain in ['all']:
                    #for kgrain in range(150 +0*ngrains):
                    #if ((kgrain != igrain) and (kgrain != jgrain)):
                    #print(kgrain, " kgrain")
                    myexpr = '(cluster != 0) && (cluster != '+str(igrain+1)+') && (cluster != '+str(jgrain+1)+')'
                    node.modifiers.append(SelectExpressionModifier(expression = myexpr))
                    node.modifiers.append(ExpandSelectionModifier(cutoff=closeby))
                    node.modifiers.append(ComputePropertyModifier(output_property="NearGrainK",expressions=["Selection"]))
                    #print(node.modifiers)
                    print("Overlap GBs")
                    node.compute()
                    NearGrainK = node.output.particle_properties['NearGrainK'].array
                    StillGood = np.all([StillGood, 1.0-NearGrainK],0) # element-wise-and
                    for imod in range(len(node.modifiers)):
                        node.modifiers.pop()
            
            #GBsize_vs_GBID[igrain*ngrains+jgrain] = np.sum(StillGood)
            GBsize_vs_GBID[iGB] = np.sum(StillGood)
            for iatom in range(natoms):
                    if StillGood[iatom]:
                        GBID_vs_atomindex[iatom] = iGB # igrain*ngrains+jgrain
            iGB += 1

atomsgrainIDsfilename = infile.replace(basenameprefix,'grainindex_vs_atomindex_from_'+basenameprefix)
atomsgrainIDsfile = open(atomsgrainIDsfilename,'w')
for item in GBID_vs_atomindex:
        atomsgrainIDsfile.write("%s\n" % str(item))
atomsgrainIDsfile.close()
print("Wrote "+atomsgrainIDsfilename)
nGBs = len(GBsize_vs_GBID)
print("nGBs",nGBs)

if 0:
    import sys
    sys.exit()
if 1:
    import subprocess
    # Compiled "insertGBIDsintodump.cpp"
    subprocess.call(["./a.out", infile, atomsgrainIDsfilename])

