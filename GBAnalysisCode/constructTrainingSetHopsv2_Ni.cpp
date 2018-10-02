#include <iostream>
#include <sstream>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Resources/SimpleGrid.h"
#include "Resources/Settings.h"
#include "Minimization/minimizer.h"
#include "Database/Database.h"
#include "Resources/OutputLAMMPS.h"
#include "Resources/LoadD2Min.h"
#include "Resources/ComputeSymmetryFunctions.h"
#include "Resources/LoadSoftness.h"
#include "Resources/NetCDFInterface.h"
#include <set>



using namespace LiuJamming;
#define DIM 3
using namespace std;
typedef Eigen::Matrix<dbl, DIM, 1> dvec;
typedef Eigen::Matrix<dbl,DIM,DIM> dmat;
typedef MatrixInterface<dbl> MI;

//#define TASDEBUG

using namespace std;

//sort based on frame
bool Sort(const pair<pair<int,int>, int> &first, const pair<pair<int,int>, int> &second)
{
	return first.first.first < second.first.first;	
}

// this function, the sort call, requires c++11.  
vector<int> sortA_typesAndB_basedOnB(vector<int>& trainingParticles, vector<int>& trainingFrames) 
{

        //cout << "oldp" ; for (int i=0; i<trainingParticles.size(); i++) cout << " " <<  trainingParticles[i] ; cout << endl;
        //cout << "oldf" ; for (int i=0; i<trainingParticles.size(); i++) cout << " " <<  trainingFrames[i] ; cout << endl;
        struct comparer {
                // This struct takes two arguments, returning true if the first is biggest. 
                // Its constructor sets thingstocompare to the trainingFrames
                comparer(vector<int>& myframes) : thingstocompare(myframes){}
                bool operator () (const int& a, const int& b) const 
                        { return thingstocompare[a] < thingstocompare[b]; }
                vector<int>& thingstocompare;
        } mycomparer(trainingFrames);
        std::vector<int> neworderindices(trainingParticles.size());
        for (int i = 0; i < trainingParticles.size(); i++) neworderindices[i]=i;
        //cout << "old" ; for (int i=0; i<trainingParticles.size(); i++) cout << " " <<  neworderindices[i] ; cout << endl;
        std::sort(neworderindices.begin(), neworderindices.end(), mycomparer);
        // Note that I changed compile to -std=c++11 to allow this
        //cout << "new" ; for (int i=0; i<trainingParticles.size(); i++) cout << " " <<  neworderindices[i] ; cout << endl;

        vector<int> sortedTrainingType(trainingParticles.size());
        vector<int> sortedTrainingParticles(trainingParticles.size());
        vector<int> sortedTrainingFrames(trainingParticles.size());

        for (int i = 0; i < trainingParticles.size(); i++) {
                sortedTrainingType[i] = neworderindices[i] < (trainingParticles.size()/2)? 1 : 0;
                sortedTrainingParticles[i] = trainingParticles[neworderindices[i]];
                sortedTrainingFrames[i] = trainingFrames[neworderindices[i]];
        }
        //for (int i = 0; i < min(25,(int)(trainingParticles.size())); i++) cout << "sorted" << sortedTrainingFrames[i] << endl;
        trainingParticles = sortedTrainingParticles;
        trainingFrames = sortedTrainingFrames;
        //cout << "newp" ; for (int i=0; i<trainingParticles.size(); i++) cout << " " <<  trainingParticles[i] ; cout << endl;
        //cout << "newf" ; for (int i=0; i<trainingParticles.size(); i++) cout << " " <<  trainingFrames[i] ; cout << endl;
        return sortedTrainingType;
}


void ComputeDisplacement(int i, CStaticState<DIM> &StateI, int j, CStaticState<DIM> StateJ, dvec &Vector)
{
	dvec PosI, PosJ;
	StateI.GetParticlePositionVirtual(PosI,i);
	StateJ.GetParticlePositionVirtual(PosJ,j);

	StateI.GetBox()->MinimumDisplacement(PosI,PosJ,Vector);
	StateI.GetBox()->Transform(Vector);
}

	template <class T>
T max(vector<T> &list)
{
	T m = -100;
	for(int i = 0 ; i < list.size() ; i++)
		if( m < list[i])
			m = list[i];

	return m;
}


int ConstructTrainingSet(Settings &settings)
{

        ////////////////////
        // Import Settings
        ////////////////////
	string DataFilePath = settings["DataFilePath"];
	string InputFileSuffix = settings["InputFileSuffix"];
	string InputFilePrefix = settings["InputFilePrefix"];
	int Number = settings["Number"];
	string OutFileTag = settings["OutFileTag"];
	int Seed = settings["Seed"];
	int Samples = settings["Samples"];
        // Interval from which to draw training set.  Inclusive [StartFrame, EndFrame].
	int startFrame = settings["StartFrame"];
	int maxFrame = settings["EndFrame"];
        // Perform a check that training sets will only come from frames where phop was able to be defined
        vector<int> Frames = settings["Frames"];
	// Number of frames of window tR when calculating phop, and becomes file name of phop
	int HopWindowNFrames = settings["HopWindowNFrames"];
        int HopWindow = HopWindowNFrames*Frames[2]; // units is the file id's
	double HopCutoff = settings["HopCutoff"];
	//the timescale in which no rearrangement should happen if a particle is to be labelled slow moving
	int SlowTime = settings["SlowTime"];
	double excitationCutoff = settings["ExcitationCutoff"];
	int IncludeCNA = settings["IncludeCNA"]; // Include CNA as a symmetry function. 0 for no, 1 for yes

        ///////////////////
        // Check settings
        ///////////////////
        if (startFrame < Frames[0]+HopWindow) {
                printf("ERROR: First training set-accessible frame is earlier than phop could have been defined\n");
                return;
        }
        if (maxFrame > Frames[1]-HopWindow+Frames[2]) {
                printf("ERROR: Last training set-accessible frame is later than phop could have been defined\n");
                return;
        }


	/////////////////////////////
	//load the excitation data.
        /////////////////////////////
	char filename[256];
        sprintf(filename,"%s/keyschandlerKA_N%i_tag_AnalysisHops_w%i_c%.2f.nc",DataFilePath,Number,HopWindow,HopCutoff);
        cout << " Reading excitations from " << filename << endl; 

	vector<int> ExcitationStart;
	vector<int> ExcitationParticle;
	//vector<int> ExcitationParticleCNA;
	vector<double> ExcitationHop;

	LoadVector(filename,"Start",ExcitationStart);
	LoadVector(filename,"Particle",ExcitationParticle);
	//LoadVector(filename,"ParticleCNA",ExcitationParticleCNA);
	LoadVector(filename,"p_hop",ExcitationHop);

	//compute the frames / particles for the rearrangements
	vector<int> trainingFrames;
	vector<int> trainingParticles;
	//vector<int> trainingParticlesCNA;


#ifdef TASDEBUG
        cout << "ExcitationStart.size " << ExcitationStart.size() << endl;
#endif
        
	// First generate the list excitationList of all excitations 
        // that had peak phop big enough
	cout << "Identifying excitations.\n";
	list<int> excitationList;
        // TAS for each excitation
	for(int i = 0 ; i < ExcitationStart.size() ; i++) {
		if( ExcitationHop[i] >= excitationCutoff) // TAS if its phop value is large enough
			excitationList.push_back(i); // TAS Add this index of the excitation to a excitationList
#ifdef TASDEBUG
                if (i < 20) cout << "i" << 1*i << " phop " << 1*ExcitationHop[i] << " new cut" << 1*excitationCutoff << endl;
#endif
        }

        // Add to trainingFrames the frame number of the windowStart of the excitation
	cout << "Accumulate a subset of the " << excitationList.size() << " large-enough excitations. \n";
        if (Samples > excitationList.size()) {printf("Error: Not enough hops bigger than excitationCutoff to sample\n"); return 1;}
	MTRand random(Seed);
	for (int i = 0 ; i < Samples ; i++)
	{
                // We should prune instead of accumulate if we want to keep most of the excitations TAS
		int index = random.randInt(excitationList.size()-1);
		list<int>::iterator excitation = excitationList.begin();
		for(int it = 0 ; it < index ; it++)
			excitation++;

#ifdef TASDEBUG
                printf("We push excitation start frame %i and particle %i \n", ExcitationStart[*excitation], ExcitationParticle[*excitation]);
#endif                
		trainingFrames.push_back(ExcitationStart[*excitation]);  
		trainingParticles.push_back(ExcitationParticle[*excitation]);
		//trainingParticlesCNA.push_back(ExcitationParticleCNA[*excitation]);
		excitationList.erase(excitation);
	}



        cout << "Now identifying slowly moving particles, of which hard particles are considered a subset\n";
        // The code assumes that N atoms are in the dump file, with ids 1-N, not necessarily sorted
        // TAS an array that is Natoms rows where each row is the lits of excitations' windowStartFrames
        // TAS These sorts of lines would have trouble if IDs were not consecutive 1-N
        // Particle 1: [frame 1528120, frame 1529020]
        // Particle 2: [frame 1529000, frame 1529200, frame 1529700, frame 1529500]
        // Particle 3: []
        // Particle 4: [frame 128140]
	vector<list<int> > smallExcitations(Number);
	for(int i = 0; i < ExcitationStart.size() ; i++)
		smallExcitations[ExcitationParticle[i]].push_back(ExcitationStart[i]);

        // currentExcitation is an Natom-array of iterators i
        // that will go over smallExcitations
	vector<list<int>::iterator > currentExcitation(Number);
	for(int i =0 ; i < Number ; i++)
		currentExcitation[i] = smallExcitations[i].begin();


	// Select particles that do not undergo a "small" rearrangement over a timescale SlowTime
	cout << "There are a total of " << (maxFrame - startFrame + Frames[2])/SlowTime << " SlowTime Intervals in which to find slow particles.\n";
        // startFrame is the first center frame id_label (femtosec)
        // maxFrame ("End_Frame" in Settings file) is the last centerframe
        // TAS SlowTime should be the time in femtoseconds to say something is not re-arrange-able.
	// Go through the list of frames in increments of SlowTime.
        // Break up [startFrame, maxFrame] interval into chunks of length SlowTime. If a chunk does not have a rearrangement
        // Add i at time fr+SlowTime/2 as a particle to train on.  Maybe should ensure that did not finish rearranging
        if (maxFrame > Frames[1] - HopWindow + Frames[2]) {printf("ERROR: Cannot train on frames where phop undefined\n"); return 1;}
        if (SlowTime > Frames[1] - Frames[0] - 2*HopWindow + Frames[2]) {printf("ERROR: Not enough simulation time to recognize slow particles\n"); return 1;}
        if (SlowTime/2 < HopWindow){printf("ERROR: HopWindow implies that hops may still be occuring SlowTime/2 into the slow interval\n"); return 1;}


        vector<int> hardParticles;
        vector<int> hardFrames;
        // For each possibly quiet interval, then for each particle
        for(int fr = startFrame ; fr <= maxFrame - SlowTime; fr+=SlowTime)
	{
		cout << "Creating a list of slow particles to train on, considering frame " << fr << endl;
		for(int j = 0 ; j < Number ; j++)	
		{

                        // if frameStart is LARGER than the current centerframe and yet LESS
                        // than this center frame + SlowTime, then the particle is not slow over this interval

                        // this fast forwards currentExcitation until it's inside or beyond the SlowTime interval
			while(currentExcitation[j]!=smallExcitations[j].end() && *currentExcitation[j] < fr) 
			        currentExcitation[j]++;

			bool didRearrange = false;
                        // This particle did not rearrange in this interval if there are no future rearrangements of the particle.
                        // This particle rearranges if there is an excitation in this block of time
                        if (currentExcitation[j]!=smallExcitations[j].end() && *currentExcitation[j] < fr + SlowTime)
			{
				didRearrange = true; // Or at least initial stages of re-arrangement are within this SlowTime interval
                                // this fast forwards currentExcitation until it's outside of the SlowTime interval.
				while(currentExcitation[j]!=smallExcitations[j].end() && *currentExcitation[j] < fr+SlowTime) 
					currentExcitation[j]++;
			}

			if(!didRearrange)
			{
                                // TAS If no rearrangement, this notes the time in the middle of that SlowTime chunk
				hardFrames.push_back(fr + SlowTime/2);
				hardParticles.push_back(j); // TAS Add to hardParticles the particle index
			}	
		}
	}

	cout << "Identified " << hardParticles.size() << " slow moving particle-intervals.\n";
#ifdef TASDEBUG
        cout << hardParticles[0] << " " << hardParticles[1] << " " << hardParticles[2] << endl;
        cout << hardFrames[0] << " " << hardFrames[1] << " " << hardFrames[2] << endl;
#endif
	if(Samples > hardParticles.size()) { cout << "Insuff. hard particles."<< endl; }//return 1;}
	cout << "Pruning list of slow moving particles.\n";	


        // Helper array to remember which particles we've already selected
	list<int> remainingParticles;
	for(int i = 0 ; i < hardParticles.size() ; i++)
		remainingParticles.push_back(i);

	set<pair<int,int> > tmpset;

	for(int i = 0 ; i < Samples ; i++)
	{       // Go to a random index and move it into trainingArrays, then it erase from remainingParticles
		int index = random.randInt(remainingParticles.size()-1);
		list<int>::iterator it = remainingParticles.begin();
		for(int j = 0 ; j < index ; j++)
			it++;
		trainingParticles.push_back(hardParticles[*it]);
		trainingFrames.push_back(hardFrames[*it]);
		tmpset.insert(pair<int,int>(hardParticles[*it],hardFrames[*it]));
		remainingParticles.erase(it);
	}
	
	cout << "Successively selected " << tmpset.size() << " hard-particle slow-intervals for training\n";
//	cout << "Soft particles = ";
///	for(int i = 0 ; i < Samples ; i++)
//		cout << "[" << trainingParticles[i] << ", " << trainingFrames[i] << "], ";
/////	cout << endl;

//	cout << "Hard particles = ";
//	for(int i = Samples ; i < Samples*2 ; i++)
//		cout << "[" << trainingParticles[i] << ", " << trainingFrames[i] << "], ";
//	cout << endl;


	cout << "Completed identifying training particles. Opening NC file.\n";	
	sprintf(filename,"%s/Training_Nt%i_N%i_tag%s.nc",DataFilePath,Number,Samples,OutFileTag.c_str());
        cout << " Will put into Training NC file " << filename << endl; // TAS

	double RadDistCutoff = settings["RadDistCutoff"];
	double AngularCutoff = settings["AngularCutoff"];
	vector<double> Mu = settings["RadialMu"];
	vector<double> Eta = settings["AngularEta"];
        //cout << "Presumably eta is supposed to be doubles" << endl;
	int Number_Symmetry = Mu.size() + Eta.size() + IncludeCNA;

	NcFile Training_Out(filename,NcFile::Replace);

	if(!Training_Out.is_valid())
	{
		cout << "ERROR Cannot open training set file.\n";
		return 1;
	}

	NcDim *Samples_Dim = Training_Out.add_dim("Samples",2*Samples);
	NcDim *Symmetry_Dim = Training_Out.add_dim("Symmetry",Number_Symmetry);
	NcVar *Type_Var = Training_Out.add_var("Type",ncInt,Samples_Dim); // TAS refers to whether movable or pinned particles
	NcVar *Symmetry_Var = Training_Out.add_var("SymmetryFunctions",ncDouble,Samples_Dim,Symmetry_Dim);


	// Now go through the training set, construct the symmetry functions and write it out
        // All we have is trainingParticles and trainingFrames.  
        // The first "Samples" entries are supposed to be soft, the second are hard.
        // We sort them by the frame number we have to open.  
        // Sorting particles by frames requires this round-about, and maybe c++11
        vector<int> sortedTrainingType = sortA_typesAndB_basedOnB(trainingParticles, trainingFrames);

        for (int i = 0; i < min(25,(int)(trainingParticles.size())); i++) 
                cout << "sorted" << sortedTrainingType[i] << endl;


	CStaticState<DIM> State(Number);
	CSimpleGrid<DIM> *Grid = NULL;
	list<list<CSimpleNeighbor<DIM> > > NeighborList;

        // (Faster might be to sort them by which frame the training particle is in
        //  Then open the file once and read out the training atoms of that frame  ) // TAS

	//compute all the symmetry functions for the current frame.
	for(int i =0 ; i < trainingParticles.size() ; i++)
	{

                // Output a sign of progress
                if ((i == trainingParticles.size()/2) || (i == trainingParticles.size()/20) || 
                    (i == trainingParticles.size()/200) || (i == trainingParticles.size()/2000) || i<4)
                cout << "On particle " << i << " of " << trainingParticles.size() << endl;
                //if (i % 1000 == 0) cout << "Structure functions from particle " << 
                //        trainingParticles[i] << endl << 
                //        " (pinned or rearranged) in frame: " << filename << endl; // TAS


                // Adding this if statement is the whole pay-off of sorting by frame
                if ((i == 0) || (trainingFrames[i] != trainingFrames[i-1])) {
                        if (trainingFrames[i] < trainingFrames[i-1]) { 
                                printf("Error: Frames supposed to be sorted \n"); return 1; 
                        }

                        sprintf(filename,"%s/%s.%07i.%s",DataFilePath,InputFilePrefix,trainingFrames[i],InputFileSuffix);
		        if (!InputLAMMPS_Dump(filename,State)) { 
                                printf("\n ERROR: Reading dump file for trainingparticle%i \n", i); return 1; 
                        }
		        if(Grid==NULL){
		                cout << "New grid needs all particles to be in box [0,1) so now PBC-wrapped at input in InputLAMMPS_Dump" << endl;
                                //printf("TAS_label cutoff %f \n", cutoff);
			        Grid = new CSimpleGrid<DIM>(State,max(AngularCutoff,RadDistCutoff)); // TAS made sure cut off includes all neighbors
			        //AngularGrid = new CSimpleGrid<DIM>(State,AngularCutoff);
		        } else {
			        Grid->Repopulate(State);
			        //AngularGrid->Repopulate(State);
		        }
                }

                // Create a list consisting of a single particle so that we can get its neighbors
		list<int> particles;
		particles.push_back(trainingParticles[i]);
		NeighborList.clear();
		Grid->PopulateNeighborList(State,RadDistCutoff,particles,NeighborList);

		list<vector<double> > symmetry;

		if (!ComputeSymmetryFunctions(State,NeighborList,settings,particles,symmetry)) return 1;

		// write data to a file.
                // NSamples  of "moving" and also NSamples of "pinned"
		Type_Var->set_cur(i);
		// int type = i < Samples ? 1 : 0;
		int type = sortedTrainingType[i];
		Type_Var->put(&type,1);

		Symmetry_Var->set_cur(i);
		Symmetry_Var->put(&symmetry.front()[0],1,Number_Symmetry);

		Training_Out.sync();
	}

        return 0;
}

int main(int argc, char* argv[])
{
	//Disable harsh errors.
	NcError nc_err(NcError::silent_nonfatal);

	cout << "Starting Analysis...\n";

	ifstream input;

	if(argc==2){
		input.open(argv[1]);
	}else{
		cout << "Not enough arguments specified.\n";
		return 0;
	}

	if(!input.is_open())
	{
		cout << "Cannot Open File\n";
		return 0;
	}

	//Add settings items and read the settings file
	Settings settings;
	//The number of particles
	settings.AddItem("Number",INTEGER);
	//The temperature
	settings.AddItem("OutFileTag",STRING);
	//The number of the seed to use
	settings.AddItem("Seed",INTEGER);
	settings.AddItem("Samples",INTEGER);
	//symmetry function stuff
	//the number of samples to put in our set
	settings.AddItem("SlowTime",INTEGER);
	settings.AddItem("ExcitationCutoff",DOUBLE);
	settings.AddItem("HopWindowNFrames",INTEGER);
	settings.AddItem("HopCutoff",DOUBLE);
	//settings.AddItem("MD_Seed",INTEGER);
	settings.AddItem("StartFrame",INTEGER);
	settings.AddItem("EndFrame",INTEGER);
        //symmetry function stuff
        settings.AddItem("RadialMu",vector<double>());
        settings.AddItem("RadialSigma",DOUBLE);
        settings.AddItem("AngularEta",vector<double>());
        settings.AddItem("Zeta",vector<double>());
        settings.AddItem("Lambda",vector<double>());
        settings.AddItem("RadDistCutoff",DOUBLE);
        settings.AddItem("AngularCutoff",DOUBLE);
	settings.AddItem("DataFilePath",STRING);	
	settings.AddItem("InputFilePrefix",STRING);	
	settings.AddItem("InputFileSuffix",STRING);	
	settings.AddItem("OutFileTag",STRING);	
        settings.AddItem("Frames",vector<int>());
	settings.AddItem("IncludeCNA",INTEGER);
	input >> settings;
	cout << "DONE IMPORTING SETTINGS\n";
	cout.flush();
	input.close();


	if (ConstructTrainingSet(settings)) return -1;

        printf("Completed ConstructTrainingSet\n");
	return 0;
}

