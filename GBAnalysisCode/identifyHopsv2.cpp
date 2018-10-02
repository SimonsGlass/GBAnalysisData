#include <iostream>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Resources/SimpleGrid.h"
#include "Resources/Settings.h"
#include "Computers/StaticComputer.h"
#include "Computers/WallComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
#include "Database/Database.h"
#include "Computers/SoftSpotComputer.h"
#include "Computers/BarrierHeight.h"
#include "Resources/OutputLAMMPS.h"
#include "Resources/LoadD2Min.h"
#include "Resources/ComputeSymmetryFunctions.h"
#include "Resources/LoadSoftness.h"

using namespace LiuJamming;
#define DIM 3
using namespace std;
typedef Eigen::Matrix<dbl, DIM, 1> dvec;
typedef Eigen::Matrix<dbl,DIM,DIM> dmat;
typedef MatrixInterface<dbl> MI;

using namespace std;


//#define TASDEBUG

//sort based on frame
bool Sort(const pair<pair<int,int>, int> &first, const pair<pair<int,int>, int> &second)
{
	return first.first.first < second.first.first;	
}


void ComputeDisplacement(int i, CStaticState<DIM> &StateI, int j, CStaticState<DIM> StateJ, dvec &Vector)
{
        dvec PosI, PosJ;
        StateI.GetParticlePositionVirtual(PosI,i);
        StateJ.GetParticlePositionVirtual(PosJ,j);

        StateI.GetBox()->MinimumDisplacement(PosI,PosJ,Vector);
        StateI.GetBox()->Transform(Vector);
}


void SumPositions(CStaticState<DIM> &State, Eigen::VectorXd &mean)
{
	Eigen::VectorXd positions,displacement;
	State.GetPositionsVirtual(positions);
	State.GetBox()->MinimumDisplacement(mean,positions,displacement);
	mean -= displacement;
}



int IdentifyExcitations(Settings &settings)
{
	char freadname[256],fwritename[256];

	int Number = settings["Number"];
	//string OutFileTag = settings["OutFileTag"];
	double HopCutoff = settings["HopCutoff"]; // phop_c
	//the instanton time
        // In the "start" and "end" time arrays, we adjust the time by shift, which must be then the t_R/2 value.
	vector<int> Frames = settings["Frames"]; // Available file name ids from Frames[0] to [1] skipping by [2] to input to phop calculation
	// Half-window for phop calculation, t_R/2, originally reads-in units of Frames[2], then femtosec
	// the plateau time (in units of Frames[2]) 
	int HopWindowNFrames = settings["HopWindowNFrames"]; // 5 for example
	string DataFilePath = settings["DataFilePath"];
	string InputFilePrefix = settings["InputFilePrefix"];
	string InputFileSuffix = settings["InputFileSuffix"];
#ifdef TASDEBUG        
        printf("TAS_label_00018 \n");fflush(stdout); 
#endif
        // Times are in units of either file spacing or file labels.  Now it is switched to file labels.
	int HopWindow = HopWindowNFrames * Frames[2]; // now units of femtoseconds (file id units)


	cout << "Loading initial state.\n";
	cout << "Beginning analysis.\n";

	//setup output file
	sprintf(fwritename,"%s/keyschandlerKA_N%i_tag_AnalysisHops_w%i_c%.2f.nc",DataFilePath,Number,HopWindow,HopCutoff);

	NcFile Output_File(fwritename,NcFile::Replace);
	
	NcDim *Excitation_Dim = Output_File.add_dim("NumberOfExcitations");
	NcDim *Spatial_Dim = Output_File.add_dim("SpatialDimension",DIM);	

        NcVar *Start_Var = Output_File.add_var("Start",ncInt,Excitation_Dim);
	NcVar *Particle_Var = Output_File.add_var("Particle",ncInt,Excitation_Dim);
	NcVar *End_Var = Output_File.add_var("End",ncInt,Excitation_Dim);
	NcVar *StartPosition_Var = Output_File.add_var("Start_Position",ncDouble,Excitation_Dim,Spatial_Dim);
	NcVar *EndPosition_Var = Output_File.add_var("End_Position",ncDouble,Excitation_Dim,Spatial_Dim);   
	NcVar *Max_Var = Output_File.add_var("p_hop",ncDouble,Excitation_Dim);
	NcVar *StartMeanPosition_Var = Output_File.add_var("Start_Mean_Position",ncDouble,Excitation_Dim,Spatial_Dim);
	NcVar *EndMeanPosition_Var = Output_File.add_var("End_Mean_Position",ncDouble,Excitation_Dim,Spatial_Dim);	
	
	//set up states
	list<Eigen::VectorXd> States_I;
	Eigen::VectorXd MeanPositionI = Eigen::ArrayXd::Zero(Number*DIM);
	list<Eigen::VectorXd> States_J;
	Eigen::VectorXd MeanPositionJ = Eigen::ArrayXd::Zero(Number*DIM);

	//note this is a kludge that will not work on non-cubic systems. should fix it at some point.
	double boxsize = 0.0;

	//go through and load all of the states unwrapping the coordinates as we go along
        // Sum positions of each atom trajectory onto MeanPositionI array and 
        // Append position arrays into the [Num window x N particles] array States_I
	for(int dt = 0 ; dt < HopWindow ; dt+=Frames[2])
	{
                // These are in the form id xs ys zs.  Sometimes before xs there is a type column
		sprintf(freadname,"%s/%s.%07i.%s",DataFilePath,InputFilePrefix,Frames[0]+dt,InputFileSuffix);
                cout << "window A file: " << freadname << endl;

                Eigen::VectorXd Positions,Displacement;
		CStaticState<DIM> InherentStateI(Number);
                InputLAMMPS_Dump(freadname,InherentStateI);
		InherentStateI.GetPositionsVirtual(Positions);
		
		//kludge kludge kludge
		if(dt ==0)
		{
			Eigen::Matrix<double,DIM,DIM> Transform;
			InherentStateI.GetBox()->GetTransformation(Transform);
			boxsize = Transform(0,0);
		}	

		if(dt > 0)
		{
			InherentStateI.GetBox()->MinimumDisplacement(Positions,States_I.front(),Displacement);
			Positions = States_I.front() + Displacement;
		}
		MeanPositionI += Positions / HopWindowNFrames;


		States_I.push_back(Positions);


		sprintf(freadname,"%s/%s.%07i.%s",DataFilePath,InputFilePrefix,Frames[0]+HopWindow+dt,InputFileSuffix);
                cout << "window B file: " << freadname << endl;

                CStaticState<DIM> InherentStateJ(Number);
                InputLAMMPS_Dump(freadname,InherentStateJ);
		InherentStateJ.GetPositionsVirtual(Positions);
		InherentStateI.GetBox()->MinimumDisplacement(Positions,States_I.front(),Displacement);
		Positions = States_I.front() + Displacement;
		MeanPositionJ += Positions / HopWindowNFrames;

		States_J.push_back(Positions);
	}

	//a list of all the particles currently experiencing an exciation. currentExcitation[i] = -1 if particle
	//i is not experiencing an excitation and currentExcitation[i] = frameStart if particle started
	// experiencing an excitation at frame frameStart.
#ifdef TASDEBUG
        printf("TAS_label_00300 Number %i\n",Number);fflush(stdout);
#endif        
	double* currentMax = new double[Number]; // TAS seg fault stack memory // double currentMax[Number];
	for(int i = 0 ; i < Number; i++)
		currentMax[i] = -1;

        // TAS If memory will exceed that allocated to the stack, need the below 
        // syntax rather than int currentLocation[Number];
	int* currentLocation = new int[Number]; 
	int* startTime = new int[Number];
	int* endTime = new int[Number];
	for(int i = 0 ; i < Number ; i++)
		currentLocation[i] = -1;

	Eigen::Matrix<double,DIM,1> *startPosition = new Eigen::Matrix<double,DIM,1>[Number];
	Eigen::Matrix<double,DIM,1> *endPosition = new Eigen::Matrix<double,DIM,1>[Number];
	Eigen::Matrix<double,DIM,1> *startMeanPosition = new Eigen::Matrix<double,DIM,1>[Number];
	Eigen::Matrix<double,DIM,1> *endMeanPosition = new Eigen::Matrix<double,DIM,1>[Number];
	for(int i = 0 ; i < Number ; i++)
	{
		startPosition[i] = Eigen::Array<double,DIM,1>::Zero();
		endPosition[i] = Eigen::Array<double,DIM,1>::Zero();
		startMeanPosition[i] = Eigen::Array<double,DIM,1>::Zero();
		endMeanPosition[i] = Eigen::Array<double,DIM,1>::Zero();
	}

	int nExcitations = 0;

        //go through all of the frames and record any enduring displacements greater in size than a
        // 0 1 2 3 4 5 6 7
        // A A A A B B B B    the first B is the location of the frame as it moves
	for(int fr = Frames[0]+HopWindow ; fr <= Frames[1]-HopWindow+Frames[2] ; fr+= Frames[2])
	{
		cout << "We've opened frames around frame " << fr << " to compute indicator function.\n";
		for(int i = 0 ; i < Number ; i++)
		{
			double var_I = 0.0;
			double var_J = 0.0;
			//compute the variance
			list<Eigen::VectorXd>::iterator it_I = States_I.begin();
			list<Eigen::VectorXd>::iterator it_J = States_J.begin();		
			
			while(it_I != States_I.end())
			{
				var_I += (*it_I - MeanPositionJ).segment<DIM>(DIM*i).squaredNorm()/States_I.size();
				var_J += (*it_J - MeanPositionI).segment<DIM>(DIM*i).squaredNorm()/States_J.size();				
				//if(i==1177)
					//cout << boxsize*(*it_I).segment<DIM>(DIM*i).transpose() << "\n";
				it_I++;
				it_J++;
			}

			
			double p_hop = boxsize*boxsize*sqrt(var_I*var_J);
			
			//check whether particle i is excited
			if(p_hop > HopCutoff && p_hop > currentMax[i])
			{
				if(currentMax[i] < 0)
				{
                                        // Seems that if pchop first exceeds phop_c, the start is assigned to first window A frame.
                                        // startMeanPosition is the average of WindowA at that time.
					startTime[i] = fr - HopWindow;  // I can say start of the event i
#ifdef TASDEBUG
                                        printf("particle %i excited at time %i \n",i,startTime[i]);
#endif                                        
					startPosition[i] = States_I.front().segment<DIM>(DIM*i);
					startMeanPosition[i] = MeanPositionI.segment<DIM>(DIM*i);
				}
				currentMax[i] = p_hop;
				currentLocation[i] = fr;
			}


			//if particle i was excited and then stopped being excited then write it to the file
			if(p_hop < HopCutoff && currentMax[i] > 0)
			{
                                //printf("nExcitations %d \n",nExcitations);fflush(stdout);
				Start_Var->set_cur(nExcitations);
				//int start = currentLocation[i] - HopWindow; // This start is never used - TAS
				Start_Var->put(&startTime[i],1);

				Particle_Var->set_cur(nExcitations);
				Particle_Var->put(&i,1);
			
				End_Var->set_cur(nExcitations);
				int end = fr + HopWindow - Frames[2]; // TAS Now this is the last frame which HopWindow B does reach at the first step that we return to phop<phop_c at window B's first frame
				End_Var->put(&end,1);
		
				StartPosition_Var->set_cur(nExcitations);
				StartPosition_Var->put(startPosition[i].data(),1,DIM);
			
				StartMeanPosition_Var->set_cur(nExcitations);
                                // This is the mean position during WindowA of when phop was first exceeded
				StartMeanPosition_Var->put(startMeanPosition[i].data(),1,DIM);
			
                                // So when windowB first frame has dropped again, last frame in window B gives the endPosition
				EndPosition_Var->set_cur(nExcitations);
				endPosition[i] = States_J.back().segment<DIM>(DIM*i);
				EndPosition_Var->put(endPosition[i].data(),1,DIM);

				EndMeanPosition_Var->set_cur(nExcitations);
				endMeanPosition[i] = MeanPositionJ.segment<DIM>(DIM*i);
				EndMeanPosition_Var->put(endMeanPosition[i].data(),1,DIM);

                                // TAS write current max to the next location in the netcdf file
				Max_Var->set_cur(nExcitations);
				Max_Var->put(&currentMax[i],1); 

				currentMax[i] = -1;
				currentLocation[i] = -1;
				nExcitations++;
			}
		}

	        if (fr+HopWindow  <= Frames[1]) // Only if there is a next frame do we load up a next frame
                {
#ifdef TASDEBUG
                        cout << "Loading new states and adjusting means.\n";
#endif                
    		        MeanPositionI -= States_I.front() / HopWindowNFrames;
    		        MeanPositionI += States_J.front() / HopWindowNFrames;
		        MeanPositionJ -= States_J.front() / HopWindowNFrames;
		
		        Eigen::VectorXd Positions,Displacement;

                        sprintf(freadname,"%s/%s.%07i.%s",DataFilePath,InputFilePrefix,fr+HopWindow,InputFileSuffix);
		        cout << "Window moving to include file " << freadname  << endl;

                        CStaticState<DIM> InherentState(Number);
                        InputLAMMPS_Dump(freadname,InherentState);
                        // Get positions of the frame that is newly in windowB 
                        InherentState.GetPositionsVirtual(Positions);
                        // Correct for wrapping by finding the atoms' unwrapped-Displacement relative to their positions in
                        // the first frame in Window A, and adding it to that reference.
                        InherentState.GetBox()->MinimumDisplacement(Positions,States_I.front(),Displacement);
                
                        Positions = States_I.front() + Displacement;
                        MeanPositionJ += Positions / HopWindowNFrames;

		        States_I.pop_front();
		        States_I.push_back(States_J.front());
		        States_J.pop_front();
		        States_J.push_back(Positions);
                }
	}
	Output_File.close();

        // Since allocated in the explicit syntax, must deallocate memory // TAS
        delete currentMax;
	delete currentLocation;
	delete startTime;
	delete endTime;
        delete startPosition;
        delete endPosition;
        delete startMeanPosition;
        delete endMeanPosition;

        cout << nExcitations << " excitations detected \n";
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
	settings.AddItem("Temperature",STRING);
	//the number of samples to put in our set
	settings.AddItem("Frames",vector<int>());
	settings.AddItem("HopWindowNFrames",INTEGER);
	settings.AddItem("HopCutoff",DOUBLE);
	settings.AddItem("Density",DOUBLE);
	settings.AddItem("DataFilePath",STRING);
	settings.AddItem("InputFilePrefix",STRING);
	settings.AddItem("InputFileSuffix",STRING);
	input >> settings;
	cout << "DONE IMPORTING SETTINGS\n";
	cout.flush();
	input.close();

	if (IdentifyExcitations(settings)) return -1;

        printf("Completed IdentifyExcitations\n");

	return 0;
}

