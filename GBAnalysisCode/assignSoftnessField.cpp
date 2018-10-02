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

//constructs symmetry functions for a balanced set of particles of a given size randomly selected from the 
int ConstructSoftnessField(Settings &settings)
{

        // Read in settings
	int Number = settings["Number"];
	string OutFileTag = settings["OutFileTag"];
	string Seed = settings["Seed"];
	int Samples = settings["Samples"];
	double RadDistCutoff = settings["RadDistCutoff"];
	double AngularCutoff = settings["AngularCutoff"];
	double C = settings["C"];
	string DataFilePath = settings["DataFilePath"];
	string InputFilePrefix = settings["InputFilePrefix"];
	string InputFileSuffix = settings["InputFileSuffix"];
	vector<int> SoftnessFrames = settings["SoftnessFrames"];


	cout << "About to load model";
	char freadname[256],fwritename[256];
	sprintf(freadname,"%s/KA_N%i_tag%s_SVMmodel_N%i_S%s_C%f.nc",DataFilePath,Number,OutFileTag.c_str(),Samples,Seed.c_str(),C);
        cout << " from file " << freadname << endl;  cout.flush();

	NcFile Model_In_A(freadname);	
	int N_Model_Symmetry_A = Model_In_A.get_dim("SymmetryFunctions")->size();
	double Mean_A[N_Model_Symmetry_A];
	Model_In_A.get_var("Mean")->get(Mean_A,N_Model_Symmetry_A);
	double Variance_A[N_Model_Symmetry_A];
	Model_In_A.get_var("Variance")->get(Variance_A,N_Model_Symmetry_A);
	double W_A[N_Model_Symmetry_A];
	Model_In_A.get_var("W")->get(W_A,N_Model_Symmetry_A);
	int Excluded_A[Model_In_A.get_dim("ExcludedSize")->size()];
	Model_In_A.get_var("Excluded")->get(Excluded_A,Model_In_A.get_dim("ExcludedSize")->size());
	double Bias_A;
	Bias_A = Model_In_A.get_att("Bias")->as_double(0);
	cout << "Loaded model with bias " << Bias_A << endl; // " and " << Bias_A << endl; 


	//construct a list of all particle ids
	list<int> particles;
	for(int i = 0 ; i < Number ; i++)
		particles.push_back(i);
	
#ifdef TASDEBUG
        cout << "temp01a"  << Number << endl; // " and " << Bias_A << endl; 
#endif

	CStaticState<DIM> State(Number);
	CSimpleGrid<DIM> *Grid = NULL;
	CSimpleGrid<DIM> *AngularGrid = NULL;

	cout << "temp02"  << endl; // " and " << Bias_A << endl; 

	//now go through the frames for which we want softness and calculate the symmetry functions
	for(int fr = SoftnessFrames[0] ; fr <= SoftnessFrames[1] ; fr+= SoftnessFrames[2])
	{
		cout << "Will calculate softness of frame: " << fr << endl;
                sprintf(freadname,"%s/%s.%07i.%s",DataFilePath,InputFilePrefix,fr,InputFileSuffix);
                printf("Reading file %s \n", freadname);
		if (!InputLAMMPS_Dump(freadname,State)) {printf("ERROR Reading %s \n",freadname); return;}

                // Construct a grid so we can consider only atoms' neighbors
		cout << "Constructing grid: " << endl;
		if(Grid==NULL)
		{
			Grid = new CSimpleGrid<DIM>(State,RadDistCutoff);
			AngularGrid = new CSimpleGrid<DIM>(State,AngularCutoff);
		} else {
			Grid->Repopulate(State);
			AngularGrid->Repopulate(State);
		}
		Grid->ConstructNeighborList(State,RadDistCutoff);
		AngularGrid->ConstructNeighborList(State,AngularCutoff);

		
		cout << "Computing symmetry functions.\n";
		list<vector<double> > symmetry;
		if (!ComputeSymmetryFunctions(State,Grid,AngularGrid,settings,particles,symmetry)) return;
		
		sprintf(fwritename,"%s/KA_N%i_tag%s_Softness.%i",DataFilePath,Number,OutFileTag.c_str(),fr);
		ofstream Softness_Out(fwritename);
		cout << "Computing softness of each particle \n";
		int i = 0;
		//compute the softness.
		for(list<vector<double> >::iterator particle_it = symmetry.begin() ; particle_it != symmetry.end() ; particle_it++)
		{
			double dp  = 0.0;
			int it = 0;
			for(int j =0 ; j< (*particle_it).size() ; j++)
			{
                                //if (j % (*particle_it).size()/2 == 0) printf("Considering a symmetry of particle %i in frame %i\n",j,fr);
				bool bExcluded = false;
				int Excluded_Size = Model_In_A.get_dim("ExcludedSize")->size();
				for(int k = 0 ; k < Excluded_Size ; k++)
					if(Excluded_A[k]==j)
						bExcluded = true;

				if(!bExcluded)
				{
				//	cout << "Symmetry = " << (*particle_it)[j] << " Mean = " << Mean_A[it] << " Std_Dev = " << sqrt(Variance_A[it]) << endl;
				//	cout << "Sum = " << (((*particle_it)[j]-Mean_A[it])/sqrt(Variance_A[it])) << " W = " << W_A[it] << " dp = " << dp << endl;
					dp +=	((*particle_it)[j]-Mean_A[it])/sqrt(Variance_A[it]) * W_A[it];
					it++; 
				}
			}

			Softness_Out << (dp + Bias_A) << endl;
			i++;
		}

		Softness_Out.close();
	}
	
	return 0;
}

int main(int argc, char* argv[])
{
        //Disable harsh errors.
        NcError nc_err(NcError::silent_nonfatal);

        cout << "Assigning estimated softness...\n";

        ifstream input;

        if(argc==2){
                input.open(argv[1]);
        } else {
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
	//settings.AddItem("Number_A",INTEGER);	
        //The temperature
        settings.AddItem("Temperature",STRING);
        //The number of seeds to use
        settings.AddItem("Seed",STRING);
	//symmetry function stuff
	settings.AddItem("RadialMu",vector<double>());
	settings.AddItem("RadialSigma",DOUBLE);
	settings.AddItem("AngularEta",vector<double>());
	settings.AddItem("Zeta",vector<double>());
	settings.AddItem("Lambda",vector<double>());
	settings.AddItem("RadDistCutoff",DOUBLE);
	settings.AddItem("AngularCutoff",DOUBLE);
	//the number of samples to put in our set
	settings.AddItem("SoftnessFrames",vector<int>());
	settings.AddItem("Samples",INTEGER);
	settings.AddItem("C",DOUBLE);
	settings.AddItem("DataFilePath",STRING);
	settings.AddItem("InputFilePrefix",STRING);
	settings.AddItem("InputFileSuffix",STRING);
	settings.AddItem("OutFileTag",STRING);
 	settings.AddItem("IncludeCNA",INTEGER);
        input >> settings;
        cout << "DONE IMPORTING SETTINGS\n";
        cout.flush();
        input.close();


       	if (ConstructSoftnessField(settings)) return -1;

        printf("Completed assignSoftnessField.\n");
        return 0;
}

