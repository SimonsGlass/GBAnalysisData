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
#include "Resources/svm.h"

using namespace LiuJamming;
#define DIM 3
using namespace std;
typedef Eigen::Matrix<dbl, DIM, 1> dvec;
typedef Eigen::Matrix<dbl,DIM,DIM> dmat;
typedef MatrixInterface<dbl> MI;

using namespace std;

//#define TASDEBUG

int trainSVM(Settings &settings)
{
	char filename[512];

	//Read in settings
	int Number = settings["Number"];
	string OutFileTag = settings["OutFileTag"];
	int Samples = settings["Samples"];
	string Seed = settings["Seed"];
	string DataFilePath = settings["DataFilePath"];
	vector<int> ExcludedSymmetry = settings["ExcludedFunctions"];

	if(ExcludedSymmetry.size() <= 1) {
                cout << "TAS No symmetry functions will be excluded" << endl;
		ExcludedSymmetry.clear();
        }

	//perform cross validation if > 0
	int CrossValidation = settings["CrossValidation"];
	//sets the error penalty, C.
	double C = settings["C"];
	double gamma = settings["Gamma"];
	string kernel = settings["Kernel"];

	//Load the symmetry function database
	sprintf(filename,"%s/Training_Nt%i_N%i_tag%s.nc",DataFilePath,Number,Samples,OutFileTag.c_str());
        cout << "Will read training set " << filename << endl;
	NcFile TrainingSet_In(filename);
        cout << "Completed reading training set " << filename << endl;

	int NumberSymmetry = TrainingSet_In.get_dim("Symmetry")->size();
        cout << NumberSymmetry << "is Number of Symmetry functions " << endl;
	NcVar *Type_Var = TrainingSet_In.get_var("Type");
	NcVar *Symmetry_Var = TrainingSet_In.get_var("SymmetryFunctions");

	//Construct the svm problem
	svm_problem problem_data;
	problem_data.l = TrainingSet_In.get_dim("Samples")->size();
	problem_data.y = new double[problem_data.l];
	problem_data.x = new svm_node*[problem_data.l];

	Type_Var->get(problem_data.y,problem_data.l);
        cout << "TAS_label " << problem_data.l << " " << problem_data.y[0] << " " << problem_data.y[1] << endl;
        cout << "Sample y data: " << endl;
        for (int i =0; i< min(10,problem_data.l); i++) // DEBUG
                cout << " " << problem_data.y[i] << endl;

	double symmetry[NumberSymmetry];

	cout << "Preparing problem.\n";

	//load the symmetry functions and normalize them to have mean zero unit variance
        cout << "size of excluded symm is " << ExcludedSymmetry.size() << endl;
	double SymmetryMean[NumberSymmetry - ExcludedSymmetry.size()];
	double SymmetryMeanSq[NumberSymmetry-ExcludedSymmetry.size()];

	for(int i = 0 ; i < NumberSymmetry-ExcludedSymmetry.size() ; i++)
	{
		SymmetryMean[i] = 0.0;
		SymmetryMeanSq[i] = 0.0;
	}

        // For each of the training-particles
	for(int i = 0 ; i < problem_data.l ; i++)
	{
		Symmetry_Var->set_cur(i);
		Symmetry_Var->get(symmetry,1,NumberSymmetry);
	
		problem_data.y[i]++;
	
                // Why is there a + 1 for the number of nodes?  Ah, to hold a -1 in the last element
		problem_data.x[i] = new svm_node[NumberSymmetry-ExcludedSymmetry.size()+1];
		int it = 0;
                // For each structure function
		for(int j = 0 ; j < NumberSymmetry ; j++)
		{
			bool acceptedFunction = true;
			for(int k = 0 ; k < ExcludedSymmetry.size() ; k++)
				if(ExcludedSymmetry[k] == j)
					acceptedFunction = false;
		
                        // if not an excluded structure function, find its mean value
			if(acceptedFunction){
				SymmetryMean[it] += symmetry[j] / problem_data.l; // this structure functino has this mean (avgd over particles)
                        	SymmetryMeanSq[it] += symmetry[j]*symmetry[j] / problem_data.l;

				problem_data.x[i][it].index = it + 1;
				problem_data.x[i][it].value = symmetry[j];
#ifdef TASDEBUG
                                printf("Giving to SVM pt %i the value %f\n",i,symmetry[j]);
#endif                                
				it++;
			}		
		}
                //cout << "means are" << SymmetryMean[0] << SymmetryMean[1] << endl;

		problem_data.x[i][NumberSymmetry-ExcludedSymmetry.size()].index = -1;
	}

        // This will give the variance of the symmetry
	for(int j =0 ; j < NumberSymmetry-ExcludedSymmetry.size() ; j++)
		SymmetryMeanSq[j] -= SymmetryMean[j]*SymmetryMean[j];

	for(int i = 0 ; i < problem_data.l ; i++)
	{
		for(int j= 0 ; j < NumberSymmetry - ExcludedSymmetry.size() ; j++)
		{
                        if (SymmetryMeanSq[j] > 0)
                    	    problem_data.x[i][j].value = (problem_data.x[i][j].value - SymmetryMean[j])/sqrt(SymmetryMeanSq[j]);
			//cout << problem_data.x[i][j].value << " : " <<  problem_data.x[i][j].index << ", ";
		}
		//cout << endl << "------------\n" << problem_data.y[i] << " : " << problem_data.x[i][0].value << endl; 
		//cout << endl;
		
	}
#ifdef TASDEBUG
        cout << "prob data x[0][0] is  " << problem_data.x[0][0].value << endl;
#endif        

	svm_parameter problem_parameters;
	problem_parameters.svm_type = C_SVC;
	problem_parameters.kernel_type = kernel.compare("LINEAR") == 0 ? LINEAR : RBF;
        cout << "Kernel:"; if (problem_parameters.kernel_type == LINEAR) cout << "LINEAR"; else cout << "RBF"; cout << endl;
	problem_parameters.cache_size = 2000;
	problem_parameters.eps = 1.0e-3;
        cout << "We are sending 2^C and s^gamma before sending to SVM" << endl;
	problem_parameters.C = pow(2.0,C);
	problem_parameters.gamma = pow(2.0,gamma);
	cout << "C = " << problem_parameters.C << endl;
	problem_parameters.nr_weight = 0;
	problem_parameters.weight_label = NULL;
	problem_parameters.weight = NULL;
	problem_parameters.probability = 0;
	problem_parameters.shrinking = 1;
	problem_parameters.nu = 0.5;
	problem_parameters.degree = 3;
	problem_parameters.coef0 = 0;
	problem_parameters.p = 0.1;

	cout << "Beginning training.\n";

	//If we are performing cross validation then just compute the accuracy probabilities, but don't bother saving the model
	if(CrossValidation > 0)
	{
		double target[problem_data.l]; // This should be filled with 1 and 2..
		
		if(!svm_check_parameter(&problem_data,&problem_parameters)==NULL)
		{
			cout << "Error: " << svm_check_parameter(&problem_data,&problem_parameters) << endl;
		}

		svm_cross_validation(&problem_data,&problem_parameters,CrossValidation,target);
	        
		double success = 0.0;
		double softSuccess = 0.0;
		double hardSuccess = 0.0;

                for(int i = 0 ; i < problem_data.l ; i++) {

#ifdef TASDEBUG
                        cout << "SVMpt " << i << "target=" << target[i];
                        cout << " with g(r=1.26)= " << problem_data.x[i][0].value << "  \tis ";
                        if(problem_data.y[i] != target[i]) cout << "in";
                        cout << "correctly identified as ";
                        if (problem_data.y[i] == 2) cout << "soft"; else cout << "hard"; 
                        cout << endl;
#endif                        

			if(problem_data.y[i] == target[i])
				success++;
		        
			if(problem_data.y[i] == 2 && problem_data.y[i] == target[i])
				softSuccess++;
			
			if(problem_data.y[i] == 1 && problem_data.y[i] == target[i])
				hardSuccess++;
		}

		cout << "Cross validation accuracy at C = " << C << " = " << success/problem_data.l << endl;

		cout << "Hard accuracy = " << hardSuccess*2.0/problem_data.l << endl;
		cout << "Soft accuracy = " << softSuccess*2.0/problem_data.l << endl;

	} else {
	//If we are not performing cross validation then construct a model and save it		
		svm_model *model = svm_train(&problem_data,&problem_parameters);

		sprintf(filename,"%s/KA_N%i_tag%s_SVMlibsvm_model_N%i_S%s_C%f",DataFilePath,Number,OutFileTag.c_str(),Samples,Seed.c_str(),C);
	        printf("We are writing model to file %s \n",filename);
		//save the model
		if(svm_save_model(filename,model))
		{
			cout << "Error: could not save svm model.\n";
			return;
		}

		sprintf(filename,"%s/KA_N%i_tag%s_SVMmodel_N%i_S%s_C%f.nc",DataFilePath,Number,OutFileTag.c_str(),Samples,Seed.c_str(),C);
	
		//compute the normal
		double w[NumberSymmetry-ExcludedSymmetry.size()];
		for(int j = 0 ; j < NumberSymmetry - ExcludedSymmetry.size() ; j++)
			w[j] = 0.0;
		
		for(int i = 0 ; i < model->l ; i ++)
			for(int j = 0 ; j < NumberSymmetry - ExcludedSymmetry.size() ; j++)
			{
				w[j] += model->sv_coef[0][i] * model->SV[i][j].value;
			}

		NcFile Model_Out(filename,NcFile::Replace);
		
		NcDim *SymmetryDim = Model_Out.add_dim("SymmetryFunctions",NumberSymmetry-ExcludedSymmetry.size());
		NcDim *ExcludedDim = Model_Out.add_dim("ExcludedSize",ExcludedSymmetry.size());		

		Model_Out.add_var("Mean",ncDouble,SymmetryDim)->put(SymmetryMean,NumberSymmetry-ExcludedSymmetry.size());
		Model_Out.add_var("Variance",ncDouble,SymmetryDim)->put(SymmetryMeanSq,NumberSymmetry-ExcludedSymmetry.size());
		Model_Out.add_var("W",ncDouble,SymmetryDim)->put(w,NumberSymmetry-ExcludedSymmetry.size());
		Model_Out.add_var("Excluded",ncInt,ExcludedDim)->put(&ExcludedSymmetry[0],ExcludedSymmetry.size());
		double bias = - model->rho[0];
		Model_Out.add_att("Bias",bias);

		Model_Out.close();
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
	settings.AddItem("Number",INTEGER);
        settings.AddItem("OutFileTag",STRING);
	settings.AddItem("Seed",STRING);
	settings.AddItem("Samples",INTEGER);
	settings.AddItem("ExcludedFunctions",vector<int>());
	settings.AddItem("CrossValidation",INTEGER);
	settings.AddItem("C",DOUBLE);
	settings.AddItem("Gamma",DOUBLE);
	settings.AddItem("Kernel",STRING);
	settings.AddItem("DataFilePath",STRING);

	input >> settings;
        cout << "DONE IMPORTING SETTINGS\n";
        cout.flush();
        input.close();


	if (trainSVM(settings)) return -1;

        printf("Completed TrainSVM\n");
        return 0;
}

