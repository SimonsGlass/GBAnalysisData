#ifndef LOAD_SOFTNESS

#define LOAD_SOFTNESS

#include "std_include.h"


void LoadSoftness(string Filename, Eigen::VectorXd &Softness)
{
	ifstream Softness_In(Filename.c_str());

	int i;

	//while(!D2Min_In.eof())
	for(int t = 0 ; t < Softness.rows() ; t++)
	{
		Softness_In >> Softness(t);	
	}
	Softness_In.close();
}

void LoadSoftness(string Filename, Eigen::VectorXd &Softness,int Number_A, double Scale_A, double Scale_B)
{
        ifstream Softness_In(Filename.c_str());
	if(!Softness_In.is_open()) {
		cout << "Cannot open file: " << Filename << endl;
		return ;
	}

        int i;

        //while(!D2Min_In.eof())
        for(int t = 0 ; t < Softness.rows() ; t++)
        {
                Softness_In >> Softness(t);
                Softness(t) *= t < Number_A ? Scale_A : Scale_B ;
        }
        Softness_In.close();
}

template <int Dim>
void LoadGradients(string Filename, vector<Eigen::Matrix<double,Dim,1> > &Gradient)
{
        ifstream Gradient_In(Filename.c_str());

	if(Gradient.size() > 0)
		Gradient.clear();

        int i;

        while(!Gradient_In.eof())
        {
		Eigen::Matrix<double,Dim,1> vector;
                Gradient_In >> vector(0) >> vector(1) >> vector(2);
        	if(Gradient_In.eof())
			return;
		Gradient.push_back(vector);
	}
        Gradient_In.close();
}

#endif
