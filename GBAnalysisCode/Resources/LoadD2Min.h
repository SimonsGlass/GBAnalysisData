#ifndef LOAD_D2MIN

#define LOAD_D2MIN

#include "std_include.h"


void LoadD2Min(string Filename, Eigen::VectorXd &D2Min)
{
	ifstream D2Min_In(Filename.c_str());

	int i;

	//while(!D2Min_In.eof())
	for(int t = 0 ; t < D2Min.rows() ; t++)
	{
		D2Min_In >> i;
		D2Min_In >> D2Min(i);	
	}
	D2Min_In.close();
}

void LoadD2MinIso(string Filename, Eigen::VectorXd &D2Min, Eigen::VectorXd &Displacement)
{
	ifstream D2Min_In(Filename.c_str());

        //while(!D2Min_In.eof())
        for(int t = 0 ; t < D2Min.rows() ; t++)
        {
                D2Min_In >> D2Min(t) >> Displacement(3*t) >> Displacement(3*t+1) >> Displacement(3*t+2);
        }
        D2Min_In.close();
}

#endif
