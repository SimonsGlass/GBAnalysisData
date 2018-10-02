#ifndef LOAD_STRESS

#define LOAD_STRESS

#include "std_include.h"

bool LoadStress(const char *filename, vector<Eigen::Matrix3d> &stress)
{
	ifstream input(filename);

	if(!input.is_open()){
		cout << "Cannot open input file.\n";
		return false;
	}

	char line[256];

	//skip 9 lines
	for(int i = 0 ; i < 9 ; i++)
		input.getline(line,256);

	double c[6];
	int index,type;

	while(!input.eof())
	{
		input >> index >> type >> c[0] >> c[1] >> c[2] >> c[3] >> c[4] >> c[5];

		index--;

		stress[index](0,0) = c[0];
		stress[index](1,1) = c[1];
		stress[index](2,2) = c[2];
		stress[index](1,0) = c[3];
		stress[index](0,1) = c[3];
		stress[index](2,0) = c[4];
		stress[index](0,2) = c[4];
		stress[index](2,1) = c[5];
		stress[index](1,2) = c[5]; 
	}

	return true;
}



#endif
