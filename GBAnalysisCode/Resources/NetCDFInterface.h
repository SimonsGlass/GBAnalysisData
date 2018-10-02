#include <netcdfcpp.h>
#include "../Resources/std_include.h"

//Load a vector of class T from a netcdf file of with name "varname". 
//It is assumed that T is the same as the format of the data saved in
//the netcdf file.
template <class T>
bool LoadVector(const char *filename, const char *varname, vector<T> &vec)
{
	NcFile file(filename,NcFile::ReadOnly);

	if(!file.is_valid())
	{
		cout << "Cannot open NetCDF file: " << filename << endl;
		return false;
	}

	NcVar *variable = file.get_var(varname);
	int size = variable->get_dim(0)->size();
	
	vec.resize(size);
	
	variable->get(vec.data(),size);

	file.close();

	return true;
}


template <class T>
bool LoadMatrix(const char *filename, const char *varname, vector<vector<T> > &mat)
{
	NcFile file(filename,NcFile::ReadOnly);

        if(!file.is_valid())
        {
                cout << "Cannot open NetCDF file: " << filename << endl;
                return false;
        }

        NcVar *variable = file.get_var(varname);
	int rows = variable->get_dim(0)->size();
	int columns = variable->get_dim(1)->size();        

        mat.resize(rows);
	for(int i = 0 ; i < rows ; i++){
		mat[i].resize(columns);
		variable->set_cur(rows);
        	variable->get(mat[i].data(),1,columns);
	}

        file.close();

        return true;
}


