#ifndef EXCEPTION

#define EXCEPTION

/////////////////////////////////////////////////////////////////////////////////
//Exception class. 
//
//Description
//		This class is provides a format for other classes to throw exceptions
//		when they encounter errors. The information that should be provided
//		is the function that the error occurs in and a brief description of the
//		error. This error is then printed out.
//
//Variables
// 		Function: string specifying the function whence the error.
//		Error: string describing the error
///
//Implements
//		Reading/Writing to netcdf files.
//		Getting/Setting particle positions.
//						particle radii.
//						box.
//						potential.
//
/////////////////////////////////////////////////////////////////////////////////

#include <ostream>
#include <string>

namespace LiuJamming
{

using namespace std;

class CException 
{
private:
	string Function;
	string Error;
	
public:
	CException();
	CException(string _Function, string _Error);
	
	const string &GetFunction() const;
	const string &GetError() const;
	
	friend std::ostream &operator<<(std::ostream &out, const CException &exception);

};

}


#endif