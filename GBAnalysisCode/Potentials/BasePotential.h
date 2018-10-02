#ifndef BASE_POTENTIAL_H
#define BASE_POTENTIAL_H

/////////////////////////////////////////////////////////////////////////////////
//Potential class. 
//
//Description
//		This is a virtual class that describes a generic isotropic potential. 
//		This class implements functions to read and write the potential and its
//		parameters formatted as a string to a NetCDF file. This class provides
//		functions to compute various derivated of the potential with respect 
//		to the distance dr between two particles. This class also has a function
//		that gives the support of the potential. It is generally assumed that all
//		potentials will have compact support.
//
//		Generically the energy will be related to the ratio of the distance between
//		particles to the sum of their radii.
//
//Global Variables
//		A map from string to potential type as an STD map.
//
//Variables
//		
//Implements
//		Reading/Writing to netcdf files.
//
//Virtual Functions
//		Calculating various orders of derivatives of the potential
//		Calculating the support of the potential as a function of radius.
//
//File Formats
//		NetCDF
//			## Note: To denote whether the NetCDF file has been populated with variables,
//			## an attribute "Potential_Populated" gets created.
//			-> One dimensions: records (unlimited).
//			-> One attribute Potential_Populated.
//			-> String of parameters as a variable.
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include "../Resources/Exception.h"
#include "../Resources/Resources.h"
#ifndef DONT_USE_NETCDF
	#include "netcdfcpp.h"
#endif

namespace LiuJamming
{

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//!Abstract base class for pair-potentials.
class CPotential
{
private:
	static std::map<string,CPotential*> PotentialTypes; //!<Global map to reference potentials by name

protected:
	static const int STRING_SIZE = 128;
	static const char STRING_DELIMITER = ':';

public:
//! @name Global functions to read and write
///@{
	static CPotential*	SetFromString(string str);	//!<From an input string, return a pointer to a newly created potential object with the given parameter values
	static string		FillString(stringstream &ss);
#ifndef DONT_USE_NETCDF
	static bool		CheckNetCDF(const NcFile &file);
	static void		PopulateNetCDF(NcFile &file);
	static void		PrepareNetCDF(NcFile &file, bool CheckOnly);
	void			NetCDFWrite(NcFile &file, int record);
	CPotential*		NetCDFRead(NcFile &file, int record);
#endif
///@}

//! @name Functions for setting, copying and saving
///@{
	virtual string		DataToString() const = 0;			//!<Return a string that contains the data necessary for the potential object (e.g. interaction strength, etc.)
 	virtual void		StringToData(string Data) = 0;		//!<From an input string, set the data variables.
 	virtual CPotential* Clone() const = 0;					//!<Return a pointer to a newly created clone of the object.

///@}

//! @name Functions to compute various derivatives
///@{
	virtual dbl  Compute(dbl dr,dbl r) const = 0;													//!<Compute the 0th derivative
	virtual dbl  ComputeFirstDerivative(dbl dr,dbl r) const = 0;									//!<Compute the 1st derivative
	virtual dbl  ComputeSecondDerivative(dbl dr, dbl r) const = 0;									//!<Compute the 2nd derivative
	virtual dbl  ComputeThirdDerivative(dbl dr, dbl r) const = 0;									//!<Compute the 3th derivative
	virtual void ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const = 0;						//!<Compute the 0th and 1st derivatives
	virtual void ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const = 0;			//!<Compute the 0th and 1st and 2nd derivatives
	virtual void ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const = 0;	//!<Compute the 0th and 1st and 2nd and 3rd derivatives
	
///@}

//! @name Misc.
///@{
	//!Get the maximum distance between two particles of unit diameter such that they interact
	virtual dbl ComputeSupport() const = 0;
	
	//!Set sigma and determine if rad1, rad2 and rlen2 are such that there is overlap
	virtual bool Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const = 0;
///@}
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


CPotential* CPotential::SetFromString(string str)
{
	vector<string> split = SplitString(str,STRING_DELIMITER);
	CPotential *pot = PotentialTypes[split[0]]->Clone();
	pot->StringToData(str);
	return pot;
}

string CPotential::FillString(stringstream &ss)
{
	string s = ss.str();
	s.append(STRING_SIZE-s.size(), STRING_DELIMITER);
	assert(s.size()==STRING_SIZE);
	return s;
}

#ifndef DONT_USE_NETCDF
bool CPotential::CheckNetCDF(const NcFile &file)
{
	NcError err(NcError::verbose_nonfatal);
	if(!file.is_valid())	throw(CException("CPotential::CheckNetCDF","NetCDF file is not valid."));

	NcAtt *a = file.get_att("Potential_Populated");
	bool populated = (a==0)?false:true;  //If file.get_att() returns 0 if attribute does not exist.
	delete a;	//attributes need to be deleted by the user to avoid memory leaks.
	return populated;
}

void CPotential::PopulateNetCDF(NcFile &file)
{
	NcError err(NcError::verbose_nonfatal);
	file.add_att("Potential_Populated", 1);

	//Load the record dimension
	NcDim *recDim = file.get_dim("rec");
	if(!recDim->is_valid())	throw(CException("CPotential::PopulateNetCDF","Record dimension is not valid."));

	//Load or create the string size dimension
	NcDim *PotStrSizeDim = file.get_dim("PotStrSize");
	if(PotStrSizeDim==NULL)			PotStrSizeDim = file.add_dim("PotStrSize",STRING_SIZE);
	if(!PotStrSizeDim->is_valid())	throw(CException("CPotential::PopulateNetCDF","PotStrSize dimension is not valid."));

	//Create the variable
	file.add_var("PotString",ncChar,recDim,PotStrSizeDim);
}

void CPotential::PrepareNetCDF(NcFile &file, bool CheckOnly)
{
	NcError err(NcError::silent_nonfatal);
	if(!file.is_valid())	throw(CException("CPotential::CheckNetCDF","NetCDF file is not valid."));

	//Attempt to load the dimensions and variables
	NcDim *recDim = file.get_dim("rec");
	NcDim *PotStrSizeDim = file.get_dim("PotStrSize");
	NcVar *PotStringVar = file.get_var("PotString");

	//Check the record dimension
	if(!recDim->is_valid())	throw(CException("CPotential::PopulateNetCDF","Record dimension is not valid."));

	//Create and check the string size dimension
	if(PotStrSizeDim==NULL && !CheckOnly)	PotStrSizeDim = file.add_dim("PotStrSize",STRING_SIZE);
	if(!PotStrSizeDim->is_valid())			throw(CException("CPotential::PrepareNetCDF","PotStrSize dimension is not valid."));
	if(PotStrSizeDim->size()!=STRING_SIZE)	throw(CException("CPotential::PrepareNetCDF","PotStrSize dimension has wrong size."));
	
	//Create and check the variable
	if(PotStringVar==NULL && !CheckOnly)	PotStringVar = file.add_var("PotString",ncChar,recDim,PotStrSizeDim);
	if(!PotStringVar->is_valid())			throw(CException("CPotential::PrepareNetCDF","PotString variable is not valid."));
}

void CPotential::NetCDFWrite(NcFile &file, int record)
{
//	if(!CheckNetCDF(file))		//Check if necessary dimensions and variables exist
//		PopulateNetCDF(file);	//Create necessary dimensions and variables

	PrepareNetCDF(file,false);
	NcVar *PotStringVar= file.get_var("PotString");
	PotStringVar->put_rec(DataToString().c_str(),   record);
}

CPotential* CPotential::NetCDFRead(NcFile &file, int record)
{
	//Perform checks
//	if(!CheckNetCDF(file))							throw(CException("CPotential::NetCDFRead","Attempting to read a potential from a file that has no appropriate data."));
	
	PrepareNetCDF(file,true);
	if(record>=file.get_dim("rec")->size())		throw(CException("CPotential::NetCDFRead","Attempting to read a potential from a record that does not exist."));
	
	char c_str[STRING_SIZE];
	NcVar *PotStringVar= file.get_var("PotString");
	PotStringVar -> set_cur(record,0);
	
	PotStringVar -> get(&c_str[0], 1, STRING_SIZE);

	string str;
	for(int i=0; i<STRING_SIZE; ++i)
		str.push_back(c_str[i]);
	// The above loop is to correct the following problem: c_str, for whatever reason, doesn't have the actual length of STRING_SIZE, so str can have random garbage at the end of the 128 characters.
//	string str = c_str;

	return SetFromString(str);
}
#endif 


}

#endif //BASE_POTENTIAL_H
