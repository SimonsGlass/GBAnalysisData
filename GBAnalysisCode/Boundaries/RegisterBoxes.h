#ifndef REGISTERBOXES

#define REGISTERBOXES

/////////////////////////////////////////////////////////////////////////////////
//Potential factory class. 
//
//Description
//		This is a class that loads different potential types from files.  
//
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

#include "Box.h"
#include "PeriodicBox.h"

namespace LiuJamming
{

template<int Dim>
void RegisterBoxes()
{
	CBox<Dim>::AddBoxType("PeriodicBox",new CPeriodicBox<Dim>());
	CBox<Dim>::AddBoxType("PeriodicBox1",new CPeriodicBox<Dim,1>());
	
}

}

#endif
