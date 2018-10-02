#ifndef BASECOMPUTER

#define BASECOMPUTER

/////////////////////////////////////////////////////////////////////////////////
//Grid class. 
//
//Description
//		This class maintains a spatial partition of particle systems so that
//		only local regions of a given particle need to be searched for neighbors.
//		The f 
//
//Variables
// 		Particle positions as an eigen vector
//		Particle radii as an eigen vector
//		Number of particles as an integer
//		Box as a CBox class
//		Potential as a CPotential class
//		
//Implements
//		Reading/Writing to netcdf files.
//		Getting/Setting particle positions.
//						particle radii.
//						box.
//						potential.
//
//
//File Formats
//		NetCDF
//			## Note: only files of a single N and Dimension may be stored in a given 
//			## file. To denote whether the NetCDF file has been populated with variables,
//			## dimensions, etc... an attribute "System_Populated" gets created.
//			-> Three dimensions: record number (unlimited), degrees of freedom (DOF), and
//			   particle number.
//			-> Dimension is stored as an attribute
//			-> Particle positions are stored as a variable.
//			-> Particle radii are stored as a variable.
//			-> Box is stored as a box (see CBox definition.)
//			-> Potential is stored as a potential (see CPotential definition.)
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include "../Potentials/Potentials.h"
#include "../Boundaries/Boxes.h"
#include "../Resources/MersenneTwister.h"
#include "BondList.h"
#include <list>


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


//!Abstract base class for Computers.

template <int Dim>
class CBaseComputer
{
public:
//! @name Get Info
///@{
	virtual dbl GetVolume() const = 0;		//!<Get the volume
	virtual int GetNdof() const = 0;		//!<Get the number of degrees of freedom
//	virtual dbl GetTimeUnit() const = 0;	//!<Get the unit of time
//	virtual dbl GetLengthUnit() const = 0;	//!<Get the unit of length
//	virtual dbl GetEnergyUnit() const = 0;	//!<Get the unit of energy
//	virtual dbl GetMassUnit() const = 0;	//!<Get the unit of mass

///@}




//! @name Stataic Computations
///@{

/*
//Compute the energy of the system
    virtual double ComputeEnergy() = 0;
    
//compute the gradient of the energy of the system.
    virtual void ComputeGradient(Eigen::VectorXd &tar) = 0;
    virtual void ComputeForce(Eigen::VectorXd &tar) = 0;
    
//Hessian stuff
    virtual void ComputeHessian(Eigen::MatrixXd &tar) = 0;

//Dynamical Matrix Stuff (all done in mass-normalized coordinates)
//computes at q = 0
    virtual void ComputeDynamicalMatrix(Eigen::MatrixXd &tar) = 0;
*/
///@}

//! @name Methods for Minimization
///@{
	virtual void Evaluate(Eigen::VectorXd &grad, dbl &fx) = 0;
	virtual bool Progress(Eigen::VectorXd const &grad, dbl fx, int iteration, int print_iter, dbl tol) = 0;
	virtual void Move(Eigen::VectorXd const &step) = 0;
	virtual void ReportHeader()
	{
		printf("         ii                       Energy            norm        max_grad \n");
	};
	virtual dbl GetMinimizationTimeScale() const = 0;

///@}
   
};

}

#endif
