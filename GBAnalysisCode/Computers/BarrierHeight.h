#ifndef BARRIERHEIGHT

#define BARRIERHEIGHT

/////////////////////////////////////////////////////////////////////////////////
//Barrier Height Class. 
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
#include "../State/StaticState.h"
#include "../Minimization/minimizer.h"
#include "BaseComputer.h"
#include "BondList.h"
#include "Grid.h"
#include "Data.h"
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

template <int Dim>
class CBarrierHeight
{
private:
	CStaticComputer<Dim> &Computer;

	double *Barriers;
	double *Displacements;
public:
	CBarrierHeight(CStaticComputer<Dim> &);
	CBarrierHeight(CStaticState<Dim> &);
	~CBarrierHeight();

	//Compute the actual barrier heights for a collection of N polarization vectors
	//The displacements are all found so that the barrier is within \pm tolerance.
	void Compute(double *,int,double);
	
	const double *GetBarriers() {return Barriers;}
	const double *GetDisplacements() {return Displacements;}	
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION    //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
CBarrierHeight<Dim>::CBarrierHeight(CStaticComputer<Dim> &c) : Computer(c)
{
	Barriers = NULL;
	Displacements = NULL;
}

template <int Dim>
CBarrierHeight<Dim>::CBarrierHeight(CStaticState<Dim> &s) : Computer(CStaticComputer<Dim>(s))
{
	Barriers = NULL;
	Displacements = NULL;
}


template <int Dim>
CBarrierHeight<Dim>::~CBarrierHeight()
{
	if(Barriers)
		delete[] Barriers;

	if(Displacements)
		delete[] Displacements;
}

template <int Dim>
void CBarrierHeight<Dim>::Compute(double *modes,int Number,double Tolerance)
{
	if(Barriers)
		delete[] Barriers;
	
	if(Displacements)
		delete[] Displacements;
	
	Barriers = new double[Number];
	Displacements = new double[Number];

	//Get the number of degrees of freedom
	int DOF = Computer.GetNdof();

	//find the energy of the initial configuration (to compare with the energy at the transition point)
	Computer.StdPrepareSystem(false);
        Computer.CalculateStdData(false,false);
        dbl initialE = Computer.Data.Energy;

	cout << "Computing barrier heights for " << Number << " modes.\n";

	//Go through each mode
	for(int i = 0 ; i < Number ; i++)
	{
		double Height = 1e10;
		double Displace = 0.0;

		cout << "Evaluating mode " << i << endl;

		//push in both the forward and backward direction
		for(double dir = -1 ; dir<2 ; dir+=2)
		{

			//cout << "Considering the " << (dir<0 ? "forward " : "backward ") << "direction.\n";
			double T_Height = 0.0;
			double T_Displacement = 1.0;

			//Get the polarization vector for the current mode
			Eigen::VectorXd Mode = Eigen::Map<Eigen::VectorXd>(modes+DOF*i,DOF);
			
			//cout << "MODE: \n" << Mode.transpose() << "\n-----\n";

			//Compute the maximum magnitude and use this to normalize the mode. (This is the most suspect part of the process
			//I wonder if there's a better way to choose the normalization.)
			double Max_Pol = 0.0;
			for(int j = 0 ; j < DOF/Dim ; j++)
			{
				double Pol = 0.0;
				for(int k = 0 ; k < Dim ; k++)
					Pol+=Mode(Dim*j+k)*Mode(Dim*j+k);

				if(Pol>Max_Pol)
					Max_Pol = Pol;
			}
			//cout << Mode.transpose() << endl; 			
			Max_Pol = sqrt(Max_Pol);
			dbl init_T_disp = 0.1/Max_Pol;


			//We assume that the particle radii are centered around 1.0
			T_Displacement = init_T_disp;
			double Delta = T_Displacement/2.0;

			//Record the initial configuration
			Eigen::VectorXd Initial_Position = Eigen::ArrayXd::Zero(DOF);
			Eigen::VectorXd Displacement = Eigen::ArrayXd::Zero(DOF);
			Eigen::VectorXd Current_Position = Eigen::ArrayXd::Zero(DOF);
			Computer.GetState().GetPositionsVirtual(Initial_Position);
			Computer.StdPrepareSystem(false);
			Computer.CalculateStdData(false,false);
			int iteration = 0;
			bool AtLeastOnceStable = false;

			double Norm = 0.0;

			while(T_Displacement==0||Delta/T_Displacement>Tolerance)
			{	
				//Move the system along the mode.
				Computer.Move(dir*Mode*T_Displacement);
				Computer.ComputeBondList_NoGrid(Computer.Bonds);
				Computer.CalculateStdData(false,false);
				T_Height = Computer.Data.Energy;

				//Reminimize
				CSimpleMinimizer<Dim> Min(Computer,CSimpleMinimizer<Dim>::FIRE);

				//Compute Displacement
				Computer.GetState().GetPositionsVirtual(Current_Position);
				Computer.GetState().GetBox()->MinimumDisplacement(Initial_Position,Current_Position,Displacement);

				//Subtract off any constant displacement
				Eigen::Matrix<double,Dim,1> CM_Displacement;
				CM_Displacement(0) = 0;
				CM_Displacement(1) = 0;
				for(int j = 0 ; j < DOF/Dim ; j++)
					CM_Displacement += Displacement.segment<Dim>(Dim*j)/DOF*Dim;


				for(int j = 0 ; j < DOF/Dim ; j++)
					Displacement.segment<Dim>(Dim*j) -= CM_Displacement;

				//if the residual displacement is less than some small amount we can increase shift
				//Otherwise reduce the delta if we've found a baseline T_Displacement.

				Norm = Displacement.norm();

				cout << "Finishing iteration " << iteration << " with Displacement = " << T_Displacement << ", Delta = " << Delta << " and Height = " << T_Height << " and Norm = " << Norm << endl;

				if(Norm < 5e-4){
					T_Displacement += Delta;
					Delta=min(Delta*1.25,init_T_disp*2.5);
					AtLeastOnceStable = true;
				}else {
					if(!AtLeastOnceStable && T_Displacement <= Delta) Delta/=2.0;

					T_Displacement -= Delta;
					if(AtLeastOnceStable){
						Delta /= 2.0;
						T_Displacement+=Delta;
					}
				}
				
				Computer.GetState().SetPositionsVirtual(Initial_Position);
				iteration++;

			}

			if(T_Height < Height)
			{
				Height = T_Height;
				Displace = T_Displacement;
			}
		}
		Barriers[i]  = Height-initialE;
		Displacements[i] = Displace;
	}
}



}

#endif
