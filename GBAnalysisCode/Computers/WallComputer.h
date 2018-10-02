#ifndef WALLCOMPUTER

#define WALLCOMPUTER

/////////////////////////////////////////////////////////////////////////////////
//Wall Computer class. 
//
//Description
//		This class maintains a spatial partition of particle systems so that
//		only local regions of a given particle need to be searched for neighbors.
//		Implements bonds between particles and the wall. 
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
class CWallComputer : public CBaseComputer<Dim>
{
private:
	typedef Eigen::Matrix<dbl,Dim,1>   dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
	CStaticState<Dim> &State;
	CGrid<Dim> Grid;

	bool FixDof;
	vector<bool> FixedDof;

public:
	CBondList<Dim> Bonds;	//!<Bond list for standard computations
	CBondList<Dim> WallBonds; //!<Bond list between particles and the x-wall
	index_map RattlerMap;	//!<Map for rattlers.
	CStdData<Dim>  Data;	//!<storage for standard data

	CWallComputer();
	CWallComputer(CStaticState<Dim> &_State);
	CWallComputer(CStaticState<Dim> &_State,bool usegrid);
	CWallComputer(const CWallComputer<Dim> &_Copy);
	
	bool UseGrid;	
		
//get and set the state
	CStaticState<Dim> &GetState() const;
	void SetState(CStaticState<Dim> & State);

//Get generic info from State
	dbl GetVolume() const;
	int GetNdof() const;

//Compute the bond list
	void ComputeBondList(CBondList<Dim> &bonds);
	void ComputeBondList_Grid(CBondList<Dim> &bonds);
	void ComputeBondList_NoGrid(CBondList<Dim> &bonds) const;
	void SetUseGrid(bool usegrid=true);


	int  StdPrepareSystem(bool verbose = true);								//!<Compute the bonds list and remove rattlers
	int  StdPrepareSystem(vector<bool> const &fixed, bool verbose = true);	//!<Compute the bonds list and remove rattlers, assuming some fixed particles

//Calculate standard data
	void CalculateStdData(bool CalcCijkl = true, bool CalcHess = true);
	void CalculateStdData(CStdData<Dim> &data, bool CalcCijkl=true, bool CalcHess = true);
	void CalculateStdData_Unstressed(bool CalcCijkl=true, bool CalcHess = true);
	void CalculateStdData_Unstressed(CStdData<Dim> &data, bool CalcCijkl=true, bool CalcHess = true);
	dbl  CalculatePressure();


//Compute Hessian
	void ComputeHessian(dbl unstress_coeff, dbl stress_coeff, dbl tether);
	void ComputeHessian(CBondList<Dim> const &bonds, Eigen::SparseMatrix<dbl> &hess) const;
	void ComputeHessianBZ(Eigen::SparseMatrix<cdbl> &hess, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const;
	void ComputeHessianBZ(CBondList<Dim> const &bonds, Eigen::SparseMatrix<cdbl> &hess, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const;
	
//Stress stuff	
	void ComputeCijkl(dbl unstress_coeff, dbl stress_coeff);
	void PrintModuli();

	
//Needed for minimization routines
	void SetFixedDof(vector<bool> const &_FixedDof);
	void UnsetFixedDof();

	void Evaluate(Eigen::VectorXd &grad, dbl &fx);
	bool Progress(Eigen::VectorXd const &grad, dbl fx, int iteration, int print_iter, dbl tol);
	void Move(Eigen::VectorXd const &step);
	dbl GetMinimizationTimeScale() const;
    

	void CalculateBZTransformation(Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const;
	void CalculateBZTransformation(CBondList<Dim> const &bonds, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const;
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION    //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
CWallComputer<Dim>::CWallComputer() : FixDof(false), UseGrid(true)
{
}


template <int Dim>
CWallComputer<Dim>::CWallComputer(CStaticState<Dim> &_State) : State(_State), Grid(&State), FixDof(false), UseGrid(true)
{	
	
	dvec MaxDistance;
	State.GetMaxDistance(MaxDistance);
	bool enough_cells = true;
	for (int dd = 0; dd< Dim; ++dd)
		if((int) floor(1.0/MaxDistance[dd])<3) enough_cells = false;
	
	if(enough_cells)
	{
		SetUseGrid(true);	
		Grid.Allocate();
	}
	if(!enough_cells)
		SetUseGrid(false);
}

template <int Dim>
CWallComputer<Dim>::CWallComputer(CStaticState<Dim> &_State,bool griduse) : State(_State), Grid(&State), FixDof(false)
{
        SetUseGrid(griduse);
        if(griduse == true) Grid.Allocate();
}



template <int Dim>
CWallComputer<Dim>::CWallComputer(const CWallComputer<Dim> &_Copy) : State(_Copy.State), Grid(&State), FixDof(_Copy.FixDof), FixedDof(_Copy.FixedDof), UseGrid(_Copy.UseGrid)
{
	dvec MaxDistance;
        State.GetMaxDistance(MaxDistance);
        bool enough_cells = true; 
        for (int dd = 0; dd< Dim; ++dd)
                if((int) floor(1.0/MaxDistance[dd])<3) enough_cells = false;

        if(enough_cells)
        {
                SetUseGrid(true);
                Grid.Allocate();
        }
        if(!enough_cells)
                SetUseGrid(false);
}


//get and set the state
template <int Dim>
CStaticState<Dim> &CWallComputer<Dim>::GetState() const
{
	return State;
}

template <int Dim>
void CWallComputer<Dim>::SetState(CStaticState<Dim> & _State)
{
	State = _State;
	Grid.SetState(&State);
}
	
template <int Dim>
inline dbl CWallComputer<Dim>::GetVolume() const
{
	return State.GetVolume();
}

template <int Dim>
inline int CWallComputer<Dim>::GetNdof() const
{
	return State.GetParticleNumber()*Dim;
}

template <int Dim>
inline dbl CWallComputer<Dim>::GetMinimizationTimeScale() const
{
	dbl avgRad = State.GetAverageRadius();
	return avgRad; //This is emperical!!! NOT GENERAL!
}

template <int Dim>
void CWallComputer<Dim>::SetUseGrid(bool usegrid)
{
	UseGrid = usegrid;
}

//Compute the bond list
template <int Dim>
void CWallComputer<Dim>::ComputeBondList(CBondList<Dim> &bonds)
{
	if(UseGrid)
		ComputeBondList_Grid(bonds);
	else
		ComputeBondList_NoGrid(bonds);

	//Compute wallbonds
        WallBonds.RemoveAllBonds();

        WallBonds.SetN(State.GetParticleNumber());
        dvec Displacement;
	dvec wallpos;
	dvec partpos;

	dbl xwall = pow(State.GetVolume(), 1./3.);
	dbl xdist, xdist2;
	
	dbl x1, x2;
	
        dbl sigma, E, g, k;
	for(int i=0; i < State.GetParticleNumber(); i++)
	{
		State.GetParticlePosition(partpos,i);
		x1 = fabs(partpos(0)-0.0);
		x2 = fabs(partpos(0)-xwall); 

		xdist = min(x1,x2);
		xdist2= xdist*xdist;
		for(int dd = 1; dd < Dim; ++dd) Displacement(dd) = 0.0;

		if(State.GetPotential()->Overlapping(State.GetRadius(i),0.5, xdist2, sigma))
		{
			Displacement(0) =  xdist;
			if(x2<x1) Displacement(0)= -xdist;
			State.GetPotential()->ComputeDerivatives012(xdist, sigma, E, g, k);
                        WallBonds.AddBond( CBond<Dim>(i, i, xdist, E, g, k, Displacement) );
		}

	}	


}

template <int Dim>
void CWallComputer<Dim>::ComputeBondList_Grid(CBondList<Dim> &bonds)
{
	//Make sure bonds is empty.
	bonds.RemoveAllBonds();

	Grid.Construct();
	//Grid.PrintGrid();

	bonds.SetN(State.GetParticleNumber());
	dvec Displacement;
	dbl sigma, rlen, rlen2, E, g, k;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		for(typename CGrid<Dim>::iterator it = Grid.GetParticleIterator(i) ; (*it)!=-1 ; it++) //(*it) is the particle index of a potential neighbor
		{
			if((*it)>i)
			{
				State.GetDisplacement(i,(*it),Displacement);
				rlen2 = Displacement.squaredNorm();
				if(State.GetPotential()->Overlapping(State.GetRadius(i),State.GetRadius((*it)), rlen2, sigma))
				{
					rlen = sqrt(rlen2);
					State.GetPotential()->ComputeDerivatives012(rlen, sigma, E, g, k);
					bonds.AddBond( CBond<Dim>(i, (*it), rlen, E, g, k, Displacement) );
				}
			}
		}
	}
	bonds.SetVolume(GetVolume());
}

template <int Dim>
void CWallComputer<Dim>::ComputeBondList_NoGrid(CBondList<Dim> &bonds) const
{
	//Make sure bonds is empty.
	bonds.RemoveAllBonds();

	dvec Displacement;
	bonds.SetN(State.GetParticleNumber());
	dbl sigma, rlen, rlen2, E, g, k;
	for(int i = 0 ; i < State.GetParticleNumber()-1 ; i++)
		for(int j = i+1 ; j < State.GetParticleNumber() ; j++)
		{
			State.GetDisplacement(i,j,Displacement);
			rlen2 = Displacement.squaredNorm();
			if(State.GetPotential()->Overlapping(State.GetRadius(i),State.GetRadius(j), rlen2, sigma))
			{
				rlen = sqrt(rlen2);
				State.GetPotential()->ComputeDerivatives012(rlen, sigma, E, g, k);
				bonds.AddBond( CBond<Dim>(i, j, rlen, E, g, k, Displacement) );
			}
		}
	bonds.SetVolume(GetVolume());
}



/**
 *	@param[in] verbose flag to control output
 *	@return 0 if the bonds list is prepared correctly, 1 if the bonds list is empty
 */
template <int Dim>
int CWallComputer<Dim>::StdPrepareSystem(bool verbose)
{
	std::vector<bool> fixed; fixed.assign(State.GetParticleNumber(),false);
	return StdPrepareSystem(fixed,verbose);
}

/**
 *	@param[in] fixed vector of bools that labels particles as being fixed
 *	@param[in] verbose flag to control output
 *	@return 0 if the bonds list is prepared correctly, 1 if the bonds list is empty
 */
template <int Dim>
int CWallComputer<Dim>::StdPrepareSystem(vector<bool> const &fixed, bool verbose)
{
	ComputeBondList(Bonds);										//Use the member variable derived from CBaseComputer
	Bonds.RemoveRattlers(RattlerMap, fixed, Dim+1, verbose);	//Remove rattlers
	return (Bonds.Empty())?1:0;
}

template <int Dim>
void CWallComputer<Dim>::CalculateStdData(bool CalcCijkl, bool CalcHess)
{
	CBondList<Dim> bonds;
        ComputeBondList(bonds);

        dbl p = CalculatePressure();
        Data.Pressure = p;
        Data.NPp = State.GetParticleNumber();
        Data.Volume = bonds.GetVolume();
        Data.Nc = (int) bonds.GetNBonds()+2*WallBonds.GetNBonds();

        Eigen::VectorXd grad,gradw;
        bonds.ComputeGradient(grad);
        WallBonds.ComputeGradient(gradw);
        Data.MaxGrad = max(max_abs_element(grad),max_abs_element(gradw));

        if(CalcHess)
        {
                ComputeHessian(1.0,1.0,1e-12);
        };

        if(CalcCijkl)
        {
                ComputeCijkl(1.0,1.0);
        };
}
template <int Dim>
void CWallComputer<Dim>::CalculateStdData(CStdData<Dim> &data, bool CalcCijkl, bool CalcHess)
{
	CBondList<Dim> bonds;
        ComputeBondList(bonds);

	dbl p = CalculatePressure();
	Data.Pressure = p;
	Data.NPp = State.GetParticleNumber();
	Data.Volume = bonds.GetVolume();
	Data.Nc = (int) bonds.GetNBonds()+2*WallBonds.GetNBonds();

	Eigen::VectorXd grad,gradw;
        bonds.ComputeGradient(grad);
	WallBonds.ComputeGradient(gradw);
	Data.MaxGrad = max(max_abs_element(grad),max_abs_element(gradw));
	
	if(CalcHess)
	{
		ComputeHessian(1.0,1.0,1e-12);
	};

	if(CalcCijkl)
	{
		ComputeCijkl(1.0,1.0);
	};

	data = Data;
};

template <int Dim>
void CWallComputer<Dim>::PrintModuli()
{
	CalculateStdData(true,true);
	dbl bulk = Data.cijkl.CalculateBulkModulus();
	dbl shear = Data.cijkl.CalculateShearModulus();
	dbl weakest = Data.cijkl.calc_weakest_deformation();

	cout << "Bulk modulus = " << bulk << endl;
	cout << "shear modulus = " << shear << endl;
	cout << "weakest response to any strain = " << weakest << endl;

}


template <int Dim>
void CWallComputer<Dim>::CalculateStdData_Unstressed(bool CalcCijkl, bool CalcHess)
{
	CalculateStdData_Unstressed(Data,CalcCijkl,CalcHess);
}


template <int Dim>
void CWallComputer<Dim>::CalculateStdData_Unstressed(CStdData<Dim> &data, bool CalcCijkl, bool CalcHess)
{
	CBondList<Dim> BondsTemp = Bonds;
	BondsTemp.MakeUnstressed();
	BondsTemp.CalculateStdData(data,CalcCijkl,CalcHess);
}











//Needed for minimization routines

template <int Dim>
void CWallComputer<Dim>::SetFixedDof(vector<bool> const &_FixedDof)
{
	FixDof = true;
	FixedDof = _FixedDof;
	assert(FixedDof.size() == Dim*State.GetParticleNumber());
}

template <int Dim>
void CWallComputer<Dim>::UnsetFixedDof()
{
	FixDof = false;
}

template <int Dim>
void CWallComputer<Dim>::Evaluate(Eigen::VectorXd &grad, dbl &fx) 
{
	CBondList<Dim> bonds;
	ComputeBondList(bonds);
	fx = FixDof?bonds.ComputeGradient(grad,FixedDof):bonds.ComputeGradient(grad);

	//now evaluate the wall bonds
        dvec temp;
        CBond<Dim> b;
	int bondnum = WallBonds.GetNBonds();
	for(int bb = 0; bb<bondnum; ++bb)
        {
                WallBonds.GetBond(bb,b);
                
		fx += b.E;
		temp = b.g/b.rlen*b.r;
                for(int dd=0; dd<Dim; ++dd)
                {
                        grad[Dim*b.i+dd] += temp[dd];
                };
        };	


};

template <int Dim>
bool CWallComputer<Dim>::Progress(Eigen::VectorXd const &grad, dbl fx, int iteration, int print_iter, dbl tol) 
{
//	const int print_iter = 1000;
	dbl gradNorm = grad.norm();
	dbl max_grad = max_abs_element(grad.size(), grad.data());
	bool converged = (max_grad < tol)?true:false;
	if(converged){ printf("converged\n"); fflush(stdout);}
	if(print_iter > 0)
	{
		if(iteration%print_iter==0 || converged)
		{
			printf("% 11i   %22.20e   % e   % e\n", iteration, fx, gradNorm, max_grad);
			fflush(stdout);
		}
	}
	return converged;
};

template <int Dim>
inline void CWallComputer<Dim>::Move(Eigen::VectorXd const &step) 
{
	State.MoveParticles(step);
};
  





template <int Dim>
dbl CWallComputer<Dim>::CalculatePressure()
{
	CBondList<Dim> bonds;
        ComputeBondList(bonds);
	dbl vol = bonds.GetVolume();
	dmat stress;
	stress.setZero();

	CBond<Dim> b;
        int bondnum = bonds.GetNBonds();
        for(int bb = 0; bb<bondnum; ++bb)
        {
                bonds.GetBond(bb,b);
		stress += b.g*((b.r)*(b.r.transpose()))/b.rlen;
        };
	cout << "internal bond pressure:" << -stress.trace()/(((dbl)Dim)*vol) <<endl;


        bondnum = WallBonds.GetNBonds();
        for(int bb = 0; bb<bondnum; ++bb)
        {
        	WallBonds.GetBond(bb,b);
		stress += b.g*((b.r)*(b.r.transpose()))/b.rlen;
	};

	stress /= vol;
	Data.Stress=stress;	
	cout << "total bond pressure:" << -stress.trace()/((dbl)Dim) <<endl;

	return -stress.trace()/((dbl)Dim);
}

template <int Dim>
void CWallComputer<Dim>::ComputeCijkl(dbl unstress_coeff, dbl stress_coeff)
{
	CBondList<Dim> bonds;
        ComputeBondList(bonds);	

	int Nvar = Dim*State.GetParticleNumber();
	Data.H.LUdecomp(); //requires Calculate Hessian to have already been called. Not ideal.
	dmat strain_tensor;
        dbl prefactor = 2./bonds.GetVolume();//The prefactor contains a 2 in it from the equation dE/V = (1/2)cijkl uij ukl.

	Eigen::VectorXd n_d2Udvdgamma = Eigen::VectorXd::Zero(Nvar);
        Eigen::VectorXd uNonAffine_node = Eigen::VectorXd::Zero(Nvar);
        Eigen::VectorXd DeltaR_bond = Eigen::VectorXd::Zero(Dim*(bonds.GetNBonds()+WallBonds.GetNBonds()));


	for(int ii=0; ii<cCIJKL<Dim>::num_constants; ++ii)
        {
                //Set the strain tensor
                Data.cijkl.set_strain_tensor(strain_tensor, ii);


		int i,j;
        	dbl g, k, DeltaRparallel;
        	dvec D2i;
		n_d2Udvdgamma = Eigen::VectorXd::Zero(Nvar);

		//Calculate the displacement vector for each bond in the metric defined by the strain tensor.
                //Also, calculate the forces on each particle due to the change in metric.
		CBond<Dim> b;
        	int bondnum = bonds.GetNBonds();
        	for(int bb = 0; bb<bondnum; ++bb)
        	{
			bonds.GetBond(bb,b);
			g=b.g*stress_coeff;
			k=b.k*unstress_coeff;
			//Calculate the affine extension of the bond.
                	DeltaR_bond.segment<Dim>(Dim*bb) = strain_tensor*b.r;

                	//Project that onto the direction of b.r;
                	DeltaRparallel = DeltaR_bond.segment<Dim>(Dim*bb).dot(b.r)/b.rlen;

                	//Calculate the induced forces on particles i and j.
                	D2i = (DeltaRparallel/b.rlen) * (k - g/b.rlen) * b.r  +  (g/b.rlen) * DeltaR_bond.segment<Dim>(Dim*bb);

	                //Add the induced force to the n_d2Udvdgamma vector
        	        n_d2Udvdgamma.segment<Dim>(Dim*b.i) += D2i;
                	n_d2Udvdgamma.segment<Dim>(Dim*b.j) -= D2i;
		};
		int wallbondnum = WallBonds.GetNBonds();
                for(int bb = 0; bb<wallbondnum; ++bb)
                {
                        WallBonds.GetBond(bb,b);
                        g=b.g*stress_coeff;
                        k=b.k*unstress_coeff;
                        DeltaR_bond.segment<Dim>(Dim*(bb+bondnum)) = strain_tensor*b.r;
                        DeltaRparallel = DeltaR_bond.segment<Dim>(Dim*(bb+bondnum)).dot(b.r)/b.rlen;
                        //Calculate the induced forces on particle i
                        D2i = (DeltaRparallel/b.rlen) * (k - g/b.rlen) * b.r  +  (g/b.rlen) * DeltaR_bond.segment<Dim>(Dim*(bb+bondnum));
                        n_d2Udvdgamma.segment<Dim>(Dim*b.i) += D2i;
                };

                //Solve for the non-affine displacement
                Data.H.solve_Mx_equals_b(uNonAffine_node, n_d2Udvdgamma);

                //Add the non-affine extension of each bond to the affine extension
		for(int bb = 0; bb<bondnum; ++bb)
		{
			bonds.GetBond(bb,b);
			DeltaR_bond.segment<Dim>(Dim*bb) += uNonAffine_node.segment<Dim>(Dim*b.j) - uNonAffine_node.segment<Dim>(Dim*b.i);
		}
		for(int bb = 0; bb<wallbondnum; ++bb)
                {
                        WallBonds.GetBond(bb,b);
                        DeltaR_bond.segment<Dim>(Dim*(bb+bondnum)) += - uNonAffine_node.segment<Dim>(Dim*b.i);
                }


                //Calculate the change in energy
		dbl dE = 0;
		dbl DeltaRparallel2, DeltaRperp2;
		for(int bb = 0; bb<bondnum; ++bb)
                {
                        bonds.GetBond(bb,b);
			g = b.g*stress_coeff;
                	k = b.k*unstress_coeff;
                	DeltaRparallel2 = POW2(DeltaR_bond.segment<Dim>(Dim*bb).dot(b.r)/b.rlen);
                	DeltaRperp2     = DeltaR_bond.segment<Dim>(Dim*bb).squaredNorm() - DeltaRparallel2;
			dE += k*DeltaRparallel2 + g*DeltaRperp2/b.rlen;
                };
                for(int bb = 0; bb<wallbondnum; ++bb)
                {
                        WallBonds.GetBond(bb,b);
			g = b.g*stress_coeff;
                        k = b.k*unstress_coeff;
                        DeltaRparallel2 = POW2(DeltaR_bond.segment<Dim>(Dim*(bb+bondnum)).dot(b.r)/b.rlen);
                        DeltaRperp2     = DeltaR_bond.segment<Dim>(Dim*(bb+bondnum)).squaredNorm() - DeltaRparallel2;
                        dE += k*DeltaRparallel2 + g*DeltaRperp2/b.rlen;                	
                };
		dE = dE *0.5;
		dE *= prefactor; //dE is now 2*dE/V. This comes from the equation dE/V = (1/2)cijkl uij ukl.

                //Set the ii'th elastic constant
                Data.cijkl.set_constant(dE, ii);
        };


};      

template <int Dim>
void CWallComputer<Dim>::ComputeHessian(dbl unstress_coeff, dbl stress_coeff, dbl tether)
{
	typedef Eigen::Triplet< dbl>  TRIP;
        vector<TRIP> coeffs;
	//first compute the hessian blocks from the particle-particle bonds
        CBondList<Dim> bonds;
        ComputeBondList(bonds);
        bonds.ComputeHessianElements(coeffs, unstress_coeff, stress_coeff, tether);

	//Next compute the hessian blocks from the particle-wall bonds
	dmat Fii, Fjj, Fij; //Fji = Fij.transpoze(); stressed block
        dmat Kii, Kjj, Kij; //Kji = Kij.transpose(); unstressed block
        dmat Bii, Bjj, Bij; //Bji = Bij.transpose(); B = unstress_coeff*K + stress_coeff*F
        int icorner;
	CBond<Dim> b;
	int bondnum = WallBonds.GetNBonds();
        for(int bb = 0; bb<bondnum; ++bb)
        {
		WallBonds.GetBond(bb,b);
		b.CalculateMatrixBlocks(Fij,Kij);
		Fii = -Fij;
                Kii = -Kij;

                Bii = unstress_coeff*Kii + stress_coeff*Fii;

                //Add the matrix blocks to the list of coefficients
                icorner = Dim*b.i;
		for(int d1=0; d1<Dim; ++d1)
                        for(int d2=0; d2<Dim; ++d2)
                        {
                                coeffs.push_back( TRIP(icorner+d1, icorner+d2, Bii(d1,d2)) );
			};

	};


        //Create a temporary matrix, and set the hessian from the acquired coefficients
	Eigen::SparseMatrix<dbl> hess;
        Eigen::SparseMatrix<dbl> temp(Dim*State.GetParticleNumber(),Dim*State.GetParticleNumber());
        hess = temp;
        hess.setFromTriplets(coeffs.begin(), coeffs.end());
	Data.H.A = hess;
};


template <int Dim>
void CWallComputer<Dim>::ComputeHessian(CBondList<Dim> const &bonds, Eigen::SparseMatrix<dbl> &hess) const
{
	bonds.ComputeHessian(hess);
}

template <int Dim>
void CWallComputer<Dim>::ComputeHessianBZ(Eigen::SparseMatrix<cdbl> &hess, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const
{
	ComputeHessianBZ(Bonds, hess, transformation, k);
}
template <int Dim>
void CWallComputer<Dim>::ComputeHessianBZ(CBondList<Dim> const &bonds, Eigen::SparseMatrix<cdbl> &hess, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const
{
	CalculateBZTransformation(bonds, transformation, k);
	Eigen::SparseMatrix<dbl> H0;
	ComputeHessian(bonds, H0);
	printf("transformation -> %i %i\n", transformation.rows(), transformation.cols());
	printf("H0             -> %i %i\n", H0.rows(), H0.cols());
	hess = transformation.adjoint() * H0 * transformation;
}

template <int Dim>
void CWallComputer<Dim>::CalculateBZTransformation(Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const
{
	CalculateBZTransformation(Bonds, transformation, k);
}

template <int Dim>
void CWallComputer<Dim>::CalculateBZTransformation(CBondList<Dim> const &bonds, Eigen::SparseMatrix<cdbl> &transformation, dvec const &k) const
{
	//Get the real positions
	Eigen::VectorXd RealPos;
	State.GetPositions(RealPos);

	//Remove Positions corresponding to rattlers
	Eigen::VectorXd RealPosNoRatt = Eigen::VectorXd::Zero(Dim*RattlerMap.size());
	for(int im=0; im<RattlerMap.size(); ++im)
		RealPosNoRatt.segment<Dim>(Dim*im) = RealPos.segment<Dim>(Dim*RattlerMap[im]);
	bonds.PlaneWaveAnsatz(transformation, RealPosNoRatt, k);
}



}

#endif
