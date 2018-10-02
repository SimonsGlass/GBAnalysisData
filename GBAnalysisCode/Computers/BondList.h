#ifndef BOND_LIST_H
#define BOND_LIST_H

#include "../Resources/std_include.h"
#include "../Resources/index_map.h"
#include <Eigen/Sparse>
#include "Bond.h"
#include "Neighbors.h"
#include "MatrixInterface.h"
#include "cijkl.h"
#include "Data.h"

namespace LiuJamming
{

//! Class to store a list of bonds.

/*!
 * This class stores a list of bonds and provides numerous functions for making static calculations,
 * such as computing the energy, gradient and the hessian matrix. It also provides functions for
 * removing rattlers.
 */
template<int Dim>
class CBondList
{
	typedef Eigen::Matrix<dbl, Dim, 1> dvec;
	typedef Eigen::Matrix<dbl, Dim, Dim> dmat;
	typedef CBond<Dim> BOND;
	typedef Eigen::Triplet< dbl>  TRIP;
	typedef Eigen::Triplet<cdbl> cTRIP;

//! @name Storage Variables
///@{
	int N;						//!<Number of possible nodes.
	dbl Volume;					//!<Volume containing the bonds.
	vector<BOND> list;			//!<A list of the bonds.

///@}

public:
//! @name Constructors and Operators
///@{
	CBondList(int _N=0, dbl V=0.);							//!<Construct a CBondsList.
	CBondList(const CBondList &src);						//!<Copy constructor.
	CBondList<Dim>& operator=(const CBondList<Dim> &src);	//!<Overloaded equals operator.

///@}

//! @name Get and Set Methods
///@{
	void AddBond(const BOND &b);	//!<Add a bond to the list
	void SetVolume(dbl _V);			//!<Set the volume.
	void SetN(int _N);				//!<Set the number of nodes.
	void RemoveAllBonds();			//!<Delete all the bonds and set N and Volume to zero.

	int GetN() const;				//!<Get the number of nodes.
	int GetNBonds() const;			//!<Get the number of bonds.
	dbl GetVolume() const;			//!<Get the volume.

	void GetBond(int i, BOND &b) const; //!<Copy a bond 

	bool Empty() const;				//!<Returns true if the bond list is empty.

	typename vector<BOND>::const_iterator begin() const;
	typename vector<BOND>::const_iterator end() const;

///@}

//! @name Bond Manipulation
///@{
	void CalculateNeighbors(vector< vector<int> > &nbrs) const;		//!<Get a list of the neighbors of each particle
	void CalculateNeighbors(vector< vector<int> > &nbrs,dbl cutoff) const;             //!<Get a list of the neighbors of each particle within cutoff

	void RemoveBonds(vector<bool> const &BondsToRemove);			//!<Remove bonds from the list
	void RemoveBondsAccordingToMap(index_map const &map);			//!<Remove any bonds that involve a node that is removed from the index_map.
	void UpdateBondIndices(index_map const &map);					//!<When some nodes are removed, as expressed by the index_map map, decrease the i and j indices of all bonds accordingly.
	int  IdentifyRattlers(vector< vector<int> > &nbrs, vector<bool> &rattlers, vector<bool> const &fixed, int c=Dim+1, bool Verbose=false) const; //!<Identify nodes that are not fixed and are involved in less than c bonds.
	int  IdentifyRattlers(vector< vector<int> > &nbrs, vector<bool> &rattlers, int c=Dim+1, bool Verbose=false) const; //!<Identify rattlers assuming no fixed nodes.
	void RemoveRattlers(index_map &map, vector<bool> const &fixed, int c=Dim+1, bool Verbose=false);	//!<Remove rattlers (i.e. nodes with less than c bonds), and return the corresponding index_map.
	void RemoveRattlers(vector<bool> const &fixed, int c=Dim+1, bool Verbose=false);					//!<Overloaded method
	void RemoveRattlers(index_map &map, int c=Dim+1, bool Verbose=false);								//!<Overloaded method
	void RemoveRattlers(int c=Dim+1, bool Verbose=false);												//!<Overloaded method

	void MakeUnstressed();				//!<Remove the prestress form every bond.
	void MultiplyForces(dbl m);			//!<Multiply all forces by a constant.
	void MultiplyStiffnesses(dbl m);	//!<Multiply all spring constants by a constant.

///@}

//! @name Computations
///@{
	dbl  ComputeEnergy() const;							//!<Compute the energy.
	dbl  ComputePressure() const;						//!<Compute the pressure.
	dbl  ComputePressure(dmat &stress) const;			//!<Compute the pressure.
	void ComputeStressTensor(dmat &stress) const;		//!<Compute the Dim by Dim stress tensor.
	dbl  ComputeGradient(Eigen::VectorXd &grad) const;	//!<Compute the Dim*N dimensional energy gradient (i.e. -Fnet), and return the energy.
	dbl  ComputeGradient(Eigen::VectorXd &grad, vector<bool> const &FixedDOF) const;

	void ComputeHessianElements(vector<TRIP> &coefficients, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.) const;				//!<Compute the elements of the hessian as a list.
	void ComputeHessianElements(vector<TRIP> &coefficients,vector<bool> const &FixedDOF) const;
	void PlaneWaveAnsatz(Eigen::SparseMatrix<cdbl> &transformation, Eigen::VectorXd const &Positions, dvec const &k) const;
	//void ComputeHessianElements(vector<cTRIP> &coefficients, dvec k, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.) const;	//!<Compute the elements of the hessian as a list.
	void ComputeHessian(Eigen::SparseMatrix<dbl> &hess, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.)  const;					//!<Compute the hessian, note: the hessian is NOT mass-normalized.
	void ComputeHessian(Eigen::SparseMatrix<dbl> &hess,vector<bool> const &FixedDOF)  const; 
	void ComputeHessian_BZ(Eigen::SparseMatrix<cdbl> &hess, dvec k, dbl unstress_coeff=1., dbl stress_coeff=1., dbl tether=0.) const;		//!<Compute the hessian at non-zero wavevector k, note: the hessian is NOT mass-normalized.

	//void ComputeEquilibriumMatrixElements(vector<cTRIP> &coefficients, dvec k) const;
	void ComputeEquilibriumMatrix   (Eigen::SparseMatrix<dbl> &Amatrix) const;
	//void ComputeEquilibriumMatrix_BZ(Eigen::SparseMatrix<dbl> &Amatrix, dvec k) const;

	void CalculateCijkl(cCIJKL<Dim> &cijkl, dbl unstress_coeff=1., dbl stress_coeff=1., bool Verbose=true) const;	//!<Compute the elastic constants.

	void CalculateStdData(CStdData<Dim> &Data, bool CalcCijkl=true, bool CalcHess=true) const;
///@}

private:
	void Calculate_n_d2Udvdgamma(dmat const &strain_tensor, Eigen::VectorXd &n_d2Udvdgamma, Eigen::VectorXd &AffineDeltaR, dbl unstress_coeff, dbl stress_coeff) const;
	dbl  CalculateEnergyChange(Eigen::VectorXd const &DeltaR_bond, dbl unstress_coeff, dbl stress_coeff) const;

public:

//! @name Misc.
///@{
	bool CheckConsistency() const;				//!<Check that the i and j index of every bond is less than N and greater than or equal to 0.
	void PrintBonds() const;					//!<Print the bond list to stdout.
	void PrintBonds_txt(char *filename) const;	//!<Print the bond list to a text file

///@}
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

		

template<int Dim>
CBondList<Dim>::CBondList(int _N, dbl V)
	:N(_N), Volume(V) 
{};

template<int Dim>
CBondList<Dim>::CBondList(const CBondList &src) 
{
	(*this) = src; 
}

template<int Dim>
CBondList<Dim> &CBondList<Dim>::operator=(const CBondList<Dim> &src)
{
	if(this != &src)
	{
		N = src.N;
		Volume = src.Volume;
		list = src.list;
	}
	return *this;
}

/**
 *	N is automatically updated so that N>b.i and N>b.j. Asserts that b.i>=0 and b.j >= 0.
 */
template<int Dim>
inline void CBondList<Dim>::AddBond(const BOND &b)
{
	assert(b.i>=0);
	assert(b.j>=0);
	if(b.i>=N) N=b.i+1;
	if(b.j>=N) N=b.j+1;
	list.push_back(b);
}

template<int Dim>
void CBondList<Dim>::SetVolume(dbl _V)
{
	Volume = _V;
}

template<int Dim>
void CBondList<Dim>::SetN(int _N)
{
	N = _N;
}

template<int Dim>
int CBondList<Dim>::GetN() const
{
	return N;
}

template<int Dim>
void CBondList<Dim>::RemoveAllBonds()
{
	list.clear();
	N=0;
	Volume=0;
}

template<int Dim>
int CBondList<Dim>::GetNBonds() const
{
	return (int)list.size();
}

template<int Dim>
dbl CBondList<Dim>::GetVolume() const
{
	return Volume;
}

template<int Dim>
void CBondList<Dim>::GetBond(int i, BOND &b) const
{
	//assert(i>=0);
	//assert(i<GetNBonds());
	b = list[i];
}

template<int Dim>
bool CBondList<Dim>::Empty() const
{
	return list.empty();
}
	
template<int Dim>
typename vector< CBond<Dim> >::const_iterator CBondList<Dim>::begin() const
{
	return list.begin();
}

template<int Dim>
typename vector< CBond<Dim> >::const_iterator CBondList<Dim>::end() const
{
	return list.end();
}

template<int Dim>
void CBondList<Dim>::CalculateNeighbors(vector< vector<int> > &nbrs) const
{
	nbrs.assign(N, vector<int>());
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
	{
		nbrs[b->i].push_back(b->j);
		nbrs[b->j].push_back(b->i);
	}
}

template<int Dim>
void CBondList<Dim>::CalculateNeighbors(vector< vector<int> > &nbrs,dbl cutoff) const
{
        nbrs.assign(N, vector<int>());
        for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
        {
		if(b->rlen < cutoff)
		{
                	nbrs[b->i].push_back(b->j);
                	nbrs[b->j].push_back(b->i);
		};
        }
}

/**
 *	Bonds that are not removed are not altered in any way, and N does not change.
 *
 *	@param[in] BondsToRemove A constant vector<bool> of length list.size() that determines
 *	which bonds to remove. Bond i is removed only if BondsToRemove[i]==true.
 */
template<int Dim>
void CBondList<Dim>::RemoveBonds(vector<bool> const &BondsToRemove)
{
	assert(BondsToRemove.size() == list.size());
	vector<BOND> tlist;
	tlist.reserve(list.size());

	for(int i=0; i<(int)list.size(); ++i)
		if(!BondsToRemove[i])
			tlist.push_back(list[i]);

	tlist.swap(list);
}
	
template<int Dim>
void CBondList<Dim>::UpdateBondIndices(index_map const &map)
{
	assert(map.full_size == N);
	for(typename vector<BOND>::iterator b=list.begin(); b!=list.end(); ++b)
	{
		//check that the bond does not involve a rattler
		assert(map.inv(b->i) != -1);
		assert(map.inv(b->j) != -1);

		//change the i and j index according to the map
		b->i = map.inv(b->i);
		b->j = map.inv(b->j);
	}
	N = map.size();
}

/**
 *	@param[in,out] nbrs A vector of vectors of ints that gives the neighbors of each node. This is updated so that rattlers are removed.
 *	@param[in,out] rattlers A vector of bools such that node i is a rattler if rattlers[i]==true. Nodes can be manually designated as a rattler by setting 
 *	rattlers[i]=true before calling this function (normally, all elements of rattlers should be initialied to false). On return, this identifies the nodes that
 *	are rattlers.
 *	@param[in] fixed A vector of bools that identifies nodes as being fixed, and thus cannot become a rattler.
 *	@param[in] c The number of bonds involving a node required for that node to not be a rattler. Default is Dim+1.
 *	@param[in] Verbose Bool to determine if number of identified rattlers, etc. should be printed to stdout.
 *	@return The total number of rattlers.
 */
template<int Dim>
int CBondList<Dim>::IdentifyRattlers(vector< vector<int> > &nbrs, vector<bool> &rattlers, vector<bool> const &fixed, int c, bool Verbose) const
{
	assert(rattlers.size() == nbrs.size());
	assert(rattlers.size() == N);
	int num_total_rattlers = 0;
	int rattlers_found = 0;
	int new_rattlers_found;
	if(Verbose)
		printf("Identifying nodes with less than %i bonds...\t", c);

	//If a particle is already designated a rattler, clear its nbrs list...
	for(int i=0; i<N; ++i) 
		if(rattlers[i])
		{
			++num_total_rattlers;
			nbrs[i].clear();
		}

	// ... and remove it from other particle's nbrs list.
	for(int i=0; i<N; ++i) 
	{    
		for(int j=(int)nbrs[i].size()-1; j>=0; --j) 
		{    
			if(rattlers[nbrs[i][j]])
				nbrs[i].erase(nbrs[i].begin()+j);
		}    
	} 
    //Look for new rattlers recursively.
	do{
		new_rattlers_found = 0;
		for(int i=0; i<N; ++i)
		{
			if((int)nbrs[i].size() < c && !rattlers[i] && !fixed[i])
			{
				RemoveFromNeighborsList(nbrs,i);
				rattlers[i] = true;
				num_total_rattlers++;
				rattlers_found++;
				new_rattlers_found++;
			}
		}
	}while(new_rattlers_found);
	if(Verbose)
		printf("    --Identified % 5i new rattlers--\t", rattlers_found);

	//Perform Checks
	for(int i=0; i<N; ++i)
	{
		if(!fixed[i])
		{
			if(rattlers[i] && ((int)nbrs[i].size())>0)
			{
				printf("Error removing rattlers: rattlers have neighbors!\n");
				exit(EXIT_FAILURE);
			}   
			if(!rattlers[i] && ((int)nbrs[i].size())<c)
			{
				printf("Error removing rattlers: non-rattlers have too few neighbors!\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	//Check the consistency of the nbrs list
	if(!CheckNeighborsConsistency(nbrs))
	{
		printf("Error: nbrs list is not consistent!\n");
		exit(EXIT_FAILURE);
	}

	if(Verbose)
		 printf("    --% 5i total rattlers--\n", num_total_rattlers);

	return num_total_rattlers;
}

template<int Dim>
int CBondList<Dim>::IdentifyRattlers(vector< vector<int> > &nbrs, vector<bool> &rattlers, int c, bool Verbose) const
{
	vector<bool> fixed;	fixed.assign(N,false);
	return IdentifyRattlers(nbrs,rattlers,fixed,c,Verbose);
}

/**
 *	This method calculates the nbrs list, identifies rattlers, sets map, removes any bonds that may involve rattlers, and updates the i and j indices of each bond (as well as N) to account for the removed rattlers.
 *	
 *	@param[out] map An index_map indicating which nodes are rattlers.
 *	@param[in] fixed A vector of bool's that indicates if a node is fixed (true) and thus should never be considered.
 *	@param[in] c The number of bonds involving a node required for that node to not be a rattler. Default is Dim+1.
 *	@param[in] Verbose Bool to determine if number of identified rattlers, etc. should be printed to stdout.
 */
template<int Dim>
void CBondList<Dim>::RemoveRattlers(index_map &map, vector<bool> const &fixed, int c, bool Verbose)
{
	//Calculate the nbrs list
	vector< vector<int> > nbrs;
	CalculateNeighbors(nbrs);
	
	//Identify the rattlers and set the map
	vector<bool> rattlers(N,false);
	IdentifyRattlers(nbrs, rattlers, fixed, c, Verbose);
	map.set_map(rattlers);

	//Remove the bonds that correspond to rattlers
	vector<bool> BondsToRemove(list.size(),false);
	for(int i=0; i<(int)list.size(); ++i)
		if(rattlers[list[i].i] || rattlers[list[i].j])
			BondsToRemove[i] = true;
	RemoveBonds(BondsToRemove);

	//Update the i and j indices of each bond (as well as N) in accordence with the map
	UpdateBondIndices(map);

	if(!CheckConsistency())
		assert(false);
}

template<int Dim>
void CBondList<Dim>::RemoveRattlers(int c, bool Verbose)
{
	index_map map;
	vector<bool> fixed;	fixed.assign(N,false);
	RemoveRattlers(map, fixed, c, Verbose);
}

template<int Dim>
void CBondList<Dim>::RemoveRattlers(vector<bool> const &fixed, int c, bool Verbose)
{
	index_map map;
	RemoveRattlers(map, fixed, c, Verbose);
}
template<int Dim>
void CBondList<Dim>::RemoveRattlers(index_map &map, int c, bool Verbose)
{
	vector<bool> fixed;	fixed.assign(N,false);
	RemoveRattlers(map, fixed, c, Verbose);
}


template<int Dim>
void CBondList<Dim>::RemoveBondsAccordingToMap(index_map const &map)
{
	assert(N==map.full_size);
	std::vector<bool> rattlers(N,false);
	for(int ii=0;ii <N;++ii)
	{
		if(map.inv(ii)==-1) rattlers[ii]=true;
	};

	//Remove the bonds that correspond to rattlers
	std::vector<bool> BondsToRemove(list.size(),false);
	for(int i=0; i<(int)list.size(); ++i)
		if(rattlers[list[i].i] || rattlers[list[i].j])
			BondsToRemove[i] = true;
	RemoveBonds(BondsToRemove);

	//Update the i and j indices of each bond (as well as N) in accordence with the map
	//UpdateBondIndices(map);
	//DMS: actually, no…if one then runs another "check for rattlers" sweep with a new
	//definition of contacts needed to be a non-rattler, it is important that the original
	//indices are preserved
}


template<int Dim>
void CBondList<Dim>::MakeUnstressed()
{
	MultiplyForces(0.);
}

template<int Dim>
void CBondList<Dim>::MultiplyForces(dbl m)
{
	for(typename vector<BOND>::iterator b=list.begin(); b!=list.end(); ++b)
		b->g *= m;
}

template<int Dim>
void CBondList<Dim>::MultiplyStiffnesses(dbl m)
{
	for(typename vector<BOND>::iterator b=list.begin(); b!=list.end(); ++b)
		b->k *= m;
}

template<int Dim>
dbl  CBondList<Dim>::ComputeEnergy() const
{
	dbl energy=0.;
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
		energy += b->E;
	return energy;
}

template<int Dim>
dbl  CBondList<Dim>::ComputePressure() const
{
	dmat stress;
	ComputeStressTensor(stress);
	return ComputePressure(stress);
}

template<int Dim>
dbl  CBondList<Dim>::ComputePressure(dmat &stress) const
{
	return -stress.trace()/((dbl)Dim);
}

template<int Dim>
void CBondList<Dim>::ComputeStressTensor(dmat &stress) const
{
	stress.setZero();
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
		stress += b->g*((b->r)*(b->r.transpose()))/(b->rlen);
	stress /= Volume;
}


template<int Dim>
dbl  CBondList<Dim>::ComputeGradient(Eigen::VectorXd &grad) const
{
	//CheckConsistency();
	grad.setZero(N*Dim);
	dbl energy = 0.;
	dvec temp;
	int dd;
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
	{
		energy += b->E;
		temp = (b->g/b->rlen)*b->r;
		for(dd=0; dd<Dim; ++dd)
		{
			grad[Dim*b->i+dd] += temp[dd];
			grad[Dim*b->j+dd] -= temp[dd];
		}
	}
	return energy;
}

template<int Dim>
dbl CBondList<Dim>::ComputeGradient(Eigen::VectorXd &grad, vector<bool> const &FixedDOF) const
{
	double energy = ComputeGradient(grad);	

	assert(grad.size() == FixedDOF.size());

	for(int i = 0 ; i < grad.size() ; i++)
		if(FixedDOF[i]==true)
			grad[i] = 0.0;
	/*//CheckConsistency();
        grad.setZero(N*Dim);
        dbl energy = 0.;
        dvec temp;
        int dd;
        for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
        {
                energy += b->E;
                temp = (b->g/b->rlen)*b->r;
                for(dd=0; dd<Dim; ++dd)
                {       
			if(!FixedDOF[Dim*b->i+dd])
                        	grad[Dim*b->i+dd] += temp[dd];
                        if(!FixedDOF[Dim*b->j+dd])
	 			grad[Dim*b->j+dd] -= temp[dd];
                }
        }
	
*/
	return energy;
}

template<int Dim>
void CBondList<Dim>::ComputeHessianElements(vector<TRIP> &coefficients, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	dmat Fii, Fjj, Fij; //Fji = Fij.transpoze(); stressed block
	dmat Kii, Kjj, Kij; //Kji = Kij.transpose(); unstressed block
	dmat Bii, Bjj, Bij; //Bji = Bij.transpose(); B = unstress_coeff*K + stress_coeff*F

	int icorner, jcorner;
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
	{
		//Calculate the matrix blocks
		b->CalculateMatrixBlocks(Fij,Kij); // equal to: (*b).CalculateMatrixBlocks(Fij,Kij);
		Fii = Fjj = -Fij;
		Kii = Kjj = -Kij;

		Bii = unstress_coeff*Kii + stress_coeff*Fii;
		Bjj = unstress_coeff*Kjj + stress_coeff*Fjj;
		Bij = unstress_coeff*Kij + stress_coeff*Fij;

		//Add the matrix blocks to the list of coefficients
		icorner = Dim*b->i;
		jcorner = Dim*b->j;
		for(int d1=0; d1<Dim; ++d1)
			for(int d2=0; d2<Dim; ++d2)
			{
				coefficients.push_back( TRIP(icorner+d1, icorner+d2, Bii(d1,d2)) );
				coefficients.push_back( TRIP(jcorner+d1, jcorner+d2, Bjj(d1,d2)) );
				coefficients.push_back( TRIP(icorner+d1, jcorner+d2, Bij(d1,d2)) );
				coefficients.push_back( TRIP(jcorner+d1, icorner+d2, (Bij.transpose())(d1,d2)) );
			}
	}
	
	//add the tether
	for(int ii=0; ii<Dim*N; ++ii)
		coefficients.push_back( TRIP(ii, ii, tether) );
}

//ONLY WORKS IF ALL THE FIXED PARTICLES ARE AT THE END
template<int Dim>
void CBondList<Dim>::ComputeHessianElements(vector<TRIP> &coefficients, vector<bool> const &FixedDOF) const
{
        dmat Fii, Fjj, Fij; //Fji = Fij.transpoze(); stressed block
        dmat Kii, Kjj, Kij; //Kji = Kij.transpose(); unstressed block
        dmat Bii, Bjj, Bij; //Bji = Bij.transpose(); B = unstress_coeff*K + stress_coeff*F

	double unstress_coeff = 1.0;
	double stress_coeff = 1.0;

        int icorner, jcorner;
        for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
        {
                //Calculate the matrix blocks
                b->CalculateMatrixBlocks(Fij,Kij); // equal to: (*b).CalculateMatrixBlocks(Fij,Kij);
                Fii = Fjj = -Fij;
                Kii = Kjj = -Kij;

                Bii = unstress_coeff*Kii + stress_coeff*Fii;
                Bjj = unstress_coeff*Kjj + stress_coeff*Fjj;
                Bij = unstress_coeff*Kij + stress_coeff*Fij;

                //Add the matrix blocks to the list of coefficients
                icorner = Dim*b->i;
                jcorner = Dim*b->j;
                for(int d1=0; d1<Dim; ++d1)
                        for(int d2=0; d2<Dim; ++d2)
                        {
				if(FixedDOF[icorner+d1]==false&&FixedDOF[icorner+d2]==false)
                                	coefficients.push_back( TRIP(icorner+d1, icorner+d2, Bii(d1,d2)) );
                                if(FixedDOF[jcorner+d1]==false&&FixedDOF[jcorner+d2]==false)
					coefficients.push_back( TRIP(jcorner+d1, jcorner+d2, Bjj(d1,d2)) );
                                if(FixedDOF[icorner+d1]==false&&FixedDOF[jcorner+d2]==false)
					coefficients.push_back( TRIP(icorner+d1, jcorner+d2, Bij(d1,d2)) );
                                if(FixedDOF[jcorner+d1]==false&&FixedDOF[icorner+d2]==false)
					coefficients.push_back( TRIP(jcorner+d1, icorner+d2, (Bij.transpose())(d1,d2)) );
                        }
        }
	
}


template<int Dim>
void CBondList<Dim>::PlaneWaveAnsatz(Eigen::SparseMatrix<cdbl> &transformation, Eigen::VectorXd const &Positions, dvec const &k) const
{
	vector<cTRIP> coefficients;
	dbl kdotR;
	cdbl eikdotR;

	for(int i=0; i<Positions.size()/Dim; ++i)
	{
		kdotR = k.dot(Positions.segment<Dim>(Dim*i));
		eikdotR = cdbl(cos(kdotR), sin(kdotR));
		for(int dd=0; dd<Dim; ++dd)
			coefficients.push_back( cTRIP( Dim*i+dd, Dim*i+dd, eikdotR) );
	}

	Eigen::SparseMatrix<cdbl> temp(Positions.size(),Positions.size());
	transformation = temp;
	transformation.setFromTriplets(coefficients.begin(), coefficients.end());
	assert(transformation.isCompressed());
}

/*
template<int Dim>
void CBondList<Dim>::ComputeHessianElements(vector<cTRIP> &coefficients, dvec k, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	dmat Fii, Fjj, Fij; //Fji = Fij.transpoze(); stressed block
	dmat Kii, Kjj, Kij; //Kji = Kij.transpose(); unstressed block
	dmat Bii, Bjj, Bij; //Bji = Bij.transpose(); B = unstress_coeff*K + stress_coeff*F

	dbl kdotr;
	cdbl eikdotr;
	int icorner, jcorner;
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
	{
		//Calculate exp{i k\cdot r} = cos(k \cdot r) + i sin(k \cdot r)
		kdotr = k.dot(b->r);
		eikdotr = cdbl(cos(kdotr), sin(kdotr));

		//Calculate the matrix blocks
		b->CalculateMatrixBlocks(Fij,Kij); // equal to: (*b).CalculateMatrixBlocks(Fij,Kij);
		Fii = Fjj = -Fij;
		Kii = Kjj = -Kij;

		Bii = unstress_coeff*Kii + stress_coeff*Fii;
		Bjj = unstress_coeff*Kjj + stress_coeff*Fjj;
		Bij = unstress_coeff*Kij + stress_coeff*Fij;

		//Add the matrix blocks to the list of coefficients
		icorner = Dim*b->i;
		jcorner = Dim*b->j;
		for(int d1=0; d1<Dim; ++d1)
			for(int d2=0; d2<Dim; ++d2)
			{
				coefficients.push_back( cTRIP(icorner+d1, icorner+d2, Bii(d1,d2)) );
				coefficients.push_back( cTRIP(jcorner+d1, jcorner+d2, Bjj(d1,d2)) );
				coefficients.push_back( cTRIP(icorner+d1, jcorner+d2, eikdotr*Bij(d1,d2)) );
				coefficients.push_back( cTRIP(jcorner+d1, icorner+d2, std::conj(eikdotr)*(Bij.transpose())(d1,d2)) );
			}
	}
	
	//add the tether
	for(int ii=0; ii<Dim*N; ++ii)
		coefficients.push_back( cTRIP(ii, ii, cdbl(tether,0)) );
}
*/

/**
 *	This is the primary method for calculating the hessian, which is returned as an Eigen::SparseMatrix<dbl> through the parameter hess.
 *
 *	@param[out] hess Eigen::SparseMatrix<dbl> representing the hessian.
 *	@param[in] unstress_coeff A dbl indicating the weight given to the unstressed component of the hessian. For most purposes, this should be set to 1 (default).
 *	@param[in] stress_coeff A dbl indicating the weight given to the stressed component of the hessian. For most purposes, this should be set to 1 (default). 
 *	Set this to 0 (and unstress_coeff to 1) to generate the hessian for the unstressed system.
 *	@param[in] tether A dbl indicating the strength of the tether. A value of tether is added to every diagonal element of the hessian. For most purposes, this should be set to 0 (default).
 */
template<int Dim>
void CBondList<Dim>::ComputeHessian(Eigen::SparseMatrix<dbl> &hess, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	//First, compute the complex matrix elements with k=0
	vector<TRIP> coeffs;
	ComputeHessianElements(coeffs, unstress_coeff, stress_coeff, tether);
	
	//Create a temporary complex matrix
	Eigen::SparseMatrix<dbl> temp(Dim*N,Dim*N);
	hess = temp;
	hess.setFromTriplets(coeffs.begin(), coeffs.end());
	assert(hess.isCompressed());
}

template<int Dim>
void CBondList<Dim>::ComputeHessian(Eigen::SparseMatrix<dbl> &hess,  vector<bool> const &FixedDOF) const
{
        //First, compute the complex matrix elements with k=0
        vector<TRIP> coeffs;
        ComputeHessianElements(coeffs, FixedDOF);
        
	//Calculate the dimensions of the matrix
	int FreeDOF = 0;
	for(int i = 0 ; i < Dim*N ; i++)
		if(!FixedDOF[i])
			FreeDOF++;

	cout << "FREE DOF = " << FreeDOF << endl;

        //Create a temporary complex matrix
        Eigen::SparseMatrix<dbl> temp(FreeDOF,FreeDOF);
        hess = temp;
        hess.setFromTriplets(coeffs.begin(), coeffs.end());
        assert(hess.isCompressed());
}


/**
 *	This method calculates the hessian at a non-zero wavevector, returning it as an Eigen::SparseMatrix<cdbl> through the parameter hess.
 *
 *	@param[out] hess Eigen::SparseMatrix<cdbl> representing the hessian.
 *	@param[in] k A dvec indicating the wavevector.
 *	@param[in] unstress_coeff A dbl indicating the weight given to the unstressed component of the hessian. For most purposes, this should be set to 1 (default).
 *	@param[in] stress_coeff A dbl indicating the weight given to the stressed component of the hessian. For most purposes, this should be set to 1 (default). 
 *	Set this to 0 (and unstress_coeff to 1) to generate the hessian for the unstressed system.
 *	@param[in] tether A dbl indicating the strength of the tether. A value of tether is added to every diagonal element of the hessian. For most purposes, this should be set to 0 (default).
 */
/*
template<int Dim>
void CBondList<Dim>::ComputeHessian_BZ(Eigen::SparseMatrix<cdbl> &hess, dvec k, dbl unstress_coeff, dbl stress_coeff, dbl tether) const
{
	vector<cTRIP> coeffs;
	ComputeHessianElements(coeffs, k, unstress_coeff, stress_coeff, tether);

	Eigen::SparseMatrix<cdbl> temp(Dim*N,Dim*N);
	hess = temp;
	//hess.setZero();
	hess.setFromTriplets(coeffs.begin(), coeffs.end());
	assert(hess.isCompressed());
}
*/

template<int Dim>
void CBondList<Dim>::ComputeEquilibriumMatrix(Eigen::SparseMatrix<dbl> &Amatrix) const
{
	vector<TRIP> coefficients;
	dbl temp;
	for(int bi=0; bi<(int)list.size(); ++bi)
	{
		BOND const &b = list[bi];
		for(int dd=0; dd<Dim; ++dd)
		{   
			temp = b.r[dd]/b.drlen;
			coefficients.push_back( TRIP(Dim*b.i+dd, bi, -temp) );
			coefficients.push_back( TRIP(Dim*b.j+dd, bi,  temp) );
		}   
	}

	Amatrix = Eigen::SparseMatrix<dbl>(Dim*N,(int)list.size());
	Amatrix.setFromTriplets(coefficients.begin(), coefficients.end());
	assert(Amatrix.isCompressed());
}

















/**
 *	This method does 2 things.
 *		1. It calculates the affine extension induced on each bond from the stratin tensor.
 *		2. It calculates the induced force on each particle from the straint tensor.
 */
template<int Dim>
void CBondList<Dim>::Calculate_n_d2Udvdgamma(dmat const &strain_tensor, Eigen::VectorXd &n_d2Udvdgamma, Eigen::VectorXd &AffineDeltaR, dbl unstress_coeff, dbl stress_coeff) const
{
	assert(AffineDeltaR.size() == Dim*list.size());
	n_d2Udvdgamma = Eigen::VectorXd::Zero(Dim*N);

	//Calculate the net force on each node
	int i,j;
	dbl g, k, DeltaRparallel;
	dvec D2i;
	for(int bi=0; bi<(int)list.size(); ++bi)
	{
		BOND const &b = list[bi];
		g = b.g*stress_coeff;
		k = b.k*unstress_coeff;

		//Calculate the affine extension of the bond.
		AffineDeltaR.segment<Dim>(Dim*bi) = strain_tensor*b.r;

		//Project that onto the direction of b.r;
		DeltaRparallel = AffineDeltaR.segment<Dim>(Dim*bi).dot(b.r)/b.rlen;

		//Calculate the induced forces on particles i and j.
		D2i = (DeltaRparallel/b.rlen) * (k - g/b.rlen) * b.r  +  (g/b.rlen) * AffineDeltaR.segment<Dim>(Dim*bi);
//		for(int d1=0; d1<TDIM; d1++)
//		{   //check the sign here!
//			df_on_i[d1] = (bnd.r[d1]/bnd.drlen)*DeltaRparallel * (k + f/bnd.drlen) - (f/bnd.drlen)*DeltaRtt.r[d1];
//			df_on_j[d1] = -df_on_i[d1];
//		} 

		//Add the induced force to the n_d2Udvdgamma vector
		n_d2Udvdgamma.segment<Dim>(Dim*b.i) += D2i;
		n_d2Udvdgamma.segment<Dim>(Dim*b.j) -= D2i;
	}
}

/**
 *	This method calculates the linear order change in energy when each bond is extended according to DeltaR_bond.
 *	The change in energy of a bond is (1/2)d2Udgamma2.
 */
template<int Dim>
dbl CBondList<Dim>::CalculateEnergyChange(Eigen::VectorXd const &DeltaR_bond, dbl unstress_coeff, dbl stress_coeff) const
{
	dbl energy = 0.;
	dbl g, k, DeltaRparallel2, DeltaRperp2;

	for(int bi=0; bi<(int)list.size(); ++bi)
	{
		BOND const &b = list[bi];
		g = b.g*stress_coeff;
		k = b.k*unstress_coeff;
		DeltaRparallel2 = POW2(DeltaR_bond.segment<Dim>(Dim*bi).dot(b.r)/b.rlen);
		DeltaRperp2     = DeltaR_bond.segment<Dim>(Dim*bi).squaredNorm() - DeltaRparallel2;

		energy += k*DeltaRparallel2 + g*DeltaRperp2/b.rlen;
	}
	energy /= 2.;
	return energy;
};

/**
 *	This method calculates all the elastic constants.
 *
 *	@param[out] cijkl An object of type cCIJKL<Dim> that stores the elastic constants.
 *	@param[in] unstress_coeff Multiplicative factor applied to all spring constants.
 *	@param[in] stress_coeff Multiplicative factor applied to all forces.
 *	@param[in] Verbose if true(default), print out the elastic constants at the end.
 */
template<int Dim>
void CBondList<Dim>::CalculateCijkl(cCIJKL<Dim> &cijkl, dbl unstress_coeff, dbl stress_coeff, bool Verbose) const
{
	int Nvar = Dim*N;

	MatrixInterface<dbl> D1;
	ComputeHessian(D1.A, unstress_coeff, stress_coeff, 1e-12);
	D1.LUdecomp();

	dmat strain_tensor;
	dbl prefactor = 2./Volume;//The prefactor contains a 2 in it from the equation dE/V = (1/2)cijkl uij ukl.

	Eigen::VectorXd n_d2Udvdgamma = Eigen::VectorXd::Zero(Nvar);
	Eigen::VectorXd uNonAffine_node = Eigen::VectorXd::Zero(Nvar);
	Eigen::VectorXd DeltaR_bond = Eigen::VectorXd::Zero(Dim*(int)list.size());

	for(int ii=0; ii<cCIJKL<Dim>::num_constants; ++ii)
	{
		//Set the strain tensor
		cijkl.set_strain_tensor(strain_tensor, ii);

//Start here...

		//Calculate the displacement vector for each bond in the metric defined by the strain tensor.
		//Also, calculate the forces on each particle due to the change in metric.
		Calculate_n_d2Udvdgamma(strain_tensor, n_d2Udvdgamma, DeltaR_bond, unstress_coeff, stress_coeff); 
		
		//Solve for the non-affine displacement
		D1.solve_Mx_equals_b(uNonAffine_node, n_d2Udvdgamma);

		//Add the non-affine extension of each bond to the affine extension
		for(int bi=0; bi<(int)list.size(); ++bi)
			DeltaR_bond.segment<Dim>(Dim*bi) += uNonAffine_node.segment<Dim>(Dim*list[bi].j) - uNonAffine_node.segment<Dim>(Dim*list[bi].i);

//End here. now look for whatever...

		//Calculate the change in energy
		dbl dE = CalculateEnergyChange(DeltaR_bond, unstress_coeff, stress_coeff);
		dE *= prefactor; //dE is now 2*dE/V. This comes from the equation dE/V = (1/2)cijkl uij ukl.

		//Set the ii'th elastic constant
		cijkl.set_constant(dE, ii);
	}
	if(Verbose)
		cijkl.print();
}

template<int Dim>
void CBondList<Dim>::CalculateStdData(CStdData<Dim> &Data, bool CalcCijkl, bool CalcHess) const
{
	Data.SetZero();
	if(!Empty())
	{
		Data.NPp = N;
		Data.Nc = (int)list.size();
		Data.Volume = Volume;

		Data.Energy = ComputeEnergy();
		ComputeStressTensor(Data.Stress);
		Data.Pressure = ComputePressure(Data.Stress);

		Eigen::VectorXd grad;
		ComputeGradient(grad);
		Data.MaxGrad = max_abs_element(grad);

		Data.cijkl.SetArtificialValues(-1e12);
		if(CalcCijkl)
		{
			CalculateCijkl(Data.cijkl, 1., 1., false);
		}
		if(CalcHess)
			ComputeHessian(Data.H.A);
	}
}









template<int Dim>
void CBondList<Dim>::PrintBonds() const
{
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
		b->print();
}

template<int Dim>
void CBondList<Dim>::PrintBonds_txt(char *filename) const
{
	FILE *outfile = fopen(filename, "w");
	for(typename vector<BOND>::const_iterator b=list.begin(); b!=list.end(); ++b)
	{
		fprintf(outfile, "%5i %5i % e % e % e % e ", b->i, b->j, b->rlen, b->E, b->g, b->k);
		for(int dd=0; dd<Dim; ++dd)
			fprintf(outfile, "% e ", b->r[dd]);
		fprintf(outfile, "\n");
	}
	fflush(outfile);
	fclose(outfile);
}

template<int Dim>
bool CBondList<Dim>::CheckConsistency() const
{
	if(Empty())
	{
		if(N==0)
			return true;
		return false;
	}
	else
	{
		//Calculate the largest and smallest particle index
		int imin, imax;
		imin = imax = list[0].i;
		for(int b=1; b<(int)list.size(); ++b)
		{
			if(imin > list[b].i) imin = list[b].i;
			if(imax < list[b].i) imax = list[b].i;
		}
		for(int b=0; b<(int)list.size(); ++b)
		{
			if(imin > list[b].j) imin = list[b].j;
			if(imax < list[b].j) imax = list[b].j;
		}

		//printf("imin = %i, imax = %i, N = %i\n", imin, imax, N);

		if(imin >= 0 && imax < N)
			return true;
		assert(false);
	}
	return false;
}





}

#endif //BOND_LIST_H

