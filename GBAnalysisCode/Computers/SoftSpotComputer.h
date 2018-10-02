#ifndef SOFTSPOTCOMPUTER

#define SOFTSPOTCOMPUTER

/////////////////////////////////////////////////////////////////////////////////
//Soft Spot Computer class.
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
#include "BaseComputer.h"
#include "BondList.h"
#include "Grid.h"
#include "Data.h"
#include <list>
#include <stack>

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
class CSoftSpotComputer : public CBaseComputer<Dim>
{
private:
	typedef Eigen::Matrix<dbl,Dim,1>   dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
	CGrid<Dim> Grid;
	bool UseGrid;
	CStaticState<Dim> &State;
	bool FixDof;
	vector<bool> FixedDof;

public:
    	int N_m, N_p; //number of particles and particles per mode used to generate soft spots
    	int ClusterSize;
	CBondList<Dim> Bonds;	//!<Bond list for standard computations
	index_map RattlerMap;	//!<Map for rattlers.
	CStdData<Dim>  Data;	//!<storage for standard data
	CStaticState<Dim> Lattice; //underlying lattice to lift planewaves
	
	CSoftSpotComputer();
	CSoftSpotComputer(CStaticState<Dim> &_State);
	CSoftSpotComputer(CStaticState<Dim> &_State, bool diagonalize);

	
	void Print();
	

        //Number of particles in the system
	int N;
	int N_lattice;
    
	void SetParameters(int nm,int np,int cluster);   
        //list of the soft modes
    	vector<Eigen::VectorXd> SoftSpots;
    	vector<int> n;
 	vector<int> softparticles; //1 if a particle is in a soft spot, 0 otherwise
   
	//a function to actually compute the soft spots given a system and the mode structure (ordered by
	//energy
	void Construct(Eigen::MatrixXd &polarizations, int zero_modes);
	
	//a function to actually compute the soft spots given a system and the mode structure (ordered by
	//energy
	void Construct(double *polarizations, int zero_modes);
	void SoftParticles(double * polarizations, int zero_modes);
    	

    //soft spot functionality
	const vector<Eigen::VectorXd> &GetSoftSpots() {return SoftSpots;}//returns a reference to a vector of soft spots.
	int Count() { return SoftSpots.size();} //returns the number of soft spots
    
    //returns the particles with the N_p largest polarization vectors for the lowest N_m modes.
    	static vector<vector<int> > GetLargestPolarizations();
    
    
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

	void BondListFull(CBondList<Dim> &bonds,dbl spring);


	int  StdPrepareSystem(bool verbose = true);								//!<Compute the bonds list and remove rattlers
	int  StdPrepareSystem(vector<bool> const &fixed, bool verbose = true);	//!<Compute the bonds list and remove rattlers, assuming some fixed particles

    //Calculate standard data
	void CalculateStdData(bool CalcCijkl = true, bool CalcHess = true);
	void CalculateStdData(CStdData<Dim> &data, bool CalcCijkl=true, bool CalcHess = true);
	void CalculateStdData_Unstressed(bool CalcCijkl=true, bool CalcHess = true);
	void CalculateStdData_Unstressed(CStdData<Dim> &data, bool CalcCijkl=true, bool CalcHess = true);

    //Compute Hessian
	void ComputeHessian(Eigen::SparseMatrix<dbl> &hess) const;
	void ComputeHessian(Eigen::SparseMatrix<dbl> &hess, bool fixed) const;
	void ComputeHessian(CBondList<Dim> const &bonds, Eigen::SparseMatrix<dbl> &hess) const;
	

    //Needed for minimization routines
	void SetFixedDof(vector<bool> const &_FixedDof);
	void UnsetFixedDof();

	void Evaluate(Eigen::VectorXd &grad, dbl &fx);
	bool Progress(Eigen::VectorXd const &grad, dbl fx, int iteration, int print_iter, dbl tol);
	void Move(Eigen::VectorXd const &step);
	dbl GetMinimizationTimeScale() const;
    
	//Things to to with the lattice
	void SetHexLattice();
	void ComputeLatticeSpring(dbl rlen2, dbl &E, dbl &g, dbl &k, dbl spring);
	dbl lattice_spacing, gauss_rad;
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION    //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
CSoftSpotComputer<Dim>::CSoftSpotComputer() : FixDof(false), UseGrid(true)
{
}


template <int Dim>
CSoftSpotComputer<Dim>::CSoftSpotComputer(CStaticState<Dim> &_State) : State(_State), Grid(&_State), FixDof(false), UseGrid(true)
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
        };
        if(!enough_cells)
                SetUseGrid(false);
 
	//First, get the dynamical matrix and count the number of zero modes
        N=State.GetParticleNumber();
        StdPrepareSystem();
        CalculateStdData();
       	int ms = (int) floor(N/8); 
	Data.H.Diagonalize(ms);
        int zero_modes = 0;
        int ii=0;
        while(Data.H.GetVal(ii) < 1e-8)
        {
            zero_modes += 1;
            ii += 1;
        };
        //construct the soft spots using default parameters
	N_m = 10;
	N_p = 10;
	ClusterSize = 5;
        Construct(Data.H.Eigenvectors,zero_modes);

};

template <int Dim>
CSoftSpotComputer<Dim>::CSoftSpotComputer(CStaticState<Dim> &_State, bool diagonalize) : State(_State), Grid(&_State), FixDof(false), UseGrid(true)
{
        //use grid or not
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
        //pre-diagonalize?
        N=State.GetParticleNumber();
	cout <<"N= " << N << endl;
	if(diagonalize)
	{
		StdPrepareSystem();
        	CalculateStdData();
		int ms = (int) floor(N/8);
        	Data.H.Diagonalize(ms);
	};
};


template <int Dim>
void CSoftSpotComputer<Dim>::SetParameters(int nm,int np,int cluster)
{
	N_m=nm;
	N_p=np;
	ClusterSize=cluster;
};

template <int Dim>
void CSoftSpotComputer<Dim>::SetHexLattice()
{
	//SPECIALIZED TO 2D
	//
	lattice_spacing = 10.0; //make variable?
	gauss_rad = lattice_spacing*0.5;
	//figure out number of rows, cols, based on lattice and box
	dmat Trans;
        State.GetBox()->GetTransformation(Trans);	
	int N_per_row, N_per_col;
	N_per_row = (int) floor((Trans(0,0)-0.5*lattice_spacing)/lattice_spacing)+1;
	dbl height = 0.5*lattice_spacing*sqrt(3.0);
	N_per_col = (int) floor((Trans(1,1))/height)+1;

	//Initialize the "lattice" state to the right box and number of particles
	N_lattice = N_per_row*N_per_col;
	CStaticState<Dim> temp_state(N_lattice);
	temp_state.GetBox()->SetTransformation(Trans);
	Lattice = temp_state;

	//arrange particles in "lattice" state in a triangular arrangement
	dvec part_pos;
	dbl vert_offset = 0.5*(Trans(1,1)-height*(N_per_col-1));
	dbl h_offset = 0.5*(Trans(0,0)-lattice_spacing*(N_per_row-.5));
	for (int rr = 0; rr < N_per_col ; ++rr)
	{	
		for (int cc = 0; cc < N_per_row; ++cc)
		{
			int index = N_per_row*rr+cc; //particle index currently being positioned
			part_pos(1) = vert_offset+height*rr;
			part_pos(0) = h_offset+lattice_spacing*cc + 0.5*lattice_spacing*(rr%2);
			Lattice.SetParticlePosition(part_pos,index);
			Lattice.SetRadius(index, lattice_spacing*0.5);
		};
	};
};

template <int Dim>
void CSoftSpotComputer<Dim>::BondListFull(CBondList<Dim> &bonds, dbl spring)
{
        //Make sure bonds is empty.
        bonds.RemoveAllBonds();

        Grid.Construct();
        //Grid.PrintGrid();

        bonds.SetN(State.GetParticleNumber()+N_lattice);//particles plus lattice points
        //first, compute standard system bonds
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
	Lattice.SetPotentialHarmonic(1./spring);
	//now add all of the virtual bonds with the lattice
	dvec part_pos1, part_pos2;
	for (int ii = 0; ii< Lattice.GetParticleNumber(); ++ii)
	{
		for (int i2 = ii+1; i2 < Lattice.GetParticleNumber(); ++i2)
		{
			Lattice.GetParticlePositionVirtual(part_pos1,ii);
			Lattice.GetParticlePositionVirtual(part_pos2,i2);
			Lattice.GetBox()->MinimumDisplacement(part_pos1,part_pos2,Displacement);
			Lattice.GetBox()->Transform(Displacement);
                        rlen2=Displacement.squaredNorm();
			if(rlen2 < (lattice_spacing*1.2*1.2*lattice_spacing)&& rlen2 > (.9*lattice_spacing*lattice_spacing))
			{
				rlen = sqrt(rlen2);
				Lattice.GetPotential()->ComputeDerivatives012(rlen, rlen, E, g, k);
                                bonds.AddBond( CBond<Dim>(ii+State.GetParticleNumber(), i2+State.GetParticleNumber(), rlen, E, g, k, Displacement) );	
			
			};
		};	

		for (int jj = 0; jj< State.GetParticleNumber(); ++jj)
		{
			Lattice.GetParticlePositionVirtual(part_pos1,ii);
			State.GetParticlePositionVirtual(part_pos2,jj);
			State.GetBox()->MinimumDisplacement(part_pos1,part_pos2,Displacement);
			Lattice.GetBox()->Transform(Displacement);
			rlen2=Displacement.squaredNorm();
			//reduce number of bonds
			if (rlen2 < 1*lattice_spacing*lattice_spacing)
			{
				ComputeLatticeSpring(rlen2,E,g,k,spring);
				bonds.AddBond( CBond<Dim>(jj,State.GetParticleNumber()+ii,sqrt(rlen2),E,g,k,Displacement) );
			};

		}; 

	};

};

template <int Dim>
void CSoftSpotComputer<Dim>::ComputeLatticeSpring(dbl rlen2, dbl &E, dbl &g, dbl &k,dbl spring)
{
	dbl sig2 = gauss_rad*gauss_rad;
	dbl arg = -rlen2/(2*sig2);
	E = exp(arg)/sqrt(2*M_PI*sig2)*spring;
	g = -sqrt(rlen2)*E/sig2;
	g = -E;
	k = E*((rlen2-sig2)/(sig2*sig2));
	k = E;
};

//get and set the state
template <int Dim>
CStaticState<Dim> &CSoftSpotComputer<Dim>::GetState() const
{
	return State;
}

template <int Dim>
void CSoftSpotComputer<Dim>::SetState(CStaticState<Dim> & _State)
{
	State = _State;
	Grid.SetState(&State);
}
	
template <int Dim>
inline dbl CSoftSpotComputer<Dim>::GetVolume() const
{
	return State.GetVolume();
}

template <int Dim>
inline int CSoftSpotComputer<Dim>::GetNdof() const
{
	return State.GetParticleNumber()*Dim;
}

template <int Dim>
inline dbl CSoftSpotComputer<Dim>::GetMinimizationTimeScale() const
{
	dbl avgRad = State.GetAverageRadius();
	return avgRad; //This is emperical!!! NOT GENERAL!
}

template <int Dim>
void CSoftSpotComputer<Dim>::SetUseGrid(bool usegrid)
{
	UseGrid = usegrid;
}

//Compute the bond list
template <int Dim>
void CSoftSpotComputer<Dim>::ComputeBondList(CBondList<Dim> &bonds)
{
	if(UseGrid)
		ComputeBondList_Grid(bonds);
	else
		ComputeBondList_NoGrid(bonds);
}

template <int Dim>
void CSoftSpotComputer<Dim>::ComputeBondList_Grid(CBondList<Dim> &bonds)
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
void CSoftSpotComputer<Dim>::ComputeBondList_NoGrid(CBondList<Dim> &bonds) const
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
int CSoftSpotComputer<Dim>::StdPrepareSystem(bool verbose)
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
int CSoftSpotComputer<Dim>::StdPrepareSystem(vector<bool> const &fixed, bool verbose)
{
	ComputeBondList(Bonds);										//Use the member variable derived from CBaseComputer
	Bonds.RemoveRattlers(RattlerMap, fixed, Dim+1, verbose);	//Remove rattlers
	return (Bonds.Empty())?1:0;
}

template <int Dim>
void CSoftSpotComputer<Dim>::CalculateStdData(bool CalcCijkl, bool CalcHess)
{
	CalculateStdData(Data,CalcCijkl,CalcHess);
}
template <int Dim>
void CSoftSpotComputer<Dim>::CalculateStdData(CStdData<Dim> &data, bool CalcCijkl, bool CalcHess)
{
	Bonds.CalculateStdData(data, CalcCijkl,CalcHess);
}


template <int Dim>
void CSoftSpotComputer<Dim>::CalculateStdData_Unstressed(bool CalcCijkl, bool CalcHess)
{
	CalculateStdData_Unstressed(Data,CalcCijkl,CalcHess);
}


template <int Dim>
void CSoftSpotComputer<Dim>::CalculateStdData_Unstressed(CStdData<Dim> &data, bool CalcCijkl, bool CalcHess)
{
	CBondList<Dim> BondsTemp = Bonds;
	BondsTemp.MakeUnstressed();
	BondsTemp.CalculateStdData(data,CalcCijkl,CalcHess);
}

template <int Dim>
void CSoftSpotComputer<Dim>::ComputeHessian(Eigen::SparseMatrix<dbl> &hess, bool fixed) const
{
	Bonds.ComputeHessian(hess,FixedDof);
}

template <int Dim>
void CSoftSpotComputer<Dim>::ComputeHessian(Eigen::SparseMatrix<dbl> &hess) const
{
	ComputeHessian(Bonds, hess);
}

template <int Dim>
void CSoftSpotComputer<Dim>::ComputeHessian(CBondList<Dim> const &bonds, Eigen::SparseMatrix<dbl> &hess) const
{
	bonds.ComputeHessian(hess);
}

    
template <int Dim>
void CSoftSpotComputer<Dim>::Print()
{
	int spots = SoftSpots.size();
	for (int ii = 0; ii< spots; ++ii)
	{
		cout << "particles in spot " <<ii<< "::"<<endl;
		for(int jj = 0; jj < SoftSpots[ii].size(); ++jj)
			if(SoftSpots[ii][jj]==1)    cout << jj<<"  ";
		cout << endl;
	};
};  
    
    
    
//Needed for minimization routines
    
template <int Dim>
void CSoftSpotComputer<Dim>::SetFixedDof(vector<bool> const &_FixedDof)
{
{
        FixDof = true;
        FixedDof = _FixedDof;
        assert(FixedDof.size() == Dim*State.GetParticleNumber());
}
};
    
template <int Dim>
void CSoftSpotComputer<Dim>::UnsetFixedDof()
{
	FixDof = false;
};

template <int Dim>
void CSoftSpotComputer<Dim>::Evaluate(Eigen::VectorXd &grad, dbl &fx)
{
};
    
template <int Dim>
bool CSoftSpotComputer<Dim>::Progress(Eigen::VectorXd const &grad, dbl fx, int iteration, int print_iter, dbl tol)
{
    return true;
};
   
template <int Dim>
inline void CSoftSpotComputer<Dim>::Move(Eigen::VectorXd const &step)
{};
   

//soft spot functions
template<int Dim>
vector<vector<int> > CSoftSpotComputer<Dim>::GetLargestPolarizations()
{
};



//particles in soft spots, no clustering
template<int Dim>
void CSoftSpotComputer<Dim>::SoftParticles(double * polarizations, int zero_modes)
{
	softparticles.resize(N);
	std::fill(softparticles.begin(), softparticles.end(), 0);
	//list of the particles in each mode with the N_p largest polarization vectors.
    	//cout << "Computing largest polarization vectors.\n";
    	vector<vector<int> > lowest;
    	for(int i = 0;i<N_m;i++)
    	{
        	lowest.push_back(vector<int>());
        	double max = 1000000000.0;

        	for(int j = 0;j<N_p;j++)
        	{
            		double local_max = 0.0;
            		int select = 0;
            		for(int k = 0;k<N;k++)
            		{
                		double mag = 0.0;
                		for(int l = 0 ; l<Dim ; l++)
                        	mag+=polarizations[i*N*Dim + k*Dim+l]*polarizations[i*N*Dim + k*Dim+l];

                		if(mag>local_max&&mag<max)
                		{
                    			select = k;
                    			local_max = mag;
                		};
            		};
            		lowest[i].push_back(select);
            		max = local_max;
        	};
    	};
	
	for(int i = 0;i<N_m;i++)
          	for(int j = 0;j<N_p;j++)
        		softparticles[lowest[i][j]]=1;
	
}; 
    
//SoftSpot construction routines
template <int Dim>
void CSoftSpotComputer<Dim>::Construct(double *polarizations, int zero_modes)
{
    SoftSpots.clear();
    n.clear();
    //first, compute the neighbor list
    vector<vector<int> > fullneighborlist;
    Bonds.CalculateNeighbors(fullneighborlist,1.5);

    //list of the particles in each mode with the N_p largest polarization vectors.
    cout << "Computing largest polarization vectors.\n";
    vector<vector<int> > lowest;
    for(int i = 0;i<N_m;i++)
    {
        lowest.push_back(vector<int>());
        double max = 1000000000.0;
        
        for(int j = 0;j<N_p;j++)
        {
            double local_max = 0.0;
            int select = 0;
            for(int k = 0;k<N;k++)
            {
            	double mag = 0.0;
                for(int l = 0 ; l<Dim ; l++)
                	mag+=polarizations[i*N*Dim + k*Dim+l]*polarizations[i*N*Dim + k*Dim+l];
                
                if(mag>local_max&&mag<max)
                {
                    select = k;
                    local_max = mag;
                }
            }
            lowest[i].push_back(select);
            max = local_max;
        }
    }
    
    //construct connected neighborhoods of particles that are among "lowest"
    cout << "Constructing connected neighborhoods.\n";
    vector<vector<int> > neighborhoods;
    for(int i = 0;i<N_m;i++)
    {
        for(int j = 0;j<N_p;j++)
        {
            //determine whether this particle is already part of a neighborhood
            bool apart = false;
            for(int k = 0;k<neighborhoods.size();k++)
            {
                for(int l = 0;l<neighborhoods[k].size();l++)
                    if(neighborhoods[k][l]==lowest[i][j])
                    {
                        apart = true;
                        k = neighborhoods.size();
                        break;
                    }
            }
            
            
            //if not, construct a neighborhood based on particle this particle
            if(!apart)
            {
                vector<int> neighborhood;
                
                stack<int> temporary;
                temporary.push(lowest[i][j]);
                while(temporary.size()>0)
                {
                    int current = temporary.top();
                    neighborhood.push_back(current);
                    
                    temporary.pop();
                    
                    vector<int> neighbors = fullneighborlist[current];
                    
                    //go through the neighbors and add any to the list that are in one of the lowest 30 modes that aren't already in the neighborhood
                    for(int n = 0;n<neighbors.size(); n++)
                    {
                        bool inneighborhood = false;
                        for(int k = 0;k<neighborhood.size();k++)
                            if(neighbors[n]==neighborhood[k])
                            {
                                inneighborhood=true;
                                break;
                            }
                        
                        if(!inneighborhood)
                        {
                            for(int k = 0;k<N_m;k++)
                                for(int l = 0;l<N_p;l++)
                                    if(neighbors[n]==lowest[k][l]) {temporary.push(neighbors[n]); k = N_m; /*cout << "B, ";*/ break;}
                        }
						
                    }
                }
                
                
                neighborhoods.push_back(neighborhood);
            }
        }
    }
    
    //finally arrange these neighborhoods into the projection operators.
    cout << "Constructing projection operators.\n";
    for(int i = 0;i<neighborhoods.size();i++)
    {
        Eigen::VectorXd neighborhood = Eigen::ArrayXd::Zero(N);
        if(neighborhoods[i].size()>=ClusterSize)
	{
		cout << "writing operator " << i << endl;
		for(int j = 0;j<neighborhoods[i].size();j++)
			neighborhood(neighborhoods[i][j]) = 1.0;
		SoftSpots.push_back(neighborhood);
		n.push_back(neighborhoods[i].size());
        };
    }
    cout << "Finished constructing soft-spots.";
    for (int ii = 0; ii < n.size(); ++ii)
	cout << n[ii]<<"  ";
    cout << endl;
};

template<int Dim>
void CSoftSpotComputer<Dim>::Construct(Eigen::MatrixXd &polarizations, int zero_modes)
{
    SoftSpots.clear();
    n.clear();	
    //first, compute the neighbor list
    vector<vector<int> > fullneighborlist;
    Bonds.CalculateNeighbors(fullneighborlist);
    
    //list of the particles in each mode with the N_p largest polarization vectors.
    cout << "Computing largest polarization vectors.\n";
    vector<vector<int> > lowest;
    for(int i = 0;i<N_m;i++)
    {
        lowest.push_back(vector<int>());
        double max = 1000000000.0;
        
        
        VectorXd pol = polarizations.col(zero_modes+i);
        System.ConvertFromMassNormalized(pol);
        for(int j = 0;j<N_p;j++)
        {
            double local_max = 0.0;
            int select = 0;
            for(int k = 0;k<N;k++)
            {
                double mag = pol.segment<Dim>(Dim*k).dot(pol.segment<Dim>(Dim*k));
                if(mag>local_max&&mag<max)
                {
                    select = k;
                    local_max = mag;
                };
            };
            lowest[i].push_back(select);
            max = local_max;
        };
    };



    //construct connected neighborhoods of particles that are among "lowest"
    cout << "Constructing connected neighborhoods.\n";
    vector<vector<int> > neighborhoods;
    for(int i = 0;i<N_m;i++)
    {
        for(int j = 0;j<N_p;j++)
        {
            //determine whether this particle is already part of a neighborhood
            bool apart = false;
            for(int k = 0;k<neighborhoods.size();k++)
            {
                for(int l = 0;l<neighborhoods[k].size();l++)
                    if(neighborhoods[k][l]==lowest[i][j])
                    {
                        apart = true;
                        k = neighborhoods.size();
                        break;
                    };
            };
            
            
            //if not, construct a neighborhood based on particle this particle
            if(!apart)
            {
                vector<int> neighborhood;
                
                stack<int> temporary;
                temporary.push(lowest[i][j]);
                while(temporary.size()>0)
                {
                    int current = temporary.top();
                    neighborhood.push_back(current);
                    
                    temporary.pop();
                    
                    vector<int> neighbors = fullneighborlist(current);
                    
                    //go through the neighbors and add any to the list that are in one of the lowest 30 modes that aren't already in the neighborhood
                    for(int n = 0;n<neighbors.size(); n++)
                    {
                        bool inneighborhood = false;
                        for(int k = 0;k<neighborhood.size();k++)
                            if(neighbors[n]==neighborhood[k])
                            {
                                inneighborhood=true;
                                break;
                            };
                        
                        if(!inneighborhood)
                        {
                            for(int k = 0;k<N_m;k++)
                                for(int l = 0;l<N_p;l++)
                                    if(neighbors[n]==lowest[k][l]) {temporary.push(neighbors[n]); k = N_m; i/*cout << "A, ";*/ break;}
                        };
						
                    };
                };

                
                neighborhoods.push_back(neighborhood);
            };
        };
    };
    
    //finally arrange these neighborhoods into the projection operators.
    cout << "Constructing projection operators.\n";
    for(int i = 0;i<neighborhoods.size();i++)
    {
//      cout << "writing operator " << i << endl;
        VectorXd neighborhood = Eigen::ArrayXd::Zero(N);
        if(neighborhoods[i].size()>=ClusterSize)
	{
		for(int j = 0;j<neighborhoods[i].size();j++)
			neighborhood(neighborhoods[i][j]) = 1.0;
		SoftSpots.push_back(neighborhood);
		n.push_back(neighborhoods[i].size());
        };
    }
    cout << "Finished constructing soft-spots.\n";

};

};

#endif
