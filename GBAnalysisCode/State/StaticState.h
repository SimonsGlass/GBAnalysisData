#ifndef STATIC_STATE

#define STATIC_STATE

/////////////////////////////////////////////////////////////////////////////////
//Static state class. 
//
//Description
//		This class stores information to represent systems of isotropic particles
//		and implements functions to manipulate this information. The data that is
//		fundamentally necessary to represent this system are particle positions
//		and radii which are stored as Eigen vectors. Other information that is 
//		required is a potential energy function and boundary conditions. These 
//		ideas are implemented through inherited classes. Writing to and reading from 
//		NetCDF files is implemented.
//
//		Particle positions are stored in the unit box [0,1]x...x[0,1] and then mapped
//		into a system of proper shape by the box class. The stored positions are
//		denoted "virtual" positions whereas the actual positions are denoted "real".
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
//			-> Three dimensions: record number (unlimited), degrees of freedom (DOF), dimension,
//			   and particle number.
//			-> Particle positions are stored as a variable.
//			-> Particle radii are stored as a variable.
//			-> Box is stored as a box (see CBox definition.)
//			-> Potential is stored as a potential (see CPotential definition.)
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include "../Potentials/Potentials.h" //This loads all the potentials
#include "../Boundaries/Boxes.h"
#include "../Boundaries/PeriodicBox.h"
//#include "netcdfcpp.h"
#include "../Resources/MersenneTwister.h"
#include "../Resources/RNG_taus.h"

namespace LiuJamming
{

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim> class CStaticDatabase;

//!Class to store data for a sphere packing, including positions, radii, box and potential
template <int Dim>
class CStaticState
{
private:
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;

	int N;							//!<Number of particles
	Eigen::VectorXd Positions;		//!<dN vector of positions (in normalized coordinates)
	Eigen::VectorXd Radii;			//!<N vector of radii (in real units)
        vector<int> CNAvector;

	CBox<Dim> *Box;					//!<Pointer to a CBox object
	CPotential *Potential;			//!<Pointer to a CPotential object
		
//Functions to check that the correct dimensions, variables, and attributes are present
//in a netCDF file.
//	bool CheckNetCDF(const NcFile &file);
//	void PopulateNetCDF(NcFile &file);
	
public:
        bool CNAisvalid;

//Constructors/Destructors and copy operators
	CStaticState();													//!<Default constructor
	CStaticState(int _N);											//!<Primary constructor
	CStaticState(CStaticState const &copy);							//!<Copy constructor
	CStaticState<Dim> &operator=(CStaticState<Dim> const &copy);	//!<Copy operator
	~CStaticState();												//!<Destructor

	void ClearSystem(int _N);

//Set potential
	void SetPotential(CPotential *t_Potential);	//!<Set the potential
	void SetPotentialHarmonic(dbl epsilon=1.);	//!<Set the potential to be harmonic
	void SetPotentialHertzian(dbl epsilon=1.);	//!<Set the potential to be hertzian
	void SetPotentialLJ(dbl epsilon = 1.0);		//!< Set the potential to be Kob-Andersen
//Set box
	void SetBox(CBox<Dim> *t_Box);				//!<Set the box
	void SetBoxPeriodic(int NonPeriodicDim=0);	//!<Set the box to be periodic
	void AssumeRectangularSymmetry();			//!<Set the assumption that the box is rectangular

//Functions to construct systems 
	void RandomizePositions(long seed = 1);		//!<Randomize the positions using a new seed
	void RandomizePositions(MTRand *random);	//!<Randomize the positions using a previously created random number generator
//Set 2d lattices
	void SetSquareLattice();					//!<Set positions to be in a square lattice (only in 2d)
	void SetSquareLattice(int Lx, int Ly);		//!<Set positions to be in a square lattice (only in 2d)
	void SetHexLattice();						//!<Set positions to be in a hexagonal lattice (only in 2d)
//Set 3d lattices
	void SetFCCLattice();						//!<Set positions to be in an FCC lattice (only in 3d)
	void SetBCCLattice();						//!<Set positions to be in a BCC lattice (only in 3d)

	void SetRadiiMono();										//!<Set all diameters to 1.
	void SetRadiiBi(dbl FracSmall = 0.5, dbl SizeRatio = 1.4);	//!<Set bidisperse size distribution, with average diameter of 1.
	void SetRadiiBiGaussian(dbl s1, dbl s2, dbl FracSmall = 0.5, dbl SizeRatio = 1.4, long seed=123123);	//!<Set bidisperse size distribution, with average diameter of 1, where there is some polydispersity in each of the two sizes.
	void SetRadiiPolyUniform(dbl SizeRatio = 1.4);				//!<Set uniform polydisperse size distribution, with average diameter of 1.
	
//Functions to read and write systems
//	void Read(const NcFile &File,int Record);
//	void Write(NcFile &File,int Record);

//Functions to set properties of the system
	void SetRadii(const Eigen::VectorXd &t_Radii);	//!<Set all the radii
	void SetRadius(int i,dbl r);					//!<Set an individual radius
	void SetPackingFraction(dbl phi);				//!<Set the packing fraction without changing the relative radii
	void SetNumberDensity(dbl density);

	void SetPositions		(const Eigen::VectorXd &t_Positions);		//!<Set particle positions
	void SetPositionsVirtual(const Eigen::VectorXd &t_Positions);
        void SetCNAVirtual(const vector<int> &t_CNAvector, int length);
	void SetParticlePosition		(const dvec &t_Position,int i);		//!<Set an individual particle's position
	void SetParticlePositionVirtual	(const dvec &t_Position,int i);
	
	void MoveParticles(const Eigen::VectorXd &t_Displacement);			//!<Move the particles
	void MoveParticlesVirtual(const Eigen::VectorXd &t_Displacement);

//Functions to get properties of the system
	void GetRadii(Eigen::VectorXd &) const;		//!<Copy the vector of radii
	dbl GetRadius(int i) const;					//!<Get an individual radius
	dbl GetMinimumRadius() const;				//!<Get the minimum radius
	dbl GetMaximumRadius() const;				//!<Get the maximum radius
	dbl GetAverageRadius() const;				//!<Get the average radius

	void GetPositions		(Eigen::VectorXd &) const;		//!<Copy the vector of positions
	void GetPositionsVirtual(Eigen::VectorXd &) const;
	void GetCNAVirtual(vector<int> &) const;
	int GetSingleCNA(int i) const;
	void GetParticlePosition		(dvec &, int i) const;	//!<Get an individual particle's position
	void GetParticlePositionVirtual	(dvec &, int i) const;
	void GetDisplacement		(int i, int j, dvec &displacement) const;	//!<Get the distence between two particles
	void GetDisplacementVirtual	(int i, int j, dvec &displacement) const;

	CBox<Dim>*  GetBox() const;			//!<Return a pointer to the box object
	CPotential* GetPotential() const;	//!<Return a pointer to the potential object

	dbl  GetSphereVolume() const;		//!<Get the combined volume of all the spheres
	dbl  GetPackingFraction() const;	//!<Get the packing fraction
	dbl  GetVolume() const;				//!<Get the volume
	int  GetParticleNumber() const;		//!<Get the number of particles

	void GetMaxDistance(dvec &dist) const;	//!<Get the maximum distance between two particles such that they may be in contact (might be an overestimate)

	void PrintParticles() const;		//!<Print info on the particles to the screen

	friend class CStaticDatabase<Dim>; 
};


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

/*
//Functions to check that the correct dimensions, variables, and attributes are present
//in a netCDF file.
template <int Dim>
bool CStaticState<Dim>::CheckNetCDF(const NcFile &file)
{
	return (file.get_att("System_Populated")!=NULL);
}

template <int Dim>
void CStaticState<Dim>::PopulateNetCDF(NcFile &file)
{
	NcDim *records = file.get_dim("Records");
	if(records==NULL);
		records = file.add_dim("Records");
		
	NcDim *DOF = file.add_dim("System_DOF",Dim*N);
	NcDim *Number = file.add_dim("System_Number",N);
	file.add_dim("System_Dimension",Dim);
	
	file.add_att("System_Populated",1);
	
	file.add_var("System_Positions",ncDouble,records, DOF);
	file.add_var("System_Radii",ncDouble,records, Number);
}
*/

//Constructors/Destructors and copy operators
//
//!!!POTENTIAL MEMORY LEAK?
template <int Dim>
CStaticState<Dim>::CStaticState()
	: N(0),
	  Box(NULL),
	  Potential(NULL)
{
	SetBoxPeriodic();
	SetPotentialHarmonic();
        CNAisvalid = false;
}

template <int Dim>
CStaticState<Dim>::CStaticState(int _N)
	: N(_N),
	  Box(NULL),
	  Potential(NULL)
{
	SetBoxPeriodic();
	SetPotentialHarmonic();
	
	Positions = Eigen::ArrayXd::Zero(Dim*N);
	Radii = Eigen::ArrayXd::Constant(N,0.5);
        CNAisvalid = false;
}

template<int Dim>
CStaticState<Dim>::CStaticState(CStaticState<Dim> const &copy)
	: N(0),
	  Box(NULL),
	  Potential(NULL)
{
	(*this) = copy;
}

template<int Dim>	
CStaticState<Dim>::~CStaticState()
{
	delete Box;
	delete Potential;
}

template<int Dim>
CStaticState<Dim> &CStaticState<Dim>::operator=(CStaticState<Dim> const &copy)
{
	if(this != &copy)
	{
		N = copy.GetParticleNumber();
	
		if(Box!=NULL)		delete Box;
		if(Potential!=NULL)	delete Potential;
		
		Box			= copy.GetBox()->Clone();		//Clone() gives a deep copy
		Potential	= copy.GetPotential()->Clone();	//Clone() gives a deep copy
		
		copy.GetPositionsVirtual(Positions);
		copy.GetRadii(Radii);
		if (copy.CNAisvalid) copy.GetCNAVirtual(CNAvector); // I might want if (CNAvector != NULL) // TAS
                CNAisvalid = copy.CNAisvalid;
	}
	return *this;
}

template<int Dim>
void CStaticState<Dim>::ClearSystem(int _N)
{
	N = _N;
	Positions = Eigen::ArrayXd::Zero(Dim*N);
	Radii = Eigen::ArrayXd::Zero(N);
        CNAisvalid = false; // fill(CNAvector.begin(), CNAvector.end(), 0);
}

template <int Dim>
void CStaticState<Dim>::SetPotential(CPotential *t_Potential)
{
	if(Potential!=NULL) delete Potential;
	Potential = t_Potential->Clone();
}

template<int Dim>
void CStaticState<Dim>::SetPotentialHarmonic(dbl epsilon)
{
	if(Potential!=NULL)	delete Potential;
	Potential = new CHarmonicPotential(epsilon);
}
template<int Dim>
void CStaticState<Dim>::SetPotentialHertzian(dbl epsilon)
{
	if(Potential!=NULL)	delete Potential;
	Potential = new CHertzianPotential(epsilon);
}
template<int Dim>
void CStaticState<Dim>::SetPotentialLJ(dbl epsilon)
{
        if(Potential!=NULL)     delete Potential;
        Potential = new CLJPotential(epsilon);
}

template <int Dim>
void CStaticState<Dim>::SetBox(CBox<Dim> *t_Box)
{
	if(Box!=NULL)	delete Box;
	Box = t_Box->Clone();
}

template<int Dim>
void CStaticState<Dim>::SetBoxPeriodic(int NonPeriodicDim)
{
	if(Box!=NULL)	delete Box;
	switch(NonPeriodicDim)
	{
		case 0:		Box = new CPeriodicBox<Dim,0>();	break;
		case 1:		Box = new CPeriodicBox<Dim,1>();	break;
		case Dim:	Box = new CPeriodicBox<Dim,Dim>();	break;
		default:	assert(false);
	}
}

template<int Dim>
void CStaticState<Dim>::AssumeRectangularSymmetry()
{
	Box->SetSymmetry(Box->RECTANGULAR);
}

//Functions to construct systems 
//Randomly positions particles in the system. 
template<int Dim>
void CStaticState<Dim>::RandomizePositions(MTRand *random)
{
        CNAisvalid = false;
	for(int i = 0; i<N ; i++)
		for(int d = 0 ; d < Dim; d++)
			Positions(Dim*i + d) = random->randExc();
}

template<int Dim>
void CStaticState<Dim>::RandomizePositions(long seed)
{
        CNAisvalid = false;
	MTRand *random = new MTRand(seed);
	RandomizePositions(random);
	delete random;
}

template<int Dim>
void CStaticState<Dim>::SetSquareLattice()
{
	int L;
	IsPerfectRoot(N,Dim,L);
	if(POW2(L)!=N)
	{
		printf("ERROR: L^2 != N, with L=%i, N=%i is not a perfect root\n", L, N);
		assert(false);
	}
	
	SetSquareLattice(L,L);
}

template<int Dim>
void CStaticState<Dim>::SetSquareLattice(int Lx, int Ly)
{
	assert(Lx*Ly==N);
	if(Dim!=2)
	{
		printf("ERROR: can only set a square lattice when Dim = 2\n");
		assert(false);
	}

	dvec spacing;
	dvec offset;
	spacing[0] = 1./((dbl)Lx);
	spacing[1] = 1./((dbl)Ly);
	for(int dd=0; dd<Dim; ++dd) offset[dd] = 0.25*spacing[dd];

	int ii;
	for(int y=0; y<Ly; ++y) 
		for(int x=0; x<Lx; ++x) 
		{    
			ii = Lx*y+x;
			Positions(Dim*ii  ) = spacing[0]*((dbl)x) + offset[0];
			Positions(Dim*ii+1) = spacing[1]*((dbl)y) + offset[1];
		}
}


template <int Dim>
void CStaticState<Dim>::SetRadiiMono()
{
	const dbl sigma = 1.;
	Radii = Eigen::VectorXd::Constant(N, sigma/2.);
}

/**
 *	@param[in] Fracsmall The fraction of small particles
 *	@param[in] SizeRatio The ratio of largest radii to the smallest
 */
template <int Dim>
void CStaticState<Dim>::SetRadiiBi(dbl FracSmall, dbl SizeRatio)
{
	const dbl sigma = 1.;
	assert(Radii.size() == N);
	int Nsmall = FracSmall*N;
	dbl size1 = 1.;
	dbl size2 = SizeRatio;
	dbl avgD = 2.*(Nsmall*size1 + (N-Nsmall)*size2)/((dbl)N);
	size1*=sigma/avgD;
	size2*=sigma/avgD;

	for(int i=0; i<Nsmall; ++i)
		Radii[i] = size1;
	for(int i=Nsmall; i<N; ++i)
		Radii[i] = size2;
}

template <int Dim>
void CStaticState<Dim>::SetRadiiBiGaussian(dbl s1, dbl s2, dbl FracSmall, dbl SizeRatio, long seed)
{
	const dbl sigma = 1.;
	assert(Radii.size() == N);
	int Nsmall = FracSmall*N;
	dbl r1 = 1.0;
	dbl r2 = r1*SizeRatio;
	dbl tempr;

	//Initialize the random number generator
	RNG_taus R(seed);

	//Set Nsmall from a Gaussian distribution with mean r1 and standard deviation s1
	for(int i=0; i<Nsmall; ++i)
	{
		do{
			tempr = s1*random_normal(R)+r1;
		}while(tempr <= 0.);
		Radii[i] = tempr;
	}

	//Set N-Nsmall from a Gaussian distribution with mean r2 and standard deviation s2
	for(int i=Nsmall; i<N; ++i)
	{
		do{
			tempr = s2*random_normal(R)+r2;
		}while(tempr <= 0.);
		Radii[i] = tempr;
	}

	dbl avgD = 2.*Radii.mean();
	Radii *= sigma/avgD;
}

/**
 *	@param[in] SizeRatio The ratio of largest radii to the smallest
 */
template <int Dim>
void CStaticState<Dim>::SetRadiiPolyUniform(dbl SizeRatio)
{
	const dbl sigma = 1.;
	dbl RelSize;
	for(int i=0; i<N; ++i)
	{
		RelSize = ((dbl)1.) + (((dbl)i)/((dbl)(N-1)))*SizeRatio;
		Radii[i] = RelSize;
	}
	dbl avgD = 2.*Radii.mean();
	Radii *= sigma/avgD;
}


		
//Functions to set properties of the system
template <int Dim>
void CStaticState<Dim>::SetRadii(const Eigen::VectorXd &t_Radii)
{
	//set the radii vector  
	if(t_Radii.rows()!=N)
		throw(CException("CStaticState<Dim>::SetRadii","Radii vector has an inconsistent size.")); 

	Radii = t_Radii;
}

template <int Dim>
void CStaticState<Dim>::SetRadius(int i, dbl r)
{
	//set the radii vector  
	if(i>N)
		throw(CException("CStaticState<Dim>::SetRadius","Attempting to set the radius of a particle that doesn't exist.")); 

	Radii(i) = r;
}

template <int Dim>
void CStaticState<Dim>::SetPositions(const Eigen::VectorXd &t_Positions)
{
	//set the position vector
	if(t_Positions.rows()!=Dim*N)
		throw(CException("CStaticState<Dim>::SetPositions","Positions vector has an inconsistent size.")); 
	
	Positions = t_Positions;
        CNAisvalid = false;

	Box->InverseTransform(Positions);
}

// Following the template of SetPositionsVirtual, this function sets a CNA value for all particles
template <int Dim>
void CStaticState<Dim>::SetCNAVirtual(const vector<int> &t_CNAvector, int length)
{
	//set the position vector
	if(length!=N)
		throw(CException("CStaticState<Dim>::SetCNAVirtual","Reported CNA length doesnt match number of particles of State.")); 
		
	CNAvector = t_CNAvector;
}

template <int Dim>
void CStaticState<Dim>::SetPositionsVirtual(const Eigen::VectorXd &t_Positions)
{
	//set the position vector
	if(t_Positions.rows()!=Dim*N)
		throw(CException("CStaticState<Dim>::SetPositionsVirtual","Positions vector has an inconsistent size.")); 
		
        CNAisvalid = false;
	Positions = t_Positions;
}

template <int Dim>
void CStaticState<Dim>::SetParticlePosition(const Eigen::Matrix<dbl,Dim,1> &t_Positions,int i)
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::SetParticlePosition","Attempting to set the position of a particle that does not exist."));
		
	Eigen::Matrix<dbl,Dim,1> new_pos = t_Positions;
	Box->InverseTransform(new_pos);
	Positions.segment<Dim>(Dim*i) = new_pos;
        CNAisvalid = false;
}

template <int Dim>
void CStaticState<Dim>::SetParticlePositionVirtual(const Eigen::Matrix<dbl,Dim,1> &t_Positions,int i)
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::SetParticlePositionVirtual","Attempting to set the position of a particle that does not exist."));
		
	Positions.segment<Dim>(Dim*i) = t_Positions;
        CNAisvalid = false;
}

template<int Dim>
void CStaticState<Dim>::MoveParticles(const Eigen::VectorXd &t_Displacement)
{
	if(t_Displacement.rows()!=Dim*N)
		throw(CException("CStaticState<Dim>::MoveParticles","Displacement vector has an inconsistent size.")); 
		
//	Box->InverseTransform(t_Displacement);
//	Box->MoveParticles(Positions,t_Displacement);
	Box->InverseTransformAndMove(Positions,t_Displacement);
        CNAisvalid = false;
}

template<int Dim>
void CStaticState<Dim>::MoveParticlesVirtual(const Eigen::VectorXd &t_Displacement)
{
	if(t_Displacement.rows()!=Dim*N)
		throw(CException("CStaticState<Dim>::MoveParticlesVirtual","Displacement vector has an inconsistent size.")); 
		
	Box->MoveParticles(Positions,t_Displacement);
        CNAisvalid = false;
}

template <int Dim>
void CStaticState<Dim>::SetPackingFraction(dbl phi)
{
	Box->SetVolume(GetSphereVolume()/phi);
	assert( fabs(GetPackingFraction()-phi) < 1e-10 );
}

template <int Dim>
void CStaticState<Dim>::SetNumberDensity(dbl den)
{
        Box->SetVolume(N/den);
}

//Functions to get properties of the system
template <int Dim>
void CStaticState<Dim>::GetRadii(Eigen::VectorXd &t_Radii) const
{
	t_Radii = Radii;
}

template <int Dim>
dbl CStaticState<Dim>::GetRadius(int i) const
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::GetRadius","Attempting to get the radius of a particle that does not exist."));
	return Radii(i);
}

template <int Dim>
dbl CStaticState<Dim>::GetMinimumRadius() const
{
	return Radii.minCoeff();
}

template <int Dim>
dbl CStaticState<Dim>::GetMaximumRadius() const
{
	return Radii.maxCoeff();
}

template <int Dim>
dbl CStaticState<Dim>::GetAverageRadius() const
{
	return Radii.mean();
}

template<int Dim>
void CStaticState<Dim>::GetPositions(Eigen::VectorXd &t_Positions) const
{
	t_Positions = Positions;
	Box->Transform(t_Positions);
}

template<int Dim>
void CStaticState<Dim>::GetPositionsVirtual(Eigen::VectorXd &t_Positions) const
{
	t_Positions = Positions;
}

template<int Dim>
void CStaticState<Dim>::GetCNAVirtual(vector<int> &t_CNAvector) const
{
	if (CNAisvalid) t_CNAvector = CNAvector;
}

template <int Dim>
int CStaticState<Dim>::GetSingleCNA(int i) const
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::GetParticlePosition","Attempting to get the CNA of a particle that does not exist."));
        return CNAvector[i];
}

template <int Dim>
void CStaticState<Dim>::GetParticlePosition(Eigen::Matrix<dbl,Dim,1> &t_Position, int i) const
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::GetParticlePosition","Attempting to get the position of a particle that does not exist."));

	t_Position = Positions.segment<Dim>(Dim*i);
	Box->Transform(t_Position);
}

template <int Dim>
void CStaticState<Dim>::GetParticlePositionVirtual(Eigen::Matrix<dbl,Dim,1> &t_Position, int i) const
{
	if(i<0||i>N)
		throw(CException("CStaticState<Dim>::GetParticlePositionVirtual","Attempting to set the position of a particle that does not exist."));

	t_Position = Positions.segment<Dim>(Dim*i);
}

template <int Dim>
void CStaticState<Dim>::GetDisplacement(int i,int j,Eigen::Matrix<dbl,Dim,1> &displacement) const
{
	if(i<0||i>N||j<0||j>N)
		throw(CException("CStaticState<Dim>::GetDisplacement","Attempting to get the displacement between particles that do not exist."));

	Box->MinimumDisplacement(Positions.segment<Dim>(Dim*i),Positions.segment<Dim>(Dim*j),displacement);
	
	Box->Transform(displacement);
}

template <int Dim>
void CStaticState<Dim>::GetDisplacementVirtual(int i,int j,Eigen::Matrix<dbl,Dim,1> &displacement) const
{
	if(i<0||i>N||j<0||j>N)
		throw(CException("CStaticState<Dim>::GetDisplacement","Attempting to get the displacement between particles that do not exist."));

	Box->MinimumDisplacement(Positions.segment<Dim>(Dim*i),Positions.segment<Dim>(Dim*j),displacement);
}

template <int Dim>
CBox<Dim> *CStaticState<Dim>::GetBox() const
{
	return Box;
}

template <int Dim>
CPotential *CStaticState<Dim>::GetPotential() const
{
	return Potential;
}

template <int Dim>
dbl CStaticState<Dim>::GetSphereVolume() const
{
	dbl SphereVolume = 0.;
	for(int i=0; i<Radii.size(); ++i)
		SphereVolume += std::pow(Radii[i],Dim);
	SphereVolume *= nSphere_Vn(Dim);
	return SphereVolume;
}

template <int Dim>
dbl CStaticState<Dim>::GetPackingFraction() const
{
	return GetSphereVolume()/GetVolume();
}

template <int Dim>
dbl CStaticState<Dim>::GetVolume() const
{
	return Box->CalculateVolume();
}
	
template <int Dim>
int CStaticState<Dim>::GetParticleNumber() const
{
	return N;
}

template <int Dim>
//Eigen::Matrix<dbl,Dim,1> CStaticState<Dim>::GetMaxDistance() const
void CStaticState<Dim>::GetMaxDistance(dvec &dist) const
{
	Box->GetMaxTransformedDistance(dist);
	dist *= 2.*Radii.maxCoeff()*Potential->ComputeSupport();
}
	
template <int Dim>
void CStaticState<Dim>::PrintParticles() const
{
	Eigen::VectorXd RealPos = Positions;
	Box->Transform(RealPos);
	printf("Printing %i particles:\n", N);
	for(int i=0; i<N; ++i)
	{
		printf("  sphere:% 5i   Position: ", i); 
		for(int dd=0; dd<Dim; dd++) printf("% 20.14f  ", RealPos[Dim*i+dd]);
		printf("  Radius: %16.14f", Radii[i]);
		printf("\n");
	}
	dmat Trans;
	Box->GetTransformation(Trans);
	cout << "Transformtion matrix:\n" << Trans << endl;
}

}

#endif
