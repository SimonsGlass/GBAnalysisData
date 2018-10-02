#ifndef GRID

#define GRID

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
//#include "../Resources/MersenneTwister.h"


namespace LiuJamming
{

//using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


//!Grid class to efficiently calculate neighbors.
/**This class maintains a spatial partition of particle systems so that
 * only local regions of a given particle need to be searched for neighbors.
 */
template <int Dim>
class CGrid
{
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
public:

	//!iterator class to allow the user to iterate over potential neighbors.
	class iterator {
	private:
		list<int>::iterator CurrentCell;
		list<int>::iterator End;
		
		int CurrentParticle;
		
		int *OccupancyList;
		int *CellList;
		
	public:
		//!Constructor
		iterator(list<int> &dual,int *occ_list, int *cell_list)
		{
			CurrentCell = dual.begin();
			End = dual.end();
			OccupancyList = occ_list;
			CellList = cell_list;
			CurrentParticle = CellList[(*CurrentCell)];
			IncrementCellIfNecessary();
		};

		//!Increment operator
		iterator &operator ++() 
		{
			if(CurrentParticle!=-1)
			{
				CurrentParticle = OccupancyList[CurrentParticle];
				IncrementCellIfNecessary();
			}
			return (*this);
		};

		//!Increment operator
		void operator ++(int dummy)
		{
			if(CurrentParticle!=-1)
			{
				CurrentParticle = OccupancyList[CurrentParticle];
				IncrementCellIfNecessary();
			}
		};

private:
		//!If CurrentParticle==-1, move on to the next cell.
		void IncrementCellIfNecessary()
		{
			while(CurrentParticle==-1&&CurrentCell!=End)
			{
				CurrentCell++;
				if(CurrentCell!=End)
					CurrentParticle = CellList[(*CurrentCell)];
			}
		};

public:		
		//!Dereference operator
		int operator *()
		{
			return CurrentParticle;
		};
	};

private:
	CStaticState<Dim> *State;	//!<Pointer to a CStaticState<Dim>.
	int N;						//!<Number of particles.

	dvec CellSize;				//!<Cell size given by twice the maximum radius.
	
	int N_Cells[Dim];			//!<The number of cells in each dimension.
	int TotalCells;				//!<Total number of cells.
	
	int *OccupancyList;			//!<From one particle, points to the next particle in the same cell (or -1 if it is the last particle).
	int *CellList;				//!<Pointer to the first particle in each cell.
	
	list<int> *DualList;		//!<The list of connections between cells.
	
	//Flags
	bool DualListUpdateNecessary;	//!<A boolean to indicate whether the DualList might have changed.
	
public:
//! @name Constructors and operators
///@{
	CGrid(CStaticState<Dim> *_s);							//!<Primary constructor.
	CGrid(const CGrid &copy);								//!<Copy constructor.
	~CGrid();												//!<Destructor.
	const CGrid<Dim> &operator=(const CGrid<Dim> &copy);	//!<Copy operator.
	bool ReallocationNecessary;             //!<A boolean to indicate whether we have to reallocate the grid.
	void ClearLists();										//!<Deallocate memory for the lists.
	void SetN(int _N);										//!<Set N and allocate memory for OccupancyList.
	void SetState(CStaticState<Dim> *s);					//!<Set the State pointer.

///@}


//! @name Grid manipulation
///@{
	void Allocate();		//!<Allocate memory for the grid.
	void Construct();		//!<Place particles in the grid.
	void UpdateDualList();	//!<Update the list of cell connections.

///@}

//! @name Functions to access the grid
///@{
	iterator GetParticleIterator(int i);				//!<Return an iterator to loop over potential neighbors.

///@}
	
//! @name Misc.
///@{
 	void CellToCoordinates(int i, dvec &coordinates);	//!<Convert a cell index into a dvec of cell coordinates.
	int CoordinatesToCell(const dvec &coordinates);		//!<Convert a dvec of cell coordinates into a cell index.
	dvec ComputeCellSize();								//!<Function to compute the grid cell size.
	void PrintGrid() const;								//!<Print the grid to stdout.

///@}

};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//Constructor and copy operators
template <int Dim>
CGrid<Dim>::CGrid(CStaticState<Dim> *s) : State(s), N(0), OccupancyList(NULL), CellList(NULL), DualList(NULL)
{
	ClearLists();
	SetN(State->GetParticleNumber());
	for(int i=0; i<N; ++i) OccupancyList[i] = -1;
//	CellList = NULL;
//	DualList = NULL;
}

/*
template <int Dim>
CGrid<Dim>::CGrid(const CGrid &copy) : State(copy.State), N(0)
{
	ClearLists();
	SetN(State->GetParticleNumber());
	for(int i=0; i<N; ++i) OccupancyList[i] = -1;
//	CellList = NULL;
//	DualList = NULL;
}
*/

template <int Dim>
CGrid<Dim>::~CGrid()
{
	ClearLists();
}

/**
 * WARNING: this is not implemented correctly.
 */
template <int Dim>
const CGrid<Dim> &CGrid<Dim>::operator=(const CGrid<Dim> &copy)
{
	if(this != &copy)
	{
		printf("WARNING: copying CGrid is not implemented correctly.\n");
		assert(false);
		ClearLists();
		State = copy.State;
		SetN(State->GetParticleNumber());
		for(int i=0; i<N; ++i) OccupancyList[i] = -1;
//		CellList = NULL;
//		DualList = NULL;
	}
	return *this;
}

template <int Dim>
void CGrid<Dim>::SetState(CStaticState<Dim> *s)
{
	State = s;
}

template <int Dim>
void CGrid<Dim>::ClearLists()
{
	if(OccupancyList	!= NULL) delete[] OccupancyList;
	if(CellList			!= NULL) delete[] CellList;
	if(DualList			!= NULL) delete[] DualList;
	OccupancyList = NULL;
	CellList = NULL;
	DualList = NULL;
}

template <int Dim>
void CGrid<Dim>::SetN(int _N)
{
	if(N!=_N)
	{
		if(OccupancyList!=NULL)
		{
			printf("about to delete[] OccupancyList.\n"); fflush(stdout);
			delete[] OccupancyList;
		}

		N = _N;
		OccupancyList = new int[N];

		ReallocationNecessary = true;
		DualListUpdateNecessary = true;
	}
}


//Functions to construct the grid
template <int Dim>
void CGrid<Dim>::Allocate()
{
	dvec MaxDistance;
	State->GetMaxDistance(MaxDistance);
	//Compute the number of grid elements in each direction. Additionally,
	//figure out the total number of cells as the product of all of these.
	TotalCells = 1;
	for(int i = 0 ; i < Dim ; i++){
		//N_Cells[i] = (int)ceil(1.0/MaxDistance(i));
		N_Cells[i] = (int)floor(1.0/MaxDistance(i));
		assert(N_Cells[i] >= 3);
		TotalCells*=N_Cells[i];
		CellSize[i] = 1./((dbl)N_Cells[i]);
		//printf("CellSize[%i] = %e\n", i, CellSize[i]);
	}

	if(!(CellList==NULL))
		delete[] CellList;

	if(!(DualList==NULL))
		delete[] DualList;

	CellList = new int[TotalCells];
	for(int i = 0; i<TotalCells; i++)
		CellList[i] = -1;
	
	DualList = new list<int>[TotalCells];

	DualListUpdateNecessary = true;
	ReallocationNecessary = false;
}

template <int Dim>
void CGrid<Dim>::Construct()
{
	//Reallocate if needed
	if(ReallocationNecessary)
		Allocate();

	//This will be taken care of later
//	for(int i = 0 ; i<N ; i++)
//		OccupancyList[i] = -1;

	//reset CellList
	for(int i = 0 ; i< TotalCells ;i++)
		CellList[i] = -1;  //This is necessary because a cell might not have any particles in it.
	

	//compute the transformed positions of the particles.
//	Eigen::VectorXd Positions;
//	State->GetPositionsVirtual(Positions);
	
	//First put the particles in the cells.
	dvec PosTemp;
	int Cell_Index;
	for(int i = 0 ; i<N ; i++)
	{
		//Calculate the cell that the particle is in.
		//We should be able to do this without copying the particle positions
		State->GetParticlePositionVirtual(PosTemp, i);
		Cell_Index = CoordinatesToCell(PosTemp);
		assert(Cell_Index >= 0);
		assert(Cell_Index <  TotalCells);
		
		//Add the particle to that cell
		OccupancyList[i] = CellList[Cell_Index];
		CellList[Cell_Index] = i;
	}

	//Now construct the dual list
	if(DualListUpdateNecessary)
		UpdateDualList();
	/*
	//go through each pair of cells and find their coordinates
	dvec Displacement;
	for(int i = 0 ; i < TotalCells ; i++)
	{
		dvec CoordinateI;
		CellToCoordinates(i,CoordinateI);
		for(int j = i ; j < TotalCells ; j++)
		{
			dvec CoordinateJ;
			CellToCoordinates(j,CoordinateJ);
			
			State->GetBox()->MinimumDisplacement(CoordinateI,CoordinateJ,Displacement);
			
			//go through the dimensions and if any dimension is closer than CellSize+\epsilon apart
			//add the cells to each others dual
			bool add = true;
			for(int k = 0 ; k < Dim ; k++)
				if(abs(Displacement(k))>CellSize(k)+CellSize(k)/100.0)	
					add = false;
	
			if(add)
			{
				//cout << "Displacement = " << max_displacement << " : CellSize = " << CellSize(0) << endl;
				DualList[i].push_back(j);
				if(j!=i)
					DualList[j].push_back(i);
			}
		}
	}
	*/
}

template<int Dim>
void CGrid<Dim>::UpdateDualList()
{
	//reset the lists
	for(int i = 0 ; i< TotalCells ;i++)
	{
		DualList[i].erase(DualList[i].begin(),DualList[i].end());
	}

	//check that there are enough cells
	vector<int> PeriodicDims;
	State->GetBox()->GetPeriodicDimensions(PeriodicDims);
	for(typename vector<int>::iterator it=PeriodicDims.begin(); it!=PeriodicDims.end(); ++it)
		assert(N_Cells[(*it)] >= 3); //this makes sure that you can never touch multiple images of the same particle

	//go through each pair of cells and find their coordinates
	dvec Displacement;
	for(int i = 0 ; i < TotalCells ; i++)
	{
		dvec CoordinateI;
		CellToCoordinates(i,CoordinateI);
		for(int j = i ; j < TotalCells ; j++)
		{
			dvec CoordinateJ;
			CellToCoordinates(j,CoordinateJ);
			
			State->GetBox()->MinimumDisplacement(CoordinateI,CoordinateJ,Displacement);
			
			//go through the dimensions and if any dimension is closer than CellSize+\epsilon apart
			//add the cells to each others dual
			bool add = true;
			for(int k = 0 ; k < Dim ; k++)
				if(abs(Displacement(k))>CellSize(k)+CellSize(k)/100.0)	
					add = false;
	
			if(add)
			{
				//cout << "Displacement = " << max_displacement << " : CellSize = " << CellSize(0) << endl;
				DualList[i].push_back(j);
				if(j!=i)
					DualList[j].push_back(i);
			}
		}
	}
//	for(int i=0; i<TotalCells; ++i)
//		printf("DualList[%i].size() = %i\n", i, (int)DualList[i].size());
	DualListUpdateNecessary = false;
}

//Functions to access the grid
template<int Dim>
void CGrid<Dim>::CellToCoordinates(int i,dvec &Coordinate)
{
	int Prod = 1;
	int LProd = 1;
	dbl Sub = 0;
	for(int k = 0 ; k < Dim ; k++){
		Prod*=N_Cells[k];
		Coordinate(k) = (i%Prod) - Sub;
		Sub+=Coordinate(k);
		Coordinate(k)/=LProd;
		LProd*=N_Cells[k];
	}

	for(int i = 0; i< Dim ; i++)
		Coordinate(i)*=CellSize(i);
}
	
template <int Dim>
int CGrid<Dim>::CoordinatesToCell(const dvec &coordinates)
{
	int Cell_Index = 0;
	int Prod = 1;
	for(int j = 0 ; j<Dim ; j++)
	{
		Cell_Index += Prod*((int)floor(coordinates(j)/CellSize(j)));
		Prod*=N_Cells[j];
	}
	
	return Cell_Index;
}

/*
template <int Dim>
void CGrid<Dim>::CellIndexToCellCoordinates(int i, ivec &Coordinate)
{
	int Prod = 1;
	int LProd = 1;
	int Sub = 0;
	for(int k = 0 ; k < Dim ; k++)
	{
		Prod*=N_Cells[k];
		Coordinate(k) = (i%Prod) - Sub;
		Sub+=Coordinate(k);
		Coordinate(k)/=LProd;
		LProd*=N_Cells[k];
	}
}

template <int Dim>
int  CGrid<Dim>::CellCoordinatesToCellIndex(const ivec &Coordinate)
{
	int Cell_Index = 0;
	int Prod = 1;
//	int Coord_temp;
	for(int j=0; j<Dim; j++)
	{
//		Coord_temp = Coordinate[j];
//		if(Coord_temp < 0) Coord_temp += N_cells[j];
//		if(Coord_temp >= N_cells[j]) Coord_temp -= N_cells[j];
		Cell_Index += Prod*Coordinate[j];
		Prod *= N_cells[j];
	}
	return Cell_Index;
}
*/

template <int Dim>
typename CGrid<Dim>::iterator CGrid<Dim>::GetParticleIterator(int i)
{
	dvec pos;
	State->GetParticlePositionVirtual(pos,i);
	int j = CoordinatesToCell(pos);
	typename CGrid<Dim>::iterator ret(DualList[j],OccupancyList,CellList);
	return ret;
}

template <int Dim>
void CGrid<Dim>::PrintGrid() const
{
	printf("Printing Grid:\n");
	printf("\tN=%i\n", N);
	printf("\tTotalCells = %i\n", TotalCells);
	printf("OccupancyList:\n");
	for(int i=0; i<N; ++i) printf("\t%5i --> %5i\n", i, OccupancyList[i]);
	printf("CellList:\n");
	for(int i=0; i<TotalCells; ++i) printf("\t head of cell %5i = %5i\n", i, CellList[i]);
	printf("DualList:\n");
	for(int i=0; i<TotalCells; ++i)
	{
		printf("\t%5i --> ", i);
		for(typename list<int>::iterator it = DualList[i].begin(); it!=DualList[i].end(); ++it)
			printf("%i ", (*it));
		printf("\n");
	}

	printf("Cell Details:\n");
	for(int i=0; i<TotalCells; ++i)
	{
		printf("\tcell %5i --> ", i);
		int curr = CellList[i];
		while(curr != -1)
		{
			printf("%i ", curr);
			curr = OccupancyList[curr];
		}
		printf("\n");
	}

}

/*
template <int Dim>
Eigen::Matrix<dbl,Dim,1> CGrid<Dim>::ComputeCellSize()
{
	dvec MaxDistance = State->GetMaxDistance();

	return MaxDistance;
}
*/

/*
//!!!!!!!!!! I THINK THIS METHOD IS BUGGY
//
//State.Radii is in REAL units, while State.Positions is in boxed coordinates.
//There should be a method in CStaticState that returns the maximum possible distance in BOXED coordinates that 
//		neighboring particles could be.
template <int Dim>
Eigen::Matrix<dbl,Dim,1> CGrid<Dim>::ComputeCellSize()
{
	dbl MaximumRadius = 0.0;
	Eigen::VectorXd Rads;
	State->GetRadii(Rads);

	for(int i = 0; i < State->GetParticleNumber() ; i++)
		if(Rads(i)>MaximumRadius)
			MaximumRadius = Rads(i);

	dbl Scale = 2*MaximumRadius;
	dbl N = floor(1.0/Scale);
	Scale = 1.0/N;

	dvec ret = dvec::Constant(Scale);
	
	return ret;
}
*/

}

#endif
