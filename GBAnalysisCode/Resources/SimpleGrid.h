#ifndef SIMPLE_GRID


#define SIMPLE GRID


namespace LiuJamming
{

template <int Dim>
class CSimpleNeighbor {
public:
	Eigen::Matrix<double,Dim,1> Displacement;
	int i,j;

	CSimpleNeighbor() {}
	CSimpleNeighbor(int _i, int _j,  Eigen::Matrix<double,Dim,1> _Displacement) { i = _i ; j = _j; Displacement = _Displacement;} 
};

template <int Dim>
class CSimpleGrid
{
public:
	list<int> *GridParticles; 
	list<int> *CellNeighbors;

	list<CSimpleNeighbor<Dim> > *NeighborList;

	Eigen::Matrix<int,Dim,1> Size; // This might be the grid size,  Nx Ny Nz in three dimensions
	Eigen::Matrix<double,Dim,1> Box; // What will Box be used for? Does this assume cubic?

	double Spacing;

public:
        // TAS Getting an error that would avoided if particles were wrapped before adding to the box
	inline int IndexFromCoordinate(Eigen::Matrix<int,Dim,1> &coord) {return coord(0) + coord(1)*Size(0) + (Dim == 3 ? coord(2)*Size(0)*Size(1) : 0); } 
	//inline int IndexFromCoordinate(Eigen::Matrix<double,Dim,1> &coord) {return floor(coord(0)/Box(0)*Size(0)) + floor(coord(1)/Box(1)*Size(1))*Size(0) + (Dim == 3 ? floor(coord(2)/Box(2)*Size(2))*Size(0)*Size(1) : 0); }

        inline int IndexFromCoordinate(Eigen::Matrix<double,Dim,1> &coord) {return floor(coord(0)/Spacing) + floor(coord(1)/Spacing)*Size(0) + (Dim == 3 ? floor(coord(2)/Spacing)*Size(0)*Size(1) : 0); }

	
	CSimpleGrid(CStaticState<Dim> &State,double Spacing = 2.5);
	~CSimpleGrid();

	//Repopulates the grid. Note: assumes the box size has not changed.
	void Repopulate(CStaticState<Dim> &State);
	//Constructs a neighbor list for a set of particles
	void PopulateNeighborList(CStaticState<Dim> &State, double Cutoff, list<int> Particles, list<list<CSimpleNeighbor<Dim> > > &Neighbors);
	
	//Construct a neighbor list
	void ConstructNeighborList(CStaticState<Dim> &State, double Cutoff);

	//Get neighboring cells and particles in a cell	
	const list<int> &GetCellParticles(int Index);
	const list<int> &GetNeighborCells(Eigen::Matrix<double,Dim,1> &Position);

	//Get the neighbors of a particle if the neighbor list has been constructed
	const list<CSimpleNeighbor<Dim> > &GetNeighbors(int particle);
};


template <int Dim>
CSimpleGrid<Dim>::CSimpleGrid(CStaticState<Dim> &State,double _Spacing)
{

        cout << "Making a new simple grid." << endl;
	Spacing = _Spacing;
	Eigen::Matrix<double,Dim,Dim> Transformation;
	State.GetBox()->GetTransformation(Transformation);
	Box = Transformation.diagonal(); // TAS This assumes orthorhombic 
        // Box is 1.0 but probably should have been read-in
        

        cout << "Found the transformation." << endl; cout.flush();

        // Prod is probably Nx * Ny * Nz for 3D
        printf("Spacing is %f \n", Spacing);fflush(stdout);
        printf("Box is %f %f %f \n", Box(0), Box(1), Box(2));fflush(stdout);
        for(int idim = 0; idim < Dim ; idim++) 
                if (Spacing > Box(idim)) 
                        cerr << "Grid spacing " << Spacing << " greater than box size " << Box(idim) << endl;

	int Prod = 1;
	for(int d = 0 ; d < Dim ; d++)
	{
		Size(d) = ceil(Box(d)/Spacing);
		Prod *= Size(d);
	}

        cout << "About to make a list of length " << Prod << endl; cout.flush();
	GridParticles = new list<int>[Prod];
	CellNeighbors = new list<int>[Prod];
        cout << "Made it.  Prod is " << Prod << endl; cout.flush();
        cout << "Size" << Size(0) << " " << Size(1) << " " << Size(2) << endl;

	//First populate the neighbor list
        // Coord is the coordinate within the 3D grid: (i,j,k)
	Eigen::Matrix<int,Dim,1> Coord = Eigen::Array<int,Dim,1>::Zero();
	while(Coord(Dim-1) < Size(Dim-1))
	{
		int Index = IndexFromCoordinate(Coord); // TAS this looks like it would work
		
                // Coord_It is a helper function, equal to -2 if near the left edge or -1 elsewhere TAS
		Eigen::Matrix<int,Dim,1> Coord_It;
		for(int d = 0 ; d < Dim ; d++)
			Coord_It(d) = Coord(d) < 2 ? -2 : -1;

                // The helper function in the last dimension remains 1 or less or else at least under 3 if also the coord in the last dimensino is at the right edge
		while(Coord_It(Dim-1) < 2 || (Coord(Dim-1) >= Size(Dim-1)-2 && Coord_It(Dim-1) < 3))
		{
                        // Coord_2 is the coordinate of the neighboring cell.  It is wrapped to the box
			Eigen::Matrix<int,Dim,1> Coord_2 = Coord + Coord_It;
			for(int d = 0 ; d < Dim ; d++){
				if(Coord_2(d) < 0) Coord_2(d) += Size(d);
				if(Coord_2(d) >= Size(d)) Coord_2(d) -= Size(d);
			}

			int Index_2 = IndexFromCoordinate(Coord_2);
                        if (Index_2 < 0) { printf("index2 %i   Coord_2 %i %i %i \n", Index_2, Coord_2(0), Coord_2(1), Coord_2(2));fflush(stdout); }
			
                        // Add index of the neighJ to the CellNeigh[IndexNeighI]
			CellNeighbors[Index].push_back(Index_2);
	
			Coord_It(0)++;
			int t_d = 0;
			while(t_d < Dim-1 && (!(Coord_It(t_d) < 2 || (Coord(t_d) >= Size(t_d)-2 && Coord_It(t_d)<3))))
			{
				Coord_It(t_d) = Coord(t_d)<2 ? -2 : -1;
				Coord_It(t_d+1)++;
				t_d++;
			}
		}

		//cout << CellNeighbors[Index].size() << endl;
/*
		Prod = 1;
		for(int d = 0 ; d < Dim ; d++) {
			Neighbors[Index].push_back(Coord(d)+1 < Size(d) ? Index + Prod : Index - (Size(d) - 1)*Prod);
			Neighbors[Index].push_back(Coord(d)-1 >= 0 ? Index - Prod : Index + (Size(d) - 1)*Prod);
			Prod *= Size(d);
		}

*/
		Coord(0)++;
		int t_d =0;
		while(t_d < Dim-1 && Coord(t_d) >= Size(t_d))
		{
			Coord(t_d) = 0;
			Coord(t_d+1)++;
			t_d++;
		}
	}
	
	//Now populate the grid 
	Eigen::Matrix<double,Dim,1> Position;
        for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		State.GetParticlePosition(Position, i);
		if ((IndexFromCoordinate(Position) < 0) || (IndexFromCoordinate(Position) >= Size(2)*Size(1)*Size(0))) printf("problem particle i %i at index %i (%f %f %f) \n", i, IndexFromCoordinate(Position), Position(0), Position(1), Position(2));
		GridParticles[IndexFromCoordinate(Position)].push_back(i);
		
	}

	NeighborList = NULL;	

}

template <int Dim>
CSimpleGrid<Dim>::~CSimpleGrid()
{
	delete []GridParticles;
	delete []CellNeighbors;
	
	if(NeighborList)
		delete[] NeighborList;
}

template <int Dim>
void CSimpleGrid<Dim>::ConstructNeighborList(CStaticState<Dim> &State, double Cutoff)
{
	if(NeighborList)
	{
		if(NeighborList->size() == State.GetParticleNumber())
			for(int i = 0 ; i < NeighborList->size() ; i++)
				NeighborList[i].clear();
		else{
			delete []NeighborList;
			NeighborList = new list<CSimpleNeighbor<Dim> >[State.GetParticleNumber()];
		}
	}else{
		NeighborList = new list<CSimpleNeighbor<Dim> >[State.GetParticleNumber()];
	}
	
	Eigen::Matrix<double,Dim,1> Position;
	Eigen::Matrix<double,Dim,1> Displacement;
	for(int i = 0 ; i < State.GetParticleNumber() ; i++)
	{
		State.GetParticlePosition(Position,i);
		
		list<int> cells = GetNeighborCells(Position);


		for(list<int>::iterator cells_it = cells.begin() ; cells_it != cells.end() ; cells_it++)
		{
			list<int> particles = GetCellParticles(*cells_it);
			
			for(list<int>::iterator particle_it =  particles.begin(); particle_it!=particles.end() ; particle_it ++){
				State.GetDisplacement(i,*particle_it,Displacement);
				
				if(i != *particle_it && Displacement.squaredNorm() < Cutoff*Cutoff){
					
					NeighborList[i].push_back(CSimpleNeighbor<Dim>(i,*particle_it,Displacement));
					
				}
			}
		}

	}
}

template <int Dim>
void CSimpleGrid<Dim>::PopulateNeighborList(CStaticState<Dim> &State, double Cutoff, list<int> Particles, list<list<CSimpleNeighbor<Dim> > > &Neighbors)
{

	Eigen::Matrix<double,Dim,1> Position;
        Eigen::Matrix<double,Dim,1> Displacement;
        for(list<int>::iterator it = Particles.begin() ; it != Particles.end() ; it++)
        {
                State.GetParticlePosition(Position,*it);
		list<CSimpleNeighbor<Dim> > TNeighbor;

                list<int> cells = GetNeighborCells(Position);
                // evidently cells includes the number -3.  If imported non-scaled coordinates.

                for(list<int>::iterator cells_it = cells.begin() ; cells_it != cells.end() ; cells_it++)
                {
                        list<int> particles = GetCellParticles(*cells_it);

                        for(list<int>::iterator particle_it =  particles.begin(); particle_it!=particles.end() ; particle_it ++){
                                State.GetDisplacement(*it,*particle_it,Displacement);

                                if(*it != *particle_it && Displacement.squaredNorm() < Cutoff*Cutoff){

                                        TNeighbor.push_back(CSimpleNeighbor<Dim>(*it,*particle_it,Displacement));\
                                }
                        }
                }
		
		Neighbors.push_back(TNeighbor);
        }
}

template <int Dim>
void CSimpleGrid<Dim>::Repopulate(CStaticState<Dim> &State)
{
	int Prod = 1;
	for(int d = 0 ; d < Dim ; d++)
		Prod *= Size(d);

	for(int i = 0 ; i < Prod ; i++)
		GridParticles[i].clear();

	Eigen::Matrix<double,Dim,1> Position;
        for(int i = 0 ; i < State.GetParticleNumber() ; i++)
        {
                State.GetParticlePosition(Position, i);
                GridParticles[IndexFromCoordinate(Position)].push_back(i);
                
        } 
}

template <int Dim>
const list<int> &CSimpleGrid<Dim>::GetCellParticles(int Index)
{
	return GridParticles[Index];
}


template <int Dim>
const list<int> &CSimpleGrid<Dim>::GetNeighborCells(Eigen::Matrix<double,Dim,1> &Position)
{
//	cout << "getting cell at index " << IndexFromCoordinate(Position) << endl; 
	return CellNeighbors[IndexFromCoordinate(Position)];
}

template <int Dim>
const list<CSimpleNeighbor<Dim> > &CSimpleGrid<Dim>::GetNeighbors(int particle)
{
	if(NeighborList)
		return NeighborList[particle];

	return list<CSimpleNeighbor<Dim> >();
}
}
#endif
