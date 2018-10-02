#ifndef BOX

#define BOX

/////////////////////////////////////////////////////////////////////////////////
//Box class. 
//
//Description
//		This is a virtual class that describes box shape. The shape of the box is
//		specified by a dxd matrix that describes the mapping from the unit box [0,1]x...x[0,1]
//		to a corresponding parallelpiped. Overloaded classes will specify the box topology
//		by specifying how the various edges of the box connect to one another. This class
//		implements functions to take points from the base space to the image space. Virtual
//		functions are provided to move points and compute minimum displacements between 
//		points. 
//
//		A framework to generically read and write box configurations to and from netcdf
//		files is also implemented. New box objects that inherit from the box class need
//		to register a string with their name to the map of box types. They must also specify
//		how to read/write any parameters (which must all be doubles) to a string. Inherited
//		objects must further implement several copy constructors. Finally, inherited classes must 
//		implement a function to create a new object of the given class.
//
//Global Variables
//		A map from string to box type as an STL map.
//
//Variables
// 		Transformation as a dxd matrix.
//		Inverse transformation as a dxd matrix.
//		
//Implements
//		Reading/Writing to netcdf files.
//		Mapping to and from the unit box to a transformed parallelopiped.
//
//Virtual Functions
//		Calculating minimal distances between points.
//		Moving particles while respecting the boundary conditions.
//		Mapping the system parameters to and from a string.
//
//File Formats
//		NetCDF
//			## Note: only files of a single Dimension may be stored in a given 
//			## file. To denote whether the NetCDF file has been populated with variables,
//			## dimensions, etc... an attribute "Box_Populated" gets created.
//			-> Two dimensions: records (unlimited) and dimension.
//			-> One attribute Box_Populated.
//			-> Transformation as a variable
//			-> String of parameters as a string.
//
/////////////////////////////////////////////////////////////////////////////////

#include "../Resources/std_include.h"
#include <Eigen/LU>
//#include "netcdfcpp.h"

namespace LiuJamming
{

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//!Abstract base class for objects that handle boundary conditions.
template <int Dim>
class CBox
{
private:
	typedef Eigen::Matrix<dbl,Dim,1> dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
	typedef Eigen::VectorBlock<Eigen::VectorXd,Dim> dvecBlock;

	static std::map<string,CBox<Dim>*> BoxTypes; //!<Global map to reference box objects by name

public:
	enum{RECTANGULAR, SYMMETRIC_QUADRILATERAL};
protected:
	int BoxSymmetry; //!<Flag that potentially increases efficiency

private:
//variables specifying the transformation
	dmat Transformation;			//!<Transformation matrix
	dmat Inverse_Transformation;	//!<Inverse transformation matrix
	
public:
//constructors and copy operators
	CBox(int symmetry=SYMMETRIC_QUADRILATERAL);						//!<Default constructor
	CBox(const dmat &Trans, int symmetry=SYMMETRIC_QUADRILATERAL);	//!<Construct from transformation matrix
	CBox(const CBox &box);											//!<Copy constructor
	const CBox &operator=(const CBox &box);							//!<Copy operator


	static CBox* SetFromStringAndMatrix(string str, dmat &trans);
	
//functions to read and write box configurations
	//static string GetName() const = 0;
	virtual string DataToString() const = 0;
 	virtual void StringToData(string Data) = 0;
 	virtual CBox* Clone() const = 0;

	void SetSymmetry(int symmetry);	//!<Set the assumed box symmetry
	int  GetSymmetry() const;		//!<Get the assumed box symmetry

//functions using the transformation matrix
	void SetTransformation(dmat &Trans);
	void GetTransformation(dmat &Trans) const;
	void GetInverseTransformation(dmat &ITrans) const;
	void Transform(dvec &Point) const;
	void Transform(Eigen::VectorXd &Points) const;
	void InverseTransform(dvec &Point) const;
	void InverseTransform(Eigen::VectorXd &Points) const;
	void InverseTransformAndMove(Eigen::VectorXd &Points, const Eigen::VectorXd &t_Displacement) const;

//set and get the volume
	virtual void SetVolume(dbl V);			//!<Set the volume of the box
	virtual dbl CalculateVolume() const;	//!<Calculate the volume of the box

	virtual void GetPeriodicDimensions(std::vector<int> &) const = 0; //!<get a list of the periodic dimensions
	
//functions involving the boundary
	virtual void MoveParticles(Eigen::VectorXd &Points, const Eigen::VectorXd &Displacements) const = 0;
	virtual void MoveParticle(dvecBlock Point, dvec const &Displacement) const = 0; //dvecBlock's are themselves references, so Point should NOT be passed as a reference.
	virtual void MinimumDisplacement(const dvec &PointA, const dvec &PointB, dvec &Displacement) const = 0;
	virtual void MinimumDisplacement(const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointA,const Eigen::VectorBlock<Eigen::VectorXd,Dim> &PointB, dvec &Displacement) const = 0;
	virtual void MinimumDisplacement(const Eigen::VectorXd &PointA, const Eigen::VectorXd &PointB, Eigen::VectorXd &Displacement) const = 0;

	void GetMaxTransformedDistance(dvec &dist) const;

};

template <int Dim>
std::map<string,CBox<Dim>*> CBox<Dim>::BoxTypes;


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
CBox<Dim>* CBox<Dim>::SetFromStringAndMatrix(string str, dmat &trans)
{
	vector<string> split = SplitString(str,":");
	CBox<Dim> *box = BoxTypes[split[0]]->Clone();
	box->StringToData(str);
	box->SetTransformation(trans);
	return box;
}



/*
//global functions to read box configurations
template <int Dim>
CBox<Dim> *CBox<Dim>::Read(const NcFile &file,int record)
{
	if(!CheckNetCDF(file))
		throw(CException("CBox<Dim>::ReadBox","Attempting to read a box from a file that has no appropriate box data."));
	
	if(Dim!=file.get_dim("System_Dimension")->size())
		throw(CException("CBox<Dim>::ReadBox","Attempting to read a box from a record has an inconsistent nubmer of dimensions."));
	
	if(record>=file.get_dim("Records")->size())
		throw(CException("CBox<Dim>::ReadBox","Attempting to read a box from a record that does not exist."));

	//read the matrix
	dmat Transformation;
	NcVar *TransVar = file.get_var("Box_Transformation");
	TransVar->set_cur(record);
	TransVar->get(Transformation.data(),1,Dim,Dim);
	
	//read the box type
	char data[256];
	NcVar *DataVar = file.get_var("Box_Data");
	DataVar->set_cur(record);
	DataVar->get(data,1,256);
	string s_data = data;
	//split up the string using colons.
	vector<string> split = SplitString(s_data,":");
	//the first element of the split string should be the box name so create that kind of box
	CBox *box = BoxTypes[split[0]]->Create();
	//go set the data using the rest of the string and set the transformation
	box->StringToData(s_data);
	box->SetTransformation(Transformation);
	
	return box;
}
*/

//constructors and copy operators
template <int Dim>
CBox<Dim>::CBox(int symmetry)
{
	Transformation = dmat::Identity();
	Inverse_Transformation = dmat::Identity();
	BoxSymmetry = symmetry;
}

template <int Dim>
CBox<Dim>::CBox(const dmat &Trans, int symmetry)
{
	Transformation = Trans;
	Inverse_Transformation = Transformation.inverse();
	BoxSymmetry = symmetry;
}
	
template <int Dim>
CBox<Dim>::CBox(const CBox &box)
{
	*this = box;
}

template <int Dim> 
const CBox<Dim>& CBox<Dim>::operator=(const CBox<Dim> &box)
{
	if(this != &box)
	{
		box.GetTransformation(Transformation);
		Inverse_Transformation = Transformation.inverse();
		BoxSymmetry = box.BoxSymmetry;
	}
	return *this;
}

template <int Dim>
void CBox<Dim>::SetSymmetry(int symmetry)
{
	BoxSymmetry = symmetry;
}

template <int Dim>
int  CBox<Dim>::GetSymmetry() const
{
	return BoxSymmetry;
}
	
//functions using the transformation matrix
template <int Dim>
void CBox<Dim>::SetTransformation(dmat &Trans)
{
	Transformation = Trans;
	Inverse_Transformation = Transformation.inverse();
}

template <int Dim>
void CBox<Dim>::GetTransformation(dmat &Trans) const
{
	Trans = Transformation;
}

template <int Dim>
void CBox<Dim>::GetInverseTransformation(dmat &ITrans) const
{
	ITrans = Inverse_Transformation;
}

template <int Dim>
inline void CBox<Dim>::Transform(dvec &Point) const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			Point = Transformation.diagonal().cwiseProduct(Point); 
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			Point = Transformation * Point;	
	}
}

template <int Dim>
inline void CBox<Dim>::Transform(Eigen::VectorXd &Points) const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			for(int i = 0 ; i<Points.size()/Dim ; i++)
				Points.segment<Dim>(Dim*i) = Transformation.diagonal().cwiseProduct(Points.segment<Dim>(Dim*i));
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			for(int i = 0 ; i<Points.size()/Dim ; i++)
				Points.segment<Dim>(Dim*i) = Transformation*Points.segment<Dim>(Dim*i);
	}
}

template <int Dim>
inline void CBox<Dim>::InverseTransform(dvec &Point) const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			Point = Inverse_Transformation.diagonal().cwiseProduct(Point); 
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			Point = Inverse_Transformation * Point;
	}
}

template <int Dim>
inline void CBox<Dim>::InverseTransform(Eigen::VectorXd &Points) const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			for(int i = 0 ; i<Points.size()/Dim ; i++)
				Points.segment<Dim>(Dim*i) = Inverse_Transformation.diagonal().cwiseProduct(Points.segment<Dim>(Dim*i));
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			for(int i = 0 ; i<Points.size()/Dim ; i++)
				Points.segment<Dim>(Dim*i) = Inverse_Transformation*Points.segment<Dim>(Dim*i);
	}
}

template <int Dim>
void CBox<Dim>::InverseTransformAndMove(Eigen::VectorXd &Points, const Eigen::VectorXd &t_Displacement) const
{
	dvec Displacement;
	int np = Points.size()/Dim;
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			for(int i=0; i<np; ++i)
			{
				Displacement = Inverse_Transformation.diagonal().cwiseProduct(t_Displacement.segment<Dim>(Dim*i));
				MoveParticle(Points.segment<Dim>(Dim*i), Displacement);
			}
			break;
		case SYMMETRIC_QUADRILATERAL:
		default:
			for(int i=0; i<np; ++i)
			{
				Displacement = Inverse_Transformation * t_Displacement.segment<Dim>(Dim*i);
				MoveParticle(Points.segment<Dim>(Dim*i), Displacement);
			}
	}
};



template <int Dim>
void CBox<Dim>::SetVolume(dbl V)
{
	dbl Vold = CalculateVolume();
	dbl Lrescale = std::pow(V/Vold,1./((dbl)Dim));
	Transformation *= Lrescale;
	Inverse_Transformation = Transformation.inverse(); //Could just divide Inverse_Transformation by Lrescale, but this is probably more stable.
}


template <int Dim>
dbl CBox<Dim>::CalculateVolume() const
{
	switch(BoxSymmetry)
	{
		case RECTANGULAR: 
			return Transformation.diagonal().prod();
		case SYMMETRIC_QUADRILATERAL:
		default:
			return fabs(Transformation.determinant());
	}
}

/**
 *	In general, a Dim-dimensional sphere of unit radius (in real units)
 *	is represented by an ellipse in transformed units. 
 *	This method calculates the largest component, in transformed units, in 
 *	each direction over the enitre ellipse.
 */
template <int Dim>
void CBox<Dim>::GetMaxTransformedDistance(dvec &dist) const
{
	for(int dd=0; dd<Dim; ++dd)
		dist[dd] = Inverse_Transformation.row(dd).norm();
}

}

#endif
