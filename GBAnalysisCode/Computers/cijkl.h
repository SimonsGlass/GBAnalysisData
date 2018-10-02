#ifndef CIJKL_H
#define CIJKL_H

#include "../Resources/std_include.h"
#include "../Resources/Resources.h"

template<int Dim> class cCIJKL;

//! Class to store elastic constants for a 2 dimensional system.
template<>
class cCIJKL<2>
{
private:
	static const int Dim = 2;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
	enum{xxxx=0, yyyy, xyxy, xxyy, xxxy, yyxy};		//!<Flags to refer to the constants.

public:
	static const int num_constants = 6;				//!<Number of independent elastic constants.

	dbl a[num_constants];	//!<Array of the elastic constants.

//! @name references to refer to the elastic constants individually.
///@{
	dbl &cxxxx;
	dbl &cyyyy;
	dbl &cxyxy;
	dbl &cxxyy;
	dbl &cxxxy;
	dbl &cyyxy;
///@}

//! @name Constructors and operators
///@{
	cCIJKL();								//!<Default constructor.
	cCIJKL(const cCIJKL &src);				//!<Copy constructor.
	cCIJKL& operator=(const cCIJKL &src);	//!<Equals operator.

///@}

//! @name Methods to help calculate the elastic constants.
///@{
	void set_constant(dbl twice_dE_over_V, int ii);					//!<Set a constant from energy of a deformation.
	void set_strain_tensor(dmat &strain_tensor, int ii) const;		//!<Get the strain tensor for the ii'th iteration.

///@}

//! @name Methods to decompose the elastic constants.
///@{
	//! Calculate some temporary values.
	void calc_temp_combos(dbl &G0, dbl &Gpi4, dbl &A2, dbl &p2, dbl &A4, dbl &p4) const;
	//! Calculate all the angle averaged elastic constants.
	void CalculateRotationalAverages(dbl &B, dbl &GDC, dbl &GAC, dbl &UDC, dbl &UAC, dbl &DDC, dbl &DAC, dbl &TDC, dbl &TAC) const;
	dbl CalculateBulkModulus() const;			//!<Get the bulk modulus
	dbl CalculateShearModulus() const;			//!<Get the shear modulus
	dbl calc_weakest_deformation() const;	//!<Get the weakest response to any strain.

///@}

//! @name Misc.
///@{
	void SetZero();								//!<Set all elastic constants to zero.
	void SetArtificialValues(dbl c);			//!<Set all elastic constants to c.
	Eigen::Matrix<dbl,3,3> get_Cik() const;		//!<Get the elastic constants as a 3 by 3 matrix.
	void print() const;							//!<Print the 3 by 3 matrix to stdout.

///@}
};


//! Class to store elastic constants for a 3 dimensional system.
template<>
class cCIJKL<3>
{
private:
	static const int Dim = 3;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
	enum{
		//Uniaxial compression
		xxxx=0, // 0
		yyyy, // 1
		zzzz, // 2
		//Shear
		yzyz, // 3
		xzxz, // 4
		xyxy, // 5
		//cross
		yyzz, // 6
		xxzz, // 7
		xxyy, // 8
		//Dilatency
		xxyz, // 9
		xxxz, //10
		xxxy, //11
		yyyz, //12
		yyxz, //13
		yyxy, //14
		zzyz, //15
		zzxz, //16
		zzxy, //17
		//Last
		yzxz, //18
		yzxy, //19
		xzxy  //20
	};
public:
	static const int num_constants = 21; //!<Number of independent elastic constants.
	
	dbl a[num_constants];	//!<Array of the elastic constants

//! @name References to refer to the elastic constants individually.
///@{
	dbl &cxxxx; // 0
	dbl &cyyyy; // 1
	dbl &czzzz; // 2
	dbl &cyzyz; // 3
	dbl &cxzxz; // 4
	dbl &cxyxy; // 5
	dbl &cyyzz; // 6
	dbl &cxxzz; // 7
	dbl &cxxyy; // 8
	dbl &cxxyz; // 9
	dbl &cxxxz; //10
	dbl &cxxxy; //11
	dbl &cyyyz; //12
	dbl &cyyxz; //13
	dbl &cyyxy; //14
	dbl &czzyz; //15
	dbl &czzxz; //16
	dbl &czzxy; //17
	dbl &cyzxz; //18
	dbl &cyzxy; //19
	dbl &cxzxy; //20
///@}

//! @name Constructors and operators.
///@{
	cCIJKL();								//!<Default constructor.
	cCIJKL(const cCIJKL &src);				//!<Copy constructor.
	cCIJKL& operator=(const cCIJKL &src);	//!<Equals operator.

///@}

//! @name Methods to help calculate the elastic constants.
///@{
	void set_constant(dbl twice_dE_over_V, int ii);			//!<Set a constant from the energy of a deformation.	
	void set_strain_tensor(dmat &strain_tensor, int ii);	//!<Get the strain tensor for the ii'th iteration.

///@}

//! @name Methods to decompose the elastic constants.
///@{
	dbl CalculateBulkModulus() const;	//!<Get the bulk modulus.
	dbl CalculateShearModulus() const;	//!<Get the shear modulus.
	//! Calculate all the angle averaged elastic constants.
	void CalculateRotationalAverages(dbl &B, dbl &GDC, dbl &GAC, dbl &UDC, dbl &UAC, dbl &DDC, dbl &DAC, dbl &TDC, dbl &TAC) const;	
	dbl calc_weakest_deformation() const; //!<Get the weakest response to any strain.
	
///@}

//! @name Misc.
///@{
	void SetZero();								//!<Set all elastic constants to zero.
	void SetArtificialValues(dbl c);			//!<Set all elastic constants to c.
	Eigen::Matrix<dbl,6,6> get_Cik() const;		//!<Get the elastic constants as a 6 by 6 matrix.
	void print() const;							//!<Print the 6 by 6 matrix to stdout.

///@}
};

#endif //cCIJKL_H

