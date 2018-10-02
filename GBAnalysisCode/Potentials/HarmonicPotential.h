#include "BasePotential.h"

#ifndef HARMONIC_POTENTIAL
#define HARMONIC_POTENTIAL

/////////////////////////////////////////////////////////////////////////////////
//Harmonic Potential class. 
//
//Description
//		This class describes particles that interact via a central force repulsion.
//		It takes a single parameter, the interaction strength \epsilon. The string
//		denoting this potential is "HarmonicPotential"
//
//Global Variables
//
//Variables
//		Interaction strength \epsilon as a dbl
//		
//Implements
//		Calculating various orders of derivatives of the potential
//		Calculating the support of the potential as a function of radius.
//		Convert from the parameter \epsilon to a string reprentation and vice versa
//
//
/////////////////////////////////////////////////////////////////////////////////

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

//! Class to implement harmonic interactions.

class CHarmonicPotential : public CPotential
{	
private:
	dbl epsilon;	//!<the interaction strength

public:
//! @name Constructors and copy operators
///@{
	CHarmonicPotential();
	CHarmonicPotential(dbl _e);
	CHarmonicPotential(const CHarmonicPotential &pot);
	const CHarmonicPotential &operator=(const CHarmonicPotential &pot);
///@}

//! @name Functions to read and write potential configurations
///@{
	static string GetName();
	virtual string DataToString() const;
 	virtual void StringToData(string Data);
 	virtual CHarmonicPotential *Clone() const;
	
///@}
	
//! @name Functions to compute various derivatives
///@{ 
	dbl  Compute(dbl dr,dbl r) const;													//!<Compute the 0th derivative
	dbl  ComputeFirstDerivative(dbl dr,dbl r) const;									//!<Compute the 1st derivative
	dbl  ComputeSecondDerivative(dbl dr, dbl r) const;									//!<Compute the 2nd derivative
	dbl  ComputeThirdDerivative(dbl dr, dbl r) const;									//!<Compute the 3th derivative
	void ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const;						//!<Compute the 0th and 1st derivatives
	void ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const;			//!<Compute the 0th and 1st and 2nd derivatives
	void ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const;	//!<Compute the 0th and 1st and 2nd and 3rd derivatives
	
///@}

//! @name Misc.
///@{
	//!Get the maximum distance between two particles of unit diameter such that they interact
	dbl ComputeSupport() const;

	//!Set sigma and determine if rad1, rad2 and rlen2 are such that there is overlap
	bool Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const;

///@}
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//constructors and copy operators
CHarmonicPotential::CHarmonicPotential()
	: epsilon(1.0)
{}

CHarmonicPotential::CHarmonicPotential(dbl _e)
	: epsilon(_e)
{}

CHarmonicPotential::CHarmonicPotential(const CHarmonicPotential &pot) 
	: epsilon(pot.epsilon)
{}
	
const CHarmonicPotential& CHarmonicPotential::operator=(const CHarmonicPotential &pot)
{
	if(this != &pot)
	{
		epsilon = pot.epsilon;
	}
	return *this;
}

string CHarmonicPotential::GetName()
{
	string s = "HarmonicPotential";
	return s;
}

//functions to write potential configurations
string CHarmonicPotential::DataToString() const
{
	stringstream ss;
	ss << GetName() << STRING_DELIMITER << ConvertDblToHexString(epsilon);
	return FillString(ss);
//	ss << GetName() << ":" << ConvertDblToHexString(epsilon);
//	return ss.str();
}
	
//functions to read potential configurations
void CHarmonicPotential::StringToData(string Data)
{
	vector<string> split = SplitString(Data,":");
	epsilon = ConvertHexStringToDbl(split[1]);
}

//Clone the potential opbject and return a pointer to the new copy.
CHarmonicPotential* CHarmonicPotential::Clone() const
{
	return new CHarmonicPotential( *this );
}
	
//functions to compute various derivatives
dbl CHarmonicPotential::Compute(dbl dr,dbl r) const
{
	return 0.5*epsilon*POW2(1-dr/r);
}

dbl CHarmonicPotential::ComputeFirstDerivative(dbl dr,dbl r) const
{
	return -epsilon/r*(1-dr/r);
}

dbl CHarmonicPotential::ComputeSecondDerivative(dbl dr,dbl r) const
{
	return epsilon/POW2(r);
}

dbl CHarmonicPotential::ComputeThirdDerivative(dbl dr,dbl r) const
{
	return 0.;
}

void CHarmonicPotential::ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const
{
	dbl delta = 1-dr/r;
	E = 0.5*epsilon*POW2(delta);
	g = -epsilon*delta/r;
}

void CHarmonicPotential::ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const
{
	dbl delta = 1-dr/r;
	E = 0.5*epsilon*POW2(delta);
	g = -epsilon*delta/r;
	k = epsilon/POW2(r);
}

void CHarmonicPotential::ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const
{
	dbl delta = 1-dr/r;
	E = 0.5*epsilon*POW2(delta);
	g = -epsilon*delta/r;
	k = epsilon/POW2(r);
	t = 0.;
}

dbl CHarmonicPotential::ComputeSupport() const
{
	return 1.;
}

bool CHarmonicPotential::Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const
{
	sigma = rad1 + rad2;
	return (rlen2<sigma*sigma);
}


}

#endif
