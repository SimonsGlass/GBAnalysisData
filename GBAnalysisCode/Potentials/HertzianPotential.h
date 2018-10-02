#include "BasePotential.h"

#ifndef HERTZIAN_POTENTIAL
#define HERTZIAN_POTENTIAL

/////////////////////////////////////////////////////////////////////////////////
//Harmonic Potential class. 
//
//Description
//		This class describes particles that interact via a central force repulsion.
//		It takes a single parameter, the interaction strength \epsilon. The string
//		denoting this potential is "HertzianPotential"
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

class CHertzianPotential : public CPotential
{	
private:
	dbl epsilon;	//!<the interaction strength

public:
//! @name Constructors and copy operators
///@{
	CHertzianPotential();
	CHertzianPotential(dbl _e);
	CHertzianPotential(const CHertzianPotential &pot);
	const CHertzianPotential &operator=(const CHertzianPotential &pot);
///@}

//! @name Functions to read and write potential configurations
///@{
	static string GetName();
	virtual string DataToString() const;
 	virtual void StringToData(string Data);
 	virtual CHertzianPotential *Clone() const;
	
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
CHertzianPotential::CHertzianPotential()
	: epsilon(1.0)
{}

CHertzianPotential::CHertzianPotential(dbl _e)
	: epsilon(_e)
{}

CHertzianPotential::CHertzianPotential(const CHertzianPotential &pot) 
	: epsilon(pot.epsilon)
{}
	
const CHertzianPotential& CHertzianPotential::operator=(const CHertzianPotential &pot)
{
	if(this != &pot)
	{
		epsilon = pot.epsilon;
	}
	return *this;
}

string CHertzianPotential::GetName()
{
	string s = "HertzianPotential";
	return s;
}

//functions to write potential configurations
string CHertzianPotential::DataToString() const
{
	stringstream ss;
	ss << GetName() << STRING_DELIMITER << ConvertDblToHexString(epsilon);
	return FillString(ss);
//	ss << GetName() << ":" << ConvertDblToHexString(epsilon);
//	return ss.str();
}
	
//functions to read potential configurations
void CHertzianPotential::StringToData(string Data)
{
	vector<string> split = SplitString(Data,":");
	epsilon = ConvertHexStringToDbl(split[1]);
}

//Clone the potential opbject and return a pointer to the new copy.
CHertzianPotential* CHertzianPotential::Clone() const
{
	return new CHertzianPotential( *this );
}
	
//functions to compute various derivatives
dbl CHertzianPotential::Compute(dbl dr,dbl r) const
{
	dbl delta = 1-dr/r;
	dbl SQRTdelta = sqrt(delta);
	return 0.4*epsilon*POW2(delta)*SQRTdelta;
}

dbl CHertzianPotential::ComputeFirstDerivative(dbl dr,dbl r) const
{
	dbl delta = 1-dr/r;
	dbl SQRTdelta = sqrt(delta);
	return -epsilon*delta*SQRTdelta/r;
}

dbl CHertzianPotential::ComputeSecondDerivative(dbl dr,dbl r) const
{
	dbl delta = 1-dr/r;
	dbl SQRTdelta = sqrt(delta);
	return 1.5*epsilon*SQRTdelta/POW2(r);
}

dbl CHertzianPotential::ComputeThirdDerivative(dbl dr,dbl r) const
{
	dbl delta = 1-dr/r;
	dbl SQRTdelta = sqrt(delta);
	return -0.75*epsilon/(POW3(r)*sqrt(1-dr/r));
}

void CHertzianPotential::ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const
{
	dbl delta = 1-dr/r;
	dbl SQRTdelta = sqrt(delta);
	E = 0.4*epsilon*POW2(delta)*SQRTdelta;
	g = -epsilon*delta*SQRTdelta/r;
}

void CHertzianPotential::ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const
{
	dbl delta = 1-dr/r;
	dbl SQRTdelta = sqrt(delta);
	E = 0.4*epsilon*POW2(delta)*SQRTdelta;
	g = -epsilon*delta*SQRTdelta/r;
	k = 1.5*epsilon*SQRTdelta/POW2(r);
}

void CHertzianPotential::ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const
{
	dbl delta = 1-dr/r;
	dbl SQRTdelta = sqrt(delta);
	E = 0.4*epsilon*POW2(delta)*SQRTdelta;
	g = -epsilon*delta*SQRTdelta/r;
	k = 1.5*epsilon*SQRTdelta/POW2(r);
	t = -0.75*epsilon/(POW3(r)*sqrt(1-dr/r));
}

dbl CHertzianPotential::ComputeSupport() const
{
	return 1.;
}

bool CHertzianPotential::Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const
{
	sigma = rad1 + rad2;
	return (rlen2<sigma*sigma);
}


}

#endif
