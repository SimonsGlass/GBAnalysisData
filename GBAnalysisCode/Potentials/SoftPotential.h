#include "BasePotential.h"

#ifndef SOFT_POTENTIAL
#define SOFT_POTENTIAL

/////////////////////////////////////////////////////////////////////////////////
//Harmonic Potential class. 
//
//Description
//		This class describes particles that interact via a central force repulsion.
//		It takes a single parameter, the interaction strength \epsilon. The string
//		denoting this potential is "SoftPotential"
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

class CSoftPotential : public CPotential
{	
private:
	dbl epsilon;	//!<the interaction strength
	dbl alpha;		//!<the power law of the interaction

public:
//! @name Constructors and copy operators
///@{
	CSoftPotential();
	CSoftPotential(dbl _e, dbl _alpha=2.);
	CSoftPotential(const CSoftPotential &pot);
	const CSoftPotential &operator=(const CSoftPotential &pot);
///@}

//! @name Functions to read and write potential configurations
///@{
	static string GetName();
	virtual string DataToString() const;
 	virtual void StringToData(string Data);
 	virtual CSoftPotential *Clone() const;
	
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
CSoftPotential::CSoftPotential()
	: epsilon(1.0), alpha(2.0)
{}

CSoftPotential::CSoftPotential(dbl _e, dbl _alpha)
	: epsilon(_e), alpha(_alpha)
{}

CSoftPotential::CSoftPotential(const CSoftPotential &pot) 
	: epsilon(pot.epsilon), alpha(pot.alpha)
{}
	
const CSoftPotential& CSoftPotential::operator=(const CSoftPotential &pot)
{
	if(this != &pot)
	{
		epsilon = pot.epsilon;
		alpha   = pot.alpha;
	}
	return *this;
}

string CSoftPotential::GetName()
{
	string s = "SoftPotential";
	return s;
}

//functions to write potential configurations
string CSoftPotential::DataToString() const
{
	stringstream ss;
	ss << GetName() << STRING_DELIMITER << ConvertDblToHexString(epsilon) << STRING_DELIMITER << ConvertDblToHexString(alpha);
	return FillString(ss);
//	ss << GetName() << ":" << ConvertDblToHexString(epsilon) << ":" << ConvertDblToHexString(alpha);
//	return ss.str();
}
	
//functions to read potential configurations
void CSoftPotential::StringToData(string Data)
{
	vector<string> split = SplitString(Data,":");
	epsilon = ConvertHexStringToDbl(split[1]);
	alpha   = ConvertHexStringToDbl(split[2]);
}

//Clone the potential opbject and return a pointer to the new copy.
CSoftPotential* CSoftPotential::Clone() const
{
	return new CSoftPotential( *this );
}
	
//functions to compute various derivatives
dbl CSoftPotential::Compute(dbl dr,dbl r) const
{
	dbl delta = 1-dr/r;
	return epsilon*std::pow(delta,alpha)/alpha;
}

dbl CSoftPotential::ComputeFirstDerivative(dbl dr,dbl r) const
{
	dbl delta = 1-dr/r;
	return -epsilon*std::pow(delta,alpha-1)/r;
}

dbl CSoftPotential::ComputeSecondDerivative(dbl dr,dbl r) const
{
	dbl delta = 1-dr/r;
	return epsilon*(alpha-1)*std::pow(delta,alpha-2)/POW2(r);
}

dbl CSoftPotential::ComputeThirdDerivative(dbl dr,dbl r) const
{
	dbl delta = 1-dr/r;
	return -epsilon*(alpha-1)*(alpha-2)*std::pow(delta,alpha-3)/POW3(r);
}

void CSoftPotential::ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const
{
	dbl delta = 1-dr/r;
	E = epsilon*std::pow(delta,alpha)/alpha;
	g =-epsilon*std::pow(delta,alpha-1)/r;
}

void CSoftPotential::ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const
{
	dbl delta = 1-dr/r;
	E = epsilon*std::pow(delta,alpha)/alpha;
	g =-epsilon*std::pow(delta,alpha-1)/r;
	k = epsilon*(alpha-1)*std::pow(delta,alpha-2)/POW2(r);
}

void CSoftPotential::ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const
{
	dbl delta = 1-dr/r;
	E = epsilon*std::pow(delta,alpha)/alpha;
	g =-epsilon*std::pow(delta,alpha-1)/r;
	k = epsilon*(alpha-1)*std::pow(delta,alpha-2)/POW2(r);
	t =-epsilon*(alpha-1)*(alpha-2)*std::pow(delta,alpha-3)/POW3(r);
}

dbl CSoftPotential::ComputeSupport() const
{
	return 1.;
}

bool CSoftPotential::Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const
{
	sigma = rad1 + rad2;
	return (rlen2<sigma*sigma);
}


}

#endif
