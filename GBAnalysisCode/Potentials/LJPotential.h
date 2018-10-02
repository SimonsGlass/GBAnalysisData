#include "BasePotential.h"

#ifndef LJ_POTENTIAL
#define LJ_POTENTIAL

/////////////////////////////////////////////////////////////////////////////////
//Lennard-Jones Potential class. 
//
//Description
//		This class describes particles that interact via a particular form
//		of a Lennard-Jones interaction, standard in the literature. (Kob-Andersen)
//		It takes a single parameter, the interaction strength \epsilon. The string
//		denoting this potential is "LJPotential"
//		
//		It has been specialized to the standard glass-forming set-up (sigma_AB = 0.8 sigma_AA
//		Sigma_BB = 0.88 sigma_AA, epsilon_AB = 1.5 epsilon, epsilon_BB = 0.5 epsilon
//		Cut off at r_ij=2.5 sigma_ij, and shifted to satisfy v(2.5 sigma_ij) = V'(2.5 sigma_ij) V''(2.5 sigma_ij) = 0.
//
//		NOTE: Current implementation specifically requires the ratio of particle sizes to be set like state.SetRadiiBi(f,1.136364) where f is
//			the fraction of small particles...i.e. the precision on the ratio of the sizes should be 1.0/1.136364
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

class CLJPotential : public CPotential
{	
private:
	dbl epsilon;	//!<the interaction strength

public:
//! @name Constructors and copy operators
///@{
	CLJPotential();
	CLJPotential(dbl _e);
	CLJPotential(const CLJPotential &box);
	const CLJPotential &operator=(const CLJPotential &box);
///@}

//! @name Functions to read and write potential configurations
///@{
	static string GetName();
	string DataToString() const;
 	void StringToData(string Data);
 	CLJPotential *Clone() const;
	
///@}
	
//! @name Functions to compute various derivatives
///@{ 
	dbl  Compute(dbl dr,dbl r) const;											//!<Compute the 0th derivative
	dbl  ComputeFirstDerivative(dbl dr,dbl r) const;							//!<Compute the 1st derivative
	dbl  ComputeSecondDerivative(dbl dr, dbl r) const;							//!<Compute the 2nd derivative
	dbl  ComputeThirdDerivative(dbl dr, dbl r) const;							//!<Compute the 3th derivative
	void ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const;				//!<Compute the 0th and 1st derivatives
	void ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const;	//!<Compute the 0th and 1st and 2nd derivatives
///@}

	void ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const;

//! @name Misc.
///@{
	//!Get the maximum distance between two particles of unit diameter such that they interact
	dbl ComputeSupport() const;
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
CLJPotential::CLJPotential()
{
	epsilon = 1.0;
}

CLJPotential::CLJPotential(dbl _e)
{
	epsilon = _e;
}

CLJPotential::CLJPotential(const CLJPotential &pot) : epsilon(pot.epsilon), CPotential(pot) //CPG comment: is is necessary to call the CPotential constructor? If so, then this should be done in CHarmonicPotential as well.
{
}
	
const CLJPotential &CLJPotential::operator=(const CLJPotential &pot)
{
	epsilon = pot.epsilon;
	CPotential::operator=(pot);
	return *this;
}

string CLJPotential::GetName()
{
	string s = "LJPotential";
	return s;
}

//functions to write potential configurations
string CLJPotential::DataToString() const
{
	stringstream ss; 
        ss << GetName() << STRING_DELIMITER << ConvertDblToHexString(epsilon);                                                  
        return FillString(ss);
}
	
//functions to read potential configurations
void CLJPotential::StringToData(string Data)
{
	vector<string> split = SplitString(Data,":");
	epsilon = atof(split[1].c_str());
	//CPG comment: I think this might be the problem. Try switching to the following:
	epsilon = ConvertHexStringToDbl(split[1]);
}

//Clone the potential opbject and return a pointer to the new copy.
CLJPotential* CLJPotential::Clone() const
{
        return new CLJPotential( *this );
} 
	
//functions to compute various derivatives
dbl CLJPotential::Compute(dbl dr,dbl r) const
{
	//CPG comment: be careful with using == in floating point comparison. 
	dbl epsij=epsilon;
	if (r < 0.84) epsij = 1.5*epsilon;
	if (r > 0.84 && r<0.9) epsij = 0.5* epsilon;
	dbl rat = r/dr;
	dbl rat3 = rat*rat*rat;
	dbl rat6 = rat3*rat3;
	dbl rat12 = rat6*rat6;

	return 0.01388888888888*epsij*(rat12-rat6)+epsij*(0.00157168-0.00107661/rat+0.000188239/(rat*rat));
}

dbl CLJPotential::ComputeFirstDerivative(dbl dr,dbl r) const
{
	dbl epsij=epsilon;
        if (r < 0.84) epsij = 1.5*epsilon;
        if (r > 0.84 && r<0.9) epsij = 0.5* epsilon;
	dbl rat = r/dr;
	dbl rat3 = rat*rat*rat;
        dbl rat6 = rat3*rat3;
        dbl rat12 = rat6*rat6;

	return epsij*(0.000376477*dr-0.00107661*r)/(r*r)+0.01388888888888*epsij*(6.0*rat6/dr-12.0*rat12/dr);
}

dbl CLJPotential::ComputeSecondDerivative(dbl dr,dbl r) const
{	
	dbl epsij=epsilon;
        if (r < 0.84) epsij = 1.5*epsilon;
        if (r > 0.84 && r<0.9) epsij = 0.5* epsilon;
	dbl rat = r/dr;
	dbl rat2 = rat*rat;
	dbl rat4 = rat2*rat2;
	dbl rat14 = rat4*rat4*rat4*rat2;
	
	return epsij*(0.000376477-0.583333*rat4*rat4+2.16667*rat14)/(r*r);
}

dbl CLJPotential::ComputeThirdDerivative(dbl dr,dbl r) const
{
	dbl epsij=epsilon;
        if (r < 0.84) epsij = 1.5*epsilon;
        if (r > 0.84 && r<0.9) epsij = 0.5* epsilon;
	dbl rat = r/dr;
	dbl rat3 = rat*rat*rat;
        dbl rat6 = rat3*rat3;
        dbl rat12 = rat6*rat6;
	dbl dri = 1.0/dr;
	dbl dri3 = dri*dri*dri;
        return 0.01388888888888*epsij*(336.0 * rat6*dri3-2184.0*rat12*dri3);
}

void CLJPotential::ComputeDerivatives01(dbl dr, dbl r, dbl &E, dbl &g) const
{
	dbl epsij=epsilon;
        if (r < 0.84) epsij = 1.5*epsilon;
        if (r > 0.84 && r<0.9) epsij = 0.5* epsilon;
	dbl rat = r/dr;
	dbl rat3 = rat*rat*rat;
        dbl rat6 = rat3*rat3;
        dbl rat12 = rat6*rat6;

	E = 0.01388888888888*epsij*(rat12-rat6)+epsij*(0.00157168-0.00107661/rat+0.000188239/(rat*rat));
	g = epsij*(0.000376477*dr-0.00107661*r)/(r*r)+0.01388888888888*epsij*(6.0*rat6/dr-12.0*rat12/dr);
}

void CLJPotential::ComputeDerivatives012(dbl dr, dbl r, dbl &E, dbl &g, dbl &k) const
{
	dbl epsij=epsilon;
        if (r < 0.84) epsij = 1.5*epsilon;
        if (r > 0.84 && r<0.9) epsij = 0.5* epsilon;
	dbl rat = r/dr;
        dbl rat3 = rat*rat*rat;
        dbl rat6 = rat3*rat3;
        dbl rat12 = rat6*rat6;
	dbl rat2 = rat*rat;
        dbl rat4 = rat2*rat2;
        dbl rat14 = rat4*rat4*rat4*rat2;

        k= epsij*(0.000376477-0.583333*rat4*rat4+2.16667*rat14)/(r*r);	

        E = 0.01388888888888*epsij*(rat12-rat6)+epsij*0.000395193-0.00013541485*epsij/rat;
        g = -0.00013541485*epsij/r+0.01388888888888*epsij*(6.0*rat6/dr-12.0*rat12/dr);
	if(g!=g)
		cout <<r<<"  "<<dr<<"  "<<E <<endl;
}

void CLJPotential::ComputeDerivatives0123(dbl dr, dbl r, dbl &E, dbl &g, dbl &k, dbl &t) const
{
	dbl epsij=epsilon;
        if (r < 0.84) epsij = 1.5*epsilon;
        if (r > 0.84 && r<0.9) epsij = 0.5* epsilon;
	dbl rat = r/dr;
        dbl rat3 = rat*rat*rat;
        dbl rat6 = rat3*rat3;
        dbl rat12 = rat6*rat6;
        dbl rat2 = rat*rat;
        dbl rat4 = rat2*rat2;
        dbl rat14 = rat4*rat4*rat4*rat2;

        k= epsij*(0.000376477-0.583333*rat4*rat4+2.16667*rat14)/(r*r);

        E = 0.01388888888888*epsij*(rat12-rat6)+epsij*0.000395193-0.00013541485*epsij/rat;
        g = -0.00013541485*epsij/r+0.01388888888888*epsij*(6.0*rat6/dr-12.0*rat12/dr);
	
	dbl dri = 1.0/dr;
        dbl dri3 = dri*dri*dri;
        t= 0.01388888888888*epsij*(336.0 * rat6*dri3-2184.0*rat12*dri3);
}


dbl CLJPotential::ComputeSupport() const
{
	return 2.5;
}

bool CLJPotential::Overlapping(dbl rad1, dbl rad2, dbl rlen2, dbl &sigma) const
{	
	dbl range = ComputeSupport();
	sigma=1.0;
	//if(rad1 == rad2 && fabs(rad1-.45929)<.001) sigma =1.0;
	//if(rad1 == rad2 && fabs(rad1-.521921)<.001) sigma = 1.13636;
	if(rad1 == rad2 && rad1 < 0.48) sigma =0.88;
        if(rad1 == rad2 && rad1 > 0.48) sigma =1.0;
	if(rad1 != rad2) sigma = 0.8; 
        return (rlen2<sigma*sigma*range*range);
}

}

#endif
