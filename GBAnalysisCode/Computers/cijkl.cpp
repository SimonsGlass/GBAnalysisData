#include "cijkl.h"
#include <Eigen/Dense>

cCIJKL<2>::cCIJKL()
	:	cxxxx (a[xxxx]),
		cyyyy (a[yyyy]),
		cxyxy (a[xyxy]),
		cxxyy (a[xxyy]),
		cxxxy (a[xxxy]),
		cyyxy (a[yyxy])
{
	for(int i=0; i<num_constants; i++)
		a[i] = 0.;
};

cCIJKL<2>::cCIJKL(const cCIJKL<2> &src)
	:	cxxxx (a[xxxx]),
		cyyyy (a[yyyy]),
		cxyxy (a[xyxy]),
		cxxyy (a[xxyy]),
		cxxxy (a[xxxy]),
		cyyxy (a[yyxy])
{
	(*this) = src;
};

cCIJKL<2>& cCIJKL<2>::operator=(const cCIJKL<2> &src)
{
	if(this != &src)
	{
		for(int i=0; i<num_constants; ++i)
			a[i] = src.a[i];
	}
	return *this;
};

void cCIJKL<2>::set_constant(dbl twice_dE_over_V, int ii)
{
	assert(ii>=0 && ii<num_constants);
	switch(ii)
	{
		case 0: //Apply ((1 0) (0 0))
			cxxxx = twice_dE_over_V;
			break;
		case 1: //Apply ((0 0) (0 1))
			cyyyy = twice_dE_over_V;
			break;
		case 2: //Apply ((0 1) (1 0))
			cxyxy = 0.25*twice_dE_over_V;
			break;
		case 3: //Apply ((1 0) (0 1))
			cxxyy = 0.5*(twice_dE_over_V - cxxxx - cyyyy);
			break;
		case 4: //Apply ((1 1) (1 0))
			cxxxy = 0.25*(twice_dE_over_V - 4.*cxyxy - cxxxx);
			break;
		case 5: //Apply ((1 1) (1 0))
			cyyxy = 0.25*(twice_dE_over_V - 4.*cxyxy - cyyyy);
			break;
	}
};

void cCIJKL<2>::set_strain_tensor(dmat &strain_tensor, int ii) const
{
	strain_tensor.setZero();
	const dbl strain = 1.;
	switch(ii)
	{
		case 0:	//Apply ((1 0) (0 0))
			strain_tensor(0,0) = strain;
			break;
		case 1: //Apply ((0 0) (0 1))
			strain_tensor(1,1) = strain;
			break;
		case 2: //Apply ((0 1) (1 0))
			strain_tensor(0,1) = strain_tensor(1,0) = strain;
			break;
		case 3: //Apply ((1 0) (0 1))
			strain_tensor(0,0) = strain_tensor(1,1) = strain;
			break;
		case 4: //Apply ((1 1) (1 0))
			strain_tensor(0,0) = strain_tensor(0,1) = strain_tensor(1,0) = strain;
			break;
		case 5: //Apply ((0 1) (1 1))
			strain_tensor(1,1) = strain_tensor(0,1) = strain_tensor(1,0) = strain;
			break;
	}
};


void cCIJKL<2>::calc_temp_combos(dbl &G0, dbl &Gpi4, dbl &A2, dbl &p2, dbl &A4, dbl &p4) const
{
	G0 = cxyxy;
	Gpi4 = 0.25*(cxxxx+cyyyy-2.0*cxxyy);

	A2 = sqrt(0.25*POW2(cxxxx-cyyyy) + POW2(cxxxy+cyyxy));
	p2 = std::atan2(cxxxx-cyyyy, -2.0*(cxxxy+cyyxy));
	A4 = -0.5*sqrt(POW2(cxxxy-cyyxy) + POW2(G0-Gpi4));
	p4 = std::atan2(G0-Gpi4, cxxxy-cyyxy);
}

void cCIJKL<2>::CalculateRotationalAverages(dbl &B, dbl &GDC, dbl &GAC, dbl &UDC, dbl &UAC, dbl &DDC, dbl &DAC, dbl &TDC, dbl &TAC) const
{
	dbl G0, Gpi4, A2, p2, A4, p4;
	calc_temp_combos(G0, Gpi4, A2, p2, A4, p4);

	B = CalculateBulkModulus();
	GDC = CalculateShearModulus();
	GAC = sqrt( (POW2(G0-Gpi4) + POW2(cxxxy-cyyxy))/8.0 );
	UDC = B + GDC;
	UAC = sqrt( 0.5*(POW2(A2) + POW2(A4)) );
	DDC = 0.0;
	DAC = sqrt( (POW2(A2) + 4.0*POW2(A4))/8.0 );
	TDC = 0.0;
	TAC = 0.0;
}

dbl cCIJKL<2>::CalculateBulkModulus() const
{
	return 0.25*(cxxxx+cyyyy+2.*cxxyy);
}

dbl cCIJKL<2>::CalculateShearModulus() const
{
	return 0.125*(4.*cxyxy + cxxxx+cyyyy-2.*cxxyy); //0.125 = 1/8
}

dbl cCIJKL<2>::calc_weakest_deformation() const
{
	//Eigen::Matrix<dbl, 3, 3> Cik = get_Cik();
	//Eigen::SelfAdjointEigenSolver< Eigen::Matrix<dbl,3,3> > eigensolver((Cik.selfadjointView<Eigen::Upper>()));
	//return (eigensolver.eigenvalues()).minCoeff();
	
	typedef Eigen::Matrix<dbl, Eigen::Dynamic, Eigen::Dynamic> MAT_TYPE;
	MAT_TYPE Cik = get_Cik();
	Eigen::SelfAdjointEigenSolver<MAT_TYPE> eigensolver((Cik.selfadjointView<Eigen::Upper>()));
	return (eigensolver.eigenvalues()).minCoeff();
}

	


void cCIJKL<2>::SetZero()
{
	SetArtificialValues(0.);
}

void cCIJKL<2>::SetArtificialValues(dbl c)
{
	for(int i=0; i<num_constants; ++i)
		a[i] = c;
};

Eigen::Matrix<dbl,3,3> cCIJKL<2>::get_Cik() const
{
	Eigen::Matrix<dbl, 3, 3> Cik;
	Cik << 	cxxxx,	cxxyy,	2.*cxxxy,
				0.,	cyyyy,	2.*cyyxy,
				0.,		0.,	4.*cxyxy;
	return Cik;
}

 
void cCIJKL<2>::print() const
{
	Eigen::Matrix<dbl, 3, 3> Cik = get_Cik();
	std::cout << "Cik:\n" << Cik << std::endl;
}


cCIJKL<3>::cCIJKL()
	:	cxxxx (a[xxxx]), // 0
		cyyyy (a[yyyy]), // 1
		czzzz (a[zzzz]), // 2
		cyzyz (a[yzyz]), // 3
		cxzxz (a[xzxz]), // 4
		cxyxy (a[xyxy]), // 5
		cyyzz (a[yyzz]), // 6
		cxxzz (a[xxzz]), // 7
		cxxyy (a[xxyy]), // 8
		cxxyz (a[xxyz]), // 9
		cxxxz (a[xxxz]), //10
		cxxxy (a[xxxy]), //11
		cyyyz (a[yyyz]), //12
		cyyxz (a[yyxz]), //13
		cyyxy (a[yyxy]), //14
		czzyz (a[zzyz]), //15
		czzxz (a[zzxz]), //16
		czzxy (a[zzxy]), //17
		cyzxz (a[yzxz]), //18
		cyzxy (a[yzxy]), //19
		cxzxy (a[xzxy])  //20
{
	for(int i=0; i<num_constants; i++)
		a[i] = 0.;
};

cCIJKL<3>::cCIJKL(const cCIJKL<3> &src)
	:	cxxxx (a[xxxx]), // 0
		cyyyy (a[yyyy]), // 1
		czzzz (a[zzzz]), // 2
		cyzyz (a[yzyz]), // 3
		cxzxz (a[xzxz]), // 4
		cxyxy (a[xyxy]), // 5
		cyyzz (a[yyzz]), // 6
		cxxzz (a[xxzz]), // 7
		cxxyy (a[xxyy]), // 8
		cxxyz (a[xxyz]), // 9
		cxxxz (a[xxxz]), //10
		cxxxy (a[xxxy]), //11
		cyyyz (a[yyyz]), //12
		cyyxz (a[yyxz]), //13
		cyyxy (a[yyxy]), //14
		czzyz (a[zzyz]), //15
		czzxz (a[zzxz]), //16
		czzxy (a[zzxy]), //17
		cyzxz (a[yzxz]), //18
		cyzxy (a[yzxy]), //19
		cxzxy (a[xzxy])  //20
{
	(*this) = src;
};

cCIJKL<3>& cCIJKL<3>::operator=(const cCIJKL<3> &src)
{
	if(this != &src)
	{
		for(int i=0; i<num_constants; ++i)
			a[i] = src.a[i];
	}
	return *this;
};



void cCIJKL<3>::set_constant(dbl twice_dE_over_V, int ii)
{
	assert(ii>=0 && ii<num_constants);
	switch(ii)
	{
	//B
		case 0: //Apply	(1 0 0) 
			//	(0 0 0) 
			//	(0 0 0)
			cxxxx = twice_dE_over_V;
			break;
		case 1: //Apply	(0 0 0) 
			//	(0 1 0) 
			//	(0 0 0)
			cyyyy = twice_dE_over_V;
			break;
		case 2: //Apply	(0 0 0) 
			//	(0 0 0) 
			//	(0 0 1)
			czzzz = twice_dE_over_V;
			break;
	//Shear
		case 3: //Apply	(0 0 0) 
			//	(0 0 1) 
			//	(0 1 0)
			cyzyz = twice_dE_over_V/4.;
			break;
		case 4: //Apply	(0 0 1) 
			//	(0 0 0) 
			//	(1 0 0)
			cxzxz = twice_dE_over_V/4.;
			break;
		case 5: //Apply	(0 1 0) 
			//	(1 0 0) 
			//	(0 0 0)
			cxyxy = twice_dE_over_V/4.;
			break;
	//cross
		case 6: //Apply	(0 0 0) 
			//	(0 1 0) 
			//	(0 0 1)
			cyyzz = (twice_dE_over_V - cyyyy - czzzz)/2.;
			break;
		case 7: //Apply	(1 0 0) 
			//	(0 0 0) 
			//	(0 0 1)
			cxxzz = (twice_dE_over_V - cxxxx - czzzz)/2.;
			break;
		case 8: //Apply	(1 0 0) 
			//	(0 1 0) 
			//	(0 0 0)
			cxxyy = (twice_dE_over_V - cxxxx - cyyyy)/2.;
			break;
	//Dilatency
		case 9: //Apply	(1 0 0) 
			//	(0 0 1) 
			//	(0 1 0)
			cxxyz = (twice_dE_over_V - cxxxx - 4.*cyzyz)/4.;
			break;
		case 10://Apply	(1 0 1) 
			//	(0 0 0) 
			//	(1 0 0)
			cxxxz = (twice_dE_over_V - cxxxx - 4.*cxzxz)/4.;
			break;
		case 11://Apply	(1 1 0) 
			//	(1 0 0) 
			//	(0 0 0)
			cxxxy = (twice_dE_over_V - cxxxx - 4.*cxyxy)/4.;
			break;
		case 12://Apply	(0 0 0) 
			//	(0 1 1) 
			//	(0 1 0)
			cyyyz = (twice_dE_over_V - cyyyy - 4.*cyzyz)/4.;
			break;
		case 13://Apply	(0 0 1) 
			//	(0 1 0) 
			//	(1 0 0)
			cyyxz = (twice_dE_over_V - cyyyy - 4.*cxzxz)/4.;
			break;
		case 14://Apply	(0 1 0) 
			//	(1 1 0) 
			//	(0 0 0)
			cyyxy = (twice_dE_over_V - cyyyy - 4.*cxyxy)/4.;
			break;
		case 15://Apply	(0 0 0) 
			//	(0 0 1) 
			//	(0 1 1)
			czzyz = (twice_dE_over_V - czzzz - 4.*cyzyz)/4.;
			break;
		case 16://Apply	(0 0 1) 
			//	(0 0 0) 
			//	(1 0 1)
			czzxz = (twice_dE_over_V - czzzz - 4.*cxzxz)/4.;
			break;
		case 17://Apply	(0 1 0) 
			//	(1 0 0) 
			//	(0 0 1)
			czzxy = (twice_dE_over_V - czzzz - 4.*cxyxy)/4.;
			break;
	//Last
		case 18://Apply	(0 0 1) 
			//	(0 0 1) 
			//	(1 1 0)
			cyzxz = (twice_dE_over_V - 4.*cyzyz - 4.*cxzxz)/8.;
			break;
		case 19://Apply	(0 1 0) 
			//	(1 0 1) 
			//	(0 1 0)
			cyzxy = (twice_dE_over_V - 4.*cyzyz - 4.*cxyxy)/8.;
			break;
		case 20://Apply	(0 1 1) 
			//	(1 0 0) 
			//	(1 0 0)
			cxzxy = (twice_dE_over_V - 4.*cxzxz - 4.*cxyxy)/8.;
			break;
	}
};

void cCIJKL<3>::set_strain_tensor(dmat &strain_tensor, int ii)
{
	strain_tensor.setZero();
	const dbl strain = 1.;
	switch(ii)
	{
	//B
		case 0: //Apply	(1 0 0) 
			//	(0 0 0) 
			//	(0 0 0)
			strain_tensor(0,0) = strain;
			break;
		case 1: //Apply	(0 0 0) 
			//	(0 1 0) 
			//	(0 0 0)
			strain_tensor(1,1) = strain;
			break;
		case 2: //Apply	(0 0 0) 
			//	(0 0 0) 
			//	(0 0 1)
			strain_tensor(2,2) = strain;
			break;
	//Shear
		case 3: //Apply	(0 0 0) 
			//	(0 0 1) 
			//	(0 1 0)
			strain_tensor(1,2) = strain_tensor(2,1) = strain;
			break;
		case 4: //Apply	(0 0 1) 
			//	(0 0 0) 
			//	(1 0 0)
			strain_tensor(0,2) = strain_tensor(2,0) = strain;
			break;
		case 5: //Apply	(0 1 0) 
			//	(1 0 0) 
			//	(0 0 0)
			strain_tensor(0,1) = strain_tensor(1,0) = strain;
			break;
		//cross
		case 6: //Apply	(0 0 0) 
			//	(0 1 0) 
			//	(0 0 1)
			strain_tensor(1,1) = strain_tensor(2,2) = strain;
			break;
		case 7: //Apply	(1 0 0) 
			//	(0 0 0) 
			//	(0 0 1)
			strain_tensor(0,0) = strain_tensor(2,2) = strain;
			break;
		case 8: //Apply	(1 0 0) 
			//	(0 1 0) 
			//	(0 0 0)
			strain_tensor(0,0) = strain_tensor(1,1) = strain;
			break;
	//Dilatency
		case 9: //Apply	(1 0 0) 
			//	(0 0 1) 
			//	(0 1 0)
			strain_tensor(0,0) = strain_tensor(1,2) = strain_tensor(2,1) = strain;
			break;
		case 10://Apply	(1 0 1) 
			//	(0 0 0) 
			//	(1 0 0)
			strain_tensor(0,0) = strain_tensor(0,2) = strain_tensor(2,0) = strain;
			break;
		case 11://Apply	(1 1 0) 
			//	(1 0 0) 
			//	(0 0 0)
			strain_tensor(0,0) = strain_tensor(0,1) = strain_tensor(1,0) = strain;
			break;
		case 12://Apply	(0 0 0) 
			//	(0 1 1) 
			//	(0 1 0)
			strain_tensor(1,1) = strain_tensor(1,2) = strain_tensor(2,1) = strain;
			break;
		case 13://Apply	(0 0 1) 
			//	(0 1 0) 
			//	(1 0 0)
			strain_tensor(1,1) = strain_tensor(0,2) = strain_tensor(2,0) = strain;
			break;
		case 14://Apply	(0 1 0) 
			//	(1 1 0) 
			//	(0 0 0)
			strain_tensor(1,1) = strain_tensor(0,1) = strain_tensor(1,0) = strain;
			break;
		case 15://Apply	(0 0 0) 
			//	(0 0 1) 
			//	(0 1 1)
			strain_tensor(2,2) = strain_tensor(1,2) = strain_tensor(2,1) = strain;
			break;
		case 16://Apply	(0 0 1) 
			//	(0 0 0) 
			//	(1 0 1)
			strain_tensor(2,2) = strain_tensor(0,2) = strain_tensor(2,0) = strain;
			break;
		case 17://Apply	(0 1 0) 
			//	(1 0 0) 
			//	(0 0 1)
			strain_tensor(2,2) = strain_tensor(0,1) = strain_tensor(1,0) = strain;
			break;
	//Last
		case 18://Apply	(0 0 1) 
			//	(0 0 1) 
			//	(1 1 0)
			strain_tensor(0,2) = strain_tensor(2,0) = strain_tensor(1,2) = strain_tensor(2,1) = strain;
			break;
		case 19://Apply	(0 1 0) 
			//	(1 0 1) 
			//	(0 1 0)
			strain_tensor(0,1) = strain_tensor(1,0) = strain_tensor(1,2) = strain_tensor(2,1) = strain;
			break;
		case 20://Apply	(0 1 1) 
			//	(1 0 0) 
			//	(1 0 0)
			strain_tensor(0,1) = strain_tensor(1,0) = strain_tensor(0,2) = strain_tensor(2,0) = strain;
			break;
	}
}






 
dbl cCIJKL<3>::CalculateBulkModulus() const
{
	return (cxxxx + 2*cxxyy + 2*cxxzz + cyyyy + 2*cyyzz + czzzz)/9.;
}

dbl cCIJKL<3>::CalculateShearModulus() const
{
	return (cxxxx - cxxyy - cxxzz + 3*cxyxy + 3*cxzxz + cyyyy - cyyzz + 3*cyzyz + czzzz)/15.;
}

void cCIJKL<3>::CalculateRotationalAverages(dbl &B, dbl &GDC, dbl &GAC, dbl &UDC, dbl &UAC, dbl &DDC, dbl &DAC, dbl &TDC, dbl &TAC) const
{
	B = CalculateBulkModulus();

	GDC = CalculateShearModulus();
	GAC = sqrt((8*std::pow(cxxxx,2) + 60*std::pow(cxxxy,2) + 60*std::pow(cxxxz,2) + 23*std::pow(cxxyy,2) + 60*std::pow(cxxyz,2) - 14*cxxyy*cxxzz + 23*std::pow(cxxzz,2) + 
		12*cxxyy*cxyxy + 12*cxxzz*cxyxy + 72*std::pow(cxyxy,2) + 180*std::pow(cxzxy,2) + 180*std::pow(cyzxy,2) + 12*cxxyy*cxzxz + 12*cxxzz*cxzxz - 
		36*cxyxy*cxzxz + 72*std::pow(cxzxz,2) + 180*std::pow(cyzxz,2) + 60*std::pow(cyyxy,2) - 60*cxxxz*cyyxz + 60*std::pow(cyyxz,2) - 16*cxxyy*cyyyy + 
		14*cxxzz*cyyyy - 12*cxyxy*cyyyy - 12*cxzxz*cyyyy + 8*std::pow(cyyyy,2) - 60*cxxyz*cyyyz + 60*std::pow(cyyyz,2) - 14*cxxyy*cyyzz - 
		14*cxxzz*cyyzz + 12*cxyxy*cyyzz + 12*cxzxz*cyyzz - 16*cyyyy*cyyzz + 23*std::pow(cyyzz,2) + 12*cxxyy*cyzyz + 12*cxxzz*cyzyz - 36*cxyxy*cyzyz - 
		36*cxzxz*cyzyz - 12*cyyyy*cyzyz + 12*cyyzz*cyzyz + 72*std::pow(cyzyz,2) - 60*cyyxy*czzxy + 60*std::pow(czzxy,2) - 60*cxxxy*(cyyxy + czzxy) - 
		60*cxxxz*czzxz - 60*cyyxz*czzxz + 60*std::pow(czzxz,2) - 60*cxxyz*czzyz - 60*cyyyz*czzyz + 60*std::pow(czzyz,2) + 14*cxxyy*czzzz - 
		16*cxxzz*czzzz - 12*cxyxy*czzzz - 12*cxzxz*czzzz + cyyyy*czzzz - 16*cyyzz*czzzz - 12*cyzyz*czzzz + 8*std::pow(czzzz,2) + 
		cxxxx*(-16*cxxyy - 16*cxxzz - 12*cxyxy - 12*cxzxz + cyyyy + 14*cyyzz - 12*cyzyz + czzzz))/1575.);

	UDC = (3*cxxxx + 2*cxxyy + 2*cxxzz + 4*cxyxy + 4*cxzxz + 3*cyyyy + 2*cyyzz + 4*cyzyz + 3*czzzz)/15.;
	UAC = sqrt((16*(7*std::pow(cxxxx,2) + 25*std::pow(cxxxy,2) + 25*std::pow(cxxxz,2) + 2*std::pow(cxxyy,2) + 5*std::pow(cxxyz,2) - cxxyy*cxxzz + 2*std::pow(cxxzz,2) + 
		8*cxxyy*cxyxy - 2*cxxzz*cxyxy + 8*std::pow(cxyxy,2) + 20*cxxyz*cxzxy + 20*std::pow(cxzxy,2) + 20*cxxxz*cyzxy + 20*std::pow(cyzxy,2) - 
		2*cxxyy*cxzxz + 8*cxxzz*cxzxz - 4*cxyxy*cxzxz + 8*std::pow(cxzxz,2) + 20*std::pow(cyzxz,2) + 20*cyzxz*cyyxy + 25*std::pow(cyyxy,2) + 
		10*cxxxz*cyyxz + 20*cyzxy*cyyxz + 5*std::pow(cyyxz,2) + cxxyy*cyyyy - 4*cxxzz*cyyyy + 2*cxyxy*cyyyy - 8*cxzxz*cyyyy + 7*std::pow(cyyyy,2) + 
		10*cxxyz*cyyyz + 20*cxzxy*cyyyz + 25*std::pow(cyyyz,2) - cxxyy*cyyzz - cxxzz*cyyzz - 2*cxyxy*cyyzz - 2*cxzxz*cyyzz + cyyyy*cyyzz + 
		2*std::pow(cyyzz,2) - 2*cxxyy*cyzyz - 2*cxxzz*cyzyz - 4*cxyxy*cyzyz - 4*cxzxz*cyzyz + 2*cyyyy*cyzyz + 8*cyyzz*cyzyz + 8*std::pow(cyzyz,2) + 
		20*cyzxz*czzxy + 10*cyyxy*czzxy + 5*std::pow(czzxy,2) + 10*cxxxy*(2*cyzxz + 3*cyyxy + czzxy) + 30*cxxxz*czzxz + 20*cyzxy*czzxz + 
		10*cyyxz*czzxz + 25*std::pow(czzxz,2) + 10*cxxyz*czzyz + 20*cxzxy*czzyz + 30*cyyyz*czzyz + 25*std::pow(czzyz,2) + 
		cxxxx*(cxxyy + cxxzz + 2*(cxyxy + cxzxz - 3*cyyyy - 2*cyyzz - 4*cyzyz - 3*czzzz)) - 4*cxxyy*czzzz + cxxzz*czzzz - 8*cxyxy*czzzz + 
		2*cxzxz*czzzz - 6*cyyyy*czzzz + cyyzz*czzzz + 2*cyzyz*czzzz + 7*std::pow(czzzz,2)))/1575.);

	DDC = 0.;
	DAC = sqrt((5*std::pow(cxxxx,2) + 23*std::pow(cxxxy,2) + 23*std::pow(cxxxz,2) + 3*std::pow(cxxyy,2) + 7*std::pow(cxxyz,2) - cxxyy*cxxzz + 3*std::pow(cxxzz,2) + 
		12*cxxyy*cxyxy - 2*cxxzz*cxyxy + 12*std::pow(cxyxy,2) + 28*cxxyz*cxzxy + 28*std::pow(cxzxy,2) + 4*cxxxz*cyzxy + 28*std::pow(cyzxy,2) - 
		2*cxxyy*cxzxz + 12*cxxzz*cxzxz - 4*cxyxy*cxzxz + 12*std::pow(cxzxz,2) + 28*std::pow(cyzxz,2) + 4*cyzxz*cyyxy + 23*std::pow(cyyxy,2) + 
		2*cxxxz*cyyxz + 28*cyzxy*cyyxz + 7*std::pow(cyyxz,2) - cxxyy*cyyyy - 2*cxxzz*cyyyy - 2*cxyxy*cyyyy - 4*cxzxz*cyyyy + 5*std::pow(cyyyy,2) + 
		2*cxxyz*cyyyz + 4*cxzxy*cyyyz + 23*std::pow(cyyyz,2) - cxxyy*cyyzz - cxxzz*cyyzz - 2*cxyxy*cyyzz - 2*cxzxz*cyyzz - cyyyy*cyyzz + 
		3*std::pow(cyyzz,2) - 2*cxxyy*cyzyz - 2*cxxzz*cyzyz - 4*cxyxy*cyzyz - 4*cxzxz*cyzyz - 2*cyyyy*cyzyz + 12*cyyzz*cyzyz + 12*std::pow(cyzyz,2) + 
		28*cyzxz*czzxy + 2*cyyxy*czzxy + 7*std::pow(czzxy,2) + 2*cxxxy*(2*cyzxz + 3*cyyxy + czzxy) + 6*cxxxz*czzxz + 4*cyzxy*czzxz + 2*cyyxz*czzxz + 
		23*std::pow(czzxz,2) + 2*cxxyz*czzyz + 4*cxzxy*czzyz + 6*cyyyz*czzyz + 23*std::pow(czzyz,2) - 2*cxxyy*czzzz - cxxzz*czzzz - 4*cxyxy*czzzz - 
		2*cxzxz*czzzz - 3*cyyyy*czzzz - cyyzz*czzzz - 2*cyzyz*czzzz + 5*std::pow(czzzz,2) - 
		cxxxx*(cxxyy + cxxzz + 2*cxyxy + 2*cxzxz + 3*cyyyy + 2*cyyzz + 4*cyzyz + 3*czzzz))/315.);

	TDC = 0.;
	TAC = sqrt((std::pow(cxxxx,2) + 7*std::pow(cxxxy,2) + 7*std::pow(cxxxz,2) + 3*std::pow(cxxyy,2) + 8*std::pow(cxxyz,2) + 3*std::pow(cxxzz,2) + 2*cxxzz*cxyxy + 
		9*std::pow(cxyxy,2) + 23*std::pow(cxzxy,2) + 23*std::pow(cyzxy,2) - 5*cxyxy*cxzxz + 9*std::pow(cxzxz,2) + 2*cxxxy*cyzxz + 23*std::pow(cyzxz,2) - 
		6*cxxxy*cyyxy + 2*cyzxz*cyyxy + 7*std::pow(cyyxy,2) + 8*std::pow(cyyxz,2) + 2*cxxzz*cyyyy - cxyxy*cyyyy - 2*cxzxz*cyyyy + std::pow(cyyyy,2) + 
		7*std::pow(cyyyz,2) - 2*cxxzz*cyyzz + 2*cxyxy*cyyzz + 2*cxzxz*cyyzz - 2*cyyyy*cyyzz + 3*std::pow(cyyzz,2) + 2*cxxzz*cyzyz - 5*cxyxy*cyzyz - 
		5*cxzxz*cyzyz - cyyyy*cyzyz + 9*std::pow(cyzyz,2) - cxxxx*(2*cxxyy + 2*cxxzz + cxyxy + cxzxz - 2*cyyzz + 2*cyzyz) - 
		4*(2*cxxxy + cyzxz + 2*cyyxy)*czzxy + 8*std::pow(czzxy,2) + 2*cxxxz*(cyzxy - 4*cyyxz - 3*czzxz) - 8*cyyxz*czzxz + 7*std::pow(czzxz,2) + 
		2*cyzxy*(-2*cyyxz + czzxz) - 6*cyyyz*czzyz + 7*std::pow(czzyz,2) + 2*cxzxy*(cyyyz + czzyz) - 4*cxxyz*(cxzxy + 2*(cyyyz + czzyz)) - 
		(2*cxxzz + 2*cxyxy + cxzxz + 2*cyyzz + cyzyz)*czzzz + std::pow(czzzz,2) + 2*cxxyy*(-cxxzz + cxzxz - cyyyy - cyyzz + cyzyz + czzzz))/315.);
}

dbl cCIJKL<3>::calc_weakest_deformation() const
{
	typedef Eigen::Matrix<dbl, Eigen::Dynamic, Eigen::Dynamic> MAT_TYPE;
	MAT_TYPE Cik = get_Cik();
	Eigen::SelfAdjointEigenSolver<MAT_TYPE> eigensolver((Cik.selfadjointView<Eigen::Upper>()));
	return (eigensolver.eigenvalues()).minCoeff();
}




void cCIJKL<3>::SetZero()
{
	SetArtificialValues(0.);
}

void cCIJKL<3>::SetArtificialValues(dbl c)
{
	for(int i=0; i<num_constants; ++i)
		a[i] = c;
};

Eigen::Matrix<dbl,6,6> cCIJKL<3>::get_Cik() const
{
	Eigen::Matrix<dbl, 6, 6> Cik;
	Cik << 	cxxxx,	cxxyy,	cxxzz,	2.*cxxyz,	2.*cxxxz,	2.*cxxxy,
				0.,	cyyyy,	cyyzz,	2.*cyyyz,	2.*cyyxz,	2.*cyyxy,
				0.,		0.,	czzzz,	2.*czzyz, 	2.*czzxz,	2.*czzxy,
				0.,		0.,		0.,	4.*cyzyz, 	4.*cyzxz,	4.*cyzxy,
				0.,		0.,		0.,			0.,	4.*cxzxz,	4.*cxzxy,
				0.,		0.,		0.,			0.,		 	0.,	4.*cxyxy;

	return Cik;
}

 
void cCIJKL<3>::print() const
{
	Eigen::Matrix<dbl, 6, 6> Cik = get_Cik();
	std::cout << "Cik:\n" << Cik << std::endl;
}

















