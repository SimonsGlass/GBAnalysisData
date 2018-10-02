#ifndef STD_DATA_H
#define STD_DATA_H

#include "../Computers/cijkl.h"
#include "../Computers/MatrixInterface.h"

namespace LiuJamming
{

template<int Dim>
class CStdData
{
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;
public:
	int NPp;
	int Nc;
	dbl Volume;

	dbl Energy;
	dbl Pressure;
	dmat Stress;
	dbl MaxGrad;

	cCIJKL<Dim> cijkl;
	MatrixInterface<dbl> H;
	
	inline int Nciso()		const { return Dim*(NPp-1); };
	inline int NcmNciso()	const { return Nc - Nciso(); };
	inline dbl Z()			const { return 2.*((dbl)Nc)/((dbl)NPp); };
	inline dbl deltaZ()		const { return 2.*((dbl)NcmNciso())/((dbl)NPp); };

	void SetZero()
	{
		NPp = Nc = 0;
		Volume = Energy = Pressure = MaxGrad = 0.;
		Stress = dmat::Zero();
		cijkl.SetZero();
		//Set the matrix to zero???
	};

	void Print()
	{
		printf("----------------------------------------------------------------\n");
		printf("Printing standard data:\n");
		printf("\tNPp         = %i\n", NPp);
		printf("\tNc          = %i\n", Nc);
		printf("\tNcmNciso    = %i\n", NcmNciso());
		printf("\tDelta Z     = % e\n", deltaZ());
		printf("\tEnergy      = % e\n", Energy);
		printf("\tPressure    = % e\n", Pressure);
		printf("\tMaxGrad     = % e\n", MaxGrad);
		printf("\tStress Tensor: \n");
		for(int d1=0; d1<Dim; ++d1)
		{
			printf("\t             ");
			for(int d2=0; d2<Dim; ++d2) printf(" % e  ", Stress(d1,d2));
			printf("\n");
		}
		printf("\tElastic Constants: \n");
		Eigen::MatrixXd Cik = cijkl.get_Cik();

		char dot[32];
		sprintf(dot,"      .      ");
		for(int d1=0; d1<Cik.rows(); ++d1)
		{
			printf("\t             ");
			for(int d2=0; d2<d1; ++d2) printf(" %s  ", dot);
			for(int d2=d1; d2<Cik.cols(); ++d2) printf(" % e  ", Cik(d1,d2));
			printf("\n");
		}
		printf("----------------------------------------------------------------\n\n");
	}
};

}

#endif //STD_DATA_H

