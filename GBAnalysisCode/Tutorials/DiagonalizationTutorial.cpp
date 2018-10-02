#include <iostream>
#include "State/StaticState.h"
#include "Potentials/Potentials.h"
#include "Boundaries/Boxes.h"
#include "Resources/Exception.h"
#include "Computers/StaticComputer.h"
#include "Computers/MatrixInterface.h"
#include "Minimization/minimizer.h"
using namespace std;
using namespace LiuJamming;

#define DIM 2
typedef CStaticState<DIM> STATE;
typedef CStaticComputer<DIM> COMP;
typedef CBondList<DIM> BLIST;
typedef CSimpleMinimizer<DIM> SMINER;
typedef MatrixInterface<dbl> MI;
typedef Eigen::Matrix<dbl, DIM, 1> dvec;

/*	Test 1: Quick Diagonalization
 *
 *	This shows how to find and print the first 100 eigenvalues in very few lines.
 */
void test1(STATE &s)
{
	COMP c(s);						//Create a Computer
	c.StdPrepareSystem();			//Prepare system for standard calculations (this sets the internal BondList and removes rattlers).
	c.CalculateStdData();			//This also creates a Matrix Interface and sets the Hessian. The Matrix Interface can be accessed in c.Data.H
	c.Data.Print();					//Print the data  (not necessary)
	MI &H = c.Data.H;				//For convenience, create a reference of the Matrix Interface (not necessary)
	H.VDiagonalize(100);			//Diagonalize the hessian to find the first 100 eigenvalues, and print the report to the screen
}

/*	Test 2: Reading the data
 *
 *	This shows how to read the eigenvalues and eigenvectors
 */
void test2(STATE &s)
{
	//Setup
	COMP c(s);						//Create a Computer
	c.StdPrepareSystem();			//Prepare system for standard calculations (this sets the internal BondList and removes rattlers).
	c.CalculateStdData();			//This also creates a Matrix Interface and sets the Hessian. The Matrix Interface can be accessed in c.Data.H
	c.Data.Print();					//Print the data
	MI &H = c.Data.H;				//For convenience, create a reference of the Matrix Interface.

	//Diagonalize the system
	H.Diagonalize(100);				//Diagonalize the hessian to find the first 100 eigenvalues

	//Print out the eigenvalues and the first 3 components of the corresponding eigenvectors
	//
	//	- The Matrix Interface has a parameter call num_converged. ALWAYS use this parameter when looping over eigenvalues/eigenvectors.
	//		Do not assume that the number of converged eigenvalues equals the number of requested eigenvalues.
	//	- The eigenvalues are stored in an array of length H.num_converged called H.Eigenvalues
	//	- The eigenvectors are stored in an array of length H.num_converged*H.A.rows() called H.Eigenvectors
	//			The first eigenvector is stored in the first H.A.rows() entries, etc.
	//
	//	- The Eigenvalues and Eigenvectors can be accessed through the functions
	//			GetVal(int m) const;			//m is the mode index
	//			GetVec(int m, int i) const;		//m is the mode index, i is the component of the eigenvector
	for(int i=0; i<H.num_converged; ++i)
	{
		printf("Eigenvalue[%5i] = % e     The first three components of the eigenvector are % e % e % e\n", i, H.GetVal(i),H.GetVec(i,0),H.GetVec(i,1),H.GetVec(i,2));
	}
}

/*	Test 3: Options
 *
 *	This shows some of the options that can be used
 */
void test3(STATE &s)
{
	//Setup
	COMP c(s);						//Create a Computer
	c.StdPrepareSystem();			//Prepare system for standard calculations (this sets the internal BondList and removes rattlers).
	c.CalculateStdData();			//This also creates a Matrix Interface and sets the Hessian. The Matrix Interface can be accessed in c.Data.H
	c.Data.Print();					//Print the data
	MI &H = c.Data.H;				//For convenience, create a reference of the Matrix Interface.

	//The default/simplest way to diagonalize the system is to call:
	//				H.Diagonalize();

	//There are currently 3 primary options
	//	1. Choose the number of eigenvalues/vectors to calculate.
	//		- By default, all eigenvalues are calculated (except the largest, for unknown reasons having to do with ARPACK being stupid)
	//		- This is easily specified by as in input argumen to H.Diagonalize. For example, to calculate the lowest 100 eigenvalues:
	//				H.Diagonalize(100);
	//		- Alternatively, one can call H.SetNumRequest(100); prior to calling the H.Diagonalize() routine
	//
	//	2. Choose to print the default diagonalization report after diagonalization is complete.
	//		- This is easily done by including a "V" (for verbose) in the function call, e.g.
	//				H.VDiagonalize(100);
	//		- After diagonalization, one can also call 
	//				H.Report(45);
	//			to print the report on the first 45 modes. If the argument is omitted, 30 modes are reported on by default.
	//
	//	3. Choose to diagonalize using the Shift and Invert functionality of ARPACK (see http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html for a simple explanation).
	//		- WARNING: While this is significantly faster (up to ~25X speedup), the accuracy of the eigenvalues is significantly decreased. Use at your own risk.
	//		- This is easily done by includeing "_SI" at the end of the function call, e.g.
	//				H.Diagonalize_SI(100);
	//		- Alternatively, one can call H.SetShiftAndInvert(); prior to calling the H.Diagonalize routine
	//
	//	In Summary, there are currently 8 primary routines that can be called:
	//				H. Diagonalize   ();
	//				H. Diagonalize   (100);
	//				H.VDiagonalize   ();
	//				H.VDiagonalize   (100);
	//				H. Diagonalize_SI();
					H. Diagonalize_SI(100);
	//				H.VDiagonalize_SI();
	//				H.VDiagonalize_SI(100);
	//
	//Additionally, one can choose to only compute the eigenvectors by calling
	//				H.SetComputeVecs(false);
	//		prior to the diagonalization routine.
	//
	//TODO:
	//	- Option for changing the region where eigenvalues are searched (either via "which" or via "sigma" in the shift and invert mode).
	//	- Easily calculate the ENTIRE spectrum... this might not be necessary bc probably only applicable to small systems when Eigen dense diagonalization would suffice
	//	- Diagonalize by first converting to an Eigen dense matrix.
	//
	//
	H.Report(100);

}

/*	Test 1: Complex Diagonalization
 *
 *	This shows how to find and print the first 100 eigenvalues in very few lines.
 */
void test4(STATE &s)
{
	COMP c(s);						//Create a Computer
	c.StdPrepareSystem();			//Prepare system for standard calculations (this sets the internal BondList and removes rattlers).
	c.CalculateStdData();			//This also creates a Matrix Interface and sets the Hessian. The Matrix Interface can be accessed in c.Data.H
	c.Data.Print();					//Print the data  (not necessary)

	MatrixInterface<cdbl> cH;
	dvec k = dvec::Zero();k[0] = M_PI;

	//This is currently borken.
//	c.Bonds.ComputeHessian_BZ(cH.A, k);
//	cH.VDiagonalize(100);			//Diagonalize the hessian to find the first 100 eigenvalues, and print the report to the screen
}

void run(int N, dbl phi, int seed, int test_num=1)
{
	//Create the system
	STATE s(N);
	s.RandomizePositions(seed);
	s.SetRadiiPolyUniform();
	s.SetPackingFraction(phi);

	//Minimize the energy:
	COMP c(s);
	SMINER miner(c, SMINER::FIRE); //Note: the miner class should be able to create the Computer (i.e. just pass it the state)
	
	//Run tests
	switch(test_num)
	{
		case 1:		test1(s); break;
		case 2:		test2(s); break;
		case 3:		test3(s); break;
		case 4:		test4(s); break;
	}
}






int main(int argc, char* argv[])
{
	int N = 256;			//n
	dbl phi = 0.9;			//f
	int r = 1;				//r
	int test_num = 1;		//t

	int c;
	while((c=getopt(argc, argv, "n:f:r:t:")) != -1)
		switch(c)
		{
			case 'n':	N = atoi(optarg); break;
			case 'f':	phi = atof(optarg); break;
			case 'r':	r = atoi(optarg); break;
			case 't':	test_num = atoi(optarg); break;
			case '?':	
				if(optopt == 'c') 
					std::cerr << "Option -" << optopt << "requires an argument.\n";
				else if(isprint(optopt)) 
					std::cerr << "Unknown opton '-" << optopt << "'.\n";
				else
					std::cerr << "Unknown option character '\\x" << optopt << "'.\n";
				return 1;
			default:
				abort();
		}
	for(int index = optind; index < argc; index++)
		std::cerr << "Non-option argument " << argv[index] << "\n";


	printf("Begin Diagonalization Tutorial with N=%i, phi=%f, r=%i\n", N, phi, r);
	run(N, phi, r, test_num);
}








