#ifndef MATRIX_INTERFACE_H
#define MATRIX_INTERFACE_H

#include "../Resources/std_include.h"

#ifndef DONT_USE_ARPACK
	#include "arch.h"
	#include "arscomp.h"
	#include "arssym.h"
#endif

#include <Eigen/Sparse>

#ifndef DONT_USE_SUITESPARSE
	#include <Eigen/UmfPackSupport>
#endif

namespace LiuJamming
{



template <class T> class MatrixInterface;

//Functions to interface with ARPACK++
template <class T> class ARWrapper {}; 


#ifndef DONT_USE_ARPACK

template <>
class ARWrapper<dbl> 
{
public: 
	ARSymStdEig<dbl,MatrixInterface<dbl> > Solver;
	static const std::string which;
	static const std::string which2;
};
const std::string ARWrapper<dbl>::which("SA");
const std::string ARWrapper<dbl>::which2("LM"); // CHANGED TO lowest magnitude...won't pick up negative modes necessarily!

template <>
class ARWrapper<cdbl> 
{
public: 
	ARCompStdEig<dbl,MatrixInterface<cdbl> > Solver;
	static const std::string which;
	static const std::string which2;
};
const std::string ARWrapper<cdbl>::which("SM");
const std::string ARWrapper<cdbl>::which2("LM");

#endif //DONT_USE_ARPACK


template <class T>
static void set_identity(Eigen::SparseMatrix<T> &Identity, int nrows)
{
	assert(Identity.rows() == nrows);
	assert(Identity.cols() == nrows);

	std::vector< Eigen::Triplet<T> > triplets;
	triplets.reserve(nrows);
	for(int i=0; i<nrows; ++i)
		triplets.push_back(Eigen::Triplet<T>(i,i,((T)1.)));
	Identity.setFromTriplets(triplets.begin(), triplets.end());
}

template<typename T> void PrintEigenvalues(T *Eigenvalues, Eigen::Matrix<T,Eigen::Dynamic,1> const &Residual, int num_print, int num_converged) {};


template<> void PrintEigenvalues(dbl *Eigenvalues, Eigen::Matrix<dbl,Eigen::Dynamic,1> const &Residual, int num_print, int num_converged)
{
	printf("     i:       eigenvalue         residual\n");
	for(int i=0; i<num_print; ++i)
		printf("%6i:    % e    % e\n", i, Eigenvalues[i], Residual[i]);
};

class cplxevdata{
public:
	cdbl ev, r;
	cplxevdata(cdbl _ev, cdbl _r)
		: ev(_ev), r(_r) {};
};
bool cplxevdata_lessthan(cplxevdata const &c1, cplxevdata const &c2)
{	return (std::real(c1.ev) < std::real(c2.ev)); }

template<> void PrintEigenvalues(cdbl *Eigenvalues, Eigen::Matrix<cdbl,Eigen::Dynamic,1> const &Residual, int num_print, int num_converged)
{
	vector<cplxevdata> evtemp;
	evtemp.reserve(num_converged);
	for(int i=0; i<num_converged; ++i)
		evtemp.push_back( cplxevdata(Eigenvalues[i], Residual[i]) );
	std::sort(evtemp.begin(), evtemp.end(), cplxevdata_lessthan);
	printf("     i:       eigenvalue                              residual\n");
	for(int i=0; i<num_print; ++i)
	{
		printf("%6i:    % e  +  % e i     % e  +  % e i\n", i, std::real(evtemp[i].ev), std::imag(evtemp[i].ev), std::real(evtemp[i].r), std::imag(evtemp[i].r));
	}
};
	








template<class T>
class MatrixInterface 
{
private:
	typedef Eigen::SparseMatrix<T> EMatrix;
	typedef Eigen::Matrix<T,Eigen::Dynamic,1> EVector;

	//internal data
	double diagonalization_time;
	long int Mv_counter;
	long int num_Mv_calls;

	//long int OPv_counter;
	//long int num_OPv_calls;

public:
	EMatrix A; //The matrix

	T *Eigenvalues;
	T *Eigenvectors;
	int num_request;
	int num_converged;
	bool compute_vecs;
	bool verbose_diagonalize;
	bool UseShiftAndInvert;

#ifndef DONT_USE_SUITESPARSE
	//Solver for UMFPACK
	Eigen::UmfPackLU<EMatrix> UMFsolver;
#endif

	MatrixInterface()
		: num_request(0),num_converged(0),compute_vecs(1),verbose_diagonalize(0), UseShiftAndInvert(0), Eigenvalues(NULL), Eigenvectors(NULL)
	{
	};

	MatrixInterface(EMatrix const &mat)
		: num_request(0),num_converged(0),compute_vecs(1),verbose_diagonalize(0), UseShiftAndInvert(0), Eigenvalues(NULL), Eigenvectors(NULL), 
		  A(mat)
	{
	};

	~MatrixInterface()
	{
		ClearEigenstuff();
	}


//	void set_default_options();
//	void set_verbose_invert     (int yesno ){	verbose_invert      = yesno;	};
//	void set_verbose_diagonalize(int yesno ){	verbose_diagonalize = yesno;	};
//	void set_diagonalize_region (int region){	diagonalize_region  = region;	};
//	void set_diagonalize_numev  (int numev ){	diagonalize_numev   = numev;	};

	void clear();
	int size() const { return A.rows(); };

private:
	void MultMv_ARPACK(T *v, T *w)
	{
		//Count the number of times this method is called. Print every 100 times
		++Mv_counter;
		if(Mv_counter%1000 == 0) { printf(" %li", Mv_counter); fflush(stdout); }

		Eigen::Map<EVector> Ev(v,A.rows());
		Eigen::Map<EVector> Ew(w,A.rows());
		Ew = A * Ev;
		//Ew = A.selfadjointView<>() * Ev; //Not sure why this doesn't work.
	};

	void MultMinvv_ARPACK(T *v, T *w)
	{
		//Count the number of times this method is called. Print every 100 times
		++Mv_counter;
		if(Mv_counter%100 == 0) { printf(" %li", Mv_counter); fflush(stdout); }

		solve_Mx_equals_b(w,v);
	};

public:
	void MultMv(EVector const &v, EVector &w) const
	{
		w = A * v;
	}

	void MultMv(T const *v, T *w) const
	{
		Eigen::Map<const EVector> Ev(v,A.rows());
		Eigen::Map<      EVector> Ew(w,A.rows());
		Ew = A * Ev;
	};

	T MultvMw(EVector const &v, EVector const &w) const
	{
		EVector Dw = A * w;
		return v.dot(Dw);
	}

	T MultvMw(T const *v, T const *w) const 
	{
		Eigen::Map<const EVector> Ev(v,A.rows());
		Eigen::Map<const EVector> Ew(w,A.rows());
		return MultvMw(Ev,Ew);
	}

	void ClearEigenstuff()
	{
		if(Eigenvalues != NULL)
			delete[] Eigenvalues;
		if(Eigenvectors != NULL)
			delete[] Eigenvectors;
	};

	void SetShiftAndInvert(bool _UseShiftAndInvert=true)	{ UseShiftAndInvert=_UseShiftAndInvert;};
	void SetNumRequest(int _num_request=0)					{ num_request = _num_request;};
	void SetComputeVecs(bool _compute_vecs=true)			{ compute_vecs = _compute_vecs;};

	void Diagonalize();
	void Diagonalize(int _num_request);
	void VDiagonalize();
	void VDiagonalize(int _num_request);
	void Diagonalize_SI();
	void Diagonalize_SI(int _num_request);
	void VDiagonalize_SI();
	void VDiagonalize_SI(int _num_request);

	void Report(int num_print = 30) const;
	void CalculateResidual(EVector &Residual) const;

	T    GetVal(int m) const;
	T    GetVec(int m, int i) const;

	void LUdecomp();
	void solve_Mx_equals_b(EVector &X, EVector const &B);
	void solve_Mx_equals_b(T *x, T const *b);

//	void save_matrix_nc(char *filename);
//	void load_matrix_nc(char *filename);

};

template<typename T>
void MatrixInterface<T>::Diagonalize()
{
#ifndef DONT_USE_ARPACK
	ClearEigenstuff();

	std::cout << "Begin Diagonalization..." << std::endl;
	time_t start,end;
	time(&start);

	if(UseShiftAndInvert)
	{
		cout << "Shift and Invert mode is on. First perform an LU decomposition... ";
		LUdecomp();
		cout << "Done." << endl;
	}

	Mv_counter = 0;
	std::cout << "calls to MultMv:" << std::endl;
	fflush(stdout);

	if(num_request <= 0)	num_request	= A.rows()-1;
	else					num_request	= std::min(A.rows()-1, num_request);
	Eigenvalues = new T[num_request];
	if(compute_vecs)
		Eigenvectors = new T[num_request*A.rows()];

	ARWrapper<T> Wrapper;
	if(UseShiftAndInvert)
	{
		Wrapper.Solver.DefineParameters(A.rows(),num_request,this,&MatrixInterface<T>::MultMinvv_ARPACK,(char*)ARWrapper<T>::which2.c_str());
		Wrapper.Solver.SetShiftInvertMode(0.,this,&MatrixInterface<T>::MultMinvv_ARPACK);
	}else
		Wrapper.Solver.DefineParameters(A.rows(),num_request,this,&MatrixInterface<T>::MultMv_ARPACK,(char*)ARWrapper<T>::which.c_str(), 0, 0.0, 100000000);//Wrapper.Solver.DefineParameters(A.rows(),num_request,this,&MatrixInterface<T>::MultMv_ARPACK,(char*)ARWrapper<T>::which.c_str());
	if(compute_vecs)
		num_converged = Wrapper.Solver.EigenValVectors(Eigenvectors,Eigenvalues);
	else
		num_converged = Wrapper.Solver.Eigenvalues(Eigenvalues);

	std::cout << std::endl;
	time(&end);
	diagonalization_time = difftime(end,start);
	num_Mv_calls = Mv_counter;

	if(verbose_diagonalize)
		Report();

#else
	
	printf("ERROR: ARPACK is turned off and so I cannot diagonalize this matrix for you.\n\
			\tComment out the line\n\
			\t\t#define DONT_USE_ARPACK\n\
			\tin Resources/std_include.h\n");
	fflush(stdout);
	assert(false);
#endif //DONT_USE_ARPACK
};

template<typename T>
void MatrixInterface<T>::Diagonalize(int _num_request)
{
	num_request = _num_request;
	Diagonalize();
}

template<typename T>
void MatrixInterface<T>::VDiagonalize()
{
	verbose_diagonalize = true;
	Diagonalize();
}

template<typename T>
void MatrixInterface<T>::VDiagonalize(int _num_request)
{
	verbose_diagonalize = true;
	Diagonalize(_num_request);
}

template<typename T>
void MatrixInterface<T>::Diagonalize_SI()
{
	SetShiftAndInvert();
	Diagonalize();
}

template<typename T>
void MatrixInterface<T>::Diagonalize_SI(int _num_request)
{
	SetShiftAndInvert();
	Diagonalize(_num_request);
}

template<typename T>
void MatrixInterface<T>::VDiagonalize_SI()
{
	SetShiftAndInvert();
	VDiagonalize();
}

template<typename T>
void MatrixInterface<T>::VDiagonalize_SI(int _num_request)
{
	SetShiftAndInvert();
	VDiagonalize(_num_request);
}


template<typename T>
void MatrixInterface<T>::Report(int num_print) const
{
	if(num_print < 0)
		num_print = num_converged;
	std::cout << "\nGenerating eigenvalue report..." << std::endl;

//	std::cout << "\tUsing ARPACK++ class ARSymStdEig" << std::endl;
//	std::cout << "\tReal symmetric eigenvalue problem: D*x - lambda*x" << std::endl;
	std::cout << "\tDimension of the system             : " << A.rows()				<< std::endl;
	std::cout << "\tNumber of 'requested' eigenvalues   : " << num_request			<< std::endl;
	std::cout << "\tNumber of 'converged' eigenvalues   : " << num_converged		<< std::endl;
	std::cout << "\tSeconds to diagonalize              : " << diagonalization_time	<< std::endl;
	std::cout << "\tNumber of calls to MultMv_ARPACK    : " << num_Mv_calls			<< std::endl << std::endl;

	EVector Residual;
	CalculateResidual(Residual);

	PrintEigenvalues(Eigenvalues, Residual, std::min(num_print,num_converged), num_converged);
}

template<typename T>
void MatrixInterface<T>::CalculateResidual(EVector &Residual) const
{
	if(!compute_vecs)
	{
		printf("Cannot calculate the residules without the eigenvectors\n");
		return;
	}

	const int Nvar = A.rows();
	Residual = EVector::Zero(Nvar);
	for(int i=0; i<num_converged; ++i)
	{
		Eigen::Map<const EVector> Evec(&Eigenvectors[Nvar*i], Nvar);
		T eMe = MultvMw(Evec, Evec);
		Residual[i] = eMe - Eigenvalues[i];
	}
}


template<typename T>
T MatrixInterface<T>::GetVal(int m) const
{
	return Eigenvalues[m];
}

template<typename T>
T MatrixInterface<T>::GetVec(int m, int i) const
{
	return Eigenvectors[A.rows()*m+i];
}

template<typename T>
void MatrixInterface<T>::LUdecomp()
{
#ifndef DONT_USE_SUITESPARSE
	UMFsolver.compute(A);
	if(UMFsolver.info()!=Eigen::Success) 
	{
		printf("ERROR: Eigen::UmfPackLU decomposition failed\n");	
		exit(EXIT_FAILURE);
	}
#else
	printf("ERROR: Umfpack is turned off.\n"); 
	fflush(stdout);
	assert(false);
#endif //DONT_USE_SUITESPARSE
}

template<typename T>
void MatrixInterface<T>::solve_Mx_equals_b(EVector &X, EVector const &B)
{
#ifndef DONT_USE_SUITESPARSE
	X = UMFsolver.solve(B);
#else
	printf("ERROR: Umfpack is turned off.\n"); 
	fflush(stdout);
	assert(false);
#endif //DONT_USE_SUITESPARSE
}

template<typename T>
void MatrixInterface<T>::solve_Mx_equals_b(T *x, T const *b)
{
#ifndef DONT_USE_SUITESPARSE
	Eigen::Map<      EVector> X(&x[0], A.rows());
	Eigen::Map<const EVector> B(&b[0], A.rows());
	X = UMFsolver.solve(B);
	//solve_Mx_equals_b(X,B);
#else
	printf("ERROR: Umfpack is turned off.\n"); 
	fflush(stdout);
	assert(false);
#endif //DONT_USE_SUITESPARSE
}


}

#endif //MATRIX_INTERFACE_H





