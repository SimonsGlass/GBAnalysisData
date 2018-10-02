#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "../Resources/std_include.h"
#include "../Computers/BaseComputer.h"
#include "minimizerFIRE.h"

namespace LiuJamming
{

//! A Class to clean up the Minimization interface for standard minimization operations.
template<int Dim>
class CSimpleMinimizer
{
public:
	//!Minimization algorithm flags.
	enum{
		NONE=0, //!<No minimization algorithm
		FAST,	//!<Fast combination of algorithms
		FIRE,	//!<FIRE algorithm
		LBFGS,	//!<L-BFGS algorithm
		CG_PR,	//!<Conjugate gradient algorithm- PR variant
		CG_FR	//!<Conjugate gradient algorithm- FR variant
	};
private:
	CBaseComputer<Dim> *pComputer;
	int N_dof;

public:
	//!Constructor
	CSimpleMinimizer(CBaseComputer<Dim> &TargetComputer, int minimize_type=NONE, dbl tol=1e-12, int max_iterations=-1, int print_iter=1000);
	//!Simple minimization with the FIRE algorithm
	void minimizeFIRE(dbl tol=1e-12, int max_iterations = -1, dbl delta_t_start=-1., int print_iter=1000);
};

/**
 *	To separate the declaration from the minimization (which increases flexibility), call the constructor with only the 
 *	first parameter. If the second parameter is set to one of the enumerated algorithm flags, then automatic minimization
 *	will occur, with the given values of tol and max_iterations.
 *
 *	@param[in] TargetComputer Computer object on which to perform minimization.
 *	@param[in] minimize_type determines algorithm used for automatic minimization. default = NONE.
 *	@param[in] tol Minimization tolerance for automatic minimization. default = 1e-12.
 *	@param[in] max_iterations Maximum number of iterations. Unlimited if less than or equal to zero. default = -1.
 */
template<int Dim>
CSimpleMinimizer<Dim>::CSimpleMinimizer(CBaseComputer<Dim> &TargetComputer, int minimize_type, dbl tol, int max_iterations, int print_iter)
	: pComputer(&TargetComputer)
{
	N_dof = pComputer->GetNdof();
	switch(minimize_type)
	{
		case NONE:
			break;
		case FAST:
		case FIRE:
			minimizeFIRE(tol, max_iterations,-1,print_iter);
			break;
	}
};

/**
 *	If automatic minimization is not used, this can be called to minimize with the FIRE algorithm.
 *
 *	@param[in] tol Minimization tolerance for automatic minimization. default = 1e-12.
 *	@param[in] max_iterations Maximum number of iterations. Unlimited if less than or equal to zero. default = -1.
 *	@param[in] delta_t_start initial MD time step. Set automatically from the Computer if less than or equal to zero. default = -1.
 */
template<int Dim>
void CSimpleMinimizer<Dim>::minimizeFIRE(dbl tol, int max_iterations, dbl delta_t_start, int print_iter)
{
	typedef CMinimizerFIRE< CBaseComputer<Dim> > MINIMIZER;
	MINIMIZER mm(N_dof);

	mm.set_functions(pComputer, &CBaseComputer<Dim>::Evaluate, &CBaseComputer<Dim>::Progress, &CBaseComputer<Dim>::Move, &CBaseComputer<Dim>::ReportHeader);

	if(max_iterations > 0)
		mm.set_max_iterations(max_iterations);

//	dbl delta_t_start = 0.04*plist.get_length_unit();
	if(delta_t_start <= 0.)
	{
		delta_t_start = 0.04*pComputer->GetMinimizationTimeScale();
	}
	mm.set_FIRE_delta_t_start(delta_t_start);
	
	mm.set_FIRE_delta_t_max(10.*delta_t_start);
	if(max_iterations>0) mm.set_FIRE_delta_t_max(100000.*delta_t_start);

	mm.set_print_iter(print_iter);

	mm.minimize(tol);
};


}

#endif //MINIMIZER
