#ifndef MINIMIZERLBFGS_H
#define MINIMIZERLBFGS_H

#include "../Resources/std_include.h"

template<class T>
class CMinimizerLBFGS
{

public:
	static void report_status_code(int status_code);

};























template<class T>
static void CMinimizerLBFGS<T>::report_status_code(int status_code)
{
	printf("Status report: \t\t");
	switch(status_code)
	{   
		case LBFGS_SUCCESS:                         printf("L-BFGS reaches convergence.\n"); break;
		case LBFGS_STOP:                            printf("L-BFGS was stopped.\n"); break;
		case LBFGS_ALREADY_MINIMIZED:               printf("The initial variables already minimize the objective function.\n"); break;
		case LBFGSERR_UNKNOWNERROR:                 printf("Unknown error.\n"); break;
		case LBFGSERR_LOGICERROR:                   printf("Logic error.\n"); break;
		case LBFGSERR_OUTOFMEMORY:                  printf("Insufficient memory.\n"); break;
		case LBFGSERR_CANCELED:                     printf("The minimization process has been canceled.\n"); break;
		case LBFGSERR_INVALID_N:                    printf("Invalid number of variables specified.\n"); break;
		case LBFGSERR_INVALID_N_SSE:                printf("Invalid number of variables (for SSE) specified.\n"); break;
		case LBFGSERR_INVALID_X_SSE:                printf("The array x must be aligned to 16 (for SSE).\n"); break;
		case LBFGSERR_INVALID_EPSILON:              printf("Invalid parameter lbfgs_parameter_t::epsilon specified.\n"); break;
		case LBFGSERR_INVALID_TESTPERIOD:           printf("Invalid parameter lbfgs_parameter_t::past specified.\n"); break;
		case LBFGSERR_INVALID_DELTA:                printf("Invalid parameter lbfgs_parameter_t::delta specified.\n"); break;
		case LBFGSERR_INVALID_LINESEARCH:           printf("Invalid parameter lbfgs_parameter_t::linesearch specified.\n"); break;
		case LBFGSERR_INVALID_MINSTEP:              printf("Invalid parameter lbfgs_parameter_t::max_step specified.\n"); break;
		case LBFGSERR_INVALID_MAXSTEP:              printf("Invalid parameter lbfgs_parameter_t::max_step specified.\n"); break;
		case LBFGSERR_INVALID_FTOL:                 printf("Invalid parameter lbfgs_parameter_t::ftol specified.\n"); break;
		case LBFGSERR_INVALID_WOLFE:                printf("Invalid parameter lbfgs_parameter_t::wolfe specified.\n"); break;
		case LBFGSERR_INVALID_GTOL:                 printf("Invalid parameter lbfgs_parameter_t::gtol specified.\n"); break;
		case LBFGSERR_INVALID_XTOL:                 printf("Invalid parameter lbfgs_parameter_t::xtol specified.\n"); break;
		case LBFGSERR_INVALID_MAXLINESEARCH:        printf("Invalid parameter lbfgs_parameter_t::max_linesearch specified.\n"); break;
		case LBFGSERR_INVALID_ORTHANTWISE:          printf("Invalid parameter lbfgs_parameter_t::orthantwise_c specified.\n"); break;
		case LBFGSERR_INVALID_ORTHANTWISE_START:    printf("Invalid parameter lbfgs_parameter_t::orthantwise_start specified.\n"); break;
		case LBFGSERR_INVALID_ORTHANTWISE_END:      printf("Invalid parameter lbfgs_parameter_t::orthantwise_end specified.\n"); break;
		case LBFGSERR_OUTOFINTERVAL:                printf("The line-search step went out of the interval of uncertainty.\n"); break;
		case LBFGSERR_INCORRECT_TMINMAX:            printf("A logic error occurred; alternatively, the interval of uncertainty became too small.\n"); break;
		case LBFGSERR_ROUNDING_ERROR:               printf("A rounding error occurred; alternatively, no line-search step satisfies the sufficient decrease and curvature conditions.\n"); break;
		case LBFGSERR_MINIMUMSTEP:                  printf("The line-search step became smaller than lbfgs_parameter_t::min_step.\n"); break;
		case LBFGSERR_MAXIMUMSTEP:                  printf("The line-search step became larger than lbfgs_parameter_t::max_step.\n"); break;
		case LBFGSERR_MAXIMUMLINESEARCH:            printf("The line-search routine reaches the maximum number of evaluations.\n"); break;
		case LBFGSERR_MAXIMUMITERATION:             printf("The algorithm routine reaches the maximum number of iterations.\n"); break;
		case LBFGSERR_WIDTHTOOSMALL:                printf("Relative width of the interval of uncertainty is at most lbfgs_parameter_t::xtol.\n"); break;
		case LBFGSERR_INVALIDPARAMETERS:            printf("A logic error (negative line-search step) occurred.\n"); break;
		case LBFGSERR_INCREASEGRADIENT:             printf("The current search direction increases the objective function value.\n"); break;
		default:                                    printf("Unknown status code.\n"); break;
	}   
};





#endif //MINIMIZERLBFGS
