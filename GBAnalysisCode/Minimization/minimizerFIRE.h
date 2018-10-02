#ifndef MINIMIZERFIRE_H
#define MINIMIZERFIRE_H

#include "../Resources/std_include.h"

namespace LiuJamming
{


template<class T>
class CMinimizerFIRE
{
public:
	enum //Minimization status codes
	{
		MINIMIZATION_SUCCESS = 1,				
		MINIMIZATION_FAILURE = -999,
		MAX_ITERATIONS_REACHED,
		OUT_OF_MEMORY,
		INVALID_PARAMETERS,
		INVALID_INITIALIZATION,
		UNKNOWN_ERROR
	};

	typedef void (T::*evaluate_t)		(Eigen::VectorXd &grad, dbl &fx);
	typedef bool (T::*progress_t)		(Eigen::VectorXd const &grad, dbl fx, int iteration, int print_iter, dbl tol);
	typedef void (T::*move_t)			(Eigen::VectorXd const &step);
	typedef void (T::*report_header_t)	();

private:
	int Nvar;		//number of variables to be minimized.

	//Parameters
	int  PARAM_max_iterations;		//Optimization process terminates after PARAM_max_iterations interations.
	bool PARAM_no_max_iterations; 		//If true, contiues the optimization process until convergence or an error (essentially sets PARAM_MAX_iterations = infinity).
	dbl  PARAM_convergence_tol;		//This is a parameter that gets passed to the is_converged function so the user can change their desired tolerance easily.
	int  FIRE_PARAM_N_min;
	dbl  FIRE_PARAM_f_inc;
	dbl  FIRE_PARAM_f_dec;
	dbl  FIRE_PARAM_alpha_start;
	dbl  FIRE_PARAM_f_alpha;
	dbl  FIRE_PARAM_delta_t_start;
	dbl  FIRE_PARAM_delta_t_max;
	int  IPARAM_print_iter;
	bool IPARAM_do_not_use_report_header;	//This is set true if the user does not pass a report_header_t function.
	bool IPARAM_functions_set;


	//Object and pointers object member functions
	T *obj;
	evaluate_t	 	evaluate_fn;
	progress_t		progress_fn;
	move_t			move_fn;
	report_header_t	report_header_fn;

	void evaluate(Eigen::VectorXd &grad, dbl &fx)
	{
		((*obj).*(evaluate_fn))(grad, fx);
	};
	bool progress(Eigen::VectorXd const &grad, dbl fx, int iteration, dbl tol)
	{
		return ((*obj).*(progress_fn))(grad, fx, iteration, IPARAM_print_iter, tol);
	};
	void move(Eigen::VectorXd const &step)
	{
		((*obj).*(move_fn))(step);
	};
	void report_header()
	{
		if(IPARAM_print_iter>0)
			((*obj).*(report_header_fn))();
	};

public:
	CMinimizerFIRE(int N) 
		: Nvar(N)
	{
		set_default_params();
	};


	//Set parameters
	void set_max_iterations(int imax)					{ PARAM_max_iterations					= imax;	
														  PARAM_no_max_iterations				= false;			};
	void set_print_iter(int _print_iter)				{ IPARAM_print_iter						= _print_iter;		};
	void set_no_max_iterations(bool b=true)				{ PARAM_no_max_iterations				= b;				};
	void set_tol(dbl t)									{ PARAM_convergence_tol					= t;				};
	void set_FIRE_delta_t_start(dbl t)					{ FIRE_PARAM_delta_t_start				= t;				};
	void set_FIRE_delta_t_max(dbl t)					{ FIRE_PARAM_delta_t_max				= t;				};
	void set_default_params()
	{
		PARAM_max_iterations = 1e8;
		PARAM_no_max_iterations = true;
		PARAM_convergence_tol = 1e-12;

		FIRE_PARAM_N_min			= 5;
		FIRE_PARAM_f_inc			= 1.1;
		FIRE_PARAM_f_dec			= 0.5;
		FIRE_PARAM_alpha_start		= 0.1;
		FIRE_PARAM_f_alpha			= 0.99;
		FIRE_PARAM_delta_t_start	= 1e-8;
		FIRE_PARAM_delta_t_max		= 1000000.*FIRE_PARAM_delta_t_start;
		
		IPARAM_print_iter			= 5000;
		IPARAM_do_not_use_report_header = false;
		IPARAM_functions_set = false;
	};


	//Function called by the user to set the RCI functions. 
	void set_functions(T *new_obj, evaluate_t new_fn1, progress_t new_fn2, move_t new_fn3=NULL, report_header_t new_fn4=NULL)
	{
		obj = new_obj;
		evaluate_fn			= new_fn1;
		progress_fn			= new_fn2;
		move_fn				= new_fn3;
		report_header_fn	= new_fn4;
		if(new_fn4 == NULL)	IPARAM_do_not_use_report_header = true;
		IPARAM_functions_set = true;
	};
	
	//step is just a dummy array
	void velocity_verlet_integrate(Eigen::VectorXd &v, Eigen::VectorXd &grad, Eigen::VectorXd &step, dbl &fx, dbl delta_t)
	{
		int i;
		dbl half_delta_t  = 0.5*delta_t;
		dbl half_delta_t2 = 0.5*delta_t*delta_t;
		
		//calculate x(t+delta_t) = x(t) + v(t)*delta_t + 0.5*a(t)*delta_t^2
		step = delta_t*v - half_delta_t2*grad; //a = f = -g
		move(step);
		
		//calculate v(t+delta_t/2) = v(t) + 0.5*a(t)*delta_t
		v -= half_delta_t*grad;
		
		//calculate the energy and the gradient
		evaluate(grad, fx);

		//calculate v(t+delta_t) = v(t+delta_t/2) + 0.5*a(t+delta_t)*delta_t
		v -= half_delta_t*grad;
	};   

	int minimize(dbl tol)
	{
		set_tol(tol);
		return minimize();
	};

#define min2(a, b)      ((a) <= (b) ? (a) : (b))
	int minimize()
	{
		printf("begin FIRE minimization with tol = %e\n", PARAM_convergence_tol);
		int ii;	//Iteration step
		int i;	//iterator for various for loops.
		int ret;

		Eigen::VectorXd v, grad, step;
		v = grad = step = Eigen::VectorXd::Zero(Nvar);

		//initial conditions; 
		dbl delta_t = FIRE_PARAM_delta_t_start;
		dbl alpha = FIRE_PARAM_alpha_start;

		//Other variables
		dbl last_fx = 1e12;
		dbl fx=0., P;
		dbl coeff1, coeff2;
		int NsinceNegP = 0;                     //Number of iterations since P was last not positive

		if(!IPARAM_do_not_use_report_header)
			report_header();

		//Algorithm loop
		ii = 0;
		while(true)
		{
			//Run MD integration
			velocity_verlet_integrate(v, grad, step, fx, delta_t); //Does not move positions on the first call.
			//Convergence test and report
			if(progress(grad, fx, ii, PARAM_convergence_tol))
				{ ret = MINIMIZATION_SUCCESS; break; }
			//If exceeded the maximum number of iterations
			if (PARAM_max_iterations != 0 && PARAM_max_iterations < ii+1)
				{ ret = MAX_ITERATIONS_REACHED; break; }
			//Calculate power
			P = -grad.dot(v);

			//set v = (1-alpha)v + alpha|v|fhat.
			coeff1 = 1.-alpha;
			coeff2 = alpha*sqrt(v.squaredNorm()/grad.squaredNorm());
			v = coeff1*v - coeff2*grad;
			
			//If P>0 and Nsinceneg>Nmin, increase delta_t and decrease alpha
			if(P>0 && NsinceNegP > FIRE_PARAM_N_min)
			{
				delta_t = min2(delta_t*FIRE_PARAM_f_inc,FIRE_PARAM_delta_t_max);
				alpha *= FIRE_PARAM_f_alpha;
			}
		

			//If P<=0, decrease delta_t and set v=0 and set alpha = alpha_start
			if(P <= 0.)
			{
				delta_t *= FIRE_PARAM_f_dec;
				v.setZero();
				alpha = FIRE_PARAM_alpha_start;
				NsinceNegP = 0;
			}       
			else
			{
				++NsinceNegP;
			}
			++ii;

			/*
			//If target function (fx) increased decrease delta_t, set v=0, and alpha = alpha_start
			if(fx>1.0001*last_fx)
			{
				delta_t *= FIRE_PARAM_f_dec;
				v.setZero();
				alpha=FIRE_PARAM_alpha_start;
				//NsinceNegP = 0; //keep?
			};
			*/
			last_fx=fx;


		}

		return ret;
	};

};

}

#endif //MINIMIZERFIRE
