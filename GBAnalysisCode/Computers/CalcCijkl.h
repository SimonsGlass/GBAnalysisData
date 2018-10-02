#include "cijkl.h"

namespace LiuJamming
{

template <int Dim, class OBJ_CLASS, class CONTAINER_CLASS>
dbl CalculateGenericElasticResponse(
{	
	//Calculate the displacement vector for each bond in the metric defined by the strain tensor.
	//Also, calculate the forces on each node due to the change in metric.
	obj.Calculate_n_d2Udvdgamma(strain_tensor, n_d2Udvdgamma, DeltaR_bond); 

	//Solve for the non-affine displacement
	D1.solve_Mx_equals_b(uNonAffine_node, n_d2Udvdgamma);

	dbl dE = obj.CalculateEnergyChange(DeltaR_bond, uNonAffine_node);

}

template <int Dim, class OBJ_CLASS, class CONTAINER_CLASS>
void CalculateGenericCijkl(OBJ_CLASS &obj, MatrixInterface<dbl> &D1, cCIJKL<Dim> &cijkl, dbl Volume, int Verbose)
{
	//First, do an LU decomposition on the matrix.
	D1.LUdecomp();

	int Nvar = D1.A.rows();
	dbl prefactor = 2./Volume;//The prefactor contains a 2 in it from the equation dE/V = (1/2)cijkl uij ukl.

	Eigen::VectorXd n_d2Udvdgamma		= Eigen::VectorXd::Zero(Nvar);
	Eigen::VectorXd uNonAffine_node		= Eigen::VectorXd::Zero(Nvar);
	CONTAINER_CLASS DeltaR_bond;
//	Eigen::VectorXd DeltaR_bond			= Eigen::VectorXd::Zero(Dim*(int)list.size());
	Eigen::Matrix<dbl, Dim, Dim> strain_tensor;

	for(int ii=0; ii<cCIJKL<Dim>::num_constants; ++ii)
	{   
		//Set the strain tensor
		cijkl.set_strain_tensor(strain_tensor, ii);

		//Calculate the displacement vector for each bond in the metric defined by the strain tensor.
		//Also, calculate the forces on each node due to the change in metric.
		obj.Calculate_n_d2Udvdgamma(strain_tensor, n_d2Udvdgamma, DeltaR_bond); 

		//Solve for the non-affine displacement
		D1.solve_Mx_equals_b(uNonAffine_node, n_d2Udvdgamma);

		dbl dE = obj.CalculateEnergyChange(DeltaR_bond, uNonAffine_node);

//		//Add the non-affine extension of each bond to the affine extension
//		for(int bi=0; bi<(int)list.size(); ++bi)
//			DeltaR_bond.segment<Dim>(Dim*bi) += uNonAffine_node.segment<Dim>(Dim*list[bi].j) - uNonAffine_node.segment<Dim>(Dim*list[bi].i);
//
//		//Calculate the change in energy
//		dbl dE = CalculateEnergyChange(DeltaR_bond, unstress_coeff, stress_coeff);

		dE *= prefactor; //dE is now 2*dE/V. This comes from the equation dE/V = (1/2)cijkl uij ukl.
		
		//Set the ii'th elastic constant
		cijkl.set_constant(dE, ii);
	}
	if(Verbose)
		cijkl.print();

}












}



