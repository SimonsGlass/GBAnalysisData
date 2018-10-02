#ifndef COMPUTE_SYMMETRY_FUNCTIONS

#define COMPUTE_SYMMETRY_FUNCTIONS

#include <list>
#include <vector>
#include "../Resources/std_include.h"
#include "../State/StaticState.h"
//#include "../Resources/Settings.h"

namespace LiuJamming {

template <int Dim>
void ComputeSymmetryFunctions(CStaticState<Dim> &State,CSimpleGrid<Dim> *Grid, Settings &settings, list<int> Particles, list<vector<double> > &Symmetry)
{
	//Load the symmetry function
	vector<double> RadialEta = settings["Radial_Eta"];
	vector<double> AngularEta = settings["Angular_Eta"];
	vector<double> Zeta = settings["Zeta"];
	vector<double> Lambda = settings["Lambda"]; 

	double Cutoff = settings["Cutoff"];

	int Number_A = settings["Number_A"];

	Eigen::Matrix<double,Dim,1> Displacement_JK;

	for(list<int>::iterator it = Particles.begin() ; it!=Particles.end() ; it++)
	{
	//	cout << "Computing symmetry functions for particle " << *it << endl;
		
		vector<double> Symmetry_Functions;

		//compute radial functions
		for(vector<double>::iterator eta_it = RadialEta.begin() ; eta_it!= RadialEta.end() ; eta_it++)
		{	
			double G1XA = 0.0;
			double G1XB = 0.0;

			const list<CSimpleNeighbor<Dim> > &neighbors = Grid->GetNeighbors(*it);

			for(typename list<CSimpleNeighbor<Dim> >::const_iterator neigh_it = neighbors.begin() ; neigh_it != neighbors.end() ; neigh_it++)
			{
				if((*neigh_it).j < Number_A)
				{
					double norm = (*neigh_it).Displacement.norm();
					G1XA += 0.5*exp(-(*eta_it)*norm*norm)*(cos(PI*norm/Cutoff)+1);
				}else{
					double norm = (*neigh_it).Displacement.norm();
					G1XB += 0.5*exp(-(*eta_it)*norm*norm)*(cos(PI*norm/Cutoff)+1);
				}
			}
			Symmetry_Functions.push_back(G1XA);
			Symmetry_Functions.push_back(G1XB);
		}

		//compute angular functions
		for(int a_it = 0 ; a_it < AngularEta.size() ; a_it++)
		{
	//		cout << "Computing angular function " << a_it << endl;
			double G2XAA = 0.0;
			double G2XAB = 0.0;
			double G2XBB = 0.0;

			double eta = AngularEta[a_it];
			double zeta = Zeta[a_it];
			double lambda = Lambda[a_it];

			list<CSimpleNeighbor<Dim> > neighbors = Grid->GetNeighbors(*it);
			
			for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_J = neighbors.begin() ; neighbor_J != neighbors.end() ; neighbor_J++)
			{
				for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_K = neighbors.begin() ; neighbor_K != neighbors.end() ; neighbor_K++)
             		        {
					if((*neighbor_K).j != (*neighbor_J).j)
					{
						State.GetDisplacement((*neighbor_J).j,(*neighbor_K).j,Displacement_JK);
									
						double norm_JK = Displacement_JK.squaredNorm();
						if(norm_JK < Cutoff*Cutoff)
						{
							norm_JK = sqrt(norm_JK);
							double norm_IJ = (*neighbor_J).Displacement.norm();
                                                	double norm_IK = (*neighbor_K).Displacement.norm();
	
							double ctheta = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK;											
							if((*neighbor_J).j < Number_A)
							{
								/*
								if((*neighbor_K).j < Number_A)
									G2XAA += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								else
									G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								*/
								
								if((*neighbor_K).j< Number_A)
                                                                        G2XAA += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/Cutoff)+1)*(cos(PI*norm_IK/Cutoff)+1)*(cos(PI*norm_JK/Cutoff)+1);
                                                                else
                                                                        G2XAB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/Cutoff)+1)*(cos(PI*norm_IK/Cutoff)+1)*(cos(PI*norm_JK/Cutoff)+1);
							}else{
								if((*neighbor_K).j > Number_A)
									G2XBB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/Cutoff)+1)*(cos(PI*norm_IK/Cutoff)+1)*(cos(PI*norm_JK/Cutoff)+1);
								else
									G2XAB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/Cutoff)+1)*(cos(PI*norm_IK/Cutoff)+1)*(cos(PI*norm_JK/Cutoff)+1);
								/*
								if((*neighbor_K).j < Number_A)
                                                                        G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
                                                                else
                                                      	                G2XBB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								*/
							}
						}							
					}
				}
			}
			
			Symmetry_Functions.push_back(G2XAA);
			Symmetry_Functions.push_back(G2XAB);
			Symmetry_Functions.push_back(G2XBB);
		}	
		
		Symmetry.push_back(Symmetry_Functions);
	}
}

template <int Dim>
void ComputeSymmetryFunctionsGradient(CStaticState<Dim> &State,CSimpleGrid<Dim> *Grid, Settings &settings, list<int> Particles, list<vector<Eigen::Matrix<double,Dim,1> > > &Symmetry)
{
	//Load the symmetry function
	vector<double> RadialEta = settings["Radial_Eta"];
	vector<double> AngularEta = settings["Angular_Eta"];
	vector<double> Zeta = settings["Zeta"];
	vector<double> Lambda = settings["Lambda"]; 

	double Cutoff = settings["Cutoff"];

	int Number_A = settings["Number_A"];

	Eigen::Matrix<double,Dim,1> Displacement_JK;

	for(list<int>::iterator it = Particles.begin() ; it!=Particles.end() ; it++)
	{
	//	cout << "Computing symmetry functions for particle " << *it << endl;
		
		vector<Eigen::Matrix<double,Dim,1> > Symmetry_Functions;

		//compute radial functions
		for(vector<double>::iterator eta_it = RadialEta.begin() ; eta_it!= RadialEta.end() ; eta_it++)
		{	
			Eigen::Matrix<double,Dim,1> G1XA = Eigen::Array<double,Dim,1>::Zero();
			Eigen::Matrix<double,Dim,1> G1XB = Eigen::Array<double,Dim,1>::Zero();

			const list<CSimpleNeighbor<Dim> > &neighbors = Grid->GetNeighbors(*it);

			for(typename list<CSimpleNeighbor<Dim> >::const_iterator neigh_it = neighbors.begin() ; neigh_it != neighbors.end() ; neigh_it++)
			{
				double norm = (*neigh_it).Displacement.norm();
				Eigen::Matrix<double,Dim,1> gradient = -2.0*( (*neigh_it).Displacement * (*eta_it) * (cos(PI*norm/Cutoff) + 1.0) + PI/2.0/Cutoff/norm*(*neigh_it).Displacement * sin(PI*norm/Cutoff))*exp(-(*eta_it)*norm*norm);
			
				if((*neigh_it).j < Number_A)
				{
					G1XA += gradient; 
				}else{
					G1XB += gradient; 
				}
			}
			Symmetry_Functions.push_back(G1XA);
			Symmetry_Functions.push_back(G1XB);
		}

		//compute angular functions
		for(int a_it = 0 ; a_it < AngularEta.size() ; a_it++)
		{
	//		cout << "Computing angular function " << a_it << endl;
			Eigen::Matrix<double,Dim,1> G2XAA = Eigen::Array<double,Dim,1>::Zero(); 
			Eigen::Matrix<double,Dim,1> G2XAB = Eigen::Array<double,Dim,1>::Zero();
			Eigen::Matrix<double,Dim,1> G2XBB = Eigen::Array<double,Dim,1>::Zero();

			double eta = AngularEta[a_it];
			double zeta = Zeta[a_it];
			double lambda = Lambda[a_it];

			list<CSimpleNeighbor<Dim> > neighbors = Grid->GetNeighbors(*it);
			
			for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_J = neighbors.begin() ; neighbor_J != neighbors.end() ; neighbor_J++)
			{
				for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_K = neighbors.begin() ; neighbor_K != neighbors.end() ; neighbor_K++)
             		        {
					if((*neighbor_K).j != (*neighbor_J).j)
					{
						State.GetDisplacement((*neighbor_J).j,(*neighbor_K).j,Displacement_JK);
									
						double norm_JK = Displacement_JK.squaredNorm();
						if(norm_JK < Cutoff*Cutoff)
						{
							norm_JK = sqrt(norm_JK);
							double norm_IJ = (*neighbor_J).Displacement.norm();
                                                	double norm_IK = (*neighbor_K).Displacement.norm();
	
							double ctheta = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK;	
							double prefactor = pow(2.0,1.0-zeta);							
					
							double costerm = pow(1+lambda*ctheta,zeta);
							double cutoff_IJ = 0.5*(cos(norm_IJ*PI/Cutoff) +1.0);
							double cutoff_IK = 0.5*(cos(norm_IK*PI/Cutoff) +1.0);
							double cutoff_JK = 0.5*(cos(norm_JK*PI/Cutoff) +1.0);
							

							Eigen::Matrix<double,Dim,1> gradient = -2.0*eta*((*neighbor_J).Displacement + (*neighbor_K).Displacement)*costerm*cutoff_IJ*cutoff_IK*cutoff_JK + 1.0/norm_IJ/norm_IK*((*neighbor_J).Displacement + (*neighbor_K).Displacement - (*neighbor_J).Displacement.dot((*neighbor_K).Displacement) * ((*neighbor_J).Displacement/norm_IJ/norm_IJ + (*neighbor_K).Displacement/norm_IK/norm_IK))*cutoff_IJ*cutoff_IK*cutoff_JK - PI/2.0/Cutoff*((*neighbor_J).Displacement/norm_IJ*sin(PI*norm_IJ/Cutoff)*cutoff_JK*cutoff_IK + (*neighbor_K).Displacement/norm_IK*sin(PI*norm_IK/Cutoff)*cutoff_IJ*cutoff_JK)*costerm;
				
							gradient *= prefactor*exp(-eta*(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK*norm_JK));

							if((*neighbor_J).j < Number_A)
							{		
								if((*neighbor_K).j< Number_A)
                                                                        G2XAA += gradient;
                                                                else
                                                                        G2XAB += gradient;
							}else{
								if((*neighbor_K).j > Number_A)
									G2XBB += gradient;
								else
									G2XAB += gradient;
							}
						}							
					}
				}
				
				list<CSimpleNeighbor<Dim> > neighbors_J = Grid->GetNeighbors((*neighbor_J).j); 
				for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_K = neighbors_J.begin() ; neighbor_K != neighbors_J.end() ; neighbor_K++)
             		        {
					if((*neighbor_K).j != (*neighbor_J).i)
					{
						State.GetDisplacement((*neighbor_J).i,(*neighbor_K).j,Displacement_JK);
									
						double norm_JK = Displacement_JK.squaredNorm();
						if(norm_JK < Cutoff*Cutoff)
						{
							norm_JK = sqrt(norm_JK);
							double norm_IJ = (*neighbor_J).Displacement.norm();
                                                	double norm_IK = (*neighbor_K).Displacement.norm();
	
							double ctheta = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK;	
							double prefactor = 2.0*pow(2.0,1.0-zeta);							
					
							double costerm = pow(1+lambda*ctheta,zeta);
							double cutoff_IJ = 0.5*(cos(norm_IJ*PI/Cutoff) +1.0);
							double cutoff_IK = 0.5*(cos(norm_IK*PI/Cutoff) +1.0);
							double cutoff_JK = 0.5*(cos(norm_JK*PI/Cutoff) +1.0);
							

							Eigen::Matrix<double,Dim,1> gradient = -2.0*eta*(*neighbor_J).Displacement*costerm*cutoff_IJ*cutoff_IK*cutoff_JK + 1.0/norm_IJ/norm_IK*((*neighbor_K).Displacement - (*neighbor_J).Displacement.dot((*neighbor_K).Displacement) * (*neighbor_J).Displacement/norm_IJ/norm_IJ)*cutoff_IJ*cutoff_IK*cutoff_JK - PI/2.0/Cutoff*(*neighbor_J).Displacement/norm_IJ*sin(PI*norm_IJ/Cutoff)*cutoff_JK*cutoff_IK*costerm;
				
							gradient *= prefactor*exp(-eta*(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK*norm_JK));

							if((*neighbor_J).i < Number_A)
							{		
								if((*neighbor_K).j< Number_A)
                                                                        G2XAA += gradient;
                                                                else
                                                                        G2XAB += gradient;
							}else{
								if((*neighbor_K).j > Number_A)
									G2XBB += gradient;
								else
									G2XAB += gradient;
							}
						}							
					}
				}
			}
			
			Symmetry_Functions.push_back(G2XAA);
			Symmetry_Functions.push_back(G2XAB);
			Symmetry_Functions.push_back(G2XBB);
		}	
		
		Symmetry.push_back(Symmetry_Functions);
	}
}


template <int Dim>
void ComputeSymmetryFunctions(CStaticState<Dim> &State,list<list<CSimpleNeighbor<Dim> > > &NeighborList, Settings &settings, list<int> Particles, list<vector<double> > &Symmetry)
{
	//Load the symmetry function
	vector<double> RadialEta = settings["Radial_Eta"];
	vector<double> AngularEta = settings["Angular_Eta"];
	vector<double> Zeta = settings["Zeta"];
	vector<double> Lambda = settings["Lambda"]; 

	double Cutoff = settings["Cutoff"];

	int Number_A = settings["Number_A"];

	Eigen::Matrix<double,Dim,1> Displacement_JK;

	typename list<list<CSimpleNeighbor<Dim> > >::iterator neighbor_it = NeighborList.begin();

	for(list<int>::iterator it = Particles.begin() ; it!=Particles.end() ; it++)
	{
		cout << "Computing symmetry functions for particle " << *it << endl;
		
		vector<double> Symmetry_Functions;

		//compute radial functions
		for(vector<double>::iterator eta_it = RadialEta.begin() ; eta_it!= RadialEta.end() ; eta_it++)
		{	
			double G1XA = 0.0;
			double G1XB = 0.0;

			list<CSimpleNeighbor<Dim> > &neighbors = *neighbor_it;

			for(typename list<CSimpleNeighbor<Dim> >::iterator neigh_it = neighbors.begin() ; neigh_it != neighbors.end() ; neigh_it++)
			{
				if((*neigh_it).j < Number_A)
				{
					double norm = (*neigh_it).Displacement.norm();
					G1XA += 0.5*exp(-(*eta_it)*norm*norm)*(cos(PI*norm/Cutoff)+1);
				}else{
					double norm = (*neigh_it).Displacement.norm();
					G1XB += 0.5*exp(-(*eta_it)*norm*norm)*(cos(PI*norm/Cutoff)+1);
				}
			}
			Symmetry_Functions.push_back(G1XA);
			Symmetry_Functions.push_back(G1XB);
		}

		//compute angular functions
		for(int a_it = 0 ; a_it < AngularEta.size() ; a_it++)
		{
			cout << "Computing angular function " << a_it << endl;
			double G2XAA = 0.0;
			double G2XAB = 0.0;
			double G2XBB = 0.0;

			double eta = AngularEta[a_it];
			double zeta = Zeta[a_it];
			double lambda = Lambda[a_it];

			list<CSimpleNeighbor<Dim> > &neighbors = *neighbor_it;	
		
			for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_J = neighbors.begin() ; neighbor_J != neighbors.end() ; neighbor_J++)
			{
				for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_K = neighbors.begin() ; neighbor_K != neighbors.end() ; neighbor_K++)
             		        {
					if((*neighbor_K).j != (*neighbor_J).j)
					{
						State.GetDisplacement((*neighbor_J).j,(*neighbor_K).j,Displacement_JK);
									
						double norm_IJ = (*neighbor_J).Displacement.norm();
						double norm_IK = (*neighbor_K).Displacement.norm();
						double norm_JK = Displacement_JK.norm();
						if(norm_JK < Cutoff && norm_IJ < Cutoff && norm_IK < Cutoff)
						{
							double ctheta = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK;											
							if((*neighbor_J).j < Number_A)
							{
								/*
								if((*neighbor_K).j < Number_A)
									G2XAA += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								else
									G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								*/
								if((*neighbor_K).j< Number_A)
                                                                        G2XAA += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/Cutoff)+1)*(cos(PI*norm_IK/Cutoff)+1)*(cos(PI*norm_JK/Cutoff)+1);
                                                                else
                                                                        G2XAB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/Cutoff)+1)*(cos(PI*norm_IK/Cutoff)+1)*(cos(PI*norm_JK/Cutoff)+1);
                                                        }else{
                                                                if((*neighbor_K).j > Number_A)
                                                                        G2XBB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/Cutoff)+1)*(cos(PI*norm_IK/Cutoff)+1)*(cos(PI*norm_JK/Cutoff)+1);
								else
									G2XAB += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/Cutoff)+1)*(cos(PI*norm_IK/Cutoff)+1)*(cos(PI*norm_JK/Cutoff)+1);/*
							}else{
								if((*neighbor_K).j < Number_A)
                                                                        G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
                                                                else
                                                      	                G2XBB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
							*/
							}
						}							
					}
				}
			}
			
			Symmetry_Functions.push_back(G2XAA);
			Symmetry_Functions.push_back(G2XAB);
			Symmetry_Functions.push_back(G2XBB);
		}	
		*neighbor_it++;
		Symmetry.push_back(Symmetry_Functions);
	}
}


}
#endif
