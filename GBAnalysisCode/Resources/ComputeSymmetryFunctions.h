#ifndef COMPUTE_SYMMETRY_FUNCTIONS

#define COMPUTE_SYMMETRY_FUNCTIONS

#include <list>
#include <vector>
#include "../Resources/std_include.h"
#include "../State/StaticState.h"
//#include "../Resources/Settings.h"

//#define TASDEBUG

namespace LiuJamming {

// Caution, small periodic system could have atoms as themselves as neighbors
// Compute the structure functions using the Grids to help with the neighbor lists.
template <int Dim>
int ComputeSymmetryFunctions(CStaticState<Dim> &State, CSimpleGrid<Dim> *Grid, CSimpleGrid<Dim> *Angular_Grid, Settings &settings, list<int> Particles, list<vector<double> > &Symmetry) // , vector<double> ParticlesCNAValue)
{
	//Load the symmetry function parameters
	vector<double> Mu = settings["RadialMu"]; 
	double Sigma = settings["RadialSigma"]; 
	vector<double> Eta = settings["AngularEta"]; 
	vector<double> Zeta = settings["Zeta"]; 
	vector<double> Lambda = settings["Lambda"];  
	int IncludeCNA = settings["IncludeCNA"]; 
        if (Lambda.size() != Zeta.size()) { printf("Error: Lambda and Zeta different lengths\n"); return 0;} 
        if (Eta.size() != Zeta.size()) { printf("Error: Eta and Zeta different lengths\n"); return 0;} 
	double AngularCutoff = settings["AngularCutoff"]; 
	Eigen::Matrix<double,Dim,1> Displacement_JK;


        long particlenum = 0;
	for(list<int>::iterator it = Particles.begin() ; it!=Particles.end() ; it++)
	{
                particlenum++;
#ifdef TASDEBUG
                if ((particlenum%1000)==1) cout << "Computing symmetry functions for particle " << *it << endl; 
#endif                
		
		vector<double> Symmetry_Functions;


                //The first symmetry function for this particle is its Ovito adaptive CNA value
		//Symmetry_Functions.push_back(ParticlesCNAValue[*it]);
                //if (*it < 10) printf("Particle index *it %i has a CNA value of %i \n", *it, ParticlesCNAvalue);


		//compute radial functions
		for(vector<double>::iterator mu_it = Mu.begin() ; mu_it!= Mu.end() ; mu_it++)
		{	
			double G1 = 0.0;
			const list<CSimpleNeighbor<Dim> > &neighbors = Grid->GetNeighbors(*it);

			for(typename list<CSimpleNeighbor<Dim> >::const_iterator neigh_it = neighbors.begin() ; neigh_it != neighbors.end() ; neigh_it++)
			{
				double norm = (*neigh_it).Displacement.norm();
				G1 += exp(-1.0/2.0/Sigma/Sigma*(norm - *mu_it)*(norm - *mu_it));
			}
			Symmetry_Functions.push_back(G1);
#ifdef TASDEBUG
                        if(particlenum%1000 == 1) printf("G1 = %f", G1);
#endif                        
		}

		//compute angular functions
		for(int a_it = 0 ; a_it < Eta.size() ; a_it++)
		{
                        int nprinted = 0;
	                // cout << "Computing angular function " << a_it << endl;
			double G2 = 0.0;
			double eta = Eta[a_it];
			double zeta = Zeta[a_it];
			double lambda = Lambda[a_it];

			list<CSimpleNeighbor<Dim> > neighbors = Angular_Grid->GetNeighbors(*it);
			
			for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_J = neighbors.begin() ; neighbor_J != neighbors.end() ; neighbor_J++)
			{
				for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_K = neighbors.begin() ; neighbor_K != neighbors.end() ; neighbor_K++)
             		        {
					if((*neighbor_K).j != (*neighbor_J).j)
					{
						State.GetDisplacement((*neighbor_J).j,(*neighbor_K).j,Displacement_JK);
									
						double norm_JK = Displacement_JK.squaredNorm();
						if(norm_JK < AngularCutoff*AngularCutoff)
						{
							norm_JK = sqrt(norm_JK);
							double norm_IJ = (*neighbor_J).Displacement.norm();
                                                        //if (neighbor_J == neighbors.begin() && neighbor_K != (neighbors.end())) printf("example norm_IJ %f \n",norm_IJ);
                                                	double norm_IK = (*neighbor_K).Displacement.norm();
	
							double ctheta = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK; 	
								
                                                        G2 += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1);
						}							
					}
				}
			}

#ifdef TASDEBUG
                        nprinted++; if (nprinted % 10 == 0) printf("G2 = %f", G2);
#endif                        
                        Symmetry_Functions.push_back(G2);
		}


                // finally add on CNA
                if (IncludeCNA == 1) {
                    double cna = (double)(State.GetSingleCNA(*it));
                    cout << "Adding to symmetry, CNA of particle " << *it << " which is " << cna << endl;
                    Symmetry_Functions.push_back(cna);
                }
		
		Symmetry.push_back(Symmetry_Functions);
	}

        return 1;
}

template <int Dim>
int ComputeSymmetryFunctionGradient(CStaticState<Dim> &State,CSimpleGrid<Dim> *Grid, CSimpleGrid<Dim> *Angular_Grid, Settings &settings, list<int> Particles, list<vector<Eigen::Matrix<double, Dim, 1> > > &Symmetry)
{
	//Load the symmetry function
	vector<double> Mu = settings["RadialMu"];
	double Sigma = settings["RadialSigma"];
	vector<double> Eta = settings["AngularEta"];
	vector<double> Zeta = settings["Zeta"];
	vector<double> Lambda = settings["Lambda"]; 
        if (Lambda.size() != Zeta.size()) { printf("Error: Lambda and Zeta different lengths\n"); return 0;}
        if (Eta.size() != Zeta.size()) { printf("Error: Eta and Zeta different lengths\n"); return 0;}
	double AngularCutoff = settings["AngularCutoff"];

	Eigen::Matrix<double,Dim,1> Displacement_JK;

	for(list<int>::iterator it = Particles.begin() ; it!=Particles.end() ; it++)
	{
	//	cout << "Computing symmetry functions for particle " << *it << endl;
		
		vector<Eigen::Matrix<double,Dim, 1> > Symmetry_Functions;

		//compute radial functions
		for(vector<double>::iterator mu_it = Mu.begin() ; mu_it!= Mu.end() ; mu_it++)
		{	
			Eigen::Matrix<double,Dim,1> G1X = Eigen::Matrix<double,Dim,1>::Zero();

			const list<CSimpleNeighbor<Dim> > &neighbors = Grid->GetNeighbors(*it);

			for(typename list<CSimpleNeighbor<Dim> >::const_iterator neigh_it = neighbors.begin() ; neigh_it != neighbors.end() ; neigh_it++)
			{
					double norm = (*neigh_it).Displacement.norm();
					G1X += -(norm-*mu_it)/norm/Sigma/Sigma*(*neigh_it).Displacement*exp(-1.0/2.0/Sigma/Sigma*(norm - *mu_it)*(norm - *mu_it));
			
			}
			Symmetry_Functions.push_back(G1);
		}

		//compute angular functions
		for(int a_it = 0 ; a_it < Eta.size() ; a_it++)
		{
                	// cout << "Computing angular function " << a_it << endl;
			Eigen::Matrix<double,Dim,1> G2X = Eigen::Matrix<double,Dim,1>::Zero();

			double eta = Eta[a_it];
			double zeta = Zeta[a_it];
			double lambda = Lambda[a_it];

			list<CSimpleNeighbor<Dim> > neighbors = Angular_Grid->GetNeighbors(*it);
			
			for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_J = neighbors.begin() ; neighbor_J != neighbors.end() ; neighbor_J++)
			{
				for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_K = neighbors.begin() ; neighbor_K != neighbors.end() ; neighbor_K++)
             		        {
					if((*neighbor_K).j != (*neighbor_J).j)
					{
						State.GetDisplacement((*neighbor_J).j,(*neighbor_K).j,Displacement_JK);
									
						double norm_JK = Displacement_JK.squaredNorm();
						if(norm_JK < AngularCutoff*AngularCutoff)
						{
							norm_JK = sqrt(norm_JK);
							double norm_IJ = (*neighbor_J).Displacement.norm();
                                                	double norm_IK = (*neighbor_K).Displacement.norm();
	
							double ctheta = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK;
							double dot = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement);							
											
							double sym = 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1);
						
							Eigen::Matrix<double,Dim,1> grad = (-4.0*(*neighbor_J).Displacement*eta + 2.0*zeta*lambda/norm_IJ/norm_IK*((*neighbor_K).Displacement - dot*(*neighbor_J).Displacement/norm_IJ/norm_IJ)/(1+lambda*ctheta) - 2.0*(*neighbor_J).Displacement/norm_IJ*PI/AngularCutoff*sin(PI*norm_IJ/AngularCutoff)/(1+cos(PI*norm_IJ/AngularCutoff)))*sym;
							
                                                        G2X += grad;
							/* bi-disperse
                                                        if((*neighbor_J).j < Number_A)
							{
								if((*neighbor_K).j< Number_A)
                                                                	G2XAA += grad;
								else
									G2XAB += grad;
							}else{
								if((*neighbor_K).j < Number_A)
									G2XAB += grad;
								else
									G2XBB += grad;
							}
                                                        */
						}							
					}
				}
				
				list<CSimpleNeighbor<Dim> > neighbors_J = Angular_Grid->GetNeighbors((*neighbor_J).j);
				
				for(typename list<CSimpleNeighbor<Dim> >::iterator neighbor_K = neighbors_J.begin() ; neighbor_K != neighbors_J.end() ; neighbor_K++)
             		        {
					if((*neighbor_K).j != (*neighbor_J).i)
					{
						State.GetDisplacement((*neighbor_J).i,(*neighbor_K).j,Displacement_JK);
									
						double norm_JK = Displacement_JK.squaredNorm();
						if(norm_JK < AngularCutoff*AngularCutoff)
						{
							norm_JK = sqrt(norm_JK);
							double norm_IJ = (*neighbor_J).Displacement.norm();
                                                	double norm_IK = (*neighbor_K).Displacement.norm();
	
							double ctheta = -(*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK;
							double dot = -(*neighbor_J).Displacement.dot((*neighbor_K).Displacement);							
											
							double sym = 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1);
						
							Eigen::Matrix<double,Dim,1> grad = (-2.0*(*neighbor_J).Displacement*eta + zeta*lambda/norm_IJ/norm_IK*((*neighbor_K).Displacement - dot*(*neighbor_J).Displacement/norm_IJ/norm_IJ)/(1+lambda*ctheta) - 2.0*(*neighbor_J).Displacement/norm_IJ*PI/AngularCutoff*sin(PI*norm_IJ/AngularCutoff)/(1+cos(PI*norm_IJ/AngularCutoff)))*sym;
		
                                                        G2X += 2*grad;
							/* bi-disperse
                                                        if((*neighbor_J).i < Number_A)
							{
								/star
								if((*neighbor_K).j < Number_A)
									G2XAA += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								else
									G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								star/
								
								if((*neighbor_K).j< Number_A)
                                                                	G2XAA += 2*grad;
								else
									G2XAB += 2*grad;
							}else{
								if((*neighbor_K).j < Number_A)
									G2XAB += 2*grad;
								else
									G2XBB += 2*grad;
								/star
								if((*neighbor_K).j < Number_A)
                                                                        G2XAB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
                                                                else
                                                      	                G2XBB += pow(2.0,1.0-zeta) * exp(-(norm_IJ*norm_IJ + norm_IK*norm_IK + norm_JK * norm_JK)/eta/eta/2.0)*pow(1+lambda*ctheta,zeta);
								star/
							}
                                                        */
						}							
					}
				}
			}
			
			Symmetry_Functions.push_back(G2X);
		}	
	
                // finally add on CNA symmetry function.
                // WHAT WOULD THE GRADIENT BE?  FIX THIS IF NEED GRADIENT
                int cna = State.GetSingleCNA(*it);
                cout << "Adding to symmetry, CNA of particle " << *it << " which is " << cna << endl;
                Symmetry_Functions.push_back(cna);

		Symmetry.push_back(Symmetry_Functions);
	}

        return 1;
}

// Compute the structure functions of the Particles list given the NeighborList
template <int Dim>
int ComputeSymmetryFunctions(CStaticState<Dim> &State,list<list<CSimpleNeighbor<Dim> > > &NeighborList, Settings &settings, list<int> Particles, list<vector<double> > &Symmetry)
{
	//Load the symmetry function
	vector<double> Mu = settings["RadialMu"];
	double Sigma = settings["RadialSigma"];
	vector<double> Eta = settings["AngularEta"];
	vector<double> Zeta = settings["Zeta"];
	vector<double> Lambda = settings["Lambda"]; 
        if (Lambda.size() != Zeta.size()) { printf("Error: Lambda and Zeta different lengths\n"); return 0;}
        if (Eta.size() != Zeta.size()) { printf("Error: Eta and Zeta different lengths\n"); return 0;}
	double AngularCutoff = settings["AngularCutoff"];
	int IncludeCNA = settings["IncludeCNA"];

        Eigen::Matrix<double,Dim,1> Displacement_JK;

        typename list<list<CSimpleNeighbor<Dim> > >::iterator neighbor_it = NeighborList.begin();

	for(list<int>::iterator it = Particles.begin() ; it!=Particles.end() ; it++)
	{       // How do we know that the neighbor list corresponds to this particle?
		//cout << "Computing symmetry functions for particle " << *it << endl;
		
		vector<double> Symmetry_Functions;

		//compute radial functions
		for(vector<double>::iterator mu_it = Mu.begin() ; mu_it!= Mu.end() ; mu_it++)
		{	
			double G1 = 0.0;
                        int nprinted =0;
                        // TAS I think this gets the list from NeighborList for the partile "it"...
			list<CSimpleNeighbor<Dim> > &neighbors = *neighbor_it;
                        // TAS I should verify this with test cases with cout << "May be worth while to " << neighbors << endl;

			for(typename list<CSimpleNeighbor<Dim> >::iterator neigh_it = neighbors.begin() ; neigh_it != neighbors.end() ; neigh_it++)
			{
				double norm = (*neigh_it).Displacement.norm();
                                //if ((int)(G1 *100000) % 1000 == 0) printf("norm of neighbors is numbers like %f, ensure x not xs \n", norm);
				G1 += exp(-1.0/2.0/Sigma/Sigma*(norm - *mu_it)*(norm - *mu_it));
                                //if (nprinted == 0) { cout << "check "<<norm<<" "<<Sigma<<" " <<*mu_it << endl;}
			}
			Symmetry_Functions.push_back(G1);
#ifdef TASDEBUG
                        nprinted++; if (nprinted  == 1) printf("G1 = %f \n",G1);
#endif                        
		}

		//compute angular functions
		for(int a_it = 0 ; a_it < Eta.size() ; a_it++)
		{
			//cout << "Computing angular function " << a_it << endl;
			double G2 = 0.0; // TAS added, deleting G2XAA, G2XAB, G2XBB

			double eta = Eta[a_it]; // I am going to say that eta is 1/exsi^2 // TAS
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
						if(norm_JK < AngularCutoff && norm_IJ < AngularCutoff && norm_IK < AngularCutoff)
						{
							double ctheta = (*neighbor_J).Displacement.dot((*neighbor_K).Displacement)/norm_IJ/norm_IK; 
                                                        G2 += 0.5*0.5*0.5*pow(2.0,1.0-zeta)*pow(1+lambda*ctheta,zeta)*exp(-eta*(norm_IJ*norm_IJ+norm_IK*norm_IK+norm_JK*norm_JK))*(cos(PI*norm_IJ/AngularCutoff)+1)*(cos(PI*norm_IK/AngularCutoff)+1)*(cos(PI*norm_JK/AngularCutoff)+1); // Why is there is a 2^(1-z) prefactor of G2?
						}
					}
				}
			}
			
			Symmetry_Functions.push_back(G2);
#ifdef TASDEBUG
                        cout << "G2 " << G2 << endl;
#endif                        
		}	

                if (IncludeCNA == 1) {
                    // finally add on CNA symmetry function
                    double cna = (double)(State.GetSingleCNA(*it));
                    cout << "Adding to symmetry, CNA of particle " << *it << " which is " << cna << endl;
                    Symmetry_Functions.push_back(cna);
                }

		*neighbor_it++;
		Symmetry.push_back(Symmetry_Functions);
	}
        return 1;
}


}
#endif
