#ifndef NEIGHBORS_H
#define NEIGHBORS_H


namespace LiuJamming
{

//Remove particle i from all of it's neighbors' nbrs list. Then clear its own nbrs list.
void RemoveFromNeighborsList(std::vector< std::vector<int> > &nbrs, int i)
{
	int j, nn; 
	for(int nni=0; nni<(int)nbrs[i].size(); nni++)
	{
		nn = nbrs[i][nni];
		for(j=0; j<(int)nbrs[nn].size(); j++)
			if(nbrs[nn][j]==i)
			{
				nbrs[nn].erase(nbrs[nn].begin()+j);
				break;
			}
	}
	nbrs[i].clear();
}


//Check that if i points to j, then j points to i
bool CheckNeighborsConsistency(std::vector< std::vector<int> > const &nbrs)
{
	int N = (int)nbrs.size();
	int nn;
	int nbr_found;
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<(int)nbrs[i].size(); j++)
		{
			nn = nbrs[i][j];
			nbr_found = 0;
			for(int k=0; k<(int)nbrs[nn].size(); k++)
				if(nbrs[nn][k] == i)
				{   nbr_found = 1; break;   }
			if(!nbr_found)
			{
				return false;
			}
		}
	}
	return true;
}










}

#endif //NEIGHBORS_H








