#include <fstream>
#include <vector>
#include <string>

using namespace std;

bool ImportTimesteps(string filename,vector<int> &timesteps)
{
	timesteps.clear();

	ifstream in(filename.c_str());
	if(!in.is_open())
	{
		cout << "Could not open timestep file: " << filename << endl;
		return false;
	}

	int i;
	while(in >> i)
		timesteps.push_back(i);

	in.close();

	return true;
}
