#include "index_map.h"

index_map::index_map()
{
	initialize(0);
}

index_map::index_map(int N)
{
	initialize(N);
}

index_map::~index_map()
{
}

void index_map::initialize(int N)
{
	clear_map();
	map.reserve(N);
	mapi.reserve(N);
	for(int i=0; i<N; i++)
	{
		map.push_back(i);
		mapi.push_back(i);
	}
	full_size = N;
}
	
index_map& index_map::operator=(const index_map &rhs) 
{
	if(this != &rhs)
	{
		clear_map();
		map = rhs.map;
		mapi = rhs.mapi;
		full_size = rhs.full_size;
	}
	return *this;
}

bool index_map::operator==(index_map const &other) const
{
	if(full_size!=(other.full_size) || (int)map.size()!=(other.size()))
		return false;
	
	for(int i=0; i<(int)map.size(); i++)
		if(map[i]!=other[i])
			return false;

	for(int i=0; i<full_size; i++)
		if(mapi[i]!=other.inv(i))
			return false;

	return true;
}

bool index_map::operator!=(index_map &other) const
{
	return !(*this==other);
}

/*
int index_map::operator[] (unsigned i) const
{
	if(i>=map.size()) //an unsigned int i must not be less than 0, so no need to check that.
	{	printf("ERROR: map out of bounds\n");	return -2;	}
	return map[i];
}

int index_map::inv(int i) const
{
	if(i<0 || i>=full_size)
	{	printf("ERROR: mapi out of bounds\n");	return -2;	}
	return mapi[i];
}

int index_map::size() const		
{	
	return (int)map.size();	
}
*/

void index_map::clear_map()
{
	map.clear();
	mapi.clear();
}

void index_map::set_mapi()
{
//	mapi.clear();
//	for(int i=0; i<full_size; i++)
//		mapi.push_back(-1);
	mapi.assign(full_size, -1);
	for(int i=0; i<(int)map.size(); i++)
		mapi[map[i]] = i;
}

void index_map::remove_from_map(int p) //p is real object #, not mapped #
{
	int mp = mapi[p];
	assert(mp>=0);
	map.erase(map.begin()+mp);
	mapi[p]=-1;
	for(int i=p+1; i<full_size; i++)
		if(mapi[i]>0)
			mapi[i]--;
}

void index_map::add_index(int i)
{
	map.push_back(i);
	full_size++;
	set_mapi();
}

void index_map::set_map(const vector<bool> &ToBeRemoved)
{
//	assert((int)ToBeRemoved.size() == full_size);
	full_size = (int)ToBeRemoved.size();
	map.clear();
	map.reserve(full_size);
	for(int i=0; i<full_size; ++i)
	{
		if(!ToBeRemoved[i])
			map.push_back(i);
	}
	set_mapi();
}

void index_map::set_map(const vector<int> &new_map)
{
	assert((int)new_map.size() <= full_size);
	
	map = new_map;
	set_mapi();
}

/*
void index_map::remove_many_from_map(vector<int> vp)
{
	int mp;
	sort(vp.begin(), vp.end());
	for(int i=(int)vp.size()-1; i>=0; --i)
	{
		map.erase(map.begin()
	}
}
*/

void index_map::PrintMap() const
{
	printf("Printing index_map:\n");
	int im;
	for(int i=0; i<full_size; ++i)
	{
		im = mapi[i];
		if(im != -1)
			assert( map[im] == i );
		printf("\tmap[%5i] = %i\n", im, i);
	}

}

template<typename T>
void index_map::vector_expand(T const *small, T *large) const
{
	int im;
	for(int i=0; i<full_size; ++i)
	{
		im = inv(i);
		if(im==-1) 
			large[i] = (T)0;
		else
			large[i] = small[im];
	}
}

template<typename T>
void index_map::vector_contract(T const *large, T *small) const
{
	for(int im=0; im<size(); ++im)
		small[im] = large[map[im]];
}












