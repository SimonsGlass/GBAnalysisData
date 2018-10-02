#ifndef INDEX_MAP_H
#define INDEX_MAP_H

#include <vector>
#include <assert.h>
#include <stdio.h>
using namespace std;
/**********************************************************************************************
 *		class:			index_map
 *		written by:		Carl P. Goodrich
 *		last edited:	4/9/11
 *
 *		purpose:		index_map is designed for removing objects from a list. If you
 *						want to remove m objects from a set of N objects, index_map
 *						will allow you to loop over only the N-m objects that have not been
 *						removed. map[i] runs from i=0 to i=N-m-1 and is the equal to the index
 *						of the ith object that has not been removed. mapi[i] runs from i=0
 *						to i=N-1 and is the inverse of map. Note: if the ith object has been
 *						removed, mapi[i] will equal -1.
 *
 *						For example, a set of 4 objects with the 2nd one removed:
 *
 *								map[0] = 0		mapi[0] = 0
 *								map[1] = 2		mapi[1] = -1
 *								map[2] = 3		mapi[2] = 1
 *												mapi[3] = 2
 *
 *		disclamer:		This code is not a black box and should not be treated as such.
 *						Please report any bugs to cpgoodri@sas.upenn.edu.
 *
 **********************************************************************************************/

class index_map
{
private:
	vector<int> map;
	vector<int> mapi;	//inverse map

public:
	int full_size;

	index_map();		//Initializes the map: trivial configuration with size 0.
	index_map(int N);	//Initializes the map: trivial configuration with size N.
	~index_map();

	void initialize(int N);

	index_map & operator=(const index_map &rhs);
	bool operator==(index_map const &other) const;
	bool operator!=(index_map &other) const;

	inline int operator[] (unsigned i) const	{ return map .at(i); };
	inline int inv(int i) const			{ return mapi.at(i); };
	inline int size() const				{ return (int)map.size(); };

	void clear_map();		//clears the map.

	void set_mapi();
	void remove_from_map(int p); //p is real object #, not mapped #. This is not very efficient.
	void add_index(int i);

	void set_map(const vector<bool> &ToBeRemoved);
	void set_map(const vector<int> &new_map);

	void PrintMap() const;

	template<typename T> void vector_expand(T const *small, T *large) const;  //!< Convert a vector from mapped indices to real indices, setting undefined elements to (T)0
	template<typename T> void vector_contract(T const *large, T *small) const;  //!< Convert a vector from real indices to mapped indices.

};

#endif //INDEX_MAP_H
