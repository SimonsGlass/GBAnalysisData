/*
 *  resources.h
 *  Percolation
 *
 *  Created by Samuel Schoenholz on 6/2/10.
 *  Copyright 2010 Swarthmore College. All rights reserved.
 *
 */

#ifndef RESOURCES
#define RESOURCES

#include "std_include.h"

/** @file */ 


// @file A file for commonly used global functions.


template <typename T> inline T POW2(T x){	return x*x;}	//!<x^2
template <typename T> inline T POW3(T x){	return x*x*x;}	//!<x^3
template <typename T> inline T POW4(T x) //!<x^4 
{
	T temp = x*x;
	return temp*temp;
}
template <typename T> inline T POW6(T x) //!<x^6
{
	T temp = x*x*x;
	return temp*temp;
}
template <typename T> inline T POW12(T x) //!<x^12
{
	T temp = x*x*x;
	temp = temp*temp;
	return temp*temp;
}


//! Get the largest absolute value from an array
template <typename T>
inline T max_abs_element(int N, T const *const a)
{
	T m = fabs(a[0]);
	for(int ii=1; ii<N; ii++)   if(fabs(a[ii])>m)   m=fabs(a[ii]);
	return m;
}
template <typename T>
inline T max_abs_element(Eigen::Matrix<T,Eigen::Dynamic,1> const &a)
{
	return max_abs_element(a.size(), a.data());
}

template <class ForwardIterator>
inline dbl mean(ForwardIterator first, ForwardIterator last)
{
	dbl m = (dbl)(*first);
	int Ndata = 1;
	while(++first!=last)
	{
		++Ndata;
		m += ( ((dbl)(*first))-m )/((dbl)Ndata);
	}
	return m;
}




//returns true if there exists an integer x such that x^pow == n
//the return value of x is 1 if the return value is false or if pow == 0;
inline bool IsPerfectRoot(int n, int power, int &x) 
{
	x=1;
	if(power<0) return false;
	if(power==0) return n==1;
	if(power%2==0 && n<0) return false;
	
	x = round( pow(((dbl)n),1./((dbl)power)) );
	int xpow(1);
	for(int i=0; i<power; ++i)
		xpow *= x;
	if(xpow==n)
		return true;
	x=1;
	return false;
}
inline bool IsPerfectRoot(int n, int power)
{
	int x;
	return IsPerfectRoot(n, power, x); 
}


inline std::string ConvertDblToHexString(dbl x)
{
	char cstring[256];
	sprintf(cstring, "%a", x);
	return std::string(cstring);
}

inline dbl ConvertHexStringToDbl(std::string s)
{
	dbl x;
	sscanf(s.c_str(), "%lf", &x);
	return x;
}








//! Check if a file exists.
inline bool FileExists(std::string strFilename)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;
	
	// Attempt to get the file attributes 
	intStat = stat(strFilename.c_str(),&stFileInfo);
	if(intStat == 0) {
		blnReturn = true;
	} else {
		blnReturn = false;
	}
	return(blnReturn);
}
















//function to output a vector<T> as a mathematica vector
template<class T>
ostream &operator<<(ostream &out, vector<T> &data)
{
	out << "{";
	typename vector<T>::iterator it;
	for(it = data.begin();it!=data.end();it++)
		out << (*it) << ((it+1!=data.end())?"," : "");
	
	out << "}";
	return out;
}

//function to input a vector<T> as a mathematica vector
template<class T>
istream &operator>>(istream &in, vector<T> &data)
{
    data.clear();
	char current;
	in >> current;
	if(data.size()==0){
		while(current!='}')
		{
			T temp;
			in >> temp;
			data.push_back(temp);
			in >> current;
		}
	}else {
		int t = 0;
		while(current!='}')
		{
			in >> data[t]; 
			t++;
			in >> current;
		}
	}
    
	return in;
}


//function to perform arithmatic on vector<T>.
template<class T>
vector<T> &operator/(vector<T> &vec, T &div)
{
	typename vector<T>::iterator it;
	for(it = vec.begin();it!=vec.end();it++)
		(*it) = (*it)/div;
	return vec;
}

template<class T>
vector<T> &operator+(vector<T> &vec, vector<T> &add)
{
	vector<T> ret;
	for(int it = 0;it<vec.size();it++)
		ret.push_back(vec[it]+add[it]);
	return ret;
}

template<class T>
vector<T> &operator/=(vector<T> &vec,const T &div)
{
	typename vector<T>::iterator it;
	for(it = vec.begin();it!=vec.end();it++)
		(*it) = (*it)/div;
	return vec;
}

//code to do linear regressions
struct LinearRegression {
    dbl slope;
    dbl intercept;
    
    LinearRegression() : slope(0), intercept(0) {}
    LinearRegression(dbl _i, dbl _s) : slope(_s), intercept(_i) {}
};

LinearRegression LinearFit(vector<dbl> x_data, vector<dbl> y_data);

//code to compute factorials as well as absolute values and sgns.
long int factorial(long int n);

dbl logfactorial(int n);

dbl logchoose(int n, int k);

dbl choose( int n,  int k);

int pow(int n, int e);

template <class T>
T sgn(T a)
{
    return (a>0) - (a<0);
}

template <class T>
T abs(T a)
{
    return (a>0) ? a : -a;
}


//code to calculate the volume and surface area of an n-dimensional unit sphere
static dbl nSphere_Sn(int n);
static dbl nSphere_Vn(int n)
{
	if(n==0) return 1.;
	return nSphere_Sn(n-1)/((dbl)n);
}
dbl nSphere_Sn(int n)
{
	if(n==0) return 2.;
	return 2.*M_PI*nSphere_Vn(n-1);
}

//code to normalize a field
// TAS added static, following stackoverflow.com/questions/6964819
static void NormalizeField(Eigen::VectorXd &field)
{
	double mean = 0.0;
	double meansq = 0.0;

	for(int i = 0 ; i < field.rows() ; i++)
	{
		mean += field(i);
		meansq += field(i)*field(i);
	}

	mean /= field.rows();
	meansq /= field.rows();

	for(int i = 0 ; i < field.rows() ; i++)
		field(i) = (field(i) - mean)/sqrt(meansq - mean*mean);
}	


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////CODE TO APPROXIMATE GAUSSIAN INVERSE CDF////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////From: http://www.johndcook.com/normal_cdf_inverse.html /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// compute log(1+x) without losing precision for small values of x
dbl LogOnePlusX(dbl x);

dbl RationalApproximation(dbl t);

dbl NormalCDFInverse(dbl p);

vector<string> SplitString(const string &target, char ctoken);
vector<string> SplitString(const string &target, const string &token);

#endif
