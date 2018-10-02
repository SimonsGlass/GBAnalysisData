
#include "Resources.h"

LinearRegression LinearFit(vector<double> x_data, vector<double> y_data)
{
    double x_avg=0.0, y_avg=0.0, xx_avg=0.0, xy_avg=0.0;
    for(int i = 0;i<x_data.size();i++)
    {
        x_avg+=x_data[i];
        y_avg+=y_data[i];
        xx_avg+=x_data[i]*x_data[i];
        xy_avg+=x_data[i]*y_data[i];
    }
    
    x_avg/=x_data.size();
    y_avg/=x_data.size();
    xx_avg/=x_data.size();
    xy_avg/=x_data.size();
    
    LinearRegression ret;
    ret.slope = (xy_avg-x_avg*y_avg)/(xx_avg-x_avg*x_avg);
    ret.intercept = y_avg - ret.slope*x_avg;
    return ret;
}

long int factorial(long int n)
{
    if(n<=1) return 1;
    long int t = n;
    for(int i = n-1;i>0;i--)
        t*=i;
    return t;
}

double logfactorial(int n)
{
    if(n<=1) return 0;
    double t = log(n);
    for(int i = n-1;i>0;i--)
        t+=log(i);
    return t;
}

double logchoose(int n, int k)
{
    return logfactorial(n)-logfactorial(k)-logfactorial(n-k);
}

double choose( int n,  int k)
{
    return exp(logchoose(n,k)-logchoose(n,n/2));
}

int pow(int n, int e)
{
    if(e==1)
        return n; 
    return e>2 ? n*pow(n,e-1) : n*n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////CODE TO APPROXIMATE GAUSSIAN INVERSE CDF////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////From: http://www.johndcook.com/normal_cdf_inverse.html /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// compute log(1+x) without losing precision for small values of x
double LogOnePlusX(double x)
{
    if (x <= -1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << x 
        << "); must be greater than -1.0";
        throw std::invalid_argument( os.str() );
    }
    
    if (fabs(x) > 1e-4)
    {
        // x is large enough that the obvious evaluation is OK
        return log(1.0 + x);
    }
    
    // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
    // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
    
    return (-0.5*x + 1.0)*x;
}

double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) / 
    (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p <= 0.0 || p >= 1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << p 
        << "); must be larger than 0 but less than 1.";
        throw std::invalid_argument( os.str() );
    }
    
    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}

vector<string> SplitString(const string &target, char ctoken)
{
	string token(&ctoken);
	return SplitString(target, token);
}

vector<string> SplitString(const string &target, const string &token)
{
	vector<string> split;
	string current;
	for(int i = 0 ; i < target.length() ; i++)
	{
		if(target.compare(i,1,token)==0)
		{
			if(current.length()>0){
				split.push_back(current);
				current = "";
			}
		}else{
			current.push_back(target[i]);
		}
	}
	
	if(current.length()>0){
		split.push_back(current);
		current = "";
	}
	return split;
}



