#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <string>

#include <vector>
#include <algorithm>
#include "methods.h"

namespace meth {

// Random Number Generators

int timeSeed() {

	time_t t;                   // usually an unsigned long int
    time(&t);                   // get seconds since Jan 1, 1900
    tm* pt = gmtime(&t);   // Universal (Greenwich Mean) Time
    return pt->tm_year + 70 * (pt->tm_mon + 12 * (pt->tm_hour +
           23 * (pt->tm_min + 59 * pt->tm_sec)));
}

double rand(int& seed, bool set) {
    if (set)
        srand(seed);
    else
        seed = std::rand();
    return seed / (RAND_MAX + 1.0);
}


double ranf(int& seed, bool set) {
    const int IA = 16807, IC = 2147483647, IQ = 127773, IR = 2836;
    if (!set) {
        int h = seed / IQ;
        int l = seed % IQ;
        seed = IA * l - IR * h;
    }
    if (seed <= 0)
        seed += IC;
    return seed / double(IC);
}

double qand(int& seed, bool set) {
    static unsigned long idum;
    const double TWO_POWER_32 = 4294967296.0;
    if (set) {
        idum = (unsigned long) seed;
        return idum / TWO_POWER_32;
    }
    idum = 1664525L * idum + 1013904223L;
    seed = int(idum);
    return idum / TWO_POWER_32;
}

double gasdev(int& seed) {
     static int iset = 0;
     static double gset;
     double fac, rsq, v1, v2;
     if (iset == 0) {
          do {
               v1 = 2.0*ranf(seed)-1.0;
               v2 = 2.0*ranf(seed)-1.0;
               rsq = v1*v1 + v2*v2;
          } while (rsq >= 1.0 || rsq == 0.0);
          fac = std::sqrt(-2.0*std::log(rsq)/rsq);
          gset = v1*fac;
          iset = 1;
          return v2*fac;
     } else {
          iset = 0;
          return gset;
     }
}
double absolute(double val)
{
	if(val<0)
		return -1.0*val;
	else
		return val;
}
 double minval(double val1, double val2)
 {
	 return (val1<val2?val1:val2);
 }

 double maxval(double val1, double val2)
{
	 return (val1>val2?val1:val2);
 }
double maxVec(OneDimDoubleVector vec)
{
    double mx = vec[0];
    for(unsigned int j=0;j<vec.size();++j)
    if(mx < vec[j])
        mx = vec[j];
    return mx;
}
double maxAbsVec(OneDimDoubleVector vec)
{
    double mx = meth::absolute(vec[0]);
    for(unsigned int j=0;j<vec.size();++j)
    {
        if(mx < meth::absolute(vec[j]));
        mx = meth::absolute(vec[j]);
        }
    return mx;
}
double meanAbsVec(OneDimDoubleVector vec)
{
    double m = 0;
    for(unsigned int j=0;j<vec.size();++j)
        m += absolute(vec[j]);
    m/=((double)vec.size());
    return m;
}
double zero(const double& time,const vector<double>& params)
{
    return 0;
}
double poly(const double& time,const vector<double>& params)
{ //param[0] + time*param[1] + time^2+param[2] etc
    double val=0;
    for(unsigned int j=0;j<params.size();++j)
        val+=pow(time,j)*params[j];
    return val;
}
double sine(const double& time,const vector<double>& params)
{  // A * sine(omega * t - phi) ... A=param[0], omega = param[1], phi = param[2]
    return params[0]*sin(params[1]*time-params[2]);
}
double cosine(const double& time,const vector<double>& params)
{
    return params[0]*cos(params[1]*time-params[2]);
}

} /* end namespace meth */

