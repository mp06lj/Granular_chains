
#ifndef METH_METHODS_H
#define METH_METHODS_H

#include<vector>
#include<sstream>
#include<string>
#include "OneDimDoubleVector.h"
namespace meth {

extern int timeSeed();
extern double rand(int&, bool =false);
extern double ranf(int&, bool =false);
extern double qand(int&, bool =false);
extern double gasdev(int&);

extern double absolute(double val);
extern double minval(double val1, double val2);
extern double maxval(double val1, double val2);
extern double maxVec(OneDimDoubleVector vec);
extern double maxAbsVec(OneDimDoubleVector vec);
extern double meanAbsVec(OneDimDoubleVector vec);

extern double zero(const double& time,const vector<double>& params);
extern double poly(const double& time,const vector<double>& params);
extern double sine(const double& time,const vector<double>& params);
extern double cosine(const double& time,const vector<double>& params);


} /* end namespace meth */

#endif

