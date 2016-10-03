#include "OneDimDoubleVector.h"
#include <vector>
//Adam Sokolow 1-31-14
//This class makes a one dimensional vector of doubles
//and then defines a number of operations that you might want for them.
using namespace std;

OneDimDoubleVector::OneDimDoubleVector()
{
    //ctor
    vec.clear();
}
OneDimDoubleVector::OneDimDoubleVector(int length)
{
    vec.resize(length);
}
void OneDimDoubleVector::resize(int newSize)
{
    vec.resize(newSize);
}
OneDimDoubleVector::~OneDimDoubleVector()
{
    //dtor
}

OneDimDoubleVector::OneDimDoubleVector(const OneDimDoubleVector& other)
{
    //copy ctor
    vec.clear();
    vec = other.vec;
}
unsigned int OneDimDoubleVector::size()
{
    return vec.size();
}
double OneDimDoubleVector::at(unsigned int indx)
{   //returns a copy
    return vec[indx];
}
void OneDimDoubleVector::setAt(unsigned int indx,const double &assignval)
{
    vec[indx]=assignval;
}
void OneDimDoubleVector::addAt(unsigned int indx,const double &assignval)
{
    vec[indx]+=assignval;
}

OneDimDoubleVector& OneDimDoubleVector::operator=(const OneDimDoubleVector& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    vec.clear();
    vec = rhs.vec;

    //assignment operator
    return *this;
}
double& OneDimDoubleVector::operator[](unsigned int indx)
{
    return vec[indx];
}

void OneDimDoubleVector::push_back(double val)
{
    vec.push_back(val);
}
///**********************************************************************

ostream& operator << (ostream &stream, const OneDimDoubleVector &v)
{
    for(unsigned int j=0;j<v.vec.size()-1;++j)
        stream<<v.vec[j]<< "\t";
    stream<<v.vec[v.vec.size()-1];
	return stream;
}
OneDimDoubleVector& OneDimDoubleVector::operator += (const double& d)
{
    for(unsigned int j=0;j<vec.size();++j)
            vec[j]+=d;
 return *this;
}

OneDimDoubleVector& OneDimDoubleVector::operator -= (const double& d)
{
for(unsigned int j=0;j<vec.size();++j)
            vec[j]-=d;
return *this;
}

OneDimDoubleVector& OneDimDoubleVector::operator += (const OneDimDoubleVector& dv)
{
    if(vec.size()==dv.vec.size())
    {
        for(unsigned int j=0;j<vec.size();++j)
            vec[j]+=dv.vec[j];
    }
    else
        cout<<"ERROR YOUR VECTORS ARE NOT THE SAME LENGTH - OneDimDoubleVector.cpp"<<endl;
    return *this;
}

OneDimDoubleVector& OneDimDoubleVector::operator -= (const OneDimDoubleVector& dv) {
    if(vec.size()==dv.vec.size())
    {
        for(unsigned int j=0;j<vec.size();++j)
            vec[j]-=dv.vec[j];
    }
    else
        cout<<"ERROR YOUR VECTORS ARE NOT THE SAME LENGTH - OneDimDoubleVector.cpp"<<endl;
    return *this;
}

OneDimDoubleVector& OneDimDoubleVector::operator *= (double d)
{
    for(unsigned int j=0;j<vec.size();++j)
            vec[j]*=d;
    return *this;

}
OneDimDoubleVector& OneDimDoubleVector::operator *= (vector<bool> m)
{
    for(unsigned int j=0;j<vec.size();++j)
            vec[j]*=(m[j]?1.0:0.0);
    return *this;
}

OneDimDoubleVector& OneDimDoubleVector::operator /= (double d) {
    for(unsigned int j=0;j<vec.size();++j)
            vec[j]/=d;
    return *this;
}

OneDimDoubleVector operator - (const OneDimDoubleVector& dv)
{ //negtaive operator should not change current vector

    OneDimDoubleVector temp;
    for(unsigned int j=0;j<dv.vec.size();++j)
        temp.push_back(-dv.vec[j]);
    return temp;
}

OneDimDoubleVector operator * (const OneDimDoubleVector& dv, double d)
{
    OneDimDoubleVector temp;
        for(unsigned int j=0;j<dv.vec.size();++j)
        temp.push_back(dv.vec[j]*d);
    return temp;
}

OneDimDoubleVector operator * (double d, const OneDimDoubleVector& dv) {
    return dv*d;
}

OneDimDoubleVector operator + (const OneDimDoubleVector& dv, double d)
{
    OneDimDoubleVector temp;
        for(unsigned int j=0;j<dv.vec.size();++j)
        temp.push_back(dv.vec[j]+d);
    return temp;
}

OneDimDoubleVector operator + (double d, const OneDimDoubleVector& dv) {
    return dv+d;
}

OneDimDoubleVector operator - (const OneDimDoubleVector& dv, double d)
{
    OneDimDoubleVector temp;
        for(unsigned int j=0;j<dv.vec.size();++j)
        temp.push_back(dv.vec[j]-d);
    return temp;
}

OneDimDoubleVector operator - (double d, const OneDimDoubleVector& dv) {
    return (-dv)+d;
}



OneDimDoubleVector operator / (const OneDimDoubleVector& dv, double d) {
    return dv*(1/d);
}
OneDimDoubleVector operator / (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2)
{
    OneDimDoubleVector temp;
        for(unsigned int j=0;j<v1.vec.size();++j)
        temp.push_back(v1.vec[j]/v2.vec[j]);
    return temp;
}
OneDimDoubleVector operator * (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2)
{
    OneDimDoubleVector temp;
        for(unsigned int j=0;j<v1.vec.size();++j)
        temp.push_back(v1.vec[j]*v2.vec[j]);
    return temp;

}

OneDimDoubleVector operator * (const OneDimDoubleVector& v1, const vector<bool> & v2)
{
    OneDimDoubleVector temp;
        for(unsigned int j=0;j<v1.vec.size();++j)
        temp.push_back(v1.vec[j]*( v2[j] ? 1.0 : 0.0 ));
    return temp;

}
OneDimDoubleVector operator + (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2)
{
    OneDimDoubleVector temp;
    if(v1.vec.size()==v2.vec.size())
    {
        for(unsigned int j=0;j<v1.vec.size();++j)
        temp.push_back(v1.vec[j]+v2.vec[j]);
    }
    else
        cout<<"ERROR YOUR VECTORS ARE NOT THE SAME LENGTH - OneDimDoubleVector.cpp"<<endl;
    return temp;
}

OneDimDoubleVector operator - (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2)
{
 OneDimDoubleVector temp;
    if(v1.vec.size()==v2.vec.size())
    {
        for(unsigned int j=0;j<v1.vec.size();++j)
        temp.push_back(v1.vec[j]-v2.vec[j]);
    }
    else
        cout<<"ERROR YOUR VECTORS ARE NOT THE SAME LENGTH - OneDimDoubleVector.cpp"<<endl;
    return temp;
}

bool operator == (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2)
{
    if(v1.vec.size()!=v2.vec.size())
        return false;
    for(unsigned int j=0;j<v1.vec.size();++j)
        if(v1.vec[j]!=v2.vec[j])
            return false;
	return true;
}
bool operator != (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2)
{
	return !(v1==v2);
}
