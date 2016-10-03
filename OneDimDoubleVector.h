#ifndef ONEDIMDOUBLEVECTOR_H
#define ONEDIMDOUBLEVECTOR_H
//Adam Sokolow 1-31-14
//This class makes a one dimensional vector of doubles
//and then defines a number of operations that you might want for them.
//This class could be made into a templated class that accepted say a 2d or 3d coordinate instead of a 1d coord (double).
//Then you could use these general vectors and operators on them in your force calculation
//but in the code everything is just using a general vector and doesn't care whether its 1d,2d or 3d
//So you could in the Simulation class convert the forces velocities, locations,etc to 3d vectors
//and leave the kinetic energy and potential energy as a vector of scalars, and as long as you
//update your calcForce function to take into account distances and perhaps consider all combinations
// of grain-grain interactions you can model more stuff than before.

#include<vector>
#include<iostream>
using namespace std;

class OneDimDoubleVector
{
    public:
        /** Default constructor */
        OneDimDoubleVector();
        OneDimDoubleVector(int length);
        /** Default destructor */
        virtual ~OneDimDoubleVector();
        /** Copy constructor
         *  \param other Object to copy from
         */
        OneDimDoubleVector(const OneDimDoubleVector& other);
        /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */

        unsigned int size();
        void resize(int newSize);
        double at(unsigned int indx);
        void setAt(unsigned int indx,const double &assignval);
        void addAt(unsigned int indx,const double &assignval);
        double& operator[](unsigned int indx);
        void push_back(double val);
        vector<double> Getvec() { return vec; }
        /** Set vec
         * \param val New value to set
         */
        void Setvec(vector<double> val) { vec = val; }

        OneDimDoubleVector& operator=(const OneDimDoubleVector& other);
        /** Access vec
         * \return The current value of vec
         */

        OneDimDoubleVector& operator += (const double& d);
        OneDimDoubleVector& operator -= (const double& d);

        OneDimDoubleVector& operator += (const OneDimDoubleVector& dv);
        OneDimDoubleVector& operator -= (const OneDimDoubleVector& dv);
        OneDimDoubleVector& operator *= (double d);
        OneDimDoubleVector& operator *= (vector<bool> m);
        OneDimDoubleVector& operator /= (double d);

	friend std::ostream& operator << (std::ostream& stream, const OneDimDoubleVector& v);
        vector<double> vec; //!< Member variable "vec" can be accessed directly which is not clean coding but easier
    protected:
    private:

};

inline OneDimDoubleVector operator + (const OneDimDoubleVector& dv) {
    return dv;
}

extern OneDimDoubleVector operator - (const OneDimDoubleVector& dv);

extern OneDimDoubleVector operator * (const OneDimDoubleVector& dv, double d);

extern OneDimDoubleVector operator * (double d, const OneDimDoubleVector& dv);

extern OneDimDoubleVector operator + (const OneDimDoubleVector& dv, double d);

extern OneDimDoubleVector operator + (double d, const OneDimDoubleVector& dv);


extern OneDimDoubleVector operator - (const OneDimDoubleVector& dv, double d);

extern OneDimDoubleVector operator - (double d, const OneDimDoubleVector& dv);


extern OneDimDoubleVector operator / (const OneDimDoubleVector& dv, double d);

extern OneDimDoubleVector operator / (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2);
extern OneDimDoubleVector operator * (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2);
extern OneDimDoubleVector operator * (const OneDimDoubleVector& v1, const vector<bool> & v2);
extern OneDimDoubleVector operator + (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2);

extern OneDimDoubleVector operator - (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2);

extern bool operator == (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2);
extern bool operator != (const OneDimDoubleVector& v1, const OneDimDoubleVector& v2);

#endif // ONEDIMDOUBLEVECTOR_H
