#ifndef GRAIN_H
#define GRAIN_H

#include<sstream>
#include<string>
#include <stdlib.h>
//Adam Sokolow - 1/31/14
//The newest revision of the taperchain code
//This class defines a grain in the chain.
//Each grain has material properties and a Shape
//The Shape has all the info use to calculate
//the coefficients that appear between grains
//that are in contact.
using namespace std;

struct Shape;
struct Shape
{
    double r_minus,r_plus,width;
    double alpha_minus,alpha_plus; //right now these are assumed equal, but this could be updated later
    double m_param;
    bool sphere;
    Shape(double rm,double w, double rp,double ap,double am,double m ,bool sp){r_minus=rm;r_plus=rp;width=w; alpha_minus=am;alpha_plus=ap;m_param=m;sphere=sp;}
    Shape(){r_minus=1;r_plus=1;width=0;alpha_minus=2;alpha_plus=2;m_param=1;sphere=false;}
    //Shape& operator=(Shape& a) { return a; }
    Shape& operator=(const Shape& a) {
    r_minus = a.r_minus;r_plus = a.r_plus; width  = a.width;alpha_minus=a.alpha_minus;alpha_plus=a.alpha_plus;m_param=a.m_param;sphere=a.sphere;
    return *this;
  }
    friend std::ostream& operator<< (std::ostream& stream, const Shape& mys)
    {
        stream<<"Shape Parameter"<<std::endl;
        if(!mys.sphere)
        {
        stream<<"\t R- : "<<mys.r_minus<<std::endl;
        stream<<"\t width : "<<mys.width<<std::endl;
        stream<<"\t R+ : "<<mys.r_plus<<std::endl;
        stream<<"\t alpha- : "<<mys.alpha_minus<<std::endl;
        stream<<"\t alpha+ : "<<mys.alpha_plus<<std::endl;
        stream<<"\t m : "<<mys.m_param<<std::endl;
        }
        else
        {
            stream<<"\t Spherical Grain R : "<<mys.r_minus<<std::endl;
            stream<<"\t m : "<<mys.m_param<<std::endl;
            stream<<"\t alpha : "<<mys.alpha_plus<<std::endl;
        }
        return stream;
    }
};


class Grain
{
    public:
        /** Default constructor */
        Grain(){poissons=0;mass=0;volume=0;rho=0;youngs=0;myShape = Shape();};
        /** Default destructor */
         ~Grain();
        /** Copy constructor
         *  \param other Object to copy from
         */
        Grain(const Grain& other);
        /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */
        Grain& operator=(const Grain& other);
        /** Access myShape
         * \return The current value of myShape
         */
        Shape GetmyShape() { return myShape; }
        /** Set myShape
         * \param val New value to set
         */
        void SetmyShape(Shape val) { myShape = val; }
        /** Access rho
         * \return The current value of rho
         */
        double Getrho() { return rho; }
        /** Set rho
         * \param val New value to set
         */
        void Setrho(double val) { rho = val; }
        /** Access youngs
         * \return The current value of youngs
         */
        double Getyoungs() { return youngs; }
        /** Set youngs
         * \param val New value to set
         */
        void Setyoungs(double val) { youngs = val; }
        /** Access poissons
         * \return The current value of poissons
         */
        double Getpoissons() { return poissons; }
        /** Set poissons
         * \param val New value to set
         */
        void Setpoissons(double val) { poissons = val; }
        /** Access volume
         * \return The current value of volume
         */
        double Getvolume() { return volume; }
        /** Set volume
         * \param val New value to set
         */
        void Setvolume(double val) { volume = val; }
        /** Access mass
         * \return The current value of mass
         */
        double Getmass() { return mass; }
        /** Set mass
         * \param val New value to set
         */
        void Setmass(double val) { mass = val; }

        bool SetbyString(string param,string val);

        void computeTheRest();

        double calculateForcePrefactor(const Grain &gtwo);
        double calculateForcePower(const Grain &gtwo);

        friend std::ostream& operator<< (std::ostream& stream, const Grain& grain)
        {
            stream<<endl<<"Grain Youngs Modulus : "<<grain.youngs<<std::endl;
            stream<<"Grain Poisson Ratio : "<<grain.poissons<<std::endl;
            stream<<"Grain Density : "<<grain.rho<<std::endl;
            stream<<"Grain Volume : "<<grain.volume<<std::endl;
            stream<<"Grain mass : "<<grain.mass<<std::endl;
            stream<<"Grain Shape : "<<grain.myShape<<std::endl;

            return stream;
        }


    protected:
    private:
        Shape myShape; //!< Member variable "myShape"
        double rho; //!< Member variable "rho"
        double youngs; //!< Member variable "youngs"
        double poissons; //!< Member variable "poissons"
        double volume; //!< Member variable "volume"
        double mass; //!< Member variable "mass"

    //shape functions
    double twoFone(double alpha);
    double simpsonIntegrate(double x,double a, double b);
    double calcGamma(double x);
    double bodyGamma(double x,double t);
    double calculateVolume(Shape &shape);

};

#endif // GRAIN_H
