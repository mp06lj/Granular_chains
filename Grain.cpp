//Adam Sokolow 1/31/14
//Grain class really is a nice bundle on the grains
//in case int eh future we change how we handle them

#include "Grain.h"
#include <cmath>
#include <iostream>
using namespace std;

Grain::~Grain()
{
    //dtor
}

Grain::Grain(const Grain& other)
{
    //copy ctor
    mass = other.mass;
    myShape = other.myShape;
    poissons = other.poissons;
    youngs = other.youngs;
    volume = other.volume;
    rho = other.rho;
}

Grain& Grain::operator=(const Grain& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator

    mass = rhs.mass;
    myShape = rhs.myShape;
    poissons = rhs.poissons;
    volume = rhs.volume;
    rho = rhs.rho;
    youngs = rhs.youngs;
    return *this;
}
bool Grain::SetbyString(string param,string val)
{
    if(param == "radiusplus")
        {myShape.r_plus = atof(val.c_str());myShape.sphere=false;}
    else if(param == "radiusminus")
        {myShape.r_minus = atof(val.c_str());myShape.sphere=false;}
    else if(param == "radius")
        {myShape.r_minus = atof(val.c_str());
        myShape.r_plus = atof(val.c_str());
        }
    else if(param == "widthparam")
        {myShape.width = atof(val.c_str());myShape.sphere=false;}

    else if(param == "alphaplus")
        {myShape.alpha_plus = atof(val.c_str());myShape.sphere=false;}
    else if(param == "alphaminus")
        {myShape.alpha_minus = atof(val.c_str());myShape.sphere=false;}
    else if(param == "alpha")
        {
            myShape.alpha_plus = atof(val.c_str());
            myShape.alpha_minus = atof(val.c_str());
            myShape.sphere=false;
        }
    else if(param == "mparam")
        {myShape.m_param = atof(val.c_str());myShape.sphere=false;}
     else if(param == "volume")
        volume = atof(val.c_str());
     else if(param == "rho")
        rho = atof(val.c_str());
    else if(param == "youngs")
        youngs = atof(val.c_str());
    else if(param == "poisson")
        poissons = atof(val.c_str());
    else if(param == "sphere")
        {//special case I guess
            myShape.sphere=true;
            myShape.r_minus = atof(val.c_str());
            myShape.r_plus = atof(val.c_str());
            myShape.width = 0;
            myShape.m_param = 1;
            myShape.alpha_minus =2;
            myShape.alpha_plus =2;
            volume = 4.0/3.0*(4 * atan(1.0))*pow(myShape.r_plus,3.0);
         //   cout<<"You are a sphere"<<endl;
        }
    else
       {
                cout<<"******************************************\nMost likely a typo in your Subunit Command!\n check spelling on : "<<param<<"\n******************************************"<<endl<<endl;
                return false;
       }
//If I was doing this better I wouldn't have output here in the Grain class but rather send back strings and bools
return true;
}


void Grain::computeTheRest()
{
    if(mass !=0 || volume !=0) //they were specified
    {
        if(mass!=0)
            volume = mass/rho;
        else
            mass = volume*rho;
    }
    else  //use the shape info to calculate volume
    {
        volume = calculateVolume(myShape);
        mass = volume*rho;
    }
}
double Grain::calculateVolume(Shape &shape)
{

   if(!shape.sphere)
   {
    ///*******************************************************
    ///Current volume formulation assumes that alpha is equal on both sides
    double alpha = (shape.alpha_minus+shape.alpha_plus)*0.5;
    ///*******************************************************

    double tempVol = 4*pow(shape.m_param,1-1.0/(alpha))*twoFone(alpha);
    double Rbig = (shape.r_minus>shape.r_plus ? shape.r_minus: shape.r_plus);

	double term1 = shape.width*pow(Rbig,2);

	double term2 = pow(shape.r_minus,2-2.0/alpha)*pow(shape.r_minus-shape.width/2.0,2.0/alpha+1)*pow(2.0/alpha+1,-1.0);
	double term3 = pow(shape.r_plus,2-2.0/alpha)*pow(shape.r_plus-shape.width/2.0,2.0/alpha+1)*pow(2.0/alpha+1,-1.0);

	return tempVol*(term1+term2+term3);
   }
   else return 4.0/3.0*(4 * atan(1.0))*pow(myShape.r_plus,3.0);
}

double Grain::twoFone(double alpha)
{
    return std::pow(calcGamma((1+alpha)/alpha),2)/calcGamma((2+alpha)/alpha);
}
double Grain::simpsonIntegrate(double x,double a, double b)
{
    return (b-a)/6.0*(bodyGamma(x,a)+4*bodyGamma(x,(a+b)/2.0)+bodyGamma(x,b));
}
double Grain::calcGamma(double x)
{
    double multfactor=1;
	while (x>2)
	{
		multfactor = multfactor*(x-1);
		x=x-1;
		//std::cout<<x<<std::endl;
	}
	while (x<1)
	{
		multfactor = multfactor/(x);
		x=x+1;
		//std::cout<<x<<std::endl;
	}
	double sum=0;
	for(int lcv=0; lcv<100000; ++lcv)
	{

		sum=sum+simpsonIntegrate(x,lcv*24/100000.0,(lcv+1)*24/100000.0);
	}
	return sum*multfactor;
}
double Grain::bodyGamma(double x,double t)
{
 return std::exp(-t)*std::pow(t,x-1);
}
double Grain::calculateForcePrefactor(const Grain &gtwo)
{
    ///assume that grain one is to the right of grain two
    ///right now average the two, but in the future this could be changed
    if(myShape.alpha_plus != gtwo.myShape.alpha_minus)
    cout<<"\n***********************************************\nWarning: your alphas are not matched and I will have to average them!\n***********************************************\n"<<endl;

    double alpha = (myShape.alpha_plus+gtwo.myShape.alpha_minus)*0.5;
    double m = (myShape.m_param+gtwo.myShape.m_param)*0.5;
    double Y1 = youngs;
    double Y2 = gtwo.youngs;
    double sig1 = poissons;
    double sig2 = gtwo.poissons;
    double R1 = myShape.r_plus;
    double R2 = gtwo.myShape.r_minus;

    double tempA = twoFone(alpha)*(((m+1)*2*Y1*Y2)/((1-sig1*sig1)*Y2+(1-sig2*sig2)*Y1))/(1.0/alpha+1);

    if(R1 == -1.0 || R2 == -1.0) //flat wall condition
    {
        R1 = max(R1,R2);
        return tempA * pow(R1,(alpha-1.0)/(alpha));
    }
    else
	return tempA * R1*R2*pow(pow(R1,alpha)*R2+pow(R2,alpha)*R1,-1.0/alpha);

}
double Grain::calculateForcePower(const Grain &gtwo)
{
    ///assume that grain one is to the right of grain two
    ///right now average the two, but in the future this could be changed
    if(myShape.alpha_plus != gtwo.myShape.alpha_minus)
    cout<<"\n***********************************************\nWarning: your alphas are not matched and I will have to average them!\n***********************************************\n"<<endl;

    double alpha = (myShape.alpha_plus+gtwo.myShape.alpha_minus)*0.5;
    return 1.0/alpha + 1.0;
}
