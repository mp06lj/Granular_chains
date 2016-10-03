#ifndef SIMULATION_H
#define SIMULATION_H
//Adam Sokolow 1-31-14
//A simulation class contains all the relevant information about the simulation
//It has the chain, the interfaces, the force law, the integration scheme
//It als handles the text parsing

#include <string>
#include <iostream>
#include <queue>
#include <vector>
#include <map>
#include <list>
#include "Grain.h"
#include "OneDimDoubleVector.h"
#include "methods.h"
using namespace std;

struct DriverFunction;
struct LoadStep;
struct State;
class Simulation;

typedef OneDimDoubleVector StateVector;//in the future one could use a different type
                            //of state vector (say 2d instead!) as long as the force
                            //calculation takes a state vector and returns updated forces
                            //which are also state vectors, everything will be consistent

typedef void (Simulation::*solver)(State&,State&,LoadStep&);
typedef void (Simulation::*timeUpdate)(double&,double&,const double&,const double&,double&);
typedef double (*drivingFunction)(const double&,const vector<double>&);

typedef void (Simulation::*outputter)(StateVector *,string ,bool);//grain or interface

struct outputInfo
{
    friend class Simulation;

  StateVector * svptr;
  string filename;
  outputter myOutputter;

  outputInfo& operator=(const outputInfo& a)
    {
        svptr = a.svptr;
        filename = a.filename;
        myOutputter=a.myOutputter;
    return *this;
    }
};

typedef void (Simulation::*command)(string);


///START OF SIMULATION CLASS ***********************************************
class Simulation
{
    friend struct State;
    friend struct LoadStep;
    friend struct DrivingFunction;
    public:
        /** Default constructor */
        Simulation();
        /** Default destructor */
        ~Simulation();
        /** Copy constructor
         *  \param other Object to copy from
         */
        Simulation(const Simulation& other);
        /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */
        Simulation& operator=(const Simulation& other);

    pair<bool,string> parseString(string original);

    void printHelp();
    bool shouldWeHault(){return haultFlag;}
    void setHault(){haultFlag=true;}
    bool isSilent(){return logFileOnly;}

    void runSimulation();
    void executeCommands();

    friend std::ostream& operator<< (std::ostream& stream, const vector<Grain> &ch);
    friend std::ostream& operator<< (std::ostream& stream, const Simulation &sim);

    void velocityVerletStep(State &curState,State &nextState,LoadStep &details);
    void equilibriumSolverStep(State &curState,State &nextState,LoadStep &details);
    void equilibriumSolverBodyOnly(State &curState,State &nextState,LoadStep &details);

    void rungeKuttaOrderFour(State &curState,State &nextState,LoadStep &details);
    void calculateForces(State &curState,LoadStep &details);


    string toLowerCase(string com);
    string removeWhiteSpace(string com);



    protected:


    private:

    State* myCurrentState;
    State* myPrevState;
    vector<outputInfo> myOutputs;
    priority_queue<pair<int,string> >commandList;

    //Parser variables
    map<string,pair<int,command> > myParser;//int is priority, command is function pointer
    map<string,string> oldCommands;

    //User Function definitions


    //Chain related
    vector<Grain> myChain;
    OneDimDoubleVector massVector;
    map<pair<int,int>,double> forcePrefactors;//given two grain indices, what is their prefactor?
    map<pair<int,int>,double> forcePowerLaw;

    //Output related
    int outputPrecision;
    double outputInterval;
    string baseFileName;
    bool forceOutput,relLocationOutput,velocityOutput,kineticEnergyOutput,
    potentialEnergyOutput,overlapOutput,interfaceForceOutput,absoluteOutput;

   //need function pointer to grain or interface, string of file name,StateVector pointer to force, or xrel or whatever it is

    bool logFileOnly;

    //Pattern related
    vector<int> outputGrainNumbers;//the values of these run from 1 -> N, but internally are 0 -> N-1
    vector<int> outputInterfaceNumbers;//the values of these run from 0->N (0 is leftwall/grain1 interface, 1 is grain1/grain2, ..., N is grain N/RightWal etc)
    map<string,Grain> unitDefinitions;
    vector<string> chainPattern;
    bool defaultDefined,haultFlag;
    double taperingRatio;
    string fileext;

    vector<LoadStep> myLoadSteps;


    ///PRIVATE MEMBER FUNCTIONS:
    //parsing
    void fillChainVector(string hw); //after decode, fill vector
    void applyTapering(string hw); //after it is filled, apply tapering
    void fillInMissingGrainInfo(string hw);//calculate the volume / mass for each store this in massvector!!!
    void calculatePowerLawAndPrefactors(string hw);//calculate the power laws

    void fillInMissingLoadStepInfo(string hw);//make sure the movables are correct for load steps and output grains?


    void splitStringByColon(const string &original, string &com, string &val);
    void splitStringBySymbol(const string &original, char symbol, string &com, string &val);
    void defaultFunction(string hw){cout<<"FILL ME IN ALREADY!!!!"<<endl;}
    void commentedLine(string hw){if(!logFileOnly)
                                    cout<<"Read and Ignoring Commented Line : "<<hw<<endl;}
    void setFilename(string hw);
    void setFilesToOutput(string hw);
    void setOutputPrecision(string hw);
    void setOutputInterval(string hw);
    void setOutputGrainNumbers(string hw);
    void addSubunit(string hw);
    void addDefaultSubunit(string hw);
    pair<string, double> getNameAndNumber(string namenum,char symbol);
    void decodeChainPattern(string hw);
    void silentModeOn(string hw){logFileOnly = true;}
    void setTapering(string hw);
    void setRestitutionFactor(string hw);
    void addLoadStep(string hw);
    void setFileExt(string hw);
    void freezeGrain(string hw);
    void unfreezeGrain(string hw);
    void addLeftWall(string hw);
    void addRightWall(string hw);
    void addWall(string hw, bool leftwall);
    void addVelocityBasedDamping(string hw);
    void addInitialCondition(string hw);
    void addBodyAcceleration(string hw);
    void decodeGrainSpec(string hw,vector<int> &gns);
    void addMCVars(string hw);
    void addDrivingCondition(string hw);
    //for driving?

    void retrieveComsAndVals(string hw, vector<pair<string,string> > &comvals);
    void whatComsWereSpecified(const vector<pair<string,string> > &comvalsUsed, const vector<string> &neededComs, vector<pair<bool,int> > &wasUsed);
    int getComLoc(string com, const vector<string> &neededComs, vector<pair<bool,int> > &wasUsed);

    //for output
    void runThroughOutputs( bool first);
    void grainOutput(StateVector *svptr,string preName,bool first);
    void interfaceOutput(StateVector *svptr,string preName,bool first); //might need to rethink StateVector vs. OneDimDoubleVector

    //for sovler incrementing time
    void incrementTimeExplicit(double &curTime, double &intervalRemaining, const double &stopTime,const double &trialDt, double &newDt);
    void incrementTimeAdaptive(double &curTime, double &intervalRemaining, const double &stopTime,const double &trialDt, double &newDt);
    void incrementTimeImplicit(double &curTime, double &intervalRemaining, const double &stopTime,const double &trialDt, double &newDt);

    void equilibriumHelperIterateStep(State &curState,State &nextState,LoadStep &details,double &scale);
    void applyDrivingConditions(State &curState,LoadStep &details);


    void throwError(string msg);
    void throwError(string msg,double val);
    void caclulateRemainingStateVars(LoadStep &details);
    bool alreadyThere(const vector<int> &vec,int val)
    {
        for(unsigned int j=0;j<vec.size();++j)
        if(vec[j]==val)
        return true;
        return false;
    }

};

/*
***********************************************************************
NOTES:
Back when Jan did his thesis he found that relative positions of the grains helped
reduce roundoff errors when determining the overlaps. This becomes confusing
very quickly when you deal with pre-compressions.
To deal with this, the relative locations are always relative to an undeformed chain.
Absolute positions are what they sound like.

let X^{abs}_i(t) be the absolute position of the i^th grain at time t

The absolute position of the grains in the undeformed chain at time 0 is
X^{abs_und}_i = Sum_j=0^{i-1}[ R^+_j + R^-_j + w_j ]  +  R^-_i + w_i/2

The distance between the Center of Mass of two neighboring grains is then
D_{i,i+1}(t) = X^{abs}_{i+1}(t) - X^{abs}_{i}(t)

Which if you consider the undeformed (each grain just in contact) distances:
D^{und}_{i,i+1}      = X^{abs_und}_{i+1} - X^{abs_und}_{i}
                        =   R^-_{i+1} + w_{i+1}/2 + R^+_{i} + w_i/2 + (Sum_j=0^{i-1}[ R^+_j + R^-_j + w_j ]  +  R^-_i + w_i/2)  -  (Sum_j=0^{i-1}[ R^+_j + R^-_j + w_j ]  +  R^-_i + w_i/2)
                        = w_{i+1}/2 + w_{i}/2 + R^-_{i+1} + R^+_{i}

The overlap between two neighboring grains is delta
which is the difference of half of the shape subtracting the distance between:
delta_{i,i+1}(t) = ( R^+_{i} + w_i/2 + R^-_{i+1} + w_{i+1}/2 ) - D_{i,i+1}(t)
                = D^{und}_{i,i+1} - D_{i,i+1}(t)
Thus delta is positive (+) when the grains overlap
and is (-) negative when they do not overlap.


 Now define the relative location of a grain as the difference between
 its current absolute position and its absolute position in the undeformed chain
 X^{rel}_i(t) = X^{abs}_i(t) - X^{abs_und}_i
=>   X^{abs}_i(t) = X^{rel}_i(t) + X^{abs_und}_i

SO that
D_{i,i+1}(t) = X^{abs}_{i+1}(t) - X^{abs}_{i}(t)
             =  X^{rel}_{i+1}(t) + X^{abs_und}_{i+1} - ( X^{rel}_{i}(t) + X^{abs_und}_{i}    )
             =  X^{rel}_{i+1}(t) -  X^{rel}_{i}(t)  + X^{abs_und}_{i+1} - X^{abs_und}_{i}
             =  X^{rel}_{i+1}(t) -  X^{rel}_{i}(t)  +  D^{und}_{i,i+1}

delta_{i,i+1}(t) = D^{und}_{i,i+1} - D_{i,i+1}(t)
                = D^{und}_{i,i+1} - (X^{rel}_{i+1}(t) -  X^{rel}_{i}(t)  +  D^{und}_{i,i+1})
                = D^{und}_{i,i+1} - X^{rel}_{i+1}(t) +  X^{rel}_{i}(t)  -  D^{und}_{i,i+1}
                =  - X^{rel}_{i+1}(t) +  X^{rel}_{i}(t)

A wall overlap will happen when  X^{rel}_{N-1} > 0 or X^{rel}_{0} < 0
delta_{left wall} =  -X^{rel}_{0}-WALLOFFSETLEFT  // (-) makes the delta positive in compression
delta_{right wall} = X^{rel}_{N-1}-WALLOFFSETRIGHT // (+) here keeps delta positive in compression

***********************************************************************
*/


struct State //state variables
{
    //friend class Simulation;
    StateVector relativeLocation,absolutePosition,velocity,force,acceleration;
    OneDimDoubleVector kineticEnergy;
    StateVector initialAbsolutePosition;
    OneDimDoubleVector overlap,interfaceForce,potentialEnergy;//index 0 is LeftWall/1stGrain(0th in code), index 1 is 1stGrain(0th in code)/2ndGrain(1st in code), last index N is NthGrain(N-1th in code)/RightWall

    double curtime;

    State(StateVector rl,StateVector ap,StateVector iap,StateVector vel,StateVector frc,
          StateVector acc,OneDimDoubleVector ke,OneDimDoubleVector ov,OneDimDoubleVector intf,OneDimDoubleVector pe,double t)
    {
        relativeLocation = rl;
        absolutePosition = ap;
        initialAbsolutePosition=iap;
        velocity = vel;
        force = frc;
        acceleration = acc;
        kineticEnergy=ke;
        interfaceForce=intf;
        potentialEnergy=pe;
        overlap = ov;
        curtime = t;

    }
    State(int length)
    {
        relativeLocation = StateVector(length);
        absolutePosition = StateVector(length);
        initialAbsolutePosition = StateVector(length);
        velocity = StateVector(length);
        force = StateVector(length);
        acceleration = StateVector(length);
        kineticEnergy=OneDimDoubleVector(length);

        potentialEnergy=OneDimDoubleVector(length+1); //N+1 interfaces
        overlap=OneDimDoubleVector(length+1); //N+1 interfaces
        interfaceForce=OneDimDoubleVector(length+1); //N+1 interfaces

        curtime = 0;

    }

    State& operator=(const State& a)
    {
        relativeLocation=a.relativeLocation;
        absolutePosition=a.absolutePosition;
        initialAbsolutePosition=a.initialAbsolutePosition;
        velocity=a.velocity;
        force = a.force;
        acceleration = a.acceleration;
        overlap = a.overlap;
        kineticEnergy=a.kineticEnergy;
        interfaceForce=a.interfaceForce;
        potentialEnergy=a.potentialEnergy;

        curtime = a.curtime;

    return *this;
    }
    void resize(int newSize)
    {
        relativeLocation.resize(newSize);
        absolutePosition.resize(newSize);
        initialAbsolutePosition.resize(newSize);
        velocity.resize(newSize);
        force.resize(newSize);
        acceleration.resize(newSize);
        kineticEnergy.resize(newSize);
        potentialEnergy.resize(newSize+1);
        interfaceForce.resize(newSize+1);
        overlap.resize(newSize+1);
    }
    friend std::ostream& operator<< (std::ostream& stream, const State& mys)
    {
        stream<<"This is the state at time : "<<mys.curtime<<endl;
        stream<<"\t relativeLocation : "<<mys.relativeLocation<<endl;
        stream<<"\t absolutePosition : "<<mys.absolutePosition<<endl;
        stream<<"\t velocity : "<<mys.velocity<<endl;
        stream<<"\t force : "<<mys.force<<endl;
        stream<<"\t acceleration : "<<mys.acceleration<<endl;
        OneDimDoubleVector ftp= mys.overlap;
        for(unsigned int j=0;j<mys.overlap.vec.size()-1;++j)
            stream<<"\t overlap: "<<"("<<j<<","<<j+1<<")\t :"<<ftp[j]<<endl;

        return stream;
    }
}; ///end of State structure



struct DriverFunction
{
    friend class Simulation;
    StateVector *svptr;//affects this statevector
    bool overRide; ///true -> set variable equal, else add value
    drivingFunction eval;
    vector<double> params;
    int grainNumber;

    DriverFunction()
    {
        svptr = NULL;
        overRide = false;
        eval = &meth::zero;
        grainNumber = -1;
    }
    DriverFunction& operator=(const DriverFunction& a)
    {
        svptr = a.svptr;
        overRide = a.overRide;
        eval = a.eval;
        grainNumber = a.grainNumber;
        return *this;
    }
    friend std::ostream& operator<< (std::ostream& stream, const DriverFunction& myd)
    {
        stream<<"Details regarding this driver function: \n";
        stream<<"Set or add : "<<(myd.overRide ? "set" : "add")<<endl;
        stream<<"Affects grain number : "<<myd.grainNumber<<endl;
        stream<<"parameters : ";
        for(unsigned int j=0;j<myd.params.size();++j)
        stream<<myd.params[j]<<'\t';
        stream<<endl;

        return stream;
    }
};

struct MonteCarloParams
{
    int Nmax;
    int countMax;
    double acceptRatio;
    double tolerance;
    double scale;
    double shift;
    int maxCutBacks;

    MonteCarloParams()
    {
            countMax = 100;
            Nmax = 1500;
            acceptRatio=0.4;
            tolerance =  1e-3;
            scale = 0.1;
            shift = 0.25;
            maxCutBacks=10;
    }
    MonteCarloParams& operator=(const MonteCarloParams& a)
    {
       Nmax = a.Nmax;
       countMax = a.countMax;
       acceptRatio=a.acceptRatio;
       tolerance =a.tolerance;
       scale = a.scale;
       shift = a.shift;
       maxCutBacks=a.maxCutBacks;
        return *this;
    }
     friend std::ostream& operator<< (std::ostream& stream, const MonteCarloParams& mc)
    {
        stream<<"Details regarding the Monte Carlo Params (if used): \n";
        stream<<"acceptance ratio : "<<mc.acceptRatio<<endl;
        stream<<"max iterations : "<<mc.countMax<<endl;
        stream<<"max cutbacks for scale : "<<mc.maxCutBacks<<endl;
        stream<<"MC avgs : "<<mc.Nmax<<endl;
        stream<<"tolerance  : "<<mc.tolerance<<endl;
        stream<<"random step scale : "<<mc.scale<<endl;
        stream<<"random step shift : "<<mc.shift<<endl;
        stream<<endl;

        return stream;
    }
};



///A load step is a structure that has the details about the load step!
//ie walls or no walls, equilibrium solution or dynamic
//other boundary condition or body forces/accel

struct LoadStep
{
    //friend class Simulation;
    ///need a function pointer to integration method that takes a current state to a new state
    solver mySolver;
    timeUpdate myUpdater;
    double startTime,stopTime,dt;

    double bodyAcceleration;
    double velocityBasedDamping;
    double restitutionFactor;
    ///driving function takes in a time, will affect some part of the state vector
    list<DriverFunction> myDriver;

    //OneDimDoubleVector kineticEnergy,potentialEnergy; //these are still vector of scalar

    double leftWallForcePrefactor,rightWallForcePrefactor;
    double leftWallOffset, rightWallOffset;
    double leftWallForcePower,rightWallForcePower; ///there will be forces added to the grain from teh wall if it is there

    Grain leftWallGrain, rightWallGrain;

    StateVector addThisInitialVelocity,addThisInitialDisplacement;

    vector<bool> movable; //in the load step is a grain movable?

    MonteCarloParams myMCparam;

    LoadStep()
    {
        LoadStep(0);
    }

    LoadStep(solver sv,timeUpdate myU,double st,double stp, double step, double body,double visc,double rf,
         list<DriverFunction> myd,
    //      OneDimDoubleVector ke,OneDimDoubleVector pe,
          double lwfpr,double rwfpr,double lws,double rws,double lwp,double rwp,Grain lg,Grain rg,vector<bool> mv,
          StateVector icv,StateVector icd,MonteCarloParams mcp)
    {
        mySolver = sv;
        myUpdater = myU;
        startTime=st;
        stopTime=stp;
        dt = step;
           bodyAcceleration=body;
           velocityBasedDamping=visc;
           restitutionFactor=rf;
           myDriver=myd;
        myMCparam=mcp;
        addThisInitialDisplacement=icd;
        addThisInitialVelocity=icv;
//        kineticEnergy = ke;
  //      potentialEnergy = pe;

        leftWallForcePrefactor = lwfpr;rightWallForcePrefactor=rwfpr;
        leftWallOffset=lws;rightWallOffset=rws;
        leftWallForcePower=lwp;rightWallForcePower=rwp;
        leftWallGrain=lg;rightWallGrain=rg;
        movable = mv;
    }
    LoadStep(int length)
    {
//        kineticEnergy = OneDimDoubleVector(length);
  //      potentialEnergy = OneDimDoubleVector(length);
        bodyAcceleration=0;
        velocityBasedDamping=0;
        restitutionFactor=1;//force unload = force duringload
        startTime=0;
        stopTime=0;
        dt = 1;
        mySolver = &Simulation::velocityVerletStep;
        myUpdater = &Simulation::incrementTimeExplicit;

        leftWallForcePrefactor = 0;rightWallForcePrefactor=0;
        leftWallOffset=0;rightWallOffset=0;
        leftWallForcePower=1;rightWallForcePower=1;
        leftWallGrain=Grain();rightWallGrain=Grain();
        for( int j=0;j<length;++j)
            movable.push_back(true);
        addThisInitialDisplacement.resize(length);
        addThisInitialVelocity.resize(length);

        myMCparam=MonteCarloParams();
    }

    LoadStep& operator=(const LoadStep& a)
    {
        startTime=a.startTime;
        stopTime=a.stopTime;
        dt=a.dt;
        velocityBasedDamping=a.velocityBasedDamping;
        restitutionFactor=a.restitutionFactor;
        mySolver = a.mySolver;
        myDriver = a.myDriver;
        myUpdater = a.myUpdater;
        myMCparam = a.myMCparam;
        bodyAcceleration=a.bodyAcceleration;
        addThisInitialDisplacement=a.addThisInitialDisplacement;
        addThisInitialVelocity=a.addThisInitialVelocity;
//        kineticEnergy = a.kineticEnergy;
  //      potentialEnergy=a.potentialEnergy;

        leftWallForcePrefactor = a.leftWallForcePrefactor;rightWallForcePrefactor=a.rightWallForcePrefactor;
        leftWallOffset=a.leftWallOffset;rightWallOffset=a.rightWallOffset;
        leftWallForcePower=a.leftWallForcePower;rightWallForcePower=a.rightWallForcePower;
        leftWallGrain=a.leftWallGrain;rightWallGrain=a.rightWallGrain;
        movable = a.movable;

    return *this;
    }
    friend std::ostream& operator<< (std::ostream& stream, const LoadStep& mys)
    {
        stream<<"Details about this LoadStep  : "<<endl;
        stream<<"\t start time : "<<mys.startTime<<endl;
        stream<<"\t stop time : "<<mys.stopTime<<endl;
        stream<<"\t time step: "<<mys.dt<<endl;
        stream<<"\t body acceleration : "<<mys.bodyAcceleration<<endl;
        stream<<"\t velocity based damping  : "<<mys.velocityBasedDamping<<endl;
        stream<<"\t force restitution (unload force is reduced by this factor)  : "<<mys.restitutionFactor<<endl;
        //stream<<"\t ADD something about the DRIVER and SOLVER"<<endl;
        list<DriverFunction>::iterator it;
        list<DriverFunction> myd = mys.myDriver;


        for(it = myd.begin();it!=myd.end();++it)
        stream<<*it<<endl;

//        stream<<"\t kinetic energy : "<<mys.kineticEnergy<<endl;
  //      stream<<"\t potential energy : "<<mys.potentialEnergy<<endl;
        stream<<"Wall info: "<<endl;
        stream<<"\t LEFT Wall : "<<(mys.leftWallForcePrefactor > 0?"on":"off")<<endl;
        if(mys.leftWallForcePrefactor > 0)
        stream<<"\t LEFT Wall Offset: "<<mys.leftWallOffset<<endl<<"\t LEFT Wall Force Prefactor: "<<mys.leftWallForcePrefactor<<endl<<"\t LEFT Wall Force Power: "<<mys.leftWallForcePower<<endl;

        stream<<"\t RIGHT Wall : "<<(mys.rightWallForcePrefactor > 0?"on":"off")<<endl;
        if(mys.rightWallForcePrefactor > 0)
        stream<<"\t RIGHT Wall Offset: "<<mys.rightWallOffset<<endl<<"\t RIGHT Wall Force Prefactor: "<<mys.rightWallForcePrefactor<<endl<<"\t RIGHT Wall Force Power: "<<mys.rightWallForcePower<<endl;

        stream<<"The following grains are movable ? "<<endl;
        for(unsigned int j=0;j<mys.movable.size();++j)
        stream<<"Grain : "<<j+1<<" is movable : "<<(mys.movable[j]?"yes":"no")<<endl;

        stream<<"At the beginning of the load step I will add this vector to the velocity:"<<endl;
        stream<<mys.addThisInitialVelocity<<endl;

        stream<<"At the beginning of the load step I will add this vector to the displacement:"<<endl;
        stream<<mys.addThisInitialDisplacement<<endl;

         stream<<mys.myMCparam<<endl;

        return stream;
    }

};///End of LoadSTep Structure



#endif // SIMULATION_H
