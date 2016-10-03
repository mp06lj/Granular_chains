#include "Simulation.h"
#include <locale>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <queue>

#define FIRST_PRIORITY 5
#define SECOND_PRIORITY 4
#define THIRD_PRIORITY 3
#define FOURTH_PRIORITY 2
#define FIFTH_PRIORITY 1

#define HIGH_PRIORITY 10
#define MEDIUM_HIGH_PRIORITY 9
#define MEDIUM_PRIORITY 8
#define MEDIUM_LOW_PRIORITY 7
#define LOW_PRIORITY 6

using namespace std;
//Adam Sokolow 1-31-14
// Rewrite of taperchain code
// SImulation class handles a lot of the over arching items

///Simulation Class Functions, specific to initialization, copy, destructor, and the most general types.
///Also has the stream operations

Simulation::Simulation()
{
    //ctor
    //these are above all else, need to happen first!
    myParser["loadstep"]=pair<int,command>(HIGH_PRIORITY,&Simulation::addLoadStep);
    myParser["defaultgrain"]=pair<int,command>(HIGH_PRIORITY,&Simulation::addDefaultSubunit);
    myParser["silentmode"]=pair<int,command>(HIGH_PRIORITY,&Simulation::silentModeOn);

    myParser["subunit"]=pair<int,command>(MEDIUM_PRIORITY,&Simulation::addSubunit);
    myParser["montecarloparams"]=pair<int,command>(MEDIUM_PRIORITY,&Simulation::addMCVars);

    //Fill the commands with references to the appropriate void functions declared above
    myParser["filename"]=pair<int,command>(LOW_PRIORITY,&Simulation::setFilename);
    myParser["outputfileextension"]=pair<int,command>(LOW_PRIORITY,&Simulation::setFileExt);
    myParser["outputfiles"]=pair<int,command>(LOW_PRIORITY,&Simulation::setFilesToOutput);
    myParser["outputprecision"]=pair<int,command>(LOW_PRIORITY,&Simulation::setOutputPrecision);
    myParser["outputinterval"]=pair<int,command>(LOW_PRIORITY,&Simulation::setOutputInterval);
    myParser["outputgrains"]=pair<int,command>(LOW_PRIORITY,&Simulation::setOutputGrainNumbers);


    myParser["chainpattern"]=pair<int,command>(LOW_PRIORITY,&Simulation::decodeChainPattern);
    myParser["taperingratio"]=pair<int,command>(HIGH_PRIORITY,&Simulation::setTapering);


    myParser["freezegrain"]=pair<int,command>(LOW_PRIORITY,&Simulation::freezeGrain);
    myParser["unfreezegrain"]=pair<int,command>(LOW_PRIORITY,&Simulation::unfreezeGrain);

    myParser["bodyacceleration"]=pair<int,command>(LOW_PRIORITY,&Simulation::addBodyAcceleration);
    myParser["velocitybaseddamping"]=pair<int,command>(LOW_PRIORITY,&Simulation::addVelocityBasedDamping);
    myParser["unloadforceratio"]=pair<int,command>(LOW_PRIORITY,&Simulation::setRestitutionFactor);

    myParser["leftwall"]=pair<int,command>(LOW_PRIORITY,&Simulation::addLeftWall);
    myParser["rightwall"]=pair<int,command>(LOW_PRIORITY,&Simulation::addRightWall);

    myParser["initialcondition"]=pair<int,command>(THIRD_PRIORITY,&Simulation::addInitialCondition);

    myParser["drivingfunction"]=pair<int,command>(FOURTH_PRIORITY,&Simulation::addDrivingCondition);
    myParser["c"]=pair<int,command>(LOW_PRIORITY,&Simulation::commentedLine);
    myParser["end"]=pair<int,command>(LOW_PRIORITY,&Simulation::commentedLine);
    myParser["end:"]=pair<int,command>(LOW_PRIORITY,&Simulation::commentedLine);

    ///These commands are added in by the parser to ensure that commands get called in the proper order
    myParser["fillchain"]=pair<int,command>(FIRST_PRIORITY,&Simulation::fillChainVector);
    myParser["applytaper"]=pair<int,command>(SECOND_PRIORITY,&Simulation::applyTapering);
    myParser["fillinmissing"]=pair<int,command>(THIRD_PRIORITY,&Simulation::fillInMissingGrainInfo);
    myParser["calcpreandpow"]=pair<int,command>(FOURTH_PRIORITY,&Simulation::calculatePowerLawAndPrefactors);
    myParser["fillinmissingparttwo"]=pair<int,command>(FIFTH_PRIORITY,&Simulation::fillInMissingLoadStepInfo);
    ///Do not change these above!!




    oldCommands["n"]="Use a chain pattern instead! This will avoid confusion in the long run.";
    oldCommands["files"]="Use the command \"output files\" instead, this is to clarify these are the output files (A)ll, (O)verlaps, (I)nterface force, (F)orce on grains, (X)relative location, a(B)solute positions,(V)elocity, (K)inetic energy grain, (P)otential energy interfaces.";
    oldCommands["nsteps"]="It\'s too confusing to specify this and worry about maxtime changing. Just use \"dt\" and \"maxtime\" to specify these";
    oldCommands["largeinitv"]="This command is gone. Replaced by some initial condition in the load step command";
    oldCommands["smallinitv"]="This command is gone. Replaced by some initial condition in the load step command";
    oldCommands["precision"]="This command is replaced by \"output precision\" to reduce confusion with the simulation precision";
    oldCommands["timemin"]="This command is gone. ";
    oldCommands["timemax"]="This command is gone. ";
    oldCommands["timespits"]="This command is replaced by \"output interval\". Sounds a little more professional...";
    oldCommands["rlarge"]="This command is replaced by the chain patterning.";
    oldCommands["d"]="This command is replaced by the chain patterning.";
    oldCommands["rho"]="This command is replaced by the chain patterning.";
    oldCommands["wall"]="This command is replaced by some boundary condition command.";
    oldCommands["preload"]="This command is replaced by some load step command.";
    oldCommands["exponential"]="This command is replaced by the chain patterning and shape commands (use alpha).";
    oldCommands["q"]="This command is replaced by \"tapering ratio\" and applies to the whole chain on top of a pattern. Using a (-) will go from right to left instead of left to right";
    oldCommands["w"]="This command is replaced by \"unload force ratio\".";
     oldCommands["dt"]="This command is replaced by \"loadstep\" definition. See manual.";
    oldCommands["grains"]="This command is replaced by \"output grains\". See manual.";
    myCurrentState = new State(0);
    myPrevState = NULL;
    baseFileName = "baseName";
    outputInterval = 1e-3;
    outputPrecision = 6;
    taperingRatio=1;
    fileext = ".dat";
    haultFlag=false;
    defaultDefined=false;
    forceOutput=false;
    relLocationOutput=false;
    velocityOutput=false;
    kineticEnergyOutput=false;
    potentialEnergyOutput=false;
    logFileOnly=false;
    overlapOutput=false;
    interfaceForceOutput=false;
    absoluteOutput=false;

}

Simulation::~Simulation()
{
    //dtor
}

Simulation::Simulation(const Simulation& other)
{
    //copy ctor
    myParser = other.myParser ;
    oldCommands = other.oldCommands;
    myCurrentState = other.myCurrentState;
    myPrevState = other.myPrevState;
    commandList=other.commandList;
    myLoadSteps = other.myLoadSteps;

    fileext=other.fileext;
    baseFileName=other.baseFileName;
    myChain=other.myChain;
    massVector=other.massVector;
    forceOutput=other.forceOutput;
    relLocationOutput=other.relLocationOutput;
    velocityOutput=other.velocityOutput;
    kineticEnergyOutput=other.kineticEnergyOutput;
    potentialEnergyOutput=other.potentialEnergyOutput;
    overlapOutput=other.overlapOutput;
    interfaceForceOutput=other.interfaceForceOutput;
    absoluteOutput=other.absoluteOutput;

    forcePrefactors = other.forcePrefactors;//given two grain indices, what is their prefactor?
    forcePowerLaw = other.forcePowerLaw;

    //Output related
    outputPrecision= other.outputPrecision;
    outputInterval = other.outputInterval;

    logFileOnly=other.logFileOnly;

    //Pattern related
    outputGrainNumbers=other.outputGrainNumbers;//these run from 1 -> N, but internally are 0 -> N-1
    outputInterfaceNumbers=other.outputInterfaceNumbers;
    unitDefinitions = other.unitDefinitions;
    chainPattern = other.chainPattern;
    defaultDefined=other.defaultDefined;
    taperingRatio = other.taperingRatio;

    haultFlag=other.haultFlag;

}

Simulation& Simulation::operator=(const Simulation& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    myParser = rhs.myParser;
    oldCommands = rhs.oldCommands;
    myCurrentState = rhs.myCurrentState;
    myPrevState = rhs.myPrevState;//do I want to copy the contents or just the pointer?
    myLoadSteps = rhs.myLoadSteps;
    fileext=rhs.fileext;
    commandList = rhs.commandList;
    baseFileName=rhs.baseFileName;
    myChain=rhs.myChain;
    massVector = rhs.massVector;
    forceOutput=rhs.forceOutput;
    relLocationOutput=rhs.relLocationOutput;
    velocityOutput=rhs.velocityOutput;
    kineticEnergyOutput=rhs.kineticEnergyOutput;
    potentialEnergyOutput=rhs.potentialEnergyOutput;
    overlapOutput=rhs.overlapOutput;
    interfaceForceOutput=rhs.interfaceForceOutput;
    absoluteOutput=rhs.absoluteOutput;

    forcePrefactors = rhs.forcePrefactors;//given two grain indices, what is their prefactor?
    forcePowerLaw = rhs.forcePowerLaw;

    cout<<"REMMVER TO COPY STATE and LOADSTEPS"<<endl;

    //Output related
    outputPrecision= rhs.outputPrecision;
    outputInterval = rhs.outputInterval;

    haultFlag=rhs.haultFlag;

    logFileOnly=rhs.logFileOnly;

    //Pattern related
    outputGrainNumbers=rhs.outputGrainNumbers;//these run from 1 -> N, but internally are 0 -> N-1
    outputInterfaceNumbers=rhs.outputInterfaceNumbers;
    unitDefinitions = rhs.unitDefinitions;
    chainPattern = rhs.chainPattern;
    defaultDefined=rhs.defaultDefined;
    taperingRatio = rhs.taperingRatio;

    return *this;
}

std::ostream& operator<< (std::ostream& stream, const vector<Grain> &ch)
    {
        for(unsigned int lcv =0; lcv < ch.size();++lcv)
        {
            stream<<"Grain "<<lcv+1<<std::endl;//Keep the user in the dark about N-1 issues
            stream<<"\t"<<ch[lcv]<<endl;
        }
        return stream;
    }
std::ostream& operator<< (std::ostream& stream, const Simulation &sim)
    { //NEED TO UPDATE THIS ********************
        stream<<"*******************************************"<<endl;
        stream<<"So I have read your input file and now find the following:"<<endl;
        stream<<"*******************************************"<<endl;
        stream<<"You want outputs to have the base name : "<<sim.baseFileName<<endl;
        stream<<"You want these output files: "<<endl;
        stream<<"\tForce output : "<<(sim.forceOutput ? " yes " : " no ")<<endl;
        stream<<"\tRelative Location output : "<<(sim.relLocationOutput ? " yes " : " no ")<<endl;
        stream<<"\tAbsolute Location output : "<<(sim.absoluteOutput ? " yes " : " no ")<<endl;
        stream<<"\tVelocity output : "<<(sim.velocityOutput ? " yes " : " no ")<<endl;
        stream<<"\tPotential Energy output : "<<(sim.potentialEnergyOutput ? " yes " : " no ")<<endl;
        stream<<"\tKinetic Energy output : "<<(sim.kineticEnergyOutput ? " yes " : " no ")<<endl;
        stream<<"\tOverlaps output : "<<(sim.overlapOutput ? " yes " : " no ")<<endl;
        stream<<"\tInterface forces output : "<<(sim.interfaceForceOutput ? " yes " : " no ")<<endl;
        stream<<"And you want these outputs for the following grains (if they exist): "<<endl<<'\t';
        for(unsigned int j=0; j<sim.outputGrainNumbers.size()-1;++j)
            stream<<sim.outputGrainNumbers[j]<<", ";
        stream<<sim.outputGrainNumbers[sim.outputGrainNumbers.size()-1];
        stream<<endl<<endl<<"These were all the subunits you defined: "<<endl;
        std::map<string,Grain>::iterator it;
        map<string,Grain> tp = sim.unitDefinitions;

        for( it=tp.begin(); it!=tp.end(); ++it)
        stream<<"\tUnit Name: "<<it->first<<endl<<it->second<<endl;

        stream<<endl<<"I read in a pattern which I decoded to be: "<<endl<<"\t";
        for( unsigned int cp=0; cp<sim.chainPattern.size();++cp)
                stream<<sim.chainPattern[cp]<<"|";
        stream<<endl;

        if(sim.taperingRatio!=1)
        stream<<endl<<"I also used a tapering ratio on top of the pattern with ratio : "<<sim.taperingRatio<<endl;

        stream<<endl<<"So that all in all these are the grains I have added to the chain: "<<endl;
        stream<<"From left to right (going in the positive direction)"<<endl<<endl;
        stream<<sim.myChain<<endl<<endl;

        stream<<endl<<"Now let's look at the interfaces between grains."<<endl;

        map<pair<int,int>,double> ftp = sim.forcePrefactors;
        map<pair<int,int>,double> ftp2 = sim.forcePowerLaw;

        for(unsigned int j=0;j<sim.myChain.size()-1;++j)
            stream<<"\tInterface: "<<"("<<j+1<<","<<j+2<<")\t has a force prefactor :"<<ftp[pair<int,int>(j,j+1)]<<" and force power law : "<<ftp2[pair<int,int>(j,j+1)]<<endl;
        stream<<"Wall interfaces (if used) will be treated as a special case"<<endl;


        Grain gs;
        double undeformedLength = 0;
        for(unsigned int j=0;j<sim.myChain.size();++j)
        {
            gs = sim.myChain[j];
            undeformedLength+=gs.GetmyShape().r_minus+gs.GetmyShape().width+gs.GetmyShape().r_plus;
        }
        stream<<"The undeformed length of the chain before any load steps is : "<<undeformedLength<<endl;

        for(unsigned int j=0;j<sim.myLoadSteps.size();++j)
            stream<<sim.myLoadSteps[j]<<endl;


        return stream;
    }
///Print help!!!!
void Simulation::printHelp()
{
    cout<<" Taper Chain Simulation Help "<<endl;
    //run through hash table and print commands.
    cout<<" To run this code, you need to provide a file name!!!!!"<<endl;
    cout<<" I used to assume the parameters were in a parameters.txt file, but now"<<endl;
    cout<<" you can decide what file to use. So run this program from the command line, or"<<endl;
    cout<<" drag and drop a text file that has the parameters onto the executable. "<<endl;
    cout<<" command line example >> Simulation.exe parameters.txt "<<endl;
    cout<<" Some commands have changed since 2004. "<<endl;
    cout<<endl<<" I will list the commands I know now: "<<endl;

    for (std::map<string,pair<int,command> >::iterator it=myParser.begin(); it!=myParser.end(); ++it)
        cout << ">"<<it->first <<endl<<endl;

    cout<<endl<<" Out of the kindness of my heart I will list the old commands and what to use intead: "<<endl;

    for (std::map<string,string>::iterator it=oldCommands.begin(); it!=oldCommands.end(); ++it)
        cout << ">"<<it->first <<" : "<<it->second<<endl<<endl;

}
void Simulation::runThroughOutputs(bool first)
{
    for(unsigned int i=0;i<myOutputs.size();++i)
        (this->*myOutputs[i].myOutputter)(myOutputs[i].svptr,myOutputs[i].filename,first);
}


void Simulation::grainOutput(StateVector *svptr,string preName,bool first)
{
    //Have bools regarding which files to output
    //Have list of grains to output
    //Have current state.
    ofstream out_file;

    if(first)
    {
        out_file.open((preName+baseFileName+fileext).c_str());
        out_file<<"Time\t";
        for(unsigned int j=0;j<outputGrainNumbers.size()-1;++j)
            out_file<<"Grain"<<outputGrainNumbers[j]<<'\t';
        out_file<<"Grain"<<outputGrainNumbers[outputGrainNumbers.size()-1]<<endl;
        out_file.close();
    }
    else
    {
    out_file.open((preName+baseFileName+fileext).c_str(),ios::app);
    out_file.precision(outputPrecision);
     out_file<<myCurrentState->curtime<<'\t';

    for(unsigned int j=0;j<outputGrainNumbers.size()-1;++j)
    out_file<<svptr->at(outputGrainNumbers[j]-1)<<'\t';
    out_file<<svptr->at(outputGrainNumbers[outputGrainNumbers.size()-1]-1)<<endl;
    }
    out_file.close();

}
void Simulation::interfaceOutput(StateVector *svptr,string preName,bool first)
{
    //use the outputGrainNumbers to get the interfaces surrounding that grain
    //output them from the StateVector
    //how to handle walls?
//Have bools regarding which files to output
    //Have list of grains to output
    //Have current state.
    ofstream out_file;

    if(first)
    {   //user thinks of grains 1->N, but code uses 0->N-1
        out_file.open((preName+baseFileName+fileext).c_str());
        out_file<<"Time\t";
        for(unsigned int j=0;j<outputInterfaceNumbers.size()-1;++j)
        {
            if(outputInterfaceNumbers[j]==0) //left wall/1st grain
                out_file<<"(LW,G1)\t";
            else //default
            out_file<<"(G"<<outputInterfaceNumbers[j]<<",G"<<outputInterfaceNumbers[j]+1<<")\t"; //interface to the left grain
        }
        if(outputInterfaceNumbers[outputInterfaceNumbers.size()-1]==(int)myChain.size()) //last grain/rightwall
                out_file<<"(G"<<myChain.size()<<",RW)";
        else
            out_file<<"(G"<<outputInterfaceNumbers[outputInterfaceNumbers.size()-1]<<",G"<<outputInterfaceNumbers[outputInterfaceNumbers.size()-1]+1<<")"; //interface to the left grain
        out_file<<endl;

    }
    else
    {
    out_file.open((preName+baseFileName+fileext).c_str(),ios::app);
    out_file.precision(outputPrecision);

      out_file<<myCurrentState->curtime<<'\t';
        for(unsigned int j=0;j<outputInterfaceNumbers.size()-1;++j)
            out_file<<svptr->at(outputInterfaceNumbers[j])<<"\t";
        out_file<<svptr->at(outputInterfaceNumbers[outputInterfaceNumbers.size()-1]);
        out_file<<endl;
    }
    out_file.close();

}
///*************************************************************************************
///Now mostly functions related to the parsing of input
/// and the functions associated with the command pointers in the hash table
void Simulation::throwError(string msg)
{
    haultFlag=true;
    cout<<"\n**************************************************************************\n";
    cout<<msg;
    cout<<"\n**************************************************************************\n";
}
void Simulation::throwError(string msg,double val)
{
    haultFlag=true;
    cout<<"\n**************************************************************************\n";
    cout<<msg<<val;
    cout<<"\n**************************************************************************\n";
}
void Simulation::retrieveComsAndVals(string hw, vector<pair<string,string> > &comvals)
{ //assume that all commands appear iwth a value, separated like this: | com1 # val1 | com2 # val2 | etc

string original = hw;
string leftover,val,com;
do{
    splitStringBySymbol(original,'#',com,leftover);
    original=leftover;
    splitStringBySymbol(original,'|',val,leftover);
    original=leftover;

    comvals.push_back(pair<string,string>(com,val));
    }while(leftover.size()>0);

}
void Simulation::whatComsWereSpecified(const vector<pair<string,string> > &comvalsUsed, const vector<string> &neededComs, vector<pair<bool,int> > &wasUsed)
{
    bool found = false;
    unsigned int i;
    for(unsigned int j=0; j<neededComs.size();++j)
    {
        found=false;
        for( i=0; i<comvalsUsed.size();++i)
            if(neededComs[j] == comvalsUsed[i].first)
             {found=true;break;}
        wasUsed.push_back(pair<bool,int>(found,i));
    }
    //check to see if there are any commands that are not understood
    for(unsigned int j=0; j<comvalsUsed.size();++j)
    {
        found = false;
        for( i =0; i< neededComs.size();++i)
        if(comvalsUsed[j].first == neededComs[i])
        {found=true; break;}
        if(!found)
        throwError("I dont know this command you have used: \""+comvalsUsed[j].first+"\" did you spell it correctly?");
    }
}
int Simulation::getComLoc(string com, const vector<string> &neededComs, vector<pair<bool,int> > &wasUsed)
{
    for(unsigned int j=0;j<neededComs.size();++j)
    if(com == neededComs[j])
    if(wasUsed[j].first)
    return wasUsed[j].second;

    return -1;
}
pair<bool,string> Simulation::parseString(string original)
{
    string myCommand;
    string myVal;
    splitStringByColon(original,myCommand,myVal);

    if(myParser.find(toLowerCase(removeWhiteSpace(myCommand)))!=myParser.end())
   {
        commandList.push(pair<int,string>(myParser[toLowerCase(removeWhiteSpace(myCommand))].first,original));
    // (this->*myParser[toLowerCase(removeWhiteSpace(myCommand))].second)(removeWhiteSpace(myVal));
     return pair<bool,string>(true,"All Good");
   }
    else if(oldCommands.find(toLowerCase(removeWhiteSpace(myCommand)))!=oldCommands.end())
   {
    return pair<bool,string>(false,oldCommands[toLowerCase(removeWhiteSpace(myCommand))]);
   }
   else if(removeWhiteSpace(myCommand).size()==0)
   return pair<bool,string>(true,("Blank line"));

   return pair<bool,string>(false,("Command Not Found, check spelling on  : "+myCommand));

}
void Simulation::executeCommands()
{
    pair<int,string> tp;
    string myCommand;
    string myVal;
    string original;
    vector<string> addCommands;
    addCommands.push_back("fillchain");addCommands.push_back("applytaper");addCommands.push_back("fillinmissing");
    addCommands.push_back("calcpreandpow");addCommands.push_back("fillinmissingparttwo");

    for(unsigned int j=0;j<addCommands.size();++j)
    {
        myCommand = addCommands[j];
        original = myCommand + ": blank";
        commandList.push(pair<int,string>(myParser[toLowerCase(removeWhiteSpace(myCommand))].first,original));
    }

    while(!commandList.empty())
    {
        if(haultFlag)
        break;
        tp = commandList.top();
        commandList.pop();
        original=tp.second;
        //cout<<"Priority: "<<tp.first<<"\t"<<tp.second<<endl;
        splitStringByColon(original,myCommand,myVal);
        (this->*myParser[toLowerCase(removeWhiteSpace(myCommand))].second)(removeWhiteSpace(myVal));
    }
}
string Simulation::toLowerCase(string com)
{
string tmp="";
locale loc;//not really sure what locale does, but examples show this to work...
for (string::size_type i=0; i<com.length(); ++i)
    tmp+=tolower(com[i],loc);
    return tmp;
}
string Simulation::removeWhiteSpace(string com)
{ //removes spaces, tabs, and endlines
    string tmp="";
    for(unsigned int lcv=0; lcv<com.length();++lcv)
    if(com[lcv]!=' ' && com[lcv]!='\n' && com[lcv]!='\t')
    tmp+=com[lcv];
    return tmp;
}
void Simulation::splitStringByColon(const string &original, string &com, string &val)
{  //what it sounds like, splits the input string by the FIRST colon it comes across
    com="";
    val ="";
    bool beforeColon=true;
    for(unsigned int lcv = 0; lcv< original.length();++lcv)
    {
        if(original[lcv]==':' && beforeColon)
            beforeColon = false;
            else (beforeColon ? com+=original[lcv] : val+=original[lcv]);
    }
}

void Simulation::splitStringBySymbol(const string &original,char symbol, string &com, string &val)
{  //what it sounds like, splits the input string by the FIRST colon it comes across
    com="";
    val ="";
    bool before=true;
    for(unsigned int lcv = 0; lcv< original.length();++lcv)
    {
        if(original[lcv]==symbol && before)
            before = false;
            else (before ? com+=original[lcv] : val+=original[lcv]);
    }
  //  cout<<"original: "<<original<<" command: "<<com<<" val: "<<val<<endl;
}


/// THese are the functions that are in the hash table
void Simulation::setFilename(string hw)
{
    baseFileName=hw;
}
void Simulation::setFileExt(string hw)
{
    fileext=hw;
}

void Simulation::setFilesToOutput(string hw)
{
locale loc;
outputInfo ti;

    for(unsigned int lcv =0; lcv< hw.size();++lcv)
    {
        switch(tolower(hw[lcv],loc))
    {
        case 'x':relLocationOutput=true;break;
        case 'f':forceOutput=true;break;
        case 'v':velocityOutput=true;break;
        case 'b':absoluteOutput=true;break;
        case 'k':kineticEnergyOutput=true;break;
        case 'o':overlapOutput=true;break;
        case 'i':interfaceForceOutput=true;break;
        case 'p':potentialEnergyOutput=true;break;
        case 'a':{relLocationOutput=true;forceOutput=true;velocityOutput=true;potentialEnergyOutput=true;
                kineticEnergyOutput=true;overlapOutput=true;interfaceForceOutput=true;absoluteOutput=true;break;}
        default : break;
    }}
    if(forceOutput)
    {
        ti.filename = "forceOnGrain";
            ti.svptr=&(myCurrentState->force);
            ti.myOutputter=&Simulation::grainOutput;
            myOutputs.push_back(ti);
    }
    if(relLocationOutput)
    {
        ti.filename = "relLocationGrain";
            ti.svptr=&(myCurrentState->relativeLocation);
            ti.myOutputter=&Simulation::grainOutput;
            myOutputs.push_back(ti);
    }
    if(absoluteOutput)
    {
            ti.filename = "absLocationOfGrain";
            ti.svptr=&(myCurrentState->absolutePosition);
            ti.myOutputter=&Simulation::grainOutput;
            myOutputs.push_back(ti);
    }
    if(velocityOutput)
    {
            ti.filename = "velocityOfGrain";
            ti.svptr=&(myCurrentState->velocity);
            ti.myOutputter=&Simulation::grainOutput;
            myOutputs.push_back(ti);
    }
    if(kineticEnergyOutput)
    {
            ti.filename = "KEOfGrain";
            ti.svptr=&(myCurrentState->kineticEnergy);
            ti.myOutputter=&Simulation::grainOutput;
            myOutputs.push_back(ti);
    }
    if(potentialEnergyOutput)
    {
            ti.filename = "PEOfInterface";
            ti.svptr=&(myCurrentState->potentialEnergy);
            ti.myOutputter=&Simulation::interfaceOutput;
            myOutputs.push_back(ti);
    }
    if(overlapOutput)
    {
            ti.filename = "Overlaps";
            ti.svptr=&(myCurrentState->overlap);
            ti.myOutputter=&Simulation::interfaceOutput;
            myOutputs.push_back(ti);
    }
    if(interfaceForceOutput)
    {
            ti.filename = "ForceAtInterface";
            ti.svptr=&(myCurrentState->interfaceForce);
            ti.myOutputter=&Simulation::interfaceOutput;
            myOutputs.push_back(ti);
    }
}
void Simulation::setOutputPrecision(string hw)
{
    outputPrecision = atoi(hw.c_str());
}
void Simulation::setOutputInterval(string hw)
{
    outputInterval = atof(hw.c_str());
}
void Simulation::setOutputGrainNumbers(string hw)
{
//parse thw string numbers get pushed back, commas go to next number
//dash index multiple grains

//assumes the whitespace was removed
    string tmp="";
    int first=0,second=0;
    bool dash=false;
    bool onFirst = true;
    hw+=",";
    for(unsigned int lcv =0; lcv<hw.length();++lcv)
    {
        if(hw[lcv]==',' && !dash) //comma after a number, but no dash command
        {
            outputGrainNumbers.push_back(atoi(tmp.c_str()));
            tmp="";
        }
        else if(hw[lcv]==',' && dash)
        { //means we were building the second and already had the first.
            second = atoi(tmp.c_str());
            for(int j=first;j<=second;++j)
                outputGrainNumbers.push_back(j);
            tmp="";
            onFirst = true;
            dash = false;
        }
        else if(hw[lcv]=='-')
        {
            dash=true;
            onFirst = false;
            first = atoi(tmp.c_str());
            tmp="";
        }
        else if(onFirst)
        tmp += hw[lcv];
        else
        tmp += hw[lcv];
    }

}

void Simulation::freezeGrain(string hw)
{
     if(myLoadSteps.size()==0)
    {
        throwError("You need to declare at least one load step!!!");return;
    }

 string tmp="";
 string tmp2;
 string original;
    int loadstepnum;
    vector<int> gns;

    splitStringBySymbol(hw,'#',tmp,tmp2);//removes the loadstep command
    splitStringBySymbol(tmp2,'|',tmp,original);//the number is now in tmp2
    loadstepnum = atoi(tmp2.c_str());

    splitStringBySymbol(original,':',tmp,tmp2);
    original=tmp2;
   // cout<<"REad a load step "<<loadstepnum<<endl;
    decodeGrainSpec(original,gns);

//gns has the grains that are frozen
    for(unsigned int j=0;j<gns.size();++j)
    {
            while((int)myLoadSteps[loadstepnum-1].movable.size()<gns[j])
            myLoadSteps[loadstepnum-1].movable.push_back(true);
            myLoadSteps[loadstepnum-1].movable[gns[j]-1]=false;
    }

}

void Simulation::unfreezeGrain(string hw)
{
     if(myLoadSteps.size()==0)
    {
        throwError("You need to declare at least one load step!!!");return;
    }

 string tmp="";
 string tmp2;
 string original;
    int loadstepnum;
    vector<int> gns;

    splitStringBySymbol(hw,'#',tmp,tmp2);//removes the loadstep command
    splitStringBySymbol(tmp2,'|',tmp,original);//the number is now in tmp2
    loadstepnum = atoi(tmp2.c_str());

    splitStringBySymbol(original,':',tmp,tmp2);
    original=tmp2;
   // cout<<"REad a load step "<<loadstepnum<<endl;
    decodeGrainSpec(original,gns);
    //gns has the grains that are frozen
    for(unsigned int j=0;j<gns.size();++j)
    {
            while((int)myLoadSteps[loadstepnum-1].movable.size()<gns[j])
            myLoadSteps[loadstepnum-1].movable.push_back(true);
            myLoadSteps[loadstepnum-1].movable[gns[j]-1]=true;
    }

}
void Simulation::addVelocityBasedDamping(string hw)
{    if(myLoadSteps.size()==0)
    {
        throwError("You need to declare at least one load step!!!");return;
    }

    //loadstep # 2 | const # 4
    string original = toLowerCase(hw);
    string com, val,leftover;
    int loadstepval = -1;
    double damp = 0;

    do{
    splitStringBySymbol(original,'#',com,leftover);
    original = leftover;
    splitStringBySymbol(original,'|',val,leftover);
    original = leftover;

    if(com == "loadstep")
        loadstepval = atoi(val.c_str());
    else if(com == "const")
        damp = atof(val.c_str());
    else
    throwError("Error, I don't know : \""+com+"\". Check spelling perhaps.");

    }while(leftover.size()>0);

    myLoadSteps[loadstepval-1].velocityBasedDamping = damp;
}

void Simulation::decodeGrainSpec(string hw,vector<int> &gns)
{
    string tmp, tmp2, original;
    int first=0, second=0;
    bool dash=false;
    bool onFirst = true;
    tmp="";
    original= hw;
    original+=",";

    for(unsigned int lcv =0; lcv<original.length();++lcv)
    {
        if(original[lcv]==',' && !dash) //comma after a number, but no dash command
        {
            gns.push_back(atoi(tmp.c_str()));
            tmp="";
        }
        else if(original[lcv]==',' && dash)
        { //means we were building the second and already had the first.
            second = atoi(tmp.c_str());
            for(int j=first;j<=second;++j)
                gns.push_back(j);
            tmp="";
            onFirst = true;
            dash = false;
        }
        else if(original[lcv]=='-')
        {
            dash=true;
            onFirst = false;
            first = atoi(tmp.c_str());
            tmp="";
        }
        else if(onFirst)
        tmp += original[lcv];
        else
        tmp += original[lcv];
    }

}

void Simulation::addInitialCondition(string hw)
{  if(myLoadSteps.size()==0)
    {
        throwError("You need to declare at least one load step!!!");return;
    }

    //load step # 3 |
    double addval = 0;
    int loadstepval=-1;
    string variable;
    string original = toLowerCase(hw);
    for(unsigned int j=0;j<original.size();++j)
        if(original[j]==':')
            original[j]='#';

    vector<string> neededComs;
    vector<int> gns;
    vector<pair<bool,int> > foundComs;
    neededComs.push_back("grains");
    neededComs.push_back("loadstep");
    neededComs.push_back("variable");
    neededComs.push_back("addvalue");
    vector<pair<string,string> > comAndVal;
    retrieveComsAndVals(original,comAndVal);
    whatComsWereSpecified(comAndVal,neededComs,foundComs);

    int comLoc=-1;
    for(unsigned int j=0;j<neededComs.size();++j)
    {
        comLoc=getComLoc(neededComs[j],neededComs,foundComs);
        if(comLoc ==-1)
            throwError("Command \""+neededComs[j]+"\" was not specified in initial condition definition. ");
        else
        {
           // cout<<"command : "<<neededComs[j]<<comAndVal[comLoc].first<<" val : "<<comAndVal[comLoc].second<<endl;
            if(comAndVal[comLoc].first == "grains")
            {
                decodeGrainSpec(comAndVal[comLoc].second,gns);
            }
            else if(comAndVal[comLoc].first == "addvalue")
            {
                addval = atof(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "variable")
            {
                variable = comAndVal[comLoc].second;
            }
            else if(comAndVal[comLoc].first == "loadstep")
            {
                loadstepval = atoi(comAndVal[comLoc].second.c_str());
            }
        }//else
    } //for loop


    if(loadstepval<0)
        throwError("Invalid load step number value specified: ",loadstepval);
    else
    {  //addThisInitialVelocity,addThisInitialDisplacement;
        if(variable == "displacement")
         for(unsigned int j=0; j<gns.size();++j)
         {
             while((int)myLoadSteps[loadstepval-1].addThisInitialDisplacement.size()< gns[j])
                myLoadSteps[loadstepval-1].addThisInitialDisplacement.push_back(0.0);
            myLoadSteps[loadstepval-1].addThisInitialDisplacement[gns[j]-1]+=addval;
         }
         else if(variable == "velocity")
        for(unsigned int j=0; j<gns.size();++j)
         {
             while((int)myLoadSteps[loadstepval-1].addThisInitialVelocity.size()< gns[j])
                myLoadSteps[loadstepval-1].addThisInitialVelocity.push_back(0.0);
            myLoadSteps[loadstepval-1].addThisInitialVelocity[gns[j]-1]+=addval;
        }
         else
          throwError("Error, initial condition is either applied to velocity or displacement.");

     }


}

void Simulation::addMCVars(string hw)
{  if(myLoadSteps.size()==0)
    {
        throwError("You need to declare at least one load step!!!");return;
    }

//mcstep: load step # 3 | nmax # 1500 | shift # 0.25 | accept ratio # 0.3 | scale # 0.1 | tolerance # 0.001
    double shift = 0;
    double scale = 0;
    int maxCount = 200;
    int maxCutBacks= 10;
    double tolerance = 1e-3;
    int loadstepval=-1;
    int maxstep=1500;
    double acceptratio =0.3;
    string variable;
    string original = toLowerCase(hw);
    for(unsigned int j=0;j<original.size();++j)
        if(original[j]==':')
            original[j]='#';

    vector<string> neededComs;
    vector<int> gns;
    vector<pair<bool,int> > foundComs;
    neededComs.push_back("mcavgs");
    neededComs.push_back("loadstep");
    neededComs.push_back("shift");
    neededComs.push_back("scale");
    neededComs.push_back("acceptratio");
    neededComs.push_back("tolerance");
    neededComs.push_back("maxiterations");
    neededComs.push_back("maxcutbacks");

    vector<pair<string,string> > comAndVal;
    retrieveComsAndVals(original,comAndVal);
    whatComsWereSpecified(comAndVal,neededComs,foundComs);

    int comLoc=-1;
    for(unsigned int j=0;j<neededComs.size();++j)
    {
        comLoc=getComLoc(neededComs[j],neededComs,foundComs);
        if(comLoc ==-1)
            throwError("Command \""+neededComs[j]+"\" was not specified in initial condition definition. ");
        else
        {
           // cout<<"command : "<<neededComs[j]<<comAndVal[comLoc].first<<" val : "<<comAndVal[comLoc].second<<endl;
            if(comAndVal[comLoc].first == "mcavgs")
            {
                maxstep = atoi(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "maxiterations")
            {
                maxCount = atoi(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "maxcutbacks")
            {
                maxCutBacks = atoi(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "shift")
            {
                shift = atof(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "scale")
            {
                scale = atof(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "tolerance")
            {
                tolerance = atof(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "loadstep")
            {
                loadstepval = atoi(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "acceptratio")
            {
                acceptratio = atof(comAndVal[comLoc].second.c_str());
            }
        }//else
    } //for loop


    if(loadstepval<0)
        throwError("Invalid load step number value specified: ",loadstepval);
    else if(acceptratio > 1 || acceptratio <= 0)
        throwError("Accept ratio has to be between 0 and 1, and probably up by 0.5! ");
    else if(maxstep <=0 )
        throwError("Max step needs to be larger than 0. And probably should be large.");
    else if(shift > 1 || shift <= 0)
        throwError("Shift has to be between 0 and 1, and probably a good number is 0.25! ");
    else if(scale <= 0)
        throwError("Scale has to be positive. And probably a good number is a fraction of the smallest radius used! ");
    else if(scale <= 0)
        throwError("Tolerance has to be positive. And probably a good number is a small one like 1e-6 (problem dependent).");
   else if(maxCount <= 0)
        throwError("Max iterations has to be positive. And probably a good number is 1000 (problem dependent).");
   else if(maxCutBacks<= 0)
        throwError("Max cutbacks has to be positive (problem dependent).");

    else
    {
        myLoadSteps[loadstepval-1].myMCparam.countMax = maxCount;
        myLoadSteps[loadstepval-1].myMCparam.maxCutBacks = maxCutBacks;
        myLoadSteps[loadstepval-1].myMCparam.shift = shift;
        myLoadSteps[loadstepval-1].myMCparam.scale = scale;
        myLoadSteps[loadstepval-1].myMCparam.acceptRatio = acceptratio;
        myLoadSteps[loadstepval-1].myMCparam.Nmax = maxstep;
        myLoadSteps[loadstepval-1].myMCparam.tolerance= tolerance;
     }

}

void Simulation::addDrivingCondition(string hw)
{
//driving function: load step # 3 | add or set : add | variable : displacement | grains : 1,5,10 |  function type : cosine | params : 1,2,3|
 if(myLoadSteps.size()==0)
    {
        throwError("You need to declare at least one load step!!!");return;
    }

    string variable;
    bool addTrueSetFalse=false;
    int loadstepval=-1;
    vector<int> gns;
    string funName;
    string paramString;
    vector<double> params;

    string original = toLowerCase(hw);
    for(unsigned int j=0;j<original.size();++j)
        if(original[j]==':')
            original[j]='#';

    vector<string> neededComs;

    vector<pair<bool,int> > foundComs;
    neededComs.push_back("loadstep");
    neededComs.push_back("addorset");
    neededComs.push_back("variable");
    neededComs.push_back("grains");
    neededComs.push_back("functiontype");
    neededComs.push_back("params");
    vector<pair<string,string> > comAndVal;
    retrieveComsAndVals(original,comAndVal);
    whatComsWereSpecified(comAndVal,neededComs,foundComs);

    int comLoc=-1;
    for(unsigned int j=0;j<neededComs.size();++j)
    {
        comLoc=getComLoc(neededComs[j],neededComs,foundComs);
        if(comLoc ==-1)
            throwError("Command \""+neededComs[j]+"\" was not specified in initial condition definition. ");
        else
        {
           // cout<<"command : "<<neededComs[j]<<comAndVal[comLoc].first<<" val : "<<comAndVal[comLoc].second<<endl;
            if(comAndVal[comLoc].first == "addorset")
            {
                if(comAndVal[comLoc].second == "add")
                    addTrueSetFalse = true;
                else if(comAndVal[comLoc].second == "set")
                    addTrueSetFalse = false;
                else
                    throwError("Not sure whether you want to add or set this driving function. Check spelling on: \""+comAndVal[comLoc].second+"\".");
            }
            else if(comAndVal[comLoc].first == "variable")
            {
                variable =comAndVal[comLoc].second;
            }
            else if(comAndVal[comLoc].first == "functiontype")
            {
                funName =comAndVal[comLoc].second;
            }
            else if(comAndVal[comLoc].first == "grains")
            {
               decodeGrainSpec(comAndVal[comLoc].second,gns);
            }
            else if(comAndVal[comLoc].first == "loadstep")
            {
                 loadstepval = atoi(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "params")
            {
                paramString =comAndVal[comLoc].second;
            }
        }//else
    } //for loop

    //paramString has double values separated by commas, push them into a vector
    paramString+=',';
    string tmp,leftOver;


    while(paramString.size()>0)
    {
        splitStringBySymbol(paramString,',',tmp,leftOver);
        paramString = leftOver;
        params.push_back(atof(tmp.c_str()));
    }


    DriverFunction tmpDriver;

    if(loadstepval<0)
        throwError("Invalid load step number value specified: ",loadstepval);
    else
    {
        tmpDriver.params=params;
        tmpDriver.overRide=!addTrueSetFalse;

        if(funName == "cosine")
        tmpDriver.eval = &meth::cosine;
        else if(funName == "sine")
        tmpDriver.eval = &meth::sine;
        else if(funName == "poly")
        tmpDriver.eval = &meth::poly;
        else
            {throwError("I do not know the function you want! Either \""+funName+"\" is spelled wrong or not implemented.");
            tmpDriver.eval = &meth::zero;
            }
        if(variable == "velocity")
        tmpDriver.svptr = &(myCurrentState->velocity);
        else if(variable == "displacement")
        tmpDriver.svptr = &(myCurrentState->relativeLocation);
        else
        {
            throwError("Right now you can only apply a driving function to velocity or displacement. Check spelling on \""+funName+"\" or try something else.");
            tmpDriver.svptr = NULL;
        }


        for(unsigned int j=0; j<gns.size();++j)
        {
            tmpDriver.grainNumber = gns[j];
            if(addTrueSetFalse)
                myLoadSteps[loadstepval-1].myDriver.push_back(tmpDriver);
            else
                myLoadSteps[loadstepval-1].myDriver.push_front(tmpDriver);
        }
     }

}

void Simulation::addBodyAcceleration(string hw)
{
     if(myLoadSteps.size()==0)
    {
        throwError("You need to declare at least one load step!!!");return;
    }

//loadstep # 2 | const # 4
    string original = toLowerCase(hw);
    string com, val,leftover;
    int loadstepval = -1;
    double body = 0;

    do{
    splitStringBySymbol(original,'#',com,leftover);
    original = leftover;
    splitStringBySymbol(original,'|',val,leftover);
    original = leftover;

    if(com == "loadstep")
        loadstepval = atoi(val.c_str());
    else if(com == "const")
        body = atof(val.c_str());
    else
    throwError("Error, I don't know : \""+com+"\". Check spelling perhaps.");

    }while(leftover.size()>0);

    myLoadSteps[loadstepval-1].bodyAcceleration = body;
}
void Simulation::setRestitutionFactor(string hw)
{
     if(myLoadSteps.size()==0)
    {
        throwError("You need to declare at least one load step!!!");return;
    }

//loadstep # 2 | const # 4
    string original = toLowerCase(hw);
    string com, val,leftover;
    int loadstepval = -1;
    double rf = 0;

    do{
    splitStringBySymbol(original,'#',com,leftover);
    original = leftover;
    splitStringBySymbol(original,'|',val,leftover);
    original = leftover;

    if(com == "loadstep")
        loadstepval = atoi(val.c_str());
    else if(com == "const")
        rf = atof(val.c_str());
    else
    throwError("Error, I don't know : \""+com+"\". Check spelling perhaps.");
    }while(leftover.size()>0);



    if(rf < 0 || rf > 1 || loadstepval > (int)myLoadSteps.size() || loadstepval ==-1)
    throwError("Error, bad ratio in restitution or invalid load step number.");
    else
        myLoadSteps[loadstepval-1].restitutionFactor = rf;

}

void Simulation::addWall(string hw, bool leftwall)
{//first part is load step, then subunit name, then offset number
    if(myLoadSteps.size()==0)
    {
        throwError("You need to declare at least one load step!!!");return;
    }

    string original = hw;
    string com,val,leftover;
    int loadstepval=-1;
    double walloffset=0;
    bool gd=false;bool infradius=false;
    Grain tg;
    map<string,Grain>::iterator it;
    for(unsigned int i=0;i<original.size();++i)
    if(original[i]==':')
    original[i]='#';

    do
    {
    splitStringBySymbol(original,'#',com,leftover);
    original=leftover;
    splitStringBySymbol(original,'|',val,leftover);
    original=leftover;

    com = toLowerCase(com);
    if(com == "radius")
    {
        if(val != "inf")
        throwError("Specify wall radius through its subunit. Or If you want it to be flat, use \"radius # inf\".");
        else
        infradius=true;
    }
    else if(com=="loadstep")
        loadstepval = atoi(val.c_str());
    else if (com =="offset")
        walloffset=atof(val.c_str());
    else if (com == "subunit")
    {
        it=unitDefinitions.find(val);
        if(it==unitDefinitions.end())
        throwError("Could not find the subunit \""+val+"\" in my subunit definitions.");
        else
            {tg = unitDefinitions[val];gd=true;}
    }
    else
    throwError("Error, I don't know : \""+com+"\". Check spelling perhaps.");
    }while(leftover.size()>0);

    if(!gd || loadstepval == -1)
    throwError("You need to define a load step when this wall applies or you forgot to include what subunit the wall inherits from.");

    if(infradius)
    tg.SetbyString("radius","-1.0");//look for this when calculating prefactors

    tg.SetbyString("volume","-1.0");//look for this when calculating prefactors


    if(leftwall)
    {
        myLoadSteps[loadstepval-1].leftWallGrain = tg;
        myLoadSteps[loadstepval-1].leftWallOffset = walloffset;
    }
    else
    {
        myLoadSteps[loadstepval-1].rightWallGrain = tg;
        myLoadSteps[loadstepval-1].rightWallOffset = walloffset;
    }


}
void Simulation::addLeftWall(string hw)
{
    addWall(hw,true);
}
void Simulation::addRightWall(string hw)
{
    addWall(hw,false);
}



std::pair<std::string, double> Simulation::getNameAndNumber(std::string namenum,char symbol)
{
	bool after=false;
	std::string nameval="",numval="";
	for (unsigned int lcv=0; lcv<namenum.length();++lcv)
	{

		if(namenum[lcv]==symbol)
			after=true;
		else
		{
			if(after)
			{
			numval = numval+namenum[lcv];
			}
			else
			{
			nameval = nameval+namenum[lcv];
			}

		}

	}

	if(numval.empty())
		numval="1";

	return pair<string,double>(nameval,atof(numval.c_str()));
}
void Simulation::decodeChainPattern(string hw)
{
	string temp="";
	string newhw="";
///First runs through and expands the 'mega' repeats []
	bool foundmega=false;
    for(unsigned int lcv =0; lcv<hw.length();++lcv)
	{
		if (hw[lcv]=='[')
		{		temp = "";
			foundmega=true;
		}
		else if(hw[lcv]==']')
		{
			pair<string,int> megarepeat = getNameAndNumber(temp,'^');
			for(int megalcv=0; megalcv<megarepeat.second;++megalcv)
			{
				newhw=newhw + megarepeat.first + '|';
			}
			foundmega=false;
		}
		else
		{
			if(foundmega)
			temp = temp+hw[lcv];
			else
			newhw=newhw+hw[lcv];
		}
	}
//now the chain code is expanded at the highest level
hw = newhw;
hw+='|';

temp="";
int rep;
string callsign;
//run through the new chain code and expand
	for(unsigned int lcv =0; lcv<hw.length();++lcv)
	{
		//std::cout<<hw[lcv]<<std::endl;
		if (hw[lcv]=='|')
		{
		//	std::cout<<"breaker"<<std::endl;
			//std::cout<<getNameAndNumber(temp,'#').first<<" "<<getNameAndNumber(temp,'#').second<<std::endl;
			if(!temp.empty())
			{
            callsign=getNameAndNumber(temp,'#').first;
            rep = getNameAndNumber(temp,'#').second;
			for(int j=0;j<rep;++j)
                chainPattern.push_back(callsign);
			}
			temp = "";
		}
		else
		temp = temp+hw[lcv];
	}

/*cout<<"decoded: ";
for(unsigned int j=0;j<chainPattern.size();++j)
cout<<chainPattern[j]<<"|";
cout<<endl;*/
}
void Simulation::addDefaultSubunit(string hw)
{
    addSubunit("DEFAULTGRAIN|"+ hw);
}
void Simulation::addSubunit(string hw)
{
//white space was removed, parse the unit command into a shape
//add this shape to the unitdefinitions map
    string name;
    Grain tempGrain;

    string leftover;
    string param,val;
    if(hw[hw.length()-1]!='|')
    hw+='|';
//a subunit is a name | params # val | repeat

    //get name first
    splitStringBySymbol(hw,'|',name,leftover);
    //this leaves us with | params # val | repeat
 hw = leftover;
 if(name != "DEFAULTGRAIN" && defaultDefined)
 tempGrain = unitDefinitions["DEFAULTGRAIN"];
else if(name == "DEFAULTGRAIN" && !defaultDefined)
defaultDefined = true;
else
throwError("You need to define the DEFAULTGRAIN first! \n add a default subunit before any others in text file. You should never read this message because the priority queue solves it.");

tempGrain.Setvolume(0);
tempGrain.Setmass(0);

    while(hw.size()>0)
    {
        //attempt to split off a parameter and value:
        splitStringBySymbol(hw,'#',param,leftover);
        hw=leftover;
        splitStringBySymbol(hw,'|',val,leftover);
        hw=leftover;

        if(!tempGrain.SetbyString(toLowerCase(param),val))
            haultFlag=true;

    }

        ///now need to calculate geometric parameters for this grain,
        ///namely volume and then mass
        ///**** IF the volume or mass are assigned in the subunit, they trump the volume calculation
        ///Is this a good idea?

    if(taperingRatio==1)
        tempGrain.computeTheRest();

    unitDefinitions.insert(pair<string,Grain>(name,tempGrain));
}
void Simulation::setTapering(string hw)
{
    taperingRatio=atof(hw.c_str());
}
void Simulation::addLoadStep(string hw)
{
    string original = toLowerCase(hw);
    for(unsigned int j=0;j<original.size();++j)
        if(original[j]==':')
            original[j]='#';
    int loadstepval = -1;
    vector<string> neededComs;
    vector<pair<bool,int> > foundComs;
    LoadStep temp(0);///NEED TO THINK THIS THROUGH MOVABLE???!!?!?!
    neededComs.push_back("number");
    neededComs.push_back("solver");
    neededComs.push_back("start");neededComs.push_back("step");neededComs.push_back("stop");

    vector<pair<string,string> > comAndVal;

    retrieveComsAndVals(original,comAndVal);
    whatComsWereSpecified(comAndVal,neededComs,foundComs);

    int comLoc=-1;
    for(unsigned int j=0;j<neededComs.size();++j)
    {
        comLoc=getComLoc(neededComs[j],neededComs,foundComs);
        if(comLoc ==-1)
            throwError("Command \""+neededComs[j]+"\" was not specified in load step definition. ");
        else
        {
           // cout<<"command : "<<neededComs[j]<<comAndVal[comLoc].first<<" val : "<<comAndVal[comLoc].second<<endl;
            if(comAndVal[comLoc].first == "solver")
            {
                 if( toLowerCase(comAndVal[comLoc].second) == "velocityverlet" )
                  {
                     temp.mySolver = &Simulation::velocityVerletStep;
                     temp.myUpdater = &Simulation::incrementTimeExplicit;
                  }
                else if(comAndVal[comLoc].second == "montecarloequilibrium")
                {
                    temp.mySolver = &Simulation::equilibriumSolverStep;
                    temp.myUpdater = &Simulation::incrementTimeImplicit;
                }
                else if(comAndVal[comLoc].second == "rungekuttafour")
                {
                    temp.mySolver = &Simulation::rungeKuttaOrderFour;
                    temp.myUpdater = &Simulation::incrementTimeExplicit;
                }
                 else if(comAndVal[comLoc].second == "equilibriumbodyonly")
                {
                    temp.mySolver = &Simulation::equilibriumSolverBodyOnly;
                    temp.myUpdater = &Simulation::incrementTimeImplicit;
                }
                else
                throwError("Check you spelling on your solver, I dont know what solver \""+comAndVal[comLoc].second+"\" is, check manual for spelling?");
            }
            else if(comAndVal[comLoc].first == "number")
            {
                loadstepval = atoi(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "start")
            {
                temp.startTime = atof(comAndVal[comLoc].second.c_str());
            }
            else if(comAndVal[comLoc].first == "step")
            {
                temp.dt = atof(comAndVal[comLoc].second.c_str());
            }else if(comAndVal[comLoc].first == "stop")
            {
                temp.stopTime = atof(comAndVal[comLoc].second.c_str());
            }

        }//else
    } //for loop

if(loadstepval<0)
    throwError("Invalid load step number value specified: ",loadstepval);
else
{
    if((int)myLoadSteps.size()<loadstepval)
        myLoadSteps.resize(loadstepval);

    myLoadSteps[loadstepval-1]=temp;
    cout<<"Assigning loadstep number : "<<loadstepval<<endl;
}


}

///*********************************************///
///More initialization of the Chain
void Simulation::fillChainVector(string hw) //after decode, fill vector
{
    for(unsigned int j=0;j<chainPattern.size();++j)
         myChain.push_back(unitDefinitions[chainPattern[j]]);
}

void Simulation::applyTapering(string hw)
{ ///assumes that myChain has stuff in it!!!!!
    Shape ts;

    if(taperingRatio!=1 && taperingRatio !=-1)
    {
            if(taperingRatio>0)
            for(unsigned int j=1;j<myChain.size();++j)
            {
                ts = myChain[j].GetmyShape();//current

                ts.r_minus=ts.r_minus*pow(taperingRatio,(double)j);
                ts.r_plus=ts.r_plus*pow(taperingRatio,(double)j);

                myChain[j].SetmyShape(ts);
            }
            else
            for(int j=myChain.size()-2;j>=0;--j)
            {

                ts = myChain[j].GetmyShape();//current

                ts.r_minus=ts.r_minus*pow(-taperingRatio,(double)(myChain.size()-1-j));
                ts.r_plus=ts.r_plus*pow(-taperingRatio,(double)(myChain.size()-1-j));

                myChain[j].SetmyShape(ts);

            }
    }

}
void Simulation::fillInMissingGrainInfo(string hw)//calculate the volume / mass for each store this in massvector!!!
{
    for(unsigned int j=0;j<myChain.size();++j)
    {
        //cout<<"I found a tapering ratio :"<<taperingRatio<<" it is equal to 1? "<<(taperingRatio==1?"yes":"no")<<endl;
        if(taperingRatio!=1)
            myChain[j].computeTheRest();
        massVector.push_back(myChain[j].Getmass());
    }


 //resize the current state
    myCurrentState->resize(myChain.size());

//X^{abs_und}_i = Sum_j=0^{i-1}[ R^+_j + R^-_j + w_j ]  +  R^-_i + w_i/2
//X^{abs}_i(t) = X^{rel}_i(t) + X^{abs_und}_i

    //calculate the undeformed absolute positions:
    Shape ts;
    ts = myChain[0].GetmyShape();
    myCurrentState->initialAbsolutePosition[0]=ts.r_minus+ts.width/2.0;
    for(unsigned int j=1; j<myChain.size();++j)
     {
         (myCurrentState->initialAbsolutePosition)[j]=(myCurrentState->initialAbsolutePosition)[j-1]+ts.r_plus+ts.width/2.0;
         ts = myChain[j].GetmyShape();
         (myCurrentState->initialAbsolutePosition)[j]+=ts.r_minus+ts.width/2.0;
     }
    //now when you need the absolute positions, it is the xrel + the initial abs pos.
        //assumes that there is something in the ChainPattern and subunits

    if(outputGrainNumbers.size()==0)
    for(unsigned int j=0;j<chainPattern.size();++j)
        outputGrainNumbers.push_back(j+1);



    while((unsigned int)outputGrainNumbers[outputGrainNumbers.size()-1]>chainPattern.size())//user doesn't know about 0->N-1
        outputGrainNumbers.pop_back();

    //now determine what interfaces these correspond to
    //the values of these run from 0->N (0 is leftwall/grain1 interface, 1 is grain1/grain2, ..., N is grain N/RightWal etc)
    vector<int> temp;
    for(unsigned int j=0;j<outputGrainNumbers.size();++j)
    {
          //first add in all,
        temp.push_back(outputGrainNumbers[j]-1);//outputGN has values of 1->N
        temp.push_back(outputGrainNumbers[j]);
    }
    //then remove duplicates
    for(unsigned int j=0;j<temp.size();++j)
    if(!alreadyThere(outputInterfaceNumbers,temp[j]))
        outputInterfaceNumbers.push_back(temp[j]);



}
///NOW WE are getting into some of the meat, need to figure out
///what the power law is between two grains, this depends on the
/// material properties, and the grain shape.
void Simulation::calculatePowerLawAndPrefactors(string hw)
{
//run through each grain pairing
    for(unsigned int j=0;j<myChain.size()-1;++j) //inbetweens
    {
        forcePrefactors[pair<int,int>(j,j+1)]= myChain[j].calculateForcePrefactor(myChain[j+1]);//calculate force prefactor between grain j and j+1
        forcePowerLaw[pair<int,int>(j,j+1)]= myChain[j].calculateForcePower(myChain[j+1]);
    }

    for(unsigned int j=0; j<myLoadSteps.size();++j)
    {   //run through each load step and calculate the wall prefactors if there is a wall
        if(myLoadSteps[j].leftWallGrain.Getvolume()==-1)
          {
                myLoadSteps[j].leftWallForcePrefactor = myLoadSteps[j].leftWallGrain.calculateForcePrefactor(myChain[0]);
                myLoadSteps[j].leftWallForcePower = myLoadSteps[j].leftWallGrain.calculateForcePower(myChain[0]);
          }
        if(myLoadSteps[j].rightWallGrain.Getvolume()==-1)
          {
                myLoadSteps[j].rightWallForcePrefactor = myChain[myChain.size()-1].calculateForcePrefactor(myLoadSteps[j].rightWallGrain);
                myLoadSteps[j].rightWallForcePower = myChain[myChain.size()-1].calculateForcePower(myLoadSteps[j].rightWallGrain);
          }
    }

}
void Simulation::fillInMissingLoadStepInfo(string hw)
{

    //Namely fill up the movable grains in each load step if they were
    //not specified earlier in some command (default is they are movable)
  for(unsigned int j=0; j<myLoadSteps.size();++j)
        while(myLoadSteps[j].movable.size() < myChain.size())
            myLoadSteps[j].movable.push_back(true);



    for(unsigned int j=0; j<myLoadSteps.size();++j)
        while(myLoadSteps[j].addThisInitialDisplacement.size() < myChain.size())
            myLoadSteps[j].addThisInitialDisplacement.push_back(0.0);

    for(unsigned int j=0; j<myLoadSteps.size();++j)
        while(myLoadSteps[j].addThisInitialVelocity.size() < myChain.size())
            myLoadSteps[j].addThisInitialVelocity.push_back(0.0);

   if(myLoadSteps.size()==0)
        throwError("You need to at least specify one load step!");
}
void Simulation::caclulateRemainingStateVars(LoadStep &details)
{
    //cout<<"Calculating absolute position, kinetic energy, etc"<<endl;

    ///assumes that the myCurrentState variable has valid info in its xrel and vel
    myCurrentState->absolutePosition = myCurrentState->initialAbsolutePosition+myCurrentState->relativeLocation;

    for(unsigned int j=0;j<myChain.size();++j)
      {
         myCurrentState->kineticEnergy[j]=pow(myCurrentState->velocity[j],2)*0.5*massVector[j];
         ///myCurrentState->potentialEnergy[j]+=massVector[j]*meth::absolute(details.bodyAcceleration)*(details.bodyAcceleration<0? myCurrentState->absolutePosition[j] :   myCurrentState->absolutePosition[myCurrentState->absolutePosition.size()-1]-myCurrentState->absolutePosition[j]);
      }


    /// / ///cout<<"Consider calculating lost energy here?"<<endl;
}
///RUN IT
void Simulation::runSimulation()
{///work through each load step and iterate them until done
    double time=myLoadSteps[0].startTime;

    double desiredDt;
    double timeLeftToOutput;

    //next state is a scratch variable that gets the updated info
    //and then copied back over to myCurrentState
    State* nextState =new State(*myCurrentState);
    myPrevState = new State(*myCurrentState);//create a copy of the current.
    myCurrentState->curtime=time;

    caclulateRemainingStateVars(myLoadSteps[0]);

    runThroughOutputs(true);//first time add titles

    for(unsigned int j=0; j<myLoadSteps.size();++j)
    {
        if(haultFlag)
            break;
        if(!isSilent())
        cout<<"At load step : "<<j+1<<endl;

        time = myLoadSteps[j].startTime;

        desiredDt = myLoadSteps[j].dt;
        timeLeftToOutput = 0;//outputInterval;
        myCurrentState->curtime=time;
        myCurrentState->relativeLocation+=myLoadSteps[j].addThisInitialDisplacement;
        myCurrentState->velocity+=myLoadSteps[j].addThisInitialVelocity;

      //  cout<<(*myCurrentState)<<endl;
      while(time<myLoadSteps[j].stopTime && !haultFlag)
        {
            //determine the next time step
            (this->*myLoadSteps[j].myUpdater)(time,timeLeftToOutput,myLoadSteps[j].stopTime,desiredDt, myLoadSteps[j].dt );
            nextState->curtime=time;
            myCurrentState->curtime=time;
            //update the driving conditions
            applyDrivingConditions(*myCurrentState,myLoadSteps[j]);

            //Calll the solver for this time step
             (this->*myLoadSteps[j].mySolver)(*myCurrentState,*nextState,myLoadSteps[j]);

         //   calculate variables
        caclulateRemainingStateVars(myLoadSteps[j]);

           // output variables
           if(timeLeftToOutput==0)
           {
               if(!logFileOnly)
                cout<<"Outputting time step "<<time<<" during load step : "<<j+1<<endl;

                runThroughOutputs(false);//not the first time
                timeLeftToOutput=outputInterval;
           }//time left to output if statement

            *myPrevState = *myCurrentState; //copy the contents, do not change pointers
            *myCurrentState = *nextState; //copy the contents, do not change pointers

        } //while loop for current state


    }
}
///Calculate forces
void Simulation::calculateForces(State &curState,LoadStep &details)
{  //given the current state (velocities, xrelative locations, and any extras in details) compute forces on grains.

    //Hertz Force between grains (less restitution)
 /*delta_{i,i+1}(t) = D^{und}_{i,i+1} - D_{i,i+1}(t)
                = D^{und}_{i,i+1} - (X^{rel}_{i+1}(t) -  X^{rel}_{i}(t)  +  D^{und}_{i,i+1})
                = D^{und}_{i,i+1} - X^{rel}_{i+1}(t) +  X^{rel}_{i}(t)  -  D^{und}_{i,i+1}
                =  - X^{rel}_{i+1}(t) +  X^{rel}_{i}(t)*/
    //clear any forces
    double delta=0;
    double prevdelta = 0;
    double tforce=0;
    double reductionFactor=1;
    pair<int,int> tpair;

    for(unsigned int i=0;i<curState.force.size();++i)
        {
            curState.force[i] = 0;
            curState.overlap[i]=0;
        }
        curState.overlap[curState.overlap.size()-1]=0;
    for(unsigned int i=0;i<curState.force.size()-1;++i)
    {
        //for interface params, PE,intForce,overlap: index 0 is LeftWall/1stGrain(0th in code), index 1 is 1stGrain(0th in code)/2ndGrain(1st in code), last index N is NthGrain(N-1th in code)/RightWall
        tpair.first=i;tpair.second=i+1;

        delta = -curState.relativeLocation[i+1] + curState.relativeLocation[i];
        prevdelta = -myPrevState->relativeLocation[i+1] + myPrevState->relativeLocation[i];
        prevdelta = (prevdelta<=0?0:prevdelta);

        curState.overlap[i+1]=(delta<=0 ? 0 : delta);

        if( curState.overlap[i+1] >= prevdelta ) //current overlap is larger -> compression
            reductionFactor=1.0;
        else
            reductionFactor=details.restitutionFactor;

        tforce = reductionFactor*forcePrefactors[tpair]*pow(curState.overlap[i+1],forcePowerLaw[tpair]);
        curState.interfaceForce[i+1]=tforce;
        curState.potentialEnergy[i+1]=1.0/(forcePowerLaw[tpair]+1)*forcePrefactors[tpair]*pow(curState.overlap[i+1],forcePowerLaw[tpair]+1);
        //subtract this force from i add it to i+1
        curState.force[i]-=tforce;
        curState.force[i+1]+=tforce;
    }

    //Wall Potentials
//    A wall overlap will happen when  X^{rel}_{N-1} > 0 or X^{rel}_{0} < 0
//delta_{left wall} =  -X^{rel}_{0}+WALLOFFSETLEFT  // (-) makes the delta positive in compression
//delta_{right wall} = X^{rel}_{N-1}-WALLOFFSETRIGHT // (+) here keeps delta positive in compression

    //left wall calculation
    delta = -curState.relativeLocation[0] + details.leftWallOffset;
    curState.overlap[0]=(delta<=0 ? 0 : delta);
    prevdelta = -myPrevState->relativeLocation[0] - details.leftWallOffset;
    prevdelta = (prevdelta<=0?0:prevdelta);
    if( curState.overlap[0] >= prevdelta ) //current overlap is larger -> compression
            reductionFactor=1.0;
        else
            reductionFactor=details.restitutionFactor;


    tforce = details.leftWallForcePrefactor*pow(curState.overlap[0],details.leftWallForcePower)*reductionFactor;
    curState.interfaceForce[0]=tforce;
    curState.potentialEnergy[0]=1.0/(details.leftWallForcePower+1.0)*details.leftWallForcePrefactor*pow(curState.overlap[0],details.leftWallForcePower+1.0);
    curState.force[0]+=tforce; //add this force to grain 0 => grain 1 by user notation

    //right wall calculation
    delta = curState.relativeLocation[curState.relativeLocation.size()-1] - details.rightWallOffset;
    curState.overlap[curState.overlap.size()-1]=(delta<=0 ? 0 : delta);


    prevdelta = myPrevState->relativeLocation[myPrevState->relativeLocation.size()-1] - details.rightWallOffset;
    prevdelta = (prevdelta<=0?0:prevdelta);
    if( curState.overlap[curState.overlap.size()-1] >= prevdelta ) //current overlap is larger -> compression
            reductionFactor=1.0;
    else
            reductionFactor=details.restitutionFactor;


    tforce = details.rightWallForcePrefactor*pow(curState.overlap[curState.overlap.size()-1],details.rightWallForcePower)*reductionFactor;
    curState.interfaceForce[curState.interfaceForce.size()-1]=tforce;
    curState.potentialEnergy[curState.potentialEnergy.size()-1]=1.0/(details.rightWallForcePower+1.0)*details.rightWallForcePrefactor*pow(curState.overlap[curState.overlap.size()-1],details.rightWallForcePower+1.0);
    curState.force[curState.force.size()-1]-=tforce; //subtract this force to grain N-1 => grain N by user notation



    //Body forces
    if(details.bodyAcceleration!=0)
        curState.force+=details.bodyAcceleration*massVector;

    //Viscous forces/Velocity Based Damping
    curState.force+=-1.0*details.velocityBasedDamping*curState.velocity;



    //If unmovable, set force = 0, and velocity =0:
    curState.force*=details.movable;
    curState.velocity*=details.movable;

}
void Simulation::applyDrivingConditions(State &curState,LoadStep &details)
{ /// run through the drivers in the load step details, add or set the state variable where appropriate

////MP cout<<"Need to update this to be compatible with freeze/unfreeze. Maybe take actual derivative (from the driving functions you already know) of displacement to get velocity?"<<endl;
   //this really only acts on myCurrentState and curSTate is used for the time. kinda bad
    list<DriverFunction>::iterator it;

    for(it = (details.myDriver).begin();it!=details.myDriver.end();++it)
    {
        if(it->overRide)
        it->svptr->setAt((it->grainNumber)-1,(it->eval)(curState.curtime,it->params));
        else
        it->svptr->addAt((it->grainNumber)-1,(it->eval)(curState.curtime,it->params));
    }

}

///***************************************************************************************
///SOLVERS IMPLEMENTED BELOW -----
void Simulation::velocityVerletStep(State &curState,State &nextState,LoadStep &details)
{
    //curstate is current, nextState is what we will work into

    ///Velocity Verlet Algorithm
/*
    v(n+.5) = v(n) + a(n)*dt/2
    x(n+1) = x(n) + v(n+.5)*dt
    update: a(n+1) using x(n+1)
    v(n+1) = v(n+.5) + a(n+1)*dt/2
*/
    nextState = curState;
    nextState.velocity+=details.dt/2.0*curState.force/massVector;
    nextState.relativeLocation+=nextState.velocity*details.dt;
    calculateForces(nextState,details);
    nextState.velocity+=details.dt/2.0*nextState.force/massVector;

}
void Simulation::rungeKuttaOrderFour(State &curState,State &nextState,LoadStep &details)
{///Apply BC conditions is called before this solver and affects the current state.
    ///the current choice is to let the BCs just be fixed during these solver steps
    ///might rethink that in the future.

        ///Runge Kutta Order 4 :
    /*
            if y' = f(t,y)
            and you want soln at t+h
            w = y(t) (at current step)

            k1 = h*f(t,w)
            k2 = h*f(t + h/2, w + k1/2)
            k3 = h*f(t + h/2, w + k2/2)
            k4 = h*f(t + h  , w + k3 )

            y(t+h) = w + (k1 + 2k2 + 2k3 + k4)/6

            our y is two first order systems

            ma = calcForces
            veldot = calcForces/m
            x-dot = vel

            so y = [vel x]
            and y' = [calcForces/m vel]
            [] = k1 = h*f(t,w) = h*[calcForces/m vel]

            need to settle on vel being evaluated at its most recent (for damping in calcForces)
            We need a State for each k1->k4 because thats how calcForces is designed.
    */
    ///details has the newer time step already but curState has the 'old' info. All time
    ///dependence of our derivatives was moved to Apply BCs so really t has no role, but the
    ///time step does (which we have in details).
            double h = details.dt; ///newest time step
            State w(curState.relativeLocation.size()),
            K1(curState.relativeLocation.size()),
            K2(curState.relativeLocation.size()),
            K3(curState.relativeLocation.size()),
            K4(curState.relativeLocation.size());

            w = curState;
            calculateForces(w,details); ///w has y'
            K1.velocity=w.force/massVector*h;
            K1.relativeLocation=(curState.velocity+K1.velocity)*h;

            w.velocity=curState.velocity+K1.velocity/2.0;
            w.relativeLocation=curState.relativeLocation+K1.relativeLocation/2.0;
            calculateForces(w,details);
            K2.velocity=w.force/massVector*h;
            K2.relativeLocation=(curState.velocity+K2.velocity)*h;

            w.velocity=curState.velocity+K2.velocity/2.0;
            w.relativeLocation=curState.relativeLocation+K2.relativeLocation/2.0;
            calculateForces(w,details);
            K3.velocity=w.force/massVector*h;
            K3.relativeLocation=(curState.velocity+K3.velocity)*h;

            w.velocity=curState.velocity+K3.velocity;
            w.relativeLocation=curState.relativeLocation+K3.relativeLocation;
            calculateForces(w,details);
            K4.velocity=w.force/massVector*h;
            K4.relativeLocation=(curState.velocity+K4.velocity)*h;

            nextState.relativeLocation= curState.relativeLocation+(K1.relativeLocation+2.0*K2.relativeLocation+2.0*K3.relativeLocation+K4.relativeLocation)/6.0;
            nextState.velocity = curState.velocity+(K1.velocity+2.0*K2.velocity+2.0*K3.velocity+K4.velocity)/6.0;
            calculateForces(nextState,details);//this function does not affect velocities


}
void Simulation::equilibriumSolverBodyOnly(State &curState,State &nextState,LoadStep &details)
{ ///if there is only a fody force, balance it with the
        if(details.bodyAcceleration < 0 && details.leftWallForcePrefactor ==0 )
            {throwError("You need to have a left wall if you are using a negative acceleration and this solver.");return;}
        if(details.bodyAcceleration > 0 && details.rightWallForcePrefactor ==0 )
            {throwError("You need to have a right wall if you are using a positive acceleration and this solver.");return;}

        ///if accel is less than 0, need to apply mass of al grains to the right and repeat
        ///if accel is greater than 0 need to apply mass of all grains to the left and repeat
        ///code one case first and then flip as necessary?
        double tforce,overlap;
       OneDimDoubleVector cumulativeSumMass = massVector;
        pair<int,int> tpair;
        nextState=curState;
        if(details.bodyAcceleration<0) ///handle left wall first then move rightward.
        {
            for( int j = cumulativeSumMass.size()-1;j>=0;--j)
                cumulativeSumMass[j]=cumulativeSumMass[j]+cumulativeSumMass[j+1];


            //now solve for overlap at left wall and 0th grain, has to supported cumulativeSUmMass[0]
            tforce = cumulativeSumMass[0]*meth::absolute(details.bodyAcceleration);
            tforce /= details.leftWallForcePrefactor;
            overlap = pow(tforce,1.0/details.leftWallForcePower);
            nextState.relativeLocation[0]=-overlap + details.leftWallOffset;

            //now repeat this procedure for teh remaining grains
            for(unsigned int j=1; j<massVector.size();++j)
            {
                tforce = cumulativeSumMass[j]*meth::absolute(details.bodyAcceleration);
                tpair.first = j-1;
                tpair.second = j;
                tforce /= forcePrefactors[tpair];
                overlap = pow(tforce,1.0/forcePowerLaw[tpair]);
                nextState.relativeLocation[j]=-overlap + nextState.relativeLocation[j-1];
            }

        }
        else  ///go the other way
        {
            for(unsigned int j = 1 ; j < cumulativeSumMass.size();++j)
                cumulativeSumMass[j]=cumulativeSumMass[j]+cumulativeSumMass[j-1];



            //now solve for overlap at left wall and 0th grain, has to supported cumulativeSUmMass[0]
            tforce = cumulativeSumMass[cumulativeSumMass.size()-1]*meth::absolute(details.bodyAcceleration);
            tforce /= details.rightWallForcePrefactor;
            overlap = pow(tforce,1.0/details.rightWallForcePower);
            nextState.relativeLocation[curState.relativeLocation.size()-1]=overlap + details.rightWallOffset;

            //now repeat this procedure for the remaining grains
            for( int j=nextState.relativeLocation.size()-2; j>=0;--j)
            {
                tforce = cumulativeSumMass[j]*meth::absolute(details.bodyAcceleration);
                tpair.first = j;
                tpair.second = j+1;
                tforce /= forcePrefactors[tpair];
                overlap = pow(tforce,1.0/forcePowerLaw[tpair]);
                nextState.relativeLocation[j]= nextState.relativeLocation[j+1]+overlap;
            }

        }
    calculateForces(nextState,details);

}

void Simulation::equilibriumSolverStep(State &curState,State &nextState,LoadStep &details)
{
///sum (abs ( forces*movable )) < thresh ?
///applyDrivingConditions(*myCurrentState,myLoadSteps[j]);
    double diff = details.myMCparam.tolerance*10;//want to get into loop!
    double scale = details.myMCparam.scale;
  //  for(unsigned int j=0; j<curState.relativeLocation.size();++j)
    //if(meth::absolute(details.addThisInitialDisplacement[j])*0.5>scale)
      //  scale = meth::absolute(details.addThisInitialDisplacement[j])*0.5;

    //cout<<"here:"<<details.myMCparam;
  //  cout<<"scale I chose: "<<scale<<endl;


    int count = 0;
    int cutback=0;
    double lastdiff=100;
  //reapply boundary conditions
do{
     applyDrivingConditions(curState,details);
    count = 0;
    while( count < details.myMCparam.countMax && diff > details.myMCparam.tolerance )
    {
    equilibriumHelperIterateStep(curState,nextState,details,scale);
    ++count;

    applyDrivingConditions(curState,details);
 //   cout<<" nextstate : "<<nextState.relativeLocation<<endl;
    //multiply by movable so we ignore them
    *myPrevState=nextState;
    calculateForces(nextState,details);

   // cout<<"guess acceel: "<<temp<<endl;
    diff=meth::meanAbsVec(nextState.force);

    //runThroughOutputs(false);
    if(diff<lastdiff)
        curState = nextState;

    lastdiff=diff;

    if(!logFileOnly)
    cout<<"count : "<<count<<" error  : "<<diff<<endl;

    }
    if(count == details.myMCparam.countMax)
    {   ++cutback;
        scale*=0.5;
        details.myMCparam.scale*=0.5;
    }
}while(count == details.myMCparam.countMax && cutback < details.myMCparam.maxCutBacks);
//after loop, fill in overlaps and compute remaining ineteresting values for output

    *myPrevState=curState;
    calculateForces(curState,details);
    if(cutback == details.myMCparam.maxCutBacks)
        throwError("Equilibrium solver failed! I reached the maximum number of cutbacks. \nYou can try adding more, but more likely you need to apply the boundary condition more slowly.");
}
void Simulation::equilibriumHelperIterateStep(State &curState,State &nextState,LoadStep &details,double &scale)
{ ///basically use Monte Carlo approach since there are N dimensions to minimize accelerations across
    int maxCases = details.myMCparam.Nmax;
    int seed = meth::timeSeed();
    nextState.relativeLocation*=0;


    double diff=0;
    double acceptratio = details.myMCparam.acceptRatio;
    double curdiff;
    double maxforce;
    curdiff=meth::meanAbsVec(curState.force);
       //test against this value
   // cout<<"curdiff: "<<curdiff<<endl;
    vector<State> testCases;
    //generate a number of test cases
      *myPrevState=curState;

      calculateForces(curState,details);
        maxforce = meth::maxval(1.0,meth::maxVec(curState.force));
     //   cout<<"found max force: "<<maxforce<<endl;

    for( int j=0; j<maxCases;++j)
        {


            nextState=curState;
            for(unsigned int i=0; i<nextState.relativeLocation.size()-1;++i)
            {
                //cout<<i<<endl;
                //cout<<"this guys mass: "<<i<<"  m: "<<massVector[i]<<endl;
                //cout<<"addme : "<<i<<" ? "<<(meth::ranf(seed,false)-details.myMCparam.shift)*2.0*scale/maxforce<<endl;
                     nextState.relativeLocation[i]+=(meth::ranf(seed,false)-details.myMCparam.shift)*2.0*scale*(curState.force[i]*(details.movable[i]?1.0:0.0))/massVector[i];
            }

            *myPrevState=nextState;
            calculateForces(nextState,details);


            diff=meth::meanAbsVec(nextState.force);

            // cout<<"test : "<<j<<" has diff val : "<<diff<<endl;
            if( curdiff > diff || meth::ranf(seed,false) < acceptratio)
                testCases.push_back(nextState);
        }

    nextState.relativeLocation*=0;
    for(unsigned int j=0; j<testCases.size();++j)
        nextState.relativeLocation+=testCases[j].relativeLocation;
    nextState.relativeLocation/=((double)testCases.size());

        diff = meth::meanAbsVec(nextState.force);

}
///DONE IMPLEMENTING SOLVERS --- IN THE FUTURE MAYBE NICE TO BREAK THIS SECTION OUT TO A METHODS CLASS
///***************************************************************************************
void Simulation::incrementTimeExplicit(double &curTime, double &intervalRemaining, const double &stopTime, const double &trialDt, double &newDt)
{
    if(curTime+trialDt > stopTime) //time step overshoots the duration we are interested in, so back it up
    {
        intervalRemaining=0;
        newDt=stopTime-curTime;
        curTime = stopTime;
    }
    else   //time step is within the stopping time
    {
        if(intervalRemaining-trialDt < 0) //but time step overshoots a time spit
        {
            curTime+=intervalRemaining;
            newDt=intervalRemaining;
            intervalRemaining=0;
        }
        else
        {
            curTime+=trialDt;
            intervalRemaining-=trialDt;
            newDt=trialDt;
        }
    }


}
void Simulation::incrementTimeAdaptive(double &curTime, double &intervalRemaining, const double &stopTime,const double &trialDt, double &newDt)
{
 cout<<"The  option of adaptive time stepping is not implemented yet!!!!"<<endl;
}
void Simulation::incrementTimeImplicit(double &curTime, double &intervalRemaining, const double &stopTime,const double &trialDt, double &newDt)
{
        if(trialDt + curTime < stopTime)
        {
            newDt = trialDt;
            curTime +=newDt;
        }
        else
        {
            newDt = stopTime - curTime;
            curTime+=newDt;
        }
    intervalRemaining=0;
}

