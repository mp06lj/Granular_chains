//TaperChain102704v3.cpp
//Rewriting of the Taperchain programs by Adam Sokolow 
//(velocity verlet algorithm taken from original code written by Jan Pfannes)
//In an attempt to upgrade the code...
//7/8/03
//Now actual date: 11/20/04
//Updated again 2.22.05
//5/19/07 Another revision for the breathing study, some preprogrammed driving functions added.
//10/13/08 adding precompression functionality for user and specification of initial compressions to individual grains

//Reverting back to using the large particle as particle N-1, and small as 0	  7/9/03
 
/***************************************************
READ THIS PLEASE:


This code has been altered so that it will run from a file...
any parameter that you would like to change should be done so
through this built in file interface... otherwise this program
may produce undesired results.. or no results at all.

Also if you have little to no programming experience, keep in mind
that every little thing does matter, and if something gets accidentally
removed, the code may still compile and run, and the results you could get
could be completely believable, showing no obvious sign of error... 

Look... but be very careful. (and save a backup)

******************************************************/

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
//#include <malloc.h>
#include <stdio.h>
#include <time.h>




//***********************CONSTANTS/VARIABLES BETWEEN RUNS

 const double PI = 4 * atan(1.0);
 int nptles=20;
                   // total number of particles


 //const bool Wall = true;
int DEFAULTPRECISION = 10;
 bool leftWall=true;
 bool rightWall=true;
//const double rho = 3.2;				//SiC
//const double D = 0.00326603139013;

//TiAlV material constants are below...(just swap which is commented out)

//const double rho=4.42,//TiAlV            3.2 /* SiC (mg/mm^3) */, 
 // D = 0.01206;//D=0.00326603139013 /* (mm^2/kN) */;

  double rho = 1.00, D = 0.1529; //rho(mg/mm^3), D(mm^/kN)
 
//const double rho = 7.780,//Stainless Steel
//D = 0.00825758;//D = 0.008327;//D = 0.0066619;
//D = 0.00196;


 double rlarge = 0.5;           // (radius of large ptle (mm))
double q = 0;      // (tapering factor (%), q=0: monodisp.)	 //changed below in code...

double xn = 2.5;                   // (exponent in potential)
 double dt = 0.00001;               // (timestepwidth (musec))
 unsigned long long int nsteps = 500000;//12500000;// (# steps integration loop)

 double timeSpits = .5;  //in the output files, how frequent do we want output... 
 double ForceSpits = timeSpits;
 int realspit = int(timeSpits/dt);
 int realFspit = int(ForceSpits/dt);


//const double noiseDump = 1.0e-30;  // (max. unconsid. signal: dump)	//these aren't used
//const double noiseInter = 1.0e-15; // (max. unconsid. signal: int.)					  //same

//const double almostZero = 1.0e-50;	  //these aren't used

   double AMP, PER;
   double konstForce=0;
   double randForceVal=0;
   int cursign=-1;

//************* Function at end of code called 	addForcesomethignsomething... upgrades this
 double SmallInitialVelocity = 0.0;       // (initial v /small/ ptle (mm/musec))
 double LargeInitialVelocity = 0.00;//-1*0.149*0.001*1;//-.01      // initial v /large/ ptle (mm/musec))
//************

double epsilon = 1; // ((1 - restitution factor) all ptles)
int EP = 100; //change this variable... do not change the one above...

//**********Compression - set the initial overlap between the wall and the large particle
//from this a force balance is made... and all the particles are at rest (except minor vibrations
//due to inprecision in the decimals cause some crazy very small energy noise)
double InitialOverlapLargeandWall=0.00;  
double loadingForce=0;
//******************

//************dont alter these...
double ChainLength;
double ForceInChain;
//************************
	int count;
 double Force=0;
 double storeForce=Force;
//Given any single particle... it has the following:
struct particle
{

	double relativeLocation;
	double currentVelocity;
	double currentAccel;
	double radius;
	double currentKE;
	double mass;
	double absolutePosition;
};

//now we have an array,(the chain) of the above particles
particle* Chain;


double* deltas;//the overlaps
double* smalla;			  //the small-a in the calcs...
double* overbefore;  //stuff for the algorithm
double pot;
double* energy_loss;
double otime = 0.0;
double* forceBefore;

//these variables keep track of the max velocities... going in the direction of the initial impulse...
//(so ripples going back that are of larger velocity are ignored...)... this is easily alterable in the code
//the time variables correspond to what time that large velocity occurred.
double maxOverlapLargeParticle=0,maxOverlapSmallParticle=0, maxOLargeTime = 0, maxOSmallTime = 0;
double maxVelocityLargeParticle=0,maxVelocitySmallParticle=0, maxVLargeTime = 0, maxVSmallTime = 0;



std::ofstream KEouts, ETot,Vouts,XRelOuts, ForceAll, AppliedForce;  //file streams



void initializeAllParticles(); //set up the chain
void makeReadmeFile();		   //make a readme file
void ChainRelPositions();	   //set up more stuff
void CompressChain(double &force, double &length); //apply a compression

void velocityVerletStep();
void computeAccelerations();


double addForceDueToWave(double ttime);		//new addition... cerca 8/30/03

void spitKEs(double ttime);	 //spit to the file the KEs of the particles at a given time
void spitVs(double ttime);	 //similar but velocities
void spitXs(double ttime);
void spitFapp(double ttime,double force);
double absolute(double val);
double sign(double val);

void loadFromFile();

//more of these output files can be added... but if its not needed, why output and slow the program down?

double recordTotE;


double TimeMin=0;
double TimeMax=nsteps*dt;


std::string fileName="_"; 
bool FF=true, KK=true, VV=true, XX=true, TE=true,FA=false;


int topGrain=nptles;
int bottomGrain=1;

struct deltafunc
{
	int grain;
	double time;
	double addedVelocity;
};

deltafunc* deltaVel;
int deltafunclength=0;

struct drivingforces
{
	char type;
	double amplitude,frequency,phase;
};

drivingforces* sourceforce;
int sourceforcelength =0;

int main()
{

srand((unsigned)(time(0))); 

	std::string FileName;
	std::ostringstream *buffer;
	
	loadFromFile();
						
	ForceSpits = timeSpits;
 realspit = int(timeSpits/dt);
 realFspit = int(ForceSpits/dt);

	{std::cout<<"October 13, 2008 Version 1"<<std::endl;
		
			 std::cout<<"number of grains: "<<nptles<<std::endl;
			 std::cout<<"nsteps: "<<std::setprecision(8)<<nsteps<<std::endl;
			 std::cout<<"dt: "<<dt<<std::endl;
			 std::cout<<"Exponential: "<<xn<<std::endl;
			 std::cout<<"q: "<<q<<std::endl;
			 std::cout<<"epsilon: "<<epsilon<<std::endl;


			 if(KK)
			 {		 
				buffer = new std::ostringstream();
 
				(*buffer) << "KE"<<fileName<<".dat";
				FileName = buffer->str();

				KEouts.open(FileName.c_str()); 
 				KEouts.precision(DEFAULTPRECISION);
			 }

			 if(TE)
			 {
	buffer = new std::ostringstream();
	(*buffer)<<"ETot"<<fileName<<".dat";
	FileName = buffer->str();

	ETot.open(FileName.c_str());
	ETot.precision(DEFAULTPRECISION);
			 }
if(VV)
{
	 buffer = new std::ostringstream();
    (*buffer) << "VEL"<<fileName<< ".dat";
    FileName = buffer->str();

    Vouts.open(FileName.c_str()); 
 	Vouts.precision(DEFAULTPRECISION);
}
if(XX)
{
	 buffer = new std::ostringstream();
    (*buffer) << "XRel"<<fileName<< ".dat";
    FileName = buffer->str();

    XRelOuts.open(FileName.c_str()); 
 	XRelOuts.precision(DEFAULTPRECISION);
}

	

if(FF)
{
	buffer = new std::ostringstream();
    (*buffer) << "Force"<<fileName<< ".dat";
    FileName = buffer->str();

    ForceAll.open(FileName.c_str()); 
 	ForceAll.precision(DEFAULTPRECISION);

}
if(FA)
{
	buffer = new std::ostringstream();
    (*buffer) << "DrivingForce"<<fileName<< ".dat";
    FileName = buffer->str();

    AppliedForce.open(FileName.c_str()); 
 	AppliedForce.precision(DEFAULTPRECISION);

}
	initializeAllParticles();

	ChainRelPositions();

	makeReadmeFile();


	computeAccelerations();
	 Chain[nptles-1].currentAccel-=addForceDueToWave(otime);

	  std::cout<<Chain[nptles-1].currentAccel<<std::endl;

	for(unsigned long long int	lcv = 0; lcv<nsteps; ++ lcv)
		{
		
		count = lcv;
		otime = lcv * dt;

		 velocityVerletStep();
	
		 
		for(int lcount=0;lcount<deltafunclength;++lcount)
		{
			if(deltaVel[lcount].time-otime<dt)
			{
				Chain[nptles-deltaVel[lcount].grain].currentVelocity+=deltaVel[lcount].addedVelocity ;
				deltaVel[lcount].addedVelocity=0;
			}
		}
		 if(lcv%realspit==0&&otime>=TimeMin&&otime<= TimeMax)
		 { 
			if(KK||TE)
				spitKEs(otime);
			if(VV)
				spitVs(otime);
			if(FF||XX||FA)
				spitXs(otime);
		 }
		
		
		if((-1*Chain[0].currentVelocity)>maxVelocitySmallParticle)	//Added 6/7/03
		{ 
			maxVelocitySmallParticle=(Chain[0].currentVelocity*-1);
			maxVSmallTime = otime;
		} //
		if((-1*Chain[nptles-1].currentVelocity)>maxVelocityLargeParticle)		  //
		{ 
			maxVelocityLargeParticle=(Chain[nptles-1].currentVelocity*-1);	
  			maxVLargeTime = otime;
		}	//

		if((overbefore[0])>maxOverlapSmallParticle)	//Added 6/7/03
		{ 
			maxOverlapSmallParticle=(overbefore[0]);
			maxOSmallTime = otime-dt;
		} //
		if(otime<100)
		if((overbefore[nptles-1])>maxOverlapLargeParticle)		  //
		{ 
			maxOverlapLargeParticle=(overbefore[nptles-1]);	
  			maxOLargeTime = otime-dt;
		}	//

		}

if(VV)
	Vouts.close();
	
if(KK)	
KEouts.close();
if(TE)
ETot.close();  
	
if(XX)
XRelOuts.close();

	
if(FF)
ForceAll.close();
if(FA)
AppliedForce.close();
	}

	  return 0;
}

void initializeAllParticles()
{

Chain=(particle *) malloc(nptles*sizeof(particle));

forceBefore = (double *)malloc((nptles+1)*sizeof(double));
deltas= (double*)malloc((nptles+1)*sizeof(double));///the overlaps
smalla= (double*)malloc((nptles+1)*sizeof(double));			  //the small-a in the calcs...
overbefore = (double *)malloc((nptles+1)*sizeof(double));  //stuff for the algorithm
energy_loss = (double*)malloc((nptles+1)*sizeof(double));


	Chain[nptles-1].radius = rlarge;
	
	maxVelocityLargeParticle=0;
	maxVelocitySmallParticle=0;
	maxVLargeTime = 0;
	maxVSmallTime = 0;

	maxOverlapLargeParticle=0;
	maxOverlapSmallParticle=0;
	maxOLargeTime = 0;
	maxOSmallTime = 0;

 for(int lcv = nptles-1;lcv>=0;lcv--)
 {
	Chain[lcv].relativeLocation=0;
	Chain[lcv].currentVelocity=0;
	Chain[lcv].currentAccel=0;
	Chain[lcv].currentKE=0;
	if(lcv!=nptles-1)
		Chain[lcv].radius = ((100.0-q)/100.0)*Chain[lcv+1].radius; 
 
	Chain[lcv].mass = Chain[lcv].radius*Chain[lcv].radius*Chain[lcv].radius*(4.0/3.0)*PI*rho;
 
	deltas[lcv]=0.0;
	smalla[lcv]=0.0;
	overbefore[lcv]=0.0;
	energy_loss[lcv]=0.0;
	forceBefore[lcv] = 0.0;
	Chain[lcv].absolutePosition = 0.0;
 }

 Chain[nptles-1].currentVelocity = LargeInitialVelocity;
 Chain[0].currentVelocity = SmallInitialVelocity;
 
 deltas[nptles]=0.0;
 smalla[nptles]=0.0;
 overbefore[nptles]=0.0;
 energy_loss[nptles]=0.0;
 
 forceBefore[nptles]=0.0;


}

void makeReadmeFile()
{
   std::string FileName;
   std::ostringstream *buffer;
   std::ofstream out_file;

	  buffer = new std::ostringstream();
	
    (*buffer) << "ReadMe"<<fileName<< ".dat";
    FileName = buffer->str();

    out_file.open(FileName.c_str()); 
 	 out_file.precision(DEFAULTPRECISION);
	out_file<<"Number of Particles: "<<nptles<<std::endl;
	out_file<<"potential exponential: "<<xn<<std::endl;
	out_file<<"rho: "<<rho<<" (mg/mm^3)"<<std::endl;
	out_file<<"D: "<<D<<" (mm^2/kN)"<<std::endl;
	out_file<<"Tapering percent: "<<q<<std::endl;
	out_file<<"Restitution: "<<1-epsilon<<std::endl;
	out_file<<"Epsilon (1-w): "<<epsilon<<std::endl;
	out_file<<"Radius of Large Particle: "<<rlarge<<" (mm)"<<std::endl;
	out_file<<"Radii of all Particles: (mm)"<<std::endl;
	for(int lcv=0; lcv<nptles;++lcv)
		out_file<<Chain[lcv].radius<<'\t';
	out_file<<std::endl;
	out_file<<"Masses of all Particles: (mg)"<<std::endl;
	for(int lcv2=0; lcv2<nptles;++lcv2)
		out_file<<Chain[lcv2].mass<<'\t';
	out_file<<std::endl;

	out_file<<"dt: "<<dt<<"musec"<<std::endl;
	out_file<<"number of those steps: "<<nsteps<<std::endl;
	out_file<<"Length of run: "<<dt*nsteps<<" (musec)"<<std::endl;

	out_file<<"Initial Velocity of Large Particle: "<<-LargeInitialVelocity<<" (mm/musec)"<<std::endl;
	out_file<<"Initial Velocity of Small Particle: "<<-SmallInitialVelocity<<" (mm/musec)"<<std::endl;
		if(rightWall)
	out_file<<"Right Wall"<<std::endl;
		else
			out_file<<"No Right Wall"<<std::endl;
		if(leftWall)
	out_file<<"Left Wall"<<std::endl; 
		else
			out_file<<"No Left Wall"<<std::endl;
if(loadingForce>0)
out_file<<"Initial Loading of Chain "<<loadingForce<<" (kN)"<<std::endl;
	
		out_file<<"Files use these conventions:"<<std::endl;
		out_file<<"Time in mu-sec"<<std::endl;
		out_file<<"Force in kN"<<std::endl;
		out_file<<"Energy in J"<<std::endl;
		out_file<<"Distances in mm"<<std::endl;
		out_file<<"Velocities in mm/musec"<<std::endl;
		out_file<<"Density in mg/mm^3"<<std::endl;


	out_file.close();
	delete buffer;


}

void ChainRelPositions()
{
 
	double force = 0;
	double length = 0;
	
		if(loadingForce>0)	//if(InitialOverlapLargeandWall!=0) old line, replaced 10/13/08
	{	
			smalla[0] = (1 / (xn * D)) * (sqrt(Chain[0].radius));
			smalla[nptles] =  (1 / (xn * D)) * (sqrt(Chain[nptles-1].radius));//in compresschain i set these to 0 if there are no walls

		 for (int i = 1; i < nptles; i++)
				smalla[i] = (1 / (xn * D)) * (sqrt((Chain[i].radius*Chain[i-1].radius)/(Chain[i].radius+Chain[i-1].radius)));
 
		CompressChain(force, length);



	}
	else
	 {	  
		 for(int lcv = 0; lcv<nptles; ++lcv)
		 { 
			 length = length + 2 * Chain[lcv].radius;
		 	 Chain[lcv].absolutePosition = length - Chain[lcv].radius;
		 }
	 	 ChainLength = length;

if(rightWall)
		 smalla[0] = (1.0 / (xn * D)) * (sqrt(Chain[0].radius));
else
		smalla[0] = 0;
if(leftWall)
		 smalla[nptles] =  (1.0 / (xn* D)) * (sqrt(Chain[nptles-1].radius));
else
		smalla[nptles] = 0;
		 for (int i = 1; i < nptles; i++)
				smalla[i] = (1.0 / (xn * D)) * (sqrt((Chain[i].radius*Chain[i-1].radius)/(Chain[i].radius+Chain[i-1].radius)));
 
	 }
	   ForceInChain = force;
		 std::cout<<"Force is: "<<force<<"kN"<<std::endl;
		 std::cout<<"Length is: "<<length<<"mm"<<std::endl;


}

void CompressChain(double &force, double &length)
{
	

		std::cout<<"loading chain..."<<std::endl;
		InitialOverlapLargeandWall = pow((loadingForce/smalla[nptles]*1/xn),(1/(xn-1)));

deltas[nptles]=InitialOverlapLargeandWall;//set up initial overlap for beginning of recursion/loop



 force = xn * smalla[nptles] *pow(deltas[nptles],(xn-1.0));
  
 
for(int lcv = nptles-1; lcv>=0;--lcv)
{	deltas[lcv] = pow((force/smalla[lcv]*1/xn),(1/(xn-1)));

}
if(!leftWall)
{	deltas[nptles]=0;
	smalla[nptles]=0;
} //added these ifs in b/c if there is no wall but the chain is precompressed
if(!rightWall)			//these overlaps must be 0, ACS 10/13/08
{	deltas[0]=0;
	smalla[0]=0;
}

for(int lcv2 = 0; lcv2<nptles; ++lcv2)
{ 
	length = length + 2*Chain[lcv2].radius - deltas[lcv2];
	Chain[lcv2].absolutePosition = length - Chain[lcv2].radius;
}
length = length - deltas[nptles];

ChainLength = length;

 for(int lcv3=0; lcv3<=nptles; ++lcv3)
	 overbefore[lcv3] = deltas[lcv3];



}



void velocityVerletStep()
{


    for (int j = 0; j < nptles; j++) 
	{
	
	Chain[j].relativeLocation += Chain[j].currentVelocity * dt + 0.5 * Chain[j].currentAccel * dt*dt;
    Chain[j].currentVelocity += 0.5 * Chain[j].currentAccel * dt;
	}

  computeAccelerations();  
  Chain[nptles-1].currentAccel-=addForceDueToWave(otime);
  
  for (int j2 = 0; j2 < nptles; j2++) 
    Chain[j2].currentVelocity += 0.5 * Chain[j2].currentAccel * dt;


}




void computeAccelerations()
{

	double over,overnm1,forceBetw,forceFactor,moved_way,energy_comp,energy_deco,forceSmall,forceLarge;

	for (int i = 0; i < nptles; i++) // zeroing all acc in every call
    Chain[i].currentAccel = 0.0;       

	pot = 0.0;
  /******* potential/force between neighboring ptles *************/
  for (int i2 = 0; i2 < nptles-1; i2++)// 
  {
    if ((-1*Chain[i2].relativeLocation + Chain[i2+1].relativeLocation - deltas[i2+1])<=0) //swap del	had -delta[i2]>
	{//added in deltas 5/27                  // only when overlap
      over = Chain[i2].relativeLocation + deltas[i2+1] - Chain[i2+1].relativeLocation;  //added in deltas 5/27 swap del
      overnm1 = pow(over, (xn - 1.0));
      pot += over * overnm1 * smalla[i2+1];
      forceBetw = smalla[i2+1] * xn * overnm1;

      if (overbefore[i2+1] < over)         // when compressing
        forceFactor = 1.0;
      else forceFactor = epsilon;         // when decompressing

      forceBetw *= forceFactor;
   
      Chain[i2].currentAccel -= forceBetw;           // sign(-): towards smaller x				swapped
      Chain[i2+1].currentAccel += forceBetw;         // sign(+): towards larger x		 swapped

	 forceBefore[i2+1] = forceBetw;

	

      /****** calculate energy loss between ptles ********/
      if (overbefore[i2+1] < over) {  // compressing
        moved_way = over - overbefore[i2+1];
        energy_comp = forceBetw * moved_way;
        energy_loss[i2+1] += energy_comp; // whenever loading
      }
      if (overbefore[i2+1] > over) {// deco.; (don't mind '='-case)
		moved_way = overbefore[i2+1] - over;
        energy_deco = forceBetw * moved_way;
        energy_loss[i2+1] -= energy_deco; // whenever unloading
      }				 
      /*****************************************************/

      overbefore[i2+1] = over;          // update for next timestep 
    }
    else
	{
		overbefore[i2+1] = 0.0;//deltas[i2+1];//switched from 0.0 to deltas[i2+1] 7/16/03
		//pot += deltas[i2+1]	* pow(deltas[i2+1], (xn - 1.0)) * smalla[i2+1];		 //added 7/16/03
	}          // reset when no overlap
  }

  /** pot./force between fixed wall (small, x=0) <-> small ptle **/	  //THIS IS REALLY LARGE PARTICLE?
  if ((-1*Chain[0].relativeLocation + deltas[0])>=0 ) {	 //swap deltas
    over = -1*Chain[0].relativeLocation + deltas[0];	//swap deltas 
    overnm1 = pow(over, (xn - 1.0));
    pot += over * overnm1 * smalla[0];
    forceSmall = smalla[0] * xn * overnm1;

    if (overbefore[0] < over)
      forceFactor = 1.0;
    else forceFactor = epsilon;

    forceSmall *= forceFactor;


	if(rightWall)
    Chain[0].currentAccel += forceSmall;  //***********************swap -+

forceBefore[0] = forceSmall;
    /****** calculate energy loss at wall (small) ********/
    if (overbefore[0] < over) {  // compressing
      moved_way = over - overbefore[0];
      energy_comp = forceSmall * moved_way;
      energy_loss[0] += energy_comp;
    }
    if (overbefore[0] > over) {
      moved_way = overbefore[0] - over;
      energy_deco = forceSmall * moved_way;
      energy_loss[0] -= energy_deco;
    }
    /*****************************************************/

    overbefore[0] = over;
  }
  else 
  {  
	  overbefore[0] = 0.0;//deltas[0];//0.0; previously just 0.0  7/16
		//pot += deltas[0]	* pow(deltas[0], (xn - 1.0)) * smalla[0]; //added 7/16
  }
  /*** pot./force between fixed wall (large) <-> large ptle ******/
  if ((1*Chain[nptles-1].relativeLocation + deltas[nptles] )>=0) {  //Added deltas 5/27 ACS  //swap delta 
    over = 1*Chain[nptles-1].relativeLocation + deltas[nptles];	//Added 5/27 ACS	//swap delta  
    overnm1 = pow(over, (xn - 1.0));
    pot += over * overnm1 * smalla[nptles];
    forceLarge = smalla[nptles] * xn * overnm1;

    
    if (overbefore[nptles] < over)
      forceFactor = 1.0;
    else forceFactor = epsilon;

    forceLarge *= forceFactor;
	
	forceBefore[nptles] = forceLarge;

if(leftWall)
    Chain[nptles-1].currentAccel -= forceLarge;	 //**************swap +-
//else
	//Chain[nptles-1].currentAccel -= 0;

    /****** calculate energy loss at wall (large) ********/
    if (overbefore[nptles] < over) {  // compressing
      moved_way = over - overbefore[nptles];
      energy_comp = forceLarge * moved_way;
      energy_loss[nptles] += energy_comp; 
    }
    if (overbefore[nptles] > over) {
      moved_way = overbefore[nptles] - over;
      energy_deco = forceLarge * moved_way;
      energy_loss[nptles] -= energy_deco;
    }
    /*****************************************************/

    overbefore[nptles] = over;
  }
  else 
  {
	  overbefore[nptles] = 0.0;//deltas[nptles];//0.0; changed by ACS 7/16
//	pot += deltas[nptles]	* pow(deltas[nptles], (xn - 1.0)) * smalla[nptles];	//added 7/16
  
  }
  /***** real dim of acc: division by mass **********/
  for (int i3 = 0; i3 < nptles; i3++)
	  Chain[i3].currentAccel = Chain[i3].currentAccel / Chain[i3].mass;
  
  }



void spitKEs(double ttime)
{	
	double totKE=0;
	double tV;
	if(KK)
	KEouts<<(ttime);
	if(TE)
		ETot<<(ttime);


	for(int lcv = 0; lcv<nptles;++lcv)
	{	
				tV = Chain[lcv].currentVelocity;
	totKE = totKE +  (tV*tV*0.5*Chain[lcv].mass);
	}


	if(KK)
for(int lcv = bottomGrain; lcv<=topGrain;++lcv)
{
	tV = Chain[nptles-lcv].currentVelocity;
	KEouts<<'\t'<<(tV*tV*0.5*Chain[nptles-lcv].mass);
}

if(KK)
KEouts<<std::endl;

	
	recordTotE = pot+totKE;
	if(TE)
	ETot<<'\t'<<pot<<'\t'<<totKE<<'\t'<<(pot+totKE)<<std::endl;

}
 void spitVs(double ttime)
{	
	
	double tV;
	if(VV)
	Vouts<<(ttime);
	for(int lcv = bottomGrain; lcv<=topGrain;++lcv)
	{	

				tV = -Chain[nptles-lcv].currentVelocity;
				

				if(VV)
				Vouts<<'\t'<<(tV);	
	}
	if(VV)
	Vouts<<std::endl;



}
 void spitXs(double ttime)
{	
	
	double tX;

	if(FA)
		spitFapp(ttime,storeForce);
	//ForceFile<<(ttime);
	if(XX)
	XRelOuts<<(ttime);
	if(FF)
	ForceAll<<(ttime);

if(FF&&leftWall&&bottomGrain==1)
ForceAll<<'\t'<<forceBefore[nptles];
else if(FF&&bottomGrain!=1)
ForceAll<<'\t'<<forceBefore[nptles-bottomGrain+1];
	for(int lcv = bottomGrain; lcv<topGrain;++lcv)
	{	

				
					tX = Chain[nptles-lcv].relativeLocation;
					if(FF)
						//if(bottomGrain!=1||leftWall==true)
							ForceAll<<'\t'<<forceBefore[nptles-lcv];//(smalla[lcv] * xn * pow(overbefore[lcv], (xn - 1.0))* epsilon);
						
					if(XX)
					XRelOuts<<'\t'<<-(tX);	
	}
//right wall corresponds to grain 0
	
	if(FF&&rightWall)
	ForceAll<<'\t'<<forceBefore[nptles-topGrain];//(smalla[nptles] * xn * pow(overbefore[nptles], (xn - 1.0))* epsilon);
	
	if(XX)
XRelOuts<<'\t'<<-Chain[nptles-topGrain].relativeLocation;
	
if(FF)
	ForceAll<<std::endl;
if(XX)
	XRelOuts<<std::endl;
	
}
 void spitFapp(double ttime,double force)
 {
	
	AppliedForce<<(ttime)<<'\t'<<force<<std::endl;

 }

 //This function adds an acceleration to the large particle as you determine you want to do it..
 //it is given a force you provide, and simply divides out the mass of the particle...

//that could easily be changed...


double addForceDueToWave(double ttime)
{
	Force = konstForce;
	
	for(int c=0;c<sourceforcelength;++c)
	{
	//		std::cout<<sourceforce[c].type<<'\t'<<sourceforce[c].amplitude<<'\t'<<sourceforce[c].frequency<<'\t'<<sourceforce[c].phase<<'\t'<<std::endl;
		switch(sourceforce[c].type)
		{
				case 's':
					{	
						Force+=sourceforce[c].amplitude*sin(sourceforce[c].frequency*ttime -sourceforce[c].phase);
						break;}
				case 'c':
					{Force+=sourceforce[c].amplitude*cos(sourceforce[c].frequency*ttime -sourceforce[c].phase);break;}
				case 't':
					{
						Force+=sourceforce[c].amplitude*(absolute(4.0*((ttime-sourceforce[c].phase)*sourceforce[c].frequency/(2.0*PI)-floor(sourceforce[c].frequency/(2.0*PI)*(ttime-sourceforce[c].phase)+0.5)))-1.0);
						break;}
				case 'w':
					{
						Force+=sourceforce[c].amplitude*2.0*((ttime-sourceforce[c].phase)*sourceforce[c].frequency/(2.0*PI)-floor(sourceforce[c].frequency/(2.0*PI)*(ttime-sourceforce[c].phase)+0.5));
						break;}
				case 'q':
					{
						Force+=sourceforce[c].amplitude*sign(sin(sourceforce[c].frequency*ttime-sourceforce[c].phase));
						break;}
					case 'r':
					{
						if(cursign==sign(sin(sourceforce[c].phase/2*ttime-PI)))
						{
							randForceVal = double(2*rand()/(double(RAND_MAX)+1.0)-1.0)*(sourceforce[c].amplitude-sourceforce[c].frequency)/2.0+(sourceforce[c].amplitude+sourceforce[c].frequency)/2.0;
							cursign=cursign*-1;
						}
						Force+=randForceVal;
						break;}
				default :
					std::cout<<"unknown force"<<std::endl;break;
	
	
		}	
	}			
	/*if(randomForce)
	{
	Force += double(rand()%10)/10
	}*/	
		//Force = 0;//+sin(ttime*PI/100);
		//initial force 0...
	  //if(ttime<358)						//while time is under 50... apply the following force
//Force = (AMP/1000)*pow((sin(pow(sin(ttime/10000),2)*(ttime*(PI*(PER/100))))),2);
	//Force = 0.5*AMP;
	//else
	//Force = 0;	//otherwise...its 0...
//Force = 0.05*cos(ttime*4)+0.3*sin(ttime/10)+5*sin(ttime/30);
	//Force = 4 + sin(ttime) + 2*cos(ttime*PI) + 3*sin(ttime/2);
	  //This function could easily be applied cotinuously... by removing the if statement
	  //it could be switched up any way you want...basically you have the freedom to define
	  //a piece wise function, and they can be anything you want, ... so its kinda nice...

	  //just be careful not to spit too much energy into the system, otherwise you'll
	  //get some results that may not be accurate due to over compression of the spheres...
storeForce= Force;
//	if((int(ttime/dt) % int(ForceSpits/dt))==0)
	/*if((count)%realFspit == 0)
	{
		spitFapp(ttime,(Force));
		
		//	Force = double(rand()%10)/10;
	
	}*/
	 return Force/Chain[nptles-1].mass;
}

double absolute(double val)
{
if(val<0)
return val*-1;
else
return val;
}
double sign(double val)
{
if(val<0)
return -1;
else
return 1;
}

void loadFromFile()
{
	bool newTimeRange=false;
	bool Specified=false;
	bool newGrains=false;
	bool grainsSpecified=false;
	std::string param;
	char cc;
	double parval;
	
	std::ifstream in_file("parameters.txt");

	in_file>>param;

	while(!in_file.eof()&&!in_file.fail())
	{
		
/*N: 20
dt: 0.002
nsteps: 20000
q: 5
w: 0.01
rho: 4.42
D: 0.01206
Wall: 11
rlarge: 5
timeSpits: 0.5
SmallInitV: 0
LargeInitV: -0.01*/

		for(unsigned int i =0;i<param.length();++i)
			param[i]=tolower(param[i]);

		if(param!="files:"&&param!="filename:"&&param!="addforce:")
		in_file>>parval;

		if(param=="n:")
		{
			nptles=int(parval);
			newGrains=true;
		}
		if(param=="dt:")
		{
			dt=parval;
			newTimeRange = true;
		}
		
		if(param=="nsteps:")
		{
			nsteps= (unsigned long long int) (parval); // ( ) added by YT. 06202009 
			newTimeRange = true;

		}
		
		if(param=="q:")
			q=parval;
		if(param=="w:")
			epsilon=1-parval;//std::cout<<"wfound"<<parval<<std::endl;
		if(param=="rho:")
			rho=parval;
		if(param=="precision:")
			DEFAULTPRECISION=int(parval);
		if(param =="exponential:")
			xn = parval;
		if(param =="preload:")
			loadingForce = parval;//kN
		

		if(param=="d:")
			D=parval;
		
		if(param=="wall:")
		{
			if(int(parval)==0)
			{
				leftWall=rightWall=false;
			}
			if(int(parval)==11)
			{
				leftWall=rightWall=true;
			}
			if(int(parval)==10)
			{
				leftWall=true;
				rightWall = false;
			}
			if(int(parval)==1)
			{
				leftWall=false;
				rightWall = true;
			}
		//	std::cout<<"Wallfound"<<parval<<std::endl;
		
		
		}
		
		if(param=="rlarge")
			rlarge=parval;
		if(param=="timespits:")
			timeSpits = parval;
		if(param=="smallinitv:")
			SmallInitialVelocity = -1*parval;
		if(param=="largeinitv:")
			LargeInitialVelocity = -1*parval;
		
		if(param=="deltav:")
		{
			 deltafunc* temp;
			 temp=(deltafunc*)malloc((deltafunclength+1)*sizeof(deltafunc));
			 
			 for(int c=0;c<deltafunclength;++c)
			 { temp[c].addedVelocity = deltaVel[c].addedVelocity;
			 temp[c].grain = deltaVel[c].grain;
			 temp[c].time = deltaVel[c].time;
			 }
			
			
			temp[deltafunclength].grain = int(parval);
			in_file>>parval;
			temp[deltafunclength].time = (parval);
			in_file>>parval;
			temp[deltafunclength].addedVelocity = -(parval);

			deltafunclength++;

			deltaVel = temp;
			temp = NULL;
			
		}

		if(param=="timemin:")
			TimeMin = parval;
		if(param=="timemax:")
		{	TimeMax = parval;
			Specified = true;
		}

		if(param == "files:")
		{

			
			FF=false;
			TE=false;
			KK=false;
			XX=false;
			VV=false;
			FA=false;
		
			
			if(!in_file.eof()&&!in_file.fail())
			in_file>>param;
			
			for(unsigned int i =0;i<param.length();++i)
				{
					param[i]=tolower(param[i]);
					
				if(param[i]=='f')
					FF=true;
				if(param[i]=='t')
					TE=true;
				if(param[i]=='k')
					KK=true;
				if(param[i]=='x')
					XX=true;
				if(param[i]=='v')
					VV=true;
				if(param[i]=='a')
					FA=true;
				}
		}
if(param == "addforce:")
		{
			cc ='d';//default
			if(!in_file.eof()&&!in_file.fail())
			in_file>>cc;
			
			//std::cout<<cc<<std::endl;

			cc=tolower(cc);
			if(cc=='k'){
				if(!in_file.eof()&&!in_file.fail())
				{
					parval =0;
					in_file>>parval;
				konstForce +=parval;
				}}
			else
			{
			drivingforces* temp;
			 temp=(drivingforces*)malloc((sourceforcelength+1)*sizeof(drivingforces));
			  
			 for(int c=0;c<sourceforcelength;++c)
			 { temp[c].type = sourceforce[c].type;
			 temp[c].amplitude = sourceforce[c].amplitude;
			 temp[c].frequency = sourceforce[c].frequency ;
			 temp[c].phase = sourceforce[c].phase;
			 }
			
			
			temp[sourceforcelength].type = cc;
			in_file>>parval;
			temp[sourceforcelength].amplitude = (parval);
			in_file>>parval;
			temp[sourceforcelength].frequency = (parval);
			in_file>>parval;
			temp[sourceforcelength].phase = (parval);
			sourceforcelength++;

			sourceforce = temp;
			temp = NULL;
				
			
			}
		}
	if(param == "grains:")
		{
				
			bottomGrain= (int(parval));
			in_file>>parval;
			topGrain= (int(parval));
			grainsSpecified=true;
			
		}

	if(param == "filename:")
	{
		in_file>>param;
		fileName = "_" + param;
	}

		in_file>>param;
	}
	in_file.close();


	
	if(newTimeRange&&!Specified)
		TimeMax = nsteps*dt;
	if(newGrains&&!grainsSpecified)
		topGrain = nptles;

	
}


