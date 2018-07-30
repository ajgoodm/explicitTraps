#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <string>
#include <sstream>
#include <unistd.h>
using namespace std;

/*------ structures ------*/
struct trap
{
  double xPos;
  double yPos;
  double energy;
  bool isOccupied;
};
// This represents a trap in MoS_2: it's physical location in
// 2D space (m), it's depth (in eV), and whether it's isOccupied
// by a trapped exciton.

struct exciton
{
  int zoneX;
  int zoneY;
  double xPos;
  double yPos;
  int parentPulse;
  int pulseTrapped;
  bool isTrapped;
  bool exists;
  bool hasTrapped;
  int trapIndex;
};
// This structure represents an exciton. It has a position in 2D
// space. To ease in exciton-exciton interaction checks, each
// exciton is also assigned to a zone, which is just a binning
// of the simulation box to reduce the number of pairwise
// interactions checked.
//
// There's also a pulse ID, which says which laser pulse formed
// the exciton, and bools that say whether the exciton is trapped,
// and whether it still exists. The trapIndex indexes the explicit
// trap the exciton is associated with.

struct point
{
  double xPos;
  double yPos;
};
/*------------------------*/

double drawTrapEnergy(double alpha, double lowBound, double highBound);
double randUniform(double minbound, double maxbound);
void initializeTrapArray(trap arrayOfTraps[][100][100], const int numberOfTrapsPerZone,
                         const double zoneSideLength, const double boxSideLength, const int nZonesPerSide,
                         double alpha, double lowBound, double highBound);
double distBetween(trap a, exciton b);
double distBetween(exciton a, exciton b);
double distBetween(point a, point b);
void defragArray(exciton excitonArray[], int nExcitons);
void resetExcitons(exciton a[], int nExcitonsToMake, int nExcitonsThatExist, double zoneSideLength, int nZonesPerSide, int nLaserPulses);
bool stepOver(point x1, point x2, point P, double R);
bool collisionCheck(point x1, point x2, point P, double R);

int main(){
  srand(time(NULL));


/* Trap Parameters */
  const double boxSideLength = 4e-6;                      // [m] This is the simulation box size
  const int    nZonesPerSide  = 100;                      // This is the number of zones per box size
                                                          // This is done to reduce the number of pairwise
                                                          // interactions to check.
  const double zoneSideLength = boxSideLength/double(nZonesPerSide);                   // [m]
  const double simulationArea = pow(boxSideLength,2);     // [m^2]
  const double zoneArea = pow(zoneSideLength,2);          // [m^2]
  const double trapDensity = 1e15;                        // [m^-2] This is the area density of traps
  const double kDetrapBase = 5e13;                        // [s^-1] This is the base detrapping attempt frequency
  const int    numberOfTraps = floor(simulationArea*trapDensity);
  const int    numberOfTrapsPerZone = floor(zoneArea*trapDensity);
  double alpha = 10; double lowBound = 0.15; double highBound = .4; // Trap depths are drawn from an exponential normal_distribution
                                                                    // This parameterizea the distribution from which they're drawn.
                                                                    // They are on the interval [lowBound, highBound] from an exponentially
                                                                    // decaying distribution. This is derived from experiment.
  cout << numberOfTrapsPerZone << endl;
/*-----------------*/

/* Exciton Parameters */
  double Diffusivity = 5e-4;            // The exciton diffusivity m^2/s
  double R = 0.4e-9;                    // Exciton interaction radius
  double trappedExciton_R = 0.4e-9;     // Trapped exciton interaction radius
  double kRad = 1.0/(900e-12);          // exciton radiative rate
  double trapDecayTime = 500e-9;        // Nonradiative decay time constant
                                        // long-lived trapped excitons
  exciton excitonArray[10000];          // Initialized exciton array
/*--------------------*/

/* Simulation Parameters */
  double deltat = 1.0e-13;                    // Simulation time step [s]
  double step = sqrt(4*Diffusivity*deltat);   // Fixed distance step in diffusive
                                              // random walk.
  double totalTime = 20e-9;                   // Total time between pulses
  int totalCountGoal = 10000;
  int nPhotonsOut;
  int numberOfSteps = floor(totalTime/deltat);
  double density = 1.0;
  int nExcitonsThatExist = 0;
  int nExcitonsToMake;
/*-----------------------*/

/* miscellaneous */
  const double pi = 3.14159265358979;
  int ii, jj, kk, ll, mm;     // loop indices;
  int  nLaserPulses=0;
  double randomRoll;
  point start, end, P;
/*---------------*/
  trap arrayOfTraps[numberOfTrapsPerZone][100][100];  // This 2D array keeps track of the traps
                                                      // in each zone for book-keeping purposes
  initializeTrapArray(arrayOfTraps, numberOfTrapsPerZone, zoneSideLength, boxSideLength,
                      nZonesPerSide, alpha, lowBound, highBound);

  int nChildProcesses = 8;    // This allows the simulation to run multiple instances
                              // in parallel.
  pid_t pid;
  for(ii=0; ii<nChildProcesses; ii++)
  {
    pid = fork();
    if(pid>0)
    {
      // parent process
      cout << "parent process id: " << getpid() << endl;
    }
    else if(pid == 0)
    {
      //child process
      cout << "child process id: " << getpid() << " with parent: " << getppid() << endl;
      break;
    }
    else if(pid<0)
    {
      cout << "fork error" << endl;
    }
  }


  ofstream output;
  stringstream ss;
  ss << "1perum_" << to_string(getpid()) << ".txt";
  string filename = ss.str();

  output.open(filename);


  nExcitonsToMake = int(2.0*density);
  nExcitonsThatExist = 0;
  nPhotonsOut=0;

  do {
    defragArray(excitonArray, nExcitonsThatExist);
    resetExcitons(excitonArray, nExcitonsToMake, nExcitonsThatExist, zoneSideLength, nZonesPerSide, nLaserPulses);
    nExcitonsThatExist=nExcitonsThatExist+nExcitonsToMake;
    nLaserPulses++;

    // cout << "we've collected " << nPhotonsOut << " photons. There are " << nExcitonsThatExist << " total excitons." << endl;
    for(ii=0; ii<numberOfSteps; ii++)
    {
      for(jj=0; jj<nExcitonsThatExist; jj++)
      {
        if(excitonArray[jj].exists)
        {
          if(!excitonArray[jj].isTrapped)
          {
            // Exciton hops. Does it hit a trap or annihilate?
            start.xPos=excitonArray[jj].xPos; start.yPos=excitonArray[jj].yPos;
            excitonArray[jj].zoneX=int(floor(excitonArray[jj].xPos/zoneSideLength))+nZonesPerSide/2;
            excitonArray[jj].zoneY=int(floor(excitonArray[jj].yPos/zoneSideLength))+nZonesPerSide/2;
            randomRoll = 2*pi*(double(rand())/double(RAND_MAX));
						excitonArray[jj].xPos = excitonArray[jj].xPos + step*cos(randomRoll);
						excitonArray[jj].yPos = excitonArray[jj].yPos + step*sin(randomRoll);
            end.xPos=excitonArray[jj].xPos; end.yPos=excitonArray[jj].yPos;
            for(mm=max(0,excitonArray[jj].zoneX-2); mm<min(excitonArray[jj].zoneX+2,nZonesPerSide); mm++)
            {
              for(ll=max(0,excitonArray[jj].zoneY-2); ll<min(excitonArray[jj].zoneY+2,nZonesPerSide); ll++)
              {
                for(kk=0; kk<numberOfTrapsPerZone; kk++)
                {
                  if((abs(excitonArray[jj].xPos-arrayOfTraps[kk][mm][ll].xPos)<step)&&(abs(excitonArray[jj].yPos-arrayOfTraps[kk][mm][ll].yPos)<step))
                  {
                    P.xPos=arrayOfTraps[kk][mm][ll].xPos; P.yPos=arrayOfTraps[kk][mm][ll].yPos;
                    if(collisionCheck(start, end, P, R)&&(!arrayOfTraps[kk][mm][ll].isOccupied))
                    {
                      excitonArray[jj].isTrapped=true;
                      excitonArray[jj].trapIndex=kk;
                      excitonArray[jj].zoneX=mm;
                      excitonArray[jj].zoneY=ll;
                      excitonArray[jj].pulseTrapped=nLaserPulses;
                      arrayOfTraps[kk][mm][ll].isOccupied=true;
                      // cout << "exciton #" << jj << " trapped at trap site #" << kk << endl;
                      break;
                    }
                    else if ((collisionCheck(start, end, P, trappedExciton_R)&&(arrayOfTraps[kk][mm][ll].isOccupied)))
                    {
                      excitonArray[jj].exists=false;
                      nExcitonsThatExist--;

                      // cout << "exciton #" << jj << " annihilated at kk: " << kk << endl;
                      break;
                    }
                  }
                }
                if(excitonArray[jj].isTrapped||!excitonArray[jj].exists)
                {break;}
              }
              if(excitonArray[jj].isTrapped||!excitonArray[jj].exists)
              {break;}
            }

            // Does the free exciton decay?
            if((excitonArray[jj].exists&&(!excitonArray[jj].isTrapped))&&excitonArray[jj].hasTrapped)
            {
              randomRoll = (double(rand())/double(RAND_MAX));
  						if(randomRoll<(1-exp(-deltat*kRad)))
              {
                output << deltat*ii << ", " << excitonArray[jj].xPos << ", " << excitonArray[jj].yPos << endl;
                nPhotonsOut++;
                excitonArray[jj].exists=false;
                nExcitonsThatExist--;
                // cout << "exciton #" << jj << " decayed." << endl;
                break;
              }
            }
          }
          else // Does the trapped exciton detrap?
          {
            randomRoll=(double(rand())/double(RAND_MAX));
            if(randomRoll<(1-exp(-deltat*kDetrapBase*exp(-arrayOfTraps[excitonArray[jj].trapIndex][excitonArray[jj].zoneX][excitonArray[jj].zoneY].energy/(300*8.6e-5)))))
            {
              excitonArray[jj].isTrapped=false;
              randomRoll = 2*pi*(double(rand())/double(RAND_MAX));
              excitonArray[jj].xPos = excitonArray[jj].xPos + 1.2*R*cos(randomRoll);
  						excitonArray[jj].yPos = excitonArray[jj].yPos + 1.2*R*sin(randomRoll);
              excitonArray[jj].hasTrapped=true;
	            arrayOfTraps[excitonArray[jj].trapIndex][excitonArray[jj].zoneX][excitonArray[jj].zoneY].isOccupied=false;
              // cout << "exciton #" << jj <<  " detrapped from site #" << excitonArray[jj].trapIndex << endl;
            }
            // Does the trapped exciton decay radiatively?
            if((double(nLaserPulses-excitonArray[jj].pulseTrapped))*totalTime>trapDecayTime)
            {
              excitonArray[jj].exists=false;
              nExcitonsThatExist--;
              arrayOfTraps[excitonArray[jj].trapIndex][excitonArray[jj].zoneX][excitonArray[jj].zoneY].isOccupied=false;
            }
          }
        }
      }
    }


  } while(nPhotonsOut<totalCountGoal);


  output.close();
  return 0;
}

double drawTrapEnergy(double alpha, double lowBound, double highBound){
    int seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    exponential_distribution<double> distribution(alpha);

    double energy;
    do {
      energy=distribution(generator);
    }while((energy<lowBound)||(energy>highBound));

    return energy;

}
double randUniform(double minbound, double maxbound){
    double result = ((double(rand())/double(RAND_MAX))*(maxbound-minbound))+minbound;
    return result;
}
void initializeTrapArray(trap arrayOfTraps[][100][100], const int numberOfTrapsPerZone,
                         const double zoneSideLength, const double boxSideLength, const int nZonesPerSide,
                         double alpha, double lowBound, double highBound)
{
  int ii, jj, kk;
  for(ii=0; ii<nZonesPerSide; ii++)
  {
    for(jj=0; jj<nZonesPerSide; jj++)
    {
      for(kk=0; kk<numberOfTrapsPerZone; kk++)
      {
        arrayOfTraps[kk][ii][jj].xPos=randUniform(-boxSideLength/2.0+(double(ii)*zoneSideLength),-boxSideLength/2.0+((double(ii)+1)*zoneSideLength));
        arrayOfTraps[kk][ii][jj].yPos=randUniform(-boxSideLength/2.0+(double(jj)*zoneSideLength),-boxSideLength/2.0+((double(jj)+1)*zoneSideLength));
        arrayOfTraps[kk][ii][jj].energy=drawTrapEnergy(alpha, lowBound, highBound);
        arrayOfTraps[kk][ii][jj].isOccupied=false;
      }
    }
  }

  return;
}
double distBetween(trap a, exciton b)
{
  return sqrt(pow(a.xPos-b.xPos,2)+pow(a.yPos-b.yPos,2));
}
double distBetween(exciton a, exciton b)
{
  return sqrt(pow(a.xPos-b.xPos,2)+pow(a.yPos-b.yPos,2));
}
double distBetween(point a, point b)
{
  return sqrt(pow(a.xPos-b.xPos,2)+pow(a.yPos-b.yPos,2));
}
void defragArray(exciton excitonArray[], int nExcitons) {
	int jj;
	exciton temp;

	if (nExcitons!=0)
	{
		for(int ii=0; ii<nExcitons; ii++)
		{
			if(!excitonArray[ii].exists)
			{
				jj=ii+1;
				do {
					if(excitonArray[jj].exists)
					{
						temp=excitonArray[ii];
						excitonArray[ii]=excitonArray[jj];
						excitonArray[jj]=temp;
					}
					jj++;
				} while(!excitonArray[ii].exists);
			}
		}
	}
	return;
}
void resetExcitons(exciton a[], int nExcitonsToMake, int nExcitonsThatExist, double zoneSideLength, int nZonesPerSide, int nLaserPulses){
  double FWHM = 0.8e-6;					// [m]
  random_device generator;
	normal_distribution<double> gaussianDist(0.0,FWHM/2.355);

	for(int ii=0; ii<nExcitonsToMake; ii++)
	{
		a[nExcitonsThatExist+ii].xPos = gaussianDist(generator); a[nExcitonsThatExist+ii].yPos = gaussianDist(generator);
    a[nExcitonsThatExist+ii].zoneX = int(floor(a[nExcitonsThatExist+ii].xPos/zoneSideLength))+nZonesPerSide/2;
    a[nExcitonsThatExist+ii].zoneY = int(floor(a[nExcitonsThatExist+ii].yPos/zoneSideLength))+nZonesPerSide/2;
    a[nExcitonsThatExist+ii].isTrapped=false; a[nExcitonsThatExist+ii].exists=true;
    a[nExcitonsThatExist+ii].hasTrapped=false;
    a[nExcitonsThatExist+ii].parentPulse=nLaserPulses;
    a[nExcitonsThatExist+ii].pulseTrapped=0;
	}
	return;
}
bool stepOver(point x1, point x2, point P, double R) {
	double m=(x1.yPos-x2.yPos)/(x1.xPos-x2.xPos);
	point intersect;
		if(m==0.0)
		{
			intersect.xPos=P.xPos;
			intersect.yPos=x1.yPos;
		}
		else if(x1.xPos==x2.xPos)
		{
			intersect.xPos=x1.xPos;
			intersect.yPos=P.yPos;
		}
		else
		{
			intersect.xPos=(P.yPos-x2.yPos+(P.xPos/m)+(x2.xPos*m))/(m+(1/m));
			intersect.yPos=m*(intersect.xPos-x2.xPos)+x2.yPos;
		}
	if (distBetween(P, intersect) <= R)
	{return true;}
	else
	{return false;}
}
bool collisionCheck(point x1, point x2, point P, double R) {
	if ((distBetween(x1,P) < 2*R)||(distBetween(x2,P) < 2*R))
	{return true;}
	else if ((P.xPos < max(x1.xPos,x2.xPos))&&(P.xPos > min(x1.xPos,x2.xPos)))
	{
		if(stepOver(x1,x2,P,R))
		{return true;}
	}
	return false;
}
