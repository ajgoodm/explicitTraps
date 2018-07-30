#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
using namespace std;

/*------ structures ------*/
struct trap
{
  double xPos;
  double yPos;
  double energy;
  bool isOccupied;
};

struct exciton
{
  int zoneX;
  int zoneY;
  double xPos;
  double yPos;
  bool isTrapped;
  bool exists;
  int trapIndex;
};

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
void resetExcitons(exciton a[], int nExcitonsToMake, int nExcitonsThatExist, double zoneSideLength, int nZonesPerSide);
bool stepOver(point x1, point x2, point P, double R);
bool collisionCheck(point x1, point x2, point P, double R);
void makeExciton( exciton excitonArray[], int excitonIndex, double FWHM, double alpha, double zoneSideLength, int nZonesPerSide);

int main(){
  srand(time(NULL));
  ofstream QYs;
  QYs.open("sto_QY.txt");

/* Trap Parameters */
  const double boxSideLength = 4e-6;                      // [m]
  const int    nZonesPerSide  = 100;
  const double zoneSideLength = boxSideLength/double(nZonesPerSide);                   // [m]
  const double simulationArea = pow(boxSideLength,2);     // [m^2]
  const double zoneArea = pow(zoneSideLength,2);          // [m^2]      35 ns lifetime conditions, density = 1e15
  const double trapDensity = 1e15;                        // [m^-2]     sto: 11, 5e13, R = 0.75
  const double kDetrapBase = 5e13;                        // [s^-1]     qtz: 6, 1e15, R = 1.0
  const int    numberOfTraps = floor(simulationArea*trapDensity);
  const int    numberOfTrapsPerZone = floor(zoneArea*trapDensity);
  double alpha = 11; double lowBound = 0.05; double highBound = .45;
  cout << numberOfTrapsPerZone << endl;
/*-----------------*/

/* Exciton Parameters */
  double Diffusivity = 5e-4;
  double R = 0.4e-9;
  double trappedExciton_R = 0.4e-9;
  double kRad = 1.0/(1200e-12);

  exciton excitonArray[10000];
/*--------------------*/

/* Simulation Parameters */
  double deltat = 1.0e-13;                                // [s]
  double step = sqrt(4*Diffusivity*deltat);
  double totalTime = 1e-4;

  int nExcitonsMade = 0;
  int nPhotonsOut;
  int numberOfSteps = floor(totalTime/deltat);
  int density = 50;
  int nExcitonsThatExist = 0;
  int nExcitonsToMake;
  double RgenQtz[11]={1.95e14, 5.98e14, 1.94e15, 6.02e15, 1.95e16, 5.9e16, 1.98e17, 5.97e17, 2.0e18, 5.94e18, 1.99e19};
  double RgenSap[14]={4.27e14, 7.22e14, 1.41e15, 3.56e15, 7.91e15, 1.48e16, 3.29e16, 7.87e16, 1.98e17, 4.09e17, 1.06e18, 2.73e18, 7.03e18, 1.9e19};
  double RgenSTO[13]={1.158e15, 4.0981e14, 3.1064e15, 7.35e15, 1.45e16, 4.24e16, 8.47e16, 1.98e17, 4.18e17,1.10e18, 2.1e18,5.13e18,1.026e19};
/*-----------------------*/

/* miscellaneous */
  const double pi = 3.14159265358979;
  int gg, ii, jj, kk, ll, mm;     // loop indices;
  double randomRoll;
  double Rgen;
  double spotFWHM;
  point start, end, P;
/*---------------*/
  trap arrayOfTraps[numberOfTrapsPerZone][100][100];
  initializeTrapArray(arrayOfTraps, numberOfTrapsPerZone, zoneSideLength, boxSideLength,
                      nZonesPerSide, alpha, lowBound, highBound);

  nExcitonsToMake = 2*density;
  nExcitonsThatExist = 0;
  nPhotonsOut=0;


  for(gg=0; gg<13; gg++)
  {
    Rgen = RgenSTO[gg];
    spotFWHM = 0.48;
    Rgen = Rgen * (spotFWHM/2.0) * (spotFWHM/2.0) * pi * 1.0e-8; // [1/s]
    nExcitonsMade=0;
    nPhotonsOut=0;

    for(ii=0; ii<numberOfSteps; ii++)
    {
      randomRoll = double(rand())/double(RAND_MAX);
      if(randomRoll<(deltat*Rgen))
      {
        makeExciton(excitonArray, nExcitonsThatExist, spotFWHM/1e6, alpha, zoneSideLength, nZonesPerSide);
        nExcitonsThatExist++;
        nExcitonsMade++;
      }
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
                      arrayOfTraps[kk][mm][ll].isOccupied=true;
                      // cout << "exciton #" << jj << " trapped at trap site #" << kk << endl;
                      break;
                    }
                    else if (collisionCheck(start, end, P, trappedExciton_R)&&(arrayOfTraps[kk][mm][ll].isOccupied))
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
            if(excitonArray[jj].exists&&(!excitonArray[jj].isTrapped))
            {
              randomRoll = (double(rand())/double(RAND_MAX));
              if(randomRoll<(1-exp(-deltat*kRad)))
              {
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
              arrayOfTraps[excitonArray[jj].trapIndex][excitonArray[jj].zoneX][excitonArray[jj].zoneY].isOccupied=false;
              // cout << "exciton #" << jj <<  " detrapped from site #" << excitonArray[jj].trapIndex << endl;
            }
          }
        }
      }
      defragArray(excitonArray, nExcitonsThatExist);
      if((nPhotonsOut>=20)&&(nExcitonsMade>=500))
      {break;}
    }
    cout << "Number of Excitations: " << nExcitonsMade << endl
		     << "Number of Photons out: " << nPhotonsOut << endl
		     << "QuantumYield: " << double(nPhotonsOut)/double(nExcitonsMade) << endl;
	  QYs  << double(nPhotonsOut)/double(nExcitonsMade) << endl;

  }
  QYs.close();
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
void resetExcitons(exciton a[], int nExcitonsToMake, int nExcitonsThatExist, double zoneSideLength, int nZonesPerSide){
  double FWHM = 0.8e-6;					// [m]
  random_device generator;
	normal_distribution<double> gaussianDist(0.0,FWHM/2.355);

	for(int ii=0; ii<nExcitonsToMake; ii++)
	{
		a[nExcitonsThatExist+ii].xPos = gaussianDist(generator); a[nExcitonsThatExist+ii].yPos = gaussianDist(generator);
    a[nExcitonsThatExist+ii].zoneX = int(floor(a[nExcitonsThatExist+ii].xPos/zoneSideLength))+nZonesPerSide/2;
    a[nExcitonsThatExist+ii].zoneY = int(floor(a[nExcitonsThatExist+ii].yPos/zoneSideLength))+nZonesPerSide/2;
    a[nExcitonsThatExist+ii].isTrapped=false; a[nExcitonsThatExist+ii].exists=true;
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
void makeExciton( exciton excitonArray[], int excitonIndex, double FWHM, double alpha, double zoneSideLength, int nZonesPerSide) {
	default_random_engine generator;
	default_random_engine expgenerator;
	normal_distribution<double> gaussianDist(0.0,FWHM/2.355);
	exponential_distribution<double> exp_distribution(alpha);

	excitonArray[excitonIndex].xPos = gaussianDist(generator); excitonArray[excitonIndex].yPos = gaussianDist(generator);
  excitonArray[excitonIndex].zoneX = int(floor(excitonArray[excitonIndex].xPos/zoneSideLength))+nZonesPerSide/2;
  excitonArray[excitonIndex].zoneY = int(floor(excitonArray[excitonIndex].yPos/zoneSideLength))+nZonesPerSide/2;
  excitonArray[excitonIndex].isTrapped=false; excitonArray[excitonIndex].exists=true;
	return;
}
