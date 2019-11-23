#ifndef ADM_SOLVER
#define ADM_SOLVER
#if defined(__APPLE__) || defined(__linux__)
# define PREDEF_PLATFORM_UNIX
#endif
int METHOD = 1;
/*
List of Methods:

Each method is an unconditionally-stable diagonally-implicit Runge-Kutta method
with the same order as the number of the method. Butcher tableaus are listed as
follows for the various methods:

Method 1 (Backwards Euler):

1 | 1
--------
  | 1

This method has the best stability of all of them, however it is 1st-order in 1 step.
It proves to be very fast, even compared to the higher-order solvers.
***********************************************************************************
Method 2 (Implicit Midpoint Formula):

0.5 | 0.5
------------
    |  1

This method has the worst stability and tends to suffer from oscillations in the
step size that make it very slow. That problem can be fixed, but I haven't done so.
Still, this is a 2nd-order method in 1 step. Not bad.
**********************************************************************************
Method 3:

0.7887 | 0.7887
0.2113 | -.5774  0.7887
--------------------------
       |   0.5    0.5

A third-order method in two steps, with good stability and good speed. Not as good
as either the 1st or 4th-order methods in practice, at least as far as I've seen.
*********************************************************************************
Method 4:

 1  |   1
0.5 | -.75  1.25
 0  |   2    -3    1
 --------------------------
    |  1/6   2/3  1/6

Relatively stable for a 4th-order method in 3 steps, and very fast. That being
said, the 1st-order method turns out to be way better for this problem.
**********************************************************************************
Method 5:

  2  |   2
 -.5 |  -2  1.5
 -------------------
     |  .4   .6

Trying a 2nd-order method here that has the same stability properties as the
backwards Euler method. This method takes two steps to achieve 2nd-order precision.
In theory it should be comparable to the backward Euler method, but we'll see.
***********************************************************************************
From running tests, it appears that the backwards Euler method is significantly
faster than all the others. In fact, the lower the order of the method, the faster
it runs across the board. This is in agreement I guess with what would be
expected of a stiff problem like this one. The effect is even more pronounced at
stricter tolerances. In general, then, stick with METHOD = 1;

To Do:

-When tolerances are strict and the particle number is high, the code stalls at
	low temperature. I think the Jacobian goes singular there. This should be
	fixed by adding a small random bump to all the k-values and trying from there,
	as it should eliminate the singularity.
-Remove nonsense dependencies that I don't have anymore
*/

double PARTICLE_MASS = 100.0; //GeV
double COUPLING_CONST = 0.02;
double Ns = 4e-12; //number density per comoving volume
double mPlanck = 1.2209e19; //Planck Mass, units of GeV/c^2 (c = 1)
double currentEntropy = 2969.5; //entropy density in units of cm^-3
double currentEnergyDensity = 1.26477e-6; //current dark matter energy density in GeV/cm^3. This is *not* the local density
int MAX_SOLUTIONS = 100; //Largest number of independent solutions to be kept track of; Everything after this treated as 0
double THRESHHOLD = 1e-200; //Any particle that has a proportion of the total particle number smaller than this is ignored;
double epsilon = 1e-300; //Treat this as the smallest number that can exist. Is basically 0.
double TOL = 1e-6; //Tolerance
char fileName[40] = "Output";

const int NUM_THREADS = 4;

#include <stdio.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include "ADMmodel.h"
#ifdef PREDEF_PLATFORM_UNIX
	#include <boost/thread.hpp>
#endif

#ifndef PRODUCTIONRATE_H
	#include "productionRate.h"
	#define PRODUCTIONRATE_H
#endif

#ifndef MASSCALC_H
	#include "massCalc.h"
	#define MASSCALC_H
#endif

#ifndef SPLINE_H
	#define SPLINE_H
	#include <gsl/gsl_spline.h>
#endif

#ifndef ERRNO_H
	#define ERRNO_H
	#include <gsl/gsl_errno.h>
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;
typedef std::vector<double> dvector;
typedef std::vector<dvector> dmatrix;
typedef std::vector<int> intvector;
typedef std::vector<ADMmodel> modvector;
typedef std::pair<int, double> point;
typedef std::pair<double, int> pInverse;
typedef std::vector<point> locus;
typedef std::vector<pInverse> locInverse;
typedef std::vector<gsl_spline *> splineVector;
typedef std::vector<splineVector> splineMatrix;
typedef std::vector<gsl_interp_accel *> accelVector;
typedef std::vector<accelVector> accelMatrix;



enum TAG {DISCRETE, SPLINE};

typedef std::vector<TAG> instruction;

void lineBreak()
{
	printf("****************************************\n");
	cout << flush;
	return;
}

bool same_sign(double a, double b)
{
    return copysign(a,b) == a;
}

double gStar(double T)
{
	if(T > 173.5) //All particles
		return 106.75;
	else if(T > 126.4) //no top
		return 96.25;
	else if(T > 80.4) //no higgs
		return 95.25;
	else if(T > 4.18) //no W or Z
		return 86.25;
	else if(T > 1.776) //no bottom
		return 75.75;
	else if(T > 1.275) //no tau
		return 72.25;
	else if(T > 0.15) //no charm
		return 61.75;
	else if(T > 0.1) //QCD phase transition (all quarks in pions)
		return 17.25;
	else if(T > 5e-4) //no pions or mu
		return 10.75;
	else
		return 3.36; //photons and neutrinos only
}

double entropy(double T)
{
	return (2*M_PI*M_PI/45)*gStar(T)*T*T*T;
}

//Define a bunch of functions relating H, x, and T. These functions are all described on pg118 in The Early Universe.
double x_of_T(double T){return PARTICLE_MASS/T;}
double T_of_x(double x){return PARTICLE_MASS/x;}
double t_of_T(double T){return 0.301*mPlanck/(T*T*sqrt(gStar(T)));}
double t_of_x(double x){return t_of_T(T_of_x(x));}
double H_of_m(double T){return 1.67*PARTICLE_MASS*PARTICLE_MASS*sqrt(gStar(T))/mPlanck;}
double H_of_x(double x){return H_of_m(T_of_x(x))/(x*x);}
double H_of_T(double T){return H_of_x(x_of_T(T));}
double T_from_t(double t){return 1.33/sqrt(t);}
double z_of_x(double x){return log(x);}
double x_of_z(double z){return exp(z);}
double z_of_T(double T){return z_of_x(x_of_T(T));}
double T_of_z(double z){return T_of_x(x_of_z(z));}



//Just for easy generation of numbers spaced pretty evenly on a logarithmic scale
double vecGen(int in)
{
	double out;
	if(in < 0 && in%5 != 0)
		in -=5;

	switch(in%5)
	{
	case 0:
		out = 1.0;
		break;
	case 1:
	case -4:
		out = 1.5;
		break;
	case 2:
	case -3:
		out = 2.5;
		break;
	case 3:
	case -2:
		out = 4.0;
		break;
	case 4:
	case -1:
		out = 6.5;
		break;
	}
	return out*pow(10.0,in/5);
}

class admSolver
{
	public:
		admSolver(double, double, double, int);
		int solve();
		double xMax;
		double xMin;
		double zMax;
		double zMin;
		ADMmodel * Model;
		double evalRateSpline(int, int, double);
		void createSplineLine(int, int, double, double);
	protected:
		void analyzeSolution(const dvector &, double, double);
		double mx;
		double gp;
		double ns;
		dvector Solutions;
		int solverStep(dvector &, const dvector &, dmatrix &, double, double, double);
		dmatrix butcherTableau;
		int rungeStepNum;
		int rungeMethod;
		splineMatrix rateSplineMat;
		accelMatrix accelMat;
		void readyRateSplines(int);
		void addSplines(int, double, double);
		void createSpline(int, int, const double[], int);
		void setMethod(int);
};

admSolver::admSolver(double m, double g, double n, int meth)
{
	mx = m;
	gp = g;
	ns = n;
	setMethod(meth);
	xMin = 0.25;
	zMin = z_of_x(xMin);
	xMax = 2e9;
	zMax = z_of_x(xMax);
	Model = new ADMmodel(mx, gp, T_of_x(xMax));
	Solutions.clear();
	Solutions.assign(10, 0.0);
	Solutions[0] = 1.0;
	readyRateSplines(MAX_SOLUTIONS);
	xMin = 0.5;
	zMin = z_of_x(xMin);
	xMax = 1e11;
	zMax = z_of_x(xMax);
}

void admSolver::setMethod(int meth)
{
	rungeMethod = meth;
	dvector tempVec;
	dmatrix tempMat;
	switch(meth)
	{
		case 1:
			tempVec.assign(1, 1.0);
			tempMat.assign(3, tempVec);
			rungeStepNum = 1;
			break;
		case 2:
			tempVec.assign(1, 0.5);
			tempMat.assign(3, tempVec);
			tempMat[1][0] = 1.0;
			rungeStepNum = 1;
			break;
		case 3:
			tempVec.assign(2, 0.5);
			tempMat.assign(4, tempVec);
			tempMat[0][0] = (3.0 + sqrt(3.0))/6.0;
			tempMat[0][1] = (3.0 - sqrt(3.0))/6.0;
			tempMat[2].assign(1, tempMat[0][0]);
			tempMat[3][0] = -1.0/sqrt(3.0);
			tempMat[3][1] = tempMat[0][0];
			rungeStepNum = 2;
			break;
		case 4:
			tempVec.assign(3, 1.0);
			tempMat.assign(5, tempVec);
			tempMat[0][0] = 1.0;
			tempMat[0][1] = 0.5;
			tempMat[0][2] = 0.0;
			tempMat[1][0] = 1.0/6.0;
			tempMat[1][1] = 2.0/3.0;
			tempMat[1][2] = 1.0/6.0;
			tempMat[2].assign(1, 1.0);
			tempMat[3].assign(2, -0.75);
			tempMat[3][1] = 1.25;
			tempMat[4][0] = 2.0;
			tempMat[4][1] = -3.0;
			tempMat[4][2] = 1.0;
			rungeStepNum = 3;
			break;
		case 5:
			tempVec.assign(2, 0.5);
			tempMat.assign(4, tempVec);
			tempMat[0][0] = 2.0;
			tempMat[0][1] = -1.0;
			tempMat[1][0] = 0.4;
			tempMat[1][1] = 0.6;
			tempMat[2].assign(1, tempMat[0][0]);
			tempMat[3][0] = -2.0;
			tempMat[3][1] = 1.5;
			rungeStepNum = 2;
			break;
	}
	butcherTableau = dmatrix(rungeStepNum+2,dvector(rungeStepNum,1.0));
	for(int k = 0; k < rungeStepNum+2; k++)
	{
		butcherTableau[k]=tempMat[k];
	}
	return;
}

void admSolver::createSpline(int n, int m, const double zVals[], int numVals)
{
	double * rates = new double[numVals];
	for(int k = 0; k < numVals; k++)
	{
		rates[k] = Model->singleRateCalculation(n, m, T_of_x(x_of_z(zVals[k])));
	}
	gsl_spline * gSpl = gsl_spline_alloc(gsl_interp_cspline, numVals);
	gsl_interp_accel * acc = gsl_interp_accel_alloc();
	gsl_spline_init(gSpl, zVals, rates, numVals);
	rateSplineMat[n-1][m-1] = gSpl;
	rateSplineMat[m-1][n-1] = gSpl;
	accelMat[n-1][m-1] = acc;
	accelMat[m-1][n-1] = acc;
	delete [] rates;
	return;
}

void admSolver::createSplineLine(int n, int m, double zStart, double zEnd)
{
	int numVals = int(5.0*log(fabs(zEnd-zStart)));
	double * zVals = new double[numVals+5];

	double zT = zStart;
	for(int k = 0; k < numVals+5; k++)
	{
		if(k ==0)
		{
			zT = zStart;
		}
		else if(k == numVals+1)
		{
			zT = zEnd;
		}
		else if(k > numVals+1)
		{
			zT = zEnd + double(k-(numVals+1))*fabs(zEnd-zStart)*0.125;
		}
		else
		{
			zT = zStart + pow(fabs(zEnd - zStart), double(k)/double(numVals+1));
		}
		zVals[k] = zT;
	}
	if(m < 11)
	{
		for(int k2 = 1; k2 <= m; k2++)
		{
			createSpline(k2, n, zVals, numVals+5);
		}
		printf("Finished calculating rates for n = %d\n", n);
		delete [] zVals;
		return;
	}

	int numRates = int(1.75*pow(log(double(fabs(m-4))),1.25));
	double * rates = new double[numRates+4];
	double alpha = pow(double(n - 2 - (numRates-1))/3.0, 1.0/double(numRates-1));;
	int * rateCalcs = new int[numRates+4];
	double * rateCalcsD = new double[numRates+4];
	for(int k3 = 0; k3 < numRates+4; k3++)
	{
		if(k3 < 3)
		{
			rateCalcs[k3] = k3+1;
		}
		else if (k3 > numRates)
		{
			rateCalcs[k3] = n - (numRates+3-k3);
		}
		else
		{
			rateCalcs[k3] = (k3-2) + int(3.0*pow(alpha, double(k3-2))+0.5);
		}
		rateCalcsD[k3] = double(rateCalcs[k3]);
	}
	gsl_spline * gTemp;
	gsl_interp_accel * accTemp;
	splineVector vecTemp;
	vecTemp.assign(numVals+5, gTemp);
	accelVector accVecTemp;
	accVecTemp.assign(numVals+5, accTemp);

	for(int k4 = 0; k4 < numVals+5; k4++)
	{
		for(int k5 = 0; k5 < numRates+4; k5++)
		{
			rates[k5] = Model->singleRateCalculation(rateCalcs[k5], n, T_of_x(x_of_z(zVals[k4])));
		}
		vecTemp[k4] = gsl_spline_alloc(gsl_interp_cspline, numRates+4);
		accVecTemp[k4] = gsl_interp_accel_alloc();
		gsl_spline_init(vecTemp[k4], rateCalcsD, rates, numRates+4);
	}

	delete [] rates;
	double * rates2 = new double[numVals+5];
	for(int k6 = 1; k6 <= m; k6++)
	{
		for(int k7 = 0; k7 < numVals+5; k7++)
		{
			rates2[k7] = gsl_spline_eval(vecTemp[k7], double(k6), accVecTemp[k7]);
		}
		rateSplineMat[n-1][k6-1] = gsl_spline_alloc(gsl_interp_cspline, numVals+5);
		accelMat[n-1][k6-1] = gsl_interp_accel_alloc();
		gsl_spline_init(rateSplineMat[n-1][k6-1], zVals, rates2, numVals+5);
		rateSplineMat[k6-1][n-1] = rateSplineMat[n-1][k6-1];
		accelMat[k6-1][n-1] = accelMat[n-1][k6-1];
	}

	for(int k8 = 0; k8 < numVals+5; k8++)
	{
		gsl_spline_free(vecTemp[k8]);
		gsl_interp_accel_free(accVecTemp[k8]);
	}
	delete [] zVals;
	delete [] rates2;
	printf("Finished calculating rates for n = %d (Interpolated from %d values)\n", n, numRates+4);
	return;
}

void parallelCreateSplineLine(admSolver * solver, int n, int m, double zStart, double zEnd)
{
	for(int k = n; k <= m; k++)
	{
		solver->createSplineLine(k, k, zStart, zEnd);
	}
	return;
}

#ifdef PREDEF_PLATFORM_UNIX
void admSolver::readyRateSplines(int n)
{
	gsl_spline * gTemp;
	gsl_interp_accel * accTemp;
	splineVector vecTemp;
	vecTemp.assign(n, gTemp);
	accelVector accVecTemp;
	accVecTemp.assign(n, accTemp);
	rateSplineMat.assign(n, vecTemp);
	accelMat.assign(n, accVecTemp);
	std::vector<boost::thread *> threadVec;
	if(n < 50 || true)
	{
		int numPerThread = n / NUM_THREADS;
		int numExtra = n - (numPerThread*NUM_THREADS);
		for(int k1 = 0; k1 < NUM_THREADS; k1++)
		{
			if(k1 == 0)
			{
				threadVec.push_back(new boost::thread(parallelCreateSplineLine, this, 1, numPerThread+numExtra, zMin, zMax));
			}
			else
			{
				threadVec.push_back(new boost::thread(parallelCreateSplineLine, this, (k1*numPerThread)+numExtra+1, ((k1+1)*numPerThread)+numExtra, zMin, zMax));
			}
		}
		for(int k2 = 0; k2 < NUM_THREADS; k2++)
		{
			threadVec[k2]->join();
		}
	}
	return;
}
#else
void admSolver::readyRateSplines(int n)
{
	gsl_spline * gTemp;
	gsl_interp_accel * accTemp;
	splineVector vecTemp;
	vecTemp.assign(n, gTemp);
	accelVector accVecTemp;
	accVecTemp.assign(n, accTemp);
	rateSplineMat.assign(n, vecTemp);
	accelMat.assign(n, accVecTemp);
	if(n < 50 || true)
	{
		for(int k = 0; k < n; k++)
		{
			createSplineLine(k, k, zMin, zMax);
		}
	}
	return;
}
#endif

void admSolver::addSplines(int n, double zStart, double zEnd)
{
	int nSiz = rateSplineMat.size();
	for(int k1 = 0; k1 < nSiz; k1++)
	{
		rateSplineMat[k1].resize(nSiz+n);
		accelMat[k1].resize(nSiz+n);
	}
	for(int k2 = 0; k2 < n; k2++)
	{
		rateSplineMat.push_back(rateSplineMat[0]);
		accelMat.push_back(accelMat[0]);
	}
	for(int k3 = nSiz+1; k3 <= nSiz+n; k3++)
	{
		createSplineLine(k3, k3, zStart, zEnd);
	}

	return;
}

double admSolver::evalRateSpline(int n1, int n2, double zT)
{
	return gsl_spline_eval(rateSplineMat[n1-1][n2-1], zT, accelMat[n1-1] [n2-1]);
}


void admSolver::analyzeSolution(const dvector & inDat, double z, double step)
{
	int numZero = 0;
	int numNeg = 0;
	int maxInd = 0;
	double maxVal = 0.0;
	double particleTotal = 0.0;
	double temp = 0.0;
	double energyDensity = 0.0;
	int siz = inDat.size();
	printf("\n");
	for(int k = 0; k < siz; k++)
	{
		temp = inDat[k]*Ns;
		printf("For n = %d, the density is %.5e\n", k+1, temp);
		if(temp <= 10*epsilon && inDat[k] >= 0.0)
			numZero++;

		if(inDat[k] < 0)
			numNeg++;

		if(temp > maxVal)
		{
			maxVal = temp;
			maxInd = k+1;
		}
		
		particleTotal += temp*(double)(k+1);
		energyDensity += temp * currentEntropy * Model->getMass(double(k+1));
	}
	printf("\nEvaluated at x = %.3e (step of %.3e):\n\n", x_of_z(z), step);
	printf("N up to size %d considered.\n", int(inDat.size()));
	printf("We encountered %d zeros and %d negative results.\n", numZero, numNeg);
	printf("The maximum value was %.5e, for n = %d.\n", maxVal, maxInd);
	printf("The current energy density is %.5e GeV/cm^3\n", energyDensity);
	printf("The total number of particles in this universe is %.5e.\n", particleTotal);

	ofstream fout;
	fout.open(fileName, ios::app);
	printf("Saving to file name %s\n", fileName);

	fout << x_of_z(z) << ",";
	for(int j = 0; j < siz; j++)
	{
		fout << inDat[j]*Ns;
		if(j < siz)
			fout << ",";

	}
	fout << endl;
	fout.close();

	cout << endl;
	return;
}

void vecAdd(dvector & mainVec, const dvector & addVec)
{
	int siz = mainVec.size();
	if(siz > addVec.size())
	{
		siz = addVec.size();
	}
	for(int k = 0; k < siz; k++)
	{
		mainVec[k] += addVec[k];
	}
	return;
}

void vecAdd(dvector & mainVec, const dvector & addVec, double alpha)
{
	int siz = mainVec.size();
	if(siz > addVec.size())
	{
		siz = addVec.size();
	}
	for(int k = 0; k < siz; k++)
	{
		mainVec[k] += alpha*addVec[k];
	}
	return;
}

void vecSubt(dvector & mainVec, const dvector & subtVec)
{
	int siz = mainVec.size();
	if(siz > subtVec.size())
	{
		siz = subtVec.size();
	}
	for(int k = 0; k < siz; k++)
	{
		mainVec[k] -= subtVec[k];
	}
	return;
}

void vecVecMult(dvector & mainVec, const dvector & multVec)
{
	int siz = mainVec.size();
	if(siz > multVec.size())
	{
		siz = multVec.size();
	}
	for(int k = 0; k < siz; k++)
	{
		mainVec[k] *= multVec[k];
	}
	return;
}

void vecConstMult(dvector & mainVec, double multCoef)
{
	for(int k = 0; k < mainVec.size(); k++)
	{
		mainVec[k] *= multCoef;
	}
	return;
}

void printVector(const dvector & mainVec)
{
	for(int k = 0; k < mainVec.size(); k++)
	{
		printf("k = %d -- %.5e\n", k, mainVec[k]);
	}
	printf("\n");
	return;
}

void printVector(const dvector & mainVec, const dvector & nextVec)
{
	for(int k = 0; k < mainVec.size(); k++)
	{
		printf("k = %d -- %.12e -- %.12e\n", k, mainVec[k], nextVec[k]);
	}
	printf("\n");
	return;
}

void printVector(const dvector & mainVec, const dvector & nextVec, const dvector & thirdVec)
{
	for(int k = 0; k < mainVec.size(); k++)
	{
		printf("k = %d : %.12e -- %.12e -- %.12e\n", k+1, mainVec[k], nextVec[k], thirdVec[k]);
	}
	printf("\n");
	return;
}

double estimateAbsError(const dvector & sol1, const dvector & sol2)
{
	int siz = sol1.size();
	if(sol2.size() < siz)
	{
		siz = sol2.size();
	}
	double error = 0.0;
	double logSum = 1.0;
	double temp;
	for(int k = 0; k < siz; k++)
	{
		temp = fabs((sol1[k] - sol2[k]));//log(10.0*(sol1[k] + sol2[k])/THRESHHOLD));
		error += temp*temp;
		//printf("k = %d, Result = %.5e, Err est = %.5e: Difference of %.5e\n", k, sol1[k], sol2[k], temp);
		//logSum += pow(log(10.0*(sol1[k] + sol2[k])/THRESHHOLD),2);
	}
	return sqrt(error)/sqrt((double)siz*logSum);
}

double estimateAbsError(const dvector & sol1)
{
	int siz = sol1.size();
	double error = 0.0;
	for(int k = 0; k < siz; k++)
	{
		error += sol1[k]*sol1[k];
		//printf("k = %d, error = %.5e\n", k, fabs(sol1[k]));
	}
	return sqrt(error)/sqrt((double)siz);
}

double estimateRelError(const dvector & sol1, const dvector & sol2)
{
	int siz = sol1.size();
	if(sol2.size() < siz)
	{
		siz = sol2.size();
	}
	double error = 0.0;
	double temp1 = 0.0;
	double temp2 = 0.0;
	for(int k = 0; k < siz; k++)
	{
		temp1 = sol1[k];
		temp2 = sol2[k];
		if(temp1 <= epsilon || temp2 <= epsilon)
		{
			continue;
		}
		if(temp1 > temp2)
		{
			error += (temp1-temp2)*(temp1-temp2)/(temp1*temp1);
		}
		else
		{
			error += (temp1-temp2)*(temp1-temp2)/(temp2*temp2);
		}
	}
	return sqrt(error);
}

double estimateRelErrorInd(const dvector & solErr, const dvector & solDat)
{
	int siz = solErr.size();
	if(siz > solDat.size())
	{
		siz = solDat.size();
	}
	double error = 0.0;
	for(int k = 0; k < siz; k++)
	{
		if(solDat[k] <= epsilon)
		{
			continue;
		}
		error += (solErr[k]*solErr[k])/(solDat[k]*solDat[k]);
	}
	return sqrt(error);
}

bool noNegVal(const dvector & inY, const dvector & inDelta)
{
	for(int k = 0; k < inY.size(); k++)
	{
		if(inY[k] + inDelta[k] < 0.0)
		{
			printf("\n\nk = %d, Y[k] = %.5e, D[k] = %.5e\n\n", k, inY[k], inDelta[k]);
			//string bleh;
			//cin >> bleh;
			return false;
		}
	}
	return true;
}

int populateKSolver(double * F, double * Jac, const dvector & yOld, const dvector & kNew, double baseStep, double aCoef, double alphaStep, double x, const dmatrix & rateMat, const dmatrix & decayMat)
{
	int n = yOld.size();
	double alpha = (x*x*entropy(T_of_x(x)))/(H_of_m(T_of_x(x)));
	double tempj1, tempj2, tempj3, tempf;
	double tempk1, tempk2, tempy1, tempy2;
	double temp;
	double alphaStepTemp;
	double step = baseStep*aCoef;
	for(int k1 = 0; k1 < n; k1++)
	{
		tempy1 = yOld[k1];
		tempk1 = kNew[k1];

		for(int k2 = 0; k2 + 0*k1 < n; k2++)
		{
			tempk2 = kNew[k2];
			tempy2 = yOld[k2];
			if(k1 == k2)
			{
				alphaStepTemp = 2.*alphaStep;
			}
			else
			{
				alphaStepTemp = 2.*alphaStep - 1;
			}
			if(alphaStepTemp < 0)
			{
				alphaStepTemp = 0;
			}
			if(alphaStepTemp > 1.)
			{
				alphaStepTemp = 1.;
			}

			//Collisional terms
			tempf = alpha*rateMat[k1][k2]*Ns;
			tempj1 = tempf * step * (tempy2 + alphaStepTemp*step*tempk2);
			tempj2 = tempf * step * (tempy1 + alphaStepTemp*step*tempk1);
			tempf *= tempy1*tempy2 + step*(tempy1*tempk2+tempy2*tempk1) + step*step*alphaStepTemp*tempk1*tempk2;

			//Decay terms
			if(k1+k2 < n-1)
			{
				tempj3 = alpha*decayMat[k1][k2]*step;
				tempf -= alpha*decayMat[k1][k2]*(yOld[k1+k2+1]+step*kNew[k1+k2+1]);
			}

			Jac[k1*n + k1] += tempj1;
			Jac[k1*n + k2] += tempj2;
			Jac[k2*n + k1] += tempj1;
			Jac[k2*n + k2] += tempj2;
			F[k1] -= tempf;
			F[k2] -= tempf;

			if(k1+k2 < n-1)
			{
				Jac[(k1+k2+1)*n + k1] -= tempj1;
				Jac[(k1+k2+1)*n + k2] -= tempj2;
				Jac[k1*n + (k1+k2+1)] -= tempj3;
				Jac[k2*n + (k1+k2+1)] -= tempj3;
				Jac[(k1+k2+1)*(n+1)] += tempj3;
				F[k1+k2+1] += tempf;
			}

		}
		Jac[k1*n + k1] += 1.;
		F[k1] -= tempk1;
	}
	return 0;
}

int stepKSolver(double * F, double * F2, double * Jac, double * Jac2, const dvector & yOld, dvector & kNew, double baseStep, double aCoef, double alphaStep, double x, const dmatrix & rateMat, const dmatrix & decayMat, double tol)
{
	int n = yOld.size();
	double temp;
	double temp2;
	double tempSum = 1.;
	dvector kTemp;
	dvector kTemp2;
	kTemp.assign(n, 0.0);
	kTemp2.assign(n, 0.0);
	gsl_matrix_view m = gsl_matrix_view_array(Jac, n, n);
	gsl_matrix_view m2 = gsl_matrix_view_array(Jac2, n, n);
	gsl_vector_view v = gsl_vector_view_array(F, n);
	gsl_vector_view v2 = gsl_vector_view_array(F2, n);
	gsl_vector * outB = gsl_vector_alloc(n);
	int s;
	gsl_permutation * p = gsl_permutation_alloc(n);
	int count = 0;
	int status = 0;
	while(abs(tempSum) > tol)
	{
		tempSum = 0.;
		for(int k1 = 0; k1 < n; k1++)
		{
			F[k1] = 0.;
			F2[k1] = 0.;
			kTemp[k1] = kNew[k1];
			kTemp2[k1] = 0.0;
		}
		for(int k2 = 0; k2 < n*n; k2++)
		{
			Jac[k2] = 0.;
			Jac2[k2] = 0.;
		}
		populateKSolver(F, Jac, yOld, kNew, baseStep, aCoef, alphaStep, x, rateMat, decayMat);
		//gsl_matrix_fprintf(stdout, &m.matrix, "%g");
		//status = gsl_linalg_HH_solve(&m.matrix, &v.vector, outB);
		gsl_linalg_LU_decomp(&m.matrix, p, &s);
		status = gsl_linalg_LU_solve(&m.matrix, p, &v.vector, outB);
		if(status!=0 && status != 21)
		{
			printf("Error: %s\n", gsl_strerror(status));
			printf("Error code = %d\n", status);
			abort();
		}
		else if(status == 21)
		{
			printf("Error: %s\n", gsl_strerror(status));
			return -1;
		}
		for(int k3 = 0; k3 < n; k3++)
		{
			temp = gsl_vector_get(outB, k3);
			kTemp[k3] = kNew[k3]+temp;
		}
		populateKSolver(F2, Jac2, yOld, kTemp, baseStep, aCoef, alphaStep, x, rateMat, decayMat);
		
		/*
		printf("Real Jac2:\n");
		gsl_matrix_fprintf(stdout, &m2.matrix, "%g");
		printf("F:\n");
		gsl_vector_fprintf(stdout, &v.vector, "%g");*/
		gsl_linalg_LU_decomp(&m2.matrix, p, &s);
		
		status = gsl_linalg_LU_solve(&m2.matrix, p, &v2.vector, outB);
		if(status!=0 && status != 21)
		{
			printf("Error: %s\n", gsl_strerror(status));
			printf("Error code = %d\n", status);
			printf("Failed on 2nd step\n");
			abort();
		}
		else if(status == 21)
		{
			return -1;
		}
		else
		{
			for(int k5 = 0; k5 < n; k5++)
			{
				temp = kTemp[k5]-kNew[k5] + gsl_vector_get(outB, k5);
				temp *= 0.5;
				kTemp2[k5] = kNew[k5];
				kNew[k5] += temp;
				temp = kNew[k5]-kTemp2[k5];
				kTemp2[k5] = temp;
				tempSum += temp*temp;
			}
		}
		tempSum = sqrt(tempSum);
		count++;
		if (count%2000 == 1000)
		{
			printf("%d iterations so far, tempSum = %e\n", count, tempSum);
			cout << flush;
		}
		
		if(count >=2500 && count <= 2505)
		{
			printf("count = %d\n", count);
			printVector(kNew, kTemp, kTemp2);
		}
		
		if (count >= 2505)
		{
			abort();
		}
	}
	gsl_vector_free(outB);
	gsl_permutation_free(p);
	return count;
}

int kSolver(dvector & kNew, const dvector & yOld, double baseStep, double aCoef, admSolver * solver, double x, double tol)
{
	int n = yOld.size();
	kNew.assign(n, 0.0);
	dvector kTemp = kNew;
	double * F;
	double * F2;
	double * Jac;
	double * Jac2;
	F = new double [n];
	F2 = new double [n];
	Jac = new double [n*n];
	Jac2 = new double [n*n];
	for(int p1 = 0; p1 < n; p1++)
	{
		F[p1] = 0.0;
	}
	for(int p2 = 0; p2 < n*n; p2++)
	{
		Jac[p2] = 0.0;
	}
	double h = 0.0;
	double hStep = 0.0;
	double tempSum = tol*10.0;
	int p, p2;
	double temp;
	double temp2;
	double tempTotal;
	int count = 0;
	dmatrix rateMat, decayMat;
	dvector tempRateVec, tempDecayVec;
	tempRateVec.assign(n, 0.0);
	rateMat.assign(n, tempRateVec);
	decayMat.assign(n, tempRateVec);
	int ind = 0;
	int numSteps = 0;
	int status;
	for(int q1 = 0; q1 < n; q1++)
	{
		for(int q2 = 0; q2 <= q1; q2++)
		{
			rateMat[q1][q2] = solver->evalRateSpline(q1+1, q2+1, z_of_x(x));
			rateMat[q2][q1] = rateMat[q1][q2];

			decayMat[q1][q2] = solver->Model->singleDecayCalculation(q1+1, q2+1, T_of_x(x));
			decayMat[q2][q1] = decayMat[q1][q2];

			ind++;
		}
	}
	temp = 0.0;
	bool bDisplay = false;
	bool bCont = true;
	do
	{
		kTemp = kNew;
		count++;
		if(count == 1)
		{
			status = stepKSolver(F, F2, Jac, Jac2, yOld, kNew, baseStep, aCoef, h, x, rateMat, decayMat, tol);
			hStep = 1e-6*tol;
			h = hStep;
			count = 1;
			continue;
		}
		h += hStep;
		if(h > 1.)
		{
			h = 1.;
			bCont = false;
		}
		status = stepKSolver(F, F2, Jac, Jac2, yOld, kNew, baseStep, aCoef, h, x, rateMat, decayMat, tol);
		//printf("h = %.4e\n", h);
		//printVector(kTemp, kNew);
		if(status < 0)
		{
			kNew = kTemp;
			h -= 0.5*hStep;
			hStep *= 0.5;
			printf("/n/nERROR HERE/n/n");
			cout << flush;
			abort();
			continue;
		}
		numSteps += status;
		tempSum = 0.0;
		for(int k1 = 0; k1 < n; k1++)
		{
			temp = kNew[k1]-kTemp[k1];
			tempSum += temp*temp;
		}
		if(status < 100)
		{
			hStep *= 2.;
		}
		else
		{
			hStep *= 200.0/((double)status);
			printf("Took %d iterations, new step size is %e\n", status, hStep);
			cout << flush;
		}
	}while(h < 1. && bCont);
	delete[] F;
	delete[] F2;
	delete[] Jac;
	delete[] Jac2;
	return numSteps;
}

int admSolver::solverStep(dvector & kResult, const dvector & yOld, dmatrix & kSpace, double baseStep, double z, double tol)
{
	/*
		For a p-step implicit runge kutta method, kSpace should be a dmatrix of size Nx(p+1).
		kSpace[0] is used to store things of the form (y+hk), while the other entries hold
		the individual k's.

		To re-enumerate it here, the Runge Kutta butcher tableau should be a
		(p+2)x(p)-dimensional dmatrix. The 0-index vector stores values of c,
		the points in time at which to evaluate F. The 1-index vector stores values
		of b, the final coefficients of all the k's. And the remaining p vectors
		store the a_ij in order.
	*/
	int status = 0;
	int countstep = 0;
	int n = yOld.size();

	for(int k1 = 0; k1 < rungeStepNum; k1++)
	{
		kSpace[0] = yOld;

		for(int k2 = 0; k2 < k1; k2++)
		{
			vecAdd(kSpace[0], kSpace[k2+1], baseStep*butcherTableau[k1+2][k2]);
		}
		status = kSolver(kSpace[k1+1], kSpace[0], baseStep, butcherTableau[k1+2][k1], this, x_of_z(z + butcherTableau[0][k1]*baseStep), tol);
		countstep += status;
		if(k1 == 0)
		{
			kResult = kSpace[k1+1];
			vecConstMult(kResult, butcherTableau[1][k1]);
		}
		else
		{
			vecAdd(kResult, kSpace[k1+1], butcherTableau[1][k1]);
		}
	}
	return countstep;
}

int admSolver::solve()
{
	printf("\n");
	printf("Solving from x = %.2e to %.2e:\n\n", xMin, xMax);
	int highestN = 1;
	double zTemp = zMin;
	const double zFin = zMax;
	double zLastRateCalc = zTemp;
	double baseStep = 1e-35;
	double baseStepMax = 0.2;
	double baseStepMin = 1e-100;
	double Tol = TOL;
	double error;
	int countstep = 0;
	int totalCountStep = 0;
	int k = 0;
	int nSiz = 0;
	string bleh;
	bool bSizeIncrease;
	int status;
	dvector kTemp;
	int tempMeth1,tempMeth2;
	if(rungeMethod == 1)
	{
		tempMeth1 = 5;
		tempMeth2 = 1;
	}
	kTemp.assign(MAX_SOLUTIONS*2,0.0);
	dmatrix kSpace;
	kSpace.assign(rungeStepNum+4, kTemp);

	dvector kResult;
	dvector kMidSol;
	clock_t clockS = clock();

	printf("k = 0\n\n");
	clock_t loopClockS;
	clock_t loopClockE;
	while(zTemp < zFin)
	{
		loopClockS = clock();
		nSiz = Solutions.size();
		status = 0;
		countstep = 0;
		//First step of size h

		status = solverStep(kTemp, Solutions, kSpace, baseStep, zTemp, Tol);
		countstep += status;
		//Second step, again of size h
		kMidSol = Solutions;
		vecAdd(kMidSol, kTemp, baseStep);

		status = solverStep(kResult, kMidSol, kSpace, baseStep, zTemp+baseStep, Tol);
		countstep += status;
		vecAdd(kResult, kTemp);

		//The step of size 2*h, to calculate the error
		status = solverStep(kTemp, Solutions, kSpace, 2.0*baseStep, zTemp, Tol);
		countstep += status;

		//Now calculate the error

		vecConstMult(kTemp, 2.0);

		error = baseStep*estimateAbsError(kResult, kTemp);
		//printVector(kResult, kTemp);
		if(error > 2.0*Tol)
		{
			printf("Error is %.5e - restarting with stepsize %.3e (old was %.3e)\n", error, baseStep*0.6, baseStep);
			baseStep *= 0.6;
			lineBreak();
			//cin >> bleh;
			continue;
		}
		vecAdd(Solutions, kResult, baseStep);
		k++;
		zTemp += baseStep;
		totalCountStep += countstep;
		bSizeIncrease = false;
		analyzeSolution(Solutions, zTemp, baseStep);
		if(fabs(Solutions[nSiz-1]) > THRESHHOLD && Solutions.size()<(MAX_SOLUTIONS-10))
		{
			for(int prw = 0; prw < 5; prw++)
			{
				Solutions.push_back(0.0);
			}
			if(rateSplineMat.size() < nSiz)
			{
				addSplines(5, zTemp, zMax);
			}
			if(nSiz > highestN)
			{
				highestN = Solutions.size();
			}
			bSizeIncrease = true;
		}
		nSiz = Solutions.size();
		if(abs(Solutions[nSiz-1]) <= 1e-250 && nSiz > 8 && bSizeIncrease == false && false)
		{
			if(abs(Solutions[nSiz-2])<=1e-250 && abs(Solutions[nSiz-3])<=1e-250 && abs(Solutions[nSiz-4])<=1e-250)
			{
				Solutions.pop_back();
				Solutions.pop_back();
				Solutions.pop_back();
				Solutions.pop_back();
			}
		}

		printf("Absolute Error = %.5e\n", error);
		if(error > (Tol * 0.05))
		{
			baseStep *= 0.95*pow(Tol / error, 1.0/(((double)rungeStepNum)+1.0));
		}
		else
		{
			baseStep *= 2.0;
		}
		if(baseStep > baseStepMax)
		{
			baseStep = baseStepMax;
		}
		if(baseStep < baseStepMin)
		{
			baseStep = baseStepMin;
		}
		if(bSizeIncrease == true && nSiz == highestN)
		{
			baseStep *= 0.01;
		}
		loopClockE = clock();
		printf("Step completed in %d iterations and %.3e seconds.\n", countstep, ((double)(loopClockE - loopClockS)/CLOCKS_PER_SEC));
		cout << flush;
		lineBreak();
	}

	clock_t clockE = clock();
	printf("\nEnded when xTemp = %.3e. Time was %.3e\n", zTemp, ((double)clockE-clockS)/CLOCKS_PER_SEC);
	printf("Took %d iterations of the solver and %d matrix inverse evaluations.\n\n", k, totalCountStep);
	return 0;
}

int main()
{
	gsl_error_handler_t * errHand = gsl_set_error_handler_off();
	unsigned int fileInd = 1;
	char buf[100];
	char outName[10];
	if(METHOD == 4)
	{
		strcpy(outName, "Out_TeM4");
	}
	else if(METHOD == 3)
	{
		strcpy(outName, "Out_TeM3");
	}
	else if(METHOD == 2)
	{
		strcpy(outName, "Out_Mid2");
	}
	else if(METHOD == 1)
	{
		strcpy(outName, "Out_BEul");
	}
	else
	{
		strcpy(outName, "Out_data");
	}
	std::fstream fs;
	while(true)
	{
		sprintf(buf, "%.8s_M%.1e_G%.1e_N%.1e_R%.3u.csv", outName, PARTICLE_MASS, COUPLING_CONST, Ns, fileInd);
		fs.open(buf, std::fstream::in);
		if(fs.good() == false)
		{
			fs.close();
			strcpy(fileName, buf);
			break;
		}
		else
		{
			fs.close();
			fileInd++;
		}
	}

	admSolver solver(PARTICLE_MASS, COUPLING_CONST, Ns, METHOD);

	solver.solve();

	return 0;
}
#endif
