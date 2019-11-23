#ifndef ADM_SOLVER
#define ADM_SOLVER

int METHOD = 1;//1 = Heun, 2 = Bogacki-Sampine, 3 = Backward Euler

double PARTICLE_MASS = 140.0; //GeV
double COUPLING_CONST = 0.5;
double Ns = 1e-4; //number density per comoving volume
double mPlanck = 1.2209e19; //Planck Mass, units of GeV/c^2 (c = 1)
double Trec = 3.08e-9; //Recombination temperature, everything should be frozen out by now. In units of GeV.
int MAX_SOLUTIONS = 100; //Largest number of independent solutions to be kept track of; Everything after this treated as 0
const int NUM_SOLUTIONS = MAX_SOLUTIONS;
double MAX_NUMBER = 1000.0; //Start the production rate interpolation with this as the maximum value of N
double THRESHHOLD = 1e-200; //Any particle that has a proportion of the total particle number smaller than this is ignored;
double THRESHHOLD_TWO = 1e-40;
double epsilon = 1e-200; //Treat this as the smallest number that can exist. Is basically 0.
char fileName[40] = "Output";

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <time.h>
#include "ADMmodel.h"
#include <assert.h>
#include <typeinfo>
#include <memory>

#ifndef PRODUCTIONRATE_H
	#include "productionRate.h"
	#define PRODUCTIONRATE_H
#endif

#ifndef MASSCALC_H
	#include "massCalc.h"
	#define MASSCALC_H
#endif

#ifndef TWO_DIM_INTERP_H
	#include "twoDimInterp.h"
	#define TWO_DIM_INTERP_H
#endif

#ifndef THREE_DIM_INTERP_APPROX_H
	#include "threeDimInterpApprox.h"
	#define THREE_DIM_INTERP_APPROX_H
#endif

#ifndef SPLINE_H
	#define SPLINE_H
	#include <gsl/gsl_spline.h>
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;
typedef std::vector<double> dvector;
typedef std::vector<int> intvector;
typedef std::vector<ADMmodel> modvector;
typedef std::pair<int, double> point;
typedef std::pair<double, int> pInverse;
typedef std::vector<point> locus;
typedef std::vector<pInverse> locInverse;

enum TAG {DISCRETE, SPLINE};

typedef std::vector<TAG> instruction;

void lineBreak()
{
	printf("****************************************\n");
	cout << flush;
	return;
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
class admSolver;
class admSolution
{
	public:
		admSolution(double, admSolver * );
		admSolution(admSolution *, double, double);
		admSolution(admSolution *);
		admSolution();
		int get_maxN();
		double get_xData();
		double getInterp(double, double);
		void solve(admSolution *);
		void solve();
		double evaluate(int);
		double evaluateDelta(int);
		double nFunc(int);
		double * outData;
		int dataSize;
		void analyzeSolution();
		double z;
		double x;
		double h;
		void set_z(double);
		
		void addSolution(double, admSolution *, double, bool);
		void addSolution(double, admSolution *, double);
		void addSolution(admSolution *, double, bool);
		void addSolution(admSolution *, double);
		void addSolution(admSolution *);
	private:
		admSolution * prev;
		bool bEvaluationReady;
		double calcSingleGain(int, int);
		admSolver * solver;
		int MAX_N;
};

class admSolver
{
	public:
		admSolver(double, double, double, double);
		void addADM(double);
		void initialize();
		void increaseMaxN();
		double getApproxInterp(double, double, double);
		void solveForward();
		dvector xvals;
		double xMax;
		double xMin;
		double zMax;
		double zMin;
		double MAX_N;
	protected:
		modvector Models;
		threeDimInterpApprox approxInterp;
		dvector nums;
		double mx;
		double gp;
		double ns;
		std::vector<admSolution> Solutions;
		void setApproxInterp();
		bool bStep;
};

admSolution::admSolution()
{
	z = -1.0;
	x = -1.0;
	h = -1.0;
	bEvaluationReady = false;
	dataSize = 0;
	MAX_N = 0;
}

admSolution::admSolution(double zIn, admSolver * solveIn)
{
	z = zIn;
	x = x_of_z(z);
	solver = solveIn;
	MAX_N = (int)(solver->MAX_N+0.1);
	h = 0.0;
	bEvaluationReady = true;
	outData = new double[1];
	outData[0] = 1;
	dataSize = 1;
}

admSolution::admSolution(admSolution * prevIn, double zIn, double hIn)
{
	z = zIn;
	x = x_of_z(z);
	h = hIn;
	prev = prevIn;
	solver = prev->solver;
	MAX_N = (int)(solver->MAX_N+0.1);
	bEvaluationReady = false;
	dataSize = 0;
	solve();
}

admSolution::admSolution(admSolution * prevIn)
{
	prev = prevIn->prev;
	x = prevIn->x;
	h = prevIn->h;
	solver = prev->solver;
	MAX_N = (int)(solver->MAX_N+0.1);
	outData = new double[prevIn->get_maxN()];
	for(int k = 0; k < prevIn->get_maxN(); k++)
	{
		outData[k] = prevIn->evaluate(k+1);
	}
	bEvaluationReady = true;
	dataSize = prevIn->get_maxN();
}

void admSolution::addSolution(double aCoef, admSolution * stepIn, double bCoef, bool bUseDelta)
{
	if(stepIn->get_maxN() > get_maxN())
	{
		double * outDataTemp = new double[stepIn->get_maxN()];
		for(int j = 0; j < stepIn->get_maxN(); j++)
		{
			outDataTemp[j] = evaluate(j+1);
		}
		outData = outDataTemp;
		dataSize = stepIn->get_maxN();
	}
	if(bUseDelta == true)
	{
		for(int k = 0; k < get_maxN(); k++)
		{
			outData[k] *= aCoef;
			outData[k] += bCoef * stepIn->evaluateDelta(k+1);
		}
	}
	else
	{
		for(int k = 0; k < get_maxN(); k++)
		{
			outData[k] *= aCoef;
			outData[k] += bCoef * stepIn->evaluate(k+1);
		}
	}
	for(int l = 0; l < get_maxN(); l++)
	{
		if(outData[l] < epsilon)
		{
			outData[l] = 0.0;
		}
	}
	return;
}

void admSolution::addSolution(double aCoef, admSolution * stepIn, double bCoef)
{
	addSolution(aCoef, stepIn, bCoef, false);
	return;
}

void admSolution::addSolution(admSolution * stepIn, double bCoef, bool bUseDelta)
{
	addSolution(1.0, stepIn, bCoef, bUseDelta);
	return;
}

void admSolution::addSolution(admSolution * stepIn, double bCoef)
{
	addSolution(1.0, stepIn, bCoef, false);
	return;
}

void admSolution::addSolution(admSolution * stepIn)
{
	addSolution(1.0, stepIn, 1.0, false);
	return;
}

void admSolution::solve(admSolution * prevIn)
{
	prev = prevIn;
	h = x - prev->x;
	solve();
	return;
}

void admSolution::analyzeSolution()
{
	int numZero = 0;
	int numNeg = 0;
	int maxInd = 0;
	double maxVal = 0.0;
	double particleTotal = 0.0;
	printf("\n");
	for(int k = 0; k < this->get_maxN(); k++)
	{
		printf("For n = %d, the density is %.5e\n", k+1, outData[k]*Ns);
		if(outData[k]*Ns <= 10*epsilon && outData[k] >= 0.0)
			numZero++;
		
		if(outData[k] < 0)
			numNeg++;
		
		if(outData[k]*Ns > maxVal)
		{
			maxVal = outData[k]*Ns;
			maxInd = k+1;
		}
		
		particleTotal += outData[k]*(double)(k+1);
	}
	printf("\nEvaluated at z = %.3e (step of %.3e):\n\n", z, h);
	printf("N up to size %d considered.\n", get_maxN());
	printf("We encountered %d zeros and %d negative results.\n", numZero, numNeg);
	printf("The maximum value was %.5e, for n = %d.\n", maxVal, maxInd);
	printf("The total number of particles in this universe is %.5e.\n", particleTotal*Ns);
	if(h > 1e-14)
	{
		ofstream fout;
		fout.open(fileName, ios::app);
		printf("Saving to file name %s\n\n", fileName);
	
		fout << x << ",";
		for(int j = 1; j <= get_maxN(); j++)
		{
			fout << evaluate(j)*Ns;
			if(j < get_maxN())
				fout << ",";
			
		}
		fout << endl;
		fout.close();
	}
	else
	{
		cout << endl;
	}
	return;
}

void admSolution::solve()
{
	dataSize = prev->get_maxN();
	dvector workingData;
	workingData.reserve(dataSize);
	for(int qot = 0; qot < dataSize; qot++)
	{
		workingData.push_back(prev->evaluate(qot+1));
	}
	double temp = 0;
	for(int j = 1; j <= dataSize; j++)
	{
		for(int k = 1; k <= dataSize - j + 1; k++)
		{
			temp = calcSingleGain(j, k);
			workingData[j-1] -= temp;
			workingData[k-1] -= temp;
			if(j+k <= dataSize)
			{
				workingData[j+k-1] += temp;
			}
		}
	}
	for(int ntk = 0; ntk < dataSize; ntk++)
	{
		if(workingData[ntk] < epsilon)
		{
			workingData[ntk] = epsilon;
		}
	}
	if(dataSize < MAX_SOLUTIONS && workingData[dataSize-1] > THRESHHOLD)
	{
		int newSize = dataSize * 2;
		if(newSize > MAX_SOLUTIONS)
		{
			newSize = MAX_SOLUTIONS;
		}
		workingData.reserve(newSize);
		for(int zvb = 0; zvb < newSize-dataSize; zvb++)
		{
			workingData.push_back(0.0);
		}
		dataSize = newSize;
	}
	if(dataSize > MAX_N/2)
	{
		solver->increaseMaxN();
	}
	
	outData = new double[dataSize];
	for(int cop = 0; cop < dataSize; cop++)
	{
		outData[cop] = workingData[cop];
	}
	
	bEvaluationReady = true;
	return;
}

int admSolution::get_maxN()
{
	//Since we solve for all mass numbers from 1 to N, our highest N is equal to the size of the solution
	return dataSize;
}

void admSolution::set_z(double zIn)
{
	z = zIn;
	x = x_of_z(z);
	h = z - prev->z;
	return;
}

double admSolution::evaluate(int n)
{
	if(bEvaluationReady==true && n >= 1 && n <= dataSize)
		return outData[n-1];
	else
		return 0.0;
	
}

double admSolution::evaluateDelta(int n)
{
	if(h < epsilon)
		return 0.0;
	
	return ((evaluate(n) - prev->evaluate(n))/h);
}

double admSolution::nFunc(int n)
{
	if(n >= 1 && n <= dataSize && outData[n-1] > epsilon)
		return outData[n-1];
	else
		return 0.0;
	
}

admSolver::admSolver(double m, double g, double n, double MAX_Nin)
{
	mx = m;
	gp = g;
	ns = n;
	MAX_N = MAX_Nin;
	xMax = 0;
	xMin = 1e200;
	double narray[35] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 100.0, 200.0, 400.0, 700.0, 1000.0, 2000.0, 4000.0, 7000.0, 10000.0, 20000.0, 40000.0,
						70000.0, 100000.0,200000.0, 4e5, 7e5, 1e6, 2e6, 4e6, 7e6, 1e7, 2e7, 4e7, 7e7, 1e8, 2e8, 4e8, 7e8};
	for(int k = 0; (k < 35 && narray[k] < 3e5); k++)
	{
		nums.push_back(narray[k]);
	}
}

void admSolver::addADM(double T)
{
	ADMmodel admTemp(mx, gp, T);
	if(admTemp.D < nums.size())
	{
		dvector::iterator temp = nums.begin();
		while(temp != nums.end())
		{
			printf("calculating for n = %d\n", (int)(*temp + 0.5));
			admTemp.getParticleProps(*temp);
			temp++;
		}
	}
	if(admTemp.prodRates.size() < (admTemp.D*(admTemp.D - 1))/2)
		admTemp.getAllProductionRates(T);
	
	Models.push_back(admTemp);
	admTemp.writeBoundProps(admTemp.T);
	xvals.push_back(x_of_T(T));
	if(x_of_T(T) > xMax)
	{
		xMax = x_of_T(T);
		zMax = z_of_x(xMax);
	}
	if(x_of_T(T) < xMin)
	{
		xMin = x_of_T(T);
		zMin = z_of_x(xMin);
	}
	return;
}

double admSolution::calcSingleGain(int n1, int n2)
{
	if(n1 < 1 || n1 > 2*MAX_N || n2 < 1 || n2 > 2*MAX_N)
		return 0.0;
	
	double result = (prev->nFunc(n1))*(prev->nFunc(n2))*getInterp(n1, n2);
	result *= h*x*x * Ns / (H_of_m(T_of_x(x)) * entropy(T_of_x(x)));
	
	return result;
}

double admSolution::getInterp(double xIn, double yIn)
{
	return solver->getApproxInterp((double)xIn, (double)yIn, prev->x);
}

void admSolver::setApproxInterp()
{
	lineBreak();
	lineBreak();
	printf("Calculating interpolation...\n");
	dvector xOut;
	dvector n1Out;
	dvector n2Out;
	dvector ratesOut;
	for(int j = 0; j < Models.size(); j++)
	{
		for(int k = 0; (k < nums.size() && nums[k]<=MAX_N); k++)
		{
			for(int l = 0; l <= k; l++)
			{
				xOut.push_back(log(xvals[j]));
				n1Out.push_back(log(nums[k]));
				n2Out.push_back(log(nums[l]));
				ratesOut.push_back(log(Models[j].getRateFromSpline(nums[k], nums[l])));
			}
		}
	}
	approxInterp.reset(n1Out, n2Out, xOut, ratesOut, 10, 0.0);
	lineBreak();
	return;
}

double admSolver::getApproxInterp(double n1, double n2, double x)
{
	if((n1 < 1.0 || n1 > MAX_N) || (n2 < 1.0 || n2 > MAX_N) || (x < xMin || x > xMax))
		return 0.0;
	
	if(n2 > n1)
	{
		double temp = n2;
		n2 = n1;
		n1 = temp;
	}
	return exp(approxInterp.eval(log(n1), log(n2), log(x)));
}

void admSolver::initialize()
{
	int q = -15;
	while(vecGen(q)<1200)
	{
		this->addADM(T_from_t(vecGen(q)));
		q++;
	}
	printf("xMin = %.3e, xMax = %.3e\n\n", xMin, xMax);
	this->setApproxInterp();
	Solutions.clear();
	int nStart = MAX_N;
	Solutions.push_back(admSolution(zMin, this));
	return;
}

void admSolver::increaseMaxN()
{
	MAX_N = 5*MAX_N;
	if(MAX_N > 2*MAX_SOLUTIONS + 10)
	{
		MAX_N = 2*MAX_SOLUTIONS + 10;
	}
	setApproxInterp();
	return;
}

double estimateAbsError(admSolution * sol1, admSolution * sol2)
{
	int siz = sol1->get_maxN();
	if(sol2->get_maxN() > siz)
	{
		siz = sol2->get_maxN();
	}
	double error = 0.0;
	double logSum = 0.0;
	double temp;
	for(int k = 1; k <= siz; k++)
	{
		if((sol1->evaluate(k) > THRESHHOLD && sol2->evaluate(k) > THRESHHOLD))
		{
			temp = fabs((sol1->evaluate(k) - sol2->evaluate(k))*log(10.0*(sol1->evaluate(k) + sol2->evaluate(k))/THRESHHOLD));
			error += temp*temp;
			logSum += pow(log(10.0*(sol1->evaluate(k) + sol2->evaluate(k))/THRESHHOLD),2);
		}
	}
	return sqrt(error)/sqrt((double)siz*logSum);
}

double estimateRelError(admSolution * sol1, admSolution * sol2)
{
	int siz = sol1->get_maxN();
	if(sol2->get_maxN() > siz)
	{
		siz = sol2->get_maxN();
	}
	double error = 0.0;
	double temp1 = 0.0;
	double temp2 = 0.0;
	for(int k = 1; k <= siz; k++)
	{
		temp1 = sol1->evaluate(k);
		temp2 = sol2->evaluate(k);
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

void admSolver::solveForward()
{
	printf("\n");
	printf("Solving from x = %.2e to %.2e:\n\n", xMin, xMax);
	double zTemp = zMin;
	const double zFin = zMax;
	double baseStep = 1e-50;
	double baseStepMax = 0.5;
	double baseStepMin = baseStep;
	double relTol = 1e-9;
	double absTol = 1e-9;
	double outStep;
	double error;
	size_t count = 0;
	int k = 0;
	int nSiz = 0;
	double numCount = 0;
	double logCount = 0;
	double temp = 0.0;
	string bleh;
	
	auto_ptr<admSolution> k1;
	auto_ptr<admSolution> k2;
	auto_ptr<admSolution> k3;
	auto_ptr<admSolution> kTemp;
	auto_ptr<admSolution> kError;
	auto_ptr<admSolution> kResult;
	clock_t clockS = clock();
	
	printf("k = 0\n\n");
	clock_t loopClockS;
	clock_t loopClockE;
	while(zTemp < zFin)
	{
		loopClockS = clock();
		if(METHOD == 1)
		{
			//The Heun Method, order 2.
			k1.reset(new admSolution(&Solutions[k], zTemp, baseStep));
			zTemp += baseStep;
			k2.reset(new admSolution(k1.get(), zTemp + baseStep, baseStep));
			k2->dataSize = k1->dataSize;
			
			kResult.reset(new admSolution(k1.get()));
			kResult->set_z(zTemp);
			kResult->addSolution(k1.get(), -0.5*baseStep, true);
			kResult->addSolution(k2.get(), 0.5*baseStep, true);
			
			
			Solutions.push_back(*(kResult.release()));
			k++;
			Solutions[k].analyzeSolution();
			
			error = estimateAbsError(&Solutions[k], k1.get());
			
			printf("Error = %.5e\n", error);
			if(error > (absTol * 0.01))
			{
				baseStep *= sqrt(absTol / error);
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
			k1.reset(new admSolution(&Solutions[k], zTemp + baseStep, baseStep));
			if(k%70 == 0 && absTol < 5e-7)
				absTol *= 2.0;
		}
		else if(METHOD == 2)
		{
			// The Bogacki-Shampine method, order 3.
			if(k == 0)
			{
				k1.reset(new admSolution(&Solutions[k], zTemp, 0.5*baseStep));
			}
			k2.reset(new admSolution(k1.get(), zTemp + 0.5*baseStep, 0.5*baseStep));
			
			kTemp.reset(new admSolution(k1.get()));
			kTemp->addSolution(k1.get(), -0.5*baseStep, true);
			kResult.reset(new admSolution(k1.get()));
			kResult->addSolution(k1.get(), -0.5*baseStep, true);
			kError.reset(new admSolution(k1.get()));
			kError->addSolution(k1.get(), -0.5*baseStep, true);
			
			kTemp->set_z(zTemp + 0.75*baseStep);
			kTemp->addSolution(k2.get(), 0.75*baseStep, true);
			k3.reset(new admSolution(kTemp.get(), kTemp->x, 0.5*baseStep));
			kTemp->addSolution(k2.get(), -0.75*baseStep, true);
			
			kResult->addSolution(k1.get(), 2.0*baseStep/9.0, true);
			kResult->addSolution(k2.get(), baseStep/3.0, true);
			kResult->addSolution(k3.get(), 4.0*baseStep/9.0, true);
			
			zTemp += baseStep;
			kResult->set_z(zTemp);
			
			Solutions.push_back(*(kResult.release()));
			k++;
			Solutions[k].analyzeSolution();
			
			kError->addSolution(k1.get(), 7.0*baseStep/24.0, true);
			kError->addSolution(k2.get(), baseStep/4.0, true);
			kError->addSolution(k3.get(), baseStep/3.0, true);
			
			k1.reset(new admSolution(&Solutions[k], zTemp, 0.5*baseStep));
			kError->addSolution(k1.get(), baseStep/8.0, true);
			kError->set_z(zTemp);
			
			error = estimateAbsError(&Solutions[k], kError.get());
			
			printf("Error = %.5e\n", error);
			if(error > (absTol * 0.01))
			{
				baseStep *= sqrt(absTol / error);
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
			k1.reset(new admSolution(&Solutions[k], zTemp + baseStep, baseStep));
			if(k%70 == 0 && absTol < 5e-7)
				absTol *= 2.0;
		}
		else if(METHOD == 3)
		{
			//Backward Euler Method, implicit method of order 1
			if(k == 0)
			{
				k1.reset(new admSolution(&Solutions[k], zTemp, baseStepMin));
			}
			else
			{
				k1.reset(new admSolution(&Solutions[k]));
			}
			k2.reset(new admSolution(k1.get(), zTemp+baseStep, baseStep));
			nSiz = k1->get_maxN();
			gsl_matrix * gradMat = gsl_matrix_alloc(nSiz, nSiz);
			for(int q1 = 0; q1 < nSiz; q1++)
			{
				temp = k1->data[q];
				k1->data[q] = temp*1.01;
				kTemp.reset(new admSolution(k1.get(), zTemp+baseStep, baseStep));
				k1->data[q] = temp;
				kTemp.dataSize = nsiz;
				kTemp->addSolution(k1.get(), -1.0, false);
				for(int q2 = 0; q2 < nSiz; q2++)
				{
					
				}
			}
			gsl_matrix_free(gradMat);
		}
		
		//if(k%500 == 0)
			//cin >> bleh;
		loopClockE = clock();
		printf("Step completed in %.3e seconds.\n", ((double)(loopClockE - loopClockS)/CLOCKS_PER_SEC));
		lineBreak();
	}
	clock_t clockE = clock();
	printf("\nEnded when xTemp = %.3e. Time was %.3e\n\n", zTemp, ((double)clockE-clockS)/CLOCKS_PER_SEC);
	return;
}

int main()
{
	srand(time(NULL));
	
	unsigned int fileInd = 1;
	char buffer[40];
	char outName[10];
	if(METHOD == 1)
	{
		strcpy(outName, "Out_Heun");
	}
	else if(METHOD == 2)
	{
		strcpy(outName, "Out_BogS");
	}
	else if(METHOD == 3)
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
		sprintf(buffer, "%.8s_M%.1e_G%.1e_N%.1e_R%.3u.csv", outName, PARTICLE_MASS, COUPLING_CONST, Ns, fileInd);
		fs.open(buffer, std::fstream::in);
		if(fs.good() == false)
		{
			fs.close();
			strcpy(fileName, buffer);
			break;
		}
		else
		{
			fs.close();
			fileInd++;
		}
	}
	
	admSolver solver(PARTICLE_MASS, COUPLING_CONST, Ns, MAX_NUMBER);
	
	solver.initialize();
	
	solver.solveForward();
	
	return 0;
}
#endif