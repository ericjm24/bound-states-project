#ifndef ADM_SOLVER
#define ADM_SOLVER

double PARTICLE_MASS = 140.0; //GeV
double COUPLING_CONST = 4.0;
double Ns = 4e-7; //number density per comoving volume
double MAX_NUMBER = 100000;
double mPlanck = 1.2209e19; //Planck Mass, units of GeV/c^2 (c = 1)
double Trec = 3.08e-9; //Recombination temperature, everything should be frozen out by now. In units of GeV.
int MAX_SOLUTIONS = 30; //Largest number of independent solutions to be kept track of; will move to splines after this;
double THRESHHOLD = 1e-16; //Any particle that has a proportion of the total particle number smaller than this is ignored;
double epsilon = 1e-60; //Treat this as the smallest number that can exist. Is basically 0.

#include <stdio.h>
#include <math.h>
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
typedef std::pair<double, double> point;
typedef std::vector<point> locus;

enum TAG {ARRAY, SPLINE};

struct splineStruct
{
	gsl_spline * spline;
	gsl_interp_accel * acc;
	double maxN;
	double minN;
	double * xData;
	double * yData;
};

struct instruction
{
	TAG tag;
	splineStruct spl;
	locus loc;
};

typedef std::vector<instruction> recipe;

class admSolution
{
	public:
		admSolution(double);
		admSolution(admSolution *, double);
		double get_max();
		double get_min();
		double get_minN();
		double get_maxN();
		double get_xData();
		void solve(admSolution *);
		void solve();
	private:
		std::vector<splineStruct> spl;
		dvector discX;
		dvector discY;
		double x;
		admSolution * prev;
		admSolution * next;
		recipe createRecipe();
};

admSolution::admSolution(double xIn)
{
	x = xIn;
}

admSolution::admSolution(admSolution * prevIn, double xIn)
{
	x = xIn;
	prev = prevIn;
	solve();
}

void admSolution::solve(admSolution * prevIn)
{
	prev = prevIn;
	solve();
	return;
}

void admSolution::solve()
{
	recipe rec;
	
	return;
}

class OutStruct
{
	public:
		TAG tag;
		double x;
		double get_max();
		double get_U_limit();
		splineStruct spl;
		dvector array;
		void set_array(dvector);
		void set_spline(splineStruct);
		double U_limit;
		double max;
};

void OutStruct::set_array(dvector vecIn)
{
	tag = ARRAY;
	array = vecIn;
	max = get_max();
	U_limit = get_U_limit();
	return;
}

void OutStruct::set_spline(splineStruct spIn)
{
	tag = SPLINE;
	spl = spIn;
	max = get_max();
	U_limit = get_U_limit();
	return;
}

double OutStruct::get_max()
{
	if(tag == ARRAY)
	{
		return *(std::max_element(array.begin(), array.end()));
	}
	else if(tag == SPLINE)
		{
			double out = 0.0;
			double temp;
			for(double k = 1.0; k <= spl.maxN; k++)
			{
				//printf("\n\nk = %.0f\n\n", k);
				cout << flush;
				temp = exp(gsl_spline_eval(spl.spline, k, spl.acc));
				if(temp > out)
					out = temp;
			}
			return out;
		}
	return 0.0;
}

double OutStruct::get_U_limit()
 {
	if(tag == ARRAY)
		return array.size();
	else if(tag == SPLINE)
		return spl.maxN;
	return 0.0;
 }

class admSolver
{
	public:
		admSolver(double, double, double, double);
		void addADM(double);
		void initialize();
		double getApproxInterp(double, double, double);
		void solveForward();
		dvector xvals;
		double xMax;
		double xMin;
		double MAX_N;
	protected:
		double calcGainIntegral(const OutStruct &, double);
		double calcLossIntegral(const OutStruct &, double);
		double calcTotalDelta(const OutStruct &, double);
		double stepForward(double);
		modvector Models;
		threeDimInterpApprox approxInterp;
		dvector nums;
		double mx;
		double gp;
		double ns;
		std::vector<OutStruct> outData;
		void setApproxInterp();
		bool bStep;
};

void lineBreak()
{
	printf("****************************************\n");
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

admSolver::admSolver(double m, double g, double n, double MAX_Nin)
{
	mx = m;
	gp = g;
	ns = n;
	MAX_N = MAX_Nin;
	xMax = 0;
	xMin = 1e200;
	int temp = 0;
	double narray[35] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 100.0, 200.0, 400.0, 700.0, 1000.0, 2000.0, 4000.0, 7000.0, 10000.0, 20000.0, 40000.0,
						70000.0, 100000.0,200000.0, 4e5, 7e5, 1e6, 2e6, 4e6, 7e6, 1e7, 2e7, 4e7, 7e7, 1e8, 2e8, 4e8, 7e8};
	do
	{
		nums.push_back(narray[temp]);
		temp++;
	}while(narray[temp] <= MAX_N);
	nums.push_back(MAX_N*2);
}

void admSolver::addADM(double T)
{
	ADMmodel admTemp(mx, gp, T);
	if(admTemp.D < nums.size())
	{
		dvector::iterator temp = nums.begin();
		while(temp != nums.end())
		{
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
		xMax = x_of_T(T);
	
	if(x_of_T(T) < xMin)
		xMin = x_of_T(T);
	
	return;
}

double nFunc(const OutStruct &tempS, double n)
{
	double nOut = 0.0;
	switch(tempS.tag)
	{
		case ARRAY:
			if(n <= tempS.array.size() && n >= 1.0)
				nOut = tempS.array[(int)(n - 0.9)];
			break;
		case SPLINE:
			if(n <= tempS.spl.maxN && n >= 1.0)
				nOut = exp(gsl_spline_eval(tempS.spl.spline, n, tempS.spl.acc));
			break;
	}
	if(nOut < 0.0)
		nOut = 0.0;
	return nOut;
}

double admSolver::calcGainIntegral(const OutStruct &tempS, double n)
{
	if(n<=1 || n > MAX_N)
		return 0;
	
	double gain_result = 0;
	for(double k = 1.0; k <= n-k; k++)
	{
		if(k == n - k)
			gain_result += 0.5*nFunc(tempS, k)*nFunc(tempS, n - k)*getApproxInterp(k, n-k, tempS.x);
		else
			gain_result += nFunc(tempS, k)*nFunc(tempS, n - k)*getApproxInterp(k, n-k, tempS.x);
	}
	return gain_result;
}

double admSolver::calcLossIntegral(const OutStruct &tempS, double n)
{
	if((n < 1 || n > MAX_N) || nFunc(tempS, n) <= epsilon)
		return 0;
	
	double loss_result = 0;
	for(double k = 1.0; k <= tempS.U_limit; k++)
	{
			loss_result += nFunc(tempS, k)*nFunc(tempS, n)*getApproxInterp(k, n, tempS.x);
	}
	return loss_result;
}

double admSolver::calcTotalDelta(const OutStruct &tempS, double n)
{
	double result = calcGainIntegral(tempS, n) - calcLossIntegral(tempS, n);
	result *= tempS.x / (H_of_m(T_of_x(tempS.x)) * entropy(T_of_x(tempS.x)));
	return result;
}

void distributeN(const OutStruct &tempS, double * baseVec, size_t numN)
{
	dvector derivVec;
	double deriveSum = 0.;
	double tempVal = 0;
	switch(tempS.tag)
	{
		case ARRAY:
			
			for(dvector::const_iterator tempIt = tempS.array.begin()+1; tempIt != tempS.array.end() - 1; tempIt++)
			{
				deriveSum += sqrt(0.5*fabs(*(tempIt+1) - *(tempIt - 1)));
				derivVec.push_back(deriveSum);
			}
			deriveSum += sqrt(fabs(*(tempS.array.end()-1) - *(tempS.array.end()-2)));
			derivVec.push_back(deriveSum);
			break;
		case SPLINE:
			for(double k = 1.0; k <= tempS.U_limit; k++)
			{
				deriveSum += sqrt(fabs(gsl_spline_eval_deriv(tempS.spl.spline, k, tempS.spl.acc)* nFunc(tempS, k)));
				derivVec.push_back(deriveSum);
			}
			break;
	}
	double baseStep = 0;//(tempS.U_limit - 1.0) / ((double)numN - 1.0);
	double varStep = *(derivVec.end()-1) / ((double)numN - 1.0);
	double j = 1.0;
	dvector::const_iterator It = derivVec.begin();
	size_t count = 1;
	double incStep = 1.0;
	double tempStep = 0;
	size_t nextStep = 0;
	baseVec[0] = 1.0;
	while(j < tempS.U_limit && count < (numN-1))
	{
		while(*(It+1) < ((double)count * varStep))
		{
			It++;
			incStep += 1.0;
		}
		j = incStep + fabs((((double)count * varStep) - *It)/(*(It+1)-*It));
		j *= 0.8;
		j += baseStep*((double)count);
		j = floor(j);
		if(j <= baseVec[count - 1])
		{
			j= baseVec[count-1]+1;
		}
		baseVec[count] = j;
		count++;
	}
	baseVec[numN-1] = tempS.U_limit;
	return;
}

double admSolver::stepForward(double baseStep)
{
	std::vector<OutStruct>::iterator tempS = outData.end() - 1;
	printf("\n\noutData is now a vector of length %d.\n", outData.size());
	cout << flush;
	OutStruct newDat;
	double outStep = 1e5;
	double index = 1.0;
	bool bContinue = true;
	double temp = 0.0;
	double temp2 = 0.0;
	if((tempS->tag == ARRAY) && (tempS->array.size() <= 2*MAX_SOLUTIONS))
	{
		//Array Path
		
		index = 1.0;
		bContinue = true;
		temp = 0.0;
		temp2 = 0.0;
		outStep = baseStep;
		dvector outArray;
		while(bContinue == true)
		{
			//Grab every eligible particle (density is above threshhold)
			temp = calcTotalDelta(*tempS, index)*baseStep;
			temp2 = temp + nFunc(*tempS, index);
			
			if(temp < epsilon)
				temp = epsilon;
			
			if(temp2 < 0.0)
				temp2 = 0.0;
				
			if((nFunc(*tempS, index)/temp) < outStep)
			{
				if(nFunc(*tempS, index) > 1e10 * epsilon)
					outStep = nFunc(*tempS, index)/temp;
			}
			if(index <= tempS->array.size())
				outArray.push_back(temp2);
			else if( temp2 > THRESHHOLD * tempS->get_max())
				outArray.push_back(temp2);
			else
				bContinue = false;
			
			index++;
		}
		newDat.x = tempS->x;
		newDat.set_array(outArray);
		printf("For x = %.3e, with a step size of %.2e, we have:\n", tempS->x, baseStep);
		for(int k2 = 0; k2 < outArray.size(); k2++)
		{
			printf("The density of n = %d is now %.5e\n", k2+1, outArray[k2]);
		}
	}
	else
	{
		//Spline Path
		newDat.tag = SPLINE;
		newDat.x = tempS->x;
		double * newX = new double[MAX_SOLUTIONS];
		double * newY = new double[MAX_SOLUTIONS];
		distributeN(*tempS, &newX[0], MAX_SOLUTIONS);
		index = 1.0;
		int count = 0;
		bContinue = true;
		temp = 0.0;
		temp2 = 0.0;
		int totalIncrease = 0;
		size_t siz = MAX_SOLUTIONS;
		while(bContinue == true && totalIncrease < 0.1*tempS->U_limit)
		{
			if(count < MAX_SOLUTIONS)
			{
				if(index < newX[count])
				{
					index ++;
					continue;
				}
			}
			
			
			//Grab every eligible particle (density is above threshhold)
			temp = calcTotalDelta(*tempS, index)*baseStep;
			temp2 = temp + nFunc(*tempS, index);
			
			if(temp < epsilon)
				temp = epsilon;
			
			if(temp2 <= epsilon)
				temp2 = epsilon;
				
			if((nFunc(*tempS, index)/temp) < outStep)
			{
				if(nFunc(*tempS, index) > 1e10 * epsilon)
					outStep = nFunc(*tempS, index)/temp;
			}
			if(index <= tempS->get_U_limit())
			{
				newY[count] = log(temp2);
				count++;
			}
			else if( temp2 > (sqrt(sqrt(THRESHHOLD)) * tempS->get_max()))
			{
				totalIncrease++;
			}
			else
				bContinue = false;
			
			index++;
		}
		
		
		if(totalIncrease > 0)
		{
			cout << flush;
			siz = MAX_SOLUTIONS + totalIncrease;
			double * newX2 = new double[siz];
			double * newY2 = new double[siz];
			double tbleh = 0.0;
			for(int k = 0; k < MAX_SOLUTIONS; k++)
			{
				newX2[k] = newX[k];
				newY2[k] = newY[k];
			}
			for(int k2 = 0; k2 < totalIncrease; k2++)
			{
				newX2[k2+MAX_SOLUTIONS] = tempS->get_U_limit() + (double)k2 + 1.0; 
				tbleh = calcTotalDelta(*tempS, newX2[k2+MAX_SOLUTIONS])*baseStep;
				if(newY2[k2+MAX_SOLUTIONS] <= epsilon)
					newY2[k2+MAX_SOLUTIONS] = epsilon;
				
				newY2[k2+MAX_SOLUTIONS] = log(tbleh);
			}
			
			splineStruct sT2;
			sT2.spline = gsl_spline_alloc(gsl_interp_cspline, siz);
			sT2.acc = gsl_interp_accel_alloc();
			sT2.maxN = tempS->U_limit + (double)totalIncrease;
			gsl_spline_init(sT2.spline, newX2, newY2, siz);
			newDat.set_spline(sT2);
			distributeN(newDat, &newX[0], MAX_SOLUTIONS);
			
			outStep = baseStep;
			for(int l = 0; l < MAX_SOLUTIONS; l++)
			{
				temp = calcTotalDelta(*tempS, newX[l])*baseStep;
				temp2 = temp + nFunc(*tempS, newX[l]);
				
				if(temp < epsilon)
					temp = epsilon;
				
				if(temp2 <= epsilon)
					temp2 = epsilon;
					
				if((nFunc(*tempS, newX[l])/temp) < outStep)
				{
					if(nFunc(*tempS, newX[l]) > 1e10 * epsilon)
						outStep = nFunc(*tempS, newX[l])/temp;
				}
				newY[l] = log(temp2);
			}
			
			gsl_spline_free(newDat.spl.spline);
			gsl_interp_accel_free(newDat.spl.acc);
		}
		printf("For x = %.3e, with a step size of %.2e, we have:\n", tempS->x, baseStep);
		for(int k2 = 0; k2 < MAX_SOLUTIONS; k2++)
		{
			printf("The density of n = %.0f is now %.5e\n", newX[k2], exp(newY[k2]));
		}
		splineStruct sT;
		sT.spline = gsl_spline_alloc(gsl_interp_cspline, MAX_SOLUTIONS);
		sT.acc = gsl_interp_accel_alloc();
		sT.maxN = tempS->U_limit + (double)totalIncrease;
		gsl_spline_init(sT.spline, newX, newY, MAX_SOLUTIONS);
		newDat.set_spline(sT);
		newDat.x = tempS->x;
		
	}
	
	outData.push_back(newDat);
	return outStep;
}

void admSolver::setApproxInterp()
{
	dvector xOut;
	dvector n1Out;
	dvector n2Out;
	dvector ratesOut;
	for(int j = 0; j < Models.size(); j++)
	{
		for(int k = 0; k < nums.size(); k++)
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
	approxInterp.reset(n1Out, n2Out, xOut, ratesOut, 20, 0.0);
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
	lineBreak();
	printf("Calculating interpolation...\n");
	this->setApproxInterp();
	outData.clear();
	OutStruct newDat;
	int nStart = 30;
	dvector tarray;
	tarray.push_back(ns);
	for(int k = 0; k < nStart; k++)
	{
		tarray.push_back(0.0);
	}
	newDat.set_array(tarray);
	newDat.x = xMin;
	outData.push_back(newDat);
	return;
}

void admSolver::solveForward()
{
	printf("\n");
	lineBreak();
	printf("\n");
	printf("Solving from x = %.2e to %.2e:\n\n", xMin, xMax);
	double xTemp = xMin;
	const double xFin = xMax;
	double baseStep = 0.2;
	double baseStepMax = 5*baseStep;
	double baseStepMin = 0.1*baseStep;
	double outStep;
	size_t count = 0;
	while(xTemp < xFin)
	{
		outStep = stepForward(baseStep);
		xTemp += outStep;
		count ++;
		outData[count].x = xTemp;
		if(outStep <= baseStepMax && outStep >= baseStepMin)
			baseStep = outStep;
		
	}
	printf("\nEnded when xTemp = %.3e\n\n", xTemp);
	return;
}

int main()
{
	srand(time(NULL));
	admSolver solver(PARTICLE_MASS, COUPLING_CONST, Ns, MAX_NUMBER);
	
	solver.initialize();
	
	solver.solveForward();
	
	return 0;
}
#endif