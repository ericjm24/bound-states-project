double COUPLING_CONST = .25;
double PARTICLE_MASS = 100.0;
double X_VAR = 1.0;
double MIN_N = 1e0;
double MAX_N = 2e20 * MIN_N;
int METHOD = 1;
const char PATH_STRING[] = "F:\\cygwin\\home\\Eric\\particles\\";
#include <stdio.h>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "massCalc.h"
#include <windows.h>

using namespace std;

typedef std::vector<double> dvector;

double singleCalc(double gp, double mx, double x0, bool bSave, char pathIn[])
{
	char buffer[150];
	sprintf(buffer, "%s\\particle_M%.1e_G%.1e_X%.1e.csv", pathIn, mx, gp, x0);
	struct outDat out = completeMassIntegral(x0, gp, mx);
	printf("\nN = %.2e, R = %.4e, M = %.4e, E = %.4e\n\n", out.N, out.radius, out.mass, out.energy);
	if(bSave)
	{
		std::fstream fs;
		fs.open(buffer, std::fstream::out|std::fstream::trunc);
		for(int k = 0; k < out.rVec.size(); k++)
		{
			fs << out.rVec[k] << ",";
			fs << out.massVec[k] << ",";
			fs << out.pFermiVec[k] << ",";
			fs << out.nVec[k] << ",";
			fs << out.eVec[k] << ",";
			fs << out.phiVec[k] << "\n";
		}
		fs.close();
	}
	return out.N;
}

void calcApprox(double mx, double gp, bool bUseConst)
{
	double n = MIN_N;
	printf("Mx = %.3e, Gp = %.4e\n\n", mx, gp);
	char buffer[100];
	sprintf(buffer, "%sG%.2e_M%.2e", PATH_STRING, gp, mx);
	SetCurrentDirectory(buffer);
	char buffer2[200];
	sprintf(buffer2, "%s\\G%.2e_M%.2e\\a__G%.2e_M%.2e.csv",PATH_STRING, gp, mx, gp, mx);
	std::fstream fs;
	while(n <= MAX_N)
	{
		
		struct outDat out = massCalc(n, mx, gp, bUseConst);
		if(out.energy > 0.0)
		{
			printf("For n = %.5e, R = %.5e, E = %.5e, M = %.5e\n", n, out.radius, out.energy, out.mass);
			fs.open(buffer2, std::fstream::out|std::fstream::app);
			fs << n << ", ";
			fs << out.radius << ", ";
			fs << out.energy << ", ";
			fs << out.mass << "\n";
			fs.close();
		}
		
		n *= 2.0;
	}
	return;
}

void calcX0(double gp, double mx)
{
	char buffer[100];
	sprintf(buffer, "%sG%.2e_M%.2e", PATH_STRING, gp, mx);
	printf("\n\n%s\n\n", buffer);
	CreateDirectory(buffer, NULL);
	SetCurrentDirectory(buffer);
	char buffer2[200];
	sprintf(buffer2, "%s\\G%.2e_M%.2e\\a__G%.2e_M%.2e.csv",PATH_STRING, gp, mx, gp, mx);
	std::fstream fs;
	//fs.open(buffer2, std::fstream::out|std::fstream::trunc);
	//fs.close();
	double x0 = 1.0;
	double temp = singleCalc(gp, mx, x0, false, buffer);
	//calcApprox(mx, gp, false);
	//calcApprox(mx, gp, true);
	if (temp > MIN_N)
	{
		while(temp > MIN_N)
		{
			x0 /= 2.0;
			temp = singleCalc(gp, mx, x0, false, buffer);
		}
	}
	if (temp < MIN_N)
	{
		while(true)
		{
			x0 *= 2.0;
			temp = singleCalc(gp, mx, x0, false, buffer);
			if(temp > MIN_N)
			{
				x0 /= 2.0;
				temp = singleCalc(gp, mx, x0, false, buffer);
				break;
			}
		}
	}
	double multFac = 1.5;
	double oldN = 0.0;
	do
	{
		x0 *= multFac;
		temp = singleCalc(gp, mx, x0, true, buffer);
		if(temp < oldN*1.005 && temp > oldN)
		{
			break;
		}
		else
		{
			oldN = temp;
		}
	}while(temp < MAX_N && x0 < 1e200);
	
	return;
}

void calcGP(double x0, double mx)
{
	char buffer[100];
	sprintf(buffer, "%sX%.2e_M%.2e", PATH_STRING, x0, mx);
	printf("\n\n%s\n\n", buffer);
	CreateDirectory(buffer, NULL);
	SetCurrentDirectory(buffer);
	double gp = 1.0;
	double temp = singleCalc(gp, mx, x0, false, buffer);
	
	if (temp > MIN_N)
	{
		while(temp > MIN_N)
		{
			gp *= 2.0;
			temp = singleCalc(gp, mx, x0, false, buffer);
		}
	}
	if (temp < MIN_N)
	{
		while(true)
		{
			gp /= 2.0;
			temp = singleCalc(gp, mx, x0, false, buffer);
			if(temp > MIN_N)
			{
				gp *= 2.0;
				temp = singleCalc(gp, mx, x0, false, buffer);
				break;
			}
		}
	}
	double multFac = 1.2;
	do
	{
		gp /= multFac;
		temp = singleCalc(gp, mx, x0, true, buffer);
	}while(temp < MAX_N);
	return;
}

void calcMX(double x0, double gp)
{
	char buffer[100];
	sprintf(buffer, "%sG%.2e_X%.2e", PATH_STRING, gp, x0);
	printf("\n\n%s\n\n", buffer);
	CreateDirectory(buffer, NULL);
	SetCurrentDirectory(buffer);
	double mx = 1e-3;
	double multFac = 2.0;
	double temp;
	do
	{
		mx *= multFac;
		temp = singleCalc(gp, mx, x0, true, buffer);
	}while(mx < 1e9);
	return;
}

int main()
{
		double gp = COUPLING_CONST;
		double mx = PARTICLE_MASS;
		double x0 = X_VAR;
	switch(METHOD)
	{
		case 1:
			calcX0(gp, mx);
			break;
		case 2:
			calcGP(x0, mx);
			break;
		case 3:
			calcMX(x0, gp);
			break;
		case 4:
			calcApprox(mx, gp, false);
			calcApprox(mx, gp, true);
			break;
	}
	return 0;
}