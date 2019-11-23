#ifndef ADM_MODEL_H
#define ADM_MODEL_H


#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <sstream>

/*#ifndef DVECTOR_H
	#include "dvector.h"
	#define DVECTOR_H
#endif*/
/*#include <dlib/threads.h>*/

#ifndef PRODUCTIONRATE_H
	#include "productionRate.h"
	#define PRODUCTIONRATE_H
#endif

#ifndef MASSCALC_H
	#include "massCalc.h"
	#define MASSCALC_H
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef TWO_DIM_INTERP_H
	#include "twoDimInterp.h"
	#define TWO_DIM_INTERP_H
#endif

using namespace std;

#ifndef VECTORS_INIT
	typedef std::vector<double> dvector;
	typedef std::vector<int> intvector;
	#define VECTORS_INIT
#endif

class ADMmodel
{
	public:
		ADMmodel(double, double, double);
		void writeBoundProps(double);
		void readBoundProps(double);
		void getParticleProps(double);
		void getAllParticleProps(double);
		void makeSplines(double);
		void getProductionRates(double);
		void getDecayRates(double);
		void getAllProductionRates(double);
		void getAllDecayRates(double);
		dvector getRatesForN(double);
		dvector getDecaysForN(double);
		dvector getAllRatesForN(double);
		double getMass(double);
		double getE(double);
		double getR(double);
		double getPhi(double);
		void display_getMass(double);
		void display_getE(double);
		void display_getR(double);
		void display_getPhi(double);
		double gp;
		double mx;
		double T;
		double x;
		dvector num;
		double numRange[2];
		int D;
		dvector prodRates;
		dvector decayRates;
		double getRateFromSpline(double, double);
		double getDecayFromSpline(double, double);
		void explicitRateCalculation(int, dvector &, dvector &, double);
		double singleRateCalculation(int, int, double);
		double singleDecayCalculation(int, int, double);
	protected:
		double lambda;
		string crosssection_path;
		string boundstate_properties_path;
		string decay_path;
		dvector masses;
		dvector bindEnergies;
		dvector radii;
		dvector phi0;
		gsl_spline *massSpline;
		gsl_spline *ESpline;
		gsl_spline *radSpline;
		gsl_spline *phiSpline;
		gsl_interp_accel *mAcc;
		gsl_interp_accel *eAcc;
		gsl_interp_accel *rAcc;
		gsl_interp_accel *pAcc;
		twoDimInterp rateSpline;
		twoDimInterp decaySpline;
		double maxN;
		bool bMadeSplines;
};
#endif