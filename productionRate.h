#ifndef PRODUCTIONRATE_H
#define PRODUCTIONRATE_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#ifndef GSLMATH_H
#define GSLMATH_H
#include <gsl/gsl_math.h>
#endif

#ifndef MONTE_H
#define MONTE_H
#include <gsl/gsl_monte.h>
#endif

#ifndef MONTEVEGAS_H
#define MONTEVEGAS_H
#include <gsl/gsl_monte_vegas.h>
#endif

#ifndef GSLBESSEL_H
#define GSLBESSEL_H
#include <gsl/gsl_sf_bessel.h>
#endif

#ifndef COMPLEX_H
	#include <gsl/gsl_complex.h>
	#include <gsl/gsl_complex_math.h>
	#define COMPLEX_H
#endif

#ifndef INTEGRATION_H
	#define INTEGRATION_H
	#include <gsl/gsl_integration.h>
#endif

#include <gsl/gsl_sf_gamma.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

struct myBoost {
  double pRel;
  double multFac;
};

struct intData {
  double mp;
  double mk;
  double m;
  double n;
  double bindEnergy;
  double gp;
  double T;
  double a;
  double b;
  double pRel;
  double M;
};

struct outDataStruct {
	double prodOut;
	double decayOut;
};

struct myBoost findCMFrame (double p, double mp, double k, double mk, double theta);

double phaseFactor(double p, double m, double T);

outDataStruct productionCrossSectionCalc(double p, double n, double m, double M_n, double M_m, double bindEnergy, double gx);

double productionCrossSection(double p, double n, double M_n, double k, double m, double M_m, double theta, double bindEnergy, double gp);

double integrandFun(double *x, size_t dim, void *fdata);

double productionRate (double n, double mp, double m, double mk, double bindEnergy, double gp, double T);

double decayRate (double n, double mp, double m, double mk, double bindEnergy, double gp, double T);

double testFun1(double gamma, void * params);
double testFun2(double v, void * params);

/*int main();*/

#endif