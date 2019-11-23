#ifndef MASSCALC_H
#define MASSCALC_H

#include <stdio.h>
#include <math.h>

#ifndef GSLMATH_H
	#define GSLMATH_H
	#include <gsl/gsl_math.h>
#endif

#ifndef SPLINE_H
	#define SPLINE_H
	#include <gsl/gsl_spline.h>
#endif

#ifndef ODEIV2_H
	#define ODEIV2_H
	#include <gsl/gsl_odeiv2.h>
#endif

#ifndef ROOTS_H
	#define ROOTS_H
	#include <gsl/gsl_roots.h>
#endif

#ifndef DERIV_H
	#define DERIV_H
	#include <gsl/gsl_deriv.h>
#endif

#ifndef ERRNO_H
	#define ERRNO_H
	#include <gsl/gsl_errno.h>
#endif

#ifndef INTEGRATION_H
	#define INTEGRATION_H
	#include <gsl/gsl_integration.h>
#endif

#ifndef GSLMIN_H
	#define GSLMIN_H
	#include <gsl/gsl_min.h>
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef MASSCALCFUNCS_H
#define MASSCALCFUNCS_H
#include <time.h>
#include <vector>

using namespace std;

struct outDat
{
	double mass;
	double energy;
	double radius;
	double phi0;
	double N;
	double epsilon;
	std::vector<double> massVec;
	std::vector<double> eVec;
	std::vector<double> rVec;
	std::vector<double> nVec;
	std::vector<double> pFermiVec;
	std::vector<double> phiVec;
};

struct fullIntegrandStruct
{
	double gp;
	double m;
	double yp;
};

class pFermiC
{
	public:
		pFermiC(int, double);
		pFermiC() {};
		~pFermiC();
		double evalWithR(double, double);
		double eval(double);
		void setR(double);
		void setConst(bool inB) {bUseConst = inB; return;};
		double getR();
		double N;
	
	protected:
		bool bUseConst;
		double R0;
		double R;
		gsl_interp_accel *acc;
		gsl_spline *spline;
		double coef;
};

struct odeFData
{
	pFermiC pF;
	double gp;
	double mx;
	double R;
	double pox;
	gsl_odeiv2_system sys;
	gsl_odeiv2_driver *d;
	gsl_root_fdfsolver *s;
	gsl_root_fsolver *sF;
	bool bIntMode;
	gsl_spline *spline2;
	gsl_interp_accel *acc2;
	gsl_integration_workspace *w;
	std::vector<double> outVec;
	std::vector<double> rVec;
};

double iFun(double);

double hFun(double);

double solveHydroFunc (double, void*);

double dSolveHydroFunc(double, void*);

void hydroFDF(double, void*, double*, double*);

int hydroFunc (double, const double[], double[], void*);

double hydroIntegralF(double, void*);

double solveE(double, void*);

double dSolveE(double, void*);

double d2SolveE(double, void*);

void dd2SolveE(double, void*, double*, double*);

struct outDat massCalc(double, double, double);

struct outDat massCalc(double, double, double, bool);

int testFullIntegralFunc(double, const double[], double[], void *);

double yIntFunc(double, void *);

double testEndCondition(double, double, double, double, double);

double testFullIntegralRun(double, double, double, double);

struct outDat testFullIntegralFinal(double, double, double, double);

struct outDat completeMassIntegral(double, double, double);

/*int main(int, char*[]);*/

#endif
#endif