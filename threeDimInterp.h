#ifndef THREE_DIM_INTERP_H
#define THREE_DIM_INTERP_H

#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <time.h>

/*#ifndef DVECTOR_H
	#include "dvector.h"
	#define DVECTOR_H
#endif*/

#ifndef GSLMATH_H
	#define GSLMATH_H
	#include <gsl/gsl_math.h>
#endif

#ifndef INTEGRATION_H
	#define INTEGRATION_H
	#include <gsl/gsl_integration.h>
#endif

#ifndef LINALG_H
	#define LINALG_H
	#include <gsl/gsl_linalg.h>
#endif

#ifndef BLAS_H
	#define BLAS_H
	#include <gsl/gsl_blas.h>
#endif

using namespace std;

#ifndef VECTORS_INIT
	typedef std::vector<double> dvector;
	typedef std::vector<int> intvector;
	#define VECTORS_INIT
#endif

class threeDimInterp
{
	public:
		threeDimInterp(dvector, dvector, dvector, dvector, double);
		threeDimInterp();
		void reset(dvector, dvector, dvector, dvector, double);
		double eval(double, double, double);
		double eval(double *);
		gsl_vector * y;
		gsl_matrix * Kmat;
		gsl_matrix * Pmat;
		gsl_matrix * Mmat;
		double * a;
		double * b;
		size_t m;
		double lambda;
	private:
		void calcKMat();
		void calcPMat();
		void calcMMat();
		void calcYVec();
		void solveABCoef();
		bool bReady;
		dvector xIn;
		dvector yIn;
		dvector zIn;
		dvector rates;
};

#endif