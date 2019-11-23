#ifndef TWO_DIM_INTERP_H
#define TWO_DIM_INTERP_H

#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>

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

class twoDimInterp
{
	public:
		twoDimInterp(dvector, dvector, dvector, double);
		twoDimInterp();
		void reset(dvector, dvector, dvector, double);
		double eval(double, double);
		double eval(double *);
		gsl_vector * y;
		gsl_matrix * Nmat;
		gsl_matrix * Mmat;
		gsl_matrix * Minv;
		gsl_vector * Acoef;
		gsl_vector * Bcoef;
		double * a;
		double * b;
		size_t m;
		double lambda;
	private:
		void calcNMat();
		void calcMMat();
		void setYvec(dvector);
		void calcB();
		void calcA();
		bool bReady;
		dvector xIn;
		dvector yIn;
};

#endif