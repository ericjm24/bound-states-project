#include "threeDimInterp.h"
#include <iostream>
using namespace std;

double getDist3(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}

threeDimInterp::threeDimInterp(dvector Xi, dvector Yi, dvector Zi, dvector ratesIn, double lambdaIn)
{
	rates = ratesIn;
	xIn = Xi;
	yIn = Yi;
	zIn = Zi;
	m = rates.size();	
	lambda = lambdaIn;
	calcPMat();
	calcKMat();
	calcYVec();
	calcMMat();
	solveABCoef();
	
	gsl_matrix_free(Kmat);
	gsl_matrix_free(Pmat);
	gsl_matrix_free(Mmat);
	gsl_vector_free(y);
	xIn.clear();
	yIn.clear();
	zIn.clear();
	rates.clear();
	bReady = true;
	return;
}

threeDimInterp::threeDimInterp()
{
	bReady = false;
	return;
}

void threeDimInterp::reset(dvector Xi, dvector Yi, dvector Zi, dvector ratesIn, double lambdaIn)
{
	rates = ratesIn;
	xIn = Xi;
	yIn = Yi;
	zIn = Zi;
	m = rates.size();
	lambda = lambdaIn;
	calcPMat();
	calcKMat();
	calcYVec();
	calcMMat();
	solveABCoef();
	
	gsl_matrix_free(Kmat);
	gsl_matrix_free(Pmat);
	gsl_matrix_free(Mmat);
	gsl_vector_free(y);
	xIn.clear();
	yIn.clear();
	zIn.clear();
	rates.clear();
	bReady = true;
	return;
}

double threeDimInterp::eval(double p, double q, double s)
{
	if(bReady == false)
		return 0.0;
	
	double t = b[0];
	t += b[1]*p;
	t += b[2]*q;
	t += b[3]*s;
	double r = 0.0;
	for(int j = 0; j < m; j++)
	{
		r = getDist3(p - xIn[j], q - yIn[j], s - zIn[j]);
		if(r > 0)
		{
			t += a[j]*r;
		}
	}
	return t;
}

double threeDimInterp::eval(double * inp)
{
	if(bReady == false)
		return 0.0;
	double p = inp[0];
	double q = inp[1];
	double s = inp[2];
	double t = b[0];
	t += b[1]*p;
	t += b[2]*q;
	t += b[3]*s;
	double r = 0.0;
	for(int j = 0; j < m; j++)
	{
		r = getDist3(p - xIn[j], q - yIn[j], s - zIn[j]);
		if(r > 0)
		{
			t += a[j]*r;
		}
	}
	return t;
}

void threeDimInterp::calcPMat()
{
	Pmat = gsl_matrix_alloc(m,4);
	for(int j = 0; j < m; j++)
	{
		gsl_matrix_set(Pmat, j, 0, 1.0);
		gsl_matrix_set(Pmat, j, 1, xIn[j]);
		gsl_matrix_set(Pmat, j, 2, yIn[j]);
		gsl_matrix_set(Pmat, j, 3, zIn[j]);
	}
	return;
}

void threeDimInterp::calcKMat()
{
	Kmat = gsl_matrix_calloc(m, m);
	double r = 0.0;
	for(int j = 0; j < m; j++)
	{
		for(int k = 0; k < m; k++)
		{
			r = getDist3(xIn[k]-xIn[j], yIn[k]-yIn[j], zIn[k] - zIn[j]);
			if(r > 0)
			{
				gsl_matrix_set(Kmat, j, k, r);
			}
		}
	}
	return;
}

void threeDimInterp::calcMMat()
{
	Mmat = gsl_matrix_calloc(m+4, m+4);
	for(int j1 = 0; j1 < m; j1++)
	{
		for(int k1 = 0; k1 < m; k1++)
		{
			gsl_matrix_set(Mmat, j1, k1, gsl_matrix_get(Kmat, j1, k1));
		}
	}
	for(int j2 = 0; j2 < m; j2++)
	{
		for(int k2 = 0; k2 < 4; k2++)
		{
			gsl_matrix_set(Mmat, j2, k2 + m, gsl_matrix_get(Pmat, j2, k2));
			gsl_matrix_set(Mmat, k2 + m, j2, gsl_matrix_get(Pmat, j2, k2));
		}
	}
	return;
}

void threeDimInterp::calcYVec()
{
	y = gsl_vector_calloc(m + 4);
	for(size_t j = 0; j < m; j++)
	{
		gsl_vector_set(y, j, rates[j]);
	}
	return;
}

void threeDimInterp::solveABCoef()
{
	gsl_permutation * p = gsl_permutation_alloc(m + 4);
	int signum = 0;
	gsl_linalg_LU_decomp(Mmat, p, &signum);
	gsl_vector * outVec = gsl_vector_calloc(m + 4);
	gsl_linalg_LU_solve(Mmat, p, y, outVec);
	a = new double[m];
	b = new double[4];
	for(int j = 0; j < m; j++)
	{
		a[j] = gsl_vector_get(outVec, j);
	}
	for(int k = 0; k < 4; k++)
	{
		b[k] = gsl_vector_get(outVec, k+m);
	}
	
	gsl_permutation_free(p);
	gsl_vector_free(outVec);
	return;
}

