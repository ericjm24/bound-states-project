#include "threeDimInterpApprox.h"
#include <iostream>
using namespace std;

double getDist(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}

threeDimInterpApprox::threeDimInterpApprox(dvector Xi, dvector Yi, dvector Zi, dvector ratesIn, int sizeFactor, double lambdaIn)
{
	rates = ratesIn;
	xIn = Xi;
	yIn = Yi;
	zIn = Zi;
	p = rates.size();
	m = p/sizeFactor;
	if(p%sizeFactor != 0)
		m++;
	
	lambda = lambdaIn;
	initializePoints();
	calcPMat();
	calcAKMat();
	calcYVec();
	calcMMat();
	solveABCoef();
	
	gsl_matrix_free(Amat);
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

threeDimInterpApprox::threeDimInterpApprox()
{
	bReady = false;
	return;
}

void threeDimInterpApprox::reset(dvector Xi, dvector Yi, dvector Zi, dvector ratesIn, int sizeFactor, double lambdaIn)
{
	rates = ratesIn;
	xIn = Xi;
	yIn = Yi;
	zIn = Zi;
	p = rates.size();
	m = p/sizeFactor;
	if(p%sizeFactor != 0)
		m++;
	
	lambda = lambdaIn;
	initializePoints();
	calcPMat();
	calcAKMat();
	calcYVec();
	calcMMat();
	solveABCoef();
	
	gsl_matrix_free(Amat);
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

double threeDimInterpApprox::eval(double p, double q, double s)
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
		r = getDist(p - xIn[j], q - yIn[j], s - zIn[j]);
		if(r > 0)
		{
			t += a[j]*r;
		}
	}
	return t;
}

double threeDimInterpApprox::eval(double * inp)
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
		r = getDist(p - xIn[j], q - yIn[j], s - zIn[j]);
		if(r > 0)
		{
			t += a[j]*r;
		}
	}
	return t;
}

void threeDimInterpApprox::calcPMat()
{
	Pmat = gsl_matrix_alloc(p,4);
	for(int j = 0; j < p; j++)
	{
		gsl_matrix_set(Pmat, j, 0, 1.0);
		gsl_matrix_set(Pmat, j, 1, xIn[j]);
		gsl_matrix_set(Pmat, j, 2, yIn[j]);
		gsl_matrix_set(Pmat, j, 3, zIn[j]);
	}
	return;
}

void threeDimInterpApprox::calcAKMat()
{
	Amat = gsl_matrix_calloc(m, m);
	Kmat = gsl_matrix_calloc(p, m);
	double r = 0.0;
	for(int j = 0; j < p; j++)
	{
		for(int k = 0; k < m; k++)
		{
			r = getDist(xIn[k]-xIn[j], yIn[k]-yIn[j], zIn[k] - zIn[j]);
			if(r > 0)
			{
				gsl_matrix_set(Kmat, j, k, r);
			}
			if(r > 0 && j < m)
			{
				gsl_matrix_set(Amat, j, k, r);
			}
		}
	}
	return;
}

void threeDimInterpApprox::calcMMat()
{
	Mmat = gsl_matrix_calloc(m+4, m+4);
	gsl_matrix * MMblock = gsl_matrix_calloc(m, m);
	gsl_matrix * M4block = gsl_matrix_calloc(m, 4);
	gsl_matrix * FourMblock = gsl_matrix_calloc(4, m);
	gsl_matrix * Four4block = gsl_matrix_calloc(4, 4);
	gsl_matrix_memcpy(MMblock, Amat);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Kmat, Kmat, lambda, MMblock);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Kmat, Pmat, 0.0, M4block);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Pmat, Kmat, 0.0, FourMblock);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Pmat, Pmat, 0.0, Four4block);
	for(int j1 = 0; j1 < m; j1++)
	{
		for(int k1 = 0; k1 < m; k1++)
		{
			gsl_matrix_set(Mmat, j1, k1, gsl_matrix_get(MMblock, j1, k1));
		}
	}
	for(int j2 = 0; j2 < m; j2++)
	{
		for(int k2 = 0; k2 < 4; k2++)
		{
			gsl_matrix_set(Mmat, j2, k2 + m, gsl_matrix_get(M4block, j2, k2));
		}
	}
	for(int j3 = 0; j3 < 4; j3++)
	{
		for(int k3 = 0; k3 < m; k3++)
		{
			gsl_matrix_set(Mmat, j3+m, k3, gsl_matrix_get(FourMblock, j3, k3));
		}
	}
	for(int j4 = 0; j4 < 4; j4++)
	{
		for(int k4 = 0; k4 < 4; k4++)
		{
			gsl_matrix_set(Mmat, j4+m, k4+m, gsl_matrix_get(Four4block, j4, k4));
		}
	}
	
	gsl_matrix_free(MMblock);
	gsl_matrix_free(M4block);
	gsl_matrix_free(FourMblock);
	gsl_matrix_free(Four4block);
	return;
}

void threeDimInterpApprox::initializePoints()
{
	srand(time(NULL));
	size_t g, h;
	double temp;
	for(size_t k = 0; k < 20*p; k++)
	{
		g = k%p;
		h = rand()%p;
		temp = xIn[g];
		xIn[g] = xIn[h];
		xIn[h] = temp;
		temp = yIn[g];
		yIn[g] = yIn[h];
		yIn[h] = temp;
		temp = zIn[g];
		zIn[g] = zIn[h];
		zIn[h] = temp;
		temp = rates[g];
		rates[g] = rates[h];
		rates[h] = temp;
	}
	return;
}

void threeDimInterpApprox::calcYVec()
{
	y = gsl_vector_calloc(m + 4);
	gsl_vector * v = gsl_vector_calloc(p);
	for(size_t j = 0; j < p; j++)
	{
		gsl_vector_set(v, j, rates[j]);
	}
	gsl_vector * Mvec = gsl_vector_calloc(m);
	gsl_vector * Fourvec = gsl_vector_calloc(4);
	gsl_blas_dgemv(CblasTrans, 1.0, Kmat, v, 0.0, Mvec);
	gsl_blas_dgemv(CblasTrans, 1.0, Pmat, v, 0.0, Fourvec);
	for(size_t j2 = 0; j2 < m; j2++)
	{
		gsl_vector_set(y, j2, gsl_vector_get(Mvec, j2));
	}
	for(size_t k = 0; k < 4; k++)
	{
		gsl_vector_set(y, k + m, gsl_vector_get(Fourvec, k));
	}
	gsl_vector_free(v);
	gsl_vector_free(Mvec);
	gsl_vector_free(Fourvec);
	return;
}

void threeDimInterpApprox::solveABCoef()
{
	gsl_permutation * perm = gsl_permutation_alloc(m + 4);
	int signum = 0;
	gsl_linalg_LU_decomp(Mmat, perm, &signum);
	gsl_vector * outVec = gsl_vector_calloc(m + 4);
	gsl_linalg_LU_solve(Mmat, perm, y, outVec);
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
	
	gsl_permutation_free(perm);
	gsl_vector_free(outVec);
	return;
}

