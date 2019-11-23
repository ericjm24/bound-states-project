#include "twoDimInterp.h"

double getDist(double x, double y)
{
	return sqrt(x*x + y*y);
}

twoDimInterp::twoDimInterp(dvector Xi, dvector Yi, dvector ratesIn, double lambdaIn)
{
	dvector rates = ratesIn;
	xIn = Xi;
	yIn = Yi;
	m = rates.size();
	lambda = lambdaIn;
	setYvec(rates);
	calcNMat();
	calcMMat();
	calcB();
	calcA();
	bReady = true;
	return;
}

twoDimInterp::twoDimInterp()
{
	bReady = false;
	return;
}

void twoDimInterp::reset(dvector Xi, dvector Yi, dvector ratesIn, double lambdaIn)
{
	dvector rates = ratesIn;
	xIn = Xi;
	yIn = Yi;
	m = rates.size();
	lambda = lambdaIn;
	setYvec(rates);
	calcNMat();
	calcMMat();
	calcB();
	calcA();
	bReady = true;
	return;
}

double twoDimInterp::eval(double p, double q)
{
	if(bReady == false)
		return 0.0;
	
	double t = b[0];
	t += b[1]*p;
	t += b[2]*q;
	double r = 0.0;
	for(int j = 0; j < m; j++)
	{
		r = getDist(p - xIn[j], q - yIn[j]);
		if(r > 0)
		{
			t += a[j]*r*r*log(r);
			//printf("\nContribution from %.0f, %.0f is %e. T = %e", xIn[j], yIn[j], a[j]*r*r*log(r), t);
		}
	}
	return t;
}

double twoDimInterp::eval(double * inp)
{
	if(bReady == false)
		return 0.0;
	double p = inp[0];
	double q = inp[1];
	double t = b[0];
	t += b[1]*p;
	t += b[2]*q;
	double r = 0.0;
	for(int j = 0; j < m; j++)
	{
		r = getDist(p - xIn[j], q - yIn[j]);
		if(r > 0)
		{
			t += a[j]*r*r*log(r);
			//printf("\nContribution from %.0f, %.0f is %e. T = %e", xIn[j], yIn[j], a[j]*r*r*log(r), t);
		}
	}
	return t;
}

void twoDimInterp::calcNMat()
{
	Nmat = gsl_matrix_alloc(m,3);
	for(int j = 0; j < m; j++)
	{
		gsl_matrix_set(Nmat, j, 0, 1.0);
		gsl_matrix_set(Nmat, j, 1, xIn[j]);
		gsl_matrix_set(Nmat, j, 2, yIn[j]);
	}
	return;
}

void twoDimInterp::calcMMat()
{
	Mmat = gsl_matrix_alloc(m, m);
	double r = 0.0;
	for(int j = 0; j < m; j++)
	{
		for(int k = 0; k < m; k++)
		{
			r = getDist(xIn[k]-xIn[j], yIn[k]-yIn[j]);
			if(r<=0)
			{
				gsl_matrix_set(Mmat, j, k, 0.0);
			}
			else
			{
				gsl_matrix_set(Mmat, j, k, r*r*log(r));
			}
			if(j == k)
			{
				gsl_matrix_set(Mmat, j, k, gsl_matrix_get(Mmat, j, k) + lambda);
			}
		}
	}
	gsl_permutation * p = gsl_permutation_alloc(m);
	Minv = gsl_matrix_alloc(m, m);
	gsl_matrix * MLU = gsl_matrix_alloc(m, m);
	gsl_matrix_memcpy(MLU, Mmat);
	int s;
	gsl_linalg_LU_decomp(MLU, p, &s);
	gsl_linalg_LU_invert(MLU, p, Minv);
	gsl_matrix_free(MLU);
	gsl_permutation_free(p);
	return;
}

void twoDimInterp::setYvec(dvector rates)
{
	y = gsl_vector_alloc(m);
	for(int k = 0; k < m; k++)
	{
		gsl_vector_set(y, k, rates[k]);
	}
	return;
}

void twoDimInterp::calcB()
{
	gsl_matrix * C3xM = gsl_matrix_calloc(3, m);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Nmat, Minv, 0.0, C3xM);
	gsl_matrix * C3x3 = gsl_matrix_calloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C3xM, Nmat, 0.0, C3x3);
	int s;
	gsl_permutation * p = gsl_permutation_alloc(3);
	gsl_linalg_LU_decomp(C3x3, p, &s);
	gsl_matrix * C3x3inv = gsl_matrix_alloc(3, 3);
	gsl_linalg_LU_invert(C3x3, p, C3x3inv);
	gsl_matrix_free(C3x3);
	gsl_permutation_free(p);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, C3x3inv, Nmat, 0.0, C3xM);
	gsl_matrix *C3xM2 = gsl_matrix_calloc(3, m);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C3xM, Minv, 0.0, C3xM2);
	Bcoef = gsl_vector_calloc(3);
	gsl_blas_dgemv(CblasNoTrans, 1.0, C3xM2, y, 0.0, Bcoef);
	gsl_matrix_free(C3xM);
	gsl_matrix_free(C3xM2);
	gsl_matrix_free(C3x3inv);
	b = new double[3];
	for(int k = 0; k < 3; k++)
	{
		b[k] = gsl_vector_get(Bcoef, k);
	}
	return;
}

void twoDimInterp::calcA()
{
	gsl_vector * Atemp = gsl_vector_alloc(m);
	Acoef = gsl_vector_calloc(m);
	gsl_vector_memcpy(Atemp, y);
	gsl_blas_dgemv(CblasNoTrans, -1.0, Nmat, Bcoef, 1.0, Atemp);
	gsl_blas_dgemv(CblasNoTrans, 1.0, Minv, Atemp, 0.0, Acoef);
	a = new double [m];
	for(int k = 0; k < m; k++)
	{
		a[k] = gsl_vector_get(Acoef, k);
	}
	return;
}