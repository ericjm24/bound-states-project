#include "massCalc.h"

double iFun(double x){
    return 0.5*x*sqrt(1 + x*x) - 0.5*asinh(x);
}

double hFun(double x){
    return 0.25*(iFun(x) + x*x*x*sqrt(1 + x*x));
}

pFermiC::pFermiC(int n, double r)
{
	N = n;
	R0 = r;
	R = r;
	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_cspline, 126);
  
	double x[126], y[126];
	for(int count = 0; count < 126; count++)
	{
		x[count] = count*R/100.0;
		if (x[count]>R)
			y[count] = 0;
		else
			y[count] = pow(3.0*M_PI*315.0*N/64.0, 1.0/3.0)*pow(1 - x[count]/R, 0.5)/R;
	}
	gsl_spline_init(spline, x, y, 126);
	coef = pow(3.0*M_PI*315.0*N/64.0, 1.0/3.0);
}

pFermiC::~pFermiC()
{
}

double pFermiC::eval(double x)
{
	if(bUseConst == true)
	{
		return pow(9.0*M_PI*N/4.0, 1.0/3.0)/R;
	}
	if (x<0)
		x = fabs(x);
	if (x > R)
		return 0.0;
	
	return coef*sqrt(1.0 - (x/R))/R;//(R0/R)*gsl_spline_eval(spline, x*R0/R, acc);
}

void pFermiC::setR(double rad)
{
	R = rad;
	return;
}

double pFermiC::getR()
{
	return R;
}


double solveHydroFunc (double pox, void *params)
{
	struct odeFData *sigma = (struct odeFData *)params;
	if(pox >(sigma->mx/sigma->gp)*(1+1e-10))
	{
		return solveHydroFunc((sigma->mx/sigma->gp)*(1-1e-10), sigma)*(pox/(sigma->mx/sigma->gp));
	}
	else if(pox < -1e-10)
	{
		return solveHydroFunc(0.0, sigma)*(1.0 - (pox/(sigma->mx/sigma->gp)));
	}
	gsl_odeiv2_driver *d = sigma->d;
	double t = 1e-6 * sigma->R;
	double	t1 = sigma->R;
	double y[2] = {pox, 0};
	gsl_odeiv2_driver_reset(d);
	double out;
	
	if (sigma->bIntMode){
		
		sigma->outVec.push_back(0.0);
		sigma->outVec.clear();
		sigma->rVec.push_back(0.0);
		sigma->rVec.clear();
		sigma->outVec.push_back(pox);
		sigma->rVec.push_back(0.0);
		
		for(double k = sigma->R/50.0; k <= sigma->R*(1.0 + 1e-8); k += sigma->R/50.0)
		{
			gsl_odeiv2_driver_apply(d, &t, k, y);
			sigma->rVec.push_back(k);
			sigma->outVec.push_back(y[0]);
		}
		out = y[0];
		
		}
	else{
		gsl_odeiv2_driver_apply(d, &t, t1, y);
		out = y[1] + y[0]/(sigma->R);
		//printf("out = %.5e\n", out);
		}
	return out;
}

double dSolveHydroFunc(double pox, void *params)
{
	struct odeFData *sigma = (struct odeFData *)params;
	gsl_function F;
	F.function = &solveHydroFunc;
	F.params = sigma;
	double result, abserr;
	gsl_deriv_central(&F, pox, 1e-2, &result, &abserr);
	return result;
}

double dSolveE(double R, void *params){
	if (R < 1e-15)
		R = 1e-15;
	gsl_function F;
	double result,abserr;
	F.function = &solveE;
	F.params = params;
	gsl_deriv_central(&F,R,1e-3,&result,&abserr);
	//printf("dSolveE = %f \n", result);
	return result;
}

double d2SolveE(double R, void *params){
	if (R < 1e-15)
		R = 1e-15;
	gsl_function F;
	double result,abserr;
	F.function = &dSolveE;
	F.params = params;
	gsl_deriv_central(&F,R,1e-3,&result,&abserr);
	//printf("d2SolveE = %f \n", result);
	return result;
}

void dd2SolveE(double R, void *params, double *y, double *dy){
	if (R<1e-15)
		R = 1e-15;
	*y = dSolveE(R,params);
	*dy = d2SolveE(R, params);
}

void hydroFDF(double pox, void *params, double *y, double *dy)
{
	*y = solveHydroFunc(pox,params);
	*dy = dSolveHydroFunc(pox, params);
}

int hydroFunc (double t, const double y[], double f[], void *params)
{
	struct odeFData *sigma = (struct odeFData *)params;
	f[0] = y[1];
	double m = sigma->mx - sigma->gp*y[0];
	if (t < 1e-15)
		t = 1e-15;
	
	f[1] = -(sigma->gp/(M_PI*M_PI))*m*m*m * iFun((sigma->pF.eval(t))/fabs(m)) - 2*y[1]/t;
	return GSL_SUCCESS;
}

double hydroIntegralF(double r, void *params)
{
	struct odeFData *sigma = (struct odeFData *)params;
	double gp = sigma->gp;
	double mx = sigma->mx;
	double phi = gsl_spline_eval(sigma->spline2,r,sigma->acc2);
	double ms = mx - gp*phi;
	return r*r*ms*ms*ms*(ms*hFun((sigma->pF.eval(r))/fabs(ms)) + 0.5*gp*phi*iFun((sigma->pF.eval(r))/fabs(ms)));
}

double solveE(double R, void *params)
{
	if(R<=0)
		R = 1e-15;
	
	struct odeFData *sigma = (struct odeFData *)params;
	int status;
	int iter = 0, max_iter = 3000;
	double pox = sigma->pox;
	double pox0;
	sigma->R = R;
	sigma->pF.setR(R);
	sigma->bIntMode = false;
	//gsl_error_handler_t * errT = gsl_set_error_handler_off();
	gsl_root_fsolver *sF = sigma->sF;
	gsl_function FF;
	FF.function = &solveHydroFunc;
	FF.params = sigma;
	gsl_root_fdfsolver *s = sigma->s;
	gsl_function_fdf FDF;
	FDF.f = &solveHydroFunc;
	FDF.df = &dSolveHydroFunc;
	FDF.fdf = &hydroFDF;
	FDF.params = sigma;
	double maxTemp = (sigma->mx/sigma->gp);
	status = gsl_root_fdfsolver_set(s, &FDF, pox);
	if(status == GSL_EINVAL)
	{
		printf("Pox status = %s at R = %.5e\n", gsl_strerror(status), R);
	}
	double x_lo, x_hi;
	//printf("status = %s\n", gsl_strerror(status));
	/*while(status == GSL_EINVAL && fabs(maxTemp) > 1e-50)
	{
		minTemp *= 0.8;
		status = gsl_root_fsolver_set(sF, &FF, minTemp, maxTemp);
		//printf("max = %.5e, %.5e = %s\n", maxTemp, solveHydroFunc(maxTemp, sigma), gsl_strerror(status));
	}*/
	do
	{
		iter++;
		pox0 = pox;
		status = gsl_root_fdfsolver_iterate(s);
		if(status != GSL_CONTINUE && status != GSL_SUCCESS)
		{
			
			printf("%s\n", gsl_strerror(status));
			
			
			double errout = -1.0;
			return errout;
		}
		pox = gsl_root_fdfsolver_root(s);
		//x_lo = gsl_root_fsolver_x_lower(sF);
		//x_hi = gsl_root_fsolver_x_upper(sF);
		status = gsl_root_test_delta(pox, pox0, 1e-6, 1e-6);
		
		if(status != GSL_CONTINUE && status != GSL_SUCCESS)
		{
			
			printf("\n\n%d\n\n", status);
			
			double errout = -1.0;
			return errout;
		}
		if(pox > (sigma->mx/sigma->gp))
		{
			pox = (sigma->mx/sigma->gp)*(1-1e-10);
		}
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	//printf("\n\nHERE\n\n");
	//gsl_set_error_handler(errT);
	//if(pox > (sigma->mx/sigma->gp))
	//	pox = (sigma->mx/sigma->gp);
	
	sigma->pox = pox;
	
	sigma->bIntMode = true;
	
	solveHydroFunc(pox, sigma);
	if(sigma->rVec.size() <= 0)
		return 1e15;
	
	sigma->bIntMode = false;
	sigma->R = R;
	
	sigma->spline2 = gsl_spline_alloc(gsl_interp_cspline, sigma->rVec.size());
	gsl_spline_init(sigma->spline2, sigma->rVec.data(), sigma->outVec.data(), sigma->rVec.size());
	
	double result, error;
	gsl_function F;
	F.function = &hydroIntegralF;
	F.params = sigma;
	
	gsl_integration_qag(&F, 0, R, 0, 1e-4, 10000, GSL_INTEG_GAUSS51, sigma->w, &result, &error);
	
	sigma->bIntMode = false;
	printf("At R = %.3e, Pox = %.4e, E = %.5e\n", R, pox, fabs(result*4.0/M_PI));
	return fabs(result*4.0/M_PI);
}

double solveEwithError(double R, void * params)
{
	struct odeFData *sigma = (struct odeFData *)params;
	double lim = 1.2*solveE(1000000.0, sigma);
	printf("R = %.2e, E = %.8e\n", R, lim + fabs(solveE(R, sigma) - lim));
	return lim + fabs(solveE(R, sigma) - lim);
}

struct outDat massCalc(double n, double mx, double gp)
{
	pFermiC pFermi = pFermiC(n, 0.2);
	double pox = -0.1;
	struct odeFData sigma;
	sigma.mx = mx;
	sigma.gp = gp;
	sigma.R = pFermi.getR();
	sigma.pF = pFermi;
	sigma.pox = pox;
	sigma.bIntMode = false;
	gsl_odeiv2_system sys;
	sys.function = &hydroFunc;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = &sigma;
	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,1e-10,0,1e-6);
	sigma.sys = sys;
	sigma.d = d;
	sigma.s = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
	sigma.sF = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	sigma.w = gsl_integration_workspace_alloc(10000);
	sigma.acc2 = gsl_interp_accel_alloc();
	gsl_error_handler_t * errT = gsl_set_error_handler_off();
	gsl_min_fminimizer *sFin = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
	gsl_function FDF;
	FDF.function = &solveE;
	FDF.params = &sigma;
	double R = sigma.R;
	double Rlow, Rhi;
	double alpha = 1.0;
	int status;
	status = gsl_min_fminimizer_set(sFin, &FDF, alpha*0.5, 1e-10, alpha);
	printf("status = %d\n", status);
	int tempInd = 0;
	
	while(status == GSL_EINVAL && tempInd < 20)
	{
		alpha *= 10.0;
		tempInd++;
		status = gsl_min_fminimizer_set(sFin, &FDF, alpha*0.5, 1e-10, alpha);
	}
	if(status == GSL_FAILURE || tempInd >=5)
	{
		struct outDat errout;
		errout.energy = -1.0;
		return errout;
		printf("\n\nHERE\n\n");
		FDF.function = &solveEwithError;
		gsl_set_error_handler(errT);
		status = gsl_min_fminimizer_set(sFin, &FDF, 318.9660, 1e-10, 10000);
		gsl_set_error_handler_off();
	}
	
	int iter = 0, max_iter = 100000;
	do
	{
		iter++;
		
		status = gsl_min_fminimizer_iterate(sFin);
		if(status != GSL_CONTINUE && status != GSL_SUCCESS)
		{
			printf("\n\n%d\n\n", status);
			
			struct outDat errout;
			errout.energy = -1.0;
			return errout;
		}
		
		R = gsl_min_fminimizer_x_minimum(sFin);
		Rlow = gsl_min_fminimizer_x_lower(sFin);
		Rhi = gsl_min_fminimizer_x_upper(sFin);
		status = gsl_min_test_interval(Rlow, Rhi, 1e-6, 1e-6);
		
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	gsl_set_error_handler(errT);
	sigma.R = R;
	double Energy = solveE(R, &sigma);
	struct outDat output;
	
	output.mass = Energy;
	if(Energy < 0.0)
		output.energy = Energy;
	else
		output.energy = mx*n - Energy;
	
	output.phi0 = -1.0*fabs(sigma.pox);
	
	output.radius = R;
	gsl_root_fdfsolver_free(sigma.s);
	gsl_odeiv2_driver_free(sigma.d);
	gsl_min_fminimizer_free(sFin);
	gsl_integration_workspace_free(sigma.w);
	gsl_interp_accel_free(sigma.acc2);
	gsl_spline_free(sigma.spline2);
	
	return output;
}

struct outDat massCalc(double n, double mx, double gp, bool bUseConst)
{
	pFermiC pFermi = pFermiC(n, .2);
	struct odeFData sigma;
	sigma.mx = mx;
	sigma.gp = gp;
	sigma.R = pFermi.getR();
	sigma.pF = pFermi;
	sigma.pF.setConst(bUseConst);
	sigma.pox = (mx/gp)*(1-1e-8);
	sigma.bIntMode = false;
	gsl_odeiv2_system sys;
	sys.function = &hydroFunc;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = &sigma;
	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,1e-2,1e-6,1e-6);
	sigma.sys = sys;
	sigma.d = d;
	sigma.s = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
	sigma.sF = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	sigma.w = gsl_integration_workspace_alloc(10000);
	sigma.acc2 = gsl_interp_accel_alloc();
	gsl_error_handler_t * errT = gsl_set_error_handler_off();
	gsl_root_fdfsolver *sFin = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_secant);
	gsl_function_fdf FDF;
	FDF.f = &dSolveE;
	FDF.df = &d2SolveE;
	FDF.fdf = &dd2SolveE;
	FDF.params = &sigma;
	double R = sigma.R;
	double Rold;
	double alpha = 1e4;
	int status;
	status = gsl_root_fdfsolver_set(sFin, &FDF, R);
	//printf("status = %s\n", gsl_strerror(status));
	int tempInd = 0;
	int iter = 0, max_iter = 10000;
	do
	{
		iter++;
		Rold = R;
		status = gsl_root_fdfsolver_iterate(sFin);
		if(status != GSL_CONTINUE && status != GSL_SUCCESS)
		{
			//printf("\n\n%d\n\n", status);
			
			struct outDat errout;
			errout.energy = -1.0;
			return errout;
		}
		
		R = gsl_root_fdfsolver_root(sFin);
		status = gsl_root_test_delta(R, Rold, 1e-6, 1e-6);
		
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	//gsl_set_error_handler(errT);
	sigma.R = R;
	double Energy = solveE(R, &sigma);
	struct outDat output;
	
	output.mass = Energy;
	if(Energy < 0.0)
		output.energy = Energy;
	else
		output.energy = mx*n - Energy;
	
	output.phi0 = -1.0*fabs(sigma.pox);
	
	output.radius = R;
	gsl_root_fdfsolver_free(sigma.s);
	gsl_odeiv2_driver_free(sigma.d);
	gsl_root_fdfsolver_free(sFin);
	gsl_integration_workspace_free(sigma.w);
	gsl_interp_accel_free(sigma.acc2);
	gsl_spline_free(sigma.spline2);
	//abort();
	return output;
}

/*int main(int argc, char *argv[]){
	int n = atoi(argv[1]);
	double gp = 4.0;
	double mx = 100.0;
	struct outDat output = massCalc(n, mx, gp);
	printf("For n = %i, mass = %e, binding energy = %e, R = %e, and Phi(0) = %e\n", n, output.mass, output.energy,output.radius, output.phi0);
	return 0;
}*/

int testFullIntegralFunc(double r, const double inF[], double df[], void * params)
{
	//f[0] = x
	//f[1] = u = dx/dr
	//f[2] = the function y(x(r))
	//f[3] = the mass function m(r)
	//f[4] = cumulative number of particles N(r)
	//f[5] = cumulative energy of the system E(r)
	//df[0] = dx/dr
	//df[1] = du/dr = d2x/dr2
	//df[2] = dy/dr
	//df[3] = dm/dr
	//df[4] = dN/dr = 4*Pi*r^2*n(r)
	//df[5] = dE/dr
	//printf("%.5e, %.5e, %.5e, %.5e, %.5e\n", f[0], f[1], f[2], f[3], f[4]);
	double x = inF[0];
	//printf("\nx = %.5e, ", x);
	if(x < 1e-100)
	{
		return GSL_EBADFUNC+10;
	}
	double u = inF[1];
	double y = inF[2];
	
	if(y > -1e-100)
	{
		return GSL_EBADFUNC+10;
	}
	struct fullIntegrandStruct inp = *(struct fullIntegrandStruct *)params;
	double g = inp.gp;
	double mx = inp.m;
	double yp = inp.yp;
	double dx, du, dy, dN, dE;
	
	double x2 = x*x;
	double x4 = x2*x2;
	double sqX = sqrt(x2+1.0);
	double i = 0.5*(x * sqX - asinh(x));
	double f = x4/(3.0*i*sqX);
	
	dx = u;
	
	dy = f*u/(x*f + 1.0);
	//printf("%.5e\n", 1.0 - r*dy);
	
	double mass = mx*exp(y)/(yp);
	//printf("mass = %.5e, mx - mass = %.5e\n", mass, mx-mass);
	double f_prime = x*sqX*(4.0 + x2);
	f_prime -= (4.0 + 3.0*x2)*asinh(x);
	f_prime *= 2.0*x*x2/(3.0*(1.0+x2)*sqX*4.0*i*i);
	
	du = -g*g*mass*mass*i*((1.0/f)+x)/(M_PI*M_PI);
	du -= f_prime*u*u/(f*(1.0 + f*x));
	du += 2.0*u*u*f/(1.0 + f*x);
	du -= 2.0*u/r;
	//printf("%.5e, %.5e, %.5e, %.5e\n",-g*g*mass*mass*i*(1.0+f*x)/(M_PI*M_PI*f), f_prime*u*u/(f*(1.0 + f*x)), 2.0*u*u*f/(1.0 + f*x), -2.0*u/r);
	
	
	dN = 4.0*r*r*x2*x*mass*mass*mass/(3.0*M_PI);
	
	dE = 0.25*mass*(i + x2*x*sqX);
	dE += 0.5*(mx-mass)*i;
	dE *= r*r*mass*mass*mass*4.0/M_PI;
	//printf("%.5e, %.5e, %.5e\n", 0.25*mass*(i + x2*x*sqX), 0.5*(mx-mass)*i, r*r*mass*mass*mass*4.0/M_PI);
	df[0] = dx;
	df[1] = du;
	df[2] = -dy;
	df[3] = dN;
	df[4] = dE;
	//printf("     %.5e, %.5e, %.5e, %.5e, %.5e, %.5e\n", df[0], df[1], df[2], df[3], df[4], mass);
	//printf("u*x = %.5e, x = %.5e\n", u*x, x);
	return GSL_SUCCESS;
}

double yIntFunc(double x, void * params)
{
	(void)(params);
	double x2 = x*x;
	double x4 = x2*x2;
	if(x < 0.01)
	{
		return x - 1.2*x*x2 + (272.0*x*x4/175.0);
	}
	else
	{
		return 2.0*x4/(3.0*x + 3.0*x*x2 + 2.0*x*x4 - 3.0*sqrt(1.0+x2)*asinh(x));
	}
	return 0.0;
}

double testEndCondition(double x, double u, double m, double C, double r)
{
	double i = 0.5*(x*sqrt(1.0+x*x) - asinh(x));
	double f = x*x*x*x/(3.0*i*sqrt(1.0+x*x));
	//printf("%.5e, %.5e, %.5e\n", m*u, C*x/(r*r*(x*f + 1.0)),(1.0 - fabs(C*x/(r*r*m*u*(x*f + 1.0)))));
	return (1.0 - fabs(C*x/(r*r*m*u*(x*f + 1.0))));
}

double testFullIntegralRun(double x0, double ypIn, double g, double mx)
{
	clock_t clockS = clock();
	/*if(x0 < 0.0)
	{
		x0 = fabs(x0);
	}*/
	size_t dim = 5;
	double eps_abs = 0.0;//1e-8;
	double eps_rel = 1e-10;
	gsl_odeiv2_step * step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, dim);
	gsl_odeiv2_control * con = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
	gsl_odeiv2_evolve * evolve = gsl_odeiv2_evolve_alloc(dim);
	struct fullIntegrandStruct params;
	
	double r_init = 1e-10;
	
	params.gp = g;
	params.m = mx;
	params.yp = ypIn;
	gsl_odeiv2_system sys;
	sys.function = &testFullIntegralFunc;
	sys.params = &params;
	sys.dimension = dim;
	
	double h = 1e-15;
	double r = r_init;
	double temp;
	double rMax = 1e60;
	//double stopTime = x0*eps_rel;
	/*
	if(stopTime > 1e-10)
	{
		stopTime = 1e-10;
	}*/
	
	double y0;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  
	double error;

	gsl_function F;
	F.function = &yIntFunc;

	gsl_integration_qag (&F, 0, x0, 0, 1e-10, 1000, GSL_INTEG_GAUSS61, w, &y0, &error); 

	gsl_integration_workspace_free (w);
	y0 = -1.0*y0;
	
	double f[5] = {x0, 0.0, y0, 0.0, 0.0};
	int count = 0;
	int status;
	while(r < rMax && /*f[0]>0.0 &&*/ f[2]<0.0)
	{
		status = gsl_odeiv2_evolve_apply(evolve, con, step, &sys, &r, rMax, &h, f);
		//printf("At r = %.5e, x = %.5e, x' = %.5e, y = %.5e, E = %.5e, N = %.5e\n", r, f[0], f[1], f[2], f[4], f[3]);
		count++;
		if(count > 5000 || status != GSL_SUCCESS)
		{
			break;
		}
	}
	double x = f[0];
	double u = f[1];
	double x2 = x*x;
	double x4 = x2*x2;
	double sqX = sqrt(x2+1.0);
	double i = 0.5*(x * sqX - asinh(x));
	double fTemp = x4/(3.0*i*sqX);
	double dy = fTemp*u/(1.0+x*fTemp);
	gsl_odeiv2_evolve_free(evolve);
	gsl_odeiv2_control_free(con);
	gsl_odeiv2_step_free(step);
	clock_t clockE = clock();
	//printf("For x0 = %.4e, N = %.5e, E = %.5e, epsilon = %.5e, R = %.5e\n", x0, f[3], f[4], mx - (f[4]/f[3]), r);
	//printf("Ended with x = %.5e and y = %.5e\n", f[0], f[2]);
	//printf("Trial ran in %d iterations and took %.5e seconds\n\n", count, ((double)(clockE - clockS)/CLOCKS_PER_SEC));
	return (1.0 - (r*dy));
}

struct outDat testFullIntegralFinal(double x0, double ypIn, double g, double mx)
{
	struct outDat out;
	out.massVec.push_back(0.0);
	out.massVec.clear();
	out.eVec.push_back(0.0);
	out.eVec.clear();
	out.rVec.push_back(0.0);
	out.rVec.clear();
	out.nVec.push_back(0.0);
	out.nVec.clear();
	out.pFermiVec.push_back(0.0);
	out.pFermiVec.clear();
	out.phiVec.push_back(0.0);
	out.phiVec.clear();
	clock_t clockS = clock();
	/*if(x0 < 0.0)
	{
		x0 = fabs(x0);
	}*/
	size_t dim = 5;
	double eps_abs = 0.0;//1e-8;
	double eps_rel = 1e-10;
	gsl_odeiv2_step * step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, dim);
	gsl_odeiv2_control * con = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
	gsl_odeiv2_evolve * evolve = gsl_odeiv2_evolve_alloc(dim);
	struct fullIntegrandStruct params;
	
	double r_init = 1e-10;
	
	params.gp = g;
	params.m = mx;
	params.yp = ypIn;
	gsl_odeiv2_system sys;
	sys.function = &testFullIntegralFunc;
	sys.params = &params;
	sys.dimension = dim;
	
	double h = 1e-15;
	double r = r_init;
	double temp;
	double rMax = 1e60;
	//double stopTime = x0*eps_rel;
	/*
	if(stopTime > 1e-10)
	{
		stopTime = 1e-10;
	}*/
	
	double y0;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  
	double error;

	gsl_function F;
	F.function = &yIntFunc;

	gsl_integration_qag (&F, 0, x0, 0, 1e-12, 1000, GSL_INTEG_GAUSS61, w, &y0, &error); 

	gsl_integration_workspace_free (w);
	y0 = -1.0*y0;
	
	double f[5] = {x0, 0.0, y0, 0.0, 0.0};
	int count = 0;
	int status;
	while(r < rMax && /*f[0]>0.0 &&*/ f[2]<0.0)
	{
		status = gsl_odeiv2_evolve_apply(evolve, con, step, &sys, &r, rMax, &h, f);
		//printf("At r = %.5e, x = %.5e, x' = %.5e, y = %.5e, E = %.5e, N = %.5e\n", r, f[0], f[1], f[2], f[4], f[3]);
		out.rVec.push_back(r);
		out.eVec.push_back(mx*f[3] - f[4]);
		out.massVec.push_back(f[4]);
		out.nVec.push_back(f[3]);
		out.pFermiVec.push_back(f[0]*mx*exp(f[2])/ypIn);
		out.phiVec.push_back(mx*exp(f[2])/ypIn);
		count++;
		if(count > 100000 || status != GSL_SUCCESS)
		{
			break;
		}
	}
	double x = f[0];
	double u = f[1];
	double x2 = x*x;
	double x4 = x2*x2;
	double sqX = sqrt(x2+1.0);
	double i = 0.5*(x * sqX - asinh(x));
	double fTemp = x4/(3.0*i*sqX);
	double dy = fTemp*u/(1.0+x*fTemp);
	gsl_odeiv2_evolve_free(evolve);
	gsl_odeiv2_control_free(con);
	gsl_odeiv2_step_free(step);
	clock_t clockE = clock();
	//printf("For x0 = %.4e, N = %.5e, E = %.5e, epsilon = %.5e, R = %.5e\n", x0, f[3], f[4], mx - (f[4]/f[3]), r);
	//printf("Ended with x = %.5e and y = %.5e\n", f[0], f[2]);
	//printf("%.5e = %.5e\n", ypIn, 1.0-(r*dy));
	//printf("Trial ran in %d iterations and took %.5e seconds\n\n", count, ((double)(clockE - clockS)/CLOCKS_PER_SEC));
	out.phi0 = mx*(1.0 - exp(y0)/(ypIn))/g;
	out.radius = r;
	
	out.N = f[3];
	out.mass = f[4];
	out.energy = mx*f[3] - f[4];
	out.epsilon = mx-(f[4]/f[3]);
	out.rVec.push_back(r);
	out.eVec.push_back(mx*f[3] - f[4]);
	out.massVec.push_back(f[4]);
	out.nVec.push_back(f[3]);
	out.pFermiVec.push_back(f[0]*mx*exp(f[2])/ypIn);
	out.phiVec.push_back(mx*exp(f[2])/ypIn);
	return out;
}

struct outDat completeMassIntegral(double x0, double g, double mx)
{
	double tempOld = 1.0;
	double tempNew = testFullIntegralRun(x0, tempOld, g, mx);
	while(fabs(1.0 - (tempNew/tempOld)) > 1e-8)
	{
		tempOld = tempNew;
		tempNew = testFullIntegralRun(x0, tempOld, g, mx);
	}
	return testFullIntegralFinal(x0, tempNew, g, mx);
}