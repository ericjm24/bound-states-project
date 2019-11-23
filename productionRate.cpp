#include "productionRate.h"

struct myBoost findCMFrame (double p, double mp, double k, double mk, double theta)
{
  double px = p*sin(theta), pz = p*cos(theta);
  double p0 = sqrt(p*p + mp*mp);
  double gammaP = sqrt(1.0 + p*p/(mp*mp)), gammaK = sqrt(1.0 + k*k/(mk*mk));
  double vx = px / (mp*gammaP + mk*gammaK), vz = (pz + k)/(mp*gammaP + mk*gammaK);
  double v = sqrt(vx*vx + vz*vz);
  if (v == 0){
    struct myBoost out;
    out.pRel = p;
    out.multFac = 1;
    return out;
    }
  double nx = vx/v, nz = vz/v, gammaCM = 1/sqrt(1.0 - v*v);
  double pRelx = -gammaCM*v*nx*p0 + (1.0 + (gammaCM - 1.0)*nx*nx)*px + (gammaCM - 1.0)*nx*nz*pz;
  double pRelz = -gammaCM*v*nz*p0 + (1.0 + (gammaCM - 1.0)*nz*nz)*pz + (gammaCM - 1.0)*nx*nz*px;
  double pRel = sqrt(pRelx*pRelx + pRelz*pRelz);
  struct myBoost out;
  out.pRel = pRel;
  out.multFac = sqrt(1.0 - vx*vx);
  return out;
}

double phaseFactor(double p, double m, double T) {
  return exp(-(sqrt(p*p + m*m) - m)/T)*p*p;
}

outDataStruct productionCrossSectionCalc(double p, double n, double m, double M_n, double M_m, double bindEnergy, double gx) {
	if(p==0)
	{
		p = 1e-15;
	}
	double mu = 2.0*M_n*M_m/(M_n+M_m);
	double alpha = gx*gx*n*m/(4.0*M_PI);
	double k = bindEnergy + p*p/mu;
	double bohrR = 2.0/(alpha*mu);
	double sigma;
	double eta = 1.0/bohrR;
	double xi = 1.0/(bohrR*p);
	double temp = p*p+eta*eta;
	double facA, facB, facC, coef;
	gsl_complex numer, denom, numerPow, denomPow;
	gsl_sf_result gamma, arg;
	GSL_SET_COMPLEX(&numer, eta, -1.0*p);
	GSL_SET_COMPLEX(&denom, temp, 0.0);
	GSL_SET_COMPLEX(&numerPow, 0.0, -2.0*xi);
	GSL_SET_COMPLEX(&denomPow, 1.0, -1.0*xi);
	gsl_sf_lngamma_complex_e(1.0, -1.0/xi, &gamma, &arg);
	
	numer = gsl_complex_pow(numer, numerPow);
	denom = gsl_complex_pow(denom, denomPow);
	numer = gsl_complex_div(numer, denom);
	coef = gsl_complex_abs(numer);
	coef = coef * coef * exp(2.0*gamma.val);
	
	facA = alpha*alpha/(mu*mu);
	
	facB = 1.0/(xi*xi);
	facB += 1.0;
	facB *= (k*k*p*p)/(temp*temp*temp);
	facB *= (2.0*alpha)/(3.0*mu*bohrR);
	
	facC = 1.0 + (1.0/(xi*xi));
	facC *= 23.0 + (7.0/(xi*xi));
	facC *= (k*k*k*k*p*p*p*p)/(pow(temp, 6.0));
	facC /= 15.0*bohrR*bohrR;
	
	sigma = 64.0*M_PI*eta*eta*eta * exp(M_PI*xi);
	sigma*= coef;
	sigma *= (facA - facB + facC);
	double sigDecay = sigma;
	
	sigma *= (alpha*k)/2.0;
	
	sigDecay *= (alpha*mu*p)/k;
	
	outDataStruct out;
	out.prodOut = sigma;
	out.decayOut = sigDecay;
  return out;
}

double productionCrossSection(double p, double n, double M_n, double k, double m, double M_m, double theta, double bindEnergy, double gp) {
  struct myBoost cmFrame = findCMFrame(p, M_n, k, M_m, theta);
  outDataStruct out = productionCrossSectionCalc(cmFrame.pRel, n, m, M_n, M_m, bindEnergy, gp);
  out.prodOut *= cmFrame.multFac;
  return out.prodOut;
}

double prodIntegrandFun(double *x, size_t dim, void *fdata){
  (void)(dim);
  struct intData sigma = *((struct intData *) fdata);
  double outVal = 8*M_PI*M_PI*sin(x[2]);
  outVal*=phaseFactor(sigma.a + x[0]/(1 - x[0]),sigma.mp,sigma.T)/((1-x[0])*(1-x[0]));
  outVal*=phaseFactor(sigma.b + x[1]/(1 - x[1]),sigma.mk,sigma.T)/((1-x[1])*(1-x[1]));
  outVal*=productionCrossSection(sigma.a + x[0]/(1 - x[0]),sigma.n,sigma.mp,sigma.b + x[1]/(1 - x[1]),sigma.m,sigma.mk,x[2],sigma.bindEnergy,sigma.gp);
  return outVal;
}

double productionRate (double n, double mp, double m, double mk, double bindEnergy, double gp, double T)
{
	double a = 0.0, b = 0.0;
	struct intData fdat;
	fdat.mp = mp;
	fdat.mk = mk;
	fdat.n = n;
	fdat.m = m;
	fdat.bindEnergy = bindEnergy;
	fdat.gp = gp;
	fdat.T = T;
	fdat.a = a;
	fdat.b = b;
	fdat.M = 0.0;
	fdat.pRel = 0.0;
	
	double result, error;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	gsl_function F;
	F.function = &testFun2;
	F.params = &fdat;
	
	//gsl_error_handler_t * errT = gsl_set_error_handler_off();
	
	int status;
	double reltol = 1e-4;
	double abstol = 1e-100;
	status = gsl_integration_qag (&F, 0.0, 1.0, abstol, reltol, 10000, GSL_INTEG_GAUSS51, w, &result, &error);
	
	/*if(status == GSL_EROUND)
	{
		while(status == GSL_EROUND && tol <= 1e-3)
		{
			tol *= 10.0;
			status = gsl_integration_qagiu (&F, 1.0, 0, tol, 1000, w, &result, &error);
		}
	}*/
	
	gsl_integration_workspace_free (w);
	//gsl_set_error_handler(errT);
	return result*1e-15;
	
	/*
	double a = 0.0, b = 0.0;
	struct intData fdat;
	fdat.mp = mp;
	fdat.mk = mk;
	fdat.n = n;
	fdat.m = m;
	fdat.bindEnergy = bindEnergy;
	fdat.gp = gp;
	fdat.T = T;
	fdat.a = a;
	fdat.b = b;
	fdat.M = 0.0;
	fdat.pRel = 0.0;
	gsl_monte_function F = {&prodIntegrandFun, 3, &fdat};
	double res, err;
	double xl[3] = {0, 0, 0};
	double xu[3] = {1, 1, M_PI};
	const gsl_rng_type *P;
	gsl_rng *r;
	size_t calls = 100000;
	P = gsl_rng_default;
	r = gsl_rng_alloc(P);
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(3);
	gsl_monte_vegas_integrate(&F,xl,xu,3,10000,r,s,&res,&err);

	do
	{
		gsl_monte_vegas_integrate(&F,xl,xu,3,calls,r,s,&res,&err);
	}
	while (fabs (gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

	res = res / (16*M_PI*M_PI*mp*mp*mk*mk*T*T*gsl_sf_bessel_Kn_scaled(2,mp/T)*gsl_sf_bessel_Kn_scaled(2,mk/T));
	gsl_monte_vegas_free(s);
	gsl_rng_free(r);
	return res;*/
}

double decayRate (double n, double mp, double m, double mk, double bindEnergy, double gp, double T)
{
	double M = mp + mk - bindEnergy;
	double M2 = (M*M) - (mp*mp) - (mk*mk);
	double a = 3.0;
	double b = (4.0*(mp*mp + mk*mk)) + 2*M2;
	double c = M2*M2;
	double p = sqrt((sqrt(b*b + 4.0*a*c) - b)/(2.0*a));
	
	double mu = 2.0*mp*mk/(M);
	double alpha = gp*gp*n*m/(4.0*M_PI);
	double k = bindEnergy;
	double bohrR = 2.0/(alpha*mu);
	
	double temp = p*p + (1.0/(bohrR*bohrR));
	double xi = 1.0/(bohrR*p);
	double facA, facB, facC, coef;
	double sigma;
	
	facA = alpha*alpha/(mu*mu);
	
	facB = 1.0/(xi*xi);
	facB += 1.0;
	facB *= (k*k*p*p)/(temp*temp*temp);
	facB *= (2.0*alpha)/(3.0*mu*bohrR);
	
	facC = 1.0 + (1.0/(xi*xi));
	facC *= 23.0 + (7.0/(xi*xi));
	facC *= (k*k*k*k*p*p*p*p)/(pow(temp, 6.0));
	facC /= 15.0*bohrR*bohrR;
	
	sigma = 128.0*M_PI*M_PI/(p*bohrR*bohrR*bohrR*bohrR);
	sigma /= temp*temp;
	sigma *= exp(-4.0*atan(p*bohrR)/(p*bohrR))/(1.0 - exp(-2.0*M_PI/(p*bohrR)));
	sigma *= (facA - facB + facC);
	
	sigma *= alpha*mu/k;
	
	double rateOut = sigma * gsl_sf_bessel_K1_scaled(M/T) / gsl_sf_bessel_Kn_scaled(2, M/T);
	rateOut *= pow(mp*mk*T/(M*2.0*M_PI), 1.5)*exp(-bindEnergy/T);
	//rateOut *= 2.0*T*T*T*exp(-bindEnergy/T)/(M_PI*M_PI);
	return rateOut;
}

double testFun1(double gamma, void * params)
{
	struct intData fdat = *(struct intData *)params;
	double m1 = fdat.mp;
	double m2 = fdat.mk;
	double n1 = fdat.n;
	double n2 = fdat.m;
	double bindEnergy = fdat.bindEnergy;
	double gp = fdat.gp;
	double T = fdat.T;
	
	
	double out = gamma*gamma - 1.0;
	if(out < 1e-300)
	{
		return 0.0;
	}
	else if(gamma > 1e3)
	{
		out = gamma - (0.5/gamma);
	}
	else
	{
		out = sqrt(out);
	}
	
	double E = m1*m1+m2*m2+4.0*m1*m2*gamma*gamma-2.0*m1*m2;
	if(E < 1e-300)
	{
		E = 1e-300;
	}
	else if (gamma*m1*m2 < (m1*m1+m2*m2 - 2.0*m1*m2)*1e-8)
	{
		E = fabs(m1-m2) + 0.5*fabs(4.0*m1*m2*gamma*gamma / (m1-m2));
	}
	else if (gamma*m1*m2 > (m1*m1 + m2*m2 - 2.0*m1*m2)*1e8)
	{
		E = 2.0*gamma*sqrt(m1*m2) + 0.5*fabs((m1-m2)*(m1-m2)/(2.0*sqrt(m1*m2)*gamma));
	}
	else
	{
		E = sqrt(E);
	}
	double p = sqrt(m1*m2)*out;
	double mu = m1*m2/E;
	double alpha = gp*gp*n1*n2/(4.0*M_PI);
	double k = bindEnergy + p*p/mu;
	double bohrR = 2.0/(alpha*mu);
	
	double temp = p*p + (1.0/(bohrR*bohrR));
	double xi = 1.0/(bohrR*p);
	double facA, facB, facC, coef;
	double sigma;
	
	facA = alpha*alpha/(mu*mu);
	
	facB = 1.0/(xi*xi);
	facB += 1.0;
	facB *= (k*k*p*p)/(temp*temp*temp);
	facB *= (2.0*alpha)/(3.0*mu*bohrR);
	
	facC = 1.0 + (1.0/(xi*xi));
	facC *= 23.0 + (7.0/(xi*xi));
	facC *= (k*k*k*k*p*p*p*p)/(pow(temp, 6.0));
	facC /= 15.0*bohrR*bohrR;
	
	sigma = 64.0*M_PI*M_PI/(p*bohrR*bohrR*bohrR*bohrR);
	sigma /= temp*temp;
	sigma *= exp(-4.0*atan(p*bohrR)/(p*bohrR))/(1.0 - exp(-2.0*M_PI/(p*bohrR)));
	sigma *= (facA - facB + facC);
	sigma *= alpha*k;
	
	double x1 = m1/T;
	double x2 = m2/T;
	double X = sqrt(x1*x2);
	double Q = (x1*x1+x2*x2)/(2.0*x1*x2);
	
	out /= sqrt((gamma+Q));
	out *= X * gamma/sqrt(2.0);
	
	double kInp = X*sqrt(2.0*(gamma+Q));
	out *= gsl_sf_bessel_K1_scaled(kInp) / (gsl_sf_bessel_Kn_scaled(2, x1)*gsl_sf_bessel_Kn_scaled(2, x2));
	out *= exp(x1 + x2 - kInp);
	
	out *= sigma;
	
	return out*1e15;
}
	
double testFun2(double v, void * params)
{
	
	double gamma;
	if(v < 1e-300)
	{
		gamma = 1.0;
	}
	else if(v < 1.0)
	{
		gamma = 1.0/sqrt(1.0 - v*v);
	}
	else
	{
		return 0.0;
	}
	struct intData fdat = *(struct intData *)params;
	double m1 = fdat.mp;
	double m2 = fdat.mk;
	double n1 = fdat.n;
	double n2 = fdat.m;
	double bindEnergy = fdat.bindEnergy;
	double gp = fdat.gp;
	double T = fdat.T;
	
	
	double out = gamma*gamma - 1.0;
	if(out < 1e-300)
	{
		return 0.0;
	}
	else if(gamma > 1e3)
	{
		out = gamma - (0.5/gamma);
	}
	else
	{
		out = sqrt(out);
	}
	
	double E = m1*m1+m2*m2+4.0*m1*m2*gamma*gamma-2.0*m1*m2;
	if(E < 1e-300)
	{
		E = 1e-300;
	}
	else if (gamma*m1*m2 < (m1*m1+m2*m2 - 2.0*m1*m2)*1e-8)
	{
		E = fabs(m1-m2) + 0.5*fabs(4.0*m1*m2*gamma*gamma / (m1-m2));
	}
	else if (gamma*m1*m2 > (m1*m1 + m2*m2 - 2.0*m1*m2)*1e8)
	{
		E = 2.0*gamma*sqrt(m1*m2) + 0.5*fabs((m1-m2)*(m1-m2)/(2.0*sqrt(m1*m2)*gamma));
	}
	else
	{
		E = sqrt(E);
	}
	double p = m1*m2*gamma/sqrt(m2*m2 + m1*m2*gamma);
	double mu = m1*m2/E;
	double alpha = gp*gp*n1*n2/(4.0*M_PI);
	double k = bindEnergy + p*p/mu;
	double bohrR = 2.0/(alpha*mu);
	
	double temp = p*p + (1.0/(bohrR*bohrR));
	double xi = 1.0/(bohrR*p);
	double facA, facB, facC, coef;
	double sigma;
	
	facA = alpha*alpha/(mu*mu);
	
	facB = 1.0/(xi*xi);
	facB += 1.0;
	facB *= (k*k*p*p)/(temp*temp*temp);
	facB *= (2.0*alpha)/(3.0*mu*bohrR);
	
	facC = 1.0 + (1.0/(xi*xi));
	facC *= 23.0 + (7.0/(xi*xi));
	facC *= (k*k*k*k*p*p*p*p)/(pow(temp, 6.0));
	facC /= 15.0*bohrR*bohrR;
	
	sigma = 64.0*M_PI*M_PI/(p*bohrR*bohrR*bohrR*bohrR);
	sigma /= temp*temp;
	sigma *= exp(-4.0*atan(p*bohrR)/(p*bohrR))/(1.0 - exp(-2.0*M_PI/(p*bohrR)));
	sigma *= (facA - facB + facC);
	sigma *= alpha*k;
	
	double x1 = m1/T;
	double x2 = m2/T;
	double X = sqrt(x1*x2);
	double Q = (x1*x1+x2*x2)/(2.0*x1*x2);
	
	out = gamma*v*v / sqrt((gamma+Q));
	out *= X /sqrt(2.0);
	
	double kInp = X*sqrt(2.0*(gamma+Q));
	out *= gsl_sf_bessel_K1_scaled(kInp) / (gsl_sf_bessel_Kn_scaled(2, x1)*gsl_sf_bessel_Kn_scaled(2, x2));
	out *= exp(x1 + x2 - kInp);
	
	out *= sigma;
	
	return out*1e15;
}
/*int main()
{
  return 0;
}*/