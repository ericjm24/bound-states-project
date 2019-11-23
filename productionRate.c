#include "productionRate.h"

struct boost findCMFrame (double p, double mp, double k, double mk, double theta)
{
  double px = p*sin(theta), pz = p*cos(theta);
  double p0 = sqrt(p*p + mp*mp);
  double gammaP = sqrt(1.0 + p*p/(mp*mp)), gammaK = sqrt(1.0 + k*k/(mk*mk));
  double vx = px / (mp*gammaP + mk*gammaK), vz = (pz + k)/(mp*gammaP + mk*gammaK);
  double v = sqrt(vx*vx + vz*vz);
  if (v == 0){
    struct boost out;
    out.pRel = p;
    out.multFac = 1;
    return out;
    }
  double nx = vx/v, nz = vz/v, gammaCM = 1/sqrt(1.0 - v*v);
  double pRelx = -gammaCM*v*nx*p0 + (1.0 + (gammaCM - 1.0)*nx*nx)*px + (gammaCM - 1.0)*nx*nz*pz;
  double pRelz = -gammaCM*v*nz*p0 + (1.0 + (gammaCM - 1.0)*nz*nz)*pz + (gammaCM - 1.0)*nx*nz*px;
  double pRel = sqrt(pRelx*pRelx + pRelz*pRelz);
  struct boost out;
  out.pRel = pRel;
  out.multFac = sqrt(1 - vx*vx);
  return out;
}

double phaseFactor(double p, double m, double T) {
  return exp(-sqrt(p*p + m*m)/T)*p*p;
}

double productionCrossSectionCalc(double p, double n, double m, double M_n, double M_m, double bindEnergy, double gx) {
  if(p==0)
    return 0.0;
  double mu = 2*M_n*M_m/(M_n+M_m);
  double gp = gx*n*m;
  double k = sqrt(bindEnergy + p*p/(mu*mu));
  double sigma = 16*gp*gp*k*M_PI*fabs(pow(gp,4)/(16*mu*mu*M_PI)-(256*k*k*M_PI*M_PI)/(3*pow(gp*mu,4))+(1835008*pow(k,4)*pow(M_PI,6)/(15*pow(gp*gp*mu,6))))/(p*exp(4));
  printf("%e  ", sigma);
  return sigma;
}

double productionCrossSection(double p, double n, double M_n, double k, double m, double M_m, double theta, double bindEnergy, double gp) {
  struct boost cmFrame = findCMFrame(p, M_n, k, M_m, theta);
  return cmFrame.multFac*productionCrossSectionCalc(cmFrame.pRel, n, m, M_n, M_m, bindEnergy, gp);
}

double integrandFun(double *x, size_t dim, void *fdata){
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
  double a = 0, b = 0;
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
  gsl_monte_function F = {&integrandFun, 3, &fdat};
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

  res = res / (16*M_PI*M_PI*mp*mp*mk*mk*T*T*gsl_sf_bessel_Kn(2,mp/T)*gsl_sf_bessel_Kn(2,mk/T));
  gsl_monte_vegas_free(s);
  gsl_rng_free(r);
  return res;
}

/*int main()
{
  return 0;
}*/