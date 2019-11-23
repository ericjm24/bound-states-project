#include "ADMmodel.h"

#ifndef PRODUCTIONRATE_H
	#include "productionRate.h"
	#define PRODUCTIONRATE_H
#endif

#ifndef MASSCALC_H
	#include "massCalc.h"
	#define MASSCALC_H
#endif

#ifndef TWO_DIM_INTERP_H
	#include "twoDimInterp.h"
	#define TWO_DIM_INTERP_H
#endif

#ifndef INTEGRATION_H
	#include <gsl/gsl_integration.h>
	#define INTEGRATION_H
#endif

#ifndef MULTIMIN_H
	#include <gsl/gsl_multimin.h>
	#define MULTIMIN_H
#endif

#ifndef MONTE_H
#define MONTE_H
#include <gsl/gsl_monte.h>
#endif

#ifndef MONTEVEGAS_H
#define MONTEVEGAS_H
#include <gsl/gsl_monte_vegas.h>
#endif

#ifndef SPLINE_H
#define SPLINE_H
#include <gsl/gsl_spline.h>
#endif

#ifndef ERRNO_H
#define ERRNO_H
#include <gsl/gsl_errno.h>
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#include <functional>

class own_double_less : public std::binary_function<double,double,bool>
{
public:
  own_double_less( double arg_ = 1e-7 ) : epsilon(arg_) {}
  bool operator()( const double &left, const double &right  ) const
  {
    // you can choose other way to make decision
    // (The original version is: return left < right;) 
    return (fabs(left - right) > epsilon) && (left < right);
  }
  double epsilon;
};
// your map:
//map<double,double,own_double_less> mymap;

double MAX_N = 100000.0;
double NUM_PARTICLES = 1e5;
using namespace std;

#ifndef VECTORS_INIT
	typedef std::vector<double> dvector;
	typedef std::vector<int> intvector;
	#define VECTORS_INIT
#endif

struct paramStruct
{
	double normC;
	double t;
	double N;
	ADMmodel * model;
	double gamma;
};

struct boltzStruct
{
	std::vector<ADMmodel> * models;
	dvector x; // temperatures
	dvector y; // log(numbers)
	dvector times;
	map<double, double, own_double_less> tMap;
	map<double, double, own_double_less> gammaMap;
	map<double, double, own_double_less> normMap;
};

struct derivStruct
{
	double n;
	map<double, double, own_double_less> tMap;
	map<double, double, own_double_less> gammaMap;
	map<double, double, own_double_less> normMap;
};

struct outStruct
{
	double tOut;
	double gammaOut;
	double normOut;
};

dvector::const_iterator findNearest(double val, const dvector& vec)
{
	dvector::const_iterator out;
	dvector::const_iterator it = vec.begin();
	while(val >= *it && it+1!=vec.end())
	{
		it++;
	}
	if (val <= *it && it != vec.begin())
		out = it;
	
	return out;
}

double adaptiveDStep(double (*f)(double, void*), dvector::const_iterator it, dvector::const_iterator itEnd, dvector::const_iterator itBegin, int numD, void * params)
{
	if (numD <= 0)
		return f(*it, params);
	
	dvector::const_iterator next, prev;
	if (*it == *(itEnd-1))
	{
		next = it;
	}
	else
		next = it+1;
	
	if(*it == *itBegin)
		prev = it;
	else
		prev = it-1;
	
	
	return (adaptiveDStep(f, next, itEnd, itBegin, numD-1, params) - adaptiveDStep(f, prev, itEnd, itBegin, numD-1, params))/(*next - *prev);
}

double adaptiveD(double (*f)(double, void*), double x0, const dvector& vec, int numD, void * params)
{
	double tol = 1e-5;
	dvector::const_iterator it = findNearest(x0, vec);
	dvector::const_iterator prev = it-1;
	double result = (f(*it, params) - f(*prev, params))/(*it - *prev);
	double error = 1.0;
	int d = 2;
	double factor = 1;
	double oldRes;
	while (error > tol && d<6)
	{
		oldRes = result;
		factor = factor / d;
		result -= factor*(pow(*it - x0, d) - pow(*prev - x0, d))*(adaptiveDStep(f, prev, vec.end(), vec.begin(), d - 1, params) - adaptiveDStep(f, it, vec.end(), vec.begin(), d-1, params))/(*it - *prev);
		error = fabs((result - oldRes)/result);
		d++;
	}
	return result;
}

double nDerivFunc(double x, void * params)
{
	struct derivStruct * inp = (struct derivStruct *) params;
	if(inp->gammaMap[x] < 0)
	{
		if(inp->n == 1)
			return inp->normMap[x];
		else
			return 0;
		
	}
	return inp->normMap[x]/((1 + ((inp->n - inp->tMap[x])*(inp->n - inp->tMap[x])/(inp->gammaMap[x]*inp->gammaMap[x])))*M_PI*inp->gammaMap[x]);
}

double Temperature(double t)
{
	return (1.33)/sqrt(t);
}

double TimeFromT(double T)
{
	return 1.33*1.33/(T*T);
}

double HubbleC(double t)
{
	return 1e-42;
}

double normFunc(double x, void * params)
{
	struct paramStruct * inp = (struct paramStruct *) params;
	return x/((1 + ((x - inp->t)*(x - inp->t)/(inp->gamma*inp->gamma)))*M_PI*inp->gamma);
}

double getNorm(ADMmodel * inp, paramStruct * params)
{
	gsl_function F;
	F.function = &normFunc;
	F.params = params;
	double result, error;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	gsl_integration_qag(&F, 1.0, MAX_N, 0, 1e-7, 1000, GSL_INTEG_GAUSS51, w, &result, &error);
	gsl_integration_workspace_free(w);
	return result;
}

double nFunc(double x, void * params)
{
	struct paramStruct * inp = (struct paramStruct *) params;
	return inp->normC/((1 + ((x - inp->t)*(x - inp->t)/(inp->gamma*inp->gamma)))*M_PI*inp->gamma);
}

double nGainFunc(double x, void * params)
{
	struct paramStruct *inp = (struct paramStruct *) params;
	return nFunc(x, inp) * nFunc(inp->N - x, inp) * inp->model->getRateFromSpline(x, inp->N - x);
}

double nLossFunc(double x, void * params)
{
	struct paramStruct *inp = (struct paramStruct *) params;
	return nFunc(x, inp)*nFunc(inp->N, inp) * inp->model->getRateFromSpline(inp->N, x);
}

double calcCollisionRate(double n, void * params)
{
	struct paramStruct *inp = (struct paramStruct *) params;
	inp->N = n;
	gsl_function Fg, Fl;
	Fg.function = &nGainFunc;
	Fg.params = inp;
	Fl.function = &nLossFunc;
	Fl.params = inp;
	double out = 0.0;
	double gain_result, loss_result, error;
	int status;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	if (n <= 1.0)
	{
		gain_result = 0.0;
	}
	else
	{
		gain_result = 0.0;
		status = gsl_integration_qag(&Fg, 1.0, n, 1e-4, 3e-3, 1000, GSL_INTEG_GAUSS15, w, &gain_result, &error);
		if(status)
		{
			if(status == GSL_EROUND)
			{
				size_t poopblah;
				gsl_integration_cquad_workspace * w2 = gsl_integration_cquad_workspace_alloc(2000);
				gsl_integration_cquad(&Fg, 1.0, n, 1e-3, 1e-2, w2, &gain_result, &error, &poopblah);
				gsl_integration_cquad_workspace_free(w2);
			}
		}
	}
	gsl_integration_workspace_free(w);
	w = gsl_integration_workspace_alloc (1000);
	loss_result = 0.0;
	status = gsl_integration_qag(&Fl, 1.0, MAX_N, 1e-4, 3e-3, 1000, GSL_INTEG_GAUSS15, w, &loss_result, &error);
	if(status)
	{
		if(status == GSL_EROUND)
		{
			size_t poopblah;
			gsl_integration_cquad_workspace * w2 = gsl_integration_cquad_workspace_alloc(2000);
			gsl_integration_cquad(&Fl, 1.0, MAX_N, 1e-3, 1e-2, w2, &loss_result, &error, &poopblah);
			gsl_integration_cquad_workspace_free(w2);
		}
	}
	gsl_integration_workspace_free(w);
	out -= gain_result;
	out += loss_result;
	if(out*out <= 0)
	{
		out = 1e-20;
	}
	return out;
}

struct interpStruct
{
	gsl_interp_accel * acc;
	gsl_spline * spline;
};

double boltzmannErrorIntFunc(double k, void * params)
{
	interpStruct * inp = (interpStruct *)params;
	return exp(gsl_spline_eval(inp->spline, log(k), inp->acc));
}

double boltzmannErrorTotal(const gsl_vector * v, void * inp)
{
	struct boltzStruct inParams = *(struct boltzStruct *) inp;
	std::vector<ADMmodel> * models = inParams.models;
	struct paramStruct params;
	params.t = 1.0;
	params.gamma = 1.0;
	params.normC = 1.0;
	params.N = 1.0;
	double errOut, res;
	dvector errorVec;
	double mint = *(min_element(inParams.times.begin(), inParams.times.end()));
	double maxt = *(max_element(inParams.times.begin(), inParams.times.end()));
	double errTerm = 0.0;
	map<double, double, own_double_less> tMap = inParams.tMap;
	map<double, double, own_double_less> gammaMap = inParams.gammaMap;
	map<double, double, own_double_less> normMap = inParams.normMap;
	derivStruct derivP;
	derivP.n = 1.0;
	
	int j = models->size()-1;
	if(gsl_vector_get(v, 0) < tMap[inParams.times[j-1]] || gsl_vector_get(v, 1) < 1.0 || gsl_vector_get(v, 0) > MAX_N)
	{
		return 1e100;
	}
	params.model = &(*models)[j];
	tMap[inParams.times[j]] = gsl_vector_get(v, 0);
	gammaMap[inParams.times[j]] = gsl_vector_get(v, 1);
	params.t = tMap[inParams.times[j]];
	params.gamma = gammaMap[inParams.times[j]];
	normMap[inParams.times[j]] = NUM_PARTICLES*exp(-3*HubbleC(inParams.times[j]))/getNorm(params.model, &params);
	params.normC = normMap[inParams.times[j]];
	derivP.tMap = tMap;
	derivP.gammaMap = gammaMap;
	derivP.normMap = normMap;
	for(int k = 0; k < params.model->num.size(); k++)
	{
		if(params.model->num[k] <= MAX_N)
		{
			derivP.n = params.model->num[k];
			errTerm = 0.0;
			errTerm += calcCollisionRate(params.model->num[k], &params);
			errTerm += 3.0*HubbleC(params.model->T)*nFunc(params.model->num[k], &params);
			errTerm += adaptiveD(&nDerivFunc, inParams.times[j], inParams.times, 1, &derivP);
			errorVec.push_back(log(errTerm*errTerm));
		}
	}
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, errorVec.size());
	
	double xin[errorVec.size()];
	double yin[errorVec.size()];
	for(int q = 0; q < errorVec.size(); q++)
	{
		xin[q] = inParams.y[q];
		yin[q] = errorVec[q];
	}
    gsl_spline_init (spline, xin, yin, errorVec.size());
	
	interpStruct intPar;
	intPar.acc = acc;
	intPar.spline = spline;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &boltzmannErrorIntFunc;
	F.params = &intPar;
	
	gsl_integration_qag (&F, 1, MAX_N, 0, 1e-5, 1000, GSL_INTEG_GAUSS51, w, &res, &errOut);
	
	gsl_integration_workspace_free (w);
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	
	/*
	twoDimInterp interp(inParams.x, inParams.y, errorVec, 0.0);
	double xl[2] = {mint, 1.0};
	double xu[2] = {maxt, MAX_N};
	
	gsl_monte_function G = { &boltzmannErrorIntFunc, 2, &interp };
	const gsl_rng_type *Typ;
	gsl_rng *r;
	gsl_rng_env_setup ();
	Typ = gsl_rng_default;
	r = gsl_rng_alloc (Typ);
	
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
	
	printf("\n\nStarting Integration\n\n");
    gsl_monte_vegas_integrate(&G, xl, xu, 2, 10000, r, s, &res, &errOut);
	
    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, 2, 10000, r, s, &res, &errOut);
        printf ("result = %e, sigma = % .6f, chisq/dof = %.1f\n", res, errOut, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);

	gsl_rng_free (r);*/
	
	return res;
}

outStruct solveBoltzmann(std::vector<ADMmodel> models, map<double, double, own_double_less> tMap, map<double, double, own_double_less> gammaMap, map<double, double, own_double_less> normMap)
{
	/*outStruct out2;
	out2.tOut = 1.0;
	out2.gammaOut = 1.0;
	out2.normOut = 1.0;
	return out2;*/
	dvector temps;
	dvector numLogs;
	dvector times;
	
	for(int qw = 0; qw < models.size(); qw++)
	{
		times.push_back(TimeFromT(models[qw].T));
	}
	for(int k = 0; k < models[0].num.size(); k++)
	{
		if(models[0].num[k] <= MAX_N)
		{
			temps.push_back(log(TimeFromT(models[1].T)));
			numLogs.push_back(log(models[0].num[k]));
		}
	}
	
	
	struct boltzStruct params;
	params.models = &models;
	params.x = temps;
	params.y = numLogs;
	params.times = times;
	params.tMap = tMap;
	params.gammaMap = gammaMap;
	params.normMap = normMap;
	
	size_t ndim = 2;
	string breakstring;
	gsl_multimin_function minex_func;
	minex_func.f = &boltzmannErrorTotal;
	minex_func.n = ndim;
	minex_func.params = &params;
	const gsl_multimin_fminimizer_type *minType = gsl_multimin_fminimizer_nmsimplex2;
	gsl_vector * paramVec = gsl_vector_alloc(ndim);
	gsl_vector * ss = gsl_vector_alloc(ndim);
	gsl_vector_set(paramVec, 0, tMap[*(times.end()-2)]+1.0);
	double woru = gammaMap[*(times.end()-2)];
	if(woru < 1)
		woru = 1;
	
	gsl_vector_set(paramVec, 1, woru);
	gsl_vector_set(ss, 0, 0.5*(tMap[*(times.end()-2)]+1.0));
	gsl_vector_set(ss, 1, 0.125*woru);
	
	
	size_t iter = 0;
	int status;
	double size;
	
	gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (minType, ndim);
	gsl_multimin_fminimizer_set (s, &minex_func, paramVec, ss);
	
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status) 
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 5e-3);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}

		printf ("\nRun %3d Completed: error = %7.3e, Nav = %e, gamma = %e, T = %.2e, size = %.3f\n", iter, s->fval, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),*(times.end()-1), size);
	}
	while (status == GSL_CONTINUE && iter < 500);
	outStruct out;
	out.tOut = gsl_vector_get(s->x, 0);
	out.gammaOut = gsl_vector_get(s->x, 1);
	
	struct paramStruct params2;
	params2.t = 1.0;
	params2.gamma = 1.0;
	params2.normC = 1.0;
	params2.N = 1.0;
	params2.model = &models[models.size()-1];
	
	out.normOut = NUM_PARTICLES*exp(-3*HubbleC(*(times.end()-1)))/getNorm(params2.model, &params2);
	return out;
}

double vecGen(int in)
{
	double out;
	if(in < 0 && in%5 != 0)
		in -=5;
	
	switch(in%5)
	{
	case 0:
		out = 1.0;
		break;
	case 1:
	case -4:
		out = 1.5;
		break;
	case 2:
	case -3:
		out = 2.5;
		break;
	case 3:
	case -2:
		out = 4.0;
		break;
	case 4:
	case -1:
		out = 6.5;
		break;
	}
	return out*pow(10.0,in/5);
}

int main()
{
	double mx = 140.0;
	double gp = 4.0;
	//gsl_set_error_handler_off();
	/*map<string, void (ADMmodel::*)(double)> functions;
	functions["getParticleProps"] = &ADMmodel::getParticleProps;
	functions["makeSplines"] = &ADMmodel::makeSplines;
	functions["getProductionRates"] = &ADMmodel::getProductionRates;
	functions["getMass"] = &ADMmodel::display_getMass;
	functions["getE"] = &ADMmodel::display_getE;
	functions["getR"] = &ADMmodel::display_getR;
	functions["getPhi"] = &ADMmodel::display_getPhi;
	functions["getAllProductionRates"] = &ADMmodel::getAllProductionRates;
	functions["save"] = &ADMmodel::writeBoundProps;
	functions["load"] = &ADMmodel::readBoundProps;
	
	size_t found;
	double arg;
	string func;*/
	double n[21] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 100.0, 200.0, 400.0, 700.0, 1000.0, 2000.0, 4000.0, 7000.0, 10000.0, 20000.0, 40000.0, 70000.0, 100000.0,200000.0};
	MAX_N = 100000.0;
	//double t[15] = {0.1, 0.5, 1.0, 2.0, 4.0, 7.0, 10.0, 20.0, 40.0, 70.0, 100.0, 200.0, 400.0, 700.0, 1000.0};
	double t = vecGen(-5);
	double minT = t;
	double maxT = t;
	ADMmodel adm = ADMmodel(mx, gp, Temperature(t));
	dvector numVec;
	for(int k2 = 0; k2<21; k2++)
	{
		adm.getParticleProps(n[k2]);
		if(n[k2]<=MAX_N)
			numVec.push_back(log(n[k2]));
	}
	//adm.getAllProductionRates(Temperature(t));
	//adm.writeBoundProps(Temperature(t));
	std::vector<ADMmodel> Models;
	Models.push_back(adm);
	
	ADMmodel admTemp = adm;
	map<double, double, own_double_less> tMap;
	map<double, double, own_double_less> gammaMap;
	map<double, double, own_double_less> normMap;
	tMap[t] = 1.0;
	gammaMap[t] = 0;
	
	struct paramStruct params;
	params.t = tMap[t];
	params.gamma = gammaMap[t];
	params.normC = 1.0;
	params.N = 1.0;
	normMap[t] = NUM_PARTICLES*exp(-3*HubbleC(t))/getNorm(params.model, &params);
	outStruct out;
	dvector Navs, Gammas, Norms, times;
	Navs.push_back(tMap[t]);
	Gammas.push_back(gammaMap[t]);
	Norms.push_back(normMap[t]);
	times.push_back(t);
	int q = -4;
	while(vecGen(q) < 1200)
	{
		t = vecGen(q);
		if(t > maxT)
			maxT = t;
		
		admTemp = ADMmodel(mx, gp, Temperature(t));
		Models.push_back(admTemp);
		admTemp.writeBoundProps(admTemp.T);
		out = solveBoltzmann(Models, tMap, gammaMap, normMap);
		tMap[t] = out.tOut;
		gammaMap[t] = out.gammaOut;
		normMap[t] = out.normOut;
		Navs.push_back(out.tOut);
		Gammas.push_back(out.gammaOut);
		Norms.push_back(out.normOut);
		times.push_back(t);
		q++;
	}
	return 0;
	derivStruct derivP;
	derivP.tMap = tMap;
	derivP.gammaMap = gammaMap;
	derivP.normMap = normMap;
	double errTerm = 0;
	double res, errOut;
	dvector errorVec;
	const int BASELINE_NAV = 1;
	const int BASELINE_GAMMA = 2;
	const int NAV_CALC = 3;
	const int NAV_STEP = 4;
	const int GAMMA_CALC = 5;
	const int GAMMA_STEP = 6;
	int CURRENT_STATE = BASELINE_NAV;
	double relErr, value1, value2, base, step, diffVal;
	double errTol = 1e-3;
	int ind = 0;
	bool bContinue = true;
	for(int p = 1; p < Models.size(); p++)
	{
		ind = 0;
		bContinue = true;
		int CURRENT_STATE = BASELINE_NAV;
		while(bContinue && ind < 50)
		{
			printf("Ind = %d, ", ind);
			params.t = Navs[p];
			params.gamma = Gammas[p];
			params.normC = Norms[p];
			params.model = &Models[p];
			for(int k = 0; k < params.model->num.size(); k++)
			{
				if(params.model->num[k] <= MAX_N)
				{
					derivP.n = params.model->num[k];
					errTerm = 0.0;
					errTerm += calcCollisionRate(params.model->num[k], &params);
					errTerm += 3.0*HubbleC(params.model->T)*nFunc(params.model->num[k], &params);
					errTerm += adaptiveD(&nDerivFunc, times[p], times, 1, &derivP);
					errorVec.push_back(log(errTerm*errTerm));
				}
			}
			gsl_interp_accel *acc = gsl_interp_accel_alloc ();
			gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, errorVec.size());
			
			double xin[errorVec.size()];
			double yin[errorVec.size()];
			for(int r = 0; r < errorVec.size(); r++)
			{
				xin[r] = numVec[r];
				yin[r] = errorVec[r];
			}
			gsl_spline_init (spline, xin, yin, errorVec.size());
			
			interpStruct intPar;
			intPar.acc = acc;
			intPar.spline = spline;
			
			gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
			gsl_function F;
			F.function = &boltzmannErrorIntFunc;
			F.params = &intPar;
			
			gsl_integration_qag (&F, 1, MAX_N, 0, 1e-5, 1000, GSL_INTEG_GAUSS51, w, &res, &errOut);
			
			gsl_integration_workspace_free (w);
			gsl_spline_free(spline);
			gsl_interp_accel_free(acc);
			errorVec.clear();
			bContinue = true;
			switch(CURRENT_STATE)
			{
				case BASELINE_NAV:
					base = 0.05*Navs[p];
					printf("Base = %e\n", base);
					value1 = res;
					CURRENT_STATE = NAV_STEP;
					step = 0.05*Navs[p];
					Navs[p] = Navs[p] + step;
					break;
				case BASELINE_GAMMA:
					base = res;
					value1 = res;
					CURRENT_STATE = GAMMA_STEP;
					step = 0.05*Gammas[p];
					Gammas[p] = Gammas[p] + step;
					break;
				case NAV_CALC:
					if(res > value2)
					{
						printf("Val1 = %e, N_av = %e\n", value1, Navs[p]);
						Navs[p] = Navs[p] - 0.5*step;
						step = step*0.5;
						break;
					}
					if(fabs(res/base) < errTol)
						bContinue = false;
					else
					{
						printf("Val1 = %e, N_av = %e\n", value1, Navs[p]);
						value1 = res;
						CURRENT_STATE = NAV_STEP;
						step = 0.05*Navs[p];
						Navs[p] = Navs[p] + step;
					}
					break;
				case NAV_STEP:
					value2 = res;
					diffVal = (value2 - value1)/step;
					printf("Val2 = %e, Diff = %e\n", value2, diffVal);
					if((value1/diffVal) > (Navs[p] - Navs[p-1]))
					{
						Navs[p] = 0.5*Navs[p];
						step = -0.5*Navs[p];
					}
					else
					{
						Navs[p] = Navs[p] - step - 0.5*(value1/diffVal);
						step = -0.5*(value1/diffVal);
					}
					CURRENT_STATE = NAV_CALC;
					break;
				case GAMMA_CALC:
					if(fabs(res/base) < errTol)
						bContinue = false;
					else
					{
						value1 = res;
						CURRENT_STATE = GAMMA_STEP;
						step = 0.05 * Gammas[p];
						Gammas[p] = Gammas[p]+step;
					}
					break;
				case GAMMA_STEP:
					value2 = res;
					diffVal = (value2 - value1)/step;
					Gammas[p] = Gammas[p] - step - (value1/diffVal);
					CURRENT_STATE = GAMMA_CALC;
					break;
			}
			derivP.tMap[times[p]] = Navs[p];
			derivP.gammaMap[times[p]] = Gammas[p];
			ind++;
		}
		printf("Error for T = %.2f is %e\n", times[p], res);
	}
	
	printf("\n\n\n\n");
	for(int cvb = 0; cvb < Models.size(); cvb++)
	{
		printf("At time T = %e, Nav = %e, Gamma = %e\n", TimeFromT(Models[cvb].T), Navs[cvb], Gammas[cvb]);
	}
	
	/*
	printf("Input: ");
	cin >> func;
	while(func != "exit")
	{
		found = func.find("(");
		if (functions.find(func.substr(0, found)) != functions.end())
		{
			arg = atof(func.substr(found+1, func.find(")")-found).c_str());
			(adm.*functions[func.substr(0,found)])(arg);
		};
		
		printf("Input: ");
		cin >> func;
	};
	*/
	return 0;
}