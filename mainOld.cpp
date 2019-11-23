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

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double MAX_N = 10000.0;
double NUM_PARTICLES = 1.0;
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
};

struct derivStruct
{
	double n;
	map<double, double> tMap;
	map<double, double> gammaMap;
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
	return 1.0/((1 + ((inp->n - inp->tMap[x])*(inp->n - inp->tMap[x])/(inp->gammaMap[x]*inp->gammaMap[x])))*M_PI*inp->gammaMap[x]);
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

double boltzmannErrorIntFunc(double * k, size_t dim, void * params)
{
	(void)(dim);
	twoDimInterp * interp = (twoDimInterp *)params;
	return exp(interp->eval(log(k[0]), log(k[1])));
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
	map<double, double> tMap, gammaMap;
	for(int p = 0; p < models->size(); p++)
	{
		if(p == 0)
		{
			tMap[inParams.times[p]] = 1.0;
			gammaMap[inParams.times[p]] = 1.0;
		}
		else
		{
			tMap[inParams.times[p]] = gsl_vector_get(v, 2*p-2);
			gammaMap[inParams.times[p]] = gsl_vector_get(v, 2*p - 1);
		}
		printf("At t = %.2e, Nav = %.2e, Gamma = %.2e\n", inParams.times[p], tMap[inParams.times[p]], gammaMap[inParams.times[p]]);
	}
	derivStruct derivP;
	derivP.n = 1.0;
	derivP.tMap = tMap;
	derivP.gammaMap = gammaMap;
	
	for(int j = 0; j < models->size(); j++)
	{
		params.model = &(*models)[j];
		if(j == 0)
		{
			params.t = 1.0;
			params.gamma = 1.0;
		}
		else
		{
			params.t = gsl_vector_get(v, 2*j-2);
			params.gamma = gsl_vector_get(v, 2*j - 1);
		}
		params.normC = NUM_PARTICLES*exp(-3*HubbleC(inParams.times[j]))/getNorm(params.model, &params);
		for(int k = 0; k < params.model->num.size(); k++)
		{
			if(params.model->num[k] < MAX_N)
			{
				derivP.n = params.model->num[k];
				errTerm = 0.0;
				errTerm += calcCollisionRate(params.model->num[k], &params);
				errTerm += 3.0*HubbleC(params.model->T)*nFunc(params.model->num[k], &params);
				errTerm += adaptiveD(&nDerivFunc, inParams.times[j], inParams.times, 1, &derivP);
				errorVec.push_back(log(errTerm*errTerm/derivP.n));
			}
		}
	}
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

	gsl_rng_free (r);
	
	printf("Total Squared Error is %e\n", res);
	return res;
}

void solveBoltzmann(std::vector<ADMmodel> models)
{
	dvector temps;
	dvector numLogs;
	dvector times;
	for(int j = 0; j < models.size(); j++)
	{
		times.push_back(TimeFromT(models[j].T));
		for(int k = 0; k < models[0].num.size(); k++)
		{
			if(models[0].num[k] <= MAX_N)
			{
				temps.push_back(log(TimeFromT(models[j].T)));
				numLogs.push_back(log(models[0].num[k]));
			}
		}
	}
	
	struct boltzStruct params;
	params.models = &models;
	params.x = temps;
	params.y = numLogs;
	params.times = times;
	
	size_t ndim = 2*(models.size()-1);
	string breakstring;
	gsl_multimin_function minex_func;
	minex_func.f = &boltzmannErrorTotal;
	minex_func.n = ndim;
	minex_func.params = &params;
	const gsl_multimin_fminimizer_type *minType = gsl_multimin_fminimizer_nmsimplex2;
	gsl_vector * paramVec = gsl_vector_alloc(ndim);
	gsl_vector * ss = gsl_vector_alloc(ndim);
	gsl_vector_set_all(paramVec, 1.0);
	gsl_vector_set_all(ss, 0.1);
	for(int kj = 0; kj < times.size()-1; kj++)
	{
		gsl_vector_set(paramVec, 2*kj, 1.0);
		gsl_vector_set(ss, 2*kj, 0.25);
	}
	
	
	size_t iter = 0;
	int status;
	double size;
	
	gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (minType, ndim);
	gsl_multimin_fminimizer_set (s, &minex_func, paramVec, ss);
	
	do
	{
		printf("\n");
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status) 
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-2);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}

		printf ("\n\nRun %3d Completed: f() = %7.3e, size = %.3f\n\n", iter, s->fval, size);
	}
	while (status == GSL_CONTINUE && iter < 500);
	
	return;
}

int main()
{
	double mx = 100.0;
	double gp = 1.0;
	
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
	double n[17] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 100.0, 200.0, 400.0, 700.0, 1000.0, 2000.0, 4000.0, 7000.0, 10000.0, 20000.0};
	MAX_N = 10000.0;
	double t[15] = {0.1, 0.5, 1.0, 2.0, 4.0, 7.0, 10.0, 20.0, 40.0, 70.0, 100.0, 200.0, 400.0, 700.0, 1000.0};
	
	ADMmodel adm = ADMmodel(mx, gp, Temperature(t[0]));
	for(int k = 0; k<17; k++)
	{
		adm.getParticleProps(n[k]);
	}
	//adm.getAllProductionRates(Temperature(t[0]));
	//adm.writeBoundProps(Temperature(t[0]));
	printf("Collision Rate: %e\n", adm.getRateFromSpline(8.0, 3.0));
	std::vector<ADMmodel> Models;
	Models.push_back(adm);
	
	ADMmodel admTemp = adm;
	for(int q = 1; q < 15; q++)
	{
		admTemp = ADMmodel(mx, gp, Temperature(t[q]));
		Models.push_back(admTemp);
		admTemp.writeBoundProps(admTemp.T);
	}
	
	solveBoltzmann(Models);
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