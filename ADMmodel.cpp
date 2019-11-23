#include "ADMmodel.h"

ADMmodel::ADMmodel(double m, double g, double Temp)
{
	bMadeSplines = false;
	
	lambda = 0.0;
	mx = m;
	gp = g;
	T = Temp;
	x = mx / T;
	char tStr[50];
	sprintf(tStr, "_M%.3e_G%.3e", m, g);
	crosssection_path = "cs" + string(tStr);
	decay_path = "de" + string(tStr);
	boundstate_properties_path = "prop" + string(tStr) + ".txt";
	D = 0;
	maxN = 1000;
	
	massSpline = gsl_spline_alloc(gsl_interp_cspline, 3);
	ESpline = gsl_spline_alloc(gsl_interp_cspline, 3);
	radSpline = gsl_spline_alloc(gsl_interp_cspline, 3);
	phiSpline = gsl_spline_alloc(gsl_interp_cspline, 3);
	
	mAcc = gsl_interp_accel_alloc();
	eAcc = gsl_interp_accel_alloc();
	rAcc = gsl_interp_accel_alloc();
	pAcc = gsl_interp_accel_alloc();
	
	getAllParticleProps(maxN);
	//readBoundProps(T);
}

void ADMmodel::getParticleProps(double n)
{
	if (n < 1)
		return;
	
	if(n <= maxN)
		return;
	
	struct outDat dat = massCalc(n, mx, gp);
	
	num.push_back(n);
	if(n > maxN)
		maxN = n;
		
	masses.push_back(dat.mass);
	bindEnergies.push_back(dat.energy);
	radii.push_back(dat.radius);
	phi0.push_back(dat.phi0);
	D+=1;
	
	printf("The particle of mass-number %i has the following properties:\n\n", int(n+0.5));
	printf("Mass: %e\nBinding Energy: %e\nRadius: %e\nPhi(0): %e\n\n", dat.mass, dat.energy, dat.radius, dat.phi0);
	bMadeSplines = false;
	if(D > 4)
		makeSplines(0.0);
	return;
}

void ADMmodel::getAllParticleProps(double n)
{
	if (n < 1.0)
	{
		return;
	}
	
	maxN = n;
	double x0 = 1.0;
	printf("Calculating particle properties.\nFinding starting value for x0...\n");
	printf("Testing x0 = %.4e...\n", x0);
	struct outDat temp = completeMassIntegral(x0, gp, mx);
	
	if (temp.N > 1.0)
	{
		while(temp.N > 1.0)
		{
			x0 /= 2.0;
			printf("Testing x0 = %.4e...\n", x0);
			temp = completeMassIntegral(x0, gp, mx);
		}
	}
	if (temp.N < 1.0)
	{
		while(true)
		{
			x0 *= 2.0;
			printf("Testing x0 = %.4e...\n", x0);
			temp = completeMassIntegral(x0, gp, mx);
			if(temp.N > 1.0)
			{
				x0 /= 2.0;
				temp = completeMassIntegral(x0, gp, mx);
				break;
			}
		}
	}
	num.push_back(1.0);
	masses.push_back(mx);
	bindEnergies.push_back(0.0);
	radii.push_back(temp.radius);
	phi0.push_back(temp.phi0);
	double multFac = 100.0;
	double x0old = x0;
	double lastN = num[0];
	printf("Calculating from x0 = %.4e...\n", x0old);
	while(num.size() < 20)
	{
		x0 = x0old;
		
		multFac = sqrt(multFac);
		temp.N = num[0];
		temp.mass = masses[0];
		temp.energy = bindEnergies[0];
		temp.radius = radii[0];
		temp.phi0 = phi0[0];
		
		num.clear();
		masses.clear();
		bindEnergies.clear();
		radii.clear();
		phi0.clear();
		
		num.push_back(temp.N);
		masses.push_back(temp.mass);
		bindEnergies.push_back(temp.energy);
		radii.push_back(temp.radius);
		phi0.push_back(temp.phi0);
		lastN = num[0];
		while(temp.N < n)
		{
			x0 *= multFac;
			temp = completeMassIntegral(x0, gp, mx);
			if(temp.N <= lastN || temp.N < 1.0)
			{
				continue;
			}
			num.push_back(temp.N);
			masses.push_back(temp.mass);
			if(temp.energy < 0.0)
			{
				bindEnergies.push_back(0.0);
			}
			else
			{
				bindEnergies.push_back(temp.energy);
			}
			radii.push_back(temp.radius);
			phi0.push_back(temp.phi0);
			lastN = temp.N;
		}
	}
	D = num.size();
	for(int k1 = 0; k1 < D; k1++)
	{
		printf("N = %.2e, Mass = %.5e, Energy = %.5e, Radius = %.5e\n", num[k1], masses[k1], bindEnergies[k1]/num[k1], radii[k1]);
	}
	makeSplines(0.0);
	return;
}

void ADMmodel::makeSplines(double x)
{
	if(D < 2)
		return;
		
	if(!bMadeSplines)
		bMadeSplines = true;
	(void)x;
	
	gsl_spline_free(massSpline);
	gsl_spline_free(ESpline);
	gsl_spline_free(radSpline);
	gsl_spline_free(phiSpline);
	dvector dnum(num.begin(), num.end());
	massSpline = gsl_spline_alloc(gsl_interp_cspline, D);
	ESpline = gsl_spline_alloc(gsl_interp_cspline, D);
	radSpline = gsl_spline_alloc(gsl_interp_cspline, D);
	phiSpline = gsl_spline_alloc(gsl_interp_cspline, D);
	
	gsl_spline_init(massSpline, &dnum[0], &masses[0], D);
	gsl_spline_init(ESpline, &dnum[0], &bindEnergies[0], D);
	gsl_spline_init(radSpline, &dnum[0], &radii[0], D);
	gsl_spline_init(phiSpline, &dnum[0], &phi0[0], D);
	
	return;
}

void ADMmodel::getProductionRates(double n)
{
	if(!bMadeSplines)
		makeSplines(0.0);
		
	dvector out;
	double mas = getMass(double(n));
	double bindEnergy = getE(double(n));
	double mas2, e2, bindE, r;
	for(int k = 0; k < D; k++)
	{
		if (num[k]+n > maxN || num[k] > n)
			continue;
		
		mas2 = getMass(double(num[k]));
		e2 = getE(double(num[k]));
		bindE = fabs(getE(double(n + num[k])) - bindEnergy - e2);
		r = productionRate(n, mas, num[k], mas2, bindE, gp, T);
		out.push_back(r);
		printf("Production rate for %i + %i -> %i : %e\n", int(n+0.5), int(num[k]+0.5), int(n+num[k]+0.5), r);
	};
	prodRates.insert(prodRates.end(), out.begin(), out.end());
	return;/*out;*/
}

void ADMmodel::getDecayRates(double n)
{
	if(!bMadeSplines)
		makeSplines(0.0);
		
	dvector out;
	double mas = getMass(double(n));
	double bindEnergy = getE(double(n));
	double mas2, e2, bindE, r;
	for(int k = 0; k < D; k++)
	{
		if (num[k]+n > maxN || num[k] > n)
			continue;
		
		mas2 = getMass(double(num[k]));
		e2 = getE(double(num[k]));
		bindE = fabs(getE(double(n + num[k])) - bindEnergy - e2);
		r = decayRate(n, mas, num[k], mas2, bindE, gp, T);
		out.push_back(r);
		printf("Decay rate for %i -> %i + %i : %e\n", int(n+num[k]+0.5), int(n+0.5), int(num[k]+0.5), r);
	};
	decayRates.insert(decayRates.end(), out.begin(), out.end());
	return;/*out;*/
}

void ADMmodel::getAllProductionRates(double n)
{
	(void)n;
	dvector numsort = num;
	sort(numsort.begin(), numsort.end(), std::greater<double>());
	if(numsort[0] < numsort[1]+numsort[2])
		getParticleProps(numsort[0]*2);
	D = num.size();
	prodRates.clear();
	sort(num.begin(), num.end());
	for(int k = 0; k < D; k++)
	{
		getProductionRates(num[k]);
	};
	numRange[0] = *std::min_element(num.begin(), num.end());
	numRange[1] = *std::max_element(num.begin(), num.end());
	dvector x;
	dvector y;
	for(int j = 0; j < D; j++)
	{
		for(int l = 0; l<=j; l++)
		{
			x.push_back((double)log(num[j]));
			y.push_back((double)log(num[l]));
		}
	}
	dvector ratesOut;
	for(int jk = 0; jk < prodRates.size(); jk++)
	{
		ratesOut.push_back((double)log(prodRates[jk]));
	}
	rateSpline.reset(x, y, ratesOut, (double)lambda);
	return;
}

void ADMmodel::getAllDecayRates(double n)
{
	(void)n;
	dvector numsort = num;
	sort(numsort.begin(), numsort.end(), std::greater<double>());
	if(numsort[0] < numsort[1]+numsort[2])
		getParticleProps(numsort[0]*2);
	
	sort(num.begin(), num.end());
	D = num.size();
	decayRates.clear();
	for(int k = 0; k < D; k++)
	{
		getDecayRates(num[k]);
	};
	numRange[0] = *std::min_element(num.begin(), num.end());
	numRange[1] = *std::max_element(num.begin(), num.end());
	dvector x;
	dvector y;
	for(int j = 0; j < D; j++)
	{
		for(int l = 0; l<=j; l++)
		{
			x.push_back((double)log(num[j]));
			y.push_back((double)log(num[l]));
		}
	}
	dvector decayOut;
	for(int jk = 0; jk < decayRates.size(); jk++)
	{
		decayOut.push_back((double)log(decayRates[jk]));
	}
	decaySpline.reset(x, y, decayOut, (double)lambda);
	return;
}

void ADMmodel::explicitRateCalculation(int max_num, dvector & production, dvector & decay, double Temper)
{
	double mas, bindEnergy, mas2, e2, bindE, jDoub, kDoub;
	int ind = (max_num*max_num + max_num)/2;
	production.assign(ind, 0.0);
	decay.assign(ind, 0.0);
	ind = 0;
	for(int j = 1; j < max_num; j++)
	{
		jDoub = double(j);
		mas = getMass(jDoub);
		bindEnergy = getE(jDoub);
		if(bindEnergy < 0.0)
		{
			bindEnergy = 0.0;
		}
		for(int k = 1; k <= j; k++)
		{
			kDoub = double(k);
			mas2 = getMass(kDoub);
			e2 = getE(kDoub);
			if(e2 < 0.0)
			{
				e2 = 0.0;
			}
			bindE = fabs(getE(jDoub+kDoub) - bindEnergy - e2);
			production[ind] = productionRate(jDoub, mas, kDoub, mas2, bindE, gp, Temper);
			decay[ind] = decayRate(jDoub, mas, kDoub, mas2, bindE, gp, Temper);
			ind++;
		}
	}
	return;
}

double ADMmodel::singleRateCalculation(int n1, int n2, double Temper)
{
	//printf("Calculating for n1 = %d, n2 = %d, m1 = %.3e, m2 = %.3e, bE = %.3e, T
	return productionRate(double(n1), getMass(double(n1)), double(n2), getMass(double(n2)), fabs(getE(double(n1+n2))-getE(double(n1))-getE(double(n2))), gp, Temper);
}

double ADMmodel::singleDecayCalculation(int n1, int n2, double Temper)
{
	//printf("Calculating for n1 = %d, n2 = %d, m1 = %.3e, m2 = %.3e, bE = %.3e, T
	return decayRate(double(n1), getMass(double(n1)), double(n2), getMass(double(n2)), fabs(getE(double(n1+n2))-getE(double(n1))-getE(double(n2))), gp, Temper);
}

dvector ADMmodel::getRatesForN(double n)
{
	if (n == maxN)
	{
		dvector out;
		return out;
	}
	dvector::iterator it = std::find(num.begin(), num.end(), n);
	size_t ind = std::distance(num.begin(), it);
	dvector out;
	if (it != num.end() && ind == 0)
	{
		dvector temp(prodRates.begin(), prodRates.begin()+1);
		out = temp;
	}
	else if (it != num.end())
	{
		dvector temp(&prodRates[ind*(ind+1)/2], &prodRates[(ind*(ind+3)/2) + 1]);
		out = temp;
	}
	
	return out;
}

dvector ADMmodel::getDecaysForN(double n)
{
	if (n == maxN)
	{
		dvector out;
		return out;
	}
	dvector::iterator it = std::find(num.begin(), num.end(), n);
	size_t ind = std::distance(num.begin(), it);
	dvector out;
	if (it != num.end() && ind == 0)
	{
		dvector temp(decayRates.begin(), decayRates.begin()+1);
		out = temp;
	}
	else if (it != num.end())
	{
		dvector temp(&decayRates[ind*(ind+1)/2], &decayRates[(ind*(ind+3)/2) + 1]);
		out = temp;
	}
	
	return out;
}

dvector ADMmodel::getAllRatesForN(double n)
{
	if (n >= maxN)
	{
		dvector out;
		return out;
	}
	dvector::iterator it = std::find(num.begin(), num.end(), n);
	size_t ind = std::distance(num.begin(), it);
	dvector out;
	if (it != num.end() && ind == 0)
	{
		for (int k = 0; k < D-1; k++)
		{
			ind = k+ind;
			out.push_back(prodRates[ind]);
		}
	}
	else if (it != num.end())
	{
		dvector temp(&prodRates[ind*(ind+1)/2], &prodRates[(ind*(ind+3)/2) + 1]);
		out = temp;
		size_t ind2 = ind;
		ind = ind*(ind+3)/2;
		for (size_t k = ind2+1; k<D-1; k++)
		{
			ind = ind + k;
			out.push_back(prodRates[ind]);
		}
	}
	return out;
}

double ADMmodel::getMass(double q)
{
	if(bMadeSplines)
		return gsl_spline_eval(massSpline, q, mAcc);
	else if (D > 1)
	{
		makeSplines(0.0);
		return gsl_spline_eval(massSpline, q, mAcc);
	}
	else
		return 0.0;
}

void ADMmodel::display_getMass(double q)
{
	printf("%f\n", getMass(q));
	return;
}

double ADMmodel::getE(double q)
{
	if(bMadeSplines)
		return gsl_spline_eval(ESpline, q, eAcc);
	else if (D > 1)
	{
		makeSplines(0.0);
		return gsl_spline_eval(ESpline, q, eAcc);
	}
	else
		return 0.0;
}

void ADMmodel::display_getE(double q)
{
	printf("%f\n", getE(q));
	return;
}

double ADMmodel::getR(double q)
{
	if(bMadeSplines)
		return gsl_spline_eval(radSpline, q, rAcc);
	else if (D > 1)
	{
		makeSplines(0.0);
		return gsl_spline_eval(radSpline, q, rAcc);
	}
	else
		return 0.0;
}

void ADMmodel::display_getR(double q)
{
	printf("%f\n", getR(q));
	return;
}

double ADMmodel::getPhi(double q)
{
	if(bMadeSplines)
		return gsl_spline_eval(phiSpline, q, pAcc);
	else if (D > 1)
	{
		makeSplines(0.0);
		return gsl_spline_eval(phiSpline, q, pAcc);
	}
	else
		return 0.0;
}

void ADMmodel::display_getPhi(double q)
{
	printf("%f\n", getPhi(q));
	return;
}

double ADMmodel::getRateFromSpline(double p, double q)
{
	if((p < numRange[0] || p > numRange[1]) || (q < numRange[0] || q > numRange[1]))
		return 0.0;
	
	if(q > p)
	{
		double temp = q;
		q = p;
		p = temp;
	}
	return exp(rateSpline.eval(log(p), log(q)));
}

double ADMmodel::getDecayFromSpline(double p, double q)
{
	if((p < numRange[0] || p > numRange[1]) || (q < numRange[0] || q > numRange[1]))
		return 0.0;
	
	if(q > p)
	{
		double temp = q;
		q = p;
		p = temp;
	}
	return exp(decaySpline.eval(log(p), log(q)));
}

void ADMmodel::writeBoundProps(double Temp)
{
	if (Temp == 0.0)
		Temp = T;
	char tStr[20];
	sprintf(tStr, "_T%.3e", Temp);
	string prop_path = boundstate_properties_path;
	string cs_path = crosssection_path + tStr + ".txt";
	string de_path = decay_path + tStr + ".txt";
	ofstream fout;
	fout.open(prop_path.c_str(), ios::trunc);
	dvector rates;
	for(int k = 0; k < D; k++)
	{
		fout << num[k] << ",";
		fout << masses[k] << ",";
		fout << bindEnergies[k] << ",";
		fout << radii[k] << ",";
		fout << phi0[k] << ",";
		fout << "\n";
	};
	fout.close();
	fout.open(cs_path.c_str(), ios::trunc);
	for(int j = 0; j < D-1; j++)
	{
		rates = getRatesForN(num[j]);
		copy(rates.begin(), rates.end(), std::ostream_iterator<double>(fout, ","));
		fout << "\n";
	};
	fout.close();
	fout.open(de_path.c_str(), ios::trunc);
	for(int j = 0; j < D-1; j++)
	{
		rates = getDecaysForN(num[j]);
		copy(rates.begin(), rates.end(), std::ostream_iterator<double>(fout, ","));
		fout << "\n";
	};
	fout.close();
	return;
}

void ADMmodel::readBoundProps(double Temp)
{
	if (Temp == 0.0)
		Temp = T;
	char tStr[20];
	sprintf(tStr, "_T%.3e", Temp);
	string cs_path = crosssection_path + tStr + ".txt";
	string de_path = decay_path + tStr + ".txt";
	string path = boundstate_properties_path;
	ifstream fin;
	fin.open(path.c_str());
	string line;
	if (fin.is_open())
	{
	istringstream ss;
	string field;
	num.clear();
	masses.clear();
	bindEnergies.clear();
	radii.clear();
	phi0.clear();
	prodRates.clear();
	decayRates.clear();
	maxN = 0;
	D = 0;
	printf("Loading bound state data...\n\n");
		while( getline(fin, line) )
		{
			ss.str(line);
			getline(ss, field, ',');
			num.push_back(atof(field.c_str()));
			if(maxN < atof(field.c_str()))
				maxN = atof(field.c_str());
			getline(ss, field, ',');
			masses.push_back(atof(field.c_str()));
			getline(ss, field, ',');
			bindEnergies.push_back(atof(field.c_str()));
			getline(ss, field, ',');
			radii.push_back(atof(field.c_str()));
			getline(ss, field, ',');
			phi0.push_back(atof(field.c_str()));
			D+=1;
		}
		fin.close();
		makeSplines(0.0);
		fin.open(cs_path.c_str());
		if (fin.is_open())
		{
			printf("Loading production rate data...\n\n");
			int poi = 1;
			while( getline(fin, line) )
			{
				ss.str(line);
				for(int l = 0; l < poi; l++)
				{
					getline(ss, field, ',');
					prodRates.push_back(atof(field.c_str()));
				}
				poi++;
			}
			fin.close();
			if (prodRates.size() == (D*(D-1))/2)
			{
				numRange[0] = *std::min_element(num.begin(), num.end());
				numRange[1] = *std::max_element(num.begin(), num.end());
				dvector x;
				dvector y;
				for(int j = 0; j < D; j++)
				{
					for(int k = 0; k<=j; k++)
					{
						x.push_back((double)log(num[j]));
						y.push_back((double)log(num[k]));
					}
				}
				dvector ratesOut;
				for(int jk = 0; jk < prodRates.size(); jk++)
				{
					ratesOut.push_back((double)log(prodRates[jk]));
				}
				rateSpline.reset(x, y, ratesOut, (double)lambda);
			}
			else
				getAllProductionRates(Temp);
			
		}
		else
		{
			printf("No production rate data found for these paramaters at this temperature.\n\n");
			printf("Temperature = %.3e, ", Temp);
			printf("Calculating...\n\n");
			getAllProductionRates(Temp);
		}
		
		fin.open(de_path.c_str());
		if (fin.is_open())
		{
			printf("Loading decay rate data...\n\n");
			int poi = 1;
			while( getline(fin, line) )
			{
				ss.str(line);
				for(int l = 0; l < poi; l++)
				{
					getline(ss, field, ',');
					decayRates.push_back(atof(field.c_str()));
				}
				poi++;
			}
			fin.close();
			if (decayRates.size() == (D*(D-1))/2)
			{
				numRange[0] = *std::min_element(num.begin(), num.end());
				numRange[1] = *std::max_element(num.begin(), num.end());
				dvector x;
				dvector y;
				for(int j = 0; j < D; j++)
				{
					for(int k = 0; k<=j; k++)
					{
						x.push_back((double)log(num[j]));
						y.push_back((double)log(num[k]));
					}
				}
				dvector decaysOut;
				for(int jk = 0; jk < decayRates.size(); jk++)
				{
					decaysOut.push_back((double)log(decayRates[jk]));
				}
				decaySpline.reset(x, y, decaysOut, (double)lambda);
			}
			else
				getAllDecayRates(Temp);
			
		}
		else
		{
			printf("No decay rate data found for these paramaters at this temperature.\n\n");
			printf("Temperature = %.3e, ", Temp);
			printf("Calculating...\n\n");
			getAllDecayRates(Temp);
		}
	}
	else
		printf("No file found. \n\n");
	
	return;
}