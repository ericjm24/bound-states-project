#ifndef DVECTOR_H
#define DVECTOR_H

#include <vector>
#include <assert.h>

class dvector : public std::vector<double>
{
	public:
		dvector(std::vector<double>::iterator itS, std::vector<double>::iterator itE) : std::vector<double>(itS, itE) {};
		dvector() : std::vector<double>() {};
		dvector(double * itS, double * itE) : std::vector<double>(itS, itE) {};
		double & operator[](int n);
		const double operator[](int n) const;
		dvector operator+(const dvector & v) const;
		dvector operator-(const dvector & v) const;
		dvector operator*(const dvector & v) const;
		dvector operator/(const dvector & v) const;
		dvector operator+(double y) const;
		friend dvector operator+(double y, const dvector & v);
		dvector operator-(double y) const;
		friend dvector operator-(double y, const dvector & v);
		dvector operator*(double y) const;
		friend dvector operator*(double y, const dvector & v);
		dvector operator/(double y) const;
		friend dvector operator/(double y, const dvector & v);
		dvector & operator+=(double y);
		dvector & operator-=(double y);
		dvector & operator*=(double y);
		dvector & operator/=(double y);
		dvector & operator+=(const dvector & v);
		dvector & operator-=(const dvector & v);
		dvector & operator*=(const dvector & v);
		dvector & operator/=(const dvector & v);
		dvector operator-() const;
	//private:
		
};

dvector operator+(double y, const dvector & v);
dvector operator-(double y, const dvector & v);
dvector operator*(double y, const dvector & v);
dvector operator/(double y, const dvector & v);

#endif