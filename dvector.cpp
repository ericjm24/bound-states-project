#include "dvector.h"

double & dvector::operator[](int n)
{
	if(n < size())
	{
		return at(n);
	}
	else
	{
		resize(n+1);
		return at(n);
	}
	return at(n);
}

const double dvector::operator[](int n) const
{
	if(n < size())
	{
		return at(n);
	}
	else
	{
		return 0.0;
	}
	return at(n);
}

dvector dvector::operator+(const dvector & v) const
{
	dvector out;
	int siz1, siz2;
	if(siz1 > siz2)
	{
		out.assign(siz1, 0.0);
	}
	else
	{
		out.assign(siz2, 0.0);
	}
	for(int k1; k1 < siz1; k1++)
	{
		out[k1] += operator[](k1);
	}
	for(int k2; k2 < siz2; k2++)
	{
		out[k2] += v[k2];
	}
	return out;
}

dvector dvector::operator-(const dvector & v) const
{
	dvector out;
	int siz1, siz2;
	if(siz1 > siz2)
	{
		out.assign(siz1, 0.0);
	}
	else
	{
		out.assign(siz2, 0.0);
	}
	for(int k1; k1 < siz1; k1++)
	{
		out[k1] += operator[](k1);
	}
	for(int k2; k2 < siz2; k2++)
	{
		out[k2] -= v[k2];
	}
	return out;
}

dvector dvector::operator*(const dvector & v) const
{
	dvector out;
	int siz1, siz2;
	if(siz1 > siz2)
	{
		out.assign(siz1, 0.0);
		for(int k = 0; k < siz2; k++)
		{
			out[k] = operator[](k)*v[k];
		}
	}
	else
	{
		out.assign(siz2, 0.0);
		for(int k = 0; k < siz1; k++)
		{
			out[k] = operator[](k)*v[k];
		}
	}
	return out;
}

dvector dvector::operator/(const dvector & v) const
{
	dvector out;
	int siz1, siz2;
	siz1 = size();
	siz2 = v.size();
	assert(siz1 == siz2);
	out.resize(siz1);
	for(int k = 0; k < siz1; k++)
	{
		out[k] = operator[](k)/v[k];
	}
	return out;
}

dvector dvector::operator+(double y) const
{
	dvector out = *this;
	for(int k = 0; k < size(); k++)
	{
		out[k] += y;
	}
	return out;
}

dvector dvector::operator-(double y) const
{
	dvector out = *this;
	for(int k = 0; k < size(); k++)
	{
		out[k] -= y;
	}
	return out;
}

dvector dvector::operator*(double y) const
{
	dvector out = *this;
	for(int k = 0; k < size(); k++)
	{
		out[k] *= y;
	}
	return out;
}

dvector dvector::operator/(double y) const
{
	dvector out = *this;
	for(int k = 0; k < size(); k++)
	{
		out[k] /= y;
	}
	return out;
}

dvector operator+(double y, const dvector & v)
{
	dvector out = v;
	for(int k = 0; k < v.size(); k++)
	{
		out[k] += y;
	}
	return out;
}

dvector operator-(double y, const dvector & v)
{
	dvector out;
	out.assign(v.size(), y);
	for(int k = 0; k < v.size(); k++)
	{
		out[k] -= v[k];
	}
	return out;
}

dvector operator*(double y, const dvector & v)
{
	dvector out = v;
	for(int k = 0; k < v.size(); k++)
	{
		out[k] *= y;
	}
	return out;
}

dvector operator/(double y, const dvector & v)
{
	dvector out;
	out.assign(v.size(), y);
	for(int k = 0; k < v.size(); k++)
	{
		out[k] /= v[k];
	}
	return out;
}

dvector & dvector::operator+=(double y)
{
	for(int k = 0; k < size(); k++)
	{
		at(k) += y;
	}
	return *this;
}

dvector & dvector::operator-=(double y)
{
	for(int k = 0; k < size(); k++)
	{
		at(k) -= y;
	}
	return *this;
}

dvector & dvector::operator*=(double y)
{
	for(int k = 0; k < size(); k++)
	{
		at(k) *= y;
	}
	return *this;
}

dvector & dvector::operator/=(double y)
{
	for(int k = 0; k < size(); k++)
	{
		at(k) *= y;
	}
	return *this;
}

dvector & dvector::operator+=(const dvector & v)
{
	if(v.size() > size())
	{
		resize(v.size());
	}
	for(int k = 0; k < v.size(); k++)
	{
		if(k < size())
		{
			at(k) += v[k];
		}
		else
		{
			at(k) = v[k];
		}
	}
	return *this;
}

dvector & dvector::operator-=(const dvector & v)
{
	if(v.size() > size())
	{
		resize(v.size());
	}
	for(int k = 0; k < v.size(); k++)
	{
		if(k < size())
		{
			at(k) -= v[k];
		}
		else
		{
			at(k) = -v[k];
		}
	}
	return *this;
}

dvector & dvector::operator*=(const dvector & v)
{
	if(v.size() > size())
	{
		resize(v.size());
	}
	for(int k = 0; k < v.size(); k++)
	{
		if(k < size())
		{
			at(k) += v[k];
		}
		else
		{
			at(k) = 0.0;
		}
	}
	return *this;
}

dvector & dvector::operator/=(const dvector & v)
{
	assert(v.size() == size());
	for(int k = 0; k < size(); k++)
	{
		at(k) /= v[k];
	}
	return *this;
}

dvector dvector::operator-() const
{
	dvector out = *this;
	for(int k = 0; k < size(); k++)
	{
		out[k] = -at(k);
	}
	return out;
}