#pragma once
#include <Feature.h>
#include <stdMatrix.h>
#include <vector>
using namespace std;

class Covariance
{
public:
	static Covariance computeCovariance(const vector<Feature>& x, Feature m);
	static Covariance computeCovariance(const vector<Feature>& x)
	{
		return computeCovariance(x, Feature::computeMean(x));
	}

	Covariance(int n=0)
	{
		ndim = n;
	}

	int ndim;
	math::matrix<float> mat;
	math::matrix<float> imat;
};
