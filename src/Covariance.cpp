#include <Covariance.h>


Covariance 
Covariance::computeCovariance(const vector<Feature>& x, Feature m)
{
	int ndim = m.length();
	math::matrix<float> c(ndim, ndim, 0);
	float vals0[] = { c(0, 0), c(0, 1), c(1, 0), c(1, 1) };
	for (int i = 0; i < x.size(); ++i)
	{
		math::matrix<float> df(m.length(), 1);
		for (int j = 0; j < ndim; ++j)
		{
			df[j][0] = x[i].vals[j] - m.vals[j];			
		}
		c = c + df * (~df);
	}
	c /= x.size();

	float vals[] = { c(0, 0), c(0, 1), c(1, 0), c(1, 1) };
	math::matrix<float> ic = c.Inv();
	Covariance cov(ndim);
	cov.mat = c;
	float vals2[] = { cov.mat(0, 0), cov.mat(0, 1), cov.mat(1, 0), cov.mat(1, 1) };
	cov.imat = ic;
	return cov;
}
