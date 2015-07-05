#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>

#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szParticleF.h>
#include <Kmeans.h>
#include <EMGaussian.h>
#include <DisjointSet.h>

vector<int>
initialClustering(vector<Feature>& P, float space)
{
	vector<Node<int>*> nodes;
	for (int i = 0; i < P.size(); ++i)
	{
		nodes.push_back(makeset(i));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			if (Feature::distance(P[i], P[j]) <= space)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<Node<int>*> c = clusters(nodes);
	vector<int> labels(nodes.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = distance(c.begin(), find(c.begin(), c.end(), findset(nodes[i])));
		labels[i] = k;
	}
	return labels;
}

/*
Divide the feature set (P) into multiple bins, according to the labeling (LABEL).
Small clusters less than (MINSIZE) are not kept.
*/
vector<vector<Feature>> divideClusters(vector<Feature>& P, vector<int>& labels, int minSize)
{
	int numC = 0;
	for (int i = 0; i < labels.size(); ++i)
	{
		numC = Max(labels[i], numC);
	}
	vector<int> counts(numC+1, 0);
	for (int i = 0; i < labels.size(); ++i)
	{
		counts[labels[i]]++;
	}
	vector<vector<Feature>> bins(numC + 1);
	for (int lb = 0; lb <= numC; ++lb)
	{
		if (counts[lb] >= minSize)
		{
			for (int j = 0; j < labels.size(); ++j)
			{
				if (labels[j] == lb)
				{
					bins[lb].push_back(P[j]);
				}
			}
		}
	}
	for (int i = bins.size() - 1; i >= 0; i--)
	{
		if (bins[i].empty())
		{
			bins.erase(bins.begin() + i);
		}
	}
	return bins;
}

float 
AIC(vector<vector<Feature>>& features)
{
	float sum = 0;
	float k = 5.0f; //# of parameters for a 2D Gaussian

	for (int i = 0; i < features.size(); ++i)
	{
		if (features[i].size() <= 6) return std::numeric_limits<float>::infinity();

		GaussianDistribution dist(features[i]);
		float sum2 = 0;
		for (int j = 0; j < features[i].size(); ++j)
		{
			sum2 -= log(dist.eval(features[i][j]));
		}
		sum2 += k + k*(k+1)/(features[i].size()-k-1);
		sum += sum2;
	}
	return sum;
}

/*
*/
vector<vector<Feature>>
recursiveSubdivision(vector<Feature>& P, int minSize)
{
	int K = 2;
	int iter = 3;
	vector<vector<Feature>> clusters(1, P);
	if (P.size() < minSize)
	{
		return clusters; //empty
	}
	else {
		vector<GaussianDistribution> dist = EMGaussianClustering(P, K, iter);
		vector<vector<Feature>> C(2);
		for (int j = 0; j < P.size(); ++j)
		{
			int sel = selectGaussianCluster(P[j], dist);
			C[sel].push_back(P[j]);
		}
		float a1 = AIC(clusters);
		float a2 = AIC(C);
		if (a1 > a2)
		{
			vector<vector<Feature>> C2;
			for (int j = 0; j < K; ++j)
			{
				vector<vector<Feature>> clusters2 = recursiveSubdivision(C[j], minSize);
				C2.insert(C2.end(), clusters2.begin(), clusters2.end());
			}
			return C2;
		}
		else
		{
			return clusters;
		}
	}
}

float
GoodnessFit(GaussianDistribution& d, vector<Feature>& P)
{
	std::valarray<float> evals;
	techsoft::matrix<float> evecs;
	if (d.sgm.mat.eigen(evals, evecs))
	{
		/*CParticleF a(d.mu.vals[0], d.mu.vals[1]);
		CParticleF b;
		if (evals[0] > evals[1])
		{
			b.m_X = evecs[0][0]; b.m_Y = evecs[1][0];
		}
		else
		{
			b.m_X = evecs[0][1]; b.m_Y = evecs[1][1];
		}
		float maxD = 0;
		for (int i = 0; i < P.size(); ++i)
		{
			CParticleF x(P[i].vals[0], P[i].vals[1]);
			float d = Distance2Line(a, b, x);
			if (d > maxD)
			{
				maxD = d;
			}
		}*/
		return max(evals[0], evals[1]) / min(evals[0], evals[1]);
	}
	else
	{
		return 0;
	}
}

/*
*/
vector<vector<Feature>>
recursiveSubdivision2(vector<Feature>& P, int minSize, float thres)
{
	int K = 2;
	int iter = 3;
	vector<vector<Feature>> clusters(1, P);
	if (P.size() < minSize)
	{
		return clusters; //empty
	}
	else {
		GaussianDistribution d(P);
		float fit = GoodnessFit(d, P);
		if (fit < thres)
		{
			vector<GaussianDistribution> dist = EMGaussianClustering(P, K, iter);
			vector<vector<Feature>> C(K);
			for (int j = 0; j < P.size(); ++j)
			{
				int sel = selectGaussianCluster(P[j], dist);
				C[sel].push_back(P[j]);
			}
			return C;
		}
		else
		{
			return clusters;
		}
	}
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	//printf("%s: This build was compiled at %s %s\n", "OverlapArea", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: L = EMGClustering(P, K)");
		return;
	}
	//Points
	vector<Feature> P;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for (int i = 0; i<dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			Feature f(2);
			f.vals[0] = x;
			f.vals[1] = y;
			P.push_back(f);
		}
	}
	float space = 10.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(space, prhs[1], classMode);
	}
	int minSize = 2;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minSize, prhs[2], classMode);
	}

	techsoft::matrix<float> m1(2, 2, 0.0f), m2;
	m1[0][0] = 1; m1[0][1] = m1[1][0] = 2; m1[1][1] = 3;
	std::valarray<float> ev;
	m1.eigen(ev, m2);
	float evec[4] = { m2[0][0], m2[0][1], m2[1][0], m2[1][1] };

	vector<int> labels = initialClustering(P, space);
	vector<vector<Feature>> clusters = divideClusters(P, labels, minSize);

	vector<vector<Feature>> clusters2;
	for (int i = 0; i < clusters.size(); ++i)
	{
		//vector<vector<Feature>> clusters0 = recursiveSubdivision(clusters[i], 10);
		vector<vector<Feature>> clusters0 = recursiveSubdivision2(clusters[i], 10, 8);
		clusters2.insert(clusters2.end(), clusters0.begin(), clusters0.end());
	}

	if (nlhs >= 1)
	{
		const int dims[] = { clusters2.size(), 4 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < clusters2.size(); ++i)
		{
			if (clusters2[i].size() >= 2)
			{
				GaussianDistribution d(clusters2[i]);
				std::valarray<float> ev;
				techsoft::matrix<float> evec;
				if (d.sgm.mat.eigen(ev, evec))
				{
					if (ev[0] > ev[1])
					{
						float z = sqrt(ev[0]);
						SetData2(F, i, 0, dims[0], dims[1], (d.mu.vals[0] - z * evec[0][0]));
						SetData2(F, i, 1, dims[0], dims[1], (d.mu.vals[1] - z * evec[1][0]));
						SetData2(F, i, 2, dims[0], dims[1], (d.mu.vals[0] + z * evec[0][0]));
						SetData2(F, i, 3, dims[0], dims[1], (d.mu.vals[1] + z * evec[1][0]));
					}
					else
					{
						float z = sqrt(ev[1]);
						SetData2(F, i, 0, dims[0], dims[1], (d.mu.vals[0] - z * evec[0][1]));
						SetData2(F, i, 1, dims[0], dims[1], (d.mu.vals[1] - z * evec[1][1]));
						SetData2(F, i, 2, dims[0], dims[1], (d.mu.vals[0] + z * evec[0][1]));
						SetData2(F, i, 3, dims[0], dims[1], (d.mu.vals[1] + z * evec[1][1]));
					}
				}
			}
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		vector<Feature> Q;
		for (int i = 0; i < clusters2.size(); ++i)
		{
			for (int j = 0; j < clusters2[i].size(); ++j)
			{
				Q.push_back(clusters2[i][j]);
			}
		}
		const int dims[] = { Q.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], Q[i].vals[0]);
			SetData2(F, i, 1, dims[0], dims[1], Q[i].vals[1]);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		vector<int> labels2;
		for (int i = 0; i < clusters2.size(); ++i)
		{
			for (int j = 0; j < clusters2[i].size(); ++j)
			{
				labels2.push_back(i);
			}
		}
		const int dims[] = { labels2.size(), 1 };
		plhs[2] = StoreData(labels2, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
}

