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

/*
Select the closest cluster center based on Euclidean distance.
*/
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
	int K = 2;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(K, prhs[1], classMode);
	}
	int iter = 0;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(iter, prhs[2], classMode);
	}

	vector<GaussianDistribution> dist = EMGaussianClustering(P, K, iter);
	vector<int> labels(P.size());
	for (int i = 0; i < labels.size(); ++i)
	{
		labels[i] = selectGaussianCluster(P[i], dist);
	}

	if (nlhs >= 1)
	{
		const int dims[] = { labels.size(), 1 };
		plhs[0] = StoreData(labels, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { dist.size(), 6 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], dist[i].mu.vals[0]);
			SetData2(F, i, 1, dims[0], dims[1], dist[i].mu.vals[1]);
			SetData2(F, i, 2, dims[0], dims[1], (float)dist[i].sgm.mat(0, 0));
			SetData2(F, i, 3, dims[0], dims[1], (float)dist[i].sgm.mat(0, 1));
			SetData2(F, i, 4, dims[0], dims[1], (float)dist[i].sgm.mat(1, 0));
			SetData2(F, i, 5, dims[0], dims[1], (float)dist[i].sgm.mat(1, 1));
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

