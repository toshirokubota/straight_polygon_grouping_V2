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
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <FragmentInfo.h>
#include <DotGroupUtility.h>
#include <Kmeans.h>

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

	vector<Feature> means = KmeansClustering(P, K);
	vector<int> labels(P.size());
	for (int i = 0; i < labels.size(); ++i)
	{
		labels[i] = selectKmeansCluster(P[i], means);
	}

	if (nlhs >= 1)
	{
		const int dims[] = { labels.size(), 1 };
		plhs[0] = StoreData(labels, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
}

