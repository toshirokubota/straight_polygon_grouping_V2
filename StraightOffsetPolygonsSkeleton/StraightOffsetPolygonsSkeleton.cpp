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
#include <set>
#include <hash_map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <DisjointSet.h>
#include <IntersectionConvexPolygons.h>
#include <GraphFactory.h>
#include <MovingParticle.h>
#include <ParticleSimulatorSkeleton.h>

//GraphFactory<StationaryParticle*>* GraphFactory<StationaryParticle*>::_instance = NULL;
//ParticleFactory* ParticleFactory::_instance = NULL;

int MovingParticle::_id = 0;
int Polygon::_id = 0;

vector<pair<int, int>>
indices2pairs(vector<int> T, const int* dims)
{
	int n = dims[0];
	vector<pair<int, int>> pairs(n);
	for (int i = 0; i<n; ++i)
	{
		int k1 = GetData2(T, i, 0, dims[0], dims[1], 0);
		int k2 = GetData2(T, i, 1, dims[0], dims[1], 0);
		pairs[i] = pair<int, int>(k1 - 1, k2 - 1);
	}
	return pairs;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "StraightMedialAxis", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = StraightMedialAxis(P, [iter delta])");
		return;
	}
	ParticleSimulatorSkeleton simulator;
	//Points
	vector<StationaryParticle*> points;
	const int* dimsP;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);

		StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
		for (int i = 0; i < dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			points.push_back(sfactory.makeParticle(CParticleF(x, y)));
		}
		//edges (indices to the points)
		vector<pair<int, int>> E;
		{
			vector<int> T0;
			mxClassID classIdT;
			int ndimT;
			const int* dimsT;
			LoadData(T0, prhs[1], classIdT, ndimT, &dimsT);
			E = indices2pairs(T0, dimsT);
		}

		float initTime = 0.01f;
		if (nrhs >= 6)
		{
			mxClassID classMode;
			ReadScalar(initTime, prhs[5], classMode);
		}
		simulator.Prepare(points, E, initTime);
	}


	float thres = 0.99f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[2], classMode);
	}
	simulator.thres = thres;
	simulator.Simulate();

	if (nlhs >= 1)
	{
		const int dims[] = { simulator.skeletons.size(), 8 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < simulator.skeletons.size(); ++i)
		{
			ParticleSkeleton ps = simulator.skeletons[i];
			SetData2(F, i, 0, dims[0], dims[1], ps.p->getP().m_X);
			SetData2(F, i, 1, dims[0], dims[1], ps.p->getP().m_Y);
			SetData2(F, i, 2, dims[0], dims[1], ps.q->getP().m_X);
			SetData2(F, i, 3, dims[0], dims[1], ps.q->getP().m_Y);
			SetData2(F, i, 4, dims[0], dims[1], ps.o.m_X);
			SetData2(F, i, 5, dims[0], dims[1], ps.o.m_Y);
			SetData2(F, i, 6, dims[0], dims[1], (float)ps.p->getId());
			SetData2(F, i, 7, dims[0], dims[1], (float)ps.q->getId());
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}

	ParticleFactory::getInstance().clean();
	mexUnlock();
}

