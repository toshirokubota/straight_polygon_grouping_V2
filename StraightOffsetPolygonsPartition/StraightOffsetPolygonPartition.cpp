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
#include <ParticleSimulatorGreedyPartition.h>

int MovingParticle::_id = 0;
int Polygon::_id = 0;


vector<pair<int,int>>
indices2pairs(vector<int> T, const int* dims)
{
	int n = dims[0];
	vector<pair<int,int>> pairs(n);
	for(int i=0; i<n; ++i)
	{
		int k1 = GetData2(T, i, 0, dims[0], dims[1], 0);
		int k2 = GetData2(T, i, 1, dims[0], dims[1], 0);
		pairs[i] = pair<int,int>(k1-1, k2-1);
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
	ParticleSimulatorGreedyPartition simulator;
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

		float initTime = 0.0f; //keep this to 0 for accurately extracting closed polygons
		/*if (nrhs >= 6)
		{
			mxClassID classMode;
			ReadScalar(initTime, prhs[5], classMode);
		}*/
		simulator.Prepare(points, E, initTime);
	}
	float endtime = 2.0f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(endtime, prhs[2], classMode);
	}
	float delta = 0.1f;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(delta, prhs[3], classMode);
	}
	bool bdebug = false;
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(bdebug, prhs[4], classMode);
	}

	try 
	{
		simulator.Simulate(endtime, delta, bdebug);
	}
	catch (exception ex)
	{
		mexErrMsgTxt("0. Exception in ParticleSimulator::Simulate");
	}
	if (nlhs >= 1) 
	{
		plhs[0] = Snapshot::StoreSnapshots(simulator.snapshots);
	}
	if (nlhs >= 2)
	{
		plhs[1] = simulator.SaveParticles();
	}
	if (nlhs >= 3)
	{
		plhs[2] = Snapshot::StoreSnapshots0(simulator.closedRegions);
	}
	if (nlhs >= 4)
	{
		//plhs[3] = Snapshot::StoreSnapshots0(simulator.closedRegions);
		plhs[3] = Snapshot::StoreSnapshots(simulator.polygons);
		//plhs[3] = simulator.SaveDoneEvents();
	}
	if (nlhs >= 5)
	{
		plhs[4] = Snapshot::StoreSnapshots(simulator.traces);
		//plhs[3] = simulator.SaveDoneEvents();
	}
	if (nlhs >= 6)
	{
		//vector<Snapshot> snapshots = simulator.closedRegions;
		//plhs[5] = Snapshot::StoreSnapshots(snapshots);
		plhs[5] = simulator.SaveDoneEvents();
	}

	ParticleFactory::getInstance().clean();
	PolygonFactory::getInstance().clean();
	mexUnlock();
}

