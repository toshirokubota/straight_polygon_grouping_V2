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
#include <StationaryParticle.h>
#include <OffsetPolygonDAGBuilder.h>

GraphFactory<StationaryParticle*>* GraphFactory<StationaryParticle*>::_instance = NULL;
ParticleFactory* ParticleFactory::_instance = NULL;

int MovingParticle::_id = 0;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "StraightMedialAxis", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = StraightMedialAxis(P, [iter delta])");
		return;
	}
	OffsetPolygonDAGBuilder builder(prhs[0], prhs[1]);
	float endtime = 2.0f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(endtime, prhs[2], classMode);
	}

	builder.Polygonify();
	builder.Expand(endtime);

	if (nlhs >= 1)
	{
		plhs[0] = builder.dag.PointsToMxArray();
	}
	if (nlhs >= 2)
	{
		plhs[1] = builder.dag.EdgesToMxArray();
	}


	ParticleFactory::getInstance().clean();
	mexUnlock();
}

