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
#include <ParticleSimulator.h>

GraphFactory<StationaryParticle*>* GraphFactory<StationaryParticle*>::_instance = NULL;
ParticleFactory* ParticleFactory::_instance = NULL;

int MovingParticle::_id = 0;

mxArray*
StoreContours(const vector<vector<CParticleF>>& polygons)
{
	const int dims[] = {polygons.size()};
	mxArray* cell = mxCreateCellArray(1, (mwSize*) dims);
	for(int i=0; i<polygons.size(); ++i)
	{
		int n = polygons[i].size() ;
		const int dimsC[] = {n, 3};
		mxArray* ar = mxCreateNumericArray(2, (mwSize*) dimsC, mxSINGLE_CLASS, mxREAL);
		float* p = (float*) mxGetData(ar);
		for(int j=0; j<n; ++j)
		{
			p[j] = polygons[i][j].m_X;
			p[n+j] = polygons[i][j].m_Y;
			p[2*n+j] = polygons[i][j].m_Z;
		}
		mxSetCell(cell, i, ar);
	}
	return cell;
}

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

vector<vector<MovingParticle*>> 
collectPolygons()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<vector<MovingParticle*>> polygons;
	set<MovingParticle*> pset;
	while (true)
	{
		bool bdone = true;
		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (pset.find(p) == pset.end())
			{
				vector<MovingParticle*> v = MovingParticle::vectorize(p);
				for (int i = 0; i < v.size(); ++i)
				{
					pset.insert(v[i]);
				}
				polygons.push_back(v);
				bdone = false;
			}
		}
		if (bdone) break;
	}
	return polygons;
}

float
collisionTime(CParticleF& o, MovingParticle* p, MovingParticle* q)
{
	ParticleDirection u = p->getFront().normal();
	float v[2];
	p->getVelocity(v[0], v[1]);
	float dq = v[0] * u.x + v[1]* u.y;

	CParticleF y = Closest2Line(p->getP(), q->getP(), o);
	float deriv = (y.m_X - o.m_X)*u.x + (y.m_Y - o.m_Y) * u.y;
	if (deriv > 0) return std::numeric_limits<float>::infinity(); //moving away

	float dval = Distance(o, y);
	return dval / dq;
}

/*bool
onSideAt(CParticleF& o, MovingParticle* q, MovingParticle* r, float t, float eps)
{
	CParticleF q2 = q->move(t);
	CParticleF r2 = r->move(t);
	float d = Distance2LineSegment(q2, r2, o);
	return d <= eps;
}*/

bool
onSideAt(CParticleF& o, MovingParticle* q, MovingParticle* r, float t, float eps)
{
	CParticleF y = Closest2Line(q->getP(), r->getP(), o);
	float df = Distance(y, q->getP()) + Distance(y, r->getP()) - Distance(q->getP(), r->getP());
	return df <= eps;
}

float
energyFromBoundary(CParticleF&  o, vector<MovingParticle*>& vp, float thres)
{
	float sum = 0;
	for (int i = 0; i < vp.size(); ++i)
	{
		MovingParticle* p = vp[i];
		MovingParticle* q = p->getNext();
		float t = collisionTime(o, p, q);
		if (onSideAt(o, p, q, t, 1.0e-3))
		{
			sum += Max(0, thres - t);
		}
	}
	return sum;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "StraightMedialAxis", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = StraightMedialAxis(P, [iter delta])");
		return;
	}
	ParticleSimulator simulator;
	//Points
	vector<StationaryParticle*> points;
	const int* dimsP;
	vector<float> P0;
	mxClassID classIdP;
	int ndimP;
	LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
	float maxX = 0, maxY = 0;
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	for (int i = 0; i < dimsP[0]; ++i)
	{
		float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
		float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
		points.push_back(sfactory.makeParticle(CParticleF(x, y)));
		maxX = Max(maxX, x);
		maxY = Max(maxY, y);
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

	float endtime = 2.0f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(endtime, prhs[2], classMode);
	}
	simulator.Prepare(points, E, 0);

	const int dimsA[] = { (int)(maxX + endtime + 10), (int)(maxY + endtime + 10) };
	vector<float> A(dimsA[0] * dimsA[1]);

	vector<vector<MovingParticle*>> polygons = collectPolygons();
	vector<vector<CParticleF>> bounds;
	for (int i = 0; i < polygons.size(); ++i)
	{
		vector<CParticleF> bound;
		for (int j = 0; j < polygons[i].size(); ++j)
		{
			bound.push_back(polygons[i][j]->getP());
		}
		bounds.push_back(bound);
	}

	for (int y = 0; y < dimsA[1]; ++y)
	{
		for (int x = 0; x < dimsA[0]; ++x)
		{
			CParticleF o(x, y);
			if (x == 50 && y == 50)
				x += 0;
			float a = GetData2(A, x, y, dimsA[0], dimsA[1], 0.0f);
			float dmin = std::numeric_limits<float>::infinity();
			for (int i = 0; i < polygons.size(); ++i)
			{
				//if (inside(o, bounds[i]))
				{
					float d = energyFromBoundary(o, polygons[i], endtime);
					//dmin = Min(d, dmin);
					a += d; // Max(0, endtime - d);
				}
			}
			//a = Max(0, endtime - dmin);
			SetData2(A, x, y, dimsA[0], dimsA[1], a);
		}
	}

	if (nlhs >= 1) 
	{
		plhs[0] = StoreData(A, mxSINGLE_CLASS, 2, dimsA);

	}

	ParticleFactory::getInstance()->clean();
	mexUnlock();
}

