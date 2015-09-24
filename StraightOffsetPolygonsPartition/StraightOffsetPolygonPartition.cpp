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

vector<StationaryParticle*>
collectStationaryParticles(Polygon* P)
{
	vector<StationaryParticle*> particles;
	for (int i = 0; i < P->getNumParticles(); ++i)
	{
		particles.push_back(P->getParticle(i)->getInitParticle());
	}
	return particles;
}

vector<CParticleF>
MovingParticles2CParticleF(vector<MovingParticle*>& sp)
{
	vector<CParticleF> cp(sp.size());
	for (int i = 0; i < sp.size(); ++i)
	{
		cp[i] = sp[i]->getP();
	}
	return cp;
}

/*
Two polygons can be merged if they share common sides and one is not completely enclosed inside the other.
*/
bool
checkMergable(Polygon* P, Polygon* Q)
{
	set<StationaryParticle*> sP = P->getStationaryParticleSet();
	set<StationaryParticle*> sQ = Q->getStationaryParticleSet();

	//check for degenerate cases [no overlap, complete inclusion of one by the other]
	set<StationaryParticle*> iset;
	set_intersection(sP.begin(), sP.end(),
		sQ.begin(), sQ.end(),
		std::inserter(iset, iset.begin()));
	if (iset.size()<2) return false;
	set<StationaryParticle*> dset1;
	set_difference(sP.begin(), sP.end(),
		sQ.begin(), sQ.end(),
		std::inserter(dset1, dset1.begin()));
	if (dset1.empty()) return false;
	set<StationaryParticle*> dset2;
	set_difference(sQ.begin(), sQ.end(),
		sP.begin(), sP.end(),
		std::inserter(dset2, dset2.begin()));
	if (dset2.empty()) return false;

	return true;
}


template<class T>
vector<Vertex<T>*>
shortestHops(Vertex<T>* u, Vertex<T>* v)
{
	vector<Vertex<T>*> Q(1, u);
	set<Vertex<T>*> S;
	S.insert(u);
	u->color = Black;
	vector<Vertex<T>*> path;
	bool bfound = false;
	while (Q.empty() == false && bfound == false)
	{
		vector<Vertex<T>*> Q2;
		for (int i = 0; i < Q.size() && !bfound; ++i)
		{
			for (int j = 0; j < Q[i]->aList.size() && !bfound; ++j)
			{
				Vertex<T>* w = Q[i]->aList[j]->v;
				if (w->color == White)
				{
					Q2.push_back(w);
					w->color = Black;
					w->pi = Q[i]->aList[j]->u;
					S.insert(w);
				}
				if (w == v)
				{
					bfound = true;
					Vertex<T>* x = v->pi;
					while (x != NULL)
					{
						path.push_back(x);
						x = x->pi;
					}
				}
			}
		}
		Q = Q2;
	}
	//reset touched vertices
	for (set<Vertex<T>*>::iterator it = S.begin(); it != S.end(); ++it)
	{
		(*it)->Reset();
	}
	return path;
}

vector<Vertex<StationaryParticle*>*>
mergedPolygonGraph(Polygon* P, Polygon* Q)
{
	GraphFactory<StationaryParticle*>& factory = GraphFactory<StationaryParticle*>::GetInstance();
	vector<Vertex<StationaryParticle*>*> G;
	Polygon* polys[2] = { P, Q };
	map<StationaryParticle*, Vertex<StationaryParticle*>*> vmap;
	for (int k = 0; k < 2; ++k)
	{
		for (int i = 0; i < polys[k]->getNumParticles(); ++i)
		{
			MovingParticle* p = polys[k]->getParticle(i);
			StationaryParticle* p0 = p->getInitParticle();
			Vertex<StationaryParticle*>* u = NULL;
			if (vmap.find(p0) == vmap.end())
			{
				u = factory.makeVertex(p0);
				vmap[p0] = u;
				G.push_back(u);
			}
			else
			{
				u = vmap[p0];
			}
		}
	}
	for (int k = 0; k < 2; ++k)
	{
		for (int i = 0; i < polys[k]->getNumParticles(); ++i)
		{
			MovingParticle* p = polys[k]->getParticle(i);
			StationaryParticle* p0 = p->getInitParticle();
			MovingParticle* q = polys[k]->getParticle((i+1) % polys[k]->getNumParticles());
			StationaryParticle* q0 = q->getInitParticle();
			Vertex<StationaryParticle*>* u = vmap[p0];
			Vertex<StationaryParticle*>* v = vmap[q0];
			u->Add(factory.makeEdge(u, v));
			v->Add(factory.makeEdge(v, u));
		}
	}

	return G;
}

vector<StationaryParticle*>
curveBoundary(vector<Vertex<StationaryParticle*>*>& vertices)
{
	vector<CParticleF> pnts;
	for (int i = 0; i < vertices.size(); ++i)
	{
		pnts.push_back(vertices[i]->key->getP());
	}
	vector<CParticleF> hull = ConvexHull2D(pnts);
	vector<Vertex<StationaryParticle*>*> sp;
	for (int i = 0; i < hull.size(); ++i)
	{
		int k = distance(pnts.begin(), find(pnts.begin(), pnts.end(), hull[i]));
		sp.push_back(vertices[k]);
	}
	vector<StationaryParticle*> boundary;
	for (int i = 0; i < hull.size(); ++i)
	{
		Vertex<StationaryParticle*>* u = sp[i];
		Vertex<StationaryParticle*>* v = sp[(i+1) % hull.size()];
		vector<Vertex<StationaryParticle*>*> path = shortestHops(u, v);
		if(path.empty())
		{
			boundary.clear();
			return boundary;
		}
		else
		{
			for (int j = path.size() - 1; j >= 0; --j) 
			{
				boundary.push_back(path[j]->key);
			}
		}
	}
	return boundary;
}

Polygon*
mergePolygons(Polygon* P, Polygon* Q)
{
	if (clockWise(P->getParticles()) > 0 || clockWise(Q->getParticles()) > 0) return NULL;
	if (checkMergable(P, Q) == false) return NULL;

	vector<Vertex<StationaryParticle*>*> G = mergedPolygonGraph(P, Q);
	vector<StationaryParticle*> S = curveBoundary(G);

	map<StationaryParticle*, MovingParticle*> pmap;
	for (int i = 0; i < Q->getNumParticles(); ++i)
	{
		pmap[Q->getParticle(i)->getInitParticle()] = Q->getParticle(i);
	}
	for (int i = 0; i < P->getNumParticles(); ++i)
	{
		pmap[P->getParticle(i)->getInitParticle()] = P->getParticle(i);
	}

	vector<MovingParticle*> vp;
	for (int i = 0; i < S.size(); ++i)
	{
		vp.push_back(pmap[S[i]]);
	}
	GraphFactory<StationaryParticle*>& gfactory = GraphFactory<StationaryParticle*>::GetInstance();
	gfactory.Clean();

	PolygonFactory& pfactory = PolygonFactory::getInstance();
	return pfactory.makePolygon(vp, 0.0f);
}

float 
fitnessMeasure(vector<MovingParticle*>& vp, float scale)
{
	float sum = 0, sum2 = 0;
	for (int i = 0; i < vp.size(); ++i)
	{
		MovingParticle* p = vp[i];
		MovingParticle* q = vp[(i + 1) % vp.size()];
		MovingParticle* r = vp[(i - 1 + vp.size()) % vp.size()];
		CParticleF p0 = p->getP();
		CParticleF q0 = q->getP();
		CParticleF r0 = r->getP();
		if (vp[i]->getNext() != q)
		{
			float d = Distance(p0, q0);
			sum += d * d;
		}
		if (vp[i]->getNext() != q || vp[i]->getPrev() != r)
		{
			float dx1 = q0.m_X - p0.m_X;
			float dx2 = p0.m_X - r0.m_X;
			float dy1 = q0.m_Y - p0.m_Y;
			float dy2 = p0.m_Y - r0.m_Y;
			sum2 += (dx1 - dx2)*(dx1 - dx2) + (dy1 - dy2)*(dy1 - dy2);
		}
	}
	return 1.0 / (0.01 + sum + scale * sum2);
}

float
fitnessMeasure0(vector<CParticleF>& vp, float scale)
{
	float area = polygonArea(vp);
	float area2 = polygonArea(ConvexHull2D(vp));
	float maxlen = 0;
	float sum = 0, sum2 = 0;
	for (int i = 0; i < vp.size(); ++i)
	{
		CParticleF p = vp[i];
		CParticleF q = vp[(i + 1) % vp.size()];
		//CParticleF r = vp[(i - 1 + vp.size()) % vp.size()];
		float d = Distance(p, q);
		maxlen = Max(d, maxlen);
		sum += d;
		sum2 += d * d;
	}
	//return sqrt(area) / maxlen;
	//return sqrt(area) / (scale + sum);
	//return area * area / (area2 * sum2);
	//return sqrt(area) / scale;*/
	return area / ((scale + sqrt(sum2)) * sqrt(area2));
}


vector<Snapshot>
clusterPolygons(vector<Snapshot> regions)
{
	struct PolygonPair {
		Polygon* P;
		Polygon* Q;
		Polygon* R;
		float fitness;
		bool operator <(PolygonPair& P)
		{
			return fitness < P.fitness;
		}
	};
	set<Polygon*> polygons;
	vector<PolygonPair> pairs;
	for (int i = 0; i < regions.size(); ++i)
	{
		polygons.insert(regions[i].getPolygon());
		for (int j = i + 1; j < regions.size(); ++j)
		{
			Polygon* m = mergePolygons(regions[i].getPolygon(), regions[j].getPolygon());
			if (m != NULL)
			{
				PolygonPair pa;
				pa.P = regions[i].getPolygon();
				pa.Q = regions[j].getPolygon();
				pa.R = m;
				//vector<CParticleF> mp = MovingParticles2CParticleF(m->getParticles());
				pa.fitness = fitnessMeasure(m->getParticles(), 1.0);
				pairs.push_back(pa);
				printf("mergePolygons: [%d, %d], [%d, %d] => [%d, %d] %f\n",
					pa.P->getId(), pa.P->getNumParticles(), 
					pa.Q->getId(), pa.Q->getNumParticles(),
					pa.R->getId(), pa.R->getNumParticles(),
					pa.fitness);
			}
		}
	}
	sort(pairs.begin(), pairs.end());

	vector<Snapshot> clusters;
	while (pairs.size()>0 && (pairs.end() - 1)->fitness > 0)
	{
		Polygon* P = (pairs.end() - 1)->P;
		Polygon* Q = (pairs.end() - 1)->Q;
		Polygon* R = (pairs.end() - 1)->R;
		/*if (P->getId() == 99 && Q->getId()==412)
		{
			P->print();
			Q->print();
			R->print();
			mexErrMsgTxt("for debugging...");
		}*/
		Snapshot shot(0.0, 0.0, R->getParticles());
		clusters.push_back(shot);
		polygons.insert(R);

		//erase pairs that involve P and Q
		for (int i = pairs.size() - 1; i >= 0; i--)
		{
			if (pairs[i].P->getId() == P->getId() || pairs[i].P->getId() == Q->getId() ||
				pairs[i].Q->getId() == P->getId() || pairs[i].Q->getId() == Q->getId())
			{
				pairs.erase(pairs.begin() + i);
			}
		}
		polygons.erase(polygons.find(P));
		polygons.erase(polygons.find(Q));
		for (set<Polygon*>::iterator it = polygons.begin(); it != polygons.end(); ++it)
		{
			Polygon* S = *it;
			Polygon* m = mergePolygons(R, S);
			if (m != NULL)
			{
				PolygonPair pa;
				pa.P = R;
				pa.Q = S;
				pa.R = m;
				vector<CParticleF> mp = MovingParticles2CParticleF(m->getParticles());
				pa.fitness = fitnessMeasure(m->getParticles(), 1.0f);
				pairs.push_back(pa);
				printf("mergePolygons: [%d, %d], [%d, %d] => [%d, %d] %f\n",
					pa.P->getId(), pa.P->getNumParticles(),
					pa.Q->getId(), pa.Q->getNumParticles(),
					pa.R->getId(), pa.R->getNumParticles(),
					pa.fitness);
			}
		}
		sort(pairs.begin(), pairs.end());
	}
	return clusters;
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
	int maxLevel = 10;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(maxLevel, prhs[2], classMode);
	}
	float thres = 0.1f;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[3], classMode);
	}
	int minLength = 3;
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(minLength, prhs[4], classMode);
	}

	simulator.Simulate(maxLevel, thres, minLength);

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
		vector<Snapshot> clusters = clusterPolygons(simulator.closedRegions);
		plhs[4] = Snapshot::StoreSnapshots(clusters);
		//plhs[4] = Snapshot::StoreSnapshots(simulator.traces);
	}
	if (nlhs >= 6)
	{
		plhs[5] = simulator.SaveDoneEvents();
	}

	ParticleFactory::getInstance().clean();
	PolygonFactory::getInstance().clean();
	mexUnlock();
}

