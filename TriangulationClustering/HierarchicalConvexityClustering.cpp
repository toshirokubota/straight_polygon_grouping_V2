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
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <Triangulation.h>
#include <FragmentInfo.h>
#include <TriangulationHelper.h>
#include <StationaryParticle.h>
#include <MovingParticle.h>
#include <Polygon.h>
#include <Snapshot.h>
#include <GraphFactory.h>

int MovingParticle::_id = 0;
int Polygon::_id = 0;

bool
sanityCheck(Vertex<Polygon*>* u)
{
	Polygon* P = u->key;
	set<StationaryParticle*> sP = P->getStationaryParticleSet();
	for (int i = 0; i < u->aList.size(); ++i)
	{
		Polygon* Q = u->aList[i]->v->key;
		set<StationaryParticle*> sQ = Q->getStationaryParticleSet();
		set<StationaryParticle*> iset;
		set_intersection(sP.begin(), sP.end(),
			sQ.begin(), sQ.end(),
			std::inserter(iset, iset.begin()));
		if (iset.size() < 2)
		{
			printf("sanityCheck:\n");
			P->print();
			Q->print();
			printf("sanityCheck: no common edges");
			return false;
		}
		/*if (iset.size() >= sQ.size() || iset.size() >= sP.size())
		{
			printf("sanityCheck:\n");
			P->print();
			Q->print();
			printf("sanityCheck: one consumes the other.");
			return false;
		}*/
	}
	return true;
}


vector<StationaryParticle*>
mergedTrace(Polygon* P, Polygon* Q)
{
	set<StationaryParticle*> sP = P->getStationaryParticleSet();
	set<StationaryParticle*> sQ = Q->getStationaryParticleSet();
	set<StationaryParticle*> iset;
	set_intersection(sP.begin(), sP.end(),
		sQ.begin(), sQ.end(),
		std::inserter(iset, iset.begin()));
	if (iset.size() < 2)
	{
		printf("Failed to merge.\n");
		P->print();
		Q->print();
		mexErrMsgTxt("mergedTrace: Two polygons do not share a side.");
	}
	MovingParticle* p1 = NULL;
	MovingParticle* q1 = NULL;
	{
		vector<MovingParticle*> vp = P->getParticles();
		for (int i = 0; i < vp.size(); ++i)
		{
			StationaryParticle* p0 = vp[i]->getInitParticle();
			StationaryParticle* q0 = vp[i]->getNext()->getInitParticle();
			if (iset.find(p0) == iset.end() && iset.find(q0) != iset.end())
			{
				q1 = vp[i];
			}
			else if (iset.find(p0) != iset.end() && iset.find(q0) == iset.end())
			{
				p1 = vp[i];
			}
		}
	}
	MovingParticle* p2 = NULL;
	MovingParticle* q2 = NULL;
	{
		vector<MovingParticle*> vp = Q->getParticles();
		for (int i = 0; i < vp.size(); ++i)
		{
			StationaryParticle* p0 = vp[i]->getInitParticle();
			StationaryParticle* q0 = vp[i]->getNext()->getInitParticle();
			if (iset.find(p0) == iset.end() && iset.find(q0) != iset.end())
			{
				q2 = vp[i];
			}
			else if (iset.find(p0) != iset.end() && iset.find(q0) == iset.end())
			{
				p2 = vp[i];
			}
		}
	}

	vector<StationaryParticle*> vs;
	vector<MovingParticle*> tr;
	if (p1 != NULL && p2 != NULL)
	{
		if (p2 != NULL && p2 != NULL)
		{
			tr = MovingParticle::extractPath(p1, q1->getNext());
			vector<MovingParticle*> tr2 = MovingParticle::extractPath(p2, q2->getNext());
			tr.insert(tr.end(), tr2.begin(), tr2.end());
		}
		else
		{
			tr = MovingParticle::extractPath(p1, q1->getNext()->getNext());
		}
	}
	else
	{
		tr = MovingParticle::extractPath(p2, q2->getNext()->getNext());
	}
	for (int i = 0; i < tr.size(); ++i)
	{
		vs.push_back(tr[i]->getInitParticle());
	}

	return vs;
}

float computeFitness(vector<StationaryParticle*>& trace, float scale = 1.0f)
{
	vector<CParticleF> vp;
	for (int i = 0; i < trace.size(); ++i)
	{
		vp.push_back(trace[i]->getP());
	}
	float area = polygonArea(vp);
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
	return sqrt(area) / (scale + sum);
	//return area / (sum2);
	//return sqrt(area) / scale;
}

Polygon*
makePolygon(vector<StationaryParticle*>& pnts)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<MovingParticle*> vp;
	for (int j = 0; j < pnts.size(); ++j)
	{
		vp.push_back(factory.makeParticle(pnts[j], Initial, 0.0f));
	}
	for (int j = 0; j < vp.size(); ++j)
	{
		MovingParticle::setNeighbors(vp[j], vp[(j - 1 + vp.size()) % vp.size()], vp[(j + 1) % vp.size()]);
	}
	for (int j = 0; j < vp.size(); ++j)
	{
		vp[j]->clearVelocity();
	}
	MovingParticle::updatePolygon(vp[0]);
	return vp[0]->getPolygon();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "TriangleClustering", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: T = TriangleCLustering(P, F)");
		return;
	}

	//Points
	vector<CParticleF> points;
	vector<StationaryParticle*> spnts;
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	ParticleFactory& factory = ParticleFactory::getInstance();
	PolygonFactory& polyfactory = PolygonFactory::getInstance();
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		points = vector2particle(P0, dimsP);
	}
	float thres = .1;	
	if(nrhs>=2) 
	{
		mxClassID classMode;
		ReadScalar(thres,prhs[1],classMode);
	} 

	Triangulation::Triangulator trmap(points);
	map<Triangulation::_Internal::_vertex*, int> pmap;
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		pmap[trmap.points[i]] = i;
		spnts.push_back(sfactory.makeParticle(trmap.points[i]->p));
	}

	map<Triangulation::_Internal::_edge*, Polygon*> emap1;
	map<Triangulation::_Internal::_edge*, Polygon*> emap2;
	for (int i = 0; i < trmap.edges.size(); ++i)
	{
		emap1[trmap.edges[i]] = NULL;
		emap2[trmap.edges[i]] = NULL;
	}

	vector<Polygon*> polygons;
	for (int i = 0; i < trmap.faces.size(); ++i)
	{
		vector<StationaryParticle*> sp;
		for (int j = 0; j < 3; ++j)
		{
			int k = pmap[trmap.faces[i]->vertices[j]];
			sp.push_back(spnts[k]);
		}
		Polygon* poly = makePolygon(sp);
		polygons.push_back(poly);
		for (int j = 0; j < 3; ++j)
		{
			Triangulation::_Internal::_edge* e = trmap.faces[i]->edges[j];
			if (emap1[e] == NULL)
			{
				emap1[e] = poly;
			}
			else
			{
				emap2[e] = poly;
			}
		}
	}

	//finally construct a graph (polygons are vertices and adjacency relations give edges)
	GraphFactory<Polygon*>& gfactory = GraphFactory<Polygon*>::GetInstance();
	vector<Vertex<Polygon*>*>  vertices;
	map<Polygon*, int> polymap;
	for (int i = 0; i < polygons.size(); ++i)
	{
		polymap[polygons[i]] = i;
		Vertex<Polygon*>* v = gfactory.makeVertex(polygons[i]);
		vertices.push_back(v);
	}
	for (int i = 0; i < trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e = trmap.edges[i];
		if (emap1[e] != NULL && emap2[e] != NULL)
		{
			Polygon* p1 = emap1[e];
			Polygon* p2 = emap2[e];
			int k1 = polymap[p1];
			int k2 = polymap[p2];
			Edge<Polygon*>* ed1 = gfactory.makeEdge(vertices[k1], vertices[k2], 1.0f);
			vertices[k1]->Add(ed1);
			Edge<Polygon*>* ed2 = gfactory.makeEdge(vertices[k2], vertices[k1], 1.0f);
			vertices[k2]->Add(ed2);
		}
	}
	vector<Snapshot> shots;
	while (true)
	{
		float bestval = -std::numeric_limits<float>::infinity();
		Edge<Polygon*>* selected = NULL;
		for (int i = 0; i < vertices.size(); ++i)
		{
			for (int j = 0; j < vertices[i]->aList.size(); ++j)
			{
				Edge<Polygon*>* ed = vertices[i]->aList[j];
				if (ed->u->key->getId() > ed->v->key->getId()) continue; //no need to do both because symetric

				vector<StationaryParticle*> pnts = mergedTrace(ed->u->key, ed->v->key);
				float fit = computeFitness(pnts);
				if (fit > bestval)
				{
					bestval = fit;
					selected = ed;
				}
			}
		}
		if (bestval > thres)
		{
			//create a new polygon and its vertex
			vector<StationaryParticle*> vp = mergedTrace(selected->u->key, selected->v->key);
			Polygon* poly = makePolygon(vp);

			printf("Selected 1:\n");
			selected->u->key->print();
			printf("Selected 2:\n");
			selected->v->key->print();
			printf("Merged:\n");
			poly->print();

			Vertex<Polygon*>* w = gfactory.makeVertex(poly);
			shots.push_back(Snapshot(0.0f, 0.0f, poly->getParticles()));

			//find vertices adjacent to the merged vertices
			set<Vertex<Polygon*>*> adj;
			{
				Vertex<Polygon*>* u = selected->u;
				Vertex<Polygon*>* v = selected->v;
				for (int j = 0; j < u->aList.size(); ++j)
				{
					if (u->aList[j]->v == v) continue;
					adj.insert(u->aList[j]->v);
				}
				for (int j = 0; j < v->aList.size(); ++j)
				{
					if (v->aList[j]->v == u) continue;
					adj.insert(v->aList[j]->v);
				}
			}

			//add new edges 
			for (set<Vertex<Polygon*>*>::iterator it = adj.begin(); it != adj.end(); ++it)
			{
				Vertex<Polygon*>* v = *it;
				w->Add(gfactory.makeEdge(w, v, 1.0f));
				v->Add(gfactory.makeEdge(v, w, 1.0f));
			}
			//remove edges from adjacent ones to the merged ones
			for (set<Vertex<Polygon*>*>::iterator it = adj.begin(); it != adj.end(); ++it)
			{
				while ((*it)->Remove(selected->u));
				while ((*it)->Remove(selected->v));
			}
			vertices.erase(find(vertices.begin(), vertices.end(), selected->u));
			vertices.erase(find(vertices.begin(), vertices.end(), selected->v));
			vertices.push_back(w);
			if (sanityCheck(w) == false) break;
		}
		else
		{
			break;
		}
	}

	if(nlhs >= 1)
	{
		plhs[0] = Snapshot::StoreSnapshots(shots);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { trmap.edges.size(), 2 };
		vector<float> F(dims[0] * dims[1], 0);
		for (int i = 0; i<trmap.edges.size(); ++i)
		{
			int j1 = pmap[trmap.edges[i]->vertices[0]];
			int j2 = pmap[trmap.edges[i]->vertices[1]];
			SetData2(F, i, 0, dims[0], dims[1], (float)j1 + 1);
			SetData2(F, i, 1, dims[0], dims[1], (float)j2 + 1);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

