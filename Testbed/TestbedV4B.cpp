#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>
#include <mex.h>
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
#include <StationaryParticle.h>
#include <DisjointSet.h>
#include <Triangulation.h>
#include <IntersectionConvexPolygons.h>

struct MAPoint
{
	MAPoint(Triangulation::_Internal::_vertex* p, Triangulation::_Internal::_vertex* q, int pid, int qid)
	{
		this->x = (p->p.m_X + q->p.m_X) / 2.0;
		this->y = (p->p.m_Y	+ q->p.m_Y) / 2.0;
		this->z = Distance(p->p, q->p) / 2.0;
		this->theta = GetVisualDirection(q->p.m_X, q->p.m_X, p->p.m_X, p->p.m_Y);
		this->p = p;
		this->q = q;
		this->pid = pid;
		this->qid = qid;
	}
	void print(char* tab="", char* end="\n")
	{
		printf("%s%d (%3.3f,%3.3f) %d (%3.3f,%3.3f) %3.3f  %3.3f %3.3f %3.3f%s",
			tab, pid, p->p.m_X, p->p.m_Y, qid, q->p.m_X, q->p.m_Y, x, y, z, theta, end);
	}
	float x;
	float y;
	float z;
	float theta;
	Triangulation::_Internal::_vertex* p;
	Triangulation::_Internal::_vertex* q;
	int pid, qid; //for debugging
};

/*
For each vertes in a triangulation map, find the 2nd shortest edge.
Among the 2nd shortest edges, return the length of the longest one.
After we remove edges that are longer than this return value, each vertex still retains two edges in the triangulation map.
*/
float calculateScaleUnit(Triangulation::Triangulator& trmap)
{
	float separation = 0;
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* p = trmap.points[i];
		vector<float> vlen;
		for (int j = 0; j < p->edges.size(); ++j)
		{
			vlen.push_back(p->edges[j]->Length());
		}
		sort(vlen.begin(), vlen.end());
		if (vlen[1]> separation)
		{
			separation = vlen[1];
		}
	}
	printf("Separation = %f\n", separation);
	return separation;
}

/*float distance(MAPoint* p, MAPoint* q, float scale=1.0f, float  wz = 1.0f, float wt = 30.0f,  bool bprint=false)
{
	float dx = (p->x - q->x)/scale;
	float dy = (p->y - q->y)/scale;
	float dz = wz * (p->z - q->z)/scale;
	float dt = wt * sin(p->theta - q->theta);
	float d = sqrt(dx*dx + dy*dy + dz*dz + dt*dt);
	if (bprint)
	{
		printf("dx=%3.3f dy=%3.3f dz=%3.3f dt=%3.3f d=%3.3f\n",
			dx, dy, dz, dt, d);
	}
	return d;
}*/

float distance(MAPoint* p, MAPoint* q, float scale=1.0f, float wt = 1.0f,  bool bprint=false)
{
	pair<float, float> param = _IntersectConvexPolygon::intersect(p->p->p, p->q->p, q->p->p, q->q->p);
	if ((param.first > 0 && param.first <1.0) && (param.second > 0 && param.second < 1.0))
	{
		return std::numeric_limits<float>::infinity();
	}
	float dx = q->x - p->x;
	float dy = q->y - p->y;
	float dz = q->z - p->z;
	float u = -sin(p->theta);
	float v = cos(p->theta);
	float w = 0;
	float t = u*dx + v*dy;
	float dd = ((dx - u*t)*(dx - u*t) + (dy - v*t)*(dy - v*t) + dz*dz)/(scale * scale);
	float dt = Abs(asin(sin(p->theta - q->theta)));
	float d = sqrt(dd + wt * dt * dt);
	if (bprint)
	{
		printf("dd=%3.3f dt=%3.3f d=%3.3f\n", dd, dt, d);
	}
	return d;
}


vector<MAPoint*>
collectMedialAxisPoints(Triangulation::Triangulator& trmap)
{
	vector<MAPoint*> mapoints;
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		for (int j = i + 1; j < trmap.points.size(); ++j)
		{
			mapoints.push_back(new MAPoint(trmap.points[i], trmap.points[j], i, j));
		}
	}
	return mapoints;
}

bool
isNeighbor(Triangulation::_Internal::_vertex* p, Triangulation::_Internal::_vertex* q)
{
	if (p == q)
	{
		return false;
	}

	for (int i = 0; i < p->edges.size(); ++i)
	{
		if (q == p->edges[i]->vertices[0] || q == p->edges[i]->vertices[1])
		{
			return true;
		}
	}
	return false;
}

/*
Collect a set of MAPoints that are neighbors to the given MAPoint.
An MAPoint y is a neighbor of x if y->p is a neighbor of x->p AND y->q is a neighbor of x->q OR
y->p is a neighbor of x->q AND y->q is a neighbor of x->p in the triangulation map.
*/
set<MAPoint*>
collectNeighbors(MAPoint* p, vector<MAPoint*>& points)
{
	set<MAPoint*> pset;
	for (int i = 0; i < points.size(); ++i)
	{
		//if two MAPoints share a vertex, they are not neighbors
		if (p->p == points[i]->p || p->p == points[i]->q) continue;
		if (p->q == points[i]->p || p->q == points[i]->q) continue;

		if ((isNeighbor(p->p, points[i]->p) && isNeighbor(p->q, points[i]->q)) ||
			(isNeighbor(p->p, points[i]->q) && isNeighbor(p->q, points[i]->p)))
		{
			pset.insert(points[i]);
		}
	}
	return pset;
}

vector<int>
clusterMAPoints(vector<MAPoint*>& P, float thres,  float scale, float weight)
{
	vector<Node<MAPoint*>*> nodes;
	map<MAPoint*, int> imap;
	for (int i = 0; i < P.size(); ++i)
	{
		nodes.push_back(makeset(P[i]));
		imap[P[i]] = i;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		MAPoint* p = nodes[i]->key;
		//collect possible neighbors
		bool bprint = false; // p->pid == 4 && p->qid == 12 || p->pid == 4 && p->qid == 28;
		//bool bprint = p->pid == 13 && p->qid == 14;
		if (bprint)
		{
			p->print();
		}
		set<MAPoint*> candidates = collectNeighbors(p, P);
		for (set<MAPoint*>::iterator it = candidates.begin(); it != candidates.end(); ++it)
		{
			MAPoint* q = *it;
			if (bprint)
			{
				q->print("\t", ": ");
			}
			Node<MAPoint*>* n = nodes[imap[q]];
			float d = distance(p, q, scale, weight, bprint);
			if(d < thres)
			{
				merge(nodes[i], n);
				printf("%d %d,  %d %d  %3.3f\n", p->pid, p->qid, q->pid, q->qid, d);
			}
		}
	}
	vector<Node<MAPoint*>*> reps = clusters(nodes);
	map<Node<MAPoint*>*, int> imap2;
	for (int i = 0; i < reps.size(); ++i)
	{
		imap2[reps[i]] = i;
	}
	vector<int> labels;
	for (int i = 0; i < nodes.size(); ++i)
	{
		labels.push_back(imap2[findset(nodes[i])]);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return labels;
}

vector<vector<Triangulation::_Internal::_vertex*> >
groupPoints(vector<MAPoint*>& points, vector<int> labels)
{
	set<int> iset;
	map<int, int> imap;
	int count = 0;
	for (int i = 0; i < labels.size(); ++i)
	{
		if (iset.find(labels[i]) == iset.end())
		{
			imap[labels[i]] = count;
			count++;
			iset.insert(labels[i]);
		}
	}
	vector<set<Triangulation::_Internal::_vertex*>> G(iset.size());
	for (int i = 0; i < points.size(); ++i)
	{
		int lb = labels[i];
		int k = imap[lb];
		G[k].insert(points[i]->p);
		G[k].insert(points[i]->q);
	}
	vector<vector<Triangulation::_Internal::_vertex*>> vG(iset.size());
	for (int i = 0; i < G.size(); ++i)
	{
		for (set<Triangulation::_Internal::_vertex*>::iterator it = G[i].begin(); it != G[i].end(); it++)
		{
			vG[i].push_back(*it);
		}
	}
	return vG;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = Testbed(P, [sigma, step])");
		return;
	}
	//Points
	vector<CParticleF> P;
	const int* dimsP;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for (int i = 0; i < dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			P.push_back(CParticleF(x, y));
		}
	}

	float thres = 5.0f;
	if (nrhs >= 2)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(thres, prhs[1], classId);
	}
	float weight = 10.0f;
	if (nrhs >= 3)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(weight, prhs[2], classId);
	}
	Triangulation::Triangulator trmap(P);
	float unit = calculateScaleUnit(trmap); //set the unit
	vector<MAPoint*> mapoints = collectMedialAxisPoints(trmap);
	vector<int> labels = clusterMAPoints(mapoints, thres, unit, weight); 
	vector<vector<Triangulation::_Internal::_vertex*>> G = groupPoints(mapoints, labels);

	if (nlhs >= 1)
	{
		const int dims[] = { G.size(), 1 };
		map<Triangulation::_Internal::_vertex*, int> imap;
		for (int i = 0; i < trmap.points.size(); ++i)
		{
			imap[trmap.points[i]] = i;
		}
		vector<vector<int>> vvint;
		for (int i = 0; i < dims[0]; ++i)
		{
			const int dims2[] = { G[i].size(), 1 };
			vector<int> vint(dims2[0] * dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				SetData2(vint, j, 0, dims2[0], dims2[1], imap[G[i][j]]+1);
			}
			vvint.push_back(vint);
		}
		plhs[0] = StoreDataCell(vvint, mxINT32_CLASS, 2, dims, 1);
	}
	if (nlhs >= 2)
	{
		int cnt = 0;
		for (int i = 0; i < labels.size(); ++i)
		{
			cnt = Max(cnt, labels[i]);
		}
		const int dims[] = { cnt+1, 1 };
		vector<vector<int>> vvint(dims[0]);
		for (int i = 0; i < labels.size(); ++i)
		{
			MAPoint* p = mapoints[i];
			int k = labels[i];
			vvint[k].push_back(p->pid);
			vvint[k].push_back(p->qid);
		}
		vector<vector<int>> vvint2(dims[0]);
		for (int i = 0; i < vvint.size(); ++i)
		{
			for (int j = 0; j < vvint[i].size() / 2; ++j)
			{
				vvint2[i].push_back(vvint[i][j]);
				vvint2[i].push_back(vvint[i][j + vvint[i].size() / 2]);
			}
		}
		plhs[1] = StoreDataCell(vvint2, mxINT32_CLASS, 2, dims, 2);
	}
	for (int i = 0; i < mapoints.size(); ++i)
	{
		delete mapoints[i];
	}
	/*for (int i = 0; i < trmap.points.size(); ++i)
	{
		links[trmap.points[i]]->clear();
		delete links[trmap.points[i]];
		links[trmap.points[i]] = NULL;
	}*/
	StationaryParticleFactory::getInstance().clean();
	mexUnlock();
}

