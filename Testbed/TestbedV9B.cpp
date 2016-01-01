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
#include <MiscGeometry.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <Dijkstra.h>
#include <Kruskal.h>

struct MAPoint
{
	MAPoint(Triangulation::_Internal::_vertex* p, Triangulation::_Internal::_vertex* q, int pid, int qid)
	{
		this->x = (p->p.m_X + q->p.m_X) / 2.0;
		this->y = (p->p.m_Y	+ q->p.m_Y) / 2.0;
		this->z = Distance(p->p, q->p) / 2.0; 
		//this->z = sqrt(this->z);
		this->theta = GetVisualDirection(q->p.m_X, q->p.m_Y, p->p.m_X, p->p.m_Y);
		this->p = p;
		this->q = q;
		this->pid = pid;
		this->qid = qid;
		this->id = _id++;
		this->x0 = this->x;
		this->y0 = this->y;
		this->z0 = this->z;
		//next = NULL;
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
	int id, pid, qid; //for debugging
	//temporal location
	float x0;
	float y0;
	float z0;
	set<MAPoint*> neighbors;
	//MAPoint* avater; //this allows to move to another point to  seek for a higher ground.
	//MAPoint* next; //used to build a DAG
	static int _id;
};

int MAPoint::_id = 0;

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

float distance0(MAPoint* p, MAPoint* q, float scale = 1.0f, float  wz = 1.0f, float wt = 30.0f, bool bprint = false)
{
	float dx = (p->x - q->x) / scale;
	float dy = (p->y - q->y) / scale;
	float dz = wz * (p->z - q->z) / scale;
	float dt = wt * sin(p->theta - q->theta);
	float d = sqrt(dx*dx + dy*dy + dz*dz + dt*dt);
	if (bprint)
	{
		printf("dx=%3.3f dy=%3.3f dz=%3.3f dt=%3.3f d=%3.3f\n",
			dx, dy, dz, dt, d);
	}
	return d;
}

float distance1(MAPoint* p, MAPoint* q, float thres = 0.9)
{
	if (Abs(cos(p->theta - q->theta)) < thres)
	{
		return std::numeric_limits<float>::infinity();
	}
	float dx = (p->x - q->x);
	float dy = (p->y - q->y);
	float dz = (p->z - q->z);
	//float dt = wt * sin(p->theta - q->theta);
	float d = sqrt(dx*dx + dy*dy + dz*dz);
	return d;
}

/*
compute the closest point from x on the line connecting p and q.
p and q have to be distinct.
This is done in 3D.
*/
CParticleF Closest2Line3d( CParticleF& p, CParticleF& q, CParticleF& x)
{
	CParticleF u = NormalizedDirection3D(q, p);
	CParticleF xp(x.m_X - p.m_X, x.m_Y - p.m_Y, x.m_Z - p.m_Z);
	float dp = xp.m_X*u.m_X + xp.m_Y*u.m_Y + xp.m_Z*u.m_Z; 
	CParticleF y(p.m_X + dp*u.m_X, p.m_Y + dp*u.m_Y, p.m_Z + dp*u.m_Z);
	return y;
}


float distanceAsymmetric(MAPoint* p, MAPoint* q)
{
	CParticleF p0(p->x, p->y, p->z);
	CParticleF q0(q->x, q->y, p->z); //disregard Z - change this to q->z to regard Z
	//return Distance(p0, q0);
	CParticleF p2(p->x - sin(p->theta), p->y + cos(p->theta), p->z);
	CParticleF y = Closest2Line3d(p0, p2, q0);
	return Distance(y, q0);
}

float distance(MAPoint* p, MAPoint* q)
{
	return distanceAsymmetric(p, q) + distanceAsymmetric(q, p);
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

int
_key(int pid, int qid, int size)
{
	int m = Min(pid, qid);
	int n = Max(pid, qid);
	return size * m + n;
}

/*
Collect a set of MAPoints that are neighbors to the given MAPoint.
An MAPoint y is a neighbor of x if y->p is a neighbor of x->p AND y->q is a neighbor of x->q OR
y->p is a neighbor of x->q A;ND y->q is a neighbor of x->p in the triangulation map.
*/
void
collectNeighbors(vector<MAPoint*>& points, Triangulation::Triangulator& trmap, float thres)
{
	int n = trmap.points.size();
	//builds data structures for efficiency
	vector<vector<int>> N(n);
	for (int i = 0; i < n; ++i)
	{
		Triangulation::_Internal::_vertex* u = trmap.points[i];
		for (int j = 0; j < n; ++j)
		{
			Triangulation::_Internal::_vertex* v = trmap.points[j];
			if (Distance(u->p, v->p) <= thres)
			{
				N[i].push_back(j);
			}
		}
	}
	int maxKey = 0;
	for (int i = 0; i < points.size(); ++i)
	{
		maxKey = Max(maxKey, _key(points[i]->pid, points[i]->qid, n));
	}
	vector<MAPoint*> pmap(maxKey + 1, NULL);
	for (int i = 0; i < points.size(); ++i)
	{
		pmap[_key(points[i]->pid, points[i]->qid, n)] = points[i];
	}
	for (int k = 0; k < points.size(); ++k)
	{
		MAPoint* p = points[k];
		for (int i = 0; i<N[p->pid].size(); ++i)
		{
			for (int j = 0; j < N[p->qid].size(); ++j)
			{
				int key = _key(N[p->pid][i], N[p->qid][j], n);
				if (key <= maxKey)
				{
					MAPoint* q = pmap[key];
					if (q != NULL)
					{
						p->neighbors.insert(q);
					}
				}
			}
		}
	}
}

vector<MAPoint*>
collectMedialAxisPoints(Triangulation::Triangulator& trmap, float thres)
{
	vector<MAPoint*> mapoints;
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		for (int j = i + 1; j < trmap.points.size(); ++j)
		{
			//if (trmap.points[i]->p.m_Life != trmap.points[j]->p.m_Life) continue; //TK!!

			if (Distance(trmap.points[i]->p, trmap.points[j]->p) <= thres) 
			{
				mapoints.push_back(new MAPoint(trmap.points[i], trmap.points[j], i, j));
			}
		}
	}
	return mapoints;
}

/*
*/
vector<Vertex<MAPoint*>*>
buildDag(vector<MAPoint*>& P, float thres)
{
	GraphFactory<MAPoint*>& factory = GraphFactory<MAPoint*>::GetInstance();
	vector<Vertex<MAPoint*>*> vertices(P.size());
	map<MAPoint*, int>  imap;
	for (int i = 0; i < P.size(); ++i)
	{
		vertices[i] = factory.makeVertex(P[i]);
		imap[P[i]] = i;
	}
	vector<Edge<MAPoint*>*> edges;

	for (int i = 0; i < P.size(); ++i)
	{
		MAPoint* p = P[i];
		vector<pair<float, MAPoint*>> pairs;
		for (set<MAPoint*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
		{
			MAPoint* q = *it;
			if (p == q) continue;
			float w = distance(p, q);
			pairs.push_back(pair<float, MAPoint*>(w*w, q));
		}
		if (pairs.empty()) continue;

		sort(pairs.begin(), pairs.end());
		float minW = pairs.size()>1 ? pairs[1].first: pairs[0].first; //should keep two choices at least
		for (int j = 0; j < pairs.size(); ++j)
		{
			if (pairs[j].first > thres) break;
			//if (pairs[j].first > minW * 2.0) break;

			MAPoint* chosen = pairs[j].second;
			if (chosen->z < p->z) continue; //only consider upward direction

			int k = imap[chosen];
			Edge<MAPoint*>* ed = factory.makeEdge(vertices[i], vertices[k], pairs[j].first);
			vertices[i]->Add(ed);
		}
	}
	return vertices;
}

vector<Vertex<MAPoint*>*>
mergeVertices(vector<Vertex<MAPoint*>*>& peaks)
{
	vector<Node<Vertex<MAPoint*>*>*> nodes;
	for (int i = 0; i < peaks.size(); ++i)
	{
		nodes.push_back(makeset(peaks[i]));
	}
	for (int i = 0; i < peaks.size(); ++i)
	{
		MAPoint* p = peaks[i]->key;
		for (int j = 0; j < peaks.size(); ++j)
		{
			if (i == j) continue;
			MAPoint* q = peaks[j]->key;
			if (p->neighbors.find(q) != p->neighbors.end())
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<Node<Vertex<MAPoint*>*>*> reps = clusters(nodes);
	GraphFactory<MAPoint*>& factory = GraphFactory<MAPoint*>::GetInstance();
	vector<Vertex<MAPoint*>*> X;
	for (int i = 0; i < reps.size(); ++i)
	{
		Vertex<MAPoint*>* u = factory.makeVertex(NULL);
		for (int j = 0; j < peaks.size(); ++j)
		{
			if (findset(nodes[j]) == reps[i])
			{
				nodes[j]->key->Add(factory.makeEdge(nodes[j]->key, u));
			}
		}
		X.push_back(u);
	}
	return X;
}

vector<vector<int>>
groupVertices(vector<Vertex<MAPoint*>*>& vertices)
{
	vector<Vertex<MAPoint*>*> peaks0;
	for (int i = 0; i < vertices.size(); ++i)
	{
		Vertex<MAPoint*>* u = vertices[i];
		if (u->aList.empty())
		{
			peaks0.push_back(u);
		}
	}
	vector<Vertex<MAPoint*>*> peaks = mergeVertices(peaks0);
	reverseEdges(vertices);

	vector<vector<int>> labels;
	for (int i = 0; i < peaks.size(); ++i)
	{
		set<Vertex<MAPoint*>*> traced;
		vector<Vertex<MAPoint*>*> Q(1, peaks[i]);
		while (Q.empty() == false)
		{
			vector<Vertex<MAPoint*>*> Q2;
			for (int j = 0; j < Q.size(); ++j)
			{
				Vertex<MAPoint*>* u = Q[j];
				for (int k = 0; k < u->aList.size(); ++k)
				{
					Vertex<MAPoint*>* v = u->aList[k]->v;
					if (traced.find(v) == traced.end())
					{
						Q2.push_back(v);
						traced.insert(v);
					}
				}
			}
			Q = Q2;
		}
		vector<int> L;
		for (set<Vertex<MAPoint*>*>::iterator it = traced.begin(); it != traced.end(); ++it)
		{
			Vertex<MAPoint*>* u = *it;
			if (u->key == NULL) continue; //artificial super peak vertex

			L.push_back(u->key->id);
		}
		labels.push_back(L);
	}
	
	GraphFactory<MAPoint*>& factory = GraphFactory<MAPoint*>::GetInstance();
	for (int i = 0; i < peaks.size(); ++i)
	{
		factory.Clean(peaks[i]);
	}
	reverseEdges(vertices);
	return labels;
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
			float val = GetData2(P0, i, 2, dimsP[0], dimsP[1], (float)0);
			P.push_back(CParticleF(x, y, 0, val));
		}
	}

	float thres = std::numeric_limits<float>::infinity();
	if (nrhs >= 2)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(thres, prhs[1], classId);
	}
	float maxDist = std::numeric_limits<float>::infinity();
	if (nrhs >= 3)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(maxDist, prhs[2], classId);
	}
	MAPoint::_id = 0;

	Triangulation::Triangulator trmap(P);
	float unit = calculateScaleUnit(trmap); //set the unit

	vector<MAPoint*> mapoints = collectMedialAxisPoints(trmap, maxDist);
	collectNeighbors(mapoints, trmap, unit*1.25);

	GraphFactory<MAPoint*>& factory = GraphFactory<MAPoint*>::GetInstance();
	vector<Vertex<MAPoint*>*> vertices = buildDag(mapoints, thres);
	vector<vector<int>> labels = groupVertices(vertices);

	if (nlhs >= 1)
	{
		const int dims[] = { labels.size(), 1 };
		vector<vector<int>> vvint;
		for (int i = 0; i < dims[0]; ++i)
		{
			const int dims2[] = { labels[i].size(), 1 };
			vector<int> vint(dims2[0] * dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				SetData2(vint, j, 0, dims2[0], dims2[1], labels[i][j]);
			}
			vvint.push_back(vint);
		}
		plhs[0] = StoreDataCell(vvint, mxINT32_CLASS, 2, dims, 1);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { factory.edges.size(), 7 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			Edge<MAPoint*>* ed = factory.edges[i];
			SetData2(F, i, 0, dims[0], dims[1], (float)ed->u->key->pid);
			SetData2(F, i, 1, dims[0], dims[1], (float)ed->u->key->qid);
			SetData2(F, i, 2, dims[0], dims[1], (float)ed->v->key->pid);
			SetData2(F, i, 3, dims[0], dims[1], (float)ed->v->key->qid);
			SetData2(F, i, 4, dims[0], dims[1], (float)ed->u->key->id);
			SetData2(F, i, 5, dims[0], dims[1], (float)ed->v->key->id);
			SetData2(F, i, 6, dims[0], dims[1], (float)ed->w);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		const int dims[] = { mapoints.size(), 7};
		vector<float> F(dims[0]*dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], (float)mapoints[i]->id);
			SetData2(F, i, 1, dims[0], dims[1], mapoints[i]->x);
			SetData2(F, i, 2, dims[0], dims[1], mapoints[i]->y);
			SetData2(F, i, 3, dims[0], dims[1], mapoints[i]->z);
			SetData2(F, i, 4, dims[0], dims[1], mapoints[i]->theta);
			SetData2(F, i, 5, dims[0], dims[1], (float)mapoints[i]->pid);
			SetData2(F, i, 6, dims[0], dims[1], (float)mapoints[i]->qid);
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	for (int i = 0; i < mapoints.size(); ++i)
	{
		delete mapoints[i];
	}
	GraphFactory<MAPoint*>::GetInstance().Clean();

	mexUnlock();
}

