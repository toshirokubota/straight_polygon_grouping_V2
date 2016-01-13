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

float distanceSymmetric(MAPoint* p, MAPoint* q)
{
	return distanceAsymmetric(p, q) + distanceAsymmetric(q, p);
}

float distance2D(MAPoint* p, MAPoint* q)
{
	return Distance(p->x, p->y, q->x, q->y);
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
	map<Triangulation::_Internal::_vertex*, int> vmap;
	for (int i = 0; i < n; ++i)
	{
		vmap[trmap.points[i]] = i;
	}
	vector<vector<int>> N(n);
	for (int i = 0; i < n; ++i)
	{
		Triangulation::_Internal::_vertex* u = trmap.points[i];
		/*for (int j = 0; j < n; ++j)
		{
		Triangulation::_Internal::_vertex* v = trmap.points[j];
		if (Distance(u->p, v->p) <= thres)
		{
		N[i].push_back(j);
		}
		}*/
		N[i].push_back(i); //push itself
		for (int j = 0; j < u->edges.size(); ++j)
		{
			Triangulation::_Internal::_vertex* v = u->edges[j]->vertices[0] == u ? u->edges[j]->vertices[1] : u->edges[j]->vertices[0];
			if (Distance(u->p, v->p) <= thres)
			{
				N[i].push_back(vmap[v]);
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
			float w = distanceSymmetric(p, q);
			if (w < thres)
			{
				pairs.push_back(pair<float, MAPoint*>(w*w, q));
			}
		}
		if (pairs.empty()) continue;

		//sort(pairs.begin(), pairs.end());
		//float minW = pairs.size()>1 ? pairs[1].first: pairs[0].first; //should keep two choices at least
		for (int j = 0; j < pairs.size(); ++j)
		{
			//if (pairs[j].first > thres) break;
			//if (pairs[j].first > minW * 2.0) break;

			MAPoint* chosen = pairs[j].second;
			if (chosen->z < p->z) continue; //only consider upward direction
			if (distance2D(p, chosen) > p->z) continue; //has to be within the radius of convexity

			int k = imap[chosen];
			Edge<MAPoint*>* ed = factory.makeEdge(vertices[i], vertices[k], pairs[j].first);
			vertices[i]->Add(ed);
		}
	}
	return vertices;
}
struct TimedPoint
{
	TimedPoint(MAPoint* p, MAPoint* q, int t)
	{
		this->p = p;
		this->q = q;
		this->time = t;
	}
	MAPoint* p;
	MAPoint* q;
	int time;
};

vector<TimedPoint>
groupVertices(vector<Vertex<MAPoint*>*>& vertices, vector<Vertex<MAPoint*>*>& peaks)
{
	reverseEdges(vertices);

	for (int i = 0; i < vertices.size(); ++i)
	{
		vertices[i]->Reset();
	}

	int iter = 0;
	vector<Vertex<MAPoint*>*> Q = peaks;
	vector<TimedPoint> G;

	for (int i = 0; i < peaks.size(); ++i)
	{
		G.push_back(TimedPoint(peaks[i]->key, NULL, iter));
		peaks[i]->color = Black;
	}
	while (Q.empty() == false)
	{
		iter++;
		vector<Vertex<MAPoint*>*> Q2;
		for (int j = 0; j < Q.size(); ++j)
		{
			Vertex<MAPoint*>* u = Q[j];
			if (u->key->id == 74840)
			{
				MAPoint* p = u->key;
				printf("%d %d %d %f %s\n", p->id, p->pid, p->qid, p->z, u->color == White ? "White" : "Black");
			}
			for (int k = 0; k < u->aList.size(); ++k)
			{
				Vertex<MAPoint*>* v = u->aList[k]->v;
				if (v->color == White)
				{
					bool bOk = true;
					MAPoint* q = v->key;
					if (u->key->id == 74840)
					{
						printf("%d %d %d %f %s\n", q->id, q->pid, q->qid, q->z, v->color == White ? "White" : "Black");
					}
					for (set<MAPoint*>::iterator it = q->neighbors.begin(); it != q->neighbors.end(); ++it)
					{
						MAPoint* r = *it;
						if (r->z > u->key->z)
						{
							if (u->key->id == 74840)
							{
								printf("=== %d %d %d %f\n", r->id, r->pid, r->qid, r->z);
							}
							bOk = false;
							break;
						}
					}
					if (bOk)
					{
						Q2.push_back(v);
						v->color = Black;
						G.push_back(TimedPoint(v->key, u->key, iter));
					}
				}
			}
			if (u->key->id == 74840)
			{
				mexErrMsgTxt("Done for debugging.");
			}
		}
		Q = Q2;
	}

	reverseEdges(vertices);
	return G;
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
	vector<int> ids;
	{
		const int* dims;
		mxClassID classIdP;
		int ndimP;
		LoadData(ids, prhs[1], classIdP, ndimP, &dims);
	}

	float thres = std::numeric_limits<float>::infinity();
	if (nrhs >= 3)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(thres, prhs[2], classId);
	}
	float maxDist = std::numeric_limits<float>::infinity();
	if (nrhs >= 4)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(maxDist, prhs[3], classId);
	}
	MAPoint::_id = 0;

	Triangulation::Triangulator trmap(P);
	float unit = calculateScaleUnit(trmap); //set the unit

	vector<MAPoint*> mapoints = collectMedialAxisPoints(trmap, maxDist);
	collectNeighbors(mapoints, trmap, unit*1.25);

	GraphFactory<MAPoint*>& factory = GraphFactory<MAPoint*>::GetInstance();
	vector<Vertex<MAPoint*>*> vertices = buildDag(mapoints, thres);
	vector<Vertex<MAPoint*>*> peaks;
	for (int i = 0; i < ids.size(); ++i)
	{
		peaks.push_back(vertices[ids[i]]);
	}

	vector<TimedPoint> group = groupVertices(vertices, peaks);

	if (nlhs >= 1)
	{
		const int dims[] = { group.size(), 5 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], group[i].p->id);
			if (group[i].q == NULL)
			{
				SetData2(F, i, 1, dims[0], dims[1], -1);
			}
			else
			{
				SetData2(F, i, 1, dims[0], dims[1], group[i].q->id);
			}
			SetData2(F, i, 2, dims[0], dims[1], group[i].time);
			SetData2(F, i, 3, dims[0], dims[1], group[i].p->pid);
			SetData2(F, i, 4, dims[0], dims[1], group[i].p->qid);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
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

