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

struct MAPoint
{
	MAPoint(Triangulation::_Internal::_vertex* p, Triangulation::_Internal::_vertex* q, int pid, int qid)
	{
		this->x = (p->p.m_X + q->p.m_X) / 2.0;
		this->y = (p->p.m_Y	+ q->p.m_Y) / 2.0;
		this->z = Distance(p->p, q->p) / 10.0; //TK!!!
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
	CParticleF q0(q->x, q->y, q->z);
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
			if (Distance(trmap.points[i]->p, trmap.points[j]->p) <= thres) 
			{
				mapoints.push_back(new MAPoint(trmap.points[i], trmap.points[j], i, j));
			}
		}
	}
	return mapoints;
}

bool
forwardLooking(MAPoint* p, MAPoint* ref)
{
	float x1 = cos(p->theta);
	float y1 = sin(p->theta);
	float x2 = cos(ref->theta);
	float y2 = sin(ref->theta);
	float a = GetVisualAngle(x1, y1, x2, y2, 0.0, 0.0);
	if (a > PI / 2.0)
	{
		x1 = -x1;
		y1 = -y1;
	}
	float a2 = GetVisualAngle2(x1, y1, x2, y2, 0.0, 0.0);
	if (a2 > 0.0) return true;
	else if (a2 < 0.0) return false;
	else {
		return p->id > ref->id;
	}
}

/*
*/
vector<Vertex<MAPoint*>*>
buildDAG(vector<MAPoint*>& P, MAPoint* ref)
{
	GraphFactory<MAPoint*>& factory = GraphFactory<MAPoint*>::GetInstance();
	vector<Vertex<MAPoint*>*> vertices(P.size()+1);
	Vertex<MAPoint*>* src = NULL;
	map<MAPoint*, int>  imap;
	for (int i = 0; i < P.size(); ++i)
	{
		vertices[i] = factory.makeVertex(P[i]);
		if (P[i] == ref) src = vertices[i];
		imap[P[i]] = i;
	}
	vertices[P.size()] = factory.makeVertex(ref);
	Vertex<MAPoint*>* sink = vertices[P.size()];

	for (int i = 0; i < P.size(); ++i)
	{
		MAPoint* p = P[i];
		for (set<MAPoint*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it) 
		{
			MAPoint* q = *it;
			if (p == q) continue;
			int k = imap[q];
			if (forwardLooking(q, p))
			{
				Vertex<MAPoint*>* v = (q == ref) ? sink : vertices[k]; //separate source and sink
				/*float dpp = Distance(p->p->p, q->p->p);
				float dpq = Distance(p->p->p, q->q->p);
				float dqp = Distance(p->q->p, q->p->p);
				float dqq = Distance(p->q->p, q->q->p);
				float w = 0;
				if (dpp + dqq < dpq + dqp)
				{
					w = Max(dpp, dqq);
				}
				else
				{
					w = Max(dpq, dqp);
				}*/
				//float w = sqrt((p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y) + (p->z - q->z)*(p->z - q->z));
				float w = distance(p, q);
				Edge<MAPoint*>* ed = factory.makeEdge(vertices[i], v, w*w);
				vertices[i]->Add(ed);
			}
		}
	}
	return vertices;
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

	Triangulation::Triangulator trmap(P);
	float unit = calculateScaleUnit(trmap); //set the unit

	vector<MAPoint*> mapoints = collectMedialAxisPoints(trmap, thres);
	collectNeighbors(mapoints, trmap, unit*1.25);
	MAPoint* ref = NULL;
	for (int i = 0; i < mapoints.size(); ++i)
	{
		if (mapoints[i]->pid == Min(ids[0], ids[1]) && mapoints[i]->qid == Max(ids[0], ids[1]))
		{
			ref = mapoints[i];
			break;
		}
	}

	GraphFactory<MAPoint*>& factory = GraphFactory<MAPoint*>::GetInstance();
	vector<Vertex<MAPoint*>*> vertices = buildDAG(mapoints, ref);
	Vertex<MAPoint*>* src = NULL;
	Vertex<MAPoint*>* sink = vertices[vertices.size()-1];
	for (int i = 0; i < mapoints.size(); ++i)
	{
		if (vertices[i]->key == ref)
		{
			src = vertices[i];
			break;
		}
	}

	Dijkstra(vertices, src);
	vector<Vertex<MAPoint*>*> path = tracePath(sink, src);

	if (nlhs >= 1)
	{
		const int dims[] = { path.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], path[i]->key->id);
			SetData2(F, i, 1, dims[0], dims[1], path[i]->key->pid);
			SetData2(F, i, 2, dims[0], dims[1], path[i]->key->qid);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { factory.edges.size(), 4 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			Edge<MAPoint*>* ed = factory.edges[i];
			SetData2(F, i, 0, dims[0], dims[1], ed->u->key->pid);
			SetData2(F, i, 1, dims[0], dims[1], ed->u->key->qid);
			SetData2(F, i, 2, dims[0], dims[1], ed->v->key->pid);
			SetData2(F, i, 3, dims[0], dims[1], ed->v->key->qid);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		const int dims[] = { mapoints.size(), 10};
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
			SetData2(F, i, 7, dims[0], dims[1], mapoints[i]->x0);
			SetData2(F, i, 8, dims[0], dims[1], mapoints[i]->y0);
			SetData2(F, i, 9, dims[0], dims[1], mapoints[i]->z0);
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	for (int i = 0; i < mapoints.size(); ++i)
	{
		delete mapoints[i];
	}
	GraphFactory<MAPoint*>::GetInstance().Clean();
	/*for (int i = 0; i < trmap.points.size(); ++i)
	{
		links[trmap.points[i]]->clear();
		delete links[trmap.points[i]];
		links[trmap.points[i]] = NULL;
	}*/
	//StationaryParticleFactory::getInstance().clean();
	mexUnlock();
}

