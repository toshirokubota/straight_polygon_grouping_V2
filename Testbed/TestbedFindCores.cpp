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
			int k1 = N[p->pid][i];
			for (int j = 0; j < N[p->qid].size(); ++j)
			{
				int k2 = N[p->qid][j];
				int key = _key(k1, k2, n);
				if (key <= maxKey)
				{
					MAPoint* q = pmap[key];
					if (q == p) continue;
					if (q != NULL)
					{
						if (q->p == p->p || q->p == p->q || q->q == p->p || q->q == p->q)
						{
							i += 0;
						}
						//else //no shared vertex
						{
							p->neighbors.insert(q);
						}
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

vector<MAPoint*> locateCores(vector<MAPoint*>& P, float thres, float thres2)
{
	vector<MAPoint*> cores;
	for (int i = 0; i < P.size(); ++i)
	{
		MAPoint* p = P[i];
		bool bCore = true;
		vector<MAPoint*> up;
		vector<MAPoint*> down;
		if (p->pid == 0 && p->qid == 25)
		{
			i += 0;
		}
		for (set<MAPoint*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
		{
			MAPoint* q = *it;
			if (Abs(cos(p->theta - q->theta)) > thres) //check for parallelness
			{
				float a = GetVisualAngle2(p->x + cos(p->theta), p->y + sin(p->theta), q->x, q->y, p->x, p->y);
				if (sin(a) > thres2)
				{
					up.push_back(q);
				}
				else if (-sin(a) > thres2)
				{
					down.push_back(q);
				}
			}
		}
		if (up.size() > 0 && down.size() > 0)
		{
			bool bOk = true;
			for (int i = 0; i < up.size(); ++i)
			{
				if (up[i]->z > p->z) bOk = false;
			}
			for (int i = 0; i < down.size(); ++i)
			{
				if (down[i]->z > p->z) bOk = false;
			}
			if (bOk)
			{
				cores.push_back(p);
			}
		}
	}

	return cores;
}

vector<MAPoint*>
mergeCores(vector<MAPoint*>& cores)
{
	vector<Node<MAPoint*>*> nodes;
	for (int i = 0; i < cores.size(); ++i)
	{
		nodes.push_back(makeset(cores[i]));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		MAPoint* p = nodes[i]->key;
		for (int j = 0; j < nodes.size(); ++j)
		{
			if (i == j) continue;
			MAPoint* q = nodes[j]->key;
			if (p->id == 0 && q->id == 25)
			{
				i += 0;
			}
			if (p->neighbors.find(q) != p->neighbors.end())
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	//cluster cores that are neighbors
	vector<Node<MAPoint*>*> reps = clusters(nodes);

	//use the one with the largest Z value as the representative.
	vector<MAPoint*> X;
	map<MAPoint*, int> pmap;
	for (int i = 0; i < reps.size(); ++i)
	{
		X.push_back(reps[i]->key);
		pmap[reps[i]->key] = i;
	}

	for (int j = 0; j < cores.size(); ++j)
	{
		MAPoint* r = findset(nodes[j])->key;
		int k = pmap[r];
		if (r->z > X[k]->z)
		{
			X[k] = r;
		}
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return X;
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

	float thres = 0.99;
	if (nrhs >= 2)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(thres, prhs[1], classId);
	}
	float thres2 = 0.95;
	if (nrhs >= 3)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(thres2, prhs[2], classId);
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

	vector<MAPoint*> cores = locateCores(mapoints, thres, thres2);
	cores = mergeCores(cores);
	for (int i = 0; i < cores.size(); ++i)
	{
	printf("Core: %d %d %d\n", cores[i]->id, cores[i]->pid, cores[i]->qid);
	}

	if (nlhs >= 1)
	{
		const int dims[] = { cores.size(), 3 };
		vector<int> F(dims[0]*dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], cores[i]->id);
			SetData2(F, i, 1, dims[0], dims[1], cores[i]->pid);
			SetData2(F, i, 2, dims[0], dims[1], cores[i]->qid);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		plhs[1] = StoreScalar(0, mxINT32_CLASS);
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

