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

struct MAPoint
{
	MAPoint(Triangulation::_Internal::_vertex* p, Triangulation::_Internal::_vertex* q, int pid, int qid)
	{
		this->x = (p->p.m_X + q->p.m_X) / 2.0;
		this->y = (p->p.m_Y	+ q->p.m_Y) / 2.0;
		this->z = 0; // Distance(p->p, q->p); // / 2.0; //TK!!!
		this->theta = GetVisualDirection(q->p.m_X, q->p.m_Y, p->p.m_X, p->p.m_Y);
		this->p = p;
		this->q = q;
		this->pid = pid;
		this->qid = qid;
		this->id = _id++;
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


float distance(MAPoint* p, MAPoint* q)
{
	CParticleF p0(p->x, p->y, p->z);
	CParticleF q0(q->x, q->y, q->z);
	CParticleF p2(p->x - sin(p->theta), p->y + cos(p->theta), p->z);
	CParticleF y = Closest2Line3d(p0, p2, q0);
	return Distance(y, q0);
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
collectMedialAxisPoints(Triangulation::Triangulator& trmap, float unit)
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

/*
Aggregation on the orientation.
*/
void
aggreagate1(vector<MAPoint*>& P, float sgm, float rate)
{
	float sgm2 = 2.*sgm*sgm;
	for (int j = 0; j < P.size(); ++j)
	{
		MAPoint* p = P[j];
		CParticleF p0(p->x, p->y, p->z);
		CParticleF r0(p->x + cos(p->theta), p->y + sin(p->theta), p->z);
		float sx = 0, sy = 0, sz = 0, sp = 0;
		for (int k = 0; k < P.size(); ++k)
		{
			MAPoint* q = P[k];
			CParticleF q0(q->x, q->y, q->z);
			CParticleF proj = Closest2Line(p0, r0, q0);
			float d2 = Distance(q0,proj);
			float wgt = exp(-d2 * d2/ sgm2);
			sx += cos(q->theta) * wgt;
			sy += sin(q->theta) * wgt;
			sp += wgt;
		}
		p->x0 = sx / sp;
		p->y0 = sy / sp;
		p->z0 = 0;
	}
	for (int j = 0; j < P.size(); ++j)
	{
		P[j]->theta = atan2(P[j]->y0, P[j]->x0);
	}
}

/*
Aggergation on position.
*/
void
aggreagate2(vector<MAPoint*>& P, float sgm, float rate)
{
	float sgm2 = 2.*sgm*sgm;
	for (int j = 0; j < P.size(); ++j)
	{
		MAPoint* p = P[j];
		float sx = 0, sy = 0, sz = 0, sp = 0;
		for (set<MAPoint*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
		{
			MAPoint* q = *it;
			float d2 = distance(p, q);
			//float d2 = Distance(p0, q0);
			float wgt =  exp(-d2 * d2 / sgm2); 
			sx += q->x * wgt;
			sy += q->y * wgt;
			sz += q->z * wgt;
			sp += wgt;
		}
		if (sp > 0)
		{
			p->x0 = sx / sp;
			p->y0 = sy / sp;
			p->z0 = sz / sp;
		}
		else
		{
			p->x0 = p->x;
			p->y0 = p->y;
			p->z0 = p->z;
		}
		if (p->id == 99 || p->id == 119)
		{
			printf("%d %f %f %f -> %f %f %f\n", p->id, sx, sy, sp, p->x0, p->y0, p->z0);
		}
	}
	for (int j = 0; j < P.size(); ++j)
	{
		MAPoint* p = P[j];
		if (p->id == 99 || p->id == 119)
		{
			printf("\t%d %f %f %f -> %f %f %f\n", p->id, p->x, p->y, p->z, (1 - rate)*p->x + rate * p->x0, (1 - rate)*p->y + rate * p->y0, p->z);
		}
		p->x = (1 - rate)*p->x + rate * p->x0;
		p->y = (1 - rate)*p->y + rate * p->y0;
		p->z = (1 - rate)*p->z + rate * p->z0;
	}
}

vector<int>
clusterMAPoints(vector<MAPoint*>& P, float thres)
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
		for (set<MAPoint*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
		{
			MAPoint* q = *it;
			Node<MAPoint*>* n = nodes[imap[q]];
			//float d = distance(p, q, scale, weight, false);
			if (p->id == 125 && q->id == 150)
			{
				i += 0;
			}
			float d1 = distance(p, q);
			float d2 = distance(q, p);
			float dm = Max(d1, d2);
			//float d = Max(1.0 - exp(-dm / thres), 1.0 - Abs(cos(p->theta - q->theta)));
			//float d = 1.0 - exp(-dm / thres);
			float d0 = Distance(p->x, p->y, q->x, q->y);
			if (d0 < 10) //0.2)
			{
				merge(nodes[i], n);
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
	float scale = 1.0f;
	if (nrhs >= 3)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(scale, prhs[2], classId);
	}
	int niter = 5;
	if (nrhs >= 4)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(niter, prhs[3], classId);
	}
	float sigma = 10.0;
	if (nrhs >= 5)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(sigma, prhs[4], classId);
	}
	float rate = 0.5;
	if (nrhs >= 6)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(rate, prhs[5], classId);
	}
	printf("sigma=%f, rate=%f, niter=%d\n", sigma, rate, niter);
	Triangulation::Triangulator trmap(P);
	float unit = calculateScaleUnit(trmap); //set the unit

	vector<MAPoint*> mapoints = collectMedialAxisPoints(trmap, unit);
	collectNeighbors(mapoints, trmap, unit*1.25);
	for (int i = 0; i < niter; ++i)
	{
		//aggreagate1(mapoints, sigma, rate);
		aggreagate2(mapoints, sigma*unit, rate);
	}
	vector<int> labels = clusterMAPoints(mapoints, thres*unit);
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
			vvint[k].push_back(i);
		}
		vector<vector<int>> vvint2;
		for (int i = 0; i < vvint.size(); ++i)
		{
			const int dims2[] = { vvint[i].size(), 3 };
			vector<int> F(dims2[0]*dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				int k = vvint[i][j];
				MAPoint* p = mapoints[k];
				SetData2(F, j, 0, dims2[0], dims2[1], p->pid);
				SetData2(F, j, 1, dims2[0], dims2[1], p->qid);
				SetData2(F, j, 2, dims2[0], dims2[1], p->id);
			}
			vvint2.push_back(F);
		}
		plhs[1] = StoreDataCell(vvint2, mxINT32_CLASS, 2, dims, 3);
	}
	if (nlhs >= 3)
	{
		const int dims[] = { mapoints.size(), 7 };
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
	/*for (int i = 0; i < trmap.points.size(); ++i)
	{
		links[trmap.points[i]]->clear();
		delete links[trmap.points[i]];
		links[trmap.points[i]] = NULL;
	}*/
	StationaryParticleFactory::getInstance().clean();
	mexUnlock();
}

