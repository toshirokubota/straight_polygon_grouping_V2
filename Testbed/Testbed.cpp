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
#include <LinkedTriple.h>
#include <StationaryParticle.h>
#include <IntersectionConvexPolygons.h>
#include <Triangulation.h>

struct RawTriple {
	RawTriple(StationaryParticle* p0, StationaryParticle* q0, StationaryParticle* r0)
	{
		p = p0;
		q = q0;
		r = r0;
		LinkedTriple::calculateVelocity(p, q, r, v[0], v[1]);
		fitness = LinkedTriple::fitnessMeasure(p, q, r);
	}

	void print()
	{
		printf("%d %d %d %3.3f %3.3f %3.3f\n", p->getId(), r->getId(), q->getId(), v[0], v[1], fitness);
	}
	bool _select(float length)
	{
		return Distance(p->getP(), r->getP()) < length && Distance(p->getP(), q->getP()) < length;
	}

	static bool identical(RawTriple* a, RawTriple* x)
	{
		return a->p == x->p && a->q == x->q && a->r == x->r;
	}

	static bool overlapping(RawTriple* a, RawTriple* x)
	{
		return a->p == x->r && a->q == x->p || a->p == x->q && a->r == x->p
			|| a->p == x->p && a->q == x->q || a->p == x->p && a->r == x->r;
	}

	static bool related(RawTriple* a, RawTriple* b, RawTriple* x, RawTriple* y)
	{
		return identical(a, x) && overlapping(b, y) || identical(a, y) && overlapping(b, x) ||
			identical(b, x) && overlapping(a, y) || identical(b, y) && overlapping(a, x)
			|| overlapping(a, x) && overlapping(b, y) || overlapping(a, y) && overlapping(b, x);
	}

	StationaryParticle* p; //center
	StationaryParticle* q; //right
	StationaryParticle* r; //left
	float v[2];
	float fitness;
};

bool identical(LinkedTriple* a, LinkedTriple* x)
{
	return a->p == x->p && a->q == x->q && a->r == x->r;
}

bool overlapping(LinkedTriple* a, LinkedTriple* x)
{
	return a->p == x->r && a->q == x->p || a->p == x->q && a->r == x->p
		|| a->p == x->p && a->q == x->q || a->p == x->p && a->r == x->r;
}

bool related(LinkedTriple* a, LinkedTriple* b, LinkedTriple* x, LinkedTriple* y)
{
	return identical(a, x) && overlapping(b, y) || identical(a, y) && overlapping(b, x) ||
		identical(b, x) && overlapping(a, y) || identical(b, y) && overlapping(a, x)
		|| overlapping(a, x) && overlapping(b, y) || overlapping(a, y) && overlapping(b, x);
}


/*
both angles are in [-PI, PI).
The function returns the angle (in radian) between the two in [0, PI).
*/
float angleDiff(float a, float b)
{
	if (a > b) swap(a, b);

	return Min(Abs(b - a), Abs(a + 2 * PI - b));
}

/*
ANGLES is a vector of angles (in radian) sorted in ascending order.
The function finds the first INDEX where ANGLES[INDEX]< a. It then returns the minimum of A-ANGLES[INDEX] and ANGLES[INDEX+1]-A
where the indexing is done in modulo.
*/
float 
angleSepration(float a, vector<float>& angles)
{
	int n = angles.size();
	if (n == 0)  return 2 * PI;
	else if (n == 1) return angleDiff(a, angles[0]);
	int idx = 0;
	while (idx < angles.size() && angles[idx] < a)
	{
		idx++;
	}
	return Min(angleDiff(a, angles[idx % n]), angleDiff(a, angles[(idx + n - 1) % n]));
}

/*
insert a value into a sorted linear array.
*/
void
insertAngle(float a, vector<float>&  angles)
{
	if (angles.empty()) angles.push_back(a);
	else
	{
		int idx = 0;
		while (idx < angles.size() && angles[idx] < a)
		{
			idx++;
		}
		angles.insert(angles.begin() + idx, a);
	}
}

vector<Triangulation::_Internal::_vertex*>
collectNeighborNodes(Triangulation::_Internal::_vertex* u, int maxHop)
{
	set<Triangulation::_Internal::_vertex*> vset;
	vset.insert(u);
	vector<Triangulation::_Internal::_vertex*> Q(1, u);
	int hop = 1;
	while (Q.empty() == false && hop <= maxHop)
	{
		vector<Triangulation::_Internal::_vertex*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			Triangulation::_Internal::_vertex* v = Q[i];
			for (int j = 0; j < v->edges.size(); ++j)
			{
				Triangulation::_Internal::_edge* ed = v->edges[j];
				Triangulation::_Internal::_vertex* w = ed->vertices[0] == v ? ed->vertices[1] : ed->vertices[0];
				if (vset.find(w) == vset.end())
				{
					vset.insert(w);
					Q2.push_back(w);
				}
			}
		}
		Q = Q2;
		hop++;
	}
	vset.erase(u); //exclude the source.
	vector<Triangulation::_Internal::_vertex*> vv;
	for (set<Triangulation::_Internal::_vertex*>::iterator it = vset.begin(); it != vset.end(); ++it)
	{
		vv.push_back(*it);
	}
	return vv;
}

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

vector<StationaryParticle*>
collectStationaryParticles(vector<CParticleF>& points)
{
	StationaryParticleFactory& factory = StationaryParticleFactory::getInstance();
	vector<StationaryParticle*> particles;
	for (int i = 0; i < points.size(); ++i)
	{
		StationaryParticle* sp = factory.makeParticle(points[i]);
		particles.push_back(sp);
	}
	return particles;
}


void
examineClusterTriples(vector<vector<LinkedTriple*>>& groups, LinkedTriple* ref)
{
	for (int i = 0; i < groups.size(); ++i)
	{
		set<int> iset;
		bool bMatch = false; //true if the group contains the reference triple
		for (int j = 0; j < groups[i].size(); ++j)
		{
			if (groups[i][j]->id == ref->id) bMatch = true;
			iset.insert(groups[i][j]->p->getId());
		}
		if (bMatch == false) continue;

		vector<int> vrt;
		for (set<int> ::iterator it = iset.begin(); it != iset.end(); ++it)
		{
			vrt.push_back((*it));
		}
		sort(vrt.begin(), vrt.end());
		printf("%d (%d): ", i, iset.size());
		for (vector<int> ::iterator it = vrt.begin(); it != vrt.end(); ++it)
		{
			printf("%d ", (*it));
		}
		printf("\n");
	}
}

vector<LinkedTriple*> 
compatiblePairs(vector<LinkedTriple*>& triples, LinkedTriple* ref, float thres)
{
	vector<LinkedTriple*> pairs;
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* t = triples[i];
		if (t == ref) continue;
		float comp = LinkedTriple::compatibility(ref->p->getP(), ref->v, t->p->getP(), t->v);
		if (comp > thres)
		{
			pairs.push_back(t);
		}
	}
	return pairs;
}

vector<vector<LinkedTriple*>> groupTriples(vector<LinkedTriple*>& pairs, LinkedTriple* ref, float thres)
{
	vector<Node<LinkedTriple*>*> nodes;
	for (int i = 0; i < pairs.size(); ++i)
	{
		nodes.push_back(makeset(pairs[i]));
	}
	for (int i = 0; i < pairs.size(); ++i)
	{
		for (int j = i + 1; j < pairs.size(); ++j)
		{
			if (related(ref, pairs[i], ref, pairs[j]))
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}

	vector<Node<LinkedTriple*>*> reps = clusters(nodes);
	map<LinkedTriple*, int> imap;
	vector<vector<LinkedTriple*>> groups(reps.size());
	for (int i = 0; i < reps.size(); ++i)
	{
		LinkedTriple* k = reps[i]->key;
		imap[k] = i;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		LinkedTriple* k = findset(nodes[i])->key;
		int j = imap[k];
		groups[j].push_back(nodes[i]->key);
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete(nodes[i]);
	}
	return groups;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "StraightMedialAxis", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = StraightMedialAxis(P, [iter delta])");
		return;
	}
	//Points
	vector<CParticleF> points;
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
			points.push_back(CParticleF(x, y));
		}

	}
	vector<StationaryParticle*> particles = collectStationaryParticles(points);
	Triangulation::Triangulator trmap(points);
	(LinkedTripleFactory::getInstance()).unit = calculateScaleUnit(trmap); //set the unit

	//Triples
	vector<LinkedTriple*> triples;
	LinkedTripleFactory&  lfactory = LinkedTripleFactory::getInstance();
	{
		mxClassID classIdP;
		const int* dimsP;
		int ndimP;
		vector<int> P0;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		for (int i = 0; i < dimsP[0]; ++i)
		{
			int p = GetData2(P0, i, 0, dimsP[0], dimsP[1], 0);
			int r = GetData2(P0, i, 1, dimsP[0], dimsP[1], 0);
			int q = GetData2(P0, i, 2, dimsP[0], dimsP[1], 0);
			triples.push_back(lfactory.makeTriple(particles[p], particles[r], particles[q]));
		}
	}
	//reference Triple
	LinkedTriple* ref = NULL;
	{
		mxClassID classIdP;
		int ndimP;
		vector<int> idx;
		LoadData(idx, prhs[2], classIdP, ndimP, &dimsP);
		for (int i = 0; i < triples.size(); ++i)
		{
			if (triples[i]->p->getId() == idx[0] && triples[i]->r->getId() == idx[1] && triples[i]->q->getId() == idx[2])
			{
				ref = triples[i];
				break;
			}
		}
	}
	float thres = 0.75;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[3], classMode);
	}

	vector<LinkedTriple*> pairs = compatiblePairs(triples, ref, thres);
	vector<vector<LinkedTriple*>> groups = groupTriples(pairs, ref, thres);

	if (nlhs >= 1)
	{
		const int dims[] = { pairs.size(), 10 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			LinkedTriple* t = pairs[i];
			t->updateFate();
			SetData2(F, i, 0, dims[0], dims[1], (float)t->id);
			SetData2(F, i, 1, dims[0], dims[1], (float)t->p->getId());
			SetData2(F, i, 2, dims[0], dims[1], (float)t->r->getId());
			SetData2(F, i, 3, dims[0], dims[1], (float)t->q->getId());
			SetData2(F, i, 4, dims[0], dims[1], t->fate.m_X);
			SetData2(F, i, 5, dims[0], dims[1], t->fate.m_Y);
			SetData2(F, i, 6, dims[0], dims[1], t->v[0]);
			SetData2(F, i, 7, dims[0], dims[1], t->v[1]);
			SetData2(F, i, 8, dims[0], dims[1], t->fitness);
			SetData2(F, i, 9, dims[0], dims[1], t->prob);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { groups.size(), 1 };
		vector<vector<int>> vvint;
		for (int i = 0; i < groups.size(); ++i)
		{
			const int dims2[] = { groups[i].size(), 3 };
			vector<int> vint(dims2[0] * dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				SetData2(vint, j, 0, dims2[0], dims2[1], groups[i][j]->p->getId());
				SetData2(vint, j, 1, dims2[0], dims2[1], groups[i][j]->r->getId());
				SetData2(vint, j, 2, dims2[0], dims2[1], groups[i][j]->q->getId());
			}
			vvint.push_back(vint);
		}
		plhs[1] = StoreDataCell(vvint, mxINT32_CLASS, 2, dims, 3);
	}

	StationaryParticleFactory::getInstance().clean();
	LinkedTripleFactory::getInstance().clean();
	mexUnlock();
}

