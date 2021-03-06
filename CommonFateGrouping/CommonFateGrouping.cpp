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

void
perturbePoints(vector<CParticleF>& points, float scale)
{
	rndm(12345);
	for (int i = 0; i < points.size(); ++i)
	{
		points[i].m_X += scale * (rndm(0) - 0.5);
		points[i].m_Y += scale * (rndm(0) - 0.5);
	}
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

vector<RawTriple*>
collectRawTriples(Triangulation::Triangulator& trmap, float thres, float separation)
{
	StationaryParticleFactory& factory = StationaryParticleFactory::getInstance();
	vector<StationaryParticle*> particles;
	map<Triangulation::_Internal::_vertex*, StationaryParticle*> pmap;
	map<StationaryParticle*,int> imap;
	vector<set<Triangulation::_Internal::_vertex*>> vpset(trmap.points.size());
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		StationaryParticle* sp = factory.makeParticle(trmap.points[i]->p);
		particles.push_back(sp);
		pmap[trmap.points[i]] = sp;
		imap[sp] = i;
	}

	for (int i = 0; i < trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* u = trmap.points[i];
		vector<Triangulation::_Internal::_vertex*> vs = collectNeighborNodes(u, 2);
		vector<pair<float, Triangulation::_Internal::_vertex*>> vpair;
		for (int j = 0; j < vs.size(); ++j)
		{
			float d = Distance(u->p, vs[j]->p);
			if (d < separation)
			{
				vpair.push_back(pair<float, Triangulation::_Internal::_vertex*>(d, vs[j]));
			}
		}
		sort(vpair.begin(), vpair.end());
		vector<float> angles;
		for (int j = 0; j < vpair.size(); ++j)
		{
			Triangulation::_Internal::_vertex* v = vpair[j].second;
			float a = GetVisualDirection(u->p.m_X, u->p.m_Y, v->p.m_X, v->p.m_Y);
			float sep = angleSepration(a, angles);
			//TK!!!
			//if (sep > PI / 4.0)
			if (sep > PI / 2.0)
			{
				int k = imap[pmap[v]];
				vpset[i].insert(trmap.points[k]);
				vpset[k].insert(trmap.points[i]);
				insertAngle(a, angles);
			}
		}
	}
	vector<RawTriple*> triples;
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		StationaryParticle* sp = pmap[trmap.points[i]];
		vector<StationaryParticle*> neighbors;
		for (set<Triangulation::_Internal::_vertex*>::iterator it = vpset[i].begin(); it != vpset[i].end(); ++it)
		{
			StationaryParticle* np = pmap[*it];
			neighbors.push_back(np);
		}
		for (int j = 0; j < neighbors.size(); ++j)
		{
			CParticleF q = neighbors[j]->getP();
			StationaryParticle* sq = neighbors[j];
			for (int k = 0; k < neighbors.size(); ++k)
			{
				if (k == j) continue;
				StationaryParticle* sr = neighbors[k];
				CParticleF r = neighbors[k]->getP();
				RawTriple* triple = new RawTriple(sp, sq, sr);
				//if (triple->fitness > thres)
				{
					triples.push_back(triple);
				}
			}
		}
	}
	for (int i = 0; i < triples.size(); ++i)
	{
		printf("%d: ", i);
		triples[i]->print();
	}
	return triples;
}

/*vector<RawTriple*>
collectRawTriples(Triangulation::Triangulator& trmap, float thres, float separation)
{
	StationaryParticleFactory& factory = StationaryParticleFactory::getInstance();
	vector<StationaryParticle*> particles;
	map<Triangulation::_Internal::_vertex*, StationaryParticle*> pmap;
	map<StationaryParticle*, int> imap;
	vector<set<Triangulation::_Internal::_vertex*>> vpset(trmap.points.size());
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		StationaryParticle* sp = factory.makeParticle(trmap.points[i]->p);
		particles.push_back(sp);
		pmap[trmap.points[i]] = sp;
		imap[sp] = i;
	}

	for (int i = 0; i < trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* u = trmap.points[i];
		StationaryParticle* sp = pmap[u];
		CParticleF p = sp->getP();
		for (int j = 0; j < trmap.points.size(); ++j)
		{
			Triangulation::_Internal::_vertex* v = trmap.points[j];
			StationaryParticle* sq = pmap[v];
			CParticleF q = sq->getP();
			if (Distance(p, q) < separation)
			{
				int k = imap[sq];
				vpset[i].insert(trmap.points[j]);
				vpset[j].insert(trmap.points[i]);
			}
		}
	}
	vector<RawTriple*> triples;
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		StationaryParticle* sp = pmap[trmap.points[i]];
		vector<StationaryParticle*> neighbors;
		for (set<Triangulation::_Internal::_vertex*>::iterator it = vpset[i].begin(); it != vpset[i].end(); ++it)
		{
			StationaryParticle* np = pmap[*it];
			neighbors.push_back(np);
		}
		for (int j = 0; j < neighbors.size(); ++j)
		{
			CParticleF q = neighbors[j]->getP();
			StationaryParticle* sq = neighbors[j];
			for (int k = 0; k < neighbors.size(); ++k)
			{
				if (k == j) continue;
				StationaryParticle* sr = neighbors[k];
				CParticleF r = neighbors[k]->getP();
				RawTriple* triple = new RawTriple(sp, sq, sr);
				if (triple->fitness > thres)
				{
					triples.push_back(triple);
				}
			}
		}
	}
	return triples;
}*/

vector<vector<RawTriple*>> 
clusterPairTriples(const vector<pair<RawTriple*,RawTriple*>>& pairs)
{
	vector<Node<int>*> nodes;
	for (int i = 0; i < pairs.size(); ++i)
	{
		nodes.push_back(makeset(i));
	}
	for (int i = 0; i < pairs.size(); ++i)
	{
		RawTriple* a = pairs[i].first;
		RawTriple* b = pairs[i].second;
		for (int j = i+1; j < pairs.size(); ++j)
		{
			RawTriple* x = pairs[j].first;
			RawTriple* y = pairs[j].second;
			if (RawTriple::related(a, b, x, y))
			{
				merge(nodes[i], nodes[j]);
				//printf("Merging %d %d %d %3.3f %3.f %d %d %d %3.3f %3.f\n",
				//	a->p->getId(), b->p->getId(), x->p->getId(), y->p->getId());
			}
		}
	}
	vector<Node<int>*> reps = clusters(nodes);
	vector<set<RawTriple*>> vset;
	map<int, int> imap;
	for (int i = 0; i < reps.size(); ++i)
	{
		int k = reps[i]->key;
		set<RawTriple*> gset;
		gset.insert(pairs[k].first);
		gset.insert(pairs[k].second);
		vset.push_back(gset);
		imap[k] = i;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = findset(nodes[i])->key;
		int j = imap[k];
		vset[j].insert(pairs[i].first);
		vset[j].insert(pairs[i].second);
	}

	vector<vector<RawTriple*>> groups;
	for (int i = 0; i < vset.size(); ++i)
	{
		vector<RawTriple*> vr;
		for (set<RawTriple*>::iterator it = vset[i].begin(); it != vset[i].end(); ++it)
		{
			vr.push_back(*it);
		}
		groups.push_back(vr);
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete(nodes[i]);
	}
	return groups;
}

vector<vector<LinkedTriple*>> 
clusterTriplesByP(vector<LinkedTriple*>& triples)
{
	vector<Node<LinkedTriple*>*> nodes;
	for (int i = 0; i < triples.size(); ++i)
	{
		nodes.push_back(makeset(triples[i]));
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		LinkedTriple* t = nodes[i]->key;
		for (int j = 0; j < nodes.size(); ++j)
		{
			if (i == j) continue;
			LinkedTriple* s = nodes[j]->key;
			if (t->p == s->p)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<vector<LinkedTriple*>> groups;
	vector<Node<LinkedTriple*>*> reps = clusters(nodes);
	vector<set<LinkedTriple*>> vset;
	map<LinkedTriple*, int> imap;
	for (int i = 0; i < reps.size(); ++i)
	{
		set<LinkedTriple*> gset;
		gset.insert(reps[i]->key);
		vset.push_back(gset);
		imap[reps[i]->key] = i;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		LinkedTriple* t = findset(nodes[i])->key;
		int j = imap[t];
		vset[j].insert(nodes[i]->key);
	}

	vector<vector<LinkedTriple*>> groups2;
	for (int i = 0; i < vset.size(); ++i)
	{
		vector<LinkedTriple*> vr;
		for (set<LinkedTriple*>::iterator it = vset[i].begin(); it != vset[i].end(); ++it)
		{
			vr.push_back(*it);
		}
		groups2.push_back(vr);
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete(nodes[i]);
	}
	return groups2;
}

void
examineClusterTriples(vector<vector<RawTriple*>> &groups)
{
	for (int i = 0; i < groups.size(); ++i)
	{
		set<int> iset;
		for (int j = 0; j < groups[i].size(); ++j)
		{
			iset.insert(groups[i][j]->p->getId());
		}
		vector<int> vrt;
		for (set<int> ::iterator it = iset.begin(); it != iset.end(); ++it)
		{
			vrt.push_back((*it));
		}
		sort(vrt.begin(), vrt.end());
		printf("%d (%d): ", i, iset.size());
		for (vector<int> ::iterator it = vrt.begin(); it != vrt.end(); ++it)
		{
			//printf("[%d %d %d] ", (*it)->r->getId(), (*it)->p->getId(), (*it)->q->getId());
			//printf("%d ", (*it)->r->getId());
			printf("%d ", (*it));
		}
		printf("\n");
	}
}


vector<vector<RawTriple*>>
groupRawTriples(vector<RawTriple*>& rt, float thres)
{
	vector<pair<RawTriple*, RawTriple*>> pairs;
	for (int i = 0; i < rt.size(); ++i)
	{
		RawTriple* t1 = rt[i];
		StationaryParticle* p = t1->p;
		for (int j = i + 1; j < rt.size(); ++j)
		{
			RawTriple* t2 = rt[j];
			StationaryParticle* q = t2->p;
			float comp = LinkedTriple::compatibility(p->getP(), t1->v, q->getP(), t2->v);
			if (comp > thres)
			{
				pairs.push_back(pair<RawTriple*, RawTriple*>(t1, t2));
			}
		}
	}
	vector<vector<RawTriple*>> groups = clusterPairTriples(pairs);
	examineClusterTriples(groups);

	return groups;
}

vector<LinkedTriple*>
collectLinkedTriples(vector<vector<RawTriple*>>& groups, float thres)
{
	LinkedTripleFactory& tfactory = LinkedTripleFactory::getInstance();
	vector<LinkedTriple*> triples;
	map<LinkedTriple*, int> gmap;
	for (int i = 0; i < groups.size(); ++i)
	{
		set<RawTriple*> rset;
		vector<LinkedTriple*> vt;
		for (int j = 0; j < groups[i].size(); ++j)
		{
			rset.insert(groups[i][j]);
		}
		for (set<RawTriple*>::iterator it = rset.begin(); it != rset.end(); ++it)
		{
			RawTriple* raw = *it;
			LinkedTriple* tr = tfactory.makeTriple(raw->p, raw->r, raw->q);
			vt.push_back(tr);
			gmap[tr] = i;
			tr->groupNumber = i;
		}
		for (int j = 0; j < vt.size(); ++j)
		{
			LinkedTriple* t1 = vt[j];
			for (int k = j + 1; k < vt.size(); ++k)
			{
				LinkedTriple* t2 = vt[k];
				float comp = LinkedTriple::compatibility(t1->p->getP(), t1->v, t2->p->getP(), t2->v);
				if (comp > thres)
				{
					t1->frontSupporters.push_back(t2);
					t2->frontSupporters.push_back(t1);
				}
			}
		}
		triples.insert(triples.end(), vt.begin(), vt.end());
	}

	//add left, right supporters
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* t1 = triples[i];
		int g1 = gmap[t1];
		for (int j = i + 1; j < triples.size(); ++j)
		{
			LinkedTriple* t2 = triples[j];
			//if (gmap[t2] == g1)
			{
				if (t1->p == t2->r && t1->q == t2->p)
				{
					t1->rightSupporters.push_back(t2);
					t2->leftSupporters.push_back(t1);
				}
				if (t1->p == t2->q && t1->r == t2->p)
				{
					t1->leftSupporters.push_back(t2);
					t2->rightSupporters.push_back(t1);
				}
			}
		}
	}

	//add competitors
	{
		vector<vector<LinkedTriple*>> groups = clusterTriplesByP(triples);
		for (int i = 0; i < groups.size(); ++i)
		{
			LinkedTriple* nil = tfactory.makeNillTriple(groups[i][0]->p);
			triples.push_back(nil);
			groups[i].push_back(nil);
			for (int j = 0; j < groups[i].size(); ++j)
			{
				LinkedTriple* t = groups[i][j];
				for (int k = 0; k < groups[i].size(); ++k)
				{
					if (j == k) continue;
					LinkedTriple* s = groups[i][k];
					t->competitors.push_back(s);
				}
			}
		}
		for (int i = 0; i < triples.size(); ++i)
		{
			triples[i]->prob = 1.0 / (triples[i]->competitors.size() + 1.0);
		}
	}
	return triples;
}

void support(vector<LinkedTriple*>& triples, float thres)
{
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* tr = triples[i];
		float sum1 = 0; // tr->fitness;
		for (int j = 0; j < tr->frontSupporters.size(); ++j)
		{
			LinkedTriple* tr2 = tr->frontSupporters[j];
			sum1 += tr2->prob * tr->compatibility(tr2);
		}
		float sum2 = 0; // tr->fitness;
		for (int j = 0; j < tr->leftSupporters.size(); ++j)
		{
			LinkedTriple* tr2 = tr->leftSupporters[j];
			sum2 += tr2->prob; // *tr->compatibility0(tr2);
		}
		float sum3 = 0; // tr->fitness;
		for (int j = 0; j < tr->rightSupporters.size(); ++j)
		{
			LinkedTriple* tr2 = tr->rightSupporters[j];
			sum3 += tr2->prob; // *tr->compatibility0(tr2);
		}
		//tr->_totalFitness = tr->prob * tr->fitness + pow(sum1*sum2*sum3, 1.0 / 3.0);
		tr->_totalFitness = tr->fitness + pow(sum1*sum2*sum3, 1.0 / 3.0);
	}
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* tr = triples[i];
		float sum = tr->prob * tr->_totalFitness;
		for (int j = 0; j < tr->competitors.size(); ++j)
		{
			LinkedTriple* tr2 = tr->competitors[j];
			sum += tr2->prob * tr2->_totalFitness;
		}
		if (sum >0)
		{
			tr->_prob0 = tr->prob * tr->_totalFitness / sum;
		}
		else
		{
			tr->_prob0 = 0.0f;
		}
	}
	for (int i = 0; i < triples.size(); ++i)
	{
		triples[i]->prob = triples[i]->_prob0;
	}
	for (int i = 0; i < triples.size(); ++i)
	{
		//if (triples[i]->p->getId() == 25)
		if (triples[i]->p->getId() == 25 && triples[i]->r->getId() == 26 && triples[i]->q->getId() == 24)
		{
			triples[i]->printDetail();
		}
	}
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
	Triangulation::Triangulator trmap(points);
	(LinkedTripleFactory::getInstance()).unit = calculateScaleUnit(trmap); //set the unit

	float perturbationScale = 0.0f;
	float thres1 = 0.5;
	float thres2 = 0.8;
	int numIter = 10;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres1, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(thres2, prhs[2], classMode);
	}
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[3], classMode);
	}
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(perturbationScale, prhs[4], classMode);
	}

	perturbePoints(points, perturbationScale);
	//vector<StationaryParticle*> particles = collectStationaryParticles(points);
	vector<RawTriple*> rawt = collectRawTriples(trmap, thres1, 1.25*(LinkedTripleFactory::getInstance()).unit);
	vector<vector<RawTriple*>> groups = groupRawTriples(rawt, thres2);
	vector<LinkedTriple*> triples = collectLinkedTriples(groups, thres2);

	for (int i = 0; i < numIter; ++i)
	{
		printf("Iteration %d:\n", i);
		support(triples, 0.01);
		//updateFate(triples);
	}
	for (int i = 0; i < triples.size(); ++i)
	{
		triples[i]->updateFate();
		triples[i]->print();
	}

	if (nlhs >= 1)
	{
		StationaryParticleFactory& spfactory = StationaryParticleFactory::getInstance();
		const int dims[] = { spfactory.particles.size(), 3 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], spfactory.particles[i]->getX());
			SetData2(F, i, 1, dims[0], dims[1], spfactory.particles[i]->getY());
			SetData2(F, i, 2, dims[0], dims[1], (float)spfactory.particles[i]->getId());
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	if (nlhs >= 2)
	{
		LinkedTripleFactory& factory = LinkedTripleFactory::getInstance();
		const int dims[] = { factory.triples.size(), 10 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			LinkedTriple* t = factory.triples[i];
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
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		const int dims[] = { triples.size(), 2 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			LinkedTriple* t = triples[i]->best();
			SetData2(F, i, 0, dims[0], dims[1], triples[i]->id);
			if (t == NULL)
			{
				SetData2(F, i, 1, dims[0], dims[1], 0);
			}
			else {
				SetData2(F, i, 1, dims[0], dims[1], t->id);
			}
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);

	}
	if (nlhs >= 4)
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
		plhs[3] = StoreDataCell(vvint, mxINT32_CLASS, 2, dims, 3);
	}
	for (int i = 0; i < rawt.size(); ++i)
	{
		delete rawt[i];
	}
	StationaryParticleFactory::getInstance().clean();
	LinkedTripleFactory::getInstance().clean();
	mexUnlock();
}

