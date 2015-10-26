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

float calculateScaleUnit(vector<CParticleF>& points)
{
	Triangulation::Triangulator trmap(points);

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

vector<LinkedTriple*>
collectTriples(vector<CParticleF>& points, float thres, float scale)
{
	StationaryParticleFactory& factory = StationaryParticleFactory::getInstance();
	LinkedTripleFactory& tfactory = LinkedTripleFactory::getInstance();
	vector<StationaryParticle*> particles;
	float separation = scale * tfactory.unit;

	for (int i = 0; i < points.size(); ++i)
	{
		StationaryParticle* sp = factory.makeParticle(points[i]);
		particles.push_back(sp);
	}

	vector<LinkedTriple*> triples;
	vector<LinkedTriple*> nils;
	for (int i = 0; i < points.size(); ++i)
	{
		CParticleF p = points[i];
		StationaryParticle* sp = particles[i];
		vector<int> neighbors;
		vector<LinkedTriple*> triples0;
		for (int j = 0; j < points.size(); ++j)
		{
			if (i == j) continue;
			if (Distance(p, points[j]) < separation)
			{
				neighbors.push_back(j);
			}
		}
		for (int j = 0; j < neighbors.size(); ++j)
		{
			CParticleF q = points[neighbors[j]];
			StationaryParticle* sq = particles[neighbors[j]];
			for (int k = 0; k < neighbors.size(); ++k)
			{
				if (k == j) continue;
				CParticleF r = points[neighbors[k]];
				StationaryParticle* sr = particles[neighbors[k]];
				if (LinkedTriple::fitnessMeasure(sp, sq, sr) > thres)
				{
					triples0.push_back(tfactory.makeTriple(sp, sq, sr));
				}
			}
		}
		LinkedTriple* nil = tfactory.makeNillTriple(sp);
		nils.push_back(nil);

		//add them as competitors
		for (int j = 0; j < triples0.size(); ++j)
		{
			for (int k = 0; k < triples0.size(); ++k)
			{
				if (j == k) continue;
				triples0[j]->competitors.push_back(triples0[k]);
			}
			triples0[j]->competitors.push_back(nil);
		}
		for (int j = 0; j < triples0.size(); ++j)
		{
			nil->competitors.push_back(triples0[j]);
		}
		triples.insert(triples.end(), triples0.begin(), triples0.end());
	}
	triples.insert(triples.end(), nils.begin(), nils.end());
	for (int i = 0; i < triples.size(); ++i)
	{
		triples[i]->prob = 1.0 / ((float) triples[i]->competitors.size()+1.0);
	}
	return triples;
}

bool
circularSupport(LinkedTriple* t, LinkedTriple* s, float time, float thres)
{
	CParticleF a = t->p->getP();
	CParticleF b = s->p->getP();
	CParticleF o(a.m_X + time*t->v[0], a.m_Y + time*t->v[1]);
	float ang = GetVisualAngle(a.m_X, a.m_Y, b.m_X, b.m_Y, o.m_X, o.m_Y);
	float sn = sin(ang / 2.0);
	return sn * sn > thres;
}

void
assignSupporters(vector<LinkedTriple*>& triples, float thres)
{
	float toobig = 1.0e4;
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* t1 = triples[i];
		if (LinkedTriple::isNil(t1)) continue; //this is a NIL

		CParticleF p1 = t1->p->getP();
		CParticleF p2(p1.m_X + t1->v[0], p1.m_Y + t1->v[1]);
		for (int j = 0; j < triples.size(); ++j)
		{
			LinkedTriple* t2 = triples[j];
			if (t1->p == t2->p) continue;
			if (LinkedTriple::isNil(t2)) continue;
			if (circularSupport(t1, t2, t1->_timeToClosestEncounter(t2), 0.5f))
			{
				CParticleF q1 = t2->p->getP();
				CParticleF q2(q1.m_X + t2->v[0], q1.m_Y + t2->v[1]);
				float comp = t1->compatibility(t2);
				/*if (t1->p->getId() == 25)
				{
				printf("[%d %d] %d (%d,%d,%d) vs. %d (%d, %d, %d) => %f\n",
				i, j, t1->id, t1->p->getId(), t1->r->getId(), t1->q->getId(), t2->id, t2->p->getId(), t2->r->getId(), t2->q->getId(), comp);
				}*/
				if (comp > thres)
				{
					t1->frontSupporters.push_back(t2);
					t1->linkWeights.push_back(1.0);
				}
			}
			if (t2->q == t1->p && t2->p == t1->r)
			{
				t1->leftSupporters.push_back(t2);
			}
			if (t2->r == t1->p && t2->p == t1->q)
			{
				t1->rightSupporters.push_back(t2);
			}
		}
	}
}

/*
void
updateFate(vector<LinkedTriple*>& triples)
{
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* t = triples[i];
		float ex = 0, ey = 0, sump = 0;
		for (int j = 0; j < t->frontSupporters.size(); ++j)
		{
			LinkedTriple* s = t->frontSupporters[j];
			ex += s->prob * (s->fate.m_X - t->fate.m_X);
			ey += s->prob * (s->fate.m_Y - t->fate.m_Y);
			sump += s->prob;
		}
		if (sump > 0)
		{
			t->_fate0.m_X = ex / sump;
			t->_fate0.m_Y = ey / sump;
		}
		else
		{
			t->_fate0 = t->fate;
		}
	}
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* t = triples[i];
		t->fate = t->_fate0;
	}
}
*/

void support(vector<LinkedTriple*>& triples, float thres)
{
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* tr = triples[i];
		float sum1 = 0; // tr->fitness;
		for (int j = 0; j < tr->frontSupporters.size(); ++j)
		{
			LinkedTriple* tr2 = tr->frontSupporters[j];
			sum1 += tr2->prob * tr->compatibility0(tr2);
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
		tr->_totalFitness = tr->prob * tr->fitness + pow(sum1*sum2*sum3, 1.0 / 3.0);
		if (tr->id == 293)
		{
			i += 0;
		}
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
		if (tr->id == 293)
		{
			i += 0;
		}
	}
	for (int i = 0; i < triples.size(); ++i)
	{
		triples[i]->prob = triples[i]->_prob0;
		if (triples[i]->p->getId() == 5)
		{
			triples[i]->print();
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
	(LinkedTripleFactory::getInstance()).unit = calculateScaleUnit(points); //set the unit

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
	StationaryParticleFactory& spfactory = StationaryParticleFactory::getInstance();
	vector<LinkedTriple*> triples = collectTriples(points, thres1, 1.25f);
	assignSupporters(triples, thres2);
	for (int i = 0; i < numIter; ++i)
	{
		support(triples, 0.01);
		//updateFate(triples);
	}

	if (nlhs >= 1)
	{
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
			if (t->p->getId() == 25)
			{
				i += 0;
			}
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
	StationaryParticleFactory::getInstance().clean();
	LinkedTripleFactory::getInstance().clean();
	mexUnlock();
}

