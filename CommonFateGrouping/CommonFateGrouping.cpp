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

bool
isNil(LinkedTriple* t)
{
	return t->q == t->p && t->r == t->p;
}

vector<LinkedTriple*>
collectTriples(vector<CParticleF>& points, float thres, float scale)
{
	Triangulation::Triangulator trmap(points);

	StationaryParticleFactory& factory = StationaryParticleFactory::getInstance();
	LinkedTripleFactory& tfactory = LinkedTripleFactory::getInstance();
	vector<StationaryParticle*> particles;
	float separation = 0;
	float margin = 1.0f;
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* p = trmap.points[i];
		StationaryParticle* sp = factory.makeParticle(p->p);
		particles.push_back(sp);
		vector<float> vlen;
		for (int j = 0; j < p->edges.size(); ++j)
		{
			vlen.push_back(p->edges[j]->Length());
		}
		sort(vlen.begin(), vlen.end());
		if (vlen[1] + margin> separation)
		{
			separation = vlen[1] + margin;
		}
	}

	vector<LinkedTriple*> triples;
	vector<LinkedTriple*> nils;
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* tp = trmap.points[i];
		/*vector<pair<float, Triangulation::_Internal::_vertex*>> pairs;
		for (int j = 0; j < tp->edges.size(); ++j)
		{
			Triangulation::_Internal::_edge* ed = tp->edges[j];
			Triangulation::_Internal::_vertex* v = ed->vertices[0] == tp ? ed->vertices[1] : ed->vertices[0];
			pairs.push_back(pair<float, Triangulation::_Internal::_vertex*>(ed->Length(), v));
		}
		sort(pairs.begin(), pairs.end());
		float separation = pairs[1].first * scale; //scale times the second smallest length*/

		CParticleF p = trmap.points[i]->p;
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
				float ang = GetVisualAngle(r.m_X, r.m_Y, q.m_X, q.m_Y, p.m_X, p.m_Y);
				if (ang > thres)
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

void
assignSupporters(vector<LinkedTriple*>& triples, float thres)
{
	float toobig = 1.0e4;
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* t1 = triples[i];
		if (isNil(t1)) continue; //this is a NIL

		CParticleF p1 = t1->p->getP();
		CParticleF p2(p1.m_X + t1->v[0], p1.m_Y + t1->v[1]);
		for (int j = 0; j < triples.size(); ++j)
		{
			if (i == j) continue;
			LinkedTriple* t2 = triples[j];
			if (isNil(t2)) continue;

			CParticleF q1 = t2->p->getP();
			CParticleF q2(q1.m_X + t2->v[0], q1.m_Y + t2->v[1]);
			pair<float, float>  param = _IntersectConvexPolygon::intersect(p1, p2, q1, q2);
			float df = Abs(param.first - param.second);
			if (param.first>0 && param.second > 0 && param.first < toobig && param.second<toobig && df < thres)
			{
				float fitness = exp(-df * df / (2.0 * thres/2 * thres/2));
				t1->supporters.push_back(t2);
				t1->compatibility.push_back(fitness);
				t1->linkWeights.push_back(1.0);
			}
		}
	}
}

void support(vector<LinkedTriple*>& triples, float thres)
{
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* tr = triples[i];
		tr->_totalFitness = tr->fitness;
		for (int j = 0; j < tr->supporters.size(); ++j)
		{
			LinkedTriple* tr2 = tr->supporters[j];
			tr->_totalFitness += tr2->prob * tr->compatibility[j] * tr->linkWeights[j];
			/*int k = distance(tr2->supporters.begin(), find(tr2->supporters.begin(), tr2->supporters.end(), tr));
			if (k < tr2->supporters.size())
			{
				tr->_totalFitness += tr2->prob * tr->compatibility[j] * tr->linkWeights[j] * tr2->linkWeights[k];
			}*/
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
	}
	for (int i = 0; i < triples.size(); ++i)
	{
		triples[i]->prob = triples[i]->_prob0;
	}
	/*for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* t = triples[i];
		t->linkWeights0.clear();
		for (int j = 0; j < t->supporters.size(); ++j)
		{
			LinkedTriple* t2 = t->supporters[j];
			int k = distance(t2->supporters.begin(), find(t2->supporters.begin(), t2->supporters.end(), t));
			if (k >= t2->supporters.size())
			{
				t->linkWeights[j] = 0.0f;
			}
			else
			{
				float dd = t->compatibility[j];
				float lnk1 = t->linkWeights[j];
				float lnk2 = t2->linkWeights[k];
				float pp = t2->prob;
				float ee = pp*lnk1*lnk2;
				float pp2 = ee*dd / (ee*dd + (1. - ee)*thres);
				t->linkWeights0.push_back(pp2);
			}
		}
	}*/
	/*for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* t = triples[i];
		t->linkWeights = t->linkWeights0;
	}*/
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
	float maxSeparation = 50.0f;
	float thres1 = PI / 2.0f;
	float thres2 = 5.0f;
	int numIter = 10;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres1, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[2], classMode);
	}
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(maxSeparation, prhs[3], classMode);
	}

	StationaryParticleFactory& spfactory = StationaryParticleFactory::getInstance();
	vector<LinkedTriple*> triples = collectTriples(points, thres1, 2.0f);
	assignSupporters(triples, thres2);
	for (int i = 0; i < numIter; ++i)
	{
		support(triples, 0.01);
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
		const int dims[] = { triples.size(), 10 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], (float)triples[i]->id);
			SetData2(F, i, 1, dims[0], dims[1], (float)triples[i]->p->getId());
			SetData2(F, i, 2, dims[0], dims[1], (float)triples[i]->r->getId());
			SetData2(F, i, 3, dims[0], dims[1], (float)triples[i]->q->getId());
			CParticleF e = triples[i]->estimate();
			SetData2(F, i, 4, dims[0], dims[1], e.m_X);
			SetData2(F, i, 5, dims[0], dims[1], e.m_Y);
			SetData2(F, i, 6, dims[0], dims[1], triples[i]->v[0]);
			SetData2(F, i, 7, dims[0], dims[1], triples[i]->v[1]);
			SetData2(F, i, 8, dims[0], dims[1], triples[i]->fitness);
			SetData2(F, i, 9, dims[0], dims[1], triples[i]->prob);
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

