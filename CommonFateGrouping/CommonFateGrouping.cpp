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


vector<LinkedTriple*> 
collectTriples(vector<CParticleF>& points, float thres, float separation)
{
	StationaryParticleFactory& factory = StationaryParticleFactory::getInstance();
	vector<StationaryParticle*> particles;
	for (int i = 0; i < points.size(); ++i)
	{
		CParticleF p = points[i];
		particles.push_back(factory.makeParticle(p));
	}

	vector<LinkedTriple*> triples;
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
				float ang = GetVisualAngle(r.m_X, r.m_Y, q.m_X, q.m_Y, p.m_X, p.m_Y);
				if (ang > thres)
				{
					triples0.push_back(new LinkedTriple(sp, sq, sr));
				}
			}
		}

		//add them as competitors
		for (int j = 0; j < triples0.size(); ++j)
		{
			for (int k = 0; k < triples0.size(); ++k)
			{
				if (j == k) continue;
				triples0[j]->competitors.push_back(triples0[k]);
			}
			triples0[j]->prob = 1.0 / (float)(triples0.size());
		}
		triples.insert(triples.end(), triples0.begin(), triples0.end());
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
		CParticleF p1 = t1->p->getP();
		CParticleF p2(p1.m_X + t1->v[0], p1.m_Y + t1->v[1]);
		for (int j = 0; j < triples.size(); ++j)
		{
			if (i == j) continue;
			LinkedTriple* t2 = triples[j];
			CParticleF q1 = t2->p->getP();
			CParticleF q2(q1.m_X + t2->v[0], q1.m_Y + t2->v[1]);
			pair<float, float>  param = _IntersectConvexPolygon::intersect(p1, p2, q1, q2);
			float df = Abs(param.first - param.second);
			if (param.first>0 && param.second > 0 && param.first < toobig && param.second<toobig && df < thres)
			{
				float fitness = exp(-df * df / (2.0 * thres/2 * thres/2));
				t1->supporters.push_back(t2);
				t1->fitness.push_back(fitness);
			}
			/*if (i == 1 && (j==15 || j==17))
			{
				printf("%d: (%f, %f)-(%f,%f)-(%f,%f): (%f, %f), (%f, %f), %f, %f\n",
					j, t2->r->getX(), t2->r->getY(), t2->p->getX(), t2->p->getY(), t2->q->getX(), t2->q->getY(), 
					t1->v[0], t1->v[1], t2->v[0], t2->v[1], param.first, df);
			}*/
		}
	}
}

void support(vector<LinkedTriple*>& triples)
{
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* tr = triples[i];
		tr->_totalFitness = 0;
		for (int j = 0; j < tr->supporters.size(); ++j)
		{
			LinkedTriple* tr2 = tr->supporters[j];
			tr->_totalFitness += tr2->prob * tr->fitness[j];
		}
	}
	for (int i = 0; i < triples.size(); ++i)
	{
		LinkedTriple* tr = triples[i];
		float sum = tr->prob * tr->_totalFitness;
		for (int j = 0; j < tr->competitors.size(); ++j)
		{
			LinkedTriple* tr2 = tr->competitors[j];
			sum += tr2->prob * tr->_totalFitness;
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
	float thres2 = 30.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres1, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(maxSeparation, prhs[2], classMode);
	}

	vector<LinkedTriple*> triples = collectTriples(points, thres1, maxSeparation);
	assignSupporters(triples, thres2);
	for (int i = 0; i < 5; ++i)
	{
		support(triples);
	}

	if (nlhs >= 1)
	{
		const int dims[] = { triples.size(), 11 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], triples[i]->p->getX());
			SetData2(F, i, 1, dims[0], dims[1], triples[i]->p->getY());
			SetData2(F, i, 2, dims[0], dims[1], triples[i]->r->getX());
			SetData2(F, i, 3, dims[0], dims[1], triples[i]->r->getY());
			SetData2(F, i, 4, dims[0], dims[1], triples[i]->q->getX());
			SetData2(F, i, 5, dims[0], dims[1], triples[i]->q->getY());
			CParticleF e = triples[i]->estimate();
			SetData2(F, i, 6, dims[0], dims[1], e.m_X);
			SetData2(F, i, 7, dims[0], dims[1], e.m_Y);
			SetData2(F, i, 8, dims[0], dims[1], triples[i]->v[0]);
			SetData2(F, i, 9, dims[0], dims[1], triples[i]->v[1]);
			SetData2(F, i, 10, dims[0], dims[1], triples[i]->prob);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	if (nlhs >= 2)
	{
		const int dims[] = { triples.size(), 11 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			LinkedTriple* t = triples[i]->best();
			if (t == NULL)
			{
				for (int j = 0; j < dims[1]; ++j)
				{
					SetData2(F, i, j, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
				}
			}
			else {
				SetData2(F, i, 0, dims[0], dims[1], t->p->getX());
				SetData2(F, i, 1, dims[0], dims[1], t->p->getY());
				SetData2(F, i, 2, dims[0], dims[1], t->r->getX());
				SetData2(F, i, 3, dims[0], dims[1], t->r->getY());
				SetData2(F, i, 4, dims[0], dims[1], t->q->getX());
				SetData2(F, i, 5, dims[0], dims[1], t->q->getY());
				CParticleF e = t->estimate();
				SetData2(F, i, 6, dims[0], dims[1], e.m_X);
				SetData2(F, i, 7, dims[0], dims[1], e.m_Y);
				SetData2(F, i, 8, dims[0], dims[1], t->v[0]);
				SetData2(F, i, 9, dims[0], dims[1], t->v[1]);
				SetData2(F, i, 10, dims[0], dims[1], t->prob);
			}
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	mexUnlock();
}

