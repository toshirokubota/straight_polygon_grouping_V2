#pragma once

#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <mex.h>
#include <szMexUtility.h>
#include <StationaryParticle.h>
#include <IntersectionConvexPolygons.h>
#include <MiscGeometry.h>

class LinkedTriple
{
public:
	LinkedTriple(StationaryParticle* p0, StationaryParticle* q0, StationaryParticle* r0, int id0)
	{
		p = p0;
		q = q0;
		r = r0;
		_setVelocity();
		fitness = fitnessMeasure();
		id = id0;
	}
	float fitnessMeasure()
	{
		float d1 = Distance(p->getP(), r->getP());
		float d2 = Distance(p->getP(), q->getP());
		return 10.0 / (d1 + d2);
	}
	bool _setVelocity()
	{
		CParticleF b = bisector(p->getP(), r->getP(), q->getP());
		v[0] = b.m_X;
		v[1] = b.m_Y;
		return true;
	}
	void print()
	{
		printf("%d (%3.3f, %3.3f) - (%3.3f, %3.3f) - (%3.3f, %3.3f) : (%3.3f, %3.3f), %f, %f\n",
			id, r->getX(), r->getY(), p->getX(), p->getY(), q->getX(), q->getY(), v[0], v[1], fitness, prob);
	}

	CParticleF estimate()
	{
		float maxFit = 0;
		int maxIndex = -1;
		CParticleF c;
		for (int i = 0; i < supporters.size(); ++i)
		{
			float f = compatibility[i] * supporters[i]->prob;
			if (f > maxFit)
			{
				maxFit = f;
				maxIndex = i;
			}
		}
		if (maxIndex >= 0)
		{
			LinkedTriple* t = supporters[maxIndex];
			CParticleF p1 = p->getP();
			CParticleF p2(p1.m_X + v[0], p1.m_Y + v[1]);
			CParticleF q1 = t->p->getP();
			CParticleF q2(q1.m_X + t->v[0], q1.m_Y + t->v[1]);
			pair<float, float>  param = _IntersectConvexPolygon::intersect(p1, p2, q1, q2);
			c.m_X = p1.m_X + v[0] * param.first;
			c.m_Y = p1.m_Y + v[1] * param.first;
		}
		else
		{
			c = p->getP();
		}
		return c;
	}
	LinkedTriple* best()
	{
		float maxFit = 0;
		LinkedTriple* b = NULL;
		for (int i = 0; i < supporters.size(); ++i)
		{
			float f = compatibility[i] * supporters[i]->prob;
			if (f > maxFit)
			{
				maxFit = f * supporters[i]->prob;
				b = supporters[i];
			}
		}
		return b;
	}
	StationaryParticle* p; //center
	StationaryParticle* q; //right
	StationaryParticle* r; //left
	float v[2];
	float fitness;
	vector<LinkedTriple*> supporters;
	vector<LinkedTriple*> competitors;
	vector<float> compatibility; //fitenss of each supporter
	vector<float> linkWeights; //link weight to each supporter
	vector<float> linkWeights0; 
	float prob;
	float _prob0; //temporary storage - this will be the prob in the next period
	float _totalFitness;

	int id;
};

struct LinkedTripleFactory
{
public:
	static LinkedTripleFactory& getInstance()
	{
		static LinkedTripleFactory instance;
		return instance;
	}
	LinkedTriple* makeTriple(StationaryParticle* sp, StationaryParticle* sr, StationaryParticle*  sq)
	{
		LinkedTriple* triple = new LinkedTriple(sp, sq, sr, _id++);
		if (_id == 605)
			_id += 0;
		triples.push_back(triple);
		return triple;
	}
	LinkedTriple* makeNillTriple(StationaryParticle* sp)
	{
		const float EPSILON = 1.0e-5;
		LinkedTriple* t = new LinkedTriple(sp, sp, sp, _id++);
		t->supporters.push_back(t);
		t->compatibility.push_back(EPSILON);
		t->linkWeights.push_back(1.0);
		t->v[0] = t->v[1] = 0.0f;
		t->fitness = 0;

		return t;
	}
	void clean()
	{
		for (int i = 0; i<triples.size(); ++i)
		{
			delete triples[i];
		}
		triples.clear();
		_id = 0;
	}
	vector<LinkedTriple*> triples;
private:
	int _id;
	LinkedTripleFactory()
	{
		_id = 0;
	}
	~LinkedTripleFactory()
	{
		clean();
	}
	LinkedTripleFactory(LinkedTripleFactory& f){}
	LinkedTripleFactory operator=(LinkedTripleFactory& f){}
};




