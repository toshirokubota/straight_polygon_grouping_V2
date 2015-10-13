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
	LinkedTriple(StationaryParticle* p0, StationaryParticle* q0, StationaryParticle* r0)
	{
		p = p0;
		q = q0;
		r = r0;
		_setVelocity();
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
		printf("(%3.3f, %3.3f) - (%3.3f, %3.3f) - (%3.3f, %3.3f) : (%3.3f, %3.3f)\n",
			r->getX(), r->getY(), p->getX(), p->getY(), q->getX(), q->getY(), v[0], v[1]);
	}

	CParticleF estimate()
	{
		float maxFit = 0;
		int maxIndex = -1;
		CParticleF c;
		for (int i = 0; i < supporters.size(); ++i)
		{
			float f = fitness[i];
			if (f > maxFit)
			{
				maxFit = f;
				maxIndex = i;
			}
		}
		if (maxIndex >= 0)
		{
			CParticleF p1 = p->getP();
			CParticleF p2(p1.m_X + v[0], p1.m_Y + v[1]);
			CParticleF q1 = supporters[maxIndex]->p->getP();
			CParticleF q2(q1.m_X + supporters[maxIndex]->v[0], q1.m_Y + supporters[maxIndex]->v[1]);
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
			float f = fitness[i] * supporters[i]->prob;
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
	vector<LinkedTriple*> supporters;
	vector<LinkedTriple*> competitors;
	vector<float> fitness; //fitenss of each supporter
	float prob;
	float _prob0; //temporary storage - this will be the prob in the next period
	float _totalFitness;
	static const float EPSILON;
};

LinkedTriple* makeNillTriple(StationaryParticle* sp)
{
	LinkedTriple* t = new LinkedTriple(sp, sp, sp);
	t->supporters.push_back(t);
	t->fitness.push_back(LinkedTriple::EPSILON);
	t->v[0] = t->v[1] = 0.0f;
}