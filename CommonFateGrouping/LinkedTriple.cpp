#include <LinkedTriple.h>
#include <szMiscOperations.h>

LinkedTriple::LinkedTriple(StationaryParticle* p0, StationaryParticle* q0, StationaryParticle* r0, int id0)
{
	p = p0;
	q = q0;
	r = r0;
	calculateVelocity(p, q, r, v[0], v[1]);
	fitness = fitnessMeasure(p, q, r);
	id = id0;
}

float 
LinkedTriple::fitnessMeasure(StationaryParticle* p, StationaryParticle* q, StationaryParticle* r)
{
	LinkedTripleFactory& factory = LinkedTripleFactory::getInstance();
	/*
	float d1 = Distance(p->getP(), r->getP()) / factory.unit;
	float d2 = Distance(p->getP(), q->getP()) / factory.unit;
	float ang = GetVisualAngle(r->getX(), r->getY(), q->getX(), q->getY(), p->getX(), p->getY());
	float sn = sin(ang/2) + 1;
	return sn * sn * exp(-d1*d2/2.0);
	*/
	float d1 = Distance(p->getP(), r->getP());
	float d2 = Distance(p->getP(), q->getP());
	float df = (d1 - d2) / (factory.unit); 
	float ang = GetVisualAngle(r->getX(), r->getY(), q->getX(), q->getY(), p->getX(), p->getY());
	float sn = sin(ang / 2);
	return sn * sn * exp(-df * df / 2.0);
}

bool 
LinkedTriple::isNil(LinkedTriple* t)
{
	return t->q == t->p && t->r == t->p;
}


bool 
LinkedTriple::calculateVelocity(StationaryParticle* p, StationaryParticle* q, StationaryParticle* r, float& vx, float& vy)
{
	CParticleF b = bisector(p->getP(), r->getP(), q->getP());
	vx = b.m_X;
	vy = b.m_Y;
	return true;
}

void 
LinkedTriple::print()
{
	printf("%d %d %d %d %d %d %d %d %3.3f %3.3f %3.3f, %3.3f\n",
		id, p->getId(), r->getId(), q->getId(), 
		competitors.size(), frontSupporters.size(), leftSupporters.size(), rightSupporters.size(), 
		v[0], v[1], fitness, prob);
}

LinkedTriple* 
LinkedTriple::best()
{
	float maxFit = 0;
	LinkedTriple* b = NULL;
	for (int i = 0; i < frontSupporters.size(); ++i)
	{
		float f = compatibility(frontSupporters[i]) * frontSupporters[i]->prob;
		if (f > maxFit)
		{
			maxFit = f * frontSupporters[i]->prob;
			b = frontSupporters[i];
		}
	}
	return b;
}

float
LinkedTriple::_timeToClosestEncounter(CParticleF& p, float u[], CParticleF& q, float v[])
{
	float dx = p.m_X - q.m_X;
	float dy = p.m_Y - q.m_Y;
	float ux = u[0] - v[0];
	float uy = u[1] - v[1];

	float lu2 = ux * ux + uy * uy;
	float s = -(dx * ux + dy * uy) / lu2;
	return s;
}

float
LinkedTriple::_timeToClosestEncounter(LinkedTriple* t)
{
	return _timeToClosestEncounter(p->getP(), v, t->p->getP(), t->v);
}

float
LinkedTriple::compatibility(CParticleF& p, float* u, CParticleF& q, float* v)
{
	LinkedTripleFactory& factory = LinkedTripleFactory::getInstance();
	float s = _timeToClosestEncounter(p, u, q, v);
	if (s < 0) return 0.0f;
	if (s >= std::numeric_limits<float>::infinity() || s != s)
	{
		return 0.0f;
	}

	float x = p.m_X;
	float y = p.m_Y;
	float x2 = q.m_X;
	float y2 = q.m_Y;
	float x0 = x + u[0] * s;
	float y0 = y + u[1] * s;

	float dx = x - x2;
	float dy = y - y2;
	float ux = u[0] - v[0];
	float uy = u[1] - v[1];


	float dd = (dx + ux * s) * (dx + ux * s) + (dy + uy * s) * (dy + uy * s);
	float sgm = factory.unit;
	float ee = exp(-dd / (2 * sgm*sgm));

	float ang = GetVisualAngle(x, y, x2, y2, x0, y0);
	float sn = sin(ang / 2.0);
	return sn * sn  * ee;
}

float
LinkedTriple::compatibility(LinkedTriple* t)
{
	return compatibility(p->getP(), v, t->p->getP(), t->v);
}

bool 
LinkedTriple::updateFate()
{
	if (frontSupporters.size() == 0) return false;

	/*float sx = 0, sy = 0, st = 0;
	for (int i = 0; i < awaySupporters.size(); ++i)
	{
		LinkedTriple* t = awaySupporters[i];
		float c = compatibility(t);
		float s = _timeToClosestEncounter(t);
		sx += s * v[0] * t->prob;
		sy += s * v[1] * t->prob;
		st += t->prob;
	}
	fate = CParticleF(p->getX() + sx / st, p->getY() + sy / st);*/
	LinkedTriple* b = best();
	if (b == NULL) return false;

	float s = this->_timeToClosestEncounter(b);
	fate = CParticleF(p->getX() + s * v[0], p->getY() + s * v[1]); 
	return true;
}
