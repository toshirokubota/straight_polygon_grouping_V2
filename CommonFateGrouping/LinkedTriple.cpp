#include <LinkedTriple.h>
#include <szMiscOperations.h>

LinkedTriple::LinkedTriple(StationaryParticle* p0, StationaryParticle* q0, StationaryParticle* r0, int id0)
{
	p = p0;
	q = q0;
	r = r0;
	_setVelocity();
	fitness = fitnessMeasure(p, q, r);
	id = id0;
}

float 
LinkedTriple::fitnessMeasure(StationaryParticle* p, StationaryParticle* q, StationaryParticle* r)
{
	LinkedTripleFactory& factory = LinkedTripleFactory::getInstance();
	float d1 = Distance(p->getP(), r->getP()) / factory.unit;
	float d2 = Distance(p->getP(), q->getP()) / factory.unit;
	float ang = GetVisualAngle(r->getX(), r->getY(), q->getX(), q->getY(), p->getX(), p->getY());
	float sn = sin(ang/2) + 1;
	return sn * sn * exp(-d1*d2/2.0);
}

bool 
LinkedTriple::isNil(LinkedTriple* t)
{
	return t->q == t->p && t->r == t->p;
}


bool 
LinkedTriple::_setVelocity()
{
	CParticleF b = bisector(p->getP(), r->getP(), q->getP());
	v[0] = b.m_X;
	v[1] = b.m_Y;
	return true;
}
void 
LinkedTriple::print()
{
	printf("%d (%3.3f, %3.3f) - (%3.3f, %3.3f) - (%3.3f, %3.3f) : (%3.3f, %3.3f), %f, %f\n",
		id, r->getX(), r->getY(), p->getX(), p->getY(), q->getX(), q->getY(), v[0], v[1], fitness, prob);
}

LinkedTriple* 
LinkedTriple::best()
{
	float maxFit = 0;
	LinkedTriple* b = NULL;
	for (int i = 0; i < supporters.size(); ++i)
	{
		float f = compatibility(supporters[i]) * supporters[i]->prob;
		if (f > maxFit)
		{
			maxFit = f * supporters[i]->prob;
			b = supporters[i];
		}
	}
	return b;
}

float
LinkedTriple::_timeToClosestEncounter(LinkedTriple* t)
{
	float dx = p->getX() - t->p->getX();
	float dy = p->getY() - t->p->getY();
	float ux = v[0] - t->v[0];
	float uy = v[1] - t->v[1];

	float lu2 = ux * ux + uy * uy;
	float s = -(dx * ux + dy * uy) / lu2;
	if (Abs(s) > 1000)
	{
		s += 0;
	}
	return s;
}

float
LinkedTriple::compatibility(LinkedTriple* t)
{
	LinkedTripleFactory& factory = LinkedTripleFactory::getInstance();
	if (isNil(t)) return 0.0f;
	float s = _timeToClosestEncounter(t);
	if (s < 0) return 0.0f;

	float x = p->getX();
	float y = p->getY();
	float x2 = t->p->getX();
	float y2 = t->p->getY();
	float x0 = x + v[0] * s;
	float y0 = y + v[1] * s;

	float dx = x - x2;
	float dy = y - y2;
	float ux = v[0] - t->v[0];
	float uy = v[1] - t->v[1];


	float dd = (dx + ux * s) * (dx + ux * s) + (dy + uy * s) + (dy + uy * s);
	float sgm = factory.unit;
	float ee = exp(-dd / (2 * sgm*sgm));

	float ang = GetVisualAngle(x, y, x2, y2, x0, y0);
	float sn = sin(ang / 2.0);
	return sn * sn  * ee;
}

bool 
LinkedTriple::updateFate()
{
	if (supporters.size() == 0) return false;

	/*float sx = 0, sy = 0, st = 0;
	for (int i = 0; i < supporters.size(); ++i)
	{
		LinkedTriple* t = supporters[i];
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
