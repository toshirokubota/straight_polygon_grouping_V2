#include <ParticleSimulatorUC.h>

/*
Tunneling simulation only considers  edge events (i.e. collision of adjacent particles)
*/
bool
ParticleSimulatorUC::Simulate(float endtime, float delta, bool bdebug)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
	{
		MovingParticle* p = *it;
		if (p->isLeaf()) continue;
		for (set<MovingParticle*>::iterator it2 = factory.activeSet.begin(); it2 != factory.activeSet.end(); ++it2)
		{
			MovingParticle* q = *it2;
			if (q->isLeaf()) continue;
			MovingParticle* r = q->getNext();
			if (p == q || p->getNext() == q || p == r || p->getPrev() == r) continue;
			float t0 = MovingParticle::intersectSideAndLine(p, q, r);
			if (t0 > 0 && t0 < std::numeric_limits<float>::infinity())
			{
				//keep if the intersection point is closer than the pair
				CParticleF x = p->move(t0);
				float dpq = Distance(p->getP(), q->getP());
				float dpx = Distance(p->getP(), x);
				float dqx = Distance(q->getP(), x);
				float ratio = dpq / (dpx + dqx);
			}
		}
	}
	return true;
}
