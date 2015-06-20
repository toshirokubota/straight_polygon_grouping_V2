//#include <EventStruct.h>
#include <MovingParticle.h>
#include <mex.h>

EventType int2EventType(int i)
{
	EventType type = UnknownEvent;
	switch (i)
	{
	case 0:
		break;
	case 1:
		type = EdgeEvent;
		break;
	case 2:
		type = SplitEvent;
		break;
	case 3:
		type = CollisionEvent;
		break;
	default:
		break;
	}
	return type;
}


EventStruct::EventStruct(float t, EventType type, const MovingParticle* p, const MovingParticle* q)
{
	this->t = t;
	this->type = type;
	this->p = p;
	this->q = q;
	if (q == NULL)
	{
		r = NULL;
	}
	else
	{
		r = q->getNext();
	}
}

void EventStruct::print()
{
	CParticleF pp = p->getP();
	if (type == SplitEvent)
	{
		CParticleF qp = q->getP();
		CParticleF rp = r->getP();
		printf("event>> Split @ %f: %d(%3.3f,%3.3f)->%d(%3.3f,%3.3f), %d(%3.3f,%3.3f)\n",
			t, p->getId(), pp.m_X, pp.m_Y, q->getId(), qp.m_X, qp.m_Y, r->getId(), rp.m_X, rp.m_Y);
	}
	else if (type == EdgeEvent)
	{
		CParticleF qp = q->getP();
		printf("event>> Edge @ %f: %d(%3.3f,%3.3f)->%d(%3.3f,%3.3f)\n",
			t, p->getId(), pp.m_X, pp.m_Y, q->getId(), qp.m_X, qp.m_Y);
	}
	else if (type == CollisionEvent)
	{
		CParticleF qp = q->getP();
		CParticleF rp = r->getP();
		printf("event>> Collision @ %f: %d(%3.3f,%3.3f)->%d(%3.3f,%3.3f), %d(%3.3f,%3.3f)\n",
			t, p->getId(), pp.m_X, pp.m_Y, q->getId(), qp.m_X, qp.m_Y, r->getId(), rp.m_X, rp.m_Y);
	}
	else
	{
		printf("event>> Unknown @ %f: %d(%3.3f,%3.3f)\n",
			t, p->getId(), pp.m_X, pp.m_Y);
	}
}
