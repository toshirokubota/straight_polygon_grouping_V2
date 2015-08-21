#include <ParticleSimulatorGreedy.h>
#include <szmexutilitytemplate.h>
#include <IntersectionConvexPolygons.h>

/*
Find a precise location where the edge should be split by a split event.
*/
CParticleF
ParticleSimulatorGreedy::makeSplitParticle(EventStruct ev)
{
	CParticleF c = ev.p->move(ev.t);
	CParticleF q = ev.q->move(ev.t);
	CParticleF r = ev.r->move(ev.t);
	pair<float, float> param = _IntersectConvexPolygon::intersect(ev.p->getP(), c, q, r);
	float t = param.second;
	CParticleF q0 = ev.q->getP();
	CParticleF r0 = ev.r->getP();
	CParticleF f((1 - t)*q0.m_X + t*r0.m_X, (1 - t)*q0.m_Y + t*r0.m_Y);
	return f;
}

/*
Find a precise location where the edge should be split by a split event.
*/
MovingParticle*
pickCollidingParticle(EventStruct ev)
{
	CParticleF c = ev.p->move(ev.t);
	CParticleF q = ev.q->move(ev.t);
	CParticleF r = ev.r->move(ev.t);
	pair<float, float> param = _IntersectConvexPolygon::intersect(ev.p->getP(), c, q, r);
	float t = param.second;

	if (t < 0.5) return (MovingParticle*) ev.q;
	else return (MovingParticle*)ev.r;
}

float
_fitnessMeasure(vector<CParticleF>& vp, float scale)
{
	float area = polygonArea(vp);
	float maxlen = 0;
	float sum = 0, sum2 = 0;
	for (int i = 0; i < vp.size(); ++i)
	{
		CParticleF p = vp[i];
		CParticleF q = vp[(i + 1) % vp.size()];
		//CParticleF r = vp[(i - 1 + vp.size()) % vp.size()];
		float d = Distance(p, q);
		maxlen = Max(d, maxlen);
		sum += d;
		sum2 += d * d;
	}
	//return sqrt(area) / maxlen;
	//return sqrt(area) / (scale + sum);
	//return area / (sum2);
	return sqrt(area) / scale;
}

/*
from a vector of particles (VP), find a simple path to pend, EXCLUDING pend.
*/
vector<MovingParticle*>
ParticleSimulatorGreedy::extractSimplePath(MovingParticle* p0, MovingParticle* pend)
{
	struct _Linker {
		MovingParticle* p;
		_Linker* prev;
		_Linker(MovingParticle* p=NULL)
		{
			this->p = p;
			prev = NULL;
		}
	};

	vector<MovingParticle*> area;
	int id = p0->getId(); //used as a unique number to detect looping

	vector<_Linker*> vl;
	_Linker* pl = new _Linker(p0);
	p0->getInitParticle()->link = (void*)pl;
	p0->getInitParticle()->info = id;
	vl.push_back(pl);
	MovingParticle* p = p0->getNext();
	while (p != pend)
	{
		int id0 = p->getInitParticle()->info;
		_Linker* pl2 = new _Linker(p);
		if (id0 == id) //a loop is found. Break it.
		{
			_Linker* pl3 = (_Linker*)p->getInitParticle()->link;
			pl2->prev = pl3->prev;
		}
		else
		{
			pl2->prev = pl;
		}

		vl.push_back(pl2);
		pl = pl2;
		p->getInitParticle()->link = (void*)pl;
		p->getInitParticle()->info = id;
		p = p->getNext();
	}

	//trace back
	while (pl != NULL)
	{
		area.insert(area.begin(), pl->p);
		pl = pl->prev;
	}
	//clean up
	for (int i = 0; i < vl.size(); ++i)
	{
		vl[i]->p->getInitParticle()->info = -1;
		vl[i]->p->getInitParticle()->link = NULL;

		delete vl[i];
		vl[i] = NULL;
	}
	vl.clear();

	return area;
}

float
ParticleSimulatorGreedy::_getCoverage(vector<MovingParticle*>& area)
{
	//return 1.0f; ///TK - just for temporary

	set<StationaryParticle*> pset;
	for (int i = 0; i < area.size(); ++i)
	{
		pset.insert(area[i]->getInitParticle());
	}
	set<StationaryParticle*> iset;
	set<StationaryParticle*> uset;
	set_intersection(pset.begin(), pset.end(),
		covered.begin(), covered.end(),
		std::inserter(iset, iset.begin()));
	/*set_union(pset.begin(), pset.end(),
		covered.begin(), covered.end(),
		std::inserter(uset, uset.begin()));*/
	//float coverage = 1 - (float)iset.size() / (float)uset.size();
	float coverage = 1 - (float)iset.size() / (float)pset.size();
	return coverage;
}

FitnessStruct
ParticleSimulatorGreedy::computeFitness(EventStruct ev, bool left)
{
	FitnessStruct fs;
	MovingParticle* pe = (MovingParticle*)ev.p;
	MovingParticle* qe = (MovingParticle*)ev.q;
	MovingParticle* re = (MovingParticle*)ev.r;
	if (pe != NULL && qe != NULL && re != NULL)
	{
		//CParticleF f = makeSplitParticle(ev);
		MovingParticle* f = pickCollidingParticle(ev);
		vector<MovingParticle*> area;
		if (left)
		{
			area = extractSimplePath(pe, qe);
			bool bFound = false;
			for (int k = 0; k < area.size(); ++k)
			{
				if (area[k]->getInitParticle() == pe->getInitParticle())
				{
					bFound = true;
					break;
				}
			}
			//for a path from PE to QE, it may be that PE is not a part of the simple path. In such case, the fitness is 0.
			if (!bFound)
			{
				area.clear();
			}
			else
			{
				area.push_back(qe);
				if (f != qe)
				{
					area.push_back(f);
				}
			}
		}
		else
		{
			area = extractSimplePath(re, pe);
			bool bFound = false;
			for (int k = 0; k < area.size(); ++k)
			{
				if (area[k]->getInitParticle() == pe->getInitParticle())
				{
					bFound = true;
					break;
				}
			}
			//for a path from RE to PE, it may be that PE is not a part of the simple path. In such case, the fitness is 0.
			//if (find(area.begin(), area.end(), pe) == area.end())
			if (!bFound)
			{
				area.clear();
			}
			else
			{
				if (f != re)
				{
					area.insert(area.begin(), f);
				}
				area.insert(area.begin(), pe);
			}
		}
		vector<CParticleF> vp0;
		for (int i = 0; i < area.size(); ++i)
		{
			vp0.push_back(area[i]->getP());
		}
		float coverage = _getCoverage(area);
		float scale = Max(1.0f, Distance(f->getP(), pe->getP()) * 2);
		float fit = _fitnessMeasure(vp0, scale);
		float value = coverage * fit;
		if (value > fs.value)
		{
			fs.value = value;
			fs.fitness = fit;
			fs.coverage = coverage;
			fs.f = f;
			fs.area = area;
			fs.bleft = left;
		}
	}
	return fs;
}

bool
ParticleSimulatorGreedy::Simulate(float endtime, float delta, bool bdebug)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	PolygonFactory& pfactory = PolygonFactory::getInstance();
	float snapTime = delta;
	bool bSuccess = true;
	{
		vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
		snapshots.insert(snapshots.end(), shots.begin(), shots.end());
	}
	float thres = delta;
	int iter = 0;
	while (iter < endtime)
	{
		iter++;
		FitnessStruct bestFitness;
		FitnessStruct bestFitness0;
		MovingParticle* psel = NULL;
		bool bLeft = true;
		int count = 0;
		for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
		{
			count++;
			MovingParticle* p = *it;
			if (p->isUnstable()) continue;
			if (p->isReflex())
			{
				p->setEvent(p->findNextSplitEvent(p->getPolygon()->getParticles()));
				EventStruct ev = p->getEvent();
				if (ev.q == NULL || ev.r == NULL)
				{
					continue;
				}
				if (ev.p->getInitParticle() == ev.q->getInitParticle() || ev.p->getInitParticle() == ev.r->getInitParticle())
				{
					continue;
				}
				FitnessStruct fitLeft = computeFitness(ev, true);
				FitnessStruct fitRight = computeFitness(ev, false);
				if (p->getId() == 89)
				{
					ev.print();
					printf("%f %f %f %f %f %f\n", 
						fitLeft.value, fitRight.value, fitLeft.fitness, fitRight.fitness, fitLeft.coverage, fitRight.coverage);
				}
				if (fitLeft.value > bestFitness.value)
				{
					bestFitness = fitLeft;
					bestFitness0 = fitRight;
					psel = p;
					bLeft = true;
				}
				if (fitRight.value > bestFitness.value)
				{
					bestFitness = fitRight;
					bestFitness0 = fitLeft;
					psel = p;
					bLeft = false;
				}
			}
		}
		MovingParticle* p = psel;
		float fvalue = bestFitness.value;
		if (fvalue <= thres || fvalue != fvalue)
			break;

		//apply the event
		//int np = p->getPolygon()->getNumParticles();
		vector<MovingParticle*> vp0 = MovingParticle::vectorize(p);
		pair<MovingParticle*,MovingParticle*> pnew = applyEventGreedy(p->getEvent());
		{
			vector<Snapshot> shots = Snapshot::TakeSnapshot((float)iter);
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			snapTime += delta;
		}
		p->getEvent().print();

		{
			MovingParticle* ap;
			MovingParticle* ap2;
			if (bLeft)
			{
				ap = pnew.first;
				ap2 = pnew.second;
			}
			else
			{
				ap = pnew.second;
				ap2 = pnew.first;
			}
			vector<MovingParticle*> area = extractSimplePath(ap, ap->getPrev());
			area.push_back(ap->getPrev());
			vector<MovingParticle*> area2 = extractSimplePath(ap2, ap2->getPrev());
			area2.push_back(ap2->getPrev());
			Snapshot shot((float)iter, 0.0f, area);
			closedRegions.push_back(shot);
			chosen.push_back(pfactory.makePolygon(area, (float)iter));
			for (int j = 0; j < area.size(); ++j)
			{
				covered.insert(area[j]->getInitParticle());
			}
			printf("%d: fitness = %f, %f, %f, %d, %d %s\n", 
				iter, bestFitness.value, bestFitness.fitness, bestFitness.coverage, bestFitness.area.size(), area.size(),
				bestFitness.bleft ? "Left": "Right");
			int dc = bestFitness.area.size() - area.size();
			if (Abs(dc)>=2)
			{
				printf("New particles:\n");
				pnew.first->print();
				pnew.first->getPrev()->print();
				pnew.second->print();
				pnew.second->getNext()->print();
				printf("Best fitness:\n");
				for (int j = 0; j < bestFitness.area.size(); ++j)
				{
					bestFitness.area[j]->print();
				}
				printf("Best fitness0:\n");
				for (int j = 0; j < bestFitness0.area.size(); ++j)
				{
					bestFitness0.area[j]->print();
				}
				printf("Area 1:\n");
				for (int j = 0; j < area.size(); ++j)
				{
					area[j]->print();
				}
				printf("Area 2:\n");
				for (int j = 0; j < area2.size(); ++j)
				{
					area2[j]->print();
				}

				vector<MovingParticle*> vp1 = MovingParticle::vectorize(ap);
				vector<MovingParticle*> vp2 = MovingParticle::vectorize(ap2);
				printf("Trace 0:\n");
				for (int j = 0; j < vp0.size(); ++j)
				{
					vp0[j]->print();
				}
				printf("Trace 1:\n");
				for (int j = 0; j < vp1.size(); ++j)
				{
					vp1[j]->print();
				}
				printf("Trace 2:\n");
				for (int j = 0; j < vp2.size(); ++j)
				{
					vp2[j]->print();
				}
				mexErrMsgTxt("Premature exist.");
			}
		}

		//chosen.push_back(p->getPolygon());
		doneEvents.push_back(p->getEvent());
	}
	return bSuccess;
}

pair<MovingParticle*,MovingParticle*>
ParticleSimulatorGreedy::applyEventGreedy(EventStruct ev)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	StationaryParticleFactory& sfactory = StationaryParticleFactory::getInstance();
	pair<MovingParticle*, MovingParticle*> particles(NULL, NULL);
	if (ev.type == SplitEvent || ev.type == CollisionEvent)
	{
		MovingParticle* p = (MovingParticle*)ev.p;
		MovingParticle* q = (MovingParticle*)ev.q;
		MovingParticle* r = (MovingParticle*)ev.r;
		MovingParticle* pn = p->getNext();
		MovingParticle* pr = p->getPrev();
		MovingParticle* sp = pickCollidingParticle(ev); // sfactory.makeParticle(makeSplitParticle(ev));
		{ //from q to p
			MovingParticle* pnew = factory.makeParticle(p->getInitParticle(), Split, 0.0f);
			MovingParticle* qnew = factory.makeParticle(sp->getInitParticle(), Split, 0.0f);
			ParticleDirection pd(sp->getP(), p->getP(), true);
			MovingFront mf(pd);
			MovingParticle::setNeighbors(qnew, q, pnew, sp==q ? q->getRear(): q->getFront(), mf);
			MovingParticle::setNeighbors(pnew, qnew, pn, mf, p->getFront());
			bool bUnstable = false;
			if (qnew->calculateVelocityR() == false) bUnstable = true;
			if (pnew->calculateVelocityR() == false) bUnstable = true;
			if (bUnstable)
			{
				MovingParticle* ptmp = pnew;
				pnew = MovingParticle::traceAndHandleUnstable(pnew, p);
				factory.inactivate(ptmp);
			}
			MovingParticle::updatePolygon(pnew);
			particles.first = pnew;
			pnew->print("pnew:");
			qnew->print("qnew:");
		}
		{ //from p to r
			MovingParticle* pnew = factory.makeParticle(p->getInitParticle(), Split, 0.0f);
			MovingParticle* qnew = factory.makeParticle(sp->getInitParticle(), Split, 0.0f);
			ParticleDirection pd(p->getP(), sp->getP(), true);
			MovingFront mf(pd);
			MovingParticle::setNeighbors(qnew, pnew, r, mf, sp==r ? r->getFront(): r->getRear());
			MovingParticle::setNeighbors(pnew, pr, qnew, p->getRear(), mf);
			bool bUnstable = false;
			if (qnew->calculateVelocityR() == false) bUnstable = true;
			if (pnew->calculateVelocityR() == false) bUnstable = true;
			if (bUnstable)
			{
				MovingParticle* ptmp = pnew;
				pnew = MovingParticle::traceAndHandleUnstable(pnew, p);
				factory.inactivate(ptmp);
			}
			MovingParticle::updatePolygon(pnew);
			particles.second = pnew;
			pnew->print("pnew:");
			qnew->print("qnew:");
		}
		factory.inactivate(p);
		if (sp == q)
		{
			factory.inactivate(q);
		}
		else
		{
			factory.inactivate(r);
		}
	}
	/*vector<MovingParticle*> vp1 = extractSimplePath(particles.first, particles.first);
	vector<MovingParticle*> vp2 = extractSimplePath(particles.second, particles.second);
	printf("# particles in resulting areas: %d %d\n", vp1.size(), vp2.size());*/
	return particles;
}

