#include <ParticleSimulatorGreedyPartition.h>
#include <szmexutilitytemplate.h>

FitnessStruct
ParticleSimulatorGreedyPartition::computeFitness(EventStruct ev, bool left)
{
	FitnessStruct fs;
	MovingParticle* pe = (MovingParticle*)ev.p;
	MovingParticle* qe = (MovingParticle*)ev.q;
	MovingParticle* re = (MovingParticle*)ev.r;
	if (pe != NULL && qe != NULL && re != NULL)
	{
		//vector<MovingParticle*> vp1 = MovingParticle::extractPath(pe, qe);
		//vector<MovingParticle*> vp2 = MovingParticle::extractPath(re, pe);
		fs.bleft = true;
		fs.coverage = 1.0f;
		fs.fitness = 1.0/(1.0+ev.t); // (float)Min(vp1.size(), vp2.size()) / (float)Max(vp1.size(), vp2.size());
		fs.p = pe;
		fs.value = fs.fitness;
	}
	return fs;
}

FitnessStruct
ParticleSimulatorGreedyPartition::findNextEventGreedy(vector<MovingParticle*>& particles)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	FitnessStruct bestFitness;
	for (int i = 0; i < particles.size(); ++i)
	{
		MovingParticle* p = particles[i];
		if (p->isUnstable()) continue;
		if (p->isReflex())
		{
			p->setEvent(p->findNextSplitEvent(particles)); // p->getPolygon()->getParticles()));
			EventStruct ev = p->getEvent();
			if (ev.q == NULL || ev.r == NULL)
			{
				continue;
			}
			if (ev.p == ev.q->getPrev() || ev.p == ev.r->getNext())
			{
				continue; //q or r is adjacent to p => if we allow this, applyEventGreedy gets complicated.
			}
			if (ev.p->getInitParticle() == ev.q->getInitParticle() || ev.p->getInitParticle() == ev.r->getInitParticle())
			{
				continue;
			}
			FitnessStruct fitLeft = computeFitness(ev, true);
			FitnessStruct fitRight = computeFitness(ev, false);
			if (fitLeft.value > bestFitness.value)
			{
				bestFitness = fitLeft;
			}
			if (fitRight.value > bestFitness.value)
			{
				bestFitness = fitRight;
			}
		}
	}
	return bestFitness;
}


bool
ParticleSimulatorGreedyPartition::_PartitionRecusive(MovingParticle* p0, int level,
													float thres, int minLength, int maxLevel)
{
	if (p0->isActive() == false) return false;
	//if (level > maxLevel) return false;

	vector<MovingParticle*> particles = MovingParticle::vectorize(p0);
	Snapshot shot(level, 0.0f, particles);
	snapshots.push_back(shot);

	vector<MovingParticle*> area = extractSimplePath(p0, p0);

	FitnessStruct bestFitness = findNextEventGreedy(particles);
	bool bLeft = bestFitness.bleft;
	MovingParticle* p = bestFitness.p;
	float fvalue = bestFitness.value;
	
	if (fvalue < thres || fvalue != fvalue || particles.size() < minLength || level > maxLevel)
	{
		//printf("fvalue = %f. Done.\n", fvalue);
		if (area.size() >= 3)
		{
			Snapshot shot2(level, 0.0f, area);
			closedRegions.push_back(shot2);
		}
		return false;
	}
	else
	{
		//print the chosen event for debugging
		bestFitness.p->getEvent().print();
	}
	//apply the event
	pair<MovingParticle*, MovingParticle*> pp = applyEventGreedy(p->getEvent());
	MovingParticle* pnew[2] = { pp.first, pp.second };
	for (int i = 0; i < 2; ++i)
	{
		if (pnew[i]->isActive())
		{
			MovingParticle* ap = i == 0 ? pnew[i] : pnew[i]->getNext();
			printf("_PartitionRecursive: %d %d %d %f\n", closedRegions.size(), level, i, fvalue);
			if (_PartitionRecusive(ap, level + 1, thres, minLength, maxLevel) == false)
			{
			}
		}
	}
	return true;
}

bool
ParticleSimulatorGreedyPartition::Simulate(int maxLevel, float thres, int minLength)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	//int minLength = 3;
	//int maxLevel = std::numeric_limits<int>::max();
	_PartitionRecusive(*(factory.activeSet.begin()), 0, thres, minLength, maxLevel);
	return true;
}

