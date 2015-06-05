#include <ParticleSimulatorTunneling.h>
#include <map>
#include <szMiscOperations.h>
#include <MiscGeometry.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <szmexutilitytemplate.h>

/*
Tunneling simulation only considers  edge events (i.e. collision of adjacent particles)
*/
bool
ParticleSimulatorTunneling::Simulate(float endtime, float delta, bool bdebug)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	float snapTime = delta;
	bool bSuccess = true;
	for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
	{
		EventStruct ev = (*it)->findNextEdgeEvent();
		(*it)->event = ev;
	}
	{
		vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
		snapshots.insert(snapshots.end(), shots.begin(), shots.end());
	}
	int iter = 0;
	while (time < endtime)
	{
		iter++;
		vector<Snapshot> shots = Snapshot::TakeSnapshot(time); //temporary
		MovingParticle* p = MovingParticle::getNextEvent();
		if (p == NULL) break;
		if (p->getEvent().t > endtime) break;

		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			(*it)->update(p->getEvent().t - time);
		}
		time = p->getEvent().t;
		if (p->applyEvent() == false) break;

		MovingParticle::removeUnstable();
		if (time >= snapTime)
		{
			vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			snapTime += delta;
		}
		MovingParticle::quickFinish();

		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			EventStruct ev = (*it)->findNextEdgeEvent();
			(*it)->event = ev;
		}
	}
	return bSuccess;
}
