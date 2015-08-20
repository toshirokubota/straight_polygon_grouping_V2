#include <ParticleSimulatorGreedyPartition.h>

bool
ParticleSimulatorGreedyPartition::Simulate(float endtime, float delta, bool bdebug)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	PolygonFactory& pfactory = PolygonFactory::getInstance();
	float snapTime = delta;
	bool bSuccess = true;
	{
		vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
		snapshots.insert(snapshots.end(), shots.begin(), shots.end());
	}
	int minsize = 5;
	int iter = 0;
	while (true)
	{
		iter++;
		EventStruct selected;
		MovingParticle* psel = NULL;
		bool bLeft = true;
		int count = 0;
		for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
		{
			count++;
			MovingParticle* p = *it;
			if (p->isUnstable()) continue;
			if (p->isReflex() && p->getType() == Initial)
			{
				p->setEvent(p->findNextSplitEvent(p->getPolygon()->getParticles()));
				EventStruct ev = p->getEvent();
				if (selected.t > ev.t)
				{
					selected = ev;
				}
			}
		}
		MovingParticle* p = (MovingParticle*)selected.p;
		if (selected.t > endtime)
			break;
		p->getEvent().print();

		//apply the event
		int np = p->getPolygon()->getNumParticles();
		pair<MovingParticle*,MovingParticle*> pnew = applyEventGreedy(p->getEvent());
		{
			vector<Snapshot> shots = Snapshot::TakeSnapshot((float)iter);
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			snapTime += delta;
		}
		{
			MovingParticle* ps[] = { pnew.first, pnew.second };
			for (int j = 0; j < 2; ++j)
			{
				vector<MovingParticle*> area = extractSimplePath(ps[j], ps[j]);
				vector<MovingParticle*> vp = MovingParticle::vectorize(ps[j]);
				//if (area.size() >= minsize)
				{
					Snapshot shot((float)iter, 0.0f, vp);
					closedRegions.push_back(shot);
					/*printf("vp from pnew\n");
					for (int k = 0; k < vp.size(); ++k)
					{
						vp[k]->print();
					}
					//mexErrMsgTxt("premature exit.");
					printf("area from pnew\n");
					for (int k = 0; k < area.size(); ++k)
					{
						area[k]->print();
					}*/
					//mexErrMsgTxt("premature exit.");
				}

				if (vp.size() < minsize)
				{
					for (int k = 0; k < vp.size(); ++k)
					{
						factory.inactivate(vp[k]);
					}
				}
			}
		}

		//chosen.push_back(p->getPolygon());
		doneEvents.push_back(p->getEvent());
	}
	return bSuccess;
}

