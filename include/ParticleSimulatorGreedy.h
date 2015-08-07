#pragma once
#include <ParticleSimulator.h>

class ParticleSimulatorGreedy : public ParticleSimulator
{
public:
	ParticleSimulatorGreedy() : ParticleSimulator()
	{
	}
	virtual bool Simulate(float endtime = 10.0f, float delta = 0.1f, bool bdebug = false);
	pair<MovingParticle*, MovingParticle*> applyEventGreedy(EventStruct ev);
	float computeFitness(EventStruct ev);
	set<MovingParticle*> covered;
};