#pragma once
#include <ParticleSimulatorGreedy.h>

class ParticleSimulatorGreedyPartition : public ParticleSimulatorGreedy
{
public:
	ParticleSimulatorGreedyPartition() : ParticleSimulatorGreedy()
	{
	}
	virtual bool Simulate(float endtime = 10.0f, float delta = 0.1f, bool bdebug = false);
};