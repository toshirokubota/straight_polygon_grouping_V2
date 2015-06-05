#pragma once;
#include <ParticleSimulator.h>

class ParticleSimulatorTunneling: public ParticleSimulator
{
public:
	virtual bool Simulate(float endtime = 10.0f, float delta = 0.1f, bool bdebug = false);
};

