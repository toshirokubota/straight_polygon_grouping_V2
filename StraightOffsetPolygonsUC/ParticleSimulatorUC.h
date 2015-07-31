#pragma once;
#include <ParticleSimulator.h>

class ParticleSimulatorUC : public ParticleSimulator
{
public:
	ParticleSimulatorUC(float thres = 0.99) : ParticleSimulator()
	{
		this->thres = thres;
	}
	virtual bool Simulate(float endtime = 10.0f, float delta = 0.1f, bool bdebug = false);
	float thres;
};

