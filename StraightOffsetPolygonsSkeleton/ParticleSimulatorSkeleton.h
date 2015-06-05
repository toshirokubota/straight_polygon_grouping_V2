#pragma once;
#include <ParticleSimulator.h>

struct ParticleSkeleton
{
public:
	ParticleSkeleton(MovingParticle* p, MovingParticle* q, CParticleF& o)
	{
		this->p = p;
		this->q = q;
		this->o = o;
	}
	MovingParticle* p;
	MovingParticle* q;
	CParticleF o;
};

class ParticleSimulatorSkeleton: public ParticleSimulator
{
public:
	ParticleSimulatorSkeleton(float thres = 0.99) : ParticleSimulator()
	{
		this->thres = thres;
	}
	virtual bool Simulate(float endtime = 10.0f, float delta = 0.1f, bool bdebug = false);
	vector<ParticleSkeleton> skeletons;
	float thres;
};

