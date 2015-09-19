#pragma once
#include <ParticleSimulatorGreedy.h>

class ParticleSimulatorGreedyPartition : public ParticleSimulatorGreedy
{
public:
	ParticleSimulatorGreedyPartition() : ParticleSimulatorGreedy()
	{
	}
	virtual bool Simulate(int maxLevel, float thres, int minLength);
	virtual bool _PartitionRecusive(MovingParticle* p, int level, float thres, int minLength, int maxLevel);
	FitnessStruct findNextEventGreedy(vector<MovingParticle*>& particles);
	virtual FitnessStruct computeFitness(EventStruct ev, bool left);
};
