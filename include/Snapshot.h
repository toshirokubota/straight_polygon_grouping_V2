#pragma once
#include <vector>
using namespace std;
#include <szParticleF.h>
#include <MovingParticle.h>
#include <mex.h>

class Snapshot
{
public:
	Snapshot(float t=0.0f, float s=0.0f)
	{
		static int _id = 0;
		created_time = t;
		projection_time = s;
		id = _id++;
	}
	Snapshot(float t, float s, vector<MovingParticle*>& particles) : Snapshot(t, s)
	{
		for (int i = 0; i < particles.size(); ++i)
		{
			add(particles[i]->getId());
		}
	}
	static vector<Snapshot> TakeSnapshot(float time);
	static Snapshot LoadSnapshot(const mxArray* ptr);
	static vector<Snapshot> LoadSnapshots(const mxArray* ptr);
	static mxArray* StoreSnapshot(Snapshot& snapshot);
	static mxArray* StoreSnapshots(vector<Snapshot>& snapshot);
	int size() const { return particles.size(); }
	float getProjectionTime() const { return projection_time; }
	float getCreationTime() const { return created_time; }
	int get(int i) const { return particles[i]; }
	bool operator==(const Snapshot& shot) const
	{
		return projection_time == shot.projection_time && shot.size() == this->size() && shot._pset == this->_pset;
	}
	int getId() const { return id; }
private:
	void add(int id)
	{
		particles.push_back(id);
		_pset.insert(id);
	}
	float projection_time;
	float created_time;
	vector<int> particles;
	set<int> _pset;
	int id;
};
