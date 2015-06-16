#pragma once;
#include <MovingParticle.h>
#include <mex.h>
#include <Snapshot.h>
#include <PolygonDAG.h>

class ParticleSimulator
{
public:
	ParticleSimulator()
	{
		time = 0.0f;
	}
	virtual bool Prepare(vector<StationaryParticle*>& points, vector<pair<int, int>>& edges, float delta0 = 0.01);
	bool LoadParticles(vector<float>& state, const int* dims);
	mxArray* SaveParticles();
	mxArray* SaveDoneEvents();
	mxArray* SaveConvexity();
	virtual bool Simulate(float endtime = 10.0f, float delta = 0.1f, bool bdebug = false);
	float getTime() const { return time; }

	static vector<MovingParticle*>	initializePolygon(vector<Edge<StationaryParticle*>*>& edges); //temporary

	vector<Snapshot> snapshots; //id-time pair.
	vector<EventStruct> doneEvents;
	vector<Snapshot> closedRegions;
	vector<Snapshot> polygons; //offset polygon
	vector<Snapshot> traces; //the trace of each offset polygon in polygons

	PolygonDAG dag;
protected:
	bool initializePolygon(vector<MovingParticle*>& particles);
	float time;
};

