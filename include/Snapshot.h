#pragma once
#include <vector>
using namespace std;
#include <szParticleF.h>
#include <MovingParticle.h>
#include <mex.h>
#include <Polygon.h>

class Snapshot
{
public:
	Snapshot(float t, float s, vector<MovingParticle*>& particles) 
	{
		PolygonFactory& factory = PolygonFactory::getInstance();
		polygon = factory.makePolygon(particles, t);
		created_time = t;
		projection_time = s;
	}
	static vector<Snapshot> TakeSnapshot(float time);
	//static Snapshot LoadSnapshot(const mxArray* ptr);
	//static vector<Snapshot> LoadSnapshots(const mxArray* ptr);
	static mxArray* StoreSnapshot(Snapshot& snapshot);
	static mxArray* StoreSnapshots(vector<Snapshot>& snapshot);
	float getProjectionTime() const { return projection_time; }
	float getCreationTime() const { return created_time; }
	bool operator==(const Snapshot& shot) const
	{
		return projection_time == shot.projection_time && *polygon == *shot.polygon;
	}
	int getId() const { return id; }
private:
	float projection_time;
	float created_time;
	Polygon* polygon;
	int id;
};
