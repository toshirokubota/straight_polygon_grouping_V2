#pragma once

#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <mex.h>
#include <szMexUtility.h>
#include <StationaryParticle.h>
#include <IntersectionConvexPolygons.h>
#include <MiscGeometry.h>

class LinkedTriple
{
public:
	LinkedTriple(StationaryParticle* p0, StationaryParticle* q0, StationaryParticle* r0, int id0);
	static float fitnessMeasure(StationaryParticle* p, StationaryParticle* q, StationaryParticle* r);
	static bool isNil(LinkedTriple* t);
	static bool calculateVelocity(StationaryParticle* p, StationaryParticle* q, StationaryParticle* r, float& vx, float& vy);
	static float compatibility(CParticleF& p, float u[], CParticleF& q, float v[]); //include the angle factor
	static float _timeToClosestEncounter(CParticleF& p, float u[], CParticleF& q, float v[]);

	void print();
	LinkedTriple* best();
	//float compatibility0(LinkedTriple* t); //only the distance factor
	float compatibility(LinkedTriple* t); //include the angle factor
	float _timeToClosestEncounter(LinkedTriple* t);
	bool updateFate();

	StationaryParticle* p; //center
	StationaryParticle* q; //right
	StationaryParticle* r; //left
	float v[2];
	float fitness;
	CParticleF fate;
	CParticleF _fate0; //temporary storage of an updated fate
	vector<LinkedTriple*> frontSupporters;
	vector<LinkedTriple*> leftSupporters;
	vector<LinkedTriple*> rightSupporters;
	vector<LinkedTriple*> competitors;
	//vector<float> compatibility; //fitenss of each supporter
	vector<float> linkWeights; //link weight to each supporter
	vector<float> linkWeights0; 
	float prob;
	float _prob0; //temporary storage - this will be the prob in the next period
	float _totalFitness;

	int id;
};

struct LinkedTripleFactory
{
public:
	static LinkedTripleFactory& getInstance()
	{
		static LinkedTripleFactory instance;
		return instance;
	}
	LinkedTriple* makeTriple(StationaryParticle* sp, StationaryParticle* sr, StationaryParticle*  sq)
	{
		LinkedTriple* triple = new LinkedTriple(sp, sq, sr, _id++);
		if (_id == 605)
			_id += 0;
		triples.push_back(triple);
		return triple;
	}
	LinkedTriple* makeNillTriple(StationaryParticle* sp)
	{
		const float EPSILON = 1.0e-5;
		LinkedTriple* t = new LinkedTriple(sp, sp, sp, _id++);
		triples.push_back(t);
		//t->compatibility.push_back(EPSILON);
		t->linkWeights.push_back(1.0);
		t->v[0] = t->v[1] = 0.0f;
		t->fitness = EPSILON;

		return t;
	}
	void clean()
	{
		for (int i = 0; i<triples.size(); ++i)
		{
			delete triples[i];
		}
		triples.clear();
		_id = 0;
	}
	vector<LinkedTriple*> triples;
	float unit;
private:
	int _id;
	LinkedTripleFactory()
	{
		_id = 0;
		unit = 1.0f;
	}
	~LinkedTripleFactory()
	{
		clean();
	}
	LinkedTripleFactory(LinkedTripleFactory& f){}
	LinkedTripleFactory operator=(LinkedTripleFactory& f){}
};





