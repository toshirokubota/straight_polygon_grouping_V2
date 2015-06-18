#include <Polygon.h>
#include <MovingParticle.h>

vector<CParticleF> 
Polygon::project(float time)
{
	vector<CParticleF> c;
	for (int i = 0; i < particles.size(); ++i)
	{
		c.push_back(particles[i]->project(time));
	}
	return c;
}

Polygon::Polygon(vector<MovingParticle*>& vp, float creation_time)
{
	this->particles = vp;
	this->creation_time = creation_time;
	id = _id++;
	level = 0;
	for (int i = 0; i < vp.size(); ++i)
	{
		//pset.insert(vp[i]);
	}
	for (int i = 0; i < vp.size(); ++i)
	{
		sset.insert(vp[i]->getInitParticle());
	}
}