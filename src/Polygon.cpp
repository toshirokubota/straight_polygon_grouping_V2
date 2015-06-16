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
