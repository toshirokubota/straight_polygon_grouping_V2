#pragma once

struct ParticleDirection {
	float x;
	float y;
	float z;
	ParticleDirection(float x = 0, float y = 0, float z = 0,  bool norm=false) 
	{
		this->x = x;
		this->y = y;
		this->z = z;
		if (norm)
		{
			float d = sqrt(x*x + y*y + z*z);
			x /= d;
			y /= d;
			z /= d;
		}
	}
	ParticleDirection(CParticleF& from, CParticleF& to, bool norm = false)
	{
		x = to.m_X - from.m_X;
		y = to.m_Y - from.m_Y;
		z = to.m_Z - from.m_Z;
		if (norm)
		{
			float d = sqrt(x*x + y*y + z*z);
			x /= d;
			y /= d;
			z /= d;
		}
	}
};