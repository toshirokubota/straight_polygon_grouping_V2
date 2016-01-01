#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>
#include <mex.h>
#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <map>
#include <set>
#include <hash_map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <LinkedTriple.h>
#include <StationaryParticle.h>
#include <IntersectionConvexPolygons.h>
#include <Triangulation.h>

struct ParticlePair {
	ParticlePair(StationaryParticle* p, StationaryParticle* q)
	{
		this->p = p;
		this->q = q;
		float d = Distance(p->getP(), q->getP()) / 2.0f;
		m = CParticleF((p->getX() + q->getX()) / 2.0, (p->getY() + q->getY()) / 2.0, d);
		id = _id++;
	}
	void print()
	{
		printf("%d - %d (%f, %f, %f)\n", p->getId(), q->getId(), m.m_X, m.m_Y, m.m_Z);
	}
	StationaryParticle* p;
	StationaryParticle* q;
	CParticleF m;
	int id;
	static int _id;
};
int ParticlePair::_id = 0;

vector<StationaryParticle*>
collectStationaryParticles(vector<CParticleF>& points)
{
	StationaryParticleFactory& factory = StationaryParticleFactory::getInstance();
	vector<StationaryParticle*> particles;
	for (int i = 0; i < points.size(); ++i)
	{
		StationaryParticle* sp = factory.makeParticle(points[i]);
		particles.push_back(sp);
	}
	return particles;
}

vector<ParticlePair*> 
collectParticlePairs(vector<StationaryParticle*>& particles)
{
	vector<ParticlePair*> pairs;
	for (int i = 0; i < particles.size(); ++i)
	{
		for (int j = 0; j < particles.size(); ++j)
		{
			pairs.push_back(new ParticlePair(particles[i], particles[j]));
		}
	}
	return pairs;
}

vector<vector<ParticlePair*>> clusterPairs(vector<ParticlePair*>& pairs, float thres, int minsize)
{
	vector<vector<ParticlePair*>> groups;
	vector<CParticleF> centroids;
	for (int i = 0; i < pairs.size(); ++i)
	{
		vector<ParticlePair*> p(1, pairs[i]);
		groups.push_back(p);
		centroids.push_back(pairs[i]->m);
	}
	while (true)
	{
		float mind = thres;
		pair<int, int> indices;
		for (int i = 0; i < groups.size(); ++i)
		{
			for (int j = i + 1; j < groups.size(); ++j)
			{
				float d = Distance(centroids[i], centroids[j]);
				if (d < mind)
				{
					mind = d;
					indices.first = i;
					indices.second = j;
				}
			}
		}
		if (mind >= thres)
		{
			break;
		}
		else
		{
			int k1 = indices.first, k2 = indices.second;
			int n1 = groups[k1].size(), n2 = groups[k2].size();
			vector<ParticlePair*> g2 = groups[k1];
			g2.insert(g2.begin(), groups[k2].begin(), groups[k2].end());
			CParticleF m((n1 * centroids[k1].m_X + n2 * centroids[k2].m_X) / (n1 + n2),
				(n1 * centroids[k1].m_Y + n2 * centroids[k2].m_Y) / (n1 + n2), (n1 * centroids[k1].m_Z + n2 * centroids[k2].m_Z) / (n1 + n2));
			groups.erase(groups.begin() + k2);
			groups.erase(groups.begin() + k1);
			groups.push_back(g2);
			centroids.erase(centroids.begin() + k2);
			centroids.erase(centroids.begin() + k1);
			centroids.push_back(m);
		}
	}
	

	for (int i = groups.size() - 1; i >= 0; i--)
	{
		if (groups[i].size() < minsize)
		{
			groups.erase(groups.begin() + i);
		}
	}

	return groups;
}

vector<StationaryParticle*>
retrieveParticle(vector<ParticlePair*>& pairs) {
	set<StationaryParticle*> pset;
	for (int i = 0; i < pairs.size(); ++i)
	{
		pset.insert(pairs[i]->p);
		pset.insert(pairs[i]->q);
	}
	vector<StationaryParticle*> P;
	for (set<StationaryParticle*>::iterator it = pset.begin(); it != pset.end(); ++it)
	{
		P.push_back(*it);
	}
	return P;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "CommonFateGrouping", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = CommonFateGrouping(P, thres)");
		return;
	}
	//Points
	vector<CParticleF> points;
	const int* dimsP;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for (int i = 0; i < dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			points.push_back(CParticleF(x, y));
		}

	}
	float thres = 1;
	int minsize = 10;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minsize, prhs[2], classMode);
	}

	vector<StationaryParticle*> particles = collectStationaryParticles(points);
	vector<ParticlePair*> pairs = collectParticlePairs(particles);
	vector<vector<ParticlePair*>> groups = clusterPairs(pairs, thres, minsize);
	vector<vector<StationaryParticle*>> groups2;
	for (int i = 0; i < groups.size(); ++i)
	{
		printf("Group %d:\n", i);
		for (int j = 0; j < groups[i].size(); ++j)
		{
			groups[i][j]->print();
		}
		groups2.push_back(retrieveParticle(groups[i]));
	}

	if (nlhs >= 1)
	{
		StationaryParticleFactory& spfactory = StationaryParticleFactory::getInstance();
		const int dims[] = { spfactory.particles.size(), 3 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], spfactory.particles[i]->getX());
			SetData2(F, i, 1, dims[0], dims[1], spfactory.particles[i]->getY());
			SetData2(F, i, 2, dims[0], dims[1], (float)spfactory.particles[i]->getId());
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	if (nlhs >= 2)
	{
		const int dims[] = { groups2.size(), 1 };
		vector<vector<int>> vvint;
		for (int i = 0; i < groups2.size(); ++i)
		{
			const int dims2[] = { groups2[i].size(), 1 };
			vector<int> vint(dims2[0] * dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				SetData2(vint, j, 0, dims2[0], dims2[1], groups2[i][j]->getId());
			}
			vvint.push_back(vint);
		}
		plhs[1] = StoreDataCell(vvint, mxINT32_CLASS, 2, dims, 1);

	}
	if (nlhs >= 3)
	{
		const int dims[] = { groups.size(), 1 };
		vector<vector<int>> vvint;
		for (int i = 0; i < groups.size(); ++i)
		{
			const int dims2[] = { groups[i].size(), 2 };
			vector<int> vint(dims2[0] * dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				SetData2(vint, j, 0, dims2[0], dims2[1], groups[i][j]->p->getId());
				SetData2(vint, j, 1, dims2[0], dims2[1], groups[i][j]->q->getId());
			}
			vvint.push_back(vint);
		}
		plhs[2] = StoreDataCell(vvint, mxINT32_CLASS, 2, dims, 2);
	}

	for (int i = 0; i < pairs.size(); ++i)
	{
		delete pairs[i];
	}
	StationaryParticleFactory::getInstance().clean();
	mexUnlock();
}

