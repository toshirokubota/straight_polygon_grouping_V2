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
		m0 = m;
		id = _id++;
	}
	void print()
	{
		printf("%d - %d (%f, %f, %f)\n", p->getId(), q->getId(), m.m_X, m.m_Y, m.m_Z);
	}
	StationaryParticle* p;
	StationaryParticle* q;
	CParticleF m;
	CParticleF m0;
	int id;
	static int _id;
};
int ParticlePair::_id = 0;

void
perturbePoints(vector<CParticleF>& points, float scale)
{
	rndm(12345);
	for (int i = 0; i < points.size(); ++i)
	{
		points[i].m_X += scale * (rndm(0) - 0.5);
		points[i].m_Y += scale * (rndm(0) - 0.5);
	}
}

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

void
agglomeratePairs(vector<ParticlePair*>& pairs, float sigma, int niter)
{
	float rate = 0.5;
	for (int n = 0; n < niter; ++n)
	{
		vector<CParticleF> M(pairs.size());
		for (int i = 0; i < pairs.size(); ++i)
		{
			float sx = 0, sy = 0, sz = 0, sd = 0;
			CParticleF m = pairs[i]->m;
			for (int j = 0; j < pairs.size(); ++j)
			{
				float dx = pairs[j]->m.m_X - m.m_X;
				float dy = pairs[j]->m.m_Y - m.m_Y;
				float dz = pairs[j]->m.m_Z - m.m_Z;
				float dd = exp(-(dx*dx + dy*dy + dz*dz)/(2*sigma*sigma));
				if (dd > 0)
				{
					sx += dx * dd;
					sy += dy * dd;
					sz += dz * dd;
					sd += dd;
				}
			}
			CParticleF m2(sx/sd, sy/sd, sz/sd);
			M[i] = CParticleF(m.m_X + rate*m2.m_X, m.m_Y + rate*m2.m_Y, m.m_Z + rate*m2.m_Z);
		}
		for (int i = 0; i < pairs.size(); ++i)
		{
			pairs[i]->m = M[i];
		}
	}
}

vector<vector<ParticlePair*>> clusterPairs(vector<ParticlePair*>& pairs, float thres, int minsize)
{
	vector<Node<ParticlePair*>*> nodes;
	for (int i = 0; i < pairs.size(); ++i)
	{
		nodes.push_back(makeset(pairs[i]));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			if (Distance(nodes[i]->key->m, nodes[j]->key->m) < thres)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<Node<ParticlePair*>*> resp = clusters(nodes);
	vector<vector<ParticlePair*>> groups(resp.size());
	map<Node<ParticlePair*>*, int> nmap;

	for (int i = 0; i < resp.size(); ++i)
	{
		nmap[resp[i]] = i;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = nmap[findset(nodes[i])];
		groups[k].push_back(nodes[i]->key);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
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
	int numIter = 5;
	float noiseLevel = 0.0f;
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
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[3], classMode);
	}
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(noiseLevel, prhs[4], classMode);
	}

	perturbePoints(points, noiseLevel);

	vector<StationaryParticle*> particles = collectStationaryParticles(points);
	vector<ParticlePair*> pairs = collectParticlePairs(particles);
	agglomeratePairs(pairs, 2 * thres, numIter);
	vector<vector<ParticlePair*>> groups = clusterPairs(pairs, thres, minsize);
	vector<vector<StationaryParticle*>> G;
	vector<vector<ParticlePair*>> H;
	for (int i = 0; i < groups.size(); ++i)
	{
		vector<StationaryParticle*> ps = retrieveParticle(groups[i]);
		if (ps.size() >= minsize)
		{
			G.push_back(ps);
			H.push_back(groups[i]);
		}
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
		const int dims[] = { G.size(), 1 };
		vector<vector<int>> vvint;
		for (int i = 0; i < G.size(); ++i)
		{
			const int dims2[] = { G[i].size(), 1 };
			vector<int> vint(dims2[0] * dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				SetData2(vint, j, 0, dims2[0], dims2[1], G[i][j]->getId());
			}
			vvint.push_back(vint);
		}
		plhs[1] = StoreDataCell(vvint, mxINT32_CLASS, 2, dims, 1);

	}
	if (nlhs >= 3)
	{
		const int dims[] = { H.size(), 1 };
		vector<vector<int>> vvint;
		for (int i = 0; i < H.size(); ++i)
		{
			const int dims2[] = { H[i].size(), 2 };
			vector<int> vint(dims2[0] * dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				SetData2(vint, j, 0, dims2[0], dims2[1], H[i][j]->p->getId());
				SetData2(vint, j, 1, dims2[0], dims2[1], H[i][j]->q->getId());
			}
			vvint.push_back(vint);
		}
		plhs[2] = StoreDataCell(vvint, mxINT32_CLASS, 2, dims, 2);
	}
	if (nlhs >= 4)
	{
		const int dims[] = { pairs.size(), 8 };
		vector<float> F(dims[0]*dims[1]);
		for (int i = 0; i < pairs.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], (float)pairs[i]->p->getId());
			SetData2(F, i, 1, dims[0], dims[1], (float)pairs[i]->q->getId());
			SetData2(F, i, 2, dims[0], dims[1], (float)pairs[i]->m.m_X);
			SetData2(F, i, 3, dims[0], dims[1], (float)pairs[i]->m.m_Y);
			SetData2(F, i, 4, dims[0], dims[1], (float)pairs[i]->m.m_Z);
			SetData2(F, i, 5, dims[0], dims[1], (float)pairs[i]->m0.m_X);
			SetData2(F, i, 6, dims[0], dims[1], (float)pairs[i]->m0.m_Y);
			SetData2(F, i, 7, dims[0], dims[1], (float)pairs[i]->m0.m_Z);
		}
		plhs[3] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}

	for (int i = 0; i < pairs.size(); ++i)
	{
		delete pairs[i];
	}
	StationaryParticleFactory::getInstance().clean();
	mexUnlock();
}

