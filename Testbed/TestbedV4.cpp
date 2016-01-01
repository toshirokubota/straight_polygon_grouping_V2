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
#include <StationaryParticle.h>
#include <DisjointSet.h>

struct MAPoint
{
	MAPoint(StationaryParticle* p, StationaryParticle* q)
	{
		this->x = (p->getX() + q->getX()) / 2.0;
		this->y = (p->getY() + q->getY()) / 2.0;
		this->z = Distance(p->getP(), q->getP()) / 2.0;
		this->theta = GetVisualDirection(q->getX(), q->getY(), p->getX(), p->getY());
		this->p = p;
		this->q = q;
	}
	void print(char* tab="", char* end="\n")
	{
		printf("%s%d(%3.3f,%3.3f) - %d(%3.3f,%3.3f) => %3.3f, %3.3f, %3.3f, %3.3f%s",
			tab, p->getId(), p->getX(), p->getY(), q->getId(), q->getX(), q->getY(), x, y, z, theta, end);
	}
	float x;
	float y;
	float z;
	float theta;
	StationaryParticle* p;
	StationaryParticle* q;
};


float distance(MAPoint* p, MAPoint* q, float  wz = 9.0f, float wt = 100.0f)
{
	float dx = p->x - q->x;
	float dy = p->y - q->y;
	float dz = p->z - q->z;
	float dt = sin(p->theta - q->theta);
	return sqrt(dx*dx + dy*dy + wz*dz*dz + wt*dt*dt);
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

vector<MAPoint*>
collectMedialAxisPoints(vector<StationaryParticle*>& points)
{
	vector<MAPoint*> mapoints;
	for (int i = 0; i < points.size(); ++i)
	{
		for (int j = i + 1; j < points.size(); ++j)
		{
			mapoints.push_back(new MAPoint(points[i], points[j]));
		}
	}
	return mapoints;
}

vector<int>
clusterMAPoints(vector<MAPoint*>& P, float thres)
{
	vector<Node<MAPoint*>*> nodes;
	for (int i = 0; i < P.size(); ++i)
	{
		nodes.push_back(makeset(P[i]));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		//printf("%d: ", i + 1);
		//nodes[i]->key->print();
		MAPoint* p = nodes[i]->key;
		if (p->p->getId() == 4 && p->q->getId() == 12)
		{
			nodes[i]->key->print("\t");
		}
		for (int j = 0; j < nodes.size(); ++j)
		{
			float d = distance(nodes[i]->key, nodes[j]->key);
			if (p->p->getId()==4 && p->q->getId()==12)
			{
				nodes[j]->key->print("", " ");
				printf("%f\n", d);
			}
			if(d < thres)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<Node<MAPoint*>*> reps = clusters(nodes);
	map<Node<MAPoint*>*, int> imap;
	for (int i = 0; i < reps.size(); ++i)
	{
		imap[reps[i]] = i;
	}
	vector<int> labels;
	for (int i = 0; i < nodes.size(); ++i)
	{
		labels.push_back(imap[findset(nodes[i])]);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return labels;
}

vector<vector<int> >
groupPoints(vector<MAPoint*>& points, vector<int> labels)
{
	set<int> iset;
	map<int, int> imap;
	int count = 0;
	for (int i = 0; i < labels.size(); ++i)
	{
		if (iset.find(labels[i]) == iset.end())
		{
			imap[labels[i]] = count;
			count++;
			iset.insert(labels[i]);
		}
	}
	vector<set<int>> G(iset.size());
	for (int i = 0; i < points.size(); ++i)
	{
		int lb = labels[i];
		int k = imap[lb];
		G[k].insert(points[i]->p->getId());
		G[k].insert(points[i]->q->getId());
	}
	vector<vector<int>> vG(iset.size());
	for (int i = 0; i < G.size(); ++i)
	{
		for (set<int>::iterator it = G[i].begin(); it != G[i].end(); it++)
		{
			vG[i].push_back(*it);
		}
	}
	return vG;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = Testbed(P, [sigma, step])");
		return;
	}
	//Points
	vector<CParticleF> P;
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
			P.push_back(CParticleF(x, y));
		}
	}

	float thres = 5.0f;
	if (nrhs >= 2)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		ReadScalar(thres, prhs[1], classId);
	}

	vector<StationaryParticle*> points = collectStationaryParticles(P);
	vector<MAPoint*> mapoints = collectMedialAxisPoints(points);
	vector<int> labels = clusterMAPoints(mapoints, thres);
	vector<vector<int>> G = groupPoints(mapoints, labels);

	if (nlhs >= 1)
	{
		const int dims[] = { G.size(), 1 };
		vector<vector<int>> vvint;
		for (int i = 0; i < dims[0]; ++i)
		{
			const int dims2[] = { G[i].size(), 1 };
			vector<int> vint(dims2[0] * dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				SetData2(vint, j, 0, dims2[0], dims2[1], G[i][j]);
			}
			vvint.push_back(vint);
		}
		plhs[0] = StoreDataCell(vvint, mxINT32_CLASS, 2, dims, 1);
	}
	/*if (nlhs >= 2)
	{
		plhs[1] = StoreData(X, mxSINGLE_CLASS, 3, dims);
	}*/

	for (int i = 0; i < mapoints.size(); ++i)
	{
		delete mapoints[i];
	}
	StationaryParticleFactory::getInstance().clean();
	mexUnlock();
}

