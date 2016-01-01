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

struct index3
{
	int x;
	int y;
	int z;
};

struct boundingBox
{
	boundingBox(vector<CParticleF>& P, vector<float>& step)
	{
		x0 = y0 = z0 = std::numeric_limits<float>::infinity();
		x1 = y1 = z1 = -std::numeric_limits<float>::infinity();
		for (int i = 0; i < P.size(); ++i)
		{
			x0 = Min(x0, P[i].m_X);
			x1 = Max(x1, P[i].m_X);
			y0 = Min(y0, P[i].m_Y);
			y1 = Max(y1, P[i].m_Y);
			z0 = Min(z0, P[i].m_Z);
			z1 = Max(z1, P[i].m_Z);
		}	
		this->step = step;
	}
	const int* getDimension()
	{
		dims[0] = (x1 - x0) / step[0];
		dims[1] = (y1 - y0) / step[1];
		dims[2] = (z1 - z0) / step[2];
		return dims;
	}
	CParticleF map2coordinate(int x, int y, int z)
	{
		CParticleF a;
		a.m_X = (x0 + x * step[0]);
		a.m_Y = (y0 + y * step[1]);
		a.m_Z = (z0 + z * step[2]);
		return a;
	}
	index3 map2index(CParticleF a)
	{
		index3 idx;
		idx.x = (a.m_X - x0) / step[0];
		idx.y = (a.m_Y - y0) / step[1];
		idx.z = (a.m_Z - z0) / step[2];
		return idx;
	}

	float x0, y0, z0, x1, y1, z1;
	vector<float> step;
	int dims[3];
};



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "StraightMedialAxis", __DATE__, __TIME__);
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
	vector<float> sigma(3, 5.0);
	if (nrhs >= 2)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		LoadData(sigma, prhs[1], classId, ndim, &dims);
	}
	vector<float> step(3, 1.0f);
	if (nrhs >= 3)
	{
		mxClassID classId;
		int ndim;
		const int* dims;
		LoadData(step, prhs[2], classId, ndim, &dims);
	}

	vector<CParticleF> Q;
	for (int i = 0; i < P.size(); ++i)
	{
		for (int j = i + 1; j < P.size(); ++j)
		{
			CParticleF q((P[i].m_X + P[j].m_X) / 2, (P[i].m_Y + P[j].m_Y) / 2.0, Distance(P[i], P[j]) / 2.0f);
			Q.push_back(q);
		}
	}

	boundingBox box(Q, step);
	const int* dims = box.getDimension();
	vector<float> X(dims[0] * dims[1] * dims[2], 0.0f);
	for (int m = 0; m < Q.size(); ++m)
	{
		for (int i = 0; i < dims[2]; i++)
		{
			for (int j = 0; j < dims[1]; j++)
			{
				for (int k = 0; k < dims[0]; k++)
				{
					CParticleF a = box.map2coordinate(k, j, i);
					float dx = a.m_X - Q[m].m_X;
					float dy = a.m_Y - Q[m].m_Y;
					float dz = a.m_Z - Q[m].m_Z;
					float val = exp(-dx*dx / (2.*sigma[0] * sigma[0]) - dy*dy / (2.*sigma[1] * sigma[1]) - dz*dz / (2.*sigma[2] * sigma[2]));
					float q = GetData3(X, k, j, i, dims[0], dims[1], dims[2], 0.0f);
					SetData3(X, k, j, i, dims[0], dims[1], dims[2], q + val);
				}
			}
		}
	}
	if (nlhs >= 1)
	{
		const int dims[] = { P.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			index3 idx = box.map2index(P[i]);
			SetData2(F, i, 0, dims[0], dims[1], idx.x);
			SetData2(F, i, 1, dims[0], dims[1], idx.y);
			SetData2(F, i, 2, dims[0], dims[1], 0);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	if (nlhs >= 2)
	{
		plhs[1] = StoreData(X, mxSINGLE_CLASS, 3, dims);
	}
	mexUnlock();
}

