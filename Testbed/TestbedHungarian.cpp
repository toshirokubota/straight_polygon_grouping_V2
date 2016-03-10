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
#include <Hungarian.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = Testbed(P, [sigma, step])");
		return;
	}
	//Points
	vector<float> C;
	const int* dimsC;
	mxClassID classC;
	int ndimC;
	LoadData(C, prhs[0], classC, ndimC, &dimsC);

	float INF = std::numeric_limits<float>::infinity();
	int n = dimsC[0];
	int m = dimsC[1];
	vector<vector<float>> a(n+1);
	a[0] = vector<float>(m + 1, 0);
	for (int i = 1; i <= n; ++i)
	{
		vector<float> a0(m+1, 0);
		for (int j = 1; j <= m; ++j)
		{
			a0[j] = -C[(i - 1)*m + j - 1]; 
		}
		a[i] = a0;
	}

	vector<float> u(n + 1);
	vector<float> v(m + 1);
	vector<int> p(m + 1);
	vector<int> way(m + 1);
	for (int i = 1; i <= n; ++i) {
		p[0] = i;
		int j0 = 0;
		vector<float> minv(m + 1, INF);
		vector<char> used(m + 1, false);
		do {
			used[j0] = true;
			int i0 = p[j0];
			float delta = INF;
			int j1;
			for (int j = 1; j <= m; ++j)
			if (!used[j]) {
				float cur = a[i0][j] - u[i0] - v[j];
				if (cur < minv[j])
					minv[j] = cur, way[j] = j0;
				if (minv[j] < delta)
					delta = minv[j], j1 = j;
			}
			for (int j = 0; j <= m; ++j)
			if (used[j])
				u[p[j]] += delta, v[j] -= delta;
			else
				minv[j] -= delta;
			j0 = j1;
		} while (p[j0] != 0);
		do {
			int j1 = way[j0];
			p[j0] = p[j1];
			j0 = j1;
		} while (j0);
	}

	vector<int> ans(n + 1);
	for (int j = 1; j <= m; ++j)
		ans[p[j]] = j;

	int cost = -v[0];

	for (int i = 1; i <= m; ++i)
	{
		printf("%d=>%d\n", i, ans[i]);
	}
	printf("cost = %d\n", cost);

	printf("\nUsing Hungarian struct.\n");
	Hungarian algo(C, dimsC);
	float cost2 = algo.solve(true);
	printf("Cost = %f\n", cost2);
	printf("Assignment: \n");
	for (int i = 0; i < algo.nJobs; ++i)
	{
		int j = algo.getWorker(i);
		int k = algo.getJob(i);
		printf("%d=>%d, %d=>%d %f\n", i, j, i, k, algo.getCost(i, j));
	}
	mexUnlock();
}

