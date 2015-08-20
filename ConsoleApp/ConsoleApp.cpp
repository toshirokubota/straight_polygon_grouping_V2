// ConsoleApp.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <set>
#include <algorithm>
#include <iterator>
#include <iostream>
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	set<int> s1;
	for (int i = 0; i < 100; i += 5)
	{
		s1.insert(i);
	}
	set<int> s2;
	for (int i = 100; i > 0; i-= 7)
	{
		s2.insert(i);
	}
	for (set<int>::iterator it = s1.begin(); it != s1.end(); ++it)
	{
		int k = *it;
		if (s2.find(k) != s2.end())
		{
			cout << k << endl;
		}
	}
	set<int> iset;
	set<int> uset;
	set_intersection(s1.begin(), s1.end(),
		s2.begin(), s2.end(),
		std::inserter(iset, iset.begin()));
	set_union(s1.begin(), s1.end(),
		s2.begin(), s2.end(),
		std::inserter(uset, uset.begin()));
	float coverage = 1 - (float)iset.size() / (float)uset.size();
	cout << iset.size() << " " << uset.size() << " " << coverage << endl;
	return 0;
}

