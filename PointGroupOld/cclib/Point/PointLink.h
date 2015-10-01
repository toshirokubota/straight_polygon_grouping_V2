#ifndef _PointLink_h_
#define _PointLink_h_

#include <iostream>
#include <fstream>
using namespace std;
#include <vector>

#include <stdlib.h>
#include <assert.h>

#include <mytype.h>
#include <Misc.h>
#include "Cpoint.h"
#include <Draw.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	vector<Cpoint> link;
	vector<int> num;
} struct_link; 

#ifdef __cplusplus
}
#endif

void
CopyPoints(Cpoint* dest, Cpoint* src, int nump);

void
ShowPoints(RealImage& map, Cpoint* points, int nump, bool closed);

bool
StartTrack(const RealImage& edge, const ByteImage& flag, int y, int x);

bool
FindNextPoint(const RealImage& edge, const ByteImage& flag, 
							int y, int x, int& y2, int& x2);

bool
FindNextPointNoCanny(const RealImage& edge, const ByteImage& flag, 
							int y, int x, int& y2, int& x2);

Cpoint*
Combine(Cpoint* p1, Cpoint* p2, int c1, int c2);

Cpoint*
Track(const RealImage& edge, ByteImage& flag, 
			int y, int x, int& cnt);

Cpoint*
TrackNoCanny(const RealImage& edge, ByteImage& flag, 
			int y, int x, int& cnt);



void
WritePointLinks(struct_link& seg, char* filename);

struct_link
ReadPointLinks(char* filename);

#endif /* _PointLink_h */
