#ifndef _Draw_h_
#define _Draw_h_

#include <Image.h>
#include <mytype.h>
#include <Misc.h>

const real LineValue=1.0;
const real CornerValue=1.0;

void
DrawLine(RealImage& image, int b, int y1, int x1, int y2, int x2, real val);

void
DrawCross(RealImage& image, int b, int y, int x, int length, real val);

void
DrawLine(RealImage& image, int b, int y1, int x1, int y2, int x2);

void
DrawCross(RealImage& image, int b, int y, int x, int length);

void
DrawLine(ByteImage& image, int b, int y1, int x1, int y2, int x2, unsigned char val);

void
DrawCross(ByteImage& image, int b, int y, int x, int length, unsigned char val);

#endif /* Draw_h */