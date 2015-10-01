#include "Squid.h"

Cpoint*
ReadSquidPoints(char* filename, int& h, int& w, int& num) {
	FILE  *fp;
	/* try to open the file */
	if ((fp = fopen(filename, "rt")) == NULL) {
		fprintf(stderr, "%s: error loading data from %s.\n", 
			"ReadSquidData()", filename);
		fclose(fp);
		abort();
	}
	fscanf(fp,"# %d", &num);
	Cpoint* pnt=new Cpoint[num];
	int i,y,x;
	h=0; w=0;
	for(i=0; i<num; ++i) {
		if(EOF==fscanf(fp, "%d %d", &y, &x))
			break;
		Cpoint pp((real)x,(real)y,0,0,0,0);
		pnt[i]=pp;
		if(y>h) h=y;
		if(x>w) w=x;
	}
	num=i;
	return pnt;
}

/*
Since the Petrakis point sets are coarsely sampled, I need to interpolate them
to get a larger number of points.
Here I used Catmull-Rom to interpolate and increase the number of points.
The parameter 'scale' provides the increase factor of the number of contour points.
*/
Cpoint*
ReadPetrakisPoints(char* filename, int& h, int& w, int& num, int scale) {
	FILE  *fp;
	/* try to open the file */
	if ((fp = fopen(filename, "rt")) == NULL) {
		fprintf(stderr, "%s: error loading data from %s.\n", 
			"ReadSquidData()", filename);
		fclose(fp);
		abort();
	}
	//Find the number of points
	for(num=0; ;++num) {
		real ty,tx;
		if(EOF==fscanf(fp, "%lf %lf", &ty, &tx))
			break;
	}
	fclose(fp);
	
	/* reopen the file */
	if ((fp = fopen(filename, "rt")) == NULL) {
		fprintf(stderr, "%s: error loading data from %s.\n", 
			"ReadSquidData()", filename);
		fclose(fp);
		abort();
	}
	int i,k;
	h=0; w=0;
	real* y=new real[num];
	real* x=new real[num];
	Cpoint* pnt=new Cpoint[num*scale];
	//Read the data
	for(i=0; i<num; ++i) {
		fscanf(fp, "%lf %lf", y+i, x+i);
		if((int)y[i]>h) h=(int)y[i];
		if((int)x[i]>w) w=(int)x[i];
	}
	//Interpolation
	real inc=1.0/scale;
	for(i=0,k=0; i<num; ++i,k+=scale) {
		//Catmull-Rom
		real x0=x[(i-1+num)%num];
		real y0=y[(i-1+num)%num];
		real x1=x[i];
		real y1=y[i];
		real x2=x[(i+1)%num];
		real y2=y[(i+1)%num];
		real x3=x[(i+2)%num];
		real y3=y[(i+2)%num];
		int n;
		real t;
		for(n=0, t=0; n<scale; ++n, t+=inc) {
			real t2=t*t;
			real t3=t2*t;
			real xx=(-t+2.*t2-t3)*x0+(2.-5.*t2+3.*t3)*x1+(t+4.*t2-3.*t3)*x2+(-t2+t3)*x3;
			real yy=(-t+2.*t2-t3)*y0+(2.-5.*t2+3.*t3)*y1+(t+4.*t2-3.*t3)*y2+(-t2+t3)*y3;
			pnt[k+n]=Cpoint(xx,yy,0,0,0,0);
		}
	}
	delete [] y;
	delete [] x;
	num*=scale;
	return pnt;
}

Cpoint*
ReadBooksteinPoints(char* filename, int& h, int& w, real scale, int& num) {
	FILE  *fp;
	/* try to open the file */
	if ((fp = fopen(filename, "rt")) == NULL) {
		fprintf(stderr, "%s: error loading data from %s.\n", 
			"ReadSquidData()", filename);
		fclose(fp);
		abort();
	}
	fscanf(fp,"# %d", &num);
	Cpoint* pnt=new Cpoint[num];
	int i;
	real* y=new real[num];
	real* x=new real[num];
	h=0; w=0;
	real miny, maxy, minx, maxx;
	for(i=0; i<num; ++i) {
		if(EOF==fscanf(fp, "%lf %lf", y+i, x+i))
			break;
		if(i==0) {
			minx=maxx=x[0];
			miny=maxy=y[0];
		}
		else {
			if(x[i]<minx) minx=x[i];
			else if(x[i]>maxx) maxx=x[i];
			if(y[i]<miny) miny=y[i];
			else if(y[i]>maxy) maxy=y[i];
		}
	}
	num=i; // there may be less points than specified

	h=(int)(scale*(maxy-miny))+4;
	w=(int)(scale*(maxx-minx))+4;
	for(i=0; i<num; ++i) {
		real yy=scale*(y[i]-miny)+2;
		real xx=scale*(x[i]-minx)+2;
		Cpoint pp(xx,yy,0,0,0,0);
		pnt[i]=pp;
		cout << pp;
	}

	delete [] y;
	delete [] x;
	return pnt;
}
