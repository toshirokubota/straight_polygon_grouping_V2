/*

 Copyright 1992, Massachusetts Institute of Technology. All Rights Reserved.
 
  --------------------------------------------------------------------
  Permission to use, copy, or modify this software and its documentation
  for educational and research purposes only and without fee is hereby
  granted, provided that this copyright notice appear on all copies and
  supporting documentation.  For any other uses of this software, in
  original or modified form, including but not limited to distribution
  in whole or in part, specific prior permission must be obtained from
  M.I.T. and the authors.  These programs shall not be used, rewritten,
  or adapted as the basis of a commercial software or hardware product
  without first obtaining appropriate licenses from M.I.T.  M.I.T. makes
  no representations about the suitability of this software for any
  purpose.  It is provided "as is" without express or implied warranty.
  
   ---------------------------------------------------------------------
   
	
	 File: jacobi.c
	 Author: Matthew Turk
	 Date: 3/22/91
	 Modified:  Thad Starner 10/6/91
	 For further information see:
	 
	  Turk, M., and Pentland, A., (1991) Eigenfaces for Recognition, 
	  Journal of Cognitive Neuroscience, Vol. 3, No. 1, pp. 71-86.
	  
*/
#include <iostream.h>
#include <math.h>

/*
Use Jacobi Transformations of a symetric matrix to obtain the 
normalized eigenvectors and eigen values 
-- Taken from NUMERICAL RECIPES IN C (pp 364-366) --
*/

#define EPS 0.000000001

#define ROTATE(a, i, j, k, l, n) g=a[i*n+j]; h=a[k*n+l]; a[i*n+j] = g-s*(h+g*tau); \
a[k*n+l] = h+s*(g-h*tau);

void jacobi(float* a, int n, float* d, float* v, int* nrot) {
	
/* 
Computes all eigenvalues and eigenvectors of a real symmetric matrix A
On output, elements of A above the diagonal are destroyed,
D returns the eigenvalues of A, and V is a matrix whose columns are the
eigenvectors of A.
	*/
	
	
	int j,iq,ip,ip_times_n,i ;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z,*vector();
	
	b = new float[n];
	z = new float[n];
	
	
	for(ip_times_n=0, ip=0; ip<n; ++ip, ip_times_n+=n)  
	{
		
		/* Initialize the identity matrix */
		for(iq=0; iq<n; ++iq)v[ip_times_n + iq] = 0.0 ;
		v[ip_times_n + ip] = 1.0 ;
		
		/* Initailize b and d to be diagonal of a */
		b[ip] = d[ip] = a[ip_times_n + ip];
		z[ip] = 0.0 ;
	}
	
	*nrot = 0 ;
	for(i=0;i<50;++i)
	{
		/* Sum off-diagonal elements */
		sm=0.0 ;
		
		for(ip_times_n=0,ip=0;ip<n-1;ip++,ip_times_n+=n)
			for(iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip_times_n + iq]);
			
			/*  If we have converged,  free the working vectors and return. */
			if(sm == 0.0)
			{
				delete [] b;
				delete [] z;
				return;
			}
			
			/* tresh is different on first three iterations...*/
			tresh=(i<3) ? 0.2*sm/(n*n) : 0.0 ;
			
			for(ip_times_n=0,ip=0;ip<n-1;ip++,ip_times_n+=n)
			{
				for(iq=ip+1;iq<n;++iq)
				{
					g=100.0*fabs(a[ip_times_n + iq]);
					
					/* After four sweeps, skip the rotation if the off-diagonal element is small */
					/* This test is taken directly from the text and looks a little suspect to me... */
					
					if(i > 3 && g < EPS)
						a[ip_times_n + iq] = 0.0 ;
					
					else if(fabs(a[ip_times_n+iq]) > tresh) 
					{
						h=d[iq]-d[ip];
						if(g < EPS)
							t = (fabs(a[ip_times_n+iq]) > EPS) ? (a[ip_times_n+iq])/h : 0.0 ; 
						else
						{ 
							theta=(fabs(h) < EPS) ? 0.0 : 0.5*h/(a[ip_times_n+iq]);
							t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
							if(theta < 0.0)
								t = -t ; 
						} 
						c=1.0/sqrt(1.0+t*t);
						s=t*c;
						tau=s/(1.0+c);
						
						h=t*a[ip_times_n+iq];
						z[ip] -= h;
						z[iq] += h;
						d[ip] -= h;
						d[iq] += h;
						a[ip_times_n+iq]=0.0;
						
						for(j=0;j<ip;j++)
						{
							ROTATE(a,j,ip,j,iq,n);
						}
						for(j=ip+1;j<iq;j++)
						{
							ROTATE(a,ip,j,j,iq,n);
						}
						for(j=iq+1;j<n;j++)
						{
							ROTATE(a,ip,j,iq,j,n);
						}
						for(j=0;j<n;j++)
						{
							ROTATE(v,j,ip,j,iq,n);
						}
						++(*nrot);
					}
				}
			}
			for(ip=0;ip<n;++ip)
			{
				b[ip] += z[ip];
				d[ip]=b[ip];
				z[ip]=0.0;
			}
	}
	
	/* Failed to converge!! What to do ??? */
	/* Well, let's at least free up memory and return without a murmur */
	
	delete [] b;
	delete [] z;
	return;	
}
