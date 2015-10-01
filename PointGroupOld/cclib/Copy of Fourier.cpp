
#include "mytype.h"
#include "Fourier.h"
#include "Complex.h"
#include <assert.h>

#define SWAP(a, b) tempr=(a); (a)=(b); (b) = tempr

/*
One dimensional Fourier Transform routine copied from Numerical Recipes.
*/

void four1(real* data, int nn, int isign)
{
	int n,mmax,m,j,istep,i;
	real wtemp,wr,wpr,wpi,wi,theta;
	real tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

/*
N dimensional Fourier Transform routine copied from Numerical Recipes.
*/

void
fourn(real *data, int *nn, int ndim, int isign)
{
	int i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
	int ibit, idim, k1, k2, n, nprev, nrem, ntot;
	real tempi, tempr;
	double theta, wi, wpi, wpr, wr, wtemp;
	
	ntot = 1;
	for (idim = 1; idim <= ndim; idim++)
		ntot *= nn[idim];
	
	nprev = 1;
	for (idim = ndim; idim >= 1; idim--)
	{
		n = nn[idim];
		nrem = ntot / (n * nprev);
		ip1 = nprev << 1;
		ip2 = ip1 * n;
		ip3 = ip2 * nrem;
		i2rev = 1;
		for (i2 = 1; i2 <= ip2; i2 += ip1)
		{
			if (i2 < i2rev)
			{
				for (i1 = i2; i1 <= i2 + ip1 - 2; i1 += 2)
				{
					for (i3 = i1; i3 <= ip3; i3 += ip2)
					{
						i3rev = i2rev + i3 - i2;
						SWAP(data[i3], data[i3rev]);
						SWAP(data[i3 + 1], data[i3rev + 1]);
					}
				}
			}
			ibit = ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit)
			{
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1 = ip1;
		while (ifp1 < ip2)
		{
			ifp2 = ifp1 << 1;
			theta = (double)(-isign * 6.28318530717959) / (ifp2 / ip1);
			wtemp = sin(.5 * theta);
			wpr = -2.0 * wtemp * wtemp;
			wpi = sin(theta);
			wr = 1.0;
			wi = 0.0;
			for (i3 = 1; i3 <= ifp1; i3 += ip1)
			{
				for (i1 = i3; i1 <= i3 + ip1 - 2; i1 += 2)
				{
					for (i2 = i1; i2 <= ip3; i2 += ifp2)
					{
						k1 = i2;
						k2 = k1 + ifp1;
						tempr = wr * data[k2] - wi * data[k2 + 1];
						tempi = wr * data[k2 + 1] + wi * data[k2];
						data[k2] = data[k1] - tempr;
						data[k2 + 1] = data[k1 + 1] - tempi;
						data[k1] += tempr;
						data[k1 + 1] += tempi;
					}
				}
				wr = (wtemp = wr) * wpr - wi * wpi + wr;
				wi = wi * wpr + wtemp * wpi + wi;
			}
			ifp1 = ifp2;
		}
		nprev *= n;
	}
}

/*
This routine performs 2-dimensional Fourier Transform.
if direction==0, it performs the foward transform.
if direction==1, it performs the inverse.
*/

ComplexImage
FourierTransform(const ComplexImage& image, int direction) {
	
	assert(direction==0 || direction==1);
	
	//first, resize the image so that its width and height are the
	//power of 2.  Zeros are padded in the expanded region.
	ComplexImage tmp_im = ResizeImagePowerTwo(image);
	
	ComplexImage result(tmp_im.NumBands(),tmp_im.NumRows(), tmp_im.NumCols());
	
	real* buffer;
	buffer = new real[2*result.NumRows()*result.NumCols()+1];
	assert(buffer);
	
	int i, j, k, index;
	int isign = (direction) ? -1 : 1;
	int nelm[3];
	nelm[1] = result.NumRows();
	nelm[2] = result.NumCols();
	
	Complex cval;
	real wgt = (direction) ? 1.0 / (result.NumRows()*result.NumCols()): 1.0;
	for(k=0; k<result.NumBands(); ++k) {
		for(i=0, index=1; i<result.NumRows(); ++i) {
			for(j=0; j<result.NumCols(); ++j, index+=2) {
				cval = tmp_im.GetPixel(k,i,j);
				buffer[index] = cval.GetReal();
				buffer[index+1] = cval.GetImag();
			}
		}
		fourn(buffer, nelm, 2, isign);
		//four1(buffer, nelm[2], isign);
		for(i=0, index=1; i<result.NumRows(); ++i) {
			for(j=0; j<result.NumCols(); ++j, index+=2) {
				cval.SetReal(buffer[index]);
				cval.SetImag(buffer[index+1]);
				result.SetPixel(k,i,j, cval*wgt);
			}
		}
	}
	delete [] buffer;
	
	return result;
}

ComplexImage
FourierTransform(const RealImage& image, int direction) {
	
	assert(direction==0 || direction==1);
	
	//first, resize the image so that its width and height are the
	//power of 2.  Zeros are padded in the expanded region.
	RealImage tmp_im = ResizeImagePowerTwo(image);
	
	ComplexImage result(tmp_im.NumBands(),tmp_im.NumRows(), tmp_im.NumCols());
	Complex cval;
	
	real* buffer;
	buffer = new real[2*result.NumRows()*result.NumCols()+1];
	assert(buffer);
	
	int i, j, k, index;
	int isign = (direction) ? -1 : 1;
	int nelm[3];
	nelm[1] = result.NumRows();
	nelm[2] = result.NumCols();
	
	real wgt = (direction) ? 1.0 / (result.NumRows()*result.NumCols()): 1.0;
	for(k=0; k<result.NumBands(); ++k) {
		for(i=0, index=1; i<result.NumRows(); ++i) {
			for(j=0; j<result.NumCols(); ++j, index+=2) {
				buffer[index] = tmp_im.GetPixel(k,i,j);
				buffer[index+1] = 0;
			}
		}
		fourn(buffer, nelm, 2, isign);
		for(i=0, index=1; i<result.NumRows(); ++i) {
			for(j=0; j<result.NumCols(); ++j, index+=2) {
				cval.SetReal(buffer[index]);
				cval.SetImag(buffer[index+1]);
				result.SetPixel(k,i,j, wgt*cval);
			}
		}
	}
	delete [] buffer;
	
	return result;
}
