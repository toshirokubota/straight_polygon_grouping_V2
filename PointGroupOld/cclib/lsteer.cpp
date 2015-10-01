/*
This file contains various routines used for Fourier Series Approximation.
*/
#include <kim.h>
#include "../magstr/dipole.h"

double (*Func_to_steer)(double);
double Steer_x, Steer_y;
int Steer_order;
int Y;

void
Bind_Func_to_steer(double (*func)(double))
{
  Func_to_steer = func;
}

double
func_cos(double x)
{
  return (Func_to_steer(x) * cos(x * Steer_order));
}

double
func_sin(double x)
{
   return (Func_to_steer(x) * sin(x * Steer_order));
}


/*
width, height: the size of the filters
order: the order of the approximation
maxorder: the number of filters generated and a subset of the filters
          are selected.
cmpl: if = 1, the complex pairs are always selected together and the
      resulted filter has the complex data type.
      if = 0, the whatever filters with big energy are selected, and
      the resulted filter has the floating point type.
*/
kimage *
steer_base_filters(int width, int height, int order, int maxorder, int cmpl)
{
   extern double Steer_x, Steer_y;
   extern int Steer_order;
   kimage *t1, *t2;
   double incx, incy, pi, itwopi, ipi, d;
   real **q1, **q2, *p1, *p2, max, *eng;
   int i, j, k, j2, mk;
   char arg[100];

   if(order > maxorder)
   {
      sprintf(arg, "order=%d, maxorder=%d", order, maxorder);
      nrwarn("steer_base_filters",
         "order has to be <= maxorder\norder set to maxorder", 
         arg);
      maxorder = order;
   }
   if(order <= 0 || maxorder <= 0)
   {
      sprintf(arg, "order=%d, maxorder=%d", order, maxorder);
      nrwarn("steer_base_filters",
         "order has to be <= maxorder\norder set to maxorder", 
         arg);
      maxorder = order;
   }

   eng= imcalloc( 2 * maxorder, real);
   t1 = alloc_kimage(width, height, maxorder * 2, KIMAGE_FLOAT);

   incx = 2.0 / (width - 1);
   incy = 2.0 / (height - 1);
   pi = (double) PI;
   itwopi = 1.0 / (2.0 * pi);
   ipi = 1.0 / pi;

   printf("Order 0\n");
   q1 = t1->data[0];
   Steer_order = 0;
   for(i = 0, Steer_y = -1.0; i < t1->col; ++i, Steer_y += incy)
   {
      p1 = q1[i];
      for(j = 0, Steer_x = -1.0; j < t1->row; ++j, Steer_x += incx)
      {
         p1[j] = qsimp(func_cos, -pi, pi) * itwopi;
         eng[0] += p1[j] * p1[j];
      }
   }

   Steer_order = 1;
   for(k = 2; k < 2 * maxorder; k += 2, Steer_order ++)
   {
      printf("Order %d\n", Steer_order);
      q1 = t1->data[k];
      q2 = t1->data[k + 1];
      for(i = 0, Steer_y = -1.0; i < t1->col; ++i, Steer_y += incy)
      {
         p1 = q1[i];
         p2 = q2[i];
         for(j = 0, Steer_x = -1.0; j < t1->row; ++j, Steer_x += incx)
         {
            p1[j] = qsimp(func_cos, -pi, pi) * ipi;
            p2[j] = qsimp(func_sin, -pi, pi) * ipi;
            eng[k] += p1[j] * p1[j];
            eng[k + 1] += p2[j] * p2[j];
         }
      }
   }

   for(k = 0; k < 2 * maxorder; ++k)
      printf("engery @ %d = %f\n", k, eng[k]);

   if(cmpl) {
     t2 = alloc_kimage(width, height, order, KIMAGE_COMPLEX);
     t2->param1 = imcalloc(t2->bands, char);
     t2->ps1 = t2->bands;
     for(k = 0; k < order; ++k){
       for(i = 0, mk = 0, max=-1.0; i < maxorder; ++i){
         d = eng[2*i] + eng[2*i+1];
         if(max < d){
           max = d;
           mk =i;
         }
       }
       printf("Select band#%d (energy=%f)\n", mk, eng[2*mk]+eng[2*mk+1]);
       eng[2*mk]= eng[2*mk+1]= -2.0;

       for(i = 0; i < t1->col; ++i) {
         for(j = 0, j2 = 0; j < t1->row; ++j, j2+=2) {
           t2->data[k][i][j2] = t1->data[2*mk][i][j];
           t2->data[k][i][j2+1] = t1->data[2*mk+1][i][j];
         }
       }
       t2->param1[k] = (unsigned char) mk; /* order */
     }
   }
   else {
     t2 = alloc_kimage(width, height, order * 2, KIMAGE_FLOAT);
     t2->param1 = imcalloc(t2->bands, char);
     t2->param2 = imcalloc(t2->bands, char);
     t2->ps1 = t2->bands;
     t2->ps2 = t2->bands;
     for(k = 0; k < 2 * order; ++k){
       for(i = 0, mk = 0, max=-1.0; i < 2 * maxorder; ++i){
         if(max < eng[i]){
           max = eng[i];
           mk =i;
         }
       }
       printf("Select band#%d (energy=%f)\n", mk, eng[mk]);
       eng[mk]= -2.0;

       replace_band_kimage(t2, t1, k, mk, 0);
       t2->param1[k] = (unsigned char) (mk >> 1); /* order */
       t2->param2[k] = (char) mk & 1;  /* sin or cosin */
     }
   }
   free(eng);
   free_kimage(&t1);
   return(t2);
}

kimage *
steer_kimage(kimage *fil, real theta)
{
   kimage *res;
   real **q1, **q2, *p1, *p2, wc1, wc2;
   int i, j, k, n, j2, cosin;

   res = alloc_kimage(fil->col, fil->row, 1, KIMAGE_FLOAT);
   
   q2 = res->data[0];
   for(k = 0; k < fil->bands; k ++)
   {
     q1 = fil->data[k];
     if(fil->type == KIMAGE_COMPLEX) {
       n = (int) fil->param1[k];
       if (!n) {
         wc1 = .5;
         wc2 = .0;
       }
       else {
         wc1 = cos(n * (double) theta);
         wc2 = sin(n * (double) theta);
       }
       for(i = 0; i < fil->col; ++i){
         p1 = q1[i];
         p2 = q2[i];
         for(j = 0, j2=0; j < fil->row; j++, j2+=2) {
           p2[j] += wc1 * p1[j2] + wc2 * p1[j2+1];
         }
       }
     }
     else {
       cosin = (int) !fil->param2[k];
       n = (int) fil->param1[k];
       if(cosin && n)
         wc1 = cos(n * (double) theta);
       else if (cosin && !n)
         wc1 = .5;
       else
         wc1 = sin(n * (double) theta);
       for(i = 0; i < fil->col; ++i){
         p1 = q1[i];
         p2 = q2[i];
         for(j = 0; j < fil->row; ++j){
           p2[j] += wc1 * p1[j];
         }
       }
     }
   }

   return(res);
}


/* the following is the functional declaration for dipole field generation.
   In order to make the routine genrate steerable filters,
   you have to bind Dipole_Func_to_steer to Func_to_steer by calling:
   Bind_Func_to_steer(Dipole_Func_to_steer);
 */

double
Dipole_Func_to_steer(double x)
{
  extern double Steer_x, Steer_y;
  double dmx, dmy, dpx, dpy, dm, dp, d, hx, hy;
  double ppx, ppy, npx, npy;
  int row;

  /* Converting into hexagonal lattice */
  row = ROUND((Steer_y + 1.0) * DIPOLE_FIELD_RANGE);
  hy = DIPOLE_FIELD_RANGE * HEIGHT_SPACE * Steer_y;
  hx = DIPOLE_FIELD_RANGE * Steer_x;
  hx = (row & 1) ? hx: hx - WIDTH_OFFSET;

  d = DIPOLE_FIELD_RANGE * DIPOLE_FIELD_RANGE * (Steer_x * Steer_x + Steer_y * Steer_y);
  if(d > .03125 * DIPOLE_LENGTH * DIPOLE_LENGTH) {
#ifndef DEBUG
    ppx = PLUS_POLE_X(.0, x, DIPOLE_LENGTH);
    ppy = PLUS_POLE_Y(.0, x, DIPOLE_LENGTH);
    npx = MINUS_POLE_X(.0, x, DIPOLE_LENGTH);
    npy = MINUS_POLE_Y(.0, x, DIPOLE_LENGTH);
    dmx = npx - hx;
    dmy = npy - hy;
    dpx = hx - ppx;
    dpy = hy - ppy;
    dm = sqrt(dmx * dmx + dmy * dmy);
    dp = sqrt(dpx * dpx + dpy * dpy);
    dm *= dm;
    dp *= dp;

    if(!Y) {
      d = dmx * pow(dm, (double) -.5*FIELD_ALPHA);
      d += dpx * pow(dp, (double) -.5*FIELD_ALPHA);
      //d = (dmx / dm + dpx / dp);
    }
    else {
      d = dmy * pow(dm, (double) -.5*FIELD_ALPHA);
      d += dpy * pow(dp, (double) -.5*FIELD_ALPHA);
      //d = (dmy / dm + dpy / dp);
    }
#else
    dm = atan2(hy, hx) - x;
    dp = (hy * hy + hx * hx);
    if(!Y) {
      d = cos(dm) / dp;
    }
    else {
      d = sin(dm) / dp;
    }
#endif
    }
    else {
      d = .0;
    }

  return(d);
}

#define ST_MAX_ORDER 10
//#define COMPLEX

kimage *
steer_filter_dipole(int order, int y)
{
   kimage *fil;
   int w, h;

   Y = y;

   w = 2*DIPOLE_FIELD_RANGE+1;
   h = 2*DIPOLE_FIELD_RANGE+1;
#ifdef COMPLEX
   fil = steer_base_filters(h, w, order, ST_MAX_ORDER, 1);
#else
   fil = steer_base_filters(h, w, order, ST_MAX_ORDER, 0);
#endif

   return(fil);
}

#undef ST_MAX_ORDER
#undef COMPLEX
