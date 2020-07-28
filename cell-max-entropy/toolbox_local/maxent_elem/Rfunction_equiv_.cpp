// Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
//
// Rfunction_equiv_(mesh,v,p,spower,mpower)
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// % Purpose
// % =======
// % Evaluate level set function using R-functions (equivalence): rho^1 ^
// % rho^2 . . . . ^ rho^n-1 where rho^i is the edge-weight for the i-th edge
// % USAGE: [D, dD] = Rfunction_equiv(v,p,in)
// %
// % Inputs:
// % v  : [x1 y1; x2 y2; ...; xn yn], n vertices of the polygon in ccw (Counter Clock Wise); 
// %      polygon is open (n vertices) or closed (n+1 vertices; repeat 1st vertex)
// % p  : [ px(1) py(1); ... ; px(m) py(m) ], the m points for which we want
// %      to compute the levet set function and its derivatives
// %
// % Outputs:
// % D  : level set function = [D_1; ...; D_m] (m x 1 array)
// % dD : gradient of level set function = [dD_1; ...; dD_m] (m x 2 array)
// %
// % N. Sukumar: UC Davis, M. Arroyo UPC-BarcelonaTech, February 2014
//
// MEX FUNCTION made by Daniel Millan, UPC-BarcelonaTech, February 2014

#include "mex.h"
#include <cmath>


// computes the absolute value of a number
#define ABS(x) ( ((x) < 0) ? -(x) : (x) )

// The  least  positive  number  that added to 1 returns a number that is greater than 1
#define EPSILON         2.220446049250313e-16
        
#define v(i,j)       v[(i)+vPts*(j)]
#define p(i,j)       p[(i)+pPts*(j)]
#define dD(i,j)     dD[(i)+pPts*(j)]

#define mesh(i,j) mesh[(i)+nEdg*(j)]
#define drho(i,j) drho[(i)+nEdg*(j)]

// template<class T> inline T sqr(T x) { return x*x; }
// template<class T> inline T dot2(T *a,T *b) { return a[0]*b[0]+a[1]*b[1]; }
// template<class T> inline T length(T x,T y) { return sqrt(sqr(x)+sqr(y)); }

// =======================================================================================
// Compute the nonnegative edge weight and its derivative of the polygon 
// using the triangle inequality
// [rho,drho] = rhoweight(mesh,v,x)
inline double rhoweight(int nEdg, double *mesh, int vPts, double *v, double *smp, double *rho, double *drho)
{
  int    i, i1, i2;
  double xx, yy;
  double v1x, v1y, v2x, v2y;
  double s1x, s1y,  sx,  sy;
  double d, id;
  double vcx, vcy;
  double f, dfx, dfy, t, dtx, dty;
  double val1, val2, ival1, f2, f3;
  
  double rho0 = 1.0; //value to return
  
  xx = smp[0]; //sample point coordinate in x-axis
  yy = smp[1]; //sample point coordinate in y-axis
  
  for ( i=0; i<nEdg ; i++ ) // i-th edge
  {
    i1 = (int)mesh(i,0)-1;  //from Fortran to C index
    i2 = (int)mesh(i,1)-1;
    
    v1x = v(i1,0);
    v1y = v(i1,1);
    v2x = v(i2,0);
    v2y = v(i2,1);

    //  s1 = smp - v1;
    //  s  = v2  - v1;
    s1x = xx  - v1x;        s1y = yy  - v1y;
    sx  = v2x - v1x;        sy  = v2y - v1y;

    //  d  = norm(s);
    //  id = 1/d;
    d   = sqrt(sx*sx + sy*sy);
    id  = 1.0/d;
    
    //  vc = (v1+v2)*0.5;
    vcx = (v1x+v2x)*0.5;
    vcy = (v1y+v2y)*0.5;
    
    //  f  = (s1(1)*s(2) - s1(2)*s(1))*id;
    //  df = [s(2), -s(1)]*id;
    f   = (s1x*sy - s1y*sx)*id;
    dfx =  sy*id;
    dfy = -sx*id;

    //  t  = ( (0.5*d)^2 - sum((x-vc).^2) )*id;
    //  dt = -2*(x-vc)*id;
    t   = ( 0.25*d*d - ((xx-vcx)*(xx-vcx) + (yy-vcy)*(yy-vcy)) )*id;
    dtx = -2.0*(xx-vcx)*id;
    dty = -2.0*(yy-vcy)*id;
    
    //   val = sqrt(t^2+f^4);
    f2   = f*f;
    f3   = f*f2;
    val1 = sqrt(t*t+f2*f2);
    ival1= 1.0/val1;
    val2 = val1-t;
    
    //   rho (i)   = sqrt(f^2 + 0.25*(val-t)^2 );
    rho[i]    = sqrt(f2 + 0.25*val2*val2 );

    if ( ABS(rho[i]) < 2.0*EPSILON ) //Spacing of floating point numbers: 2^(-52) ~ 2.22e-16
    {
      rho[i]    =  0.0;
      drho(i,0) = -dfx;
      drho(i,1) = -dfy;
      
      rho0      = 0.0;
    }
    else
    {
      //   drho(i,:) = ( f*df + 0.25*(val-t) * ((t*dt+2*f^3*df)/val - dt) )/rho(i);
      drho(i,0) = ( f*dfx + 0.25*val2 * ((t*dtx + 2.0*f3*dfx)*ival1 - dtx) )/rho[i];
      drho(i,1) = ( f*dfy + 0.25*val2 * ((t*dty + 2.0*f3*dfy)*ival1 - dty) )/rho[i];
    }
  }

  return rho0;
}

// =======================================================================================
// Compute the R function and its gradient
// function [w,dw] = formR(rho,drho,p)
inline void formR(int nEdg, double mpow, double rho0, double *rho, double *drho, double *w, double *dw)
{
  int    i;
  int    nrho0;
  double iw, dwx, dwy;
  double irho, irhop;

  *w  = 0.0;
  dwx = 0.0;
  dwy = 0.0;

  if ( rho0 > 0.0 )
  {
    for (i=0; i<nEdg; i++ )
    {
      irho  = 1.0/rho[i];
      irhop = pow(irho,mpow);

      //  w  =  w + (1/rho(i))^p;
      *w   =  *w  + irhop;

      //  dw = dw + rho(i)^(-p-1)*drho(i,:);
      dwx = dwx + irho*irhop*drho(i,0);
      dwy = dwy + irho*irhop*drho(i,1);
    }

    // ws  = sum(1./(rho.^p));
    // w   = ws^(-1/p);
    iw    = 1.0/(*w);
    *w    = pow(iw,1.0/mpow);

    // dw  = (rho.^(-p-1))' * drho * w/ws;
    dw[0] = dwx * (*w)*iw;
    dw[1] = dwy * (*w)*iw;   
  }
  else
  {
    nrho0 = 0;
    for (i=0; i<nEdg; i++ )
    {
      if ( rho[i] < 2.0*EPSILON )
      {
        dwx += drho(i,0);
        dwy += drho(i,1);
        nrho0++;
      }
    }
    dw[0] = dwx / ((double)nrho0);
    dw[1] = dwy / ((double)nrho0);
  }//end if
  
  return;
}

// ================================MAIN FUNCTION==========================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //Input
  int    nEdg  = mxGetM(prhs[0]);  //number of edges
  int    vPts  = mxGetM(prhs[1]);  //number of vertices of the mesh
  int    pPts  = mxGetM(prhs[2]);  //number of points where Rfunctions are evaluated
  
  double *mesh = mxGetPr(prhs[0]);
  double *v    = mxGetPr(prhs[1]);
  double *p    = mxGetPr(prhs[2]);
  double  spow = mxGetPr(prhs[3])[0];
  double  mpow = mxGetPr(prhs[4])[0];
  
  //printf("\tRfunction (equivalence)  spow=%4.2f  mpow=%4.2f\n", spow, mpow);
  //printf("\tNumber of edges  nEdg=%d\n", nEdg);
  //printf("\tNumber of verts  vPts=%d\n", vPts);
  //printf("\tNumber of points pPts=%d\n", pPts);
  
  //Output
  plhs[0]      = mxCreateDoubleMatrix(pPts,1,mxREAL);
  plhs[1]      = mxCreateDoubleMatrix(pPts,2,mxREAL);
  double *D    = mxGetPr(plhs[0]);
  double *dD   = mxGetPr(plhs[1]);
  
  mxArray *rho_ = mxCreateDoubleMatrix(nEdg,1,mxREAL); //N edges
  mxArray *drho_= mxCreateDoubleMatrix(nEdg,2,mxREAL); //N edges
  double *rho   = mxGetPr(rho_);
  double *drho  = mxGetPr(drho_);
  
  int     i, j, i1, i2;
  double  w_pow, w;
  double  dw[2];
  double  x[2], *v1, *v2;
  double rho0;
  
  //main loop on input points
  for ( j=0; j<pPts; j++ )
  {
    x[0] = p(j,0);
    x[1] = p(j,1);
    
    //Compute the nonnegative edge weight and its derivative
    rho0 = rhoweight(nEdg, mesh, vPts, v, x, rho, drho);
    
    // Compute the R function and its gradient
    formR(nEdg, mpow, rho0, rho, drho, &w, dw);
    
    w_pow  = pow(w,spow-1.0);
    D[j]   = w_pow*w;
    
    dD(j,0) = spow*w_pow*dw[0];
    dD(j,1) = spow*w_pow*dw[1];
    
    
    // checking NaN values at the vertices of the polygon
    if ( std::isnan(dD(j,0)) > 0)
    {
       //printf("\t--Message dDx is NaN\n");
       dD(j,0) = 0.0;
     }
    if ( std::isnan(dD(j,1)) > 0)
    {
       //printf("\t--Message dDy is NaN\n");
       dD(j,1) = 0.0;
     }
  }
  
  //free temp memory to avoid memory leaks
  mxDestroyArray(rho_);
  mxDestroyArray(drho_);
  
  return;
}
