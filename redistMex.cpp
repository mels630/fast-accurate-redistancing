#include "mex.h"
#include "redist.hpp"
#include "redist_unsigned.hpp"
#include "redist3.hpp"
#include "defs.h"
#include "array2d.hpp"
#include "array3d.hpp"
#include "toolbox.hpp"
#include "toolbox3d.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Matlab interface:
       [v,cpx,cpy] = redist(u,width,flag)
       where u is m x n 2D level set function to be redistanced,
             width is half-width of tubular region for redistancing,
  */

  double *u, *v, *cpx, *cpy, *cpz;
  double dx, dy, dz;
  int width, flag, m, n, k, N, ii;
  int *bnd;
  /* u holds the input array */
  u = mxGetPr(prhs[0]);
  width = static_cast<int>(mxGetScalar(prhs[1]));
  flag = static_cast<int>(mxGetScalar(prhs[2]));
  k = static_cast<int>(mxGetScalar(prhs[3]));
  dx = mxGetScalar(prhs[4]);
  dy = mxGetScalar(prhs[5]);
  dz = mxGetScalar(prhs[6]);
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  N = m*n;
  
  if(k == 1)
  {
    plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m,n,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(m,n,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
    v = mxGetPr(plhs[0]);
    cpx = mxGetPr(plhs[1]);
    cpy = mxGetPr(plhs[2]);
    // copy data into Array2D<double> type
    Array2D<double> uu(u,m,n,dx,dy);
    if(flag != 4)
    {
      Redist r(uu,dx,dy,width,flag);
      r.redistance();
      // handle outputs
      r.dump_u(v);
      if((flag==2) || (flag==3))
        r.dump_cp(cpx,cpy);
    }
    else // flag == 4 -- unsigned redistancing
    {
      // make uu a (quadratic) unsigned level set function with
      // 0-level set giving the interface
      for(int ii=0; ii<uu.getN(); ++ii)
        uu[ii] *= uu[ii];
      RedistU r(uu,dx,dy,width,3);
      r.redistance();
      // handle outputs
      r.dump_u(v);
      r.dump_cp(cpx,cpy);
    }
  }
  else
  {
    n /= k;
    plhs[0] = mxCreateDoubleMatrix(m,n*k,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m,n*k,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(m,n*k,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(m,n*k,mxREAL);
    v = mxGetPr(plhs[0]);
    cpx = mxGetPr(plhs[1]);
    cpy = mxGetPr(plhs[2]);
    cpz = mxGetPr(plhs[3]);
    Array3D<double> uu(m,n,k);
    for(int ii=0; ii<N; ++ii)
      uu.put(u[ii],ii);
    if( (fabs(dx-1.0f/static_cast<double>(n)) > 5.e-16) || (fabs(dy-1.0f/static_cast<double>(m)) > 5.e-16) || (fabs(dz-1.0f/static_cast<double>(k)) > 5.e-16) )
      mexPrintf(" 3D code is not enabled for dx/dy/dz other than 1/grid length yet.\n");
    Redist3 r(uu,width,flag);
    r.redistance();
    // handle outputs
    r.dump_u(v);
    if((flag==2) || (flag==3))
      r.dump_cp(cpx,cpy,cpz);
  }
}

