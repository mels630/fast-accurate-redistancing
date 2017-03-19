#include "mex.h"
#include "reinit.hpp"
#include "defs.h"
#include "array2d.hpp"
#include "toolbox.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Matlab interface:
  //   [v,cpx,cpy] = redist(tpx,tpy,width,flag)
  //   where (tpx,tpy) is the track-point function to be reinitialized
  //   width is half-width of tubular region for reinitialization (in grid cells)

  const double *_tpx = mxGetPr(prhs[0]);
  const double *_tpy = mxGetPr(prhs[1]);
  const int width = static_cast<int>(mxGetScalar(prhs[2]));
  const int flag = static_cast<int>(mxGetScalar(prhs[3]));
  const int m = mxGetM(prhs[0]);
  const int n = mxGetN(prhs[0]);
  const int N = m*n;
  
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m,n,mxREAL);
  double *cpx = mxGetPr(plhs[0]);
  double *cpy = mxGetPr(plhs[1]);
  Reinit r(Array2D<double>(_tpx,m,n),Array2D<double>(_tpy,m,n),width,flag);
  r.reinit();
  // handle outputs
  r.dump_cp(cpx,cpy);
}

