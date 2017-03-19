#include "redist_unsigned.hpp"

RedistU::RedistU(const Array2D<double> &_u, const int _width, const int _flag) :
  width(_width),
  flag(_flag),
  m(_u.getm()),
  n(_u.getn()),
  N(m*n),
  dx(1.0f/static_cast<double>(n)),
  dy(1.0f/static_cast<double>(m)),
  thres(static_cast<double>(width+1)*mymax(dx,dy)),
  cpflag((_flag==2) || (_flag==3)),
  h(10*width,cpflag),
  state(m,n),
  u0(_u,_flag),
  u(_u),
  cpx(m,n),
  cpy(m,n)
{ }

RedistU::RedistU(const Array2D<double> &_u, const double _dx, const double _dy, const int _width, const int _flag) :
  width(_width),
  flag(_flag),
  m(_u.getm()),
  n(_u.getn()),
  N(m*n),
  dx(_dx),
  dy(_dy),
  thres(static_cast<double>(width+1)*mymax(dx,dy)),
  cpflag((_flag==2) || (_flag==3)),
  h(10*width,cpflag),
  state(m,n),
  u0(_u,_flag),
  u(_u),
  cpx(m,n),
  cpy(m,n)
{ }

void RedistU::redistance()
{
  if((flag == 0) || (flag == 1))
  {
    fastMarchingRedist();
    //if(flag == 0)
    //  secondOrderIterations();
  }
  else if((flag == 2) || (flag == 3))
    directionalOptimization();
  else
    cout << "Do not recognize the flag. Cowardly refusing to do anything." << endl;
}

void RedistU::fastMarchingRedist()
{
  setInterfaceValues();
  thresholdAwayFromInterface();
  // For now, cheat: set interface values as usual; only work on the outer ring
  for(size_t ii=0;ii<bnd.size();++ii)
    updateAndAddNeighborsToHeap(bnd[ii]);

  // while heap is non-empty, fix its top value, and update the neighbors
  while(h.numberOfElements() > 0)
  {
    struct helt h1 = h.popFromHeap();
    if((state.get(h1.i) == false) && (fabs(h1.d) < thres)) // don't update if outside threshold cutoff
    {
      applyResult(h1);
      state.put(true,h1.i);
      updateAndAddNeighborsToHeap(h1.i);
    }
  }
}

void RedistU::updateAndAddNeighborsToHeap(const int idx)
{ // fix the value at u[idx], update its 4 neighbors if they are not already fixed, and make sure they're in the heap at the appropriate location
  int idx2arr[4];
  idx2arr[0] = u.xp(idx);
  idx2arr[1] = u.xm(idx); 
  idx2arr[2] = u.yp(idx);
  idx2arr[3] = u.ym(idx);
  for(int ii=0;ii<4;++ii)
  {
    if(state.get(idx2arr[ii]) == false) // if true, value is already fixed
    {
      double dtemp = estimateUpdate(idx2arr[ii]);
      if(fabs(dtemp) < fabs(u[idx2arr[ii]]))
      {
        u[idx2arr[ii]] = dtemp; 
        h.addToHeap(idx2arr[ii],fabs(dtemp));
      }
    }
  }
}

void RedistU::updateAndAddNeighborsToHeapDO(const int idx)
{ // fix the value at u[idx], update its 4 neighbors if they are not already fixed, and make sure they're in the heap at the appropriate location
  int idx2arr[4];
  idx2arr[0] = u0.xp(idx);
  idx2arr[1] = u0.xm(idx); 
  idx2arr[2] = u0.yp(idx);
  idx2arr[3] = u0.ym(idx);
  for(int ii=0;ii<4;++ii)
  {
    if(state.get(idx2arr[ii]) == false) // if true, value is already fixed
    {
      struct helt htemp = performDO_unsigned(idx2arr[ii],cpx[idx],cpy[idx]);
      if(fabs(htemp.d) < fabs(u[idx2arr[ii]]))
      {
        applyResult(htemp);
        h.addToHeap(idx2arr[ii],fabs(htemp.d),htemp.aux);
      }
    }
  }
}

double RedistU::estimateUpdate(const int idx)
{
  const double a = mymin(fabs(u.getxm(idx)),fabs(u.getxp(idx)));
  const double b = mymin(fabs(u.getym(idx)),fabs(u.getyp(idx)));
  if(((a-b)*(a-b))>=(dx*dx+dy*dy))
    if((a+dx) < (b+dy))
      if(fabs(u.getxm(idx)) < fabs(u.getxp(idx)))
	return(u.getxm(idx) + dx*mysign(u.getxm(idx)));
      else
	return(u.getxp(idx) + dx*mysign(u.getxp(idx)));
    else
      if(fabs(u.getym(idx)) < fabs(u.getyp(idx)))
	return(u.getym(idx) + dy*mysign(u.getym(idx)));
      else
	return(u.getyp(idx) + dy*mysign(u.getyp(idx)));
  else
    return((dy*dy*a+dx*dx*b+dx*dy*sqrt(dx*dx+dy*dy-(a-b)*(a-b)))/(dx*dx+dy*dy));
}

void RedistU::setInterfaceValues()
{ 
  Array2D<int> sgn(u.getm(),u.getn());
  for(int ii=0; ii<sgn.getN(); ++ii)
    sgn[ii] = static_cast<int>(mysign(u[ii]));
  for(int ii=0; ii<sgn.getN(); ++ii)
    if((abs(sgn.getxp(ii)-sgn[ii]) + abs(sgn.getxm(ii)-sgn[ii]) +
        abs(sgn.getyp(ii)-sgn[ii]) + abs(sgn.getym(ii)-sgn[ii]))>0)
      bnd.push_back(ii);
 
  vector<double> dr(bnd.size());

  for(size_t ii=0; ii<bnd.size(); ++ii)
  {
    // compute norm(grad u) with centered differences
    double rx = (u.getxp(bnd[ii])-u.getxm(bnd[ii]))/dx/2.0f;
    double ry = (u.getyp(bnd[ii])-u.getym(bnd[ii]))/dy/2.0f;
    dr[ii] = sqrt(rx*rx+ry*ry);

    // compute norm(grad u) with one-sided differences
    rx = mymax(fabs(u.getxp(bnd[ii])-u[bnd[ii]]),fabs(u[bnd[ii]]-u.getxm(bnd[ii])))/dx;
    ry = mymax(fabs(u.getyp(bnd[ii])-u[bnd[ii]]),fabs(u[bnd[ii]]-u.getym(bnd[ii])))/dy;
    const double dr2 = sqrt(rx*rx+ry*ry);

    // Accept one-sided difference is much different than centered difference
    if((dr[ii] < (0.5*dr2)) || (dr[ii] > (2.0*dr2)))
      dr[ii] = dr2;
  }
  for(size_t ii=0; ii<bnd.size();++ii)
    u[bnd[ii]] = u[bnd[ii]]/dr[ii];
}

void RedistU::setInterfaceValuesDO_unsigned(const double cutoff)
{ 
  double meanabsbndval = 0.0f;
  for(int ii=0; ii<N; ++ii)
    if(u[ii] <= cutoff)
    {
      bnd.push_back(ii);
      meanabsbndval += fabs(u0[ii]);
    }
  meanabsbndval /= static_cast<double>(bnd.size());
  // rescale u0 based on meanabsbndval
  u0 *= (mymax(dx,dy) / meanabsbndval);

  // compute interface values and update
  struct helt bndval;
  for(size_t ii=0; ii<bnd.size(); ++ii)
  {
    bndval = performDOSurf_unsigned(bnd[ii]);
    applyResult(bndval);
    state.put(true,bndval.i);
  }
}

void RedistU::thresholdAwayFromInterface()
{
  for(int ii=0;ii<state.getN();++ii)
    state.put(false,ii);
  for(vector<int>::iterator it=bnd.begin(); it != bnd.end(); ++it)
  {
    state.put(true,*it); // true indicates value is fixed (false otherwise)
  }
  for(int ii=0;ii<state.getN();++ii)
    if(!state.get(ii))
      u[ii] = mysign(u[ii])*thres;
}

void RedistU::directionalOptimization()
{
  setInterfaceValuesDO_unsigned(4.0f*mymax(dx,dy)*mymax(dx,dy)); // (2dx)^2 
  thresholdAwayFromInterface();

  for(size_t ii=0;ii<bnd.size();++ii)
    updateAndAddNeighborsToHeapDO(bnd[ii]);
  while(h.numberOfElements() > 0)
  {
    struct helt h1 = h.popFromHeap();
    if((state.get(h1.i)==false) && (fabs(h1.d) < thres)) // don't update if outside threshold cutoff, or already fixed
    {
      applyResult(h1);
      state.put(true,h1.i);
      updateAndAddNeighborsToHeapDO(h1.i);
    }
  }
}

void RedistU::applyResult(const struct helt &h)
{
  u[h.i] = mysign(u[h.i])*h.d;
  if(flag > 1)
  {
    cpx[h.i] = h.aux[0];
    cpy[h.i] = h.aux[1];
  }
}

struct helt RedistU::performDOSurf_unsigned(const int idx)
{ // make a guess for interface location and call performDO(const int ind, const double *cpguess)
  struct helt h; 
  // make initial, gradient-based guess.
  double grad0[2], cpguess[2], grad[2];
  // compute "downwind" gradient (take informatiom from *larger* values of u to avoid crossing the minimum)
  if(u0.getxp(idx) > u0.getxm(idx))
    grad0[0] = pdl(u0.getxp(idx),u0[idx],u0.lenx()) / dx;
  else
    grad0[0] = pdl(u0[idx],u0.getxm(idx),u0.lenx()) / dx;
  if(u0.getyp(idx) > u0.getym(idx))
    grad0[1] = pdl(u0.getyp(idx),u0[idx],u0.leny()) / dy;
  else
    grad0[1] = pdl(u0[idx],u0.getym(idx),u0.leny()) / dy;

  grad[0] = grad0[0];
  grad[1] = grad0[1]; // grad is overwritten in lineSearch

  lineSearch_unsigned(idx,grad,cpguess);
  h = performDO_unsigned(idx,cpguess[0],cpguess[1]);
  double grad2all[16] = {1.0f,0.0f, 0.0f,1.0f, -1.0f,0.0f, 0.0f,-1.0f, 1.0f,1.0f, -1.0f,1.0f, -1.0f,-1.0f, 1.0f,-1.0f};
  for(int ii = 0; ii<8; ++ii)
  {
    double grad2[2]; grad2[0] = grad2all[2*ii]; grad2[1] = grad2all[2*ii+1];
    if(myDot(grad2[0],grad2[1],grad0[0],grad0[1]) >= 0.0f)
    {
      double cpguess2[2];
      lineSearch_unsigned(idx,grad2,cpguess2);
      struct helt h1 = performDO_unsigned(idx,cpguess2[0],cpguess2[1]);
      if( fabs(h1.d) < fabs(h.d) )
        h = h1;
    }
  }
  return(h);
}
struct helt RedistU::performDO_unsigned(const int idx, const double cpxguess, const double cpyguess)
{
  // optimize over directions via line search
  const double x0[2] = {u0.getX(idx), u0.getY(idx)};
  const int MAXLOOPS = 5;
  const int MAXSEARCHES = 5;
  const double TOL = mymax( (0.01 / pow(mymax(m,n),flag+1)), 5.0e-16) ;

  double xl[2], xr[2], xc[2], xt[2];
  double dl,dc,dr,dt;
  xc[0] = cpxguess;
  xc[1] = cpyguess;

  dc = dist(x0,xc,u0.lenx(),u0.leny());
  // set angle increment
  double delta;
  if(dc > mymax(dx,dy))
    delta = asin(mymax(dx,dy)/dc);
  else // interface is nearby, look over a wide range
    delta = PI2;

  findNborDirections(x0,xc,xl,xr,delta);
  // bracket the interface in the search directions
  dl = search1D_unsigned(idx,xl);
  dr = search1D_unsigned(idx,xr);
  
  for(int count0=0; count0<MAXLOOPS; ++count0)
  { // (1) center on minimum; (2) perform Newton step
    int count = 0;
    while( (dc > mymin(dl,dr)) && (count < MAXSEARCHES))
    {
      if(delta < PI4-1e-4)
        delta *= 2.0f;

      if(dl < dr)
      {
        dc = dl; xc[0] = xl[0]; xc[1] = xl[1];
      }
      else // dl >= dr
      {
        dc = dr; xc[0] = xr[0]; xc[1] = xr[1];
      }
      findNborDirections(x0,xc,xl,xr,delta);
      dl = search1D_unsigned(idx,xl);
      dr = search1D_unsigned(idx,xr);
      count++;
    }
    // perform Newton step only if dl/dr are successfully computed
    if(mymax(dl,dr) < 0.9f*mymax(m,n)*mymax(dx,dy))
    {
      double ds = (dr-dl)/2.0f;
      double dss = (dr-2.0f*dc+dl);
      if(fabs(dss) < 1e-16)
        dss = 1e-16;
      double dir[2];
      dir[0] = pdl(xc[0],x0[0],u0.lenx()); dir[1] = pdl(xc[1],x0[1],u0.leny());
      double delta2 = -ds/dss * delta;
      xt[0] = dinrange2l(x0[0] + cos(delta2)*dir[0] - sin(delta2)*dir[1],u0.lenx());
      xt[1] = dinrange2l(x0[1] + sin(delta2)*dir[0] + cos(delta2)*dir[1],u0.leny());
      dt = search1D_unsigned(idx,xt);
      if(fabs(dt-dc) < TOL)
        count0 = MAXLOOPS; // no further improvement available
      else if(dt<dc)
      {
        xc[0] = xt[0]; xc[1] = xt[1]; dc = dt;
      }
    }
    if(count0 < MAXLOOPS)
    {
      delta /= 2.0f;
      findNborDirections(x0,xc,xl,xr,delta);
      dl = search1D_unsigned(idx,xl);
      dr = search1D_unsigned(idx,xr);
    }
  }
  struct helt lval;
  double val1 = u0.interpolate(xc[0],xc[1]);
  lval.i = idx;
  lval.d = dc + val1;
  lval.aux[0] = xc[0];
  lval.aux[1] = xc[1];
  return(lval);
}

void RedistU::lineSearch_unsigned(const int idx, double (&grad)[2], double (&cpguess)[2])
{ // takes the given index at estimate of gradient direction (unnormalized)
  // and performs a line search along this direction for the interface.
  // Interface is expected to lie within a few grid cells of idx
  // Return location goes in cpguess
  
  normalize(grad);
  const double dt = 1.0f; // If u0 ~= d^2, I think this is appropriate time step.
  double interpRes[3];
  // now: look along line defined by x - t*grad (t>0)
  double t = SQRT2*mymax(dx,dy);
  cpguess[0] = dinrange2l(u0.getX(idx)-t*grad[0],u0.lenx()); // initial guess
  cpguess[1] = dinrange2l(u0.getY(idx)-t*grad[1],u0.leny());
  double x1[3], x2[3];
  bool gl = bracket_unsigned(idx,cpguess,x1,x2);
  if(gl)
    double uval = bisect_unsigned(cpguess,x1,x2);
  else
  {
    for(int ii=0; ii<5; ++ii) // plain gradient descent to get near zero
    {
      u0.interpolate(cpguess[0],cpguess[1],interpRes);
      cpguess[0] = dinrange2l(cpguess[0]-dt*interpRes[1],u0.lenx()); // initial guess
      cpguess[1] = dinrange2l(cpguess[1]-dt*interpRes[2],u0.leny());
    }
  }
}

bool RedistU::bracket_unsigned(const int idx, double (&cpguess)[2], double (&x1)[3], double (&x2)[3])
{ // Look about location cpguess along line from idx and find small negative
  // and positive derivatives 
  // input to bisection.
  
  // Returns false if it's much further than initially guessed to 
  // interface in this direction, and true otherwise.
  // Assumed that x1, x2, both have room for two doubles

  const double x0[2] = {u0.getX(idx), u0.getY(idx)};
  double iRes0[3], iRes1[3], iRes2[3];
  double dir[2];
  dir[0] = pdl(cpguess[0],x0[0],u0.lenx());
  dir[1] = pdl(cpguess[1],x0[1],u0.leny());
  normalize(dir);
 
  u0.interpolate(x0[0],x0[1],iRes0);
  const double dd0 = dirDerivative(&(iRes0[1]),dir);
  if( (dd0 == 0) && (iRes0[0] <= mymax(dx,dy)/3.0f) )
    // already at local extrema; accept only if value is also small
  {
    x1[0] = x0[0]; x1[1] = x1[1]; x1[2] = dd0;
    x2[0] = x0[0]; x2[1] = x2[1]; x2[2] = dd0;
    return(true);
  }
  
  u0.interpolate(cpguess[0],cpguess[1],iRes1);
  const double dd1 = dirDerivative(&(iRes1[1]),dir);
  x1[0] = cpguess[0]; x1[1] = cpguess[1]; x1[2] = dd1;
  if(mysign(dd1) == 0) // local extrema at proposed closest point, accept.
  {
    x2[0] = cpguess[0]; x2[1] = cpguess[1]; x2[2] = dd1;
    return(true);
  }
  else
  { // work in direction of decreasing d
    int tries = 0;
    double dd2;
    do
    {
      ++tries;
      const double xguess = dinrange2l(cpguess[0]-dir[0]*mysign(dd1)*tries*mymax(dx,dy),u0.lenx());
      const double yguess = dinrange2l(cpguess[1]-dir[1]*mysign(dd1)*tries*mymax(dx,dy),u0.leny());

      u0.interpolate(xguess,yguess,iRes2);
      dd2 = dirDerivative(&(iRes2[1]),dir);
    }
    while((mysign(dd1 * dd2) == 1) && (tries < 10));

    if(mysign(dd1 * dd2) == -1)
    {
      x2[0] = dinrange2l(cpguess[0]-dir[0]*mysign(dd1)*tries*mymax(dx,dy),u0.lenx());
      x2[1] = dinrange2l(cpguess[1]-dir[1]*mysign(dd1)*tries*mymax(dx,dy),u0.leny());
      x2[2] = dd2;
      return(true);
    }
    else
      return(false);
  }
}

double RedistU::bisect_unsigned(double (&result)[2], double (&x1)[3], double (&x2)[3])
{ // uses tangent line approximation to directional derivative to try to 
  // find near-extrema of u between x1 and x2
  // Assumes that the directional derivatives at x1 and x2 along their common line have different signs.
  double u1 = x1[2]; // directional derivative of u at (x1[0],x1[1])
  double u2 = x2[2]; // directional derivative of u at (x2[0],x2[1])
  double xt[2];
  assert(mysign(u1*u2) <= 0);
  if(u1==u2)
  {
    result[0] = x1[0]; result[1] = x1[1];
    return(u1);
  }
  ccomb(x1,x2,xt,u2/(u2-u1),u0.lenx(),u0.leny());
  double ut = u0.interpolate(xt[0],xt[1]); 
  result[0] = xt[0]; result[1] = xt[1]; 
  return(ut);
}

void RedistU::findNborDirections(const double (&x0)[2], const double (&xc)[2], double (&xl)[2], double (&xr)[2], const double delta)
{
  double s,c;
  double rx,ry;
  s = sin(delta);
  c = cos(delta);
  rx = pdl(xc[0],x0[0],u0.lenx());
  ry = pdl(xc[1],x0[1],u0.leny());
  xl[0] = dinrange2l(x0[0] + c*rx + s*ry,u0.lenx());
  xl[1] = dinrange2l(x0[1] - s*rx + c*ry,u0.leny());
  xr[0] = dinrange2l(x0[0] + c*rx - s*ry,u0.lenx());
  xr[1] = dinrange2l(x0[1] + s*rx + c*ry,u0.leny());
}

double RedistU::search1D_unsigned(const int idx, double (&x)[2])
{
  // bracket the interface in the search directions
  double  x0[2], xm[3], xp[3]; // xm and xp also return values/gradients from u0
  bool gl = bracket_unsigned(idx,x,xm,xp);
  double d;

  x0[0] = u0.getX(idx); x0[1] = u0.getY(idx);
  // once bracketed successfully, bisect
  if(gl)
  {
    double val1 = bisect_unsigned(x,xm,xp);
    d = dist(x0,x,u0.lenx(),u0.leny());
    d += val1;
  }
  else
    d = mymax(m,n)*mymax(dx,dy);
  return(d);
}

void RedistU::dump_u(double *v)
{ // assumes sufficient memory is allocated into v
  for(int ii=0; ii<N; ++ii)
    v[ii] = u[ii];
}

void RedistU::dump_cp(double *cpxArr, double *cpyArr)
{ // assume sufficient memory is allocated into cpxArr, cpyArr
  for(int ii=0; ii<N; ++ii)
    cpxArr[ii] = cpx[ii];
  for(int ii=0; ii<N; ++ii)
    cpyArr[ii] = cpy[ii];
}

