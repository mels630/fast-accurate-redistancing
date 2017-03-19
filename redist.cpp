#include "redist.hpp"

Redist::Redist(const Array2D<double> &_u, const int _width, const int _flag) :
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

Redist::Redist(const Array2D<double> &_u, const double _dx, const double _dy, const int _width, const int _flag) :
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
  u0(_u.duplicateArray2D(),_flag),   // should be more careful in memory
  u(_u.duplicateArray2D()),           // at some point!!!
  cpx(m,n),
  cpy(m,n)
{ }

void Redist::redistance()
{
  if((flag == 0) || (flag == 1))
  {
    fastMarchingRedist();
    if(flag == 0)
      secondOrderIterations();
  }
  else if((flag == 2) || (flag == 3))
    directionalOptimization();
  else
    cout << "Do not recognize the flag. Cowardly refusing to do anything." << endl;
}

void Redist::fastMarchingRedist()
{
  setInterfaceValues();
  thresholdAwayFromInterface();

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

void Redist::updateAndAddNeighborsToHeap(const int idx)
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

void Redist::updateAndAddNeighborsToHeapDO(const int idx)
{ // fix the value at u[idx], update its 4 neighbors if they are not already fixed, and make sure they're in the heap at the appropriate location
  int idx2arr[4];
  idx2arr[0] = u.xp(idx);
  idx2arr[1] = u.xm(idx); 
  idx2arr[2] = u.yp(idx);
  idx2arr[3] = u.ym(idx);
  //idx2arr[4] = u.xp(u.yp(idx));
  //idx2arr[5] = u.xm(u.yp(idx));
  //idx2arr[6] = u.xp(u.ym(idx));
  //idx2arr[7] = u.xm(u.ym(idx));
  for(int ii=0;ii<4;++ii)
  {
    if(state.get(idx2arr[ii]) == false) // if true, value is already fixed
    {
      struct helt htemp = performDO(idx2arr[ii],cpx[idx],cpy[idx]);
      if(fabs(htemp.d) < fabs(u[idx2arr[ii]]))
      {
        applyResult(htemp);
        h.addToHeap(idx2arr[ii],fabs(htemp.d),htemp.aux);
      }
    }
  }
}

double Redist::estimateUpdate(const int idx)
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
    return(mysign(u[idx]) * (dy*dy*a+dx*dx*b+dx*dy*sqrt(dx*dx+dy*dy-(a-b)*(a-b)))/(dx*dx+dy*dy));
}

void Redist::setInterfaceValues()
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
    double ry = (u.getyp(bnd[ii])-u.getym(bnd[ii]))/dy/2.0f;
    double rx = (u.getxp(bnd[ii])-u.getxm(bnd[ii]))/dx/2.0f;
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

void Redist::setInterfaceValuesDO()
{ 
  double meanabsbndval = 0.0f;
  for(int ii=0; ii<N; ++ii)
    if(diffSign(ii))
    {
      bnd.push_back(ii);
      meanabsbndval += fabs(u0[ii]);
    }
  meanabsbndval /= static_cast<double>(bnd.size());
  // rescale u0 based on meanabsbndval
  u0 *= 0.5f * (mymax(dx,dy) / meanabsbndval); // scale u0 (rationale: mean distance to interface at a grid cell should be ~ 0.5 dx)

  // compute interface values and update
  struct helt bndval;
  for(size_t ii=0; ii<bnd.size(); ++ii)
  {
    bndval = performDOSurf(bnd[ii]);
    applyResult(bndval);
    state.put(true,bndval.i);
  }
}

void Redist::thresholdAwayFromInterface()
{
  for(int ii=0;ii<state.getN();++ii)
    state.put(false,ii);
  for(vector<int>::iterator it=bnd.begin(); it != bnd.end(); ++it)
  {
    state.put(true,*it); // true indicates value is fixed (false otherwise)
  }
  for(int ii=0;ii<state.getN();++ii)
    if(!state.get(ii))
      u[ii] = mysign(u0[ii])*thres;
}

void Redist::directionalOptimization()
{
  setInterfaceValuesDO(); 
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

void Redist::applyResult(const struct helt &h)
{
  u[h.i] = mysign(u0[h.i])*h.d;
  if(flag > 1)
  {
    cpx[h.i] = h.aux[0];
    cpy[h.i] = h.aux[1];
  }
}

struct helt Redist::performDO(const int idx)
{ // make a guess for interface location and call performDO(const int ind, const double *cpguess)
  double grad[2], cpguess[2];
  grad[0] = pdl(u0.getxp(idx),u0.getxm(idx),u0.lenx()) / (2.0f * dx);
  grad[1] = pdl(u0.getyp(idx),u0.getym(idx),u0.leny()) / (2.0f * dy);
  lineSearch(idx,grad,cpguess);
  return(performDO(idx,cpguess[0],cpguess[1]));
}

struct helt Redist::performDOSurf(const int idx)
{ // make a guess for interface location and call performDO(const int ind, const double *cpguess)
  struct helt h; 
  // make initial, gradient-based guess.
  double grad0[2], cpguess[2], grad[2];
  grad0[0] = pdl(u0.getxp(idx),u0.getxm(idx),u0.lenx()) / (2.0f * dx);
  grad0[1] = pdl(u0.getyp(idx),u0.getym(idx),u0.leny()) / (2.0f * dy);

  grad[0] = grad0[0];
  grad[1] = grad0[1]; // grad is overwritten in lineSearch

  lineSearch(idx,grad,cpguess);
  h = performDO(idx,cpguess[0],cpguess[1]);

  double grad2all[16] = {1.0f,0.0f, 0.0f,1.0f, -1.0f,0.0f, 0.0f,-1.0f, 1.0f,1.0f, -1.0f,1.0f, -1.0f,-1.0f, 1.0f,-1.0f};
  for(int ii = 0; ii<8; ++ii)
  {
    double grad2[2]; grad2[0] = grad2all[2*ii]; grad2[1] = grad2all[2*ii+1];
    if(myDot(grad2[0],grad2[1],grad0[0],grad0[1]) >= 0.0f)
    {
      double cpguess2[2];
      lineSearch(idx,grad2,cpguess2);
      struct helt h1 = performDO(idx,cpguess2[0],cpguess2[1]);
      if( fabs(h1.d) < fabs(h.d) )
        h = h1;
    }
  }
  return(h);
}

struct helt Redist::performDO(const int idx, const double cpxguess, const double cpyguess)
{
  // optimize over directions via line search
  const double x0[2] = {u.getX(idx), u.getY(idx)};
  const int MAXLOOPS = 5;
  const int MAXSEARCHES = 5;
  // still not sure about this TOL!
  const double TOL = mymax( (0.01 * mymin(dx,dy) / pow(static_cast<double>(mymax(m,n)),flag)), 5.0e-16);
  // const double TOL = mymax( (0.01 / pow(mymax(m,n),flag+1)), 5.0e-16);

  double xl[2], xr[2], xc[2], xt[2];
  double dl,dc,dr,dt;
  xc[0] = cpxguess;
  xc[1] = cpyguess;

  dc = dist(x0,xc,u.lenx(),u.leny());
  // set angle increment
  double delta;
  if(dc > mymax(dx,dy))
    delta = mymax(asin(mymax(dx,dy)/dc),PI/32.0f);
  else // interface is nearby, look over a wide range
    delta = PI2;

  findNborDirections(x0,xc,xl,xr,delta);
  // bracket the interface in the search directions
  dl = search1D(idx,xl);
  dr = search1D(idx,xr);

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
      dl = search1D(idx,xl);
      dr = search1D(idx,xr);
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
      dir[0] = pdl(xc[0],x0[0],u.lenx()); dir[1] = pdl(xc[1],x0[1],u.leny());
      double delta2 = -ds/dss * delta;
      xt[0] = dinrange2l(x0[0] + cos(delta2)*dir[0] - sin(delta2)*dir[1],u.lenx());
      xt[1] = dinrange2l(x0[1] + sin(delta2)*dir[0] + cos(delta2)*dir[1],u.leny());
      dt = search1D(idx,xt);

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
      dl = search1D(idx,xl);
      dr = search1D(idx,xr);
    }
  }
  struct helt lval;
  double val1 = u0.interpolate(xc[0],xc[1]);
  lval.i = idx;
  lval.d = dc + fabs(val1);

  double grad[2];
  grad[0] = pdl(xc[0],x0[0],u0.lenx());
  grad[1] = pdl(xc[1],x0[1],u0.leny());
  normalize(grad);
  lval.aux[0] = xc[0] + fabs(val1) * grad[0];
  lval.aux[1] = xc[1] + fabs(val1) * grad[1];
  return(lval);
}

void Redist::lineSearch(const int idx, double (&grad)[2], double (&cpguess)[2])
{ // takes the given index at estimate of gradient direction (unnormalized)
  // and performs a line search along this direction for the interface.
  // Interface is expected to lie within one grid cell of idx

  normalize(grad);
  // get the correct direction to search in; 
  if(u[idx] < 0.0f) // then we go uphill, keep gradient sign
  {
    grad[0] *= SQRT2*mymax(dx,dy);
    grad[1] *= SQRT2*mymax(dx,dy);
  }
  else // go downhill
  {
    grad[0] *= -SQRT2*mymax(dx,dy);
    grad[1] *= -SQRT2*mymax(dx,dy);
  }
  if(findOppSign(idx,grad))
  {
    cpguess[0] = grad[0]; cpguess[1] = grad[1];
    bisect(idx,cpguess);
  }
}

bool Redist::findOppSign(const int idx, double (&guess)[2])
{ // searches in the direction initially given by "guess" to find x s.t. u(x) has opposite sign
  // from from u(idx)
  // input suggests to look at u[y+guess], where y is spatial coordinates of idx. 
  // returns location x in guess.
  const int sgn = mysign(u0[idx]);
  double inguess[2];
  inguess[0] = guess[0];
  inguess[1] = guess[1];

  for(int ii=1; ii<=5; ++ii)
  {
    guess[0] = dinrange2l(u0.getX(idx)+static_cast<double>(ii)*inguess[0],u0.lenx());
    guess[1] = dinrange2l(u0.getY(idx)+static_cast<double>(ii)*inguess[1],u0.leny());

    double val1 = u0.interpolate(guess[0],guess[1]);
    if(sgn != mysign(val1))
      return(true);
  }
  return(false);
}

bool Redist::bracket(const int idx, double (&cpguess)[2], double (&xm)[3], double (&xp)[3])
{ // Look about location cpguess along line from idx and find (small) positive and negative values as
  // input to bisection.
  
  // Returns false if it's much further than initially guessed to interface in this direction (e.g. more
  // than one grid cell further), and true otherwise.
  // Assumed that xm, xp, both have room for two doubles

  const double x0[2] = {u.getX(idx), u.getY(idx)};
  const double uv0 = u0[idx];
  if(mysign(uv0) == 0) // already on interface
  {
    xm[0] = x0[0]; xm[1] = x0[1]; xm[2] = uv0;
    xp[0] = x0[0]; xp[1] = x0[1]; xp[2] = uv0;
    return(true);
  }

  double dir[2];
  dir[0] = pdl(cpguess[0],x0[0],u.lenx());
  dir[1] = pdl(cpguess[1],x0[1],u.leny());
  normalize(dir);

  double ug = u0.interpolate(cpguess[0],cpguess[1]);
  if(mysign(uv0 * ug) == 0) // sign(ug) == 0; 
  {
    xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = ug;
    xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = ug;
    return(true);
  }
  if(mysign(uv0 * ug) == 1) // same side of interface
  {
    double ug2 = u0.interpolate(cpguess[0]+dir[0]*mymax(dx,dy),cpguess[1]+dir[1]*mymax(dx,dy));
    if(mysign(uv0 * ug2) == 1) // still on same side of interface
      return(false);
    else
    {
      if(mysign(ug) == 0)
      {
        if(mysign(ug2) == 1)
        {
          xp[0] = cpguess[0]+dir[0]*mymax(dx,dy); xp[1] = cpguess[1]+dir[1]*mymax(dx,dy); xp[2] = ug2;
          xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = ug;
        }
        else
        {
          xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = ug;
          xm[0] = cpguess[0]+dir[0]*mymax(dx,dy); xm[1] = cpguess[1]+dir[1]*mymax(dx,dy); xm[2] = ug2;
        }
      }
      else if(mysign(ug) == 1)
      {
        xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = ug;
        xm[0] = cpguess[0]+dir[0]*mymax(dx,dy); xm[1] = cpguess[1]+dir[1]*mymax(dx,dy); xm[2] = ug2;
      }
      else // mysign(ug) == -1
      {
        xp[0] = cpguess[0]+dir[0]*mymax(dx,dy); xp[1] = cpguess[1]+dir[1]*mymax(dx,dy); xp[2] = ug2;
        xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = ug;
      }
      return(true);
    }
  }
  else // x0 and cpguess are on opposite sides of interface, work backwards
  {
    double xx,yy;
    xx = cpguess[0]-dir[0]*mymax(dx,dy);
    yy = cpguess[1]-dir[1]*mymax(dx,dy);
    double ug2 = u0.interpolate(xx,yy);
    int tries = 0;
    while((mysign(ug*ug2) == 1) && (tries < 10))
    {
      xx -= dir[0]*mymax(dx,dy);
      yy -= dir[1]*mymax(dx,dy);
      ug2 = u0.interpolate(xx,yy);
      tries++;
    }
    if(tries == 10)
      return(false);

    if(mysign(ug) == 1)
    {
      xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = ug;
      xm[0] = xx; xm[1] = yy; xm[2] = ug2;
    }
    else if(mysign(ug) == -1)
    {
      xp[0] = xx; xp[1] = yy; xp[2] = ug2;
      xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = ug;
    }
    else // mysign(ug) == 0
    {
      if(mysign(ug2) == 1)
      {
        xp[0] = xx; xp[1] = yy; xp[2] = ug2;
        xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = ug;
      }
      else
      {
        xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = ug;
        xm[0] = xx; xm[1] = yy; xm[2] = ug2;
      }
    }
    return(true);
  }
}

void Redist::bisect(const int idx, double (&guess)[2])
{ // Assumes that u0[idx] and u[guess] have different signs.
  // Tries to find the zero between [idx] and (guess) by simple linear approximation.  
  double xm[2], xp[2], xt[2];
  double um, up;
  assert(mysign(u0[idx]*u0.interpolate(guess[0],guess[1])) != 1);

  if(mysign(u0[idx]) == 0)
  {
    guess[0] = u.getX(idx);
    guess[1] = u.getY(idx);
    return;
  }
  else if(mysign(u0[idx]) == 1)
  {
    xp[0] = u.getX(idx); xp[1] = u.getY(idx); up = u0[idx];
    xm[0] = guess[0];    xm[1] = guess[1];    um = u0.interpolate(xm[0],xm[1]);
  }
  else
  {
    xm[0] = u.getX(idx); xm[1] = u.getY(idx); um = u0[idx];
    xp[0] = guess[0];    xp[1] = guess[1];    up = u0.interpolate(xp[0],xp[1]);
    if(up <= -1.0e-15)
      printf(" up = %e, um = %e\n",up,um);
  }

  if(up==um)
  {
    guess[0] = xp[0]; guess[1] = xp[1];
    return;
  }
  if(isnan(up/(up-um)))
  {
    printf(" up = %.3e, um = %.3e. Setting um = -1.0f.\n",up,um);
    um = -1.0f;
  }
  ccomb(xm,xp,xt,up/(up-um),u0.lenx(),u0.leny());
  assert(up>=0.0f); 
  assert(um<=0.0f);
  guess[0] = xt[0];
  guess[1] = xt[1];
  return;
}

double Redist::bisect(double (&result)[2], double (&xm)[3], double (&xp)[3])
{ // bisection with linear approximation
  // Assumes that u0[idx] and u[guess] have different signs.
  // Tries to find the zero between [idx] and (guess) by single linear approx.
  double up = xp[2]; //u0.interpolate(xp[0],xp[1]);
  double um = xm[2]; //u0.interpolate(xm[0],xm[1]);
  double xt[2];
  assert(up>=0.0f); 
  assert(um<=0.0f);
  if(up==um)
  {
    result[0] = xp[0]; result[1] = xp[1];
    return(up);
  }
  ccomb(xm,xp,xt,up/(up-um),u0.lenx(),u0.leny());
  double ut = u0.interpolate(xt[0],xt[1]); 
  result[0] = xt[0]; result[1] = xt[1]; 
  return(ut);
}

void Redist::findNborDirections(const double (&x0)[2], const double (&xc)[2], double (&xl)[2], double (&xr)[2], const double delta)
{
  double s,c;
  double rx,ry;
  s = sin(delta);
  c = cos(delta);
  rx = pdl(xc[0],x0[0],u.lenx());
  ry = pdl(xc[1],x0[1],u.leny());
  xl[0] = dinrange2l(x0[0] + c*rx + s*ry,u.lenx());
  xl[1] = dinrange2l(x0[1] - s*rx + c*ry,u.leny());
  xr[0] = dinrange2l(x0[0] + c*rx - s*ry,u.lenx());
  xr[1] = dinrange2l(x0[1] + s*rx + c*ry,u.leny());
}

double Redist::search1D(const int idx, double (&x)[2])
{
  // bracket the interface in the search directions
  double  x0[2], xm[3], xp[3]; // xm and xp also return values from u0
  bool gl = bracket(idx,x,xm,xp);
  double d;

  x0[0] = u.getX(idx); x0[1] = u.getY(idx);
  // once bracketed successfully, bisect
  if(gl)
  {
    double val1 = bisect(x,xm,xp);
    d = dist(x0,x,u.lenx(),u.leny());
    d += val1 * mysign(u0[idx]);
  }
  else
    d = mymax(m,n)*mymax(dx,dy);
  return(d);
}

void Redist::dump_u(double *v)
{ // assumes sufficient memory is allocated into v
  for(int ii=0; ii<N; ++ii)
    v[ii] = u[ii];
}

void Redist::dump_cp(double *cpxArr, double *cpyArr)
{ // assume sufficient memory is allocated into cpxArr, cpyArr
  for(int ii=0; ii<N; ++ii)
    cpxArr[ii] = cpx[ii];
  for(int ii=0; ii<N; ++ii)
    cpyArr[ii] = cpy[ii];
}

void Redist::secondOrderIterations()
{
  const double me = 1.0f/ 10000.0f * mymax(dx,dy);
  const double thres2 = thres - dx;
  const double dt = 1.0f / 5.0f * mymax(dx,dy);
  // Step 1: Figure out which indices to do the iterations on
  vector<int> idxc;
  for(int ii=0; ii<N; ++ii)
    if(fabs(u[ii]) < thres2)
      idxc.push_back(ii);
  
  double *G = new double[idxc.size()];
  bool *df = new bool[idxc.size()];
  for(size_t ii=0; ii<idxc.size(); ++ii)
    df[ii] = diffSign(idxc[ii]);

  const int iter2o = 3*width;
  for(int iter=1; iter<=iter2o; ++iter)
  {
    for(size_t ii=0; ii<idxc.size(); ++ii)
      if(!(df[ii]))
      {
        const int jj = idxc[ii];

        double t0[5],t1[4],t2[4],t3[4];
        /* x */
        t0[0] = u.getxm(u.xm(jj));
        t0[1] = u.getxm(jj);
        t0[2] = u[jj];
        t0[3] = u.getxp(jj);
        t0[4] = u.getxp(u.xp(jj));
        t1[0] = dx; t1[1] = dx; t1[2] = dx; t1[3] = dx;
        if((mysigntol(t0[0],me)*mysigntol(t0[1],me))==-1)
        {
          t1[0] = dx * t0[1] / (t0[1]-t0[0]);
          t0[0] = 0.0f;
        }
        if((mysigntol(t0[3],me)*mysigntol(t0[4],me))==-1)
        {
          t1[3] = -dx * t0[3] / (t0[4]-t0[3]);
          t0[4] = 0.0f;
        }
        
        t2[0] = (t0[1]-t0[0]) / t1[0];
        t2[1] = (t0[2]-t0[1]) / t1[1];
        t2[2] = (t0[3]-t0[2]) / t1[2];
        t2[3] = (t0[4]-t0[3]) / t1[3];
        t3[0] = (t2[1]-t2[0]) / (t1[1]+t1[0]);
        t3[1] = (t2[2]-t2[1]) / (t1[2]+t1[1]);
        t3[2] = (t2[3]-t2[2]) / (t1[3]+t1[2]);
        double a = t2[1]+mm(t3[0],t3[1])*t1[1];
        double b = t2[2]-mm(t3[1],t3[2])*t1[2];          
        
        /* y */
        t0[0] = u.getym(u.ym(jj));
        t0[1] = u.getym(jj);
        t0[2] = u[jj];
        t0[3] = u.getyp(jj);
        t0[4] = u.getyp(u.yp(jj));
        t1[0] = dy; t1[1] = dy; t1[2] = dy; t1[3] = dy; 
        if((mysigntol(t0[0],me)*mysigntol(t0[1],me))==-1)
        {
          t1[0] = dy * t0[1] / (t0[1]-t0[0]);
          t0[0] = 0.0f;
        }
        if((mysigntol(t0[3],me)*mysigntol(t0[4],me))==-1)
        {
          t1[3] = -dy * t0[3] / (t0[4]-t0[3]);
          t0[4] = 0.0f;
        }
        
        t2[0] = (t0[1]-t0[0]) / t1[0];
        t2[1] = (t0[2]-t0[1]) / t1[1];
        t2[2] = (t0[3]-t0[2]) / t1[2];
        t2[3] = (t0[4]-t0[3]) / t1[3];
        t3[0] = (t2[1]-t2[0]) / (t1[1]+t1[0]);
        t3[1] = (t2[2]-t2[1]) / (t1[2]+t1[1]);
        t3[2] = (t2[3]-t2[2]) / (t1[3]+t1[2]);
        double c = t2[1]+mm(t3[0],t3[1])*t1[1];
        double d = t2[2]-mm(t3[1],t3[2])*t1[2];          
        
        if(u[jj] > me)
          G[ii] = sqrt(mymax((static_cast<double>(a>0.0f))*a*a,(static_cast<double>(b<0.0f))*b*b) +
                       mymax((static_cast<double>(c>0.0f))*c*c,(static_cast<double>(d<0.0f))*d*d)) - 1.0f;
        else if(u[jj] < -me)
          G[ii] = sqrt(mymax((static_cast<double>(a<0.0f))*a*a,(static_cast<double>(b>0.0f))*b*b) +
                       mymax((static_cast<double>(c<0.0f))*c*c,(static_cast<double>(d>0.0f))*d*d)) - 1.0f;
        else
          G[ii] = 0.0f;
      } // if non-interface cell - compute values to update
    for(size_t ii=0; ii<idxc.size(); ++ii)
      if(!(df[ii])) // update all non-interface cells simultaneously
        u[idxc[ii]] = u[idxc[ii]]-dt*mysign(u[idxc[ii]])*G[ii];
  } // 2nd order iterations
  delete[] G;
}

bool Redist::diffSign(int idx)
{ // return true if sign of u[idx] differs from sign of any of its 4 neighbors
  const double ui = u0[idx];
  if(mysign(ui) != mysign(u0.getxp(idx)))
    return(true);
  if(mysign(ui) != mysign(u0.getyp(idx)))
    return(true);
  if(mysign(ui) != mysign(u0.getxm(idx)))
    return(true);
  if(mysign(ui) != mysign(u0.getym(idx)))
    return(true);
  return(false);
}
