#include "redist3.hpp"

void Redist3::redistance()
{
  if((u.getn() != u.getm()) || (u.getn() != u.getk()))
    cout << " Warning: m = " << u.getm() << ", n = " << u.getn() << ", and k = " << u.getk() << ". Code needs to be checked carefully for non-cubic grid! Trying to continue ..." << endl;

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

void Redist3::fastMarchingRedist()
{
  setInterfaceValues();
  thresholdAwayFromInterface(bnd);

  for(size_t ii=0;ii<bnd.size();++ii)
    updateAndAddNeighborsToHeap(bnd[ii]);

  // while heap is non-empty, fix its top value, and update the neighbors
  while(1)
  {
    struct helt h1 = h.popFromHeap();
    //printf(" h1.i = %d, hi.d = %f\n",h1.i,h1.d);
    if(h1.i == HEAPDONE) // if heap is empty, break from loop
      break; 
    if((state.get(h1.i) == false) && (fabs(h1.d) < thres)) // don't update if outside threshold cutoff
    {
      applyResult(h1);
      state.put(true,h1.i);
      updateAndAddNeighborsToHeap(h1.i);
    }
  }
}

void Redist3::updateAndAddNeighborsToHeap(const int idx)
{ // fix the value at u[idx], update its 6 neighbors if they are not already fixed, and make sure they're in the heap at the appropriate location
  int idx2arr[6];
  idx2arr[0] = u.xp(idx);
  idx2arr[1] = u.xm(idx); 
  idx2arr[2] = u.yp(idx);
  idx2arr[3] = u.ym(idx);
  idx2arr[4] = u.zp(idx);
  idx2arr[5] = u.zm(idx);
  for(int ii=0;ii<6;++ii)
  {
    if(state.get(idx2arr[ii]) == false) // if true, value is already fixed
    {
      double dtemp = estimateUpdate(idx2arr[ii]);
      if(fabs(dtemp) < fabs(u.get(idx2arr[ii])))
      {
        u.put(dtemp,idx2arr[ii]);
        h.addToHeap(idx2arr[ii],fabs(dtemp));
      }
    }
  }
}

void Redist3::updateAndAddNeighborsToHeapDO(const int idx)
{ // fix the value at u[idx], update its 4 neighbors if they are not already fixed, and make sure they're in the heap at the appropriate location
  int idx2arr[6];
  idx2arr[0] = u.xp(idx);
  idx2arr[1] = u.xm(idx); 
  idx2arr[2] = u.yp(idx);
  idx2arr[3] = u.ym(idx);
  idx2arr[4] = u.zp(idx);
  idx2arr[5] = u.zm(idx);
  for(int ii=0;ii<6;++ii)
  {
    if(state.get(idx2arr[ii]) == false) // if true, value is already fixed
    {
      struct helt htemp = performDO(idx2arr[ii],cpx.get(idx),cpy.get(idx),cpz.get(idx));
      if(fabs(htemp.d) < fabs(u.get(idx2arr[ii])))
      {
        applyResult(htemp);
        h.addToHeap(idx2arr[ii],fabs(htemp.d),htemp.aux);
      }
    }
  }
}

double Redist3::estimateUpdate(const int idx)
{
  double a = mymin(fabs(u.getxm(idx)),fabs(u.getxp(idx)));
  double b = mymin(fabs(u.getym(idx)),fabs(u.getyp(idx)));
  double c = mymin(fabs(u.getzm(idx)),fabs(u.getzp(idx)));
  double d;
  sort(a,b,c);
  assert(a<=b);
  assert(b<=c);

  if((b-a)>=dx)
    d = a+dx;
  else
  {
    d = 0.5f * (a+b+sqrt(2.0f*dx*dx-(a-b)*(a-b)));
    if(d > c)
      d = 1.0f/3.0f * (a+b+c+sqrt(3.0f*dx*dx+2.0f*(a*b+a*c+b*c-a*a-b*b-c*c)));
  }
  return(d*mysign(u.get(idx)));
}

void Redist3::sort(double &a1, double &a2, double &a3)
{
  /* Sort the values in a1, a2, and a3 so that
     we return with a1<=a2<=a3 using c as temp storage */
  double c;
  if(a2 < a1)
  {  c = a2; a2 = a1; a1 = c;}
  if(a3 < a1)
  { c = a3; a3 = a2; a2 = a1; a1 = c;}
  else if(a3 < a2)
  { c = a3; a3 = a2; a2 = c;}
}

void Redist3::setInterfaceValues()
{ 
  Array3D<int> sgn(u.getm(),u.getn(),u.getk());
  for(size_t ii=0; ii<sgn.getN(); ++ii)
    sgn.put(static_cast<int>(mysign(u.get(ii))),ii);
  for(size_t ii=0; ii<sgn.getN(); ++ii)
    if((abs(sgn.getxp(ii)-sgn.get(ii)) + abs(sgn.getxm(ii)-sgn.get(ii)) +
        abs(sgn.getyp(ii)-sgn.get(ii)) + abs(sgn.getym(ii)-sgn.get(ii)) +
        abs(sgn.getzp(ii)-sgn.get(ii)) + abs(sgn.getzm(ii)-sgn.get(ii)) )>0)
      bnd.push_back(ii);

  vector<double> dr(bnd.size());

  for(size_t ii=0; ii<bnd.size(); ++ii)
  {
    // compute norm(grad u) with centered differences
    double rx = (u.getxp(bnd[ii])-u.getxm(bnd[ii]))/dx/2.0f;
    double ry = (u.getyp(bnd[ii])-u.getym(bnd[ii]))/dy/2.0f;
    double rz = (u.getzp(bnd[ii])-u.getzm(bnd[ii]))/dz/2.0f;
    dr[ii] = sqrt(rx*rx+ry*ry+rz*rz);

    // compute norm(grad u) with one-sided differences
    rx = mymax(fabs(u.getxp(bnd[ii])-u.get(bnd[ii])),fabs(u.get(bnd[ii])-u.getxm(bnd[ii])))/dx;
    ry = mymax(fabs(u.getyp(bnd[ii])-u.get(bnd[ii])),fabs(u.get(bnd[ii])-u.getym(bnd[ii])))/dy;
    rz = mymax(fabs(u.getzp(bnd[ii])-u.get(bnd[ii])),fabs(u.get(bnd[ii])-u.getzm(bnd[ii])))/dz;
    const double dr2 = sqrt(rx*rx+ry*ry+rz*rz);

    // Accept one-sided difference is much different than centered difference
    if((dr[ii] < (0.5*dr2)) || (dr[ii] > (2.0*dr2)))
      dr[ii] = dr2;
  }
  for(size_t ii=0; ii<bnd.size();++ii)
    u.put(u.get(bnd[ii])/dr[ii],bnd[ii]);
}

void Redist3::setInterfaceValuesDO()
{ 
  Array3D<int> sgn(m,n,k);

  double meanabsbndval = 0.0f;

  for(size_t ii=0; ii<sgn.getN(); ++ii)
    sgn.put(static_cast<int>(mysigntol(u0.get(ii),1e-10)),ii);
  for(size_t ii=0; ii<sgn.getN(); ++ii)
    if((abs(sgn.getxp(ii)-sgn.get(ii)) + abs(sgn.getxm(ii)-sgn.get(ii)) +
        abs(sgn.getyp(ii)-sgn.get(ii)) + abs(sgn.getym(ii)-sgn.get(ii)) +
        abs(sgn.getzp(ii)-sgn.get(ii)) + abs(sgn.getzm(ii)-sgn.get(ii)))>0)
    {
      bnd.push_back(ii);
      meanabsbndval += fabs(u0.get(ii));
    }

  const double dr = mymax3(dx,dy,dz);
  meanabsbndval /= static_cast<double>(bnd.size());
  if(meanabsbndval > dr)
    for(size_t ii=0; ii<u0.getN(); ++ii)
      u0.put(u0.get(ii)/meanabsbndval*dr,ii);

  vector<struct helt> bndval;
  bndval.resize(bnd.size());
  for(size_t ii=0; ii<bnd.size(); ++ii)
  {
    bndval[ii] = performDO(bnd[ii]);
    cpx.put(bndval[ii].aux[0],bndval[ii].i);
    cpy.put(bndval[ii].aux[1],bndval[ii].i);
    cpz.put(bndval[ii].aux[2],bndval[ii].i);
  }

  for(size_t ii=0; ii<bnd.size();++ii)
  {
    applyResult(bndval[ii]);
    state.put(true,bndval[ii].i);
  }

  //printf("u[15597] = %f\n",u.get(15597));
  //printf(" 0.7\n");fflush(stdout);

}

void Redist3::thresholdAwayFromInterface(vector<int> &bndry)
{
  for(size_t ii=0;ii<state.getN();++ii)
    state.put(false,ii);
  for(vector<int>::iterator it=bndry.begin(); it != bndry.end(); ++it)
  {
    state.put(true,*it); // true indicates value is fixed (false otherwise)
  }
  for(size_t ii=0;ii<state.getN();++ii)
    if(!state.get(ii))
      u.put(mysign(u.get(ii))*thres,ii);
}

Redist3::Redist3(const Array3D<double> &_u, const int _width, const int _flag) :
  width(_width),
  flag(_flag),
  m(_u.getm()),
  n(_u.getn()),
  k(_u.getk()),
  N(m*n*k),
  dx(1.0f/static_cast<double>(n)),
  dy(1.0f/static_cast<double>(m)),
  dz(1.0f/static_cast<double>(k)),
  thres(static_cast<double>(width+1)*mymax3(dx,dy,dz)),
  cpflag((_flag==2) || (_flag==3)),
  h(10*width,cpflag),
  state(m,n,k),
  u0(_u,_flag),
  u(_u),
  cpx(m,n,k),
  cpy(m,n,k),
  cpz(m,n,k)
{ }

void Redist3::directionalOptimization()
{
  //printf("\n0\n"); fflush(stdout);
  setInterfaceValuesDO(); 
  //printf("1\n"); fflush(stdout);
  thresholdAwayFromInterface(bnd);
  //printf("2\n"); fflush(stdout);

  for(size_t ii=0;ii<bnd.size();++ii)
    updateAndAddNeighborsToHeapDO(bnd[ii]);
  //printf("3\n"); fflush(stdout);

  while(1)
  {
    struct helt h1 = h.popFromHeap();
    if(h1.i == HEAPDONE) // if heap is empty, break from loop
      break; 
    if((state.get(h1.i)==false) && (fabs(h1.d) < thres)) // don't update if outside threshold cutoff, or already fixed
    {
      applyResult(h1);
      cpx.put(h1.aux[0],h1.i);
      cpy.put(h1.aux[1],h1.i);
      cpz.put(h1.aux[2],h1.i);
      state.put(true,h1.i);
      updateAndAddNeighborsToHeapDO(h1.i);
    }
  }
  //printf("4\n"); fflush(stdout);
}

void Redist3::applyResult(const struct helt &h)
{
  u.put(mysign(u.get(h.i))*h.d,h.i);
  if(flag > 1)
  {
    cpx.put(h.aux[0],h.i);
    cpy.put(h.aux[1],h.i);
    cpz.put(h.aux[2],h.i);
  }
}

struct helt Redist3::performDO(const int idx)
{ // make a guess for interface location and call performDO(const int ind, const double *cpguess)
  double grad[3], cpguess[3];
  grad[0] = pd(u.getxp(idx),u.getxm(idx)) / (2.0f * dx);
  grad[1] = pd(u.getyp(idx),u.getym(idx)) / (2.0f * dy);
  grad[2] = pd(u.getzp(idx),u.getzm(idx)) / (2.0f * dz);
  //printf("  0.5.0\n"); fflush(stdout);
  lineSearch(idx,grad,cpguess);

  //printf(" idx = %d, grad=(%f,%f,%f), cpguess = (%f,%f,%f)\n",idx,grad[0],grad[1],grad[2],cpguess[0],cpguess[1],cpguess[2]);
  //printf("    u(idx) = %f, u(cpguess) = %f.\n",u0.get(idx),u0.interpolate(cpguess[0],cpguess[1],cpguess[2]));

  //printf("  0.5.1\n"); fflush(stdout);
  return(performDO(idx,cpguess[0],cpguess[1],cpguess[2]));
}

struct helt Redist3::performDO(const int idx, const double cpxguess, const double cpyguess, const double cpzguess) 
{
  //bool VERBOSITY = false;
  //if(idx == 10634)
  //  VERBOSITY = true;
    
  // optimize over directions via line search
  const double x0[3] = {u.getX(idx), u.getY(idx), u.getZ(idx)};
  const int MAXLOOPS = 5;
  const int MAXSEARCHES = 5;
  const double TOL = mymax( (0.01 / pow(mymax3(m,n,k),flag+1)), 5.0e-16) ;
  const int NDIR = 5;

  double xt[3], xx[5][3], rvec[3], v1[3], v2[3];
  double dd[5], dt;
  xx[0][0] = cpxguess;
  xx[0][1] = cpyguess;
  xx[0][2] = cpzguess;

  dd[0] = dist3(x0,xx[0]);

  //if(VERBOSITY)
  //  printf(" xx[0] = (%f,%f,%f), dd[0] = %f\n",xx[0][0],xx[0][1],xx[0][2],dd[0]);

  // set angle increment
  double delta;
  if(dd[0] > mymax3(dx,dy,dz))
    delta = asin(mymax3(dx,dy,dz)/dd[0]);
  else // interface is nearby, look over a wide range
    delta = PI2;

  //printf("  0.5.2\n"); fflush(stdout);

  double nn = findNborDirections(x0,xx,delta,rvec,v1,v2);

  //if(VERBOSITY)
  //  printf(" after findNborDirections: xx[0] = (%f,%f,%f), dd[0] = %f\n",xx[0][0],xx[0][1],xx[0][2],dd[0]);

  //printf("  0.5.3\n"); fflush(stdout);
  
  // bracket the interface in the search directions
  for(int ii=0; ii<NDIR; ++ii)
    dd[ii] = search1D(idx,xx[ii]);

  //if(VERBOSITY)
  //  printf(" after search1D: xx[0] = (%f,%f,%f), dd[0] = %f\n",xx[0][0],xx[0][1],xx[0][2],dd[0]);

  int minidx = minabs(dd,NDIR);

  //if(VERBOSITY)
  //  printf(" dd = [%f %f %f %f %f], minidx = %d\n",dd[0],dd[1],dd[2],dd[3],dd[4],minidx);

  //printf("  0.5.4\n"); fflush(stdout);

  for(int count0=0; count0<MAXLOOPS; ++count0)
  { // (1) center on minimum; (2) perform Newton step
    int count = 0;
    while( (minidx != 0) && (count < MAXSEARCHES))
    {
      if(delta < PI4-1e-4)
        delta *= 2.0f;
      if(minidx == 1)
      { xx[0][0] = xx[1][0]; xx[0][1] = xx[1][1]; xx[0][2] = xx[1][2]; dd[0] = dd[1];}
      else if(minidx == 2) 
      { xx[0][0] = xx[2][0]; xx[0][1] = xx[2][1]; xx[0][2] = xx[2][2]; dd[0] = dd[2];}
      else if(minidx == 3) 
      { xx[0][0] = xx[3][0]; xx[0][1] = xx[3][1]; xx[0][2] = xx[3][2]; dd[0] = dd[3];}
      else /* (minidx == 4) */ 
      { xx[0][0] = xx[4][0]; xx[0][1] = xx[4][1]; xx[0][2] = xx[4][2]; dd[0] = dd[4];}
      /* update directions around new center location */
      nn = findNborDirections(x0,xx,delta,rvec,v1,v2);
      for(int ii=1; ii<NDIR; ++ii)
        dd[ii] = search1D(idx,xx[ii]);
      minidx = minabs(dd,NDIR);

      //if(VERBOSITY)
      //  printf(" count = %d, dd = [%f %f %f %f %f], minidx = %d\n",count,dd[0],dd[1],dd[2],dd[3],dd[4],minidx);
      count++;
    }
    // perform Newton step only if dd values are successfully computed
    if(mymax(mymax(dd[1],dd[2]),mymax(dd[3],dd[4])) < 0.9f)
    {
      double us = (dd[1]-dd[3])/2.0;
      double uss= (dd[1]-2.0*dd[0]+dd[3]);
      double ut = (dd[2]-dd[4])/2.0;
      double utt= (dd[2]-2.0*dd[0]+dd[4]);
      if(fabs(uss) < 1e-16)
        uss = 1e-16;
      if(fabs(utt) < 1e-16)
        utt = 1e-16;
      double deltas = -us/uss*delta;
      double deltat = -ut/utt*delta;
      xt[0] = dinrange2(x0[0] + nn/2.0*(rvec[0]*cos(deltas)+v1[0]*sin(deltas)+rvec[0]*cos(deltat)+v2[0]*sin(deltat)));
      xt[1] = dinrange2(x0[1] + nn/2.0*(rvec[1]*cos(deltas)+v1[1]*sin(deltas)+rvec[1]*cos(deltat)+v2[1]*sin(deltat)));
      xt[2] = dinrange2(x0[2] + nn/2.0*(rvec[2]*cos(deltas)+v1[2]*sin(deltas)+rvec[2]*cos(deltat)+v2[2]*sin(deltat)));
      dt = search1D(idx,xt);

      if(fabs(dt-dd[0]) < TOL)
        count0 = MAXLOOPS; // no further improvement available
      else if(dt<dd[0])
      {
        xx[0][0] = xt[0]; xx[0][1] = xt[1]; xx[0][2] = xt[2]; dd[0] = dt;
      }
    }
    if(count0 < MAXLOOPS)
    {
      delta /= 2.0f;
      nn = findNborDirections(x0,xx,delta,rvec,v1,v2);
      for(int ii=0;ii<NDIR;ii++)
        dd[ii] = search1D(idx,xx[ii]);
      //if(VERBOSITY)
      //  printf(" count0 = %d, dd = [%f %f %f %f %f], minidx = %d\n",count0,dd[0],dd[1],dd[2],dd[3],dd[4],minidx);

    }
  }
  struct helt lval;

  //printf("  0.5.5\n"); fflush(stdout);
  double val1 = u0.interpolate(xx[0][0],xx[0][1],xx[0][2]);
  lval.i = idx;
  lval.d = dd[0] + fabs(val1);
  //if(VERBOSITY)
  //  printf("   dd[0] = %f, val1 = %f\n",dd[0],val1);
  lval.aux[0] = xx[0][0];
  lval.aux[1] = xx[0][1];
  lval.aux[2] = xx[0][2];
  //printf("  0.5.6\n"); fflush(stdout);

  return(lval);
}

void Redist3::lineSearch(const int idx, double (&grad)[3], double (&cpguess)[3])
{ // takes the given index at estimate of gradient direction (unnormalized)
  // and performs a line search along this direction for the interface.
  // Interface is expected to lie within one grid cell of idx

  normalize(grad);
  // get the correct direction to search in; 
  if(u.get(idx) < 0.0f) // then we go uphill, keep gradient sign
  {
    grad[0] *= SQRT3/static_cast<double>(n);
    grad[1] *= SQRT3/static_cast<double>(m); 
    grad[2] *= SQRT3/static_cast<double>(k);
  }
  else // go downhill
  {
    grad[0] *= -SQRT3/static_cast<double>(n);
    grad[1] *= -SQRT3/static_cast<double>(m);
    grad[2] *= -SQRT3/static_cast<double>(k);
  }

  if(findOppSign(idx,grad))
  {
    cpguess[0] = grad[0]; cpguess[1] = grad[1]; cpguess[2] = grad[2];
    bisect(idx,cpguess);
  }
  else // need to figure out where to look
  {
    const int sgn = mysign(u.get(idx));
    int nborind[18];
    nborind[0] = u.xp(idx);
    nborind[1] = u.xm(idx);
    nborind[2] = u.yp(idx);
    nborind[3] = u.ym(idx);
    nborind[4] = u.zp(idx);
    nborind[5] = u.zm(idx);
    nborind[6] = u.yp(u.xp(idx));
    nborind[7] = u.yp(u.xm(idx));
    nborind[8] = u.ym(u.xp(idx));
    nborind[9] = u.ym(u.xm(idx));
    nborind[10] = u.zp(u.xp(idx));
    nborind[11] = u.zp(u.xm(idx));
    nborind[12] = u.zm(u.xp(idx));
    nborind[13] = u.zm(u.xm(idx));
    nborind[14] = u.zp(u.yp(idx));
    nborind[15] = u.zp(u.ym(idx));
    nborind[16] = u.zm(u.yp(idx));
    nborind[17] = u.zm(u.ym(idx));

    for(int ii=0;ii<18;++ii)
    {
      if(mysign(u.get(nborind[ii])) != sgn)
      {
        cpguess[0] = u.getX(nborind[ii]);
        cpguess[1] = u.getY(nborind[ii]);
        cpguess[2] = u.getZ(nborind[ii]);
        bisect(idx,cpguess);
        return;
      }
    }
  }
}

bool Redist3::findOppSign(const int idx, double (&guess)[3])
{ // searches in the direction initially given by "guess" to find x s.t. u(x) has opposite sign
  // from from u(idx)
  // input suggests to look at u[y+guess], where y is spatial coordinates of idx. 
  // returns location x in guess.
  const int sgn = mysign(u.get(idx));

  double val1 = u0.interpolate(u.getX(idx)+guess[0],u.getY(idx)+guess[1],u.getZ(idx)+guess[2]);
  if(sgn == mysign(val1)) // need to look harder!
  {
    double val2 = u0.interpolate(u.getX(idx)+2.0f*guess[0],u.getY(idx)+2.0f*guess[1],u.getZ(idx)+2.0f*guess[2]);
    if(sgn == mysign(val2)) // still need to do more work
      return(false);
    else
    {
      guess[0] = u.getX(idx)+2.0f*guess[0];
      guess[1] = u.getY(idx)+2.0f*guess[1];
      guess[2] = u.getZ(idx)+2.0f*guess[2];
      return(true);
    }
  }
  else
  {
    guess[0] = u.getX(idx)+guess[0];
    guess[1] = u.getY(idx)+guess[1];
    guess[2] = u.getZ(idx)+guess[2];
    return(true);
  }
}

bool Redist3::bracket(const int idx, double (&cpguess)[3], double (&xm)[4], double (&xp)[4])
{ // Look about location cpguess along line from idx and find (small) positive and negative values as
  // input to bisection.
  
  // Returns false if it's much further than initially guessed to interface in this direction (e.g. more
  // than one grid cell further), and true otherwise.
  // Assumed that xm, xp, both have room for two doubles

  //bool VERBOSITY = false;
  //if(idx == 10634)
  //  VERBOSITY = true;

  const double x0[3] = {u.getX(idx), u.getY(idx), u.getZ(idx)};
  const double uv0 = u.get(idx);
  if(mysign(uv0) == 0) // already on interface
  {
    xm[0] = x0[0]; xm[1] = x0[1]; xm[2] = x0[2]; xm[3] = uv0;
    xp[0] = x0[0]; xp[1] = x0[1]; xp[2] = x0[2]; xp[3] = uv0;
    return(true);
  }

  double dir[3];
  dir[0] = pd(cpguess[0],x0[0]);
  dir[1] = pd(cpguess[1],x0[1]);
  dir[2] = pd(cpguess[2],x0[2]);
  normalize(dir);

  double ug = u0.interpolate(cpguess[0],cpguess[1],cpguess[2]);
  if(mysign(uv0 * ug) == 0) // sign(ug) == 0; 
  {
    xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = cpguess[2]; xm[3] = ug;
    xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = cpguess[2]; xp[3] = ug;
    return(true);
  }
  if(mysign(uv0 * ug) == 1) // same side of interface
  {
    double dr = mymax3(dx,dy,dz);
    double ug2 = u0.interpolate(cpguess[0]+dir[0]*dr,cpguess[1]+dir[1]*dr,cpguess[2]+dir[2]*dr);
    if(mysign(uv0 * ug2) == 1) // still on same side of interface
      return(false);
    else
    {
      if(mysign(ug) == 0)
      {
        if(mysign(ug2) == 1)
        {
          xp[0] = cpguess[0]+dir[0]*dr; xp[1] = cpguess[1]+dir[1]*dr; xp[2] = cpguess[2]+dir[2]*dr; xp[3] = ug2;
          xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = cpguess[2]; xm[3] = ug;
        }
        else
        {
          xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = cpguess[2]; xp[3] = ug;
          xm[0] = cpguess[0]+dir[0]*dr; xm[1] = cpguess[1]+dir[1]*dr; xm[2] = cpguess[2]+dir[2]*dr; xm[3] = ug2;
        }
      }
      else if(mysign(ug) == 1)
      {
        xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = cpguess[2]; xp[3] = ug;
        xm[0] = cpguess[0]+dir[0]*dr; xm[1] = cpguess[1]+dir[1]*dr; xm[2] = cpguess[2]+dir[2]*dr; xm[3] = ug2;
      }
      else // mysign(ug) == -1
      {
        xp[0] = cpguess[0]+dir[0]*dr; xp[1] = cpguess[1]+dir[1]*dr; xp[2] = cpguess[2]+dir[2]*dr; xp[3] = ug2;
        xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = cpguess[2]; xm[3] = ug;
      }
      return(true);
    }
  }
  else // x0 and cpguess are on opposite sides of interface, work backwards
  {
    double dr = mymax3(dx,dy,dz);
    double xx,yy,zz;
    xx = cpguess[0]-dir[0]*dr;
    yy = cpguess[1]-dir[1]*dr;
    zz = cpguess[2]-dir[2]*dr;
    double ug2 = u0.interpolate(xx,yy,zz);
    int tries = 0;
    while((mysign(ug*ug2) == 1) && (tries < 10))
    {
      xx -= dir[0]*dr;
      yy -= dir[1]*dr;
      zz -= dir[2]*dr;
      ug2 = u0.interpolate(xx,yy,zz);
      tries++;
    }
    if(tries == 10)
      return(false);

    if(mysign(ug) == 1)
    {
      xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = cpguess[2]; xp[3] = ug;
      xm[0] = xx; xm[1] = yy; xm[2] = zz; xm[3] = ug2;
    }
    else if(mysign(ug) == -1)
    {
      xp[0] = xx; xp[1] = yy; xp[2] = zz; xp[3] = ug2;
      xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = cpguess[2]; xm[3] = ug;
    }
    else // mysign(ug) == 0
    {
      if(mysign(ug2) == 1)
      {
        xp[0] = xx; xp[1] = yy; xp[2] = zz; xp[3] = ug2;
        xm[0] = cpguess[0]; xm[1] = cpguess[1]; xm[2] = cpguess[2]; xm[3] = ug;
      }
      else
      {
        xp[0] = cpguess[0]; xp[1] = cpguess[1]; xp[2] = cpguess[2]; xp[3] = ug;
        xm[0] = xx; xm[1] = yy; xm[2] = zz; xm[3] = ug2;
      }
    }
    return(true);
  }
}

void Redist3::bisect(const int idx, double (&guess)[3])
{ // bisection without Newton
  const int NTRIES = 10;
  // Assumes that u[idx] and u[guess] have different signs.
  // Tries to find the zero between u[idx] and u[guess] by bisection.
  
  double xm[3], xp[3], xt[3];
  double um, up;
  assert(mysign(u.get(idx)*u0.interpolate(guess[0],guess[1],guess[2])) != 1);

  if(mysign(u.get(idx)) == 0)
  {
    guess[0] = u.getX(idx);
    guess[1] = u.getY(idx);
    guess[2] = u.getZ(idx);
    return;
  }
  else if(mysign(u.get(idx)) == 1)
  {
    xp[0] = u.getX(idx); xp[1] = u.getY(idx); xp[2] = u.getZ(idx); up = u.get(idx);
    xm[0] = guess[0];    xm[1] = guess[1];    xm[2] = guess[2];    um = u0.interpolate(xm[0],xm[1],xm[2]);
  }
  else
  {
    xm[0] = u.getX(idx); xm[1] = u.getY(idx); xm[2] = u.getZ(idx); um = u.get(idx);
    xp[0] = guess[0];    xp[1] = guess[1];    xp[2] = guess[2];    up = u0.interpolate(xp[0],xp[1],xp[2]);

    if(up <= 0.0f)
      printf(" up = %e, um = %e\n",up,um);
  }
  if(up==um)
  {
    guess[0] = xp[0]; guess[1] = xp[1]; guess[2] = xp[2];
    return;
  }
  if(std::isnan(up/(up-um)))
  {
    printf(" up = %.3e, um = %.3e. Setting um = -1.0f.\n",up,um);
    um = -1.0f;
  }
  ccomb(xm,xp,xt,up/(up-um));
  assert(up>=0.0f); 
  assert(um<=0.0f);
  for(int ii=0; ii<NTRIES; ++ii)
  {
    double ut = u0.interpolate(xt[0],xt[1],xt[2]);
    if(mysign(ut) == 1)
    {
      xp[0] = xt[0]; xp[1] = xt[1]; xp[2] = xt[2]; up = ut;
    }
    else
    {
      xm[0] = xt[0]; xm[1] = xt[1]; xm[2] = xt[2]; um = ut;
    }
    if(std::isnan(up/(up-um)))
    {
      printf(" up = %.4e, um = %.4e, setting up = 1.0f\n",up,um);
      up = 1.0f;
    }
    ccomb(xm,xp,xt,up/(up-um));
    assert(up>=0.0f); 
    assert(um<=0.0f);
  }
  guess[0] = xt[0];
  guess[1] = xt[1];
  guess[2] = xt[2];
  return;
}

double Redist3::bisect(double (&result)[3], double (&xm)[4], double (&xp)[4])
{ // bisection with linear approximation
  const int NTRIES = 10;
  const double UTOL = 1e-14;
  const double XTOL = 5e-16;
  // Assumes that u[idx] and u[guess] have different signs.
  // Tries to find the zero between u[idx] and u[guess] by bisection.
  double up = xp[3]; //u0.interpolate(xp[0],xp[1],xp[2]);
  double um = xm[3]; //u0.interpolate(xm[0],xm[1],xm[2]);
  double xt[3];
  assert(up>=0.0f); 
  assert(um<=0.0f);
  if(up==um)
  {
    result[0] = xp[0]; result[1] = xp[1]; result[2] = xp[2];
    return(up);
  }
  ccomb(xm,xp,xt,up/(up-um));
  double ut = u0.interpolate(xt[0],xt[1],xt[2]); 
  for(int ii=0; ii<NTRIES; ++ii)
  {
    if(mysign(ut) == 0)
    {
      result[0] = xt[0]; result[1] = xt[1]; result[2] = xt[2];
      return(ut);
    }
    else if(mysign(ut) == 1)
    {
      xp[0] = xt[0]; xp[1] = xt[1]; xp[2] = xt[2]; up = ut;
    }
    else
    {
      xm[0] = xt[0]; xm[1] = xt[1]; xm[2] = xt[2]; um = ut;
    }
    ccomb(xm,xp,xt,up/(up-um));
    assert(up>=0.0f); 
    assert(um<=0.0f);
    if(dist3(xp,xm) < XTOL)
      break;
    ut = u0.interpolate(xt[0],xt[1],xt[2]); 
    if(mysigntol(ut,UTOL) == 0)
      break;
  }
  result[0] = xt[0]; result[1] = xt[1]; result[2] = xt[2];
  return(ut);
}

double Redist3::findNborDirections(const double (&x0)[3], double (&xx)[5][3], const double delta, double (&rvec)[3], double (&v1)[3], double (&v2)[3])
{
  double s,c,nn;
  s = sin(delta);
  c = cos(delta);

  rvec[0] = pd(xx[0][0],x0[0]);
  rvec[1] = pd(xx[0][1],x0[1]);
  rvec[2] = pd(xx[0][2],x0[2]);
  nn = mymax(sqrt(rvec[0]*rvec[0]+rvec[1]*rvec[1]+rvec[2]*rvec[2]), 1e-16);
  normalize(rvec);
  orthoVecs(rvec,v1,v2);

  xx[1][0] = dinrange2(x0[0] + nn*(rvec[0]*c+v1[0]*s));
  xx[1][1] = dinrange2(x0[1] + nn*(rvec[1]*c+v1[1]*s));
  xx[1][2] = dinrange2(x0[2] + nn*(rvec[2]*c+v1[2]*s));
  xx[2][0] = dinrange2(x0[0] + nn*(rvec[0]*c+v2[0]*s));
  xx[2][1] = dinrange2(x0[1] + nn*(rvec[1]*c+v2[1]*s));
  xx[2][2] = dinrange2(x0[2] + nn*(rvec[2]*c+v2[2]*s));
  xx[3][0] = dinrange2(x0[0] + nn*(rvec[0]*c-v1[0]*s));
  xx[3][1] = dinrange2(x0[1] + nn*(rvec[1]*c-v1[1]*s));
  xx[3][2] = dinrange2(x0[2] + nn*(rvec[2]*c-v1[2]*s));
  xx[4][0] = dinrange2(x0[0] + nn*(rvec[0]*c-v2[0]*s));
  xx[4][1] = dinrange2(x0[1] + nn*(rvec[1]*c-v2[1]*s));
  xx[4][2] = dinrange2(x0[2] + nn*(rvec[2]*c-v2[2]*s));
  return(nn);
}

double Redist3::search1D(const int idx, double (&x)[3])
{
  //bool VERBOSITY = false;
  //if(idx == 10634)
  //  VERBOSITY = true;

  // bracket the interface in the search directions
  double  x0[3], xm[4], xp[4]; // xm and xp also return values from u0
  bool gl = bracket(idx,x,xm,xp);

  //if(VERBOSITY)
  //  printf("gl = %d\n",static_cast<int>(gl));

  double d;

  x0[0] = u.getX(idx); x0[1] = u.getY(idx); x0[2] = u.getZ(idx);
  // once bracketed successfully, bisect
  if(gl)
  {
    double val1 = bisect(x,xm,xp);
    d = dist3(x0,x);
    d += val1 * mysign(u.get(idx));
  }
  else
    d = 1.0f;
  return(d);
}

void Redist3::dump_u(double *v)
{ // assumes sufficient memory is allocated into v
  for(int ii=0; ii<N; ++ii)
    v[ii] = u.get(ii);
}

void Redist3::dump_cp(double *cpx_d, double *cpy_d, double *cpz_d)
{ // assume sufficient memory is allocated into cpx, cpy, cpz
  for(int ii=0; ii<N; ++ii)
    cpx_d[ii] = cpx.get(ii);
  for(int ii=0; ii<N; ++ii)
    cpy_d[ii] = cpy.get(ii);
  for(int ii=0; ii<N; ++ii)
    cpz_d[ii] = cpz.get(ii);
}

void Redist3::secondOrderIterations()
{
  // currently, do nothing.
/*
  const double me = 1.0f/ 10000.0f / static_cast<double>(n);
  const double thres2 = thres - dx;
  const double dt = 1.0f / 5.0f / static_cast<double>(n);
  // Step 1: Figure out which indices to do the iterations on
  vector<int> idxc;
  for(int ii=0; ii<N; ++ii)
    if(fabs(u.get(ii)) < thres2)
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
        const double ujj = u.get(jj);

        double t0[5],t1[4],t2[4],t3[4];
        // x
        t0[0] = u.getxm(u.xm(jj));
        t0[1] = u.getxm(jj);
        t0[2] = ujj;
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
        
        // y
        t0[0] = u.getym(u.ym(jj));
        t0[1] = u.getym(jj);
        t0[2] = ujj;
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
        
        if(ujj > me)
          G[ii] = sqrt(mymax((static_cast<double>(a>0.0f))*a*a,(static_cast<double>(b<0.0f))*b*b) +
                       mymax((static_cast<double>(c>0.0f))*c*c,(static_cast<double>(d<0.0f))*d*d)) - 1.0f;
        else if(ujj < -me)
          G[ii] = sqrt(mymax((static_cast<double>(a<0.0f))*a*a,(static_cast<double>(b>0.0f))*b*b) +
                       mymax((static_cast<double>(c<0.0f))*c*c,(static_cast<double>(d>0.0f))*d*d)) - 1.0f;
        else
          G[ii] = 0.0f;
      } // if non-interface cell - compute values to update
    for(size_t ii=0; ii<idxc.size(); ++ii)
      if(!(df[ii])) // update all non-interface cells simultaneously
        u.put(u.get(idxc[ii])-dt*mysign(u.get(idxc[ii]))*G[ii],idxc[ii]);
  } // 2nd order iterations
  delete[] G;
*/
}

bool Redist3::diffSign(int idx)
{ // return true if sign of u[idx] differs from sign of any of its 6 first neighbors
  const double ui = u.get(idx);
  if(mysign(ui) != mysign(u.getxp(idx)))
    return(true);
  if(mysign(ui) != mysign(u.getyp(idx)))
    return(true);
  if(mysign(ui) != mysign(u.getzp(idx)))
    return(true);
  if(mysign(ui) != mysign(u.getxm(idx)))
    return(true);
  if(mysign(ui) != mysign(u.getym(idx)))
    return(true);
  if(mysign(ui) != mysign(u.getzm(idx)))
    return(true);
  return(false);
}
