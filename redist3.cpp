#include "redist3.hpp"
#include "toolbox3d.hpp"
#include <iostream>
#include <cstring>
#include <cassert>

namespace {
  void sort(double &a1, double &a2, double &a3)
  {
    /* Sort the values in a1, a2, and a3 so that
       we return with a1<=a2<=a3 */
    if(a2 < a1)
      std::swap(a1,a2);
    if(a3 < a2) {
      std::swap(a2,a3);
      if(a2 < a1)
        std::swap(a1,a2);
    }
  }
}

using BracketPair3 = std::pair<Point3, double>;

/// Constructor with full arguments
/// \param[in] _i   : Index value
/// \param[in] _aux : Auxilary double array values
Aux3::Aux3(idx_t const _i, Point3 const &_aux) :
  i(_i),
  aux(_aux)
{ }

/// Constructor from index only
/// \param[in] _i   : Index value
Aux3::Aux3(idx_t const _i) :
  i(_i),
  aux()
{ } 

/// Constructor
/// \param[in] _u     : Input level set function
/// \param[in] _width : Thresholding width (pixels)
/// \param[in] _flag  : Interpolation order flag
Redist3::Redist3(const Array3D<double> &_u, const idx_t _width, const int _flag) :
  width(_width),
  flag(_flag),
  m(_u.getm()),
  n(_u.getn()),
  k(_u.getk()),
  N(m*n*k),
  dx(1./static_cast<double>(n)),
  dy(1./static_cast<double>(m)),
  dz(1./static_cast<double>(k)),
  thres(static_cast<double>(width+1)*mymax3(dx,dy,dz)),
  cpflag((_flag==2) || (_flag==3)),
  h(10*width,cpflag),
  state(m,n,k),
  u0(_u,_flag),
  u(_u),
  cpx(m,n,k),
  cpy(m,n,k),
  cpz(m,n,k),
  bWarn(false)
{ }

/// Perform the redistancing
void Redist3::redistance()
{
  if((u.getn() != u.getm()) || (u.getn() != u.getk()))
    std::cout << " Warning: m = " << u.getm() << ", n = " << u.getn() << ", and k = " << u.getk() << ". Code needs to be checked carefully for non-cubic grid! Trying to continue ..." << std::endl;
  
  if((flag == 0) || (flag == 1)) {
    fastMarchingRedist();
    if(flag == 0)
      secondOrderIterations();
  } else if((flag == 2) || (flag == 3))
    directionalOptimization();
  else
    std::cout << "Do not recognize the flag. Cowardly refusing to do anything." << std::endl;
}

/// Internal function for fast marching redistancing
void Redist3::fastMarchingRedist()
{
  setInterfaceValues();
  thresholdAwayFromInterface();

  for (idx_t const ix : bnd)
    updateAndAddNeighborsToHeap(ix);

  // while heap is non-empty, fix its top value, and update the neighbors
  Helt3 helt;
  while(h.popFromHeap(helt)) {
    if((state.get(helt.second.i) == false) && (std::abs(helt.first) < thres)) { // don't update if outside threshold cutoff
      applyResult(helt);
      state.put(true, helt.second.i);
      updateAndAddNeighborsToHeap(helt.second.i);
    }
  }
}

/// Fix the value at u[idx], update its 4 neighbors if they are not already fixed, and make sure they're in the heap at the appropriate location
/// \param[in] idx : Location to work from
void Redist3::updateAndAddNeighborsToHeap(idx_t const idx)
{
  std::array<idx_t, 6> idx2arr = {u.xp(idx), u.xm(idx), u.yp(idx), u.ym(idx), u.zp(idx), u.zm(idx)};
  for(idx_t const ind : idx2arr) {
    if(!state.get(ind)) { // if true, value is already fixed
      double const dtemp = estimateUpdate(ind);
      if(std::abs(dtemp) < std::abs(u[ind])) {
        u[ind] = dtemp;
        h.addToHeap(Aux3(ind), std::abs(dtemp));
      }
    }
  }
}

/// Fix the value at u[idx], update its 6 neighbors if they are not already fixed, and make sure they're in the heap at the appropriate location -- for directional optimization
/// \param[in] idx : Location to work from
void Redist3::updateAndAddNeighborsToHeapDO(idx_t const idx)
{
  std::array<idx_t, 6> const idx2arr = {u.xp(idx), u.xm(idx), u.yp(idx), u.ym(idx), u.zp(idx), u.zm(idx)};
  for (idx_t const ind : idx2arr) {
    if(!state.get(ind)) { // if true, value is already fixed
      Helt3 const htemp = performDO(ind, Point3({cpx[idx],cpy[idx],cpz[idx]}));
      if(std::abs(htemp.first) < std::abs(u[ind])) {
        applyResult(htemp);
        h.addToHeap(htemp);
      }
    }
  }
}

/// Estimate update value at idx
/// \param[in] idx : location to estimate at
/// \return          Estimated distance to interface
double Redist3::estimateUpdate(idx_t const idx)
{
  double a = std::min(std::abs(u.getxm(idx)), std::abs(u.getxp(idx)));
  double b = std::min(std::abs(u.getym(idx)), std::abs(u.getyp(idx)));
  double c = std::min(std::abs(u.getzm(idx)), std::abs(u.getzp(idx)));
  ::sort(a,b,c);
  assert((a<=b) && (b<=c));
  double d;
  if((b-a)>=dx)
    d = a+dx;
  else {
    d = 0.5f * (a+b+sqrt(2.*dx*dx-(a-b)*(a-b)));
    if(d > c)
      d = 1./3. * (a+b+c+sqrt(3.*dx*dx+2.*(a*b+a*c+b*c-a*a-b*b-c*c)));
  }
  return(d*mysign(u[idx]));
}

/// Set values of output signed distance function at the interface
void Redist3::setInterfaceValues()
{
  Array3D<int> sgn(u.getm(),u.getn(),u.getk());
  std::vector<int> &vSgn = sgn.returnData();
  std::vector<double> const &vu = u.returnData();
  assert(vSgn.size() == vu.size());
  std::transform(vu.begin(), vu.end(), vSgn.begin(), [](double const d)->int{ return static_cast<int>(mysign(d)); });
  for(idx_t ii=0; ii<sgn.getN(); ++ii)
    if ((sgn.getxp(ii)-sgn[ii]) || (sgn.getxm(ii)-sgn[ii]) ||
        (sgn.getyp(ii)-sgn[ii]) || (sgn.getym(ii)-sgn[ii]) ||
        (sgn.getzp(ii)-sgn[ii]) || (sgn.getzm(ii)-sgn[ii]))
      bnd.push_back(ii);

  std::vector<double> dr(bnd.size());
  
  std::transform(bnd.begin(), bnd.end(), dr.begin(), [&](idx_t const ii)->double {
      // compute norm(grad u) with centered differences
      double rx = (u.getxp(ii)-u.getxm(ii))/dx/2.;
      double ry = (u.getyp(ii)-u.getym(ii))/dy/2.;
      double rz = (u.getzp(ii)-u.getzm(ii))/dz/2.;
      double dr = sqrt(rx*rx+ry*ry+rz*rz);

      // compute norm(grad u) with one-sided differences
      rx = std::max(fabs(u.getxp(ii)-u.get(ii)),fabs(u.get(ii)-u.getxm(ii)))/dx;
      ry = std::max(fabs(u.getyp(ii)-u.get(ii)),fabs(u.get(ii)-u.getym(ii)))/dy;
      rz = std::max(fabs(u.getzp(ii)-u.get(ii)),fabs(u.get(ii)-u.getzm(ii)))/dz;
      double const dr2 = sqrt(rx*rx+ry*ry+rz*rz);

      // Accept one-sided difference is much different than centered difference
      if((dr < (0.5*dr2)) || (dr > (2.0*dr2)))
        dr = dr2;
      return dr;
    });
  for(size_t ii=0; ii<bnd.size();++ii)
    u[bnd[ii]] /= dr[ii];
}

/// Set the interface values for directional optimization (must have correct order of accuracy)
void Redist3::setInterfaceValuesDO()
{ 
  double meanabsbndval = 0.;
  for(idx_t ii=0; ii<N; ++ii) 
    if(diffSign(ii)) {
      bnd.push_back(ii);
      meanabsbndval += std::abs(u0[ii]);
    }

  meanabsbndval /= static_cast<double>(bnd.size());
  // rescale u0 based on meanabsbndval
  u0 *= 0.5 * mymax3(dx,dy,dz) / meanabsbndval; // scale u0 (rationale: mean distance to interface at a grid cell should be ~ 0.5 dx)
  
  for (idx_t ind : bnd) {
    Helt3 const bndval = performDOSurf(ind);
    applyResult(bndval);
    state.put(true, ind);
  }
}

/// Initialize the values of output u away from the interface
void Redist3::thresholdAwayFromInterface()
{
  state.fillWithValue(false);
  for (idx_t ind : bnd)
    state.put(true, ind);
  for(size_t ii=0;ii<state.getN();++ii)
    if(!state.get(ii))
      u[ii] = mysign(u0[ii])*thres;
}

/// Perform the directional optimization routine
void Redist3::directionalOptimization()
{
  setInterfaceValuesDO(); 
  thresholdAwayFromInterface();

  for (idx_t ind : bnd)
    updateAndAddNeighborsToHeapDO(ind);

  Helt3 helt;
  while(h.popFromHeap(helt)) {
    if(!state.get(helt.second.i) && (std::abs(helt.first) < thres)) {// don't update if outside threshold cutoff, or already fixed
      applyResult(helt);
      state.put(true, helt.second.i);
      updateAndAddNeighborsToHeapDO(helt.second.i);
    }
  }
}

/// Update a single location with updated signed distance function value and closest point information if applicable
/// \param[in] h : The heap element to update with
void Redist3::applyResult(Helt3 const &h)
{
  u[h.second.i] = mysign(u0[h.second.i])*h.first;
  if(flag > 1) {
    cpx[h.second.i] = h.second.aux[0];
    cpy[h.second.i] = h.second.aux[1];
    cpz[h.second.i] = h.second.aux[2];
  }
}

/// Perform directional optimization at a single location
/// \param[in] idx : Location to perform directional optimization at
/// \return          Heap element containing the directional optimization result
Helt3 Redist3::performDO(idx_t const idx)
{ // make a guess for interface location and call performDO(const int ind, const double *cpguess)
  Point3 grad({pd(u.getxp(idx),u.getxm(idx)) / (2. * dx),
        pd(u.getyp(idx),u.getym(idx)) / (2. * dy), 
        pd(u.getzp(idx),u.getzm(idx)) / (2. * dz)});
  Point3 cpguess = lineSearch(idx, grad);
  return performDO(idx, cpguess);
}

/// Perform directional optimization near the interface
/// \param[in] idx : Location to perform directional optimization at
/// \return          Heap element containing the directional optimization result
Helt3 Redist3::performDOSurf(idx_t const idx)
{ // make a guess for interface location and call performDO
  // make initial, gradient-based guess.
  Point3 grad0({pdl(u0.getxp(idx),u0.getxm(idx),u0.lenx()) / (2. * dx),
        pdl(u0.getyp(idx),u0.getym(idx),u0.leny()) / (2. * dy),
        pdl(u0.getzp(idx),u0.getzm(idx),u0.lenz()) / (2. * dz)});
  Point3 grad = grad0;
  Point3 cpguess = lineSearch(idx,grad);
  Helt3 h = performDO(idx, cpguess);

  std::array<Point3, 26> grad2all({{{1., 0., 0.},
      {-1., 0., 0.},
      {0., 1., 0.},
      {0., -1., 0.},
      {0., 0., 1.},
      {0., 0., -1.},
      {1., 1., 0.},
      {1., -1., 0.},
      {1., 0., 1.},
      {1., 0., -1.},
      {-1., 1., 0.},
      {-1., -1., 0.},
      {-1., 0., 1.},
      {-1., 0., -1.},
      {0., 1., 1.},
      {0., 1., -1.},
      {0., -1., 1.},
      {0., -1., -1.},
      {1., 1., 1.},
      {1., 1., -1.},
      {1., -1., 1.},
      {1., -1., -1.},
      {-1., 1., 1.},
      {-1., 1., -1.},
      {-1., -1., 1.},
      {-1., -1., -1.}}});
  for (Point3 &grad2 : grad2all) {
    if(std::inner_product(grad2.begin(), grad2.end(), grad0.begin(), 0.) >= 0.) {
      Point3 cpguess2 = lineSearch(idx,grad2);
      Helt3 const h1 = performDO(idx, cpguess2);
      if( std::abs(h1.first) < std::abs(h.first) )
        h = h1;
    } 
  }
  return h;
}

/// Perform directional optimization at a single location
/// \param[in] idx     : Location to perform directional optimization at
/// \param[in] cpguess : Initial guess for closest point on the interface
/// \return              Heap element containing the directional optimization result
Helt3 Redist3::performDO(idx_t const idx, Point3 const &cpguess)
{
  // optimize over directions via line search
  Point3 const x0({u.getX(idx), u.getY(idx), u.getZ(idx)});
  unsigned int const MAXLOOPS = 5;
  unsigned int const MAXSEARCHES = 5;
  double const TOL = std::max( (0.01 / std::pow(mymax3(m,n,k),flag)), 5.0e-16);
  
  Point3 rvec, v1, v2;
  PrPtD xdt;
  std::array<PrPtD, 5> xxdd;

  xxdd[0].first = cpguess;
  xxdd[0].second = dist(x0, xxdd[0].first, u.lenx(), u.leny(), u.lenz());

  // set angle increment
  double delta;
  if(xxdd[0].second > mymax3(dx,dy,dz))
    delta = std::asin(mymax3(dx,dy,dz)/xxdd[0].second);
  else // interface is nearby, look over a wide range
    delta = PI2;
  
  double nn = findNborDirections(x0,xxdd,delta,rvec,v1,v2);
  // bracket the interface in the search directions
  for (PrPtD &elt : xxdd)
    elt.second = search1D(idx, elt.first);

  auto const MinElt([](std::array<PrPtD, 5> const &xxdd)->std::array<PrPtD, 5>::const_iterator {
      return std::min_element(xxdd.begin(), xxdd.end(),  [](std::pair<Point3, double> const &a, std::pair<Point3, double> const &b)->bool { return a.second < b.second; }); });

  auto it = MinElt(xxdd);
  
  for(unsigned int count0=0; count0<MAXLOOPS; ++count0)
  { // (1) center on minimum; (2) perform Newton step
    unsigned int count = 0;
    while( (it != xxdd.begin()) && (count++ < MAXSEARCHES)) {
      if(delta < PI4-1e-4)
        delta *= 2.;
      xxdd[0] = *it;
      /* update directions around new center location */
      nn = findNborDirections(x0,xxdd,delta,rvec,v1,v2);
      for (auto it2 = xxdd.begin()+1; it2 != xxdd.end(); ++it2)
        it2->second = search1D(idx,it2->first);
      it = MinElt(xxdd);
    }
    // perform Newton step only if dd values are successfully computed
    if(std::max(std::max(xxdd[1].second,xxdd[2].second),std::max(xxdd[3].second,xxdd[4].second)) < 0.9f) {
      double const us = (xxdd[1].second-xxdd[3].second)/2.0;
      double const ut = (xxdd[2].second-xxdd[4].second)/2.0;
      double uss= (xxdd[1].second-2.0*xxdd[0].second+xxdd[3].second);
      double utt= (xxdd[2].second-2.0*xxdd[0].second+xxdd[4].second);
      if(std::abs(uss) < 1e-16)
        uss = 1e-16;
      if(std::abs(utt) < 1e-16)
        utt = 1e-16;
      double const deltas = -us/uss*delta;
      double const deltat = -ut/utt*delta;
      xdt.first[0] = dinrange2(x0[0] + nn/2.0*(rvec[0]*cos(deltas)+v1[0]*sin(deltas)+rvec[0]*cos(deltat)+v2[0]*sin(deltat)));
      xdt.first[1] = dinrange2(x0[1] + nn/2.0*(rvec[1]*cos(deltas)+v1[1]*sin(deltas)+rvec[1]*cos(deltat)+v2[1]*sin(deltat)));
      xdt.first[2] = dinrange2(x0[2] + nn/2.0*(rvec[2]*cos(deltas)+v1[2]*sin(deltas)+rvec[2]*cos(deltat)+v2[2]*sin(deltat)));
      xdt.second = search1D(idx,xdt.first);
      
      if(std::abs(xdt.second-xxdd[0].second) < TOL)
        count0 = MAXLOOPS; // no further improvement available
      else if(xdt.second<xxdd[0].second) 
        xxdd[0] = xdt;
    }
    if(count0 < MAXLOOPS) {
      delta /= 2.;
      nn = findNborDirections(x0,xxdd,delta,rvec,v1,v2);
      for (PrPtD &elt : xxdd)
        elt.second = search1D(idx, elt.first);
    }
  }
  return std::make_pair(xxdd[0].second + std::abs(u0.interpolate(xxdd[0].first[0],xxdd[0].first[1],xxdd[0].first[2])), Aux3(idx, xxdd[0].first));
}

/// Perform a line search along the direction defined by grad for the inteface
/// \param[in] idx      : Location of the element being worked on
/// \param[in,out] grad : Direction to search in; normalized in-place
/// \return               Guess at closest point location along this line
Point3 Redist3::lineSearch(idx_t const idx, Point3 &grad)
{ // takes the given index at estimate of gradient direction (unnormalized)
  // and performs a line search along this direction for the interface.
  // Interface is expected to lie within one grid cell of idx

  normalize(grad);
  // get the correct direction to search in; 
  if(u.get(idx) < 0.) { // then we go uphill, keep gradient sign
    grad[0] *= SQRT3/static_cast<double>(n);
    grad[1] *= SQRT3/static_cast<double>(m); 
    grad[2] *= SQRT3/static_cast<double>(k);
  } else { // go downhill
    grad[0] *= -SQRT3/static_cast<double>(n);
    grad[1] *= -SQRT3/static_cast<double>(m);
    grad[2] *= -SQRT3/static_cast<double>(k);
  }
  if(findOppSign(idx,grad))
    bisect(idx,grad);
  else // need to figure out where to look
  {
    int const sgn = mysign(u.get(idx));
    std::array<idx_t,18> const nborind({u.xp(idx), u.xm(idx),
          u.yp(idx), u.ym(idx),
          u.zp(idx), u.zm(idx),
          u.yp(u.xp(idx)), u.yp(u.xm(idx)),
          u.ym(u.xp(idx)), u.ym(u.xm(idx)),
          u.zp(u.yp(idx)), u.zp(u.ym(idx)),
          u.zm(u.yp(idx)), u.zm(u.ym(idx)),
          u.xp(u.zp(idx)), u.xp(u.zm(idx)),
          u.xm(u.zp(idx)), u.xm(u.zm(idx))});
    for (idx_t const ind : nborind) {
      if(mysign(u[ind]) != sgn) {
        Point3 cpguess({u.getX(ind), u.getY(ind), u.getZ(ind)});
        bisect(idx,cpguess);
        return cpguess;
      }
    }
  }
  return grad;
}

/// Find a point across the interface in direction of input guess
/// \param[in]     idx   : Point to search for the interface from
/// \param[in/out] guess : Input as search direction with distance estimate, return the point across the interface along this line
/// \return                True if search succeeds
bool Redist3::findOppSign(idx_t const idx, Point3 &guess)
{ // searches in the direction initially given by "guess" to find x s.t. u(x) has opposite sign
  // from from u(idx)
  // input suggests to look at u[y+guess], where y is spatial coordinates of idx. 
  // returns location x in guess.
  int const sgn = mysign(u0[idx]);
  Point3 inguess = guess;

  for(int ii=1; ii<=5; ++ii)
  {
    double const dd = static_cast<double>(ii);
    guess[0] = dinrange2l(u0.getX(idx)+dd*inguess[0],u0.lenx());
    guess[1] = dinrange2l(u0.getY(idx)+dd*inguess[1],u0.leny());
    guess[2] = dinrange2l(u0.getZ(idx)+dd*inguess[2],u0.lenz());
    double const val1 = u0.interpolate(guess[0],guess[1],guess[2]);
    if(sgn != mysign(val1))
      return true;
  }
  return false;
}

/// Look about location cpguess along line from idx and find (small) positive and negative values as input to bisection.
/// \param[in]  idx     : Element to search from 
/// \param[in]  cpguess : Point across interface from this location
/// \param[out] xm      : Location with negative distance value near the interface along the line between point corersponding to idx and cpguess
/// \param[out] xp      : Location with positive distance value near the interface along the line between point corersponding to idx and cpguess
/// \return               True if the bracketing procedure succeeds
bool Redist3::bracket(idx_t const idx, Point3 const &cpguess, BracketPair3 &xm, BracketPair3 &xp)
{
  // Returns false if it's much further than initially guessed to interface in this direction (e.g. more
  // than one grid cell further), and true otherwise.
  // Assumed that xm, xp, both have room for two doubles

  Point3 const x0({u.getX(idx), u.getY(idx), u.getZ(idx)});
  double const uv0 = u0[idx];
  if(mysign(uv0) == 0) { // already on interface
    xm = std::make_pair(x0, uv0);
    xp = xm;
    return true;
  }

  Point3 dir({pdl(cpguess[0],x0[0],u.lenx()), pdl(cpguess[1],x0[1],u.leny()), pdl(cpguess[2],x0[2],u.lenz())});
  normalize(dir);

  double ug = u0.interpolate(cpguess[0],cpguess[1],cpguess[2]);
  if(mysign(uv0 * ug) == 0) { // sign(ug) == 0;
    xm = std::make_pair(cpguess, ug);
    xp = xm;
    return true;
  }
  double const dr = mymax3(dx,dy,dz);
  if(mysign(uv0 * ug) == 1) { // same side of interface
    double const ug2 = u0.interpolate(cpguess[0]+dir[0]*dr,cpguess[1]+dir[1]*dr,cpguess[2]+dir[2]*dr);
    if(mysign(uv0 * ug2) == 1) // still on same side of interface
      return false;
    else {
      xp = std::make_pair(cpguess, ug);
      xm = std::make_pair(Point3({cpguess[0]+dir[0]*dr, cpguess[1]+dir[1]*dr, cpguess[2]+dir[2]*dr}), ug2);
      if ( (mysign(ug)==-1) || ((mysign(ug)==0) && (mysign(ug2)==-1)) )
        std::swap(xp,xm);
      return true;
    }
  } else { // x0 and cpguess are on opposite sides of interface, work backwards
    Point3 xyz({cpguess[0]-dir[0]*dr, cpguess[1]-dir[1]*dr, cpguess[2]-dir[2]*dr});
    double ug2 = u0.interpolate(xyz[0],xyz[1],xyz[2]);
    unsigned int tries = 0;
    while((mysign(ug*ug2) == 1) && (tries++ < 10)) {
      xyz[0] -= dir[0]*dr;
      xyz[1] -= dir[1]*dr;
      xyz[2] -= dir[2]*dr;
      ug2 = u0.interpolate(xyz[0],xyz[1],xyz[2]);
    }
    if(tries == 10)
      return(false);

    xp = std::make_pair(cpguess, ug);
    xm = std::make_pair(xyz, ug2);
    if ( (mysign(ug)==-1) || ((mysign(ug)==0) && (mysign(ug2)==1)) )
      std::swap(xp,xm);
    return true;
  }
}

///  Tries to find the zero between [idx] and (guess) by simple linear approximation. Assumes that u0[idx] and u[guess] have different signs.
/// \param[in]     idx   : Location to search from
/// \param[in,out] guess : Initial location across the interface, returns as location of the interface along this line
void Redist3::bisect(idx_t const idx, Point3 &guess)
{ 
  Point3 xm, xp;
  double um, up;
  assert(mysign(u0[idx]*u0.interpolate(guess[0],guess[1],guess[2])) != 1);

  if(mysign(u0[idx]) == 0) {
    guess = {u.getX(idx), u.getY(idx), u.getZ(idx)};
    return;
  } else if(mysign(u0[idx]) == 1) {
    xp = {u.getX(idx), u.getY(idx), u.getZ(idx)};
    up = u0[idx];
    xm = guess;
    um = u0.interpolate(xm[0],xm[1],xm[2]);
  } else {
    xm = {u.getX(idx), u.getY(idx), u.getZ(idx)};
    um = u0[idx];
    xp = guess;
    up = u0.interpolate(xp[0],xp[1],xp[2]);
    if(up <= -1.0e-15)
      std::cout << " up = " << up << ", um = " << um << std::endl;
  }

  if(up==um) {
    guess = xp;
    return;
  }
  if(std::isnan(up/(up-um)))
  {
    std::cout << " up = " << up << ", um = " << um << ". Setting um = -1" << std::endl;
    um = -1.;
  }
  Point3 xt = ccomb(xm,xp,up/(up-um),u0.lenx(),u0.leny(),u0.lenz());
  assert(up>=0.); 
  assert(um<=0.);
  guess = xt;
  return;
}

///  Tries to find the zero between xp and xm by simple linear approximation. Assumes that u(xp) and u(xm) have different signs.
/// \param[out] result : Closest point location approximation
/// \param[in]  xm     : Negative side location and distance function value
/// \param[in]  xp     : Positive side location and distance function value
/// \return              Interpolation value at result location
double Redist3::bisect(Point3 &result, BracketPair3 const &xm, BracketPair3 const &xp)
{ // bisection with linear approximation
  double up = xp.second; //u0.interpolate(xp[0],xp[1]);
  double um = xm.second; //u0.interpolate(xm[0],xm[1]);
  assert(up>=0.); 
  assert(um<=0.);
  if(up==um) {
    result = xp.first;
    return(up);
  }
  result = ccomb(xm.first, xp.first, up/(up-um), u0.lenx(), u0.leny(), u0.lenz());
  return u0.interpolate(result[0],result[1],result[2]); 
}

/// Find directions at angles +/- delta from vector xc-x0 with same length
/// \param[in]     x0    : Initial point
/// \param[in,out] xxdd  : Initial search direction, and array of search points as return
/// \param[in]     delta : Angle to sweep through
/// \param[out]    rvec  : ONB vector in direction of xxdd[0].first - x0
/// \param[out]    v1    : ONB vector
/// \param[out]    v2    : ONB vector
/// \return Length of xxdd[0].first - x0
double Redist3::findNborDirections(Point3 const &x0, std::array<PrPtD, 5> &xxdd, double const delta, Point3 &rvec, Point3 &v1, Point3 &v2)
{
  double const s = std::sin(delta);
  double const c = std::cos(delta);

  rvec[0] = pdl(xxdd[0].first[0],x0[0],u.lenx());
  rvec[1] = pdl(xxdd[0].first[1],x0[1],u.leny());
  rvec[2] = pdl(xxdd[0].first[2],x0[2],u.lenx());
  double const nn = std::max(sqrt(rvec[0]*rvec[0]+rvec[1]*rvec[1]+rvec[2]*rvec[2]), 1e-16);
  normalize(rvec);
  orthoVecs(rvec,v1,v2);
  
  xxdd[1].first[0] = dinrange2l(x0[0] + nn*(rvec[0]*c+v1[0]*s), u.lenx());
  xxdd[1].first[1] = dinrange2l(x0[1] + nn*(rvec[1]*c+v1[1]*s), u.leny());
  xxdd[1].first[2] = dinrange2l(x0[2] + nn*(rvec[2]*c+v1[2]*s), u.lenz());
  xxdd[2].first[0] = dinrange2l(x0[0] + nn*(rvec[0]*c+v2[0]*s), u.lenx());
  xxdd[2].first[1] = dinrange2l(x0[1] + nn*(rvec[1]*c+v2[1]*s), u.leny());
  xxdd[2].first[2] = dinrange2l(x0[2] + nn*(rvec[2]*c+v2[2]*s), u.lenz());
  xxdd[3].first[0] = dinrange2l(x0[0] + nn*(rvec[0]*c-v1[0]*s), u.lenx());
  xxdd[3].first[1] = dinrange2l(x0[1] + nn*(rvec[1]*c-v1[1]*s), u.leny());
  xxdd[3].first[2] = dinrange2l(x0[2] + nn*(rvec[2]*c-v1[2]*s), u.lenz());
  xxdd[4].first[0] = dinrange2l(x0[0] + nn*(rvec[0]*c-v2[0]*s), u.lenx());
  xxdd[4].first[1] = dinrange2l(x0[1] + nn*(rvec[1]*c-v2[1]*s), u.leny());
  xxdd[4].first[2] = dinrange2l(x0[2] + nn*(rvec[2]*c-v2[2]*s), u.lenz());
  return(nn);
}

/// Perform the 1D search for the interface from point at idx in direction of x
/// \param[in]     idx : Point to compute distance from
/// \param[in,out] x   : Initial search point, returns the interface point in this direction
/// \return              Distance to interface point in this direction
double Redist3::search1D(idx_t const idx, Point3 &x)
{
  // bracket the interface in the search directions
  BracketPair3 xm, xp; // xm and xp also return values from u0
  if (bracket(idx,x,xm,xp)) {
    double const val1 = bisect(x,xm,xp);
    Point3 const x0({u.getX(idx), u.getY(idx), u.getZ(idx)});
    return dist(x0,x,u.lenx(),u.leny(),u.lenz()) + val1 * mysign(u0[idx]);
  }
  return mymax3(u.lenx(), u.leny(), u.lenz());
}

/// Write u out to a pointer (e.g., for MATLAB)
/// \param[out] v : Double array to write to. Memory assumed to be pre-allocated
void Redist3::dump_u(double *v)
{ // assumes sufficient memory is allocated into v
  std::memcpy(v, u.returnData().data(), N * sizeof(double));
}

/// Write cpx, cpy, cpz out to pointers (e.g., for MATLAB)
/// \param[out] cpx_d : Double array to write cpx to. Memory assumed to be pre-allocated
/// \param[out] cpy_d : Double array to write cpy to. Memory assumed to be pre-allocated
/// \param[out] cpz_d : Double array to write cpz to. Memory assumed to be pre-allocated
void Redist3::dump_cp(double *cpx_d, double *cpy_d, double *cpz_d)
{ // assume sufficient memory is allocated into cpx, cpy, cpz
  std::memcpy(cpx_d, cpx.returnData().data(), N * sizeof(double));
  std::memcpy(cpy_d, cpy.returnData().data(), N * sizeof(double));
  std::memcpy(cpz_d, cpz.returnData().data(), N * sizeof(double));
}

/// Perform second order iterations for non-directional optimization -- NOT enabled
void Redist3::secondOrderIterations()
{
  if (!bWarn) {
    std::cout << "WARNING: secondOrderIterations are not enabled in 3D." << std::endl;
    bWarn = true;
  }
}

/// Return true if sign of u[idx] differs from sign of any of its 4 neighbors
/// \param[in] idx : Location to compare signs around
/// \return          True if any four-neighbor of idx has different signs in u0 feom u0[idx]
bool Redist3::diffSign(idx_t const idx)
{
  int const sui = mysign(u[idx]);
  return (sui != mysign(u.getxp(idx))) || (sui != mysign(u.getxm(idx))) ||
    (sui != mysign(u.getym(idx))) || (sui != mysign(u.getym(idx))) ||
    (sui != mysign(u.getzp(idx))) || (sui != mysign(u.getzm(idx)));
}

/// Return a const reference to the signed distance function
/// \return Const reference to the signed distance function
Array3D<double> const& Redist3::dump_u()
{
  return u;
}

