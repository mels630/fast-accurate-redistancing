#include "redist.hpp"
#include "toolbox.hpp"
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cstring>

using BracketPair = std::pair<Point, double>;

/// Constructor with full arguments
/// \param[in] _i   : Index value
/// \param[in] _aux : Auxilary double array values
Aux::Aux(idx_t const _i, Point const &_aux) :
  i(_i),
  aux(_aux)
{ }

/// Constructor from index only
/// \param[in] _i   : Index value
Aux::Aux(idx_t const _i) :
  i(_i),
  aux()
{ } 

/// Constructor
/// \param[in] _u     : Input level set function
/// \param[in] _width : Thresholding width (pixels)
/// \param[in] _flag  : Interpolation order flag
Redist::Redist(Array2D<double> const &_u, idx_t const _width, int const _flag) :
  Redist(_u, 1./static_cast<double>(n), 1./static_cast<double>(m), _width, _flag)
{ }

/// Constructor
/// \param[in] _u     : Input level set function
/// \param[in] _dx    : x grid spacing
/// \param[in] _dy    : y grid spacing
/// \param[in] _width : Thresholding width (pixels)
/// \param[in] _flag  : Interpolation order flag
Redist::Redist(Array2D<double> const &_u, double const _dx, double const _dy, idx_t const _width, int const _flag) :
  width(_width),
  flag(_flag),
  m(_u.getm()),
  n(_u.getn()),
  N(m*n),
  dx(_dx),
  dy(_dy),
  thres(static_cast<double>(width+1)*std::max(dx,dy)),
  cpflag((_flag==2) || (_flag==3)),
  h(10*width,cpflag),
  state(m,n),
  u0(_u,_flag),
  u(_u),
  cpx(m,n),
  cpy(m,n)
{ }

/// Perform the redistancing
void Redist::redistance()
{
  if((flag == 0) || (flag == 1)) {
    fastMarchingRedist();
    if(flag == 0)
      secondOrderIterations();
  }
  else if((flag == 2) || (flag == 3))
    directionalOptimization();
  else
    std::cout << "Do not recognize the flag. Cowardly refusing to do anything." << std::endl;
}

/// Internal function for fast marching redistancing
void Redist::fastMarchingRedist()
{
  setInterfaceValues();
  thresholdAwayFromInterface();

  for (idx_t const ix : bnd)
    updateAndAddNeighborsToHeap(ix);

  // while heap is non-empty, fix its top value, and update the neighbors
  Helt helt;
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
void Redist::updateAndAddNeighborsToHeap(idx_t const idx)
{
  std::array<idx_t, 4> const idx2arr = {u.xp(idx), u.xm(idx), u.yp(idx), u.ym(idx)};
  for(idx_t const ind : idx2arr) {
    if(!state.get(ind)) { // if true, value is already fixed
      double const dtemp = estimateUpdate(ind);
      if(std::abs(dtemp) < std::abs(u[ind])) {
        u[ind] = dtemp;
        h.addToHeap(Aux(ind), std::abs(dtemp));
      }
    }
  }
}

/// Fix the value at u[idx], update its 4 neighbors if they are not already fixed, and make sure they're in the heap at the appropriate location -- for directional optimization
/// \param[in] idx : Location to work from
void Redist::updateAndAddNeighborsToHeapDO(idx_t const idx)
{ 
  std::array<idx_t, 4> const idx2arr = {u.xp(idx), u.xm(idx), u.yp(idx), u.ym(idx)};
  for (idx_t const ind : idx2arr) {
    if(!state.get(ind)) { // if true, value is already fixed
      Helt const htemp = performDO(ind, Point({cpx[idx],cpy[idx]}));
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
double Redist::estimateUpdate(idx_t const idx)
{
  double const a = std::min(std::abs(u.getxm(idx)), std::abs(u.getxp(idx)));
  double const b = std::min(std::abs(u.getym(idx)), std::abs(u.getyp(idx)));
  if(((a-b)*(a-b))>=(dx*dx+dy*dy)) {
    if((a+dx) < (b+dy))
      if(std::abs(u.getxm(idx)) < std::abs(u.getxp(idx)))
	return(u.getxm(idx) + dx*mysign(u.getxm(idx)));
      else
	return(u.getxp(idx) + dx*mysign(u.getxp(idx)));
    else
      if(std::abs(u.getym(idx)) < std::abs(u.getyp(idx)))
	return(u.getym(idx) + dy*mysign(u.getym(idx)));
      else
	return(u.getyp(idx) + dy*mysign(u.getyp(idx)));
  } else
    return(mysign(u[idx]) * (dy*dy*a+dx*dx*b+dx*dy*sqrt(dx*dx+dy*dy-(a-b)*(a-b)))/(dx*dx+dy*dy));
}

/// Set values of output signed distance function at the interface
void Redist::setInterfaceValues()
{ 
  Array2D<int> sgn(u.getm(),u.getn());
  std::vector<int> &vSgn = sgn.returnData();
  std::vector<double> const &vu = u.returnData();
  assert(vSgn.size() == vu.size());
  std::transform(vu.begin(), vu.end(), vSgn.begin(), [](double const d)->int{ return static_cast<int>(mysign(d)); });
  for(idx_t ii=0; ii<sgn.getN(); ++ii)
    if ((sgn.getxp(ii)-sgn[ii]) || (sgn.getxm(ii)-sgn[ii]) ||
        (sgn.getyp(ii)-sgn[ii]) || (sgn.getym(ii)-sgn[ii]))
      bnd.push_back(ii);
 
  std::vector<double> dr(bnd.size());
  
  std::transform(bnd.begin(), bnd.end(), dr.begin(), [&](idx_t const ii)->double {
      // compute norm(grad u) with centered differences
      double ry = (u.getyp(ii)-u.getym(ii))/dy/2.;
      double rx = (u.getxp(ii)-u.getxm(ii))/dx/2.;
      double dr = sqrt(rx*rx+ry*ry);

      // compute norm(grad u) with one-sided differences
      rx = std::max(std::abs(u.getxp(ii)-u[ii]),std::abs(u[ii]-u.getxm(ii)))/dx;
      ry = std::max(std::abs(u.getyp(ii)-u[ii]),std::abs(u[ii]-u.getym(ii)))/dy;
      double const dr2 = sqrt(rx*rx+ry*ry);
      
      // Accept one-sided difference is much different than centered difference
      if((dr < (0.5*dr2)) || (dr > (2.0*dr2)))
	dr = dr2;
      return dr;
    });
  for(size_t ii = 0; ii<bnd.size(); ++ii)    
    u[bnd[ii]] /= dr[ii];
}

/// Set the interface values for directional optimization (must have correct order of accuracy)
void Redist::setInterfaceValuesDO()
{ 
  double meanabsbndval = 0.;
  for(idx_t ii=0; ii<N; ++ii)
    if(diffSign(ii)) {
      bnd.push_back(ii);
      meanabsbndval += std::abs(u0[ii]);
    }
  
  meanabsbndval /= static_cast<double>(bnd.size());
  // rescale u0 based on meanabsbndval
  u0 *= 0.5 * std::max(dx,dy) / meanabsbndval; // scale u0 (rationale: mean distance to interface at a grid cell should be ~ 0.5 dx)

  // compute interface values and update
  for (idx_t ind : bnd) {
    Helt const bndval = performDOSurf(ind);
    applyResult(bndval);
    state.put(true, ind);
  }
}

/// Initialize the values of output u away from the interface
void Redist::thresholdAwayFromInterface()
{
  state.fillWithValue(false);
  for (idx_t ind : bnd)
    state.put(true, ind);
  for(idx_t ii=0;ii<state.getN();++ii)
    if(!state.get(ii))
      u[ii] = mysign(u0[ii])*thres;
}

/// Perform the directional optimization routine
void Redist::directionalOptimization()
{
  setInterfaceValuesDO(); 
  thresholdAwayFromInterface();

  for (idx_t ind : bnd) 
    updateAndAddNeighborsToHeapDO(ind);
  Helt helt;
  while(h.popFromHeap(helt)) {
    if(!state.get(helt.second.i) && (std::abs(helt.first) < thres)) // don't update if outside threshold cutoff, or already fixed
    {
      applyResult(helt);
      state.put(true, helt.second.i);
      updateAndAddNeighborsToHeapDO(helt.second.i);
    }
  }
}

/// Update a single location with updated signed distance function value and closest point information if applicable
/// \param[in] h : The heap element to update with
void Redist::applyResult(Helt const &h)
{
  u[h.second.i] = mysign(u0[h.second.i])*h.first;
  if(flag > 1) {
    cpx[h.second.i] = h.second.aux[0];
    cpy[h.second.i] = h.second.aux[1];
  }
}

/// Perform directional optimization at a single location
/// \param[in] idx : Location to perform directional optimization at
/// \return          Heap element containing the directional optimization result
Helt Redist::performDO(idx_t const idx)
{ // make a guess for interface location and call performDO(const int ind, const double *cpguess)
  Point grad({pdl(u0.getxp(idx),u0.getxm(idx),u0.lenx()) / (2. * dx),
        pdl(u0.getyp(idx),u0.getym(idx),u0.leny()) / (2. * dy)}); 
  Point cpguess = lineSearch(idx, grad);
  return performDO(idx, cpguess);
}

/// Perform directional optimization near the interface
/// \param[in] idx : Location to perform directional optimization at
/// \return          Heap element containing the directional optimization result
Helt Redist::performDOSurf(idx_t const idx)
{ // make a guess for interface location and call performDO
  // make initial, gradient-based guess.
  Point grad0({pdl(u0.getxp(idx),u0.getxm(idx),u0.lenx()) / (2. * dx), pdl(u0.getyp(idx),u0.getym(idx),u0.leny()) / (2. * dy)});
  Point grad = grad0;
  Point cpguess = lineSearch(idx,grad);
  Helt h = performDO(idx, cpguess);

  std::array<Point, 8> grad2all({{{1.,0.},
        {0.,1.},
        {-1.,0.},
        {0.,-1.},
        {1.,1.},
        {-1.,1.},
        {-1.,-1.},
        {1.,-1}}});
  for(Point &grad2 : grad2all) {
    if(std::inner_product(grad2.begin(), grad2.end(), grad0.begin(), 0.) >= 0.) {
      Point cpguess2 = lineSearch(idx,grad2);
      Helt const h1 = performDO(idx, cpguess2);
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
Helt Redist::performDO(idx_t const idx, Point const &cpguess)
{
  // optimize over directions via line search
  Point const x0({u.getX(idx), u.getY(idx)});
  unsigned int const MAXLOOPS = 5;
  unsigned int const MAXSEARCHES = 5;
  // still not sure about this TOL!
  double const TOL = std::max( (0.01 * std::min(dx,dy) / std::pow(static_cast<double>(std::max(m,n)),flag)), 5.0e-16);
  
  Point xl, xr, xc, xt;
  double dl, dc, dr, dt;
  xc = cpguess;

  dc = dist(x0,xc,u.lenx(),u.leny());
  // set angle increment
  double delta;
  if(dc > std::max(dx,dy))
    delta = std::max(asin(std::max(dx,dy)/dc),PI/32.);
  else // interface is nearby, look over a wide range
    delta = PI2;

  findNborDirections(x0,xc,xl,xr,delta);
  // bracket the interface in the search directions
  dl = search1D(idx,xl);
  dr = search1D(idx,xr);

  for(unsigned int count0=0; count0<MAXLOOPS; ++count0)
  { // (1) center on minimum; (2) perform Newton step
    unsigned int count = 0;
    while( (dc > std::min(dl,dr)) && (count++ < MAXSEARCHES))
    {
      if(delta < PI4-1e-4)
        delta *= 2.;

      if(dl < dr) {
        dc = dl; xc = xl;
      } else { // dl >= dr
        dc = dr; xc = xr;
      }
      findNborDirections(x0,xc,xl,xr,delta);
      dl = search1D(idx,xl);
      dr = search1D(idx,xr);
    }
    // perform Newton step only if dl/dr are successfully computed
    if(std::max(dl,dr) < 0.9f*std::max(m,n)*std::max(dx,dy))
    {
      double ds = (dr-dl)/2.;
      double dss = (dr-2.*dc+dl);
      if(std::abs(dss) < 1e-16)
        dss = 1e-16;
      Point dir = {pdl(xc[0],x0[0],u.lenx()), pdl(xc[1],x0[1],u.leny())};
      double const delta2 = -ds/dss * delta;
      xt[0] = dinrange2l(x0[0] + cos(delta2)*dir[0] - sin(delta2)*dir[1],u.lenx());
      xt[1] = dinrange2l(x0[1] + sin(delta2)*dir[0] + cos(delta2)*dir[1],u.leny());
      dt = search1D(idx,xt);

      if(std::abs(dt-dc) < TOL)
        count0 = MAXLOOPS; // no further improvement available
      else if(dt<dc)
        dc = dt; xc = xt;
    }
    if(count0 < MAXLOOPS)
    {
      delta /= 2.;
      findNborDirections(x0,xc,xl,xr,delta);
      dl = search1D(idx,xl);
      dr = search1D(idx,xr);
    }
  }
  double val1 = u0.interpolate(xc[0],xc[1]);
  Helt lval(dc + std::abs(val1), Aux(idx));
  
  Point grad = {pdl(xc[0],x0[0],u0.lenx()), pdl(xc[1],x0[1],u0.leny())};
  normalize(grad);
  lval.second.aux[0] = xc[0] + std::abs(val1) * grad[0];
  lval.second.aux[1] = xc[1] + std::abs(val1) * grad[1];
  return(lval);
}

/// Perform a line search along the direction defined by grad for the inteface
/// \param[in] idx      : Location of the element being worked on
/// \param[in,out] grad : Direction to search in; normalized in-place
/// \return               Guess at closest point location along this line
Point Redist::lineSearch(idx_t const idx, Point &grad)
{ // takes the given index at estimate of gradient direction (unnormalized)
  // and performs a line search along this direction for the interface.
  // Interface is expected to lie within one grid cell of idx

  normalize(grad);
  // get the correct direction to search in; 
  double const val = SQRT2*std::max(dx,dy);
  if(u[idx] < 0.) { // then we go uphill, keep gradient sign
    grad[0] *= val;
    grad[1] *= val;
  } else { // go downhill 
    grad[0] *= -val;
    grad[1] *= -val;
  }
  if(findOppSign(idx,grad))
    bisect(idx,grad);
  return grad;
}

/// Find a point across the interface in direction of input guess
/// \param[in]     idx   : Point to search for the interface from
/// \param[in/out] guess : Input as search direction with distance estimate, return the point across the interface along this line
/// \return                True if search succeeds
bool Redist::findOppSign(idx_t const idx, Point &guess)
{ // searches in the direction initially given by "guess" to find x s.t. u(x) has opposite sign
  // from from u(idx)
  // input suggests to look at u[y+guess], where y is spatial coordinates of idx. 
  // returns location x in guess.
  int const sgn = mysign(u0[idx]);
  Point inguess = guess;

  for(int ii=1; ii<=5; ++ii)
  {
    double const dd = static_cast<double>(ii);
    guess[0] = dinrange2l(u0.getX(idx)+dd*inguess[0],u0.lenx());
    guess[1] = dinrange2l(u0.getY(idx)+dd*inguess[1],u0.leny());
    double const val1 = u0.interpolate(guess[0],guess[1]);
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
bool Redist::bracket(idx_t const idx, Point const &cpguess, BracketPair &xm, BracketPair &xp)
{
  // Returns false if it's much further than initially guessed to interface in this direction (e.g. more
  // than one grid cell further), and true otherwise.
  // Assumed that xm, xp, both have room for two doubles
  
  Point const x0({u.getX(idx), u.getY(idx)});
  double const uv0 = u0[idx];
  if(mysign(uv0) == 0) { // already on interface
    xm = std::make_pair(x0, uv0);
    xp = xm;
    return true;
  }

  Point dir({pdl(cpguess[0],x0[0],u.lenx()), pdl(cpguess[1],x0[1],u.leny())});
  normalize(dir);

  double const ug = u0.interpolate(cpguess[0],cpguess[1]);
  if(mysign(uv0 * ug) == 0) { // sign(ug) == 0; 
    xm = std::make_pair(cpguess, ug);
    xp = xm;
    return true;
  }
  double const dr = std::max(dx,dy);
  if(mysign(uv0 * ug) == 1) { // same side of interface
    double const ug2 = u0.interpolate(cpguess[0]+dir[0]*dr,cpguess[1]+dir[1]*dr);
    if(mysign(uv0 * ug2) == 1) // still on same side of interface
      return false;
    else {
      xp = std::make_pair(cpguess, ug);
      xm = std::make_pair(Point({cpguess[0]+dir[0]*dr, cpguess[1]+dir[1]*dr}), ug2);
      if ( (mysign(ug)==-1) || ((mysign(ug)==0) && (mysign(ug2)==-1)) )
        std::swap(xp,xm);
      return true;
    }
  } else { // x0 and cpguess are on opposite sides of interface, work backwards
    Point xy({cpguess[0]-dir[0]*dr, cpguess[1]-dir[1]*dr});
    double ug2 = u0.interpolate(xy[0],xy[1]);
    unsigned int tries = 0;
    while((mysign(ug*ug2) == 1) && (tries++ < 10)) {
      xy[0] -= dir[0]*dr;
      xy[1] -= dir[1]*dr;
      ug2 = u0.interpolate(xy[0],xy[1]);
    }
    if(tries == 10)
      return false;

    xp = std::make_pair(cpguess, ug);
    xm = std::make_pair(xy, ug2);
    if ( (mysign(ug)==-1) || ((mysign(ug)==0) && (mysign(ug2)==1)) )
      std::swap(xp,xm);
    return true;
  }
}

///  Tries to find the zero between [idx] and (guess) by simple linear approximation. Assumes that u0[idx] and u[guess] have different signs.
/// \param[in]     idx   : Location to search from
/// \param[in,out] guess : Initial location across the interface, returns as location of the interface along this line
void Redist::bisect(idx_t const idx, Point &guess)
{
  Point xm, xp;
  double um, up;
  assert(mysign(u0[idx]*u0.interpolate(guess[0],guess[1])) != 1);

  if(mysign(u0[idx]) == 0) {
    guess = {u.getX(idx), u.getY(idx)};
    return;
  } else if(mysign(u0[idx]) == 1) {
    xp = {u.getX(idx), u.getY(idx)};
    up = u0[idx];
    xm = guess;
    um = u0.interpolate(xm[0],xm[1]);
  } else {
    xm = {u.getX(idx), u.getY(idx)};
    um = u0[idx];
    xp = guess;
    up = u0.interpolate(xp[0],xp[1]);
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
  Point xt = ccomb(xm,xp,up/(up-um),u0.lenx(),u0.leny());
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
double Redist::bisect(Point &result, BracketPair const &xm, BracketPair const &xp)
{ // bisection with linear approximation
  // Assumes that u0[idx] and u[guess] have different signs.
  // Tries to find the zero between [idx] and (guess) by single linear approx.
  double up = xp.second; //u0.interpolate(xp[0],xp[1]);
  double um = xm.second; //u0.interpolate(xm[0],xm[1]);
  assert(up>=0.); 
  assert(um<=0.);
  if(up==um) {
    result = xp.first;
    return(up);
  }
  result = ccomb(xm.first, xp.first, up/(up-um), u0.lenx(), u0.leny());
  return u0.interpolate(result[0],result[1]); 
}

/// Find directions at angles +/- delta from vector xc-x0 with same length
/// \param[in]  x0    : Initial point
/// \param[in]  xc    : Point forming vector with initial point
/// \param[out] xl    : Point delta degrees "left" of xc from x0
/// \param[out] xr    : Point delta degrees "right" of xc from x0
/// \param[in]  delta : Angle to sweep through
void Redist::findNborDirections(Point const &x0,
                                Point const &xc,
                                Point       &xl,
                                Point       &xr,
                                double                const  delta)
{
  double const s = sin(delta);
  double const c = cos(delta);
  double const rx = pdl(xc[0],x0[0],u.lenx());
  double const ry = pdl(xc[1],x0[1],u.leny());
  xl[0] = dinrange2l(x0[0] + c*rx + s*ry,u.lenx());
  xl[1] = dinrange2l(x0[1] - s*rx + c*ry,u.leny());
  xr[0] = dinrange2l(x0[0] + c*rx - s*ry,u.lenx());
  xr[1] = dinrange2l(x0[1] + s*rx + c*ry,u.leny());
}

/// Perform the 1D search for the interface from point at idx in direction of x
/// \param[in]     idx : Point to compute distance from
/// \param[in,out] x   : Initial search point, returns the interface point in this direction
/// \return              Distance to interface point in this direction
double Redist::search1D(idx_t const idx, Point &x)
{
  // bracket the interface in the search directions
  BracketPair xm, xp; // xm and xp also return values from u0
  // once bracketed successfully, bisect
  if(bracket(idx,x,xm,xp)) {
    double const val1 = bisect(x,xm,xp);
    Point const x0({u.getX(idx), u.getY(idx)});
    return dist(x0,x,u.lenx(),u.leny()) + val1 * mysign(u0[idx]);
  }
  return std::max(u.lenx(), u.leny());
}

/// Write u out to a pointer (e.g., for MATLAB)
/// \param[out] v : Double array to write to. Memory assumed to be pre-allocated
/// \todo Fix this function / memcpy
void Redist::dump_u(double *v)
{ // assumes sufficient memory is allocated into v
  std::memcpy(v, u.returnData().data(), N * sizeof(double));
}

/// Write cpx, cpy out to pointers (e.g., for MATLAB)
/// \param[out] cpxArr : Double array to write cpx to. Memory assumed to be pre-allocated
/// \param[out] cpyArr : Double array to write cpy to. Memory assumed to be pre-allocated
/// \todo Fix this function / memcpy
void Redist::dump_cp(double *cpxArr, double *cpyArr)
{ // assume sufficient memory is allocated into cpxArr, cpyArr
  std::memcpy(cpxArr, cpx.returnData().data(), N * sizeof(double));
  std::memcpy(cpyArr, cpy.returnData().data(), N * sizeof(double));
}

/// Perform second order iterations for non-directional optimization
void Redist::secondOrderIterations()
{
  double const me = 1./ 10000. * std::max(dx,dy);
  double const thres2 = thres - dx;
  double const dt = 1. / 5. * std::max(dx,dy);
  // Step 1: Figure out which indices to do the iterations on
  std::vector<int> idxc;
  for(idx_t ii=0; ii<N; ++ii)
    if(std::abs(u[ii]) < thres2)
      idxc.emplace_back(ii);

  std::vector<double> G(idxc.size(), 0.);
  std::vector<bool> df;
  df.reserve(idxc.size());
  for (idx_t const ix : idxc)
    df.emplace_back(diffSign(ix));

  idx_t const iter2o = 3*width;
  std::array<double, 5> t0; // temporaries
  std::array<double, 4> t1, t2, t3; // temporaries
  for(idx_t iter=1; iter<=iter2o; ++iter)
  {
    for(size_t ii=0; ii<idxc.size(); ++ii)
      if(!df[ii])
      {
        idx_t const jj = idxc[ii];
        /* x */
        t0[0] = u.getxm(u.xm(jj));
        t0[1] = u.getxm(jj);
        t0[2] = u[jj];
        t0[3] = u.getxp(jj);
        t0[4] = u.getxp(u.xp(jj));
        t1[0] = dx; t1[1] = dx; t1[2] = dx; t1[3] = dx;
        if((mysigntol(t0[0],me)*mysigntol(t0[1],me))==-1) {
          t1[0] = dx * t0[1] / (t0[1]-t0[0]);
          t0[0] = 0.;
        }
        if((mysigntol(t0[3],me)*mysigntol(t0[4],me))==-1) {
          t1[3] = -dx * t0[3] / (t0[4]-t0[3]);
          t0[4] = 0.;
        }
        
        t2[0] = (t0[1]-t0[0]) / t1[0];
        t2[1] = (t0[2]-t0[1]) / t1[1];
        t2[2] = (t0[3]-t0[2]) / t1[2];
        t2[3] = (t0[4]-t0[3]) / t1[3];
        t3[0] = (t2[1]-t2[0]) / (t1[1]+t1[0]);
        t3[1] = (t2[2]-t2[1]) / (t1[2]+t1[1]);
        t3[2] = (t2[3]-t2[2]) / (t1[3]+t1[2]);
        double const a = t2[1]+mm(t3[0],t3[1])*t1[1];
        double const b = t2[2]-mm(t3[1],t3[2])*t1[2];          
        
        /* y */
        t0[0] = u.getym(u.ym(jj));
        t0[1] = u.getym(jj);
        t0[2] = u[jj];
        t0[3] = u.getyp(jj);
        t0[4] = u.getyp(u.yp(jj));
        t1[0] = dy; t1[1] = dy; t1[2] = dy; t1[3] = dy; 
        if((mysigntol(t0[0],me)*mysigntol(t0[1],me))==-1) {
          t1[0] = dy * t0[1] / (t0[1]-t0[0]);
          t0[0] = 0.;
        }
        if((mysigntol(t0[3],me)*mysigntol(t0[4],me))==-1) {
          t1[3] = -dy * t0[3] / (t0[4]-t0[3]);
          t0[4] = 0.;
        }
        
        t2[0] = (t0[1]-t0[0]) / t1[0];
        t2[1] = (t0[2]-t0[1]) / t1[1];
        t2[2] = (t0[3]-t0[2]) / t1[2];
        t2[3] = (t0[4]-t0[3]) / t1[3];
        t3[0] = (t2[1]-t2[0]) / (t1[1]+t1[0]);
        t3[1] = (t2[2]-t2[1]) / (t1[2]+t1[1]);
        t3[2] = (t2[3]-t2[2]) / (t1[3]+t1[2]);
        double const c = t2[1]+mm(t3[0],t3[1])*t1[1];
        double const d = t2[2]-mm(t3[1],t3[2])*t1[2];          
        
        if(u[jj] > me) {
          double gg = 0.;
          if (a > 0.)
            gg = a*a;
          if (b > 0.)
            gg = std::max(gg, b*b);
          if (c > 0.)
            gg = std::max(gg, c*c);
          if (d > 0.)
            gg = std::max(gg, d*d);
          G[ii] = std::sqrt(gg) - 1.;
        } else if(u[jj] < -me) {
          double gg = 0.;
          if (a < 0.)
            gg = a*a;
          if (b < 0.)
            gg = std::max(gg, b*b);
          if (c < 0.)
            gg = std::max(gg, c*c);
          if (d < 0.)
            gg = std::max(gg, d*d);
          G[ii] = std::sqrt(gg) - 1.;
        } else
          G[ii] = 0.;
      } // if non-interface cell - compute values to update
    for(size_t ii=0; ii<idxc.size(); ++ii)
      if(!df[ii]) // update all non-interface cells simultaneously
        u[idxc[ii]] = u[idxc[ii]]-dt*mysign(u[idxc[ii]])*G[ii];
  } // 2nd order iterations
}

/// Return true if sign of u[idx] differs from sign of any of its 4 neighbors
/// \param[in] idx : Location to compare signs around
/// \return          True if any four-neighbor of idx has different signs in u0 feom u0[idx]
bool Redist::diffSign(idx_t const idx)
{ // 
  int const sui = mysign(u0[idx]);
  return ((sui != mysign(u0.getxp(idx))) || (sui != mysign(u0.getyp(idx))) || (sui != mysign(u0.getxm(idx))) || (sui != mysign(u0.getym(idx))));
}

/// Return a const reference to the signed distance function
/// \return Const reference to the signed distance function
Array2D<double> const& Redist::dump_u()
{
  return u;
}

