#include "reinit.hpp"

bool VERBOSITY = false;

Reinit::Reinit(const Array2D<double> &_tpx, const Array2D<double> &_tpy, const int _width, const int _flag) :
  width(_width),
  flag(_flag),
  m(_tpx.getm()),
  n(_tpy.getn()),
  N(m*n),
  dx(1.0f/static_cast<double>(n)),
  dy(1.0f/static_cast<double>(m)),
  thres(static_cast<double>(width+1)*mymax(dx,dy)),
  h(10*width,true),
  state(m,n),
  //tpx(_tpx,_flag),
  //tpy(_tpy,_flag),
  cpx(_tpx,_flag),
  cpy(_tpy,_flag),
  dd(m,n)
{
  // set dd according to input values tpx/tpy
  for(int ii=0; ii<N; ++ii)
  {
    struct point xx; xx.x = cpx.getX(ii); xx.y = cpx.getY(ii);
    struct point cp; cp.x = cpx[ii]; cp.y = cpy[ii];
    dd[ii] = dist(xx,cp);
  }
}


void Reinit::reinit()
{
  /*
  for(int ii=0; ii<5; ++ii)
  {
    struct point x;
    x.x = static_cast<double>(ii)/5.0f;
    x.y = static_cast<double>(ii)/6.0f;
    struct point d;
    d.x = static_cast<double>(ii)/11.3151f;
    d.y = static_cast<double>(ii)/7.4319581f;
    const double t = 3.3f / static_cast<double>(n);
    double ht, htder;
    computeHtAndDer(x,d,t,ht,htder);
  }
  */
  //mexPrintf(" cpy.maxval() = %f, cpy.minval() = %f\n",cpy.maxval(),cpy.minval());

//void Reinit::computeHtAndDer(const struct point x, const struct point d, const double t, double &ht, double &htder) const
  struct point x; x.x = 0.25f; x.y = 0.5f;
  struct point d; d.x = -1.0f; d.y = 1.0f;
  double t0 = -5.f*mymax(dx,dy);
  mexPrintf("data1 = [\n");
  for(int ii=0; ii<10; ++ii)
  {
    double tt = t0 + static_cast<double>(ii) * mymax(dx,dy);
    double ht,htder;
    computeHtAndDer(x,d,tt,ht,htder);
    mexPrintf(" %f %f %f %f\n",tt,ht,htder,dd.interpolate(x.x+tt*d.x,x.y+tt*d.y));
  }
  mexPrintf("];\n");

  state.fillWithValue(false);
  tempInterfaceValues();

  thresholdAwayFromInterface();

  mexPrintf("data2 = [\n");
  for(int ii=0; ii<10; ++ii)
  {
    double tt = t0 + static_cast<double>(ii) * mymax(dx,dy);
    double ht,htder;
    computeHtAndDer(x,d,tt,ht,htder);
    mexPrintf(" %f %f %f %f\n",tt,ht,htder,dd.interpolate(x.x+tt*d.x,x.y+tt*d.y));
  }
  mexPrintf("];\n");



  for(size_t ii=0;ii<bnd.size();++ii)
    updateAndAddNeighborsToHeapDO(bnd[ii]);

  while(h.numberOfElements() > 0)
  {
    struct helt h1 = h.popFromHeap();
    if((state.get(h1.i)==false) && (fabs(h1.d) < thres)) // don't update if outside threshold cutoff, or already fixed
    {
      //mexPrintf(" idx = %d, h1.d = %f, cpx = %f, cpy = %f, dist to circ = %f\n",h1.i,h1.d,h1.aux[0],h1.aux[1],.25f-sqrt((h1.aux[0]-.46875f)*(h1.aux[0]-.46875f) + (h1.aux[1]-.46875f)*(h1.aux[1]-.46875f)));
      applyResult(h1);
      state.put(true,h1.i);
      updateAndAddNeighborsToHeapDO(h1.i);
    }
  }
}

void Reinit::updateAndAddNeighborsToHeapDO(const int idx)
{ // fix the value at u[idx], update its 4 neighbors if they are not already fixed, and make sure they're in the heap at the appropriate location
  int idx2arr[4];
  idx2arr[0] = cpx.xp(idx);
  idx2arr[1] = cpx.xm(idx); 
  idx2arr[2] = cpx.yp(idx);
  idx2arr[3] = cpx.ym(idx);
  for(int ii=0;ii<4;++ii)
  {
    if(state.get(idx2arr[ii]) == false) // if true, value is already fixed
    { 
      // perform DO starting from cp[idx], store guess in htemp
      struct point cpidx;
      cpidx.x = cpx[idx]; cpidx.y = cpy[idx];
      struct helt htemp = performDO(idx2arr[ii],cpidx);
      if(fabs(htemp.d) < fabs(distToCP(idx2arr[ii],cpx[idx2arr[ii]],cpy[idx2arr[ii]])))
      {
        applyResult(htemp);
        h.addToHeap(htemp.i,fabs(htemp.d),htemp.aux);
      }
    }
  }
}

void Reinit::applyResult(const struct helt &h)
{
  cpx[h.i] = h.aux[0];
  cpy[h.i] = h.aux[1];
  double xx[2]; xx[0] = cpx.getX(h.i); xx[1] = cpx.getY(h.i);
  dd[h.i] = dist(h.aux,xx);
  cpx.freeQeIdx(h.i);
  cpy.freeQeIdx(h.i); // reset interpolation coefficients for cpx/cpy once
                      // their values have been changed
  dd.freeQeIdx(h.i);
}
  
void Reinit::thresholdAwayFromInterface()
{
  for(int ii=0;ii<state.getN();++ii)
    if(!state.get(ii))
    {
      cpx[ii] = 10.0f;
      cpy[ii] = 10.0f;
      dd[ii] = 10.0f;
      cpx.freeQeIdx(ii);
      cpy.freeQeIdx(ii);
      dd.freeQeIdx(ii);
    }
}
  

void Reinit::tempInterfaceValues()
{
  for(int ii=0; ii<N; ++ii)
  {
    if(myDist(cpx[ii],cpy[ii],static_cast<double>(ii/m)*dx,static_cast<double>(ii%m)*dy) < 3.0f*mymax(dx,dy))
    {
      //mexPrintf(" ii = %d, cp[ii] = (%f,%f)\n",ii,cpx[ii],cpy[ii]);
      bnd.push_back(ii);
      state.put(true,ii);
      struct point p1; p1.x = cpx[ii]; p1.y = cpy[ii];
      struct point p2; p2.x = static_cast<double>(ii/m)*dx; p2.y = static_cast<double>(ii%m)*dy;
      dd[ii] = dist(p1,p2);
      dd.freeQeIdx(ii);
      //state[ii] = true; // doesn't work for some reason ... why???
    }
  }
}

void Reinit::dump_cp(double *dcpx, double *dcpy) const
{ // assume sufficient memory is allocated into dcpx, dcpy
  for(int ii=0;ii<N;++ii)
    dcpx[ii] = cpx[ii];
  for(int ii=0;ii<N;++ii)
    dcpy[ii] = cpy[ii];
}

double Reinit::distToCP(const int idx, const double cpxidx, const double cpyidx) const
{
  const double x = static_cast<double>(idx/m) * dx;
  const double y = static_cast<double>(idx%m) * dy;
  return(myDist(x,y,cpxidx,cpyidx));
}

struct helt Reinit::performDO(const int idx, const struct point cpguess)
{
  if(idx == 231)
    VERBOSITY = true;
  else
    VERBOSITY = false;

  // optimize over directions via line search
  struct point x0;
  x0.x = cpx.getX(idx); x0.y = cpx.getY(idx);
  const int MAXLOOPS = 1;
  const int MAXSEARCHES = 1;
  //const int MAXLOOPS = 5;
  //const int MAXSEARCHES = 5;
  const double TOL = mymax( (0.01 * pow(mymax(dx,dy),flag+1)), 5.0e-16) ;

  struct point xl, xr, xc, xt;
  double dl,dc,dr,dt;
  xc = cpguess;
  if(VERBOSITY)
    printf(" x0 = (%f,%f), xc = (%f,%f).\n",x0.x,x0.y,xc.x,xc.y);
  //dc = dist(x0,xc);
  dc = search1D(idx,xc);
  // set angle increment
  double delta;
  if(dc > mymax(dx,dy))
    delta = asin(mymax(dx,dy)/dc);
  else // interface is nearby, look over a wide range
    delta = PI2;

  findNborDirections(x0,xc,xl,xr,delta);
  
  // bracket the interface in the search directions
  dl = search1D(idx,xl);
  dr = search1D(idx,xr);

  if(VERBOSITY)
  {
    mexPrintf(" xl = (%f,%f), xc = (%f,%f), xr = (%f,%f)\n",xl.x,xl.y,xc.x,xc.y,xr.x,xr.y);
    mexPrintf(" dl = %f, dc = %f, dr = %f\n",dl,dc,dr);
  }

  for(int count0=0; count0<MAXLOOPS; ++count0)
  { // (1) center on minimum; (2) perform Newton step
    int count = 0;
    while( (dc > mymin(dl,dr)) && (count < MAXSEARCHES))
    {
      if(delta < PI4-1e-4)
        delta *= 2.0f;

      if(dl < dr)
      {
        dc = dl; xc = xl;
      }
      else // dl >= dr
      {
        dc = dr; xc = xr;
      }
      findNborDirections(x0,xc,xl,xr,delta);
      dl = search1D(idx,xl);
      dr = search1D(idx,xr);
      count++;
    }
    // perform Newton step only if dl/dr are successfully computed
    if(mymax(dl,dr) < 0.9f)
    {
      double ds = (dr-dl)/2.0f;
      double dss = (dr-2.0f*dc+dl);
      if(fabs(dss) < 1e-16)
        dss = 1e-16;
      struct point dir;
      dir.x = pd(xc.x,x0.x); dir.y = pd(xc.y,x0.y);
      double delta2 = -ds/dss * delta;
      xt.x = dinrange2(x0.x + cos(delta2)*dir.x - sin(delta2)*dir.y);
      xt.y = dinrange2(x0.y + sin(delta2)*dir.x + cos(delta2)*dir.y);
      dt = search1D(idx,xt);
      if(fabs(dt-dc) < TOL)
        count0 = MAXLOOPS; // no further improvement available
      else if(dt<dc)
      {
        xc = xt; dc = dt;
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

  lval.i = idx;
  lval.d = dc;
  lval.aux[0] = xc.x;
  lval.aux[1] = xc.y;

  if(VERBOSITY)
    mexPrintf(" After DO, d = %f, cp = (%f,%f)\n",lval.d,lval.aux[0],lval.aux[1]);

  return(lval);
}

void Reinit::findNborDirections(const struct point x0, const struct point xc, struct point &xl, struct point &xr, const double delta) const
{
  double s,c;
  double rx,ry;
  s = sin(delta);
  c = cos(delta);
  rx = pd(xc.x,x0.x);
  ry = pd(xc.y,x0.y);
  xl.x = dinrange2(x0.x + c*rx + s*ry);
  xl.y = dinrange2(x0.y - s*rx + c*ry);
  xr.x = dinrange2(x0.x + c*rx - s*ry);
  xr.y = dinrange2(x0.y + s*rx + c*ry);
}

double Reinit::search1D(const int idx, struct point &x) const
{
  // bracket the interface in the search direction
  // want to find minimizer of h(t) = dist(x0,x0+t*d) + dist(x0+t*d,cp(x0+t*d)), with d = (x-x0)/dist(x-x0)
  struct point x0; x0.x = cpx.getX(idx); x0.y = cpx.getY(idx);
  struct point d; d.x = pd(x.x,x0.x); d.y = pd(x.y,x0.y); 
  double initT = normalize(d);
  // mexPrintf("x = (%f,%f), d = (%f,%f), t = %f, x0 = (%f,%f), x0+t*d = (%f,%f)\n",x.x,x.y,d.x,d.y,initT,x0.x,x0.y,dinrange2(x0.x+initT*d.x),dinrange2(x0.y+initT*d.y));
  double t[2];
  double ht[2];
  double htder[2];
  bool gl = bracket(x0,d,initT,t,ht,htder);

  //if(VERBOSITY)
  //  mexPrintf(" initT = %f\n",initT);

  double htnew;
  if(gl)
    htnew = bisect(x0,d,t,ht,htder,x);
  else
  {
    htnew = 1.0f;
    x.x = 10.0f;
    x.y = 10.0f;
  }

  if(VERBOSITY)
    mexPrintf(" After search1D, htnew = %f, gl = %d\n",htnew,gl);
  return(htnew);
}

void Reinit::computeHtAndDerOld(const struct point x, const struct point d, const double t, double &ht, double &htder) const
{ // h(t) = dist(x,x+t*d) + dist(x+t*d,cp(x+t*d)). Return h(t) in ht and h'(t) in htder

  // xt => x(t) = x + t*d
  struct point xt;
  xt.x = dinrange2(x.x + t*d.x);
  xt.y = dinrange2(x.y + t*d.y);
  double cpxres[3], cpyres[3];
  cpx.interpolate(xt.x,xt.y,cpxres);
  cpy.interpolate(xt.x,xt.y,cpyres);
  // cpres => cp(x(t))
  struct point cpres;
  cpres.x = cpxres[0];
  cpres.y = cpyres[0];

  double d1 = dist(x,xt);
  double d2 = dist(xt,cpres);

  ht = d1+d2;
  if( (d1>0) && (d2>0) )
  {
    struct point derterm;
    derterm.x = cpxres[1]*d.x+cpxres[2]*d.y;
    derterm.y = cpyres[1]*d.x+cpyres[2]*d.y;
    htder = (d.x*pd(xt.x,x.x)+d.y*pd(xt.y,x.y))/dist(x,xt) + ((derterm.x-d.x)*pd(cpres.x,xt.x)+(derterm.y-d.y)*pd(cpres.y,xt.y))/dist(xt,cpres);
  }
  else if( (d1>0) && (d2==0) )
    htder = (d.x*pd(xt.x,x.x)+d.y*pd(xt.y,x.y))/dist(x,xt);
  else if( (d1==0) && (d2>0) )
  {
    struct point derterm;
    derterm.x = cpxres[1]*d.x+cpxres[2]*d.y;
    derterm.y = cpyres[1]*d.x+cpyres[2]*d.y;
    htder = ((derterm.x-d.x)*pd(cpres.x,xt.x)+(derterm.y-d.y)*pd(cpres.y,xt.y))/dist(xt,cpres);
  }
  else
    htder = 0.0f;

  //if(VERBOSITY)
  //  mexPrintf(" ht = %f, htder =%f, d1 = %f, d2 = %f\n",ht,htder,d1,d2);

/*
  // sanity check: PASSED
  double delta = 1.e-4;
  double tt = t + delta;
  xt.x = x.x + tt*d.x;
  xt.y = x.y + tt*d.y;
  double cpxres2[3], cpyres2[3];
  cpx.interpolate(xt.x,xt.y,cpxres2);
  cpy.interpolate(xt.x,xt.y,cpyres2);
  cpres.x = cpxres2[0];
  cpres.y = cpyres2[0];
  double ht2 = dist(x,xt) + dist(xt,cpres);
  double htdercheck = (ht2-ht)/delta;
  mexPrintf("ht = %f, htder = %f, htdercheck = %f, (ht2 = %f)\n",ht,htder,htdercheck,ht2); 
*/
}

void Reinit::computeHtAndDer(const struct point x, const struct point d, const double t, double &ht, double &htder) const
{ // h(t) = dist(x,x+t*d) + dd(x+t*d). Return h(t) in ht and h'(t) in htder

  // xt => x(t) = x + t*d
  struct point xt;
  xt.x = dinrange2(x.x + t*d.x);
  xt.y = dinrange2(x.y + t*d.y);
  double ddres[3];
  dd.interpolate(xt.x,xt.y,ddres);

  double d1 = dist(x,xt);
  double d2 = ddres[0];
  if(d2<0)
    mexPrintf("d2 = %f -- should never be negative!!\n",d2);

  ht = d1+d2;
  if( (d1>0) && (d2>0) )
    htder = (d.x*pd(xt.x,x.x)+d.y*pd(xt.y,x.y))/d1 + ddres[1]*d.x+ddres[2]*d.y;
  else if( (d1>0) && (d2==0) )
    htder = (d.x*pd(xt.x,x.x)+d.y*pd(xt.y,x.y))/d1;
  else if( (d1==0) && (d2>0) )
    htder = ddres[1]*d.x+ddres[2]*d.y;
  else
    htder = 0.0f;
  
  //if(VERBOSITY)
  //  mexPrintf(" ht = %f, htder =%f, d1 = %f, d2 = %f\n",ht,htder,d1,d2);
}

bool Reinit::bracket(const struct point x0, const struct point d, const double initT, double (&t)[2], double (&htout)[2], double (&hderout)[2]) const
{
  const int IMAX = 10;
  const double TOL = 5.e-16;
  double dxy = mymax(dx,dy);
  t[0] = initT;
  computeHtAndDer(x0,d,t[0],htout[0],hderout[0]);
  bool lval = false;
  //if(VERBOSITY)
  //{
  //  mexPrintf(" data = [");
  //  mexPrintf("                            %f %f %f\n",t[0],htout[0],hderout[0]);
  //}
  for(int ii=0; ii<IMAX; ++ii)
  {
    int prev = ii%2;
    int curr = (ii+1)%2;
    if(hderout[prev] < 0.0f) // look further along the line
      t[curr] = t[prev]+dxy;
    else
      t[curr] = t[prev]-dxy;
    computeHtAndDer(x0,d,t[curr],htout[curr],hderout[curr]);

    //if(VERBOSITY)
    //  mexPrintf("                            %f %f %f\n",t[curr],htout[curr],hderout[curr]);

    if(hderout[prev]*hderout[curr] <= 0.0f)
    {
      lval = true;
      break;
    }
    else
      if((htout[prev] < htout[curr]) && (ii != (IMAX-1)))
      {
        htout[curr] = htout[prev];
        hderout[curr] = hderout[prev];
        t[curr] = t[prev];
        dxy /= 2.0f;
      }
  } 
  //if(VERBOSITY)
  //  mexPrintf(" ];\n");
  //mexPrintf("bracket : t = (%f,%f), htout = (%f,%f), hderout = (%f,%f) \n",t[0],t[1],htout[0],htout[1],hderout[0],hderout[1]);
  return(lval);
}

bool Reinit::bracketOld(const struct point x0, const struct point d, const double initT, double (&t)[2], double (&htout)[2], double (&hderout)[2]) const
{
  // Look about line x0+t*d near t=initT to find t's with positive and negative derivatives as input to bisection
  // Returns false if it's much further than initially guessed to minimum in this direction (e.g. more than one grid cell further), and true otherwise.
  // stores values of t, h and dh/dt at t's in t, htout, and hderout, resp.

  const double TOL = 5.e-16;
  const double dxy = mymax(dx,dy);
  double ht[2], htder[2];
  computeHtAndDer(x0,d,0,ht[0],htder[0]);
  if(fabs(htder[0]) < TOL) // x0 is a minimum along this line
  {
    t[0] = 0.0f; t[1] = 0.0f;
    hderout[0] = htder[0]; hderout[1] = htder[0];
    htout[0] = ht[0]; htout[1] = ht[0];
      if(VERBOSITY)
        mexPrintf( " End of bracket: ht = (%f,%f), htder = (%f,%f), t = (%f,%f)\n",ht[0],ht[1],htder[0],htder[1],t[0],t[1]);
    return(true);
  }
  computeHtAndDer(x0,d,initT,ht[1],htder[1]);
  if(fabs(htder[1]) < TOL) // x0+initT*d is a minimum along this line
  {
    t[0] = initT; t[1] = initT;
    hderout[0] = htder[1]; hderout[1] = htder[1];
    htout[0] = ht[1]; htout[1] = ht[1];

      if(VERBOSITY)
        mexPrintf( " End of bracket: ht = (%f,%f), htder = (%f,%f), t = (%f,%f)\n",ht[0],ht[1],htder[0],htder[1],t[0],t[1]);
    return(true);
  }

  // else, both htder[0] and htder[1] have well-signed derivatives
  if(mysign(htder[0] * htder[1]) == 1) // same side of minimum; look in direction of increasing t
  {
    double tt = initT+dxy;
    computeHtAndDer(x0,d,tt,ht[0],htder[0]);
    if(mysign(htder[0] * htder[1]) == 1) // still same side of minimum; fail.
      return(false);
    else
    {
      if(mysign(htder[0]) == 1)
      {
        t[0] = initT; t[1] = tt;
        hderout[0] = htder[1]; hderout[1] = htder[0];
        htout[0] = ht[1]; htout[1] = ht[0];
      }
      else if(mysign(htder[0]) == -1)
      {
        t[0] = tt; t[1] = initT;
        hderout[0] = htder[0]; hderout[1] = htder[1];
        htout[0] = ht[0]; htout[1] = ht[1];
      }
      else // mysign(htder[0]) == 0
      {
        if(mysign(htder[1]) == 1)
        {
          t[0] = tt; t[1] = initT;
          hderout[0] = htder[0]; hderout[1] = htder[1];
          htout[0] = ht[0]; htout[1] = ht[1];
        }
        else
        {
          t[0] = initT; t[1] = tt;
          hderout[0] = htder[1]; hderout[1] = htder[0];
          htout[0] = ht[1]; htout[1] = ht[0];
        }
      }
      
      if(VERBOSITY)
        mexPrintf( " End of bracket: ht = (%f,%f), htder = (%f,%f), t = (%f,%f)\n",ht[0],ht[1],htder[0],htder[1],t[0],t[1]);
      return(true);
    }
  }
  else // htder[0] and htder[1] have opposite signs, look back towards initial point
  {
    int tries = 0;
    double tt = initT-dxy;
    computeHtAndDer(x0,d,tt,ht[0],htder[0]);
    if(mysign(htder[0] * htder[1]) == 1) // still same side of minimum; fail.
    {
      while((mysign(htder[0] * htder[1]) == 1) && (tries < 10))
      {
        tt -= dxy;
        computeHtAndDer(x0,d,tt,ht[0],htder[0]);
        tries++;
      }
      if(mysign(htder[0] * htder[1]) == 1)
        return(false);
      else
      {
        t[0] = initT; t[1] = tt;
        hderout[0] = htder[1]; hderout[1] = htder[0];
        htout[0] = ht[1]; htout[1] = ht[0];
      } // I really don't think that order matters!
    }
    else
    {
      if(mysign(htder[0]) == 1)
      {
        t[0] = initT; t[1] = tt;
        hderout[0] = htder[1]; hderout[1] = htder[0];
        htout[0] = ht[1]; htout[1] = ht[0];
      }
      else if(mysign(htder[0]) == -1)
      {
        t[0] = tt; t[1] = initT;
        hderout[0] = htder[0]; hderout[1] = htder[1];
        htout[0] = ht[0]; htout[1] = ht[1];
      }
      else // mysign(htder[0]) == 0
      {
        if(mysign(htder[1]) == 1)
        {
          t[0] = tt; t[1] = initT;
          hderout[0] = htder[0]; hderout[1] = htder[1];
          htout[0] = ht[0]; htout[1] = ht[1];
        }
        else
        {
          t[0] = initT; t[1] = tt;
          hderout[0] = htder[1]; hderout[1] = htder[0];
          htout[0] = ht[1]; htout[1] = ht[0];
        }
      }
      
      if(VERBOSITY)
        mexPrintf( " End of bracket: ht = (%f,%f), htder = (%f,%f), t = (%f,%f)\n",ht[0],ht[1],htder[0],htder[1],t[0],t[1]);
      return(true);
    }
  }
}

double Reinit::bisect(const struct point &x0, const struct point &d, const double (&t)[2], const double (&ht)[2], const double (&htder)[2], struct point &result) const
{ // bisection with linear approximation
  // Tries to find the minimum of h(t) between t[0] and t[1] by linear
  // interpolation of derivatives. 
  // Returns the value of h (result) and the location on the line (l-value)
  assert(htder[1] >= 0.0f);
  assert(htder[0] <= 0.0f);
  double htnew, htdernew;
  if(htder[0] == htder[1])
  {
    result.x = dinrange2(x0.x + t[0]*d.x);
    result.y = dinrange2(x0.y + t[0]*d.y);
    computeHtAndDer(x0,d,t[0],htnew,htdernew);
    return(htnew);
  }

  double tstar = (t[0]*htder[1]-htder[0]*t[1]) / (htder[1]-htder[0]);
  result.x = dinrange2(x0.x + tstar*d.x);
  result.y = dinrange2(x0.y + tstar*d.y);

  computeHtAndDer(x0,d,tstar,htnew,htdernew);
  struct point x = result;
  result.x = dinrange2(x0.x + tstar*d.x);
  result.y = dinrange2(x0.y + tstar*d.y);

  if(VERBOSITY)
    mexPrintf(" BISECT: t = (%f,%f), ht = (%f,%f), htder = (%f,%f), tstar = %f\n",t[0],t[1],ht[0],ht[1],htder[0],htder[1],tstar);
  if(VERBOSITY)
    mexPrintf("       x = (%f,%f), result = (%f,%f)\n",x.x,x.y,result.x,result.y);
  if(VERBOSITY)
    printf(" After bisect, ht = %f\n",htnew);

  if(ht[0] < htnew)
  {
    htnew = ht[0];
    result.x = dinrange2(x0.x + t[0]*d.x);
    result.y = dinrange2(x0.y + t[0]*d.y);
  }
  if(ht[1] < htnew)
  {
    htnew = ht[1];
    result.x = dinrange2(x0.x + t[1]*d.x);
    result.y = dinrange2(x0.y + t[1]*d.y);
  }
  
  // replace result with cp(result)
  struct point realresult;
  realresult.x = cpx.interpolate(result.x,result.y);
  realresult.y = cpy.interpolate(result.x,result.y);
  result = realresult;
  if(VERBOSITY) 
    mexPrintf(" Final bisect result: ht = %f, res = (%f,%f)\n",htnew,result.x,result.y);
  return(htnew);
}



/*
void Reinit::fastMarchingReinit()
{
  setInterfaceValues();
  thresholdAwayFromInterface(bnd);

  for(size_t ii=0;ii<bnd.size();++ii)
    updateAndAddNeighborsToHeap(bnd[ii]);

  // while heap is non-empty, fix its top value, and update the neighbors
  while(1)
  {
    struct helt h1 = h.popFromHeap();
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

void Reinit::updateAndAddNeighborsToHeap(const int idx)
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
      if(fabs(dtemp) < fabs(u.get(idx2arr[ii])))
      {
        u.put(dtemp,idx2arr[ii]);
        h.addToHeap(idx2arr[ii],fabs(dtemp));
      }
    }
  }
}

void Reinit::updateAndAddNeighborsToHeapDO(const int idx)
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
      struct helt htemp = performDO(idx2arr[ii],cp[0].get(idx),cp[1].get(idx));
      if(fabs(htemp.d) < fabs(u.get(idx2arr[ii])))
      {
        applyResult(htemp);
        h.addToHeap(idx2arr[ii],fabs(htemp.d),htemp.aux);
      }
    }
  }
}

double Reinit::estimateUpdate(const int idx)
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
    return(mysign(u.get(idx)) * (dy*dy*a+dx*dx*b+dx*dy*sqrt(dx*dx+dy*dy-(a-b)*(a-b)))/(dx*dx+dy*dy));
}

void Reinit::setInterfaceValues()
{ 
  Array2D<int> sgn(u.getm(),u.getn());
  for(int ii=0; ii<sgn.getN(); ++ii)
    sgn.put(static_cast<int>(mysign(u.get(ii))),ii);
  for(int ii=0; ii<sgn.getN(); ++ii)
    if((abs(sgn.getxp(ii)-sgn.get(ii)) + abs(sgn.getxm(ii)-sgn.get(ii)) +
        abs(sgn.getyp(ii)-sgn.get(ii)) + abs(sgn.getym(ii)-sgn.get(ii)))>0)
      bnd.push_back(ii);
 
  vector<double> dr(bnd.size());

  for(size_t ii=0; ii<bnd.size(); ++ii)
  {
    // compute norm(grad u) with centered differences
    double rx = (u.getxp(bnd[ii])-u.getxm(bnd[ii]))/dx/2.0f;
    double ry = (u.getyp(bnd[ii])-u.getym(bnd[ii]))/dy/2.0f;
    dr[ii] = sqrt(rx*rx+ry*ry);

    // compute norm(grad u) with one-sided differences
    rx = mymax(fabs(u.getxp(bnd[ii])-u.get(bnd[ii])),fabs(u.get(bnd[ii])-u.getxm(bnd[ii])))/dx;
    ry = mymax(fabs(u.getyp(bnd[ii])-u.get(bnd[ii])),fabs(u.get(bnd[ii])-u.getym(bnd[ii])))/dy;
    const double dr2 = sqrt(rx*rx+ry*ry);

    // Accept one-sided difference is much different than centered difference
    if((dr[ii] < (0.5*dr2)) || (dr[ii] > (2.0*dr2)))
      dr[ii] = dr2;
  }
  for(size_t ii=0; ii<bnd.size();++ii)
    u.put(u.get(bnd[ii])/dr[ii],bnd[ii]);
}

void Reinit::setInterfaceValuesDO()
{ 
  Array2D<int> sgn(m,n);

  double meanabsbndval = 0.0f;

  for(int ii=0; ii<sgn.getN(); ++ii)
    sgn.put(static_cast<int>(mysigntol(u0.get(ii),1e-10)),ii);
  for(int ii=0; ii<sgn.getN(); ++ii)
    if((abs(sgn.getxp(ii)-sgn.get(ii)) + abs(sgn.getxm(ii)-sgn.get(ii)) +
        abs(sgn.getyp(ii)-sgn.get(ii)) + abs(sgn.getym(ii)-sgn.get(ii)))>0)
    {
      bnd.push_back(ii);
      meanabsbndval += fabs(u0.get(ii));
    }

  // rescale u0 based on meanabsbndval
  meanabsbndval /= static_cast<double>(bnd.size());
  if(meanabsbndval > mymax(dx,dy))
    for(int ii=0; ii<u0.getN(); ++ii)
      u0.put(u0.get(ii)/meanabsbndval*mymax(dx,dy),ii);

  vector<struct helt> bndval;
  bndval.resize(bnd.size());
  for(size_t ii=0; ii<bnd.size(); ++ii)
  {
    bndval[ii] = performDO(bnd[ii]);
    // update cp guesses for later interface cells to use
    cp[0].put(bndval[ii].aux[0],bndval[ii].i);
    cp[1].put(bndval[ii].aux[1],bndval[ii].i);
  }

  for(size_t ii=0; ii<bnd.size();++ii)
  {
    applyResult(bndval[ii]);
    state.put(true,bndval[ii].i);
  }
}


Array2D<double>* Reinit::directionalOptimization()
{
  setInterfaceValuesDO(); 
  thresholdAwayFromInterface(bnd);

  for(size_t ii=0;ii<bnd.size();++ii)
    updateAndAddNeighborsToHeapDO(bnd[ii]);

  while(1)
  {
    struct helt h1 = h.popFromHeap();
    if(h1.i == HEAPDONE) // if heap is empty, break from loop
      break; 
    if((state.get(h1.i)==false) && (fabs(h1.d) < thres)) // don't update if outside threshold cutoff, or already fixed
    {
      applyResult(h1);
      cp[0].put(h1.aux[0],h1.i);
      cp[1].put(h1.aux[1],h1.i);
      state.put(true,h1.i);
      updateAndAddNeighborsToHeapDO(h1.i);
    }
  }
  return(cp);
}


struct helt Reinit::performDO(const int idx)
{ // make a guess for interface location and call performDO(const int ind, const double *cpguess)
  double grad[2], cpguess[2];
  grad[0] = pd(u.getxp(idx),u.getxm(idx)) / (2.0f * dx);
  grad[1] = pd(u.getyp(idx),u.getym(idx)) / (2.0f * dy);
  lineSearch(idx,grad,cpguess);
  return(performDO(idx,cpguess[0],cpguess[1]));
}

void Reinit::lineSearch(const int idx, double (&grad)[2], double (&cpguess)[2])
{ // takes the given index at estimate of gradient direction (unnormalized)
  // and performs a line search along this direction for the interface.
  // Interface is expected to lie within one grid cell of idx

  normalize(grad);
  // get the correct direction to search in; 
  if(u.get(idx) < 0.0f) // then we go uphill, keep gradient sign
  {
    grad[0] *= SQRT2*dx;
    grad[1] *= SQRT2*dy;
  }
  else // go downhill
  {
    grad[0] *= -SQRT2*dx;
    grad[1] *= -SQRT2*dy;
  }

  if(findOppSign(idx,grad))
  {
    cpguess[0] = grad[0]; cpguess[1] = grad[1];
    bisect(idx,cpguess);
  }
  else // need to figure out where to look
  {
    const int sgn = mysign(u.get(idx));
    int nborind[8];
    nborind[0] = u.xp(idx);
    nborind[1] = u.xm(idx);
    nborind[2] = u.yp(idx);
    nborind[3] = u.ym(idx);
    nborind[4] = u.yp(u.xp(idx));
    nborind[5] = u.yp(u.xm(idx));
    nborind[6] = u.ym(u.xp(idx));
    nborind[7] = u.ym(u.xm(idx));
    for(int ii=0;ii<8;++ii)
    {
      if(mysign(u.get(nborind[ii])) != sgn)
      {
        cpguess[0] = u.getX(nborind[ii]);
        cpguess[1] = u.getY(nborind[ii]);
        bisect(idx,cpguess);
        return;
      }
    }
  }
}

bool Reinit::findOppSign(const int idx, double (&guess)[2])
{ // searches in the direction initially given by "guess" to find x s.t. u(x) has opposite sign
  // from from u(idx)
  // input suggests to look at u[y+guess], where y is spatial coordinates of idx. 
  // returns location x in guess.
  const int sgn = mysign(u.get(idx));

  double val1 = u0.interpolate(u.getX(idx)+guess[0],u.getY(idx)+guess[1]);
  if(sgn == mysign(val1)) // need to look harder!
  {
    double val2 = u0.interpolate(u.getX(idx)+2.0f*guess[0],u.getY(idx)+2.0f*guess[1]);
    if(sgn == mysign(val2)) // still need to do more work
      return(false);
    else
    {
      guess[0] = u.getX(idx)+2.0f*guess[0];
      guess[1] = u.getY(idx)+2.0f*guess[1];
      return(true);
    }
  }
  else
  {
    guess[0] = u.getX(idx)+guess[0];
    guess[1] = u.getY(idx)+guess[1];
    return(true);
  } 
}

void Reinit::bisect(const int idx, double (&guess)[2])
{ // bisection without Newton
  const int NTRIES = 10;
  // Assumes that u[idx] and u[guess] have different signs.
  // Tries to find the zero between u[idx] and u[guess] by bisection.
  
  double xm[2], xp[2], xt[2];
  double um, up;
  assert(mysign(u.get(idx)*u0.interpolate(guess[0],guess[1])) != 1);

  if(mysign(u.get(idx)) == 0)
  {
    guess[0] = u.getX(idx);
    guess[1] = u.getY(idx);
    return;
  }
  else if(mysign(u.get(idx)) == 1)
  {
    xp[0] = u.getX(idx); xp[1] = u.getY(idx); up = u.get(idx);
    xm[0] = guess[0];    xm[1] = guess[1];    um = u0.interpolate(xm[0],xm[1]);
  }
  else
  {
    xm[0] = u.getX(idx); xm[1] = u.getY(idx); um = u.get(idx);
    xp[0] = guess[0];    xp[1] = guess[1];    up = u0.interpolate(xp[0],xp[1]);

    if(up <= 0.0f)
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
  ccomb(xm,xp,xt,up/(up-um));
  assert(up>=0.0f); 
  assert(um<=0.0f);
  for(int ii=0; ii<NTRIES; ++ii)
  {
    double ut = u0.interpolate(xt[0],xt[1]);
    if(mysign(ut) == 1)
    {
      xp[0] = xt[0]; xp[1] = xt[1]; up = ut;
    }
    else
    {
      xm[0] = xt[0]; xm[1] = xt[1]; um = ut;
    }
    if(isnan(up/(up-um)))
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
  return;
}

void Reinit::secondOrderIterations()
{
  const double me = 1.0f/ 10000.0f * mymax(dx,dy);
  const double thres2 = thres - dx;
  const double dt = 1.0f / 5.0f * mymax(dx,dy);
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
}

bool Reinit::diffSign(int idx)
{ // return true if sign of u[idx] differs from sign of any of its 4 neighbors
  const double ui = u.get(idx);
  if(mysign(ui) != mysign(u.getxp(idx)))
    return(true);
  if(mysign(ui) != mysign(u.getyp(idx)))
    return(true);
  if(mysign(ui) != mysign(u.getxm(idx)))
    return(true);
  if(mysign(ui) != mysign(u.getym(idx)))
    return(true);
  return(false);
}
*/
