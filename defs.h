#ifndef _DEFS_H_
#define _DEFS_H_

//#include <vector>
//using std::vector;

/*********** Define Pi *********/
#ifndef PI
#define TWOPI (6.283185307179586)
#define PI (3.14159265358979323846)
#define PI2 (1.57079632679489761923)
#define PI4 (0.785398163397448)
#define SQRT2 (1.41421356237310)
#define SQRT3 (1.73205080756888)
#endif
/*************************************************/

/********** Define maximum and minimum and sign ***********/
#ifndef MAXMINSIGN
  #define MAXMINSIGN
  #define mymax(x,y) ( (x)>(y) ? (x) : (y) )
  #define mymin(x,y) ( (x)>(y) ? (y) : (x) )
  #define mymax3(a,b,c) (mymax(mymax((a),(b)),(c)))
  #define mymin3(a,b,c) (mymin(mymin((a),(b)),(c)))
  #define mysign(x)  ( (x)<(0) ? (-1) : ( (x)==(0) ? (0):(1) ))
  //#define mysign0m(x)  ( (x)>(0) ? (1) : (-1))
  #define mysigntol(x,tol) ( (x)<(-(tol)) ? (-1) : ( (x)>(tol) ? (1):(0) ))
  //#define df(sn,i,n)   ((abs(sn[xp(i,n,n)]-sn[xm(i,n,n)])+abs(sn[yp(i,n,n)]-sn[ym(i,n,n)])) > 0)
  //#define dfMN(sn,i,m,n)   ((abs(sn[xp(i,m,n)]-sn[xm(i,m,n)])+abs(sn[yp(i,m,n)]-sn[ym(i,m,n)])) > 0)
  //#define dfMN2(sn,i,m,n,myep)   ((abs(sn[xp(i,m,n)]-sn[xm(i,m,n)])+abs(sn[yp(i,m,n)]-sn[ym(i,m,n)])+(abs(sn[xp(i,m,n)])<myep)+(abs(sn[xm(i,m,n)])<myep)+(abs(sn[yp(i,m,n)])<myep)+(abs(sn[ym(i,m,n)])<0)) > 0)
  #define pos(x) ( (x)>(0) ? (x) : (0))
#endif
/*************************************************/

/****** Define in-range periodicity ***********/
#ifndef INRANGE
  #define INRANGE
  #define inrange(x,n) ( ((x)<0) ? ((x)+(n)): ((x) < (n) ? (x) : ((x)-(n))) )
  #define dinrange(x) ( ((x)<0) ? ((x)+1) : ((x) < 1 ? (x) : ((x)-1)) )
  #define dinrange2(x) ( (x)-floor(x) )
  #define dinrange2l(x,len) ( (x)-(len)*floor((x)/(len)) )
  #define pd(x,y) ( (fabs((x)-(y)) < 0.5) ? ((x)-(y)) : ( ((x)-(y) >= 0.5) ? ((x)-1.0-(y)) : ((x)+1.0-(y))) )
  #define pdl(x,y,len) ( (fabs((x)-(y)) < len/2.0f) ? ((x)-(y)) : ( ((x)-(y) >= len/2.0f) ? ((x)-len-(y)) : ((x)+len-(y))) )
  #define myDist(x0,y0,x1,y1) ( sqrt(pd((x0),(x1))*pd((x0),(x1)) + pd((y0),(y1))*pd((y0),(y1))) )
  #define myDot(x0,y0,x1,y1) ( (x0)*(x1) + (y0)*(y1) )
#endif
/*************************************************/

/****** Define minmod ****************************/
#ifndef MINMOD
  #define MINMOD
  #define mm(x,y) ( ((x)*(y))>0 ? ( fabs(x) > fabs(y) ? (y):(x)) : (0))
#endif
/*************************************************/

struct point
{
  double x;
  double y;
};

struct point3
{
  double x;
  double y;
  double z;
};

#endif
