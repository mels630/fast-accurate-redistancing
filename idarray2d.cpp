#include "idarray2d.hpp"

/// Constructor
/// \param[in] nn    : Number of elements per row and column
/// \param[in] _flag : Interpolation flag
IDArray2D::IDArray2D(int const nn, int const _flag) :
  IDArray2D(nn, nn, _flag)
{ }

/// Constructor
/// \param[in] mm    : Number of elements per row
/// \param[in] nn    : Number of elements per column
/// \param[in] _flag : Interpolation flag
IDArray2D::IDArray2D(int const mm, int const nn, int const _flag) :
  IDArray2D(mm, nn, 1./static_cast<double>(nn), 1./static_cast<double>(mm), _flag)
{ }

/// Constructor
/// \param[in] mm    : Number of elements per row
/// \param[in] nn    : Number of elements per column  
/// \param[in] _dx   : Element spacing in x
/// \param[in] _dy   : Element spacing in y
/// \param[in] _flag : Interpolation flag
IDArray2D::IDArray2D(int const mm, int const nn, double const _dx, double const _dy, int const _flag) :
  Array2D<double>(mm,nn,_dx,_dy),
  qe(),
  flag(_flag)
{
  qe.assign(N, nullptr);
}

/// Copy constructor
/// \param[in] input : Array2D<T> object to copy
/// \param[in] _flag : Interpolation flag
IDArray2D::IDArray2D(const Array2D<double> &input, int const _flag) :
  Array2D<double>(input,0),
  qe(),
  flag(_flag)
{
  qe.assign(N, nullptr);
}

/// Destructor
IDArray2D::~IDArray2D()
{ // clear all the interpolation coefficients
  for (auto &elt : qe)
    if (elt != nullptr) {
      delete [] elt;
      elt = nullptr;
    }
}

/// Compute and store interpolation coefficients
/// \param[in] idx : Flat index to create interpolation coefficients at
void IDArray2D::SetQe(idx_t const idx) const
{
  if (qe[idx] != nullptr)
    return;
  if (flag == 1) {
    qe[idx] = new double [4];
    qe[idx][0] = get(idx);
    qe[idx][1] = getyp(idx);
    qe[idx][2] = getxp(idx);
    qe[idx][3] = getyp(xp(idx));
  } else if (flag == 2) {
    qe[idx] = new double [9];
    double A[9];
    A[0] = getxm(yp(idx)); A[1] = getxm(idx); A[2] = getxm(ym(idx));
    A[3] = getyp(idx); A[4] = get(idx); A[5] = getym(idx);
    A[6] = getxp(yp(idx)); A[7] = getxp(idx); A[8] = getxp(ym(idx));
    qe[idx][0] = ((A[0]+A[2]+A[6]+A[8])-2.0*(A[1]+A[3]+A[5]+A[7])+4.0*A[4])/4.0;
    qe[idx][1] = ((A[0]-A[2]+A[6]-A[8])+2.0*(A[5]-A[3]))/4.0;
    qe[idx][2] = (A[1]-2.0*A[4]+A[7])/2.0;
    qe[idx][3] = ((-A[0]-A[2]+A[6]+A[8])+2.0*(A[1]-A[7]))/4.0;
    qe[idx][4] = (-A[0]+A[2]+A[6]-A[8])/4.0;
    qe[idx][5] = (-A[1]+A[7])/2.0;
    qe[idx][6] = (A[3]-2.0*A[4]+A[5])/2.0;
    qe[idx][7] = (A[3]-A[5])/2.0;
    qe[idx][8] = A[4];
  } else { // flag == 3
    qe[idx] = new double [16];
    double A[16];
    A[ 0] = getxm(ym(idx));
    A[ 1] = getxm(idx);
    A[ 2] = getxm(yp(idx));
    A[ 3] = getxm(yp(yp(idx)));

    A[ 4] = getym(idx);
    A[ 5] = get(idx);
    A[ 6] = getyp(idx);
    A[ 7] = getyp(yp(idx));

    A[ 8] = getxp(ym(idx));
    A[ 9] = getxp(idx);
    A[10] = getxp(yp(idx));
    A[11] = getxp(yp(yp(idx)));

    A[12] = getxp(xp(ym(idx)));
    A[13] = getxp(xp(idx));
    A[14] = getxp(xp(yp(idx)));
    A[15] = getxp(xp(yp(yp(idx))));

    qe[idx][0] = A[5];
    qe[idx][1] = -A[7]/6.0+A[6]-A[5]/2.0-A[4]/3.0;
    qe[idx][2] = (A[6]+A[4])/2.0 - A[5];
    qe[idx][3] = (A[7]-A[4])/6.0 - (A[6]-A[5])/2.0;
    qe[idx][4] = -A[13]/6.0+A[9]-A[5]/2.0-A[1]/3.0;
    qe[idx][5] = A[15]/36.0 + (A[12]+A[3])/18.0 + (A[13]+A[7])/12.0 + A[0]/9.0 + (-A[14]-A[11]+A[4]+A[1])/6.0 + A[5]/4.0 - (A[8]+A[2])/3.0 - (A[9]+A[6])/2.0 + A[10];
    qe[idx][6] = -(A[14]+A[12])/12.0 - (A[2]-A[13]+A[0])/6.0 - (A[6]+A[4])/4.0 + A[1]/3.0 + (A[10]+A[5]+A[8])/2.0 - A[9]; 
    qe[idx][7] = (-A[15]+A[12])/36.0 + (-A[3]+A[0])/18.0 + (-A[7]+A[14]-A[13]+A[4])/12.0 + (A[11]+A[2]-A[1]-A[8])/6.0 + (A[6]-A[5])/4.0 + (-A[10]+A[9])/2.0;
    qe[idx][8] = (A[9]+A[1])/2.0 - A[5];
    qe[idx][9] = -(A[11]+A[3])/12.0 - (A[8]-A[7]+A[0])/6.0 - (A[9]+A[1])/4.0 + A[4]/3.0 + (A[10]+A[5]+A[2])/2.0 - A[6];
    qe[idx][10] = (A[10]+A[8]+A[2]+A[0])/4.0 - (A[9]+A[6]+A[4]+A[1])/2.0 + A[5];
    qe[idx][11] = (A[11]-A[8]+A[3]-A[0])/12.0 - (A[10]-A[9]+A[2]-A[1])/4.0 - (A[7]-A[4])/6.0 + (A[6]-A[5])/2.0;
    qe[idx][12] = (A[13]-A[1])/6.0 - (A[9]-A[5])/2.0;
    qe[idx][13] = (-A[15]+A[3])/36.0 + (-A[12]+A[0])/18.0 + (-A[13]+A[11]-A[7]+A[1])/12.0 + (A[14]+A[8]-A[4]-A[2])/6.0 + (A[9]-A[5])/4.0 + (-A[10]+A[6])/2.0;
    qe[idx][14] = (A[14]-A[2]+A[12]-A[0])/12.0 - (A[10]-A[6]+A[8]-A[4])/4.0 - (A[13]-A[1])/6.0 + (A[9]-A[5])/2.0;
    qe[idx][15] = (A[15]-A[12]-A[3]+A[0])/36.0 - (A[14]-A[13]+A[11]-A[8]-A[7]+A[4]-A[2]+A[1])/12.0 + (A[10]-A[9]-A[6]+A[5])/4.0;
  }
}

/// Compute interpolated value at (x,y) ([0,0] element is origin)
/// \param[in] x : x location to interpolate at
/// \param[in] y : y location to interpolate at
/// \return        Interpolated value
double IDArray2D::interpolate(double x, double y) const
{
  assert((flag==1) || (flag==2) || (flag==3));
  if((x < 0.) || (x >= lenx()))
    x = dinrange2l(x,lenx());
  if((y < 0.) || (y >= leny()))
    y = dinrange2l(y,leny());

  if(flag == 1)
  {
    idx_t const ii = static_cast<idx_t>(floor(y*static_cast<double>(m)));
    idx_t const jj = static_cast<idx_t>(floor(x*static_cast<double>(n)));
    idx_t const idx = inrange(ii,m) + inrange(jj,n)*m;

    if (qe[idx] == nullptr)
      SetQe(idx);
    double const rx = x*static_cast<double>(n)-static_cast<double>(jj); // in range [-1/2,1/2], not scaled by dx -- cancels with terms in qe matrix   
    double const ry = y*static_cast<double>(m)-static_cast<double>(ii);
    double const xi1 = qe[idx][0]*(1.0-rx) + qe[idx][2]*rx;
    double const xi2 = qe[idx][1]*(1.0-rx) + qe[idx][3]*rx;
    return(xi1*(1.0-ry) + xi2*ry);
  }
  else if(flag == 2)
  {
    idx_t const ii = static_cast<idx_t>(round(y/dy));    // [0,m]
    idx_t const jj = static_cast<idx_t>(round(x/dx));    // [0,n] 
    idx_t const idx = inrange(ii,m) + inrange(jj,n)*m;                     // nearest grid point 
    double const rx = x/dx-static_cast<double>(jj);; // in range [-1/2,1/2], not scaled by dx -- cancels with terms in qe matrix 
    double const ry = y/dy-static_cast<double>(ii);  // note: this relies on ii being able to take values 0 OR m and only being cast into range for idx 
    if (qe[idx] == nullptr)
      SetQe(idx);
    return( qe[idx][0]*rx*rx*ry*ry + qe[idx][1]*rx*rx*ry + qe[idx][2]*rx*rx+
            qe[idx][3]*rx*ry*ry + qe[idx][4]*rx*ry + qe[idx][5]*rx+
            qe[idx][6]*ry*ry + qe[idx][7]*ry + qe[idx][8] );
  }
  else // flag == 3
  {
    idx_t const ii = static_cast<idx_t>(floor(y/dy)); 
    idx_t const jj = static_cast<idx_t>(floor(x/dx));
    idx_t const idx = inrange(ii,m) + inrange(jj,n)*m;                      // nearest grid point 
    double const rx = (x/dx-static_cast<double>(jj)); // in range (-1,1), not scaled by dx/dy (picked up in coefficients) 
    double const rx2 = rx*rx;
    double const rx3 = rx2*rx;
    double const ry = (y/dy-static_cast<double>(ii)); // in range (-1,1), not scaled by dx/dy (picked up in coefficients) 
    double const ry2 = ry*ry;
    double const ry3 = ry2*ry;
    if (qe[idx] == nullptr)
      SetQe(idx);
    return ((qe[idx][0]+qe[idx][1]*ry+qe[idx][2]*ry2+qe[idx][3]*ry3)+rx*(qe[idx][4]+qe[idx][5]*ry+qe[idx][6]*ry2+qe[idx][7]*ry3)+rx2*(qe[idx][8]+qe[idx][9]*ry+qe[idx][10]*ry2+qe[idx][11]*ry3)+rx3*(qe[idx][12]+qe[idx][13]*ry+qe[idx][14]*ry2+qe[idx][15]*ry3));
  }
}

/// Compute interpolated value at (x,y) ([0,0] element is origin)
/// \param[in]  x      : x location to interpolate at
/// \param[in]  y      : y location to interpolate at
/// \param[out] result : [interpolated value, (gradient of data at x,y)_x (gradient of data at x,y)_y]
void IDArray2D::interpolate(double x, double y, double(&result)[3]) const // returns interpolated value of data(x,y) in result[0] and gradient of data(x,y) in result[1:2]
{
  assert((flag==2) || (flag==3));
  if((x < 0.) || (x >= lenx()))
    x = dinrange2l(x,lenx());
  if((y < 0.) || (y >= leny()))
    y = dinrange2l(y,leny());

  if(flag == 2)
  {
    idx_t const ii = static_cast<idx_t>(round(y/dy));    // [0,m]
    idx_t const jj = static_cast<idx_t>(round(x/dx));    // [0,n] 
    idx_t const idx = inrange(ii,m) + inrange(jj,n)*m;                     // nearest grid point 
    double const rx = x/dx-static_cast<double>(jj);; // in range [-1/2,1/2], not scaled by dx -- cancels with terms in qe matrix 
    double const ry = y/dy-static_cast<double>(ii);  // note: this relies on ii being able to take values 0 OR m and only being cast into range for idx 
    double const rx2 = rx*rx;
    double const ry2 = ry*ry;
    if (qe[idx] == nullptr)
      SetQe(idx);

    result[0] = ( qe[idx][0]*rx2*ry2 + qe[idx][1]*rx2*ry + qe[idx][2]*rx2+
                  qe[idx][3]*rx*ry2 + qe[idx][4]*rx*ry + qe[idx][5]*rx+
                  qe[idx][6]*ry2 + qe[idx][7]*ry + qe[idx][8] );
    result[1] = ( qe[idx][0]*2.0f*rx*ry2 + qe[idx][1]*2.0f*rx*ry + qe[idx][2]*2.0f*rx+
                  qe[idx][3]*ry2 + qe[idx][4]*ry + qe[idx][5] ) * static_cast<double>(n);
    result[2] = ( qe[idx][0]*rx2*2.0f*ry + qe[idx][1]*rx2 +
                  qe[idx][3]*rx*2.0f*ry + qe[idx][4]*rx +
                  qe[idx][6]*2.0f*ry + qe[idx][7] )*static_cast<double>(m);
    return;
  }
  else // flag == 3
  {
    idx_t const ii = static_cast<idx_t>(floor(y*static_cast<double>(m))); 
    idx_t const jj = static_cast<idx_t>(floor(x*static_cast<double>(n)));
    idx_t const idx = inrange(ii,m) + inrange(jj,n)*m;                      // nearest grid point 
    double const rx = (x*static_cast<double>(n)-static_cast<double>(jj)); // in range (-1,1), not scaled by dx/dy (picked up in coefficients) 
    double const rx2 = rx*rx; double const rx3 = rx2*rx;
    double const ry = (y*static_cast<double>(m)-static_cast<double>(ii)); // in range (-1,1), not scaled by dx/dy (picked up in coefficients) 
    double const ry2 = ry*ry; double const ry3 = ry2*ry;
    if (qe[idx] == nullptr)
      SetQe(idx);
    result[0] = ((qe[idx][0]+qe[idx][1]*ry+qe[idx][2]*ry2+qe[idx][3]*ry3)+rx*(qe[idx][4]+qe[idx][5]*ry+qe[idx][6]*ry2+qe[idx][7]*ry3)+rx2*(qe[idx][8]+qe[idx][9]*ry+qe[idx][10]*ry2+qe[idx][11]*ry3)+rx3*(qe[idx][12]+qe[idx][13]*ry+qe[idx][14]*ry2+qe[idx][15]*ry3));
    result[1] = ((qe[idx][4]+qe[idx][5]*ry+qe[idx][6]*ry2+qe[idx][7]*ry3)+2.0f*rx*(qe[idx][8]+qe[idx][9]*ry+qe[idx][10]*ry2+qe[idx][11]*ry3)+3.0f*rx2*(qe[idx][12]+qe[idx][13]*ry+qe[idx][14]*ry2+qe[idx][15]*ry3))*dx;
    result[2] = ((qe[idx][ 1]+qe[idx][ 2]*2.0f*ry+qe[idx][ 3]*3.0f*ry2)
             +rx*(qe[idx][ 5]+qe[idx][ 6]*2.0f*ry+qe[idx][ 7]*3.0f*ry2)
            +rx2*(qe[idx][ 9]+qe[idx][10]*2.0f*ry+qe[idx][11]*3.0f*ry2)
            +rx3*(qe[idx][13]+qe[idx][14]*2.0f*ry+qe[idx][15]*3.0f*ry2))*dy;
    return;
  }
}

