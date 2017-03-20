#include "idarray2d.hpp"

IDArray2D::IDArray2D(int const nn, int const _flag) :
  IDArray2D(nn, nn, _flag)
{ }

IDArray2D::IDArray2D(int const mm, int const nn, int const _flag) :
  IDArray2D(mm, nn, 1./static_cast<double>(nn), 1./static_cast<double>(mm), _flag)
{ }

IDArray2D::IDArray2D(int const mm, int const nn, double const _dx, double const _dy, int const _flag) :
  Array2D<double>(mm,nn,_dx,_dy),
  flag(_flag)
{
  qe.assign(N, nullptr);
}

IDArray2D::IDArray2D(const Array2D<double> &input, int const _flag) :
  Array2D<double>(input,0),
  flag(_flag)
{
  qe.resize(N);
  for(idx_t ii=0; ii<N; ++ii)
    qe[ii] = static_cast<double*>(nullptr);
}

IDArray2D::~IDArray2D()
{
  freeQeIdxAll();
}

double IDArray2D::interpolate(double xx, double yy) const
/*
{
  double x = xx-0.5f+dx;
  double y = yy-0.5f+dy;
  return(0.25 - sqrt(x*x+y*y));
}
*/
{
  assert((flag==1) || (flag==2) || (flag==3));
  double x = xx;
  double y = yy;
  if((mymin(x,y) < 0.0f) || (mymax(x,y) >= 1.0f))
  {
    x = dinrange2l(x,lenx());
    y = dinrange2l(y,leny());
  }

  if(flag == 1)
  {
    const idx_t ii = static_cast<idx_t>(floor(y*static_cast<double>(m)));
    const idx_t jj = static_cast<idx_t>(floor(x*static_cast<double>(n)));
    const idx_t idx = inrange(ii,m) + inrange(jj,n)*m;

    const double rx = x*static_cast<double>(n)-static_cast<double>(jj); // in range [-1/2,1/2], not scaled by dx -- cancels with terms in qe matrix   
    const double ry = y*static_cast<double>(m)-static_cast<double>(ii);
    if(qe[idx] == NULL)
    { /* need to compute coefficients */
      qe[idx] = new double [4];
      qe[idx][0] = this->get(idx);
      qe[idx][1] = this->getyp(idx);
      qe[idx][2] = this->getxp(idx);
      qe[idx][3] = this->getyp(this->xp(idx));
    }
    const double xi1 = qe[idx][0]*(1.0-rx) + qe[idx][2]*rx;
    const double xi2 = qe[idx][1]*(1.0-rx) + qe[idx][3]*rx;
    return(xi1*(1.0-ry) + xi2*ry);
  }
  else if(flag == 2)
  {
    idx_t ii = static_cast<idx_t>(round(y/dy));    // [0,m]
    idx_t jj = static_cast<idx_t>(round(x/dx));    // [0,n] 
    idx_t idx = inrange(ii,m) + inrange(jj,n)*m;                     // nearest grid point 
    double rx = x/dx-static_cast<double>(jj);; // in range [-1/2,1/2], not scaled by dx -- cancels with terms in qe matrix 
    double ry = y/dy-static_cast<double>(ii);  // note: this relies on ii being able to take values 0 OR m and only being cast into range for idx 

    if(qe[idx] == NULL)
    { // need to compute coefficients 
      qe[idx] = new double [9];
      double A[9];
      A[0] = this->getxm(this->yp(idx)); A[1] = this->getxm(idx); A[2] = this->getxm(this->ym(idx));
      A[3] = this->getyp(idx); A[4] = this->get(idx); A[5] = this->getym(idx);
      A[6] = this->getxp(this->yp(idx)); A[7] = this->getxp(idx); A[8] = this->getxp(this->ym(idx));
      qe[idx][0] = ((A[0]+A[2]+A[6]+A[8])-2.0*(A[1]+A[3]+A[5]+A[7])+4.0*A[4])/4.0;
      qe[idx][1] = ((A[0]-A[2]+A[6]-A[8])+2.0*(A[5]-A[3]))/4.0;
      qe[idx][2] = (A[1]-2.0*A[4]+A[7])/2.0;
      qe[idx][3] = ((-A[0]-A[2]+A[6]+A[8])+2.0*(A[1]-A[7]))/4.0;
      qe[idx][4] = (-A[0]+A[2]+A[6]-A[8])/4.0;
      qe[idx][5] = (-A[1]+A[7])/2.0;
      qe[idx][6] = (A[3]-2.0*A[4]+A[5])/2.0;
      qe[idx][7] = (A[3]-A[5])/2.0;
      qe[idx][8] = A[4];
    }
    return( qe[idx][0]*rx*rx*ry*ry + qe[idx][1]*rx*rx*ry + qe[idx][2]*rx*rx+
            qe[idx][3]*rx*ry*ry + qe[idx][4]*rx*ry + qe[idx][5]*rx+
            qe[idx][6]*ry*ry + qe[idx][7]*ry + qe[idx][8] );
  }
  else // flag == 3
  {
    idx_t ii = static_cast<idx_t>(floor(y/dy)); 
    idx_t jj = static_cast<idx_t>(floor(x/dx));
    //if((fabs(x-0.734375) < 1.e-3) && (fabs(y-1.093750) < 1.e-3))
    //  mexPrintf(" x = %f, y = %f, dx = %f, dy = %f, ii = %d, jj = %d\n",x,y,dx,dy,ii,jj);

    idx_t idx = inrange(ii,m) + inrange(jj,n)*m;                      // nearest grid point 
    double rx = (x/dx-static_cast<double>(jj)); // in range (-1,1), not scaled by dx/dy (picked up in coefficients) 
    double rx2 = rx*rx; double rx3 = rx2*rx;
    double ry = (y/dy-static_cast<double>(ii)); // in range (-1,1), not scaled by dx/dy (picked up in coefficients) 
    double ry2 = ry*ry; double ry3 = ry2*ry;
    if(qe[idx] == NULL)
    {
      qe[idx] = new double [16];
      double A[16];
      A[ 0] = this->getxm(this->ym(idx));
      A[ 1] = this->getxm(idx);
      A[ 2] = this->getxm(this->yp(idx));
      A[ 3] = this->getxm(this->yp(this->yp(idx)));

      A[ 4] = this->getym(idx);
      A[ 5] = this->get(idx);
      A[ 6] = this->getyp(idx);
      A[ 7] = this->getyp(this->yp(idx));

      A[ 8] = this->getxp(this->ym(idx));
      A[ 9] = this->getxp(idx);
      A[10] = this->getxp(this->yp(idx));
      A[11] = this->getxp(this->yp(this->yp(idx)));

      A[12] = this->getxp(this->xp(this->ym(idx)));
      A[13] = this->getxp(this->xp(idx));
      A[14] = this->getxp(this->xp(this->yp(idx)));
      A[15] = this->getxp(this->xp(this->yp(this->yp(idx))));

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
    return ((qe[idx][0]+qe[idx][1]*ry+qe[idx][2]*ry2+qe[idx][3]*ry3)+rx*(qe[idx][4]+qe[idx][5]*ry+qe[idx][6]*ry2+qe[idx][7]*ry3)+rx2*(qe[idx][8]+qe[idx][9]*ry+qe[idx][10]*ry2+qe[idx][11]*ry3)+rx3*(qe[idx][12]+qe[idx][13]*ry+qe[idx][14]*ry2+qe[idx][15]*ry3));
  }
}

void IDArray2D::interpolate(double xx, double yy, double(&result)[3]) const // returns interpolated value of data(x,y) in result[0] and gradient of data(x,y) in result[1:2]
/*
{ // d = 0.25 - sqrt(x^2 + y^2)
  double x = xx-0.5f+dx;
  double y = yy-0.5f+dy;
  result[0] = 0.25 - sqrt(x*x+y*y);
  result[1] = -x/sqrt(x*x+y*y);
  result[2] = -y/sqrt(x*x+y*y);
}
*/
{
  assert((flag==2) || (flag==3));
  double x = xx;
  double y = yy;
  if((mymin(x,y) < 0.0f) || (mymax(x,y) >= 1.0f))
  {
    x = dinrange2l(x,lenx());
    y = dinrange2l(y,leny());
  }

  if(flag == 2)
  {
    const idx_t ii = static_cast<idx_t>(round(y/dy));    // [0,m]
    const idx_t jj = static_cast<idx_t>(round(x/dx));    // [0,n] 
    const idx_t idx = inrange(ii,m) + inrange(jj,n)*m;                     // nearest grid point 
    const double rx = x/dx-static_cast<double>(jj);; // in range [-1/2,1/2], not scaled by dx -- cancels with terms in qe matrix 
    const double ry = y/dy-static_cast<double>(ii);  // note: this relies on ii being able to take values 0 OR m and only being cast into range for idx 
    const double rx2 = rx*rx;
    const double ry2 = ry*ry;

    if(qe[idx] == NULL)
    { // need to compute coefficients 
      qe[idx] = new double [9];
      double A[9];
      A[0] = this->getxm(this->yp(idx)); A[1] = this->getxm(idx); A[2] = this->getxm(this->ym(idx));
      A[3] = this->getyp(idx); A[4] = this->get(idx); A[5] = this->getym(idx);
      A[6] = this->getxp(this->yp(idx)); A[7] = this->getxp(idx); A[8] = this->getxp(this->ym(idx));
      qe[idx][0] = ((A[0]+A[2]+A[6]+A[8])-2.0*(A[1]+A[3]+A[5]+A[7])+4.0*A[4])/4.0;
      qe[idx][1] = ((A[0]-A[2]+A[6]-A[8])+2.0*(A[5]-A[3]))/4.0;
      qe[idx][2] = (A[1]-2.0*A[4]+A[7])/2.0;
      qe[idx][3] = ((-A[0]-A[2]+A[6]+A[8])+2.0*(A[1]-A[7]))/4.0;
      qe[idx][4] = (-A[0]+A[2]+A[6]-A[8])/4.0;
      qe[idx][5] = (-A[1]+A[7])/2.0;
      qe[idx][6] = (A[3]-2.0*A[4]+A[5])/2.0;
      qe[idx][7] = (A[3]-A[5])/2.0;
      qe[idx][8] = A[4];
    }
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
    const idx_t ii = static_cast<idx_t>(floor(y*static_cast<double>(m))); 
    const idx_t jj = static_cast<idx_t>(floor(x*static_cast<double>(n)));
    const idx_t idx = inrange(ii,m) + inrange(jj,n)*m;                      // nearest grid point 
    const double rx = (x*static_cast<double>(n)-static_cast<double>(jj)); // in range (-1,1), not scaled by dx/dy (picked up in coefficients) 
    const double rx2 = rx*rx; const double rx3 = rx2*rx;
    const double ry = (y*static_cast<double>(m)-static_cast<double>(ii)); // in range (-1,1), not scaled by dx/dy (picked up in coefficients) 
    const double ry2 = ry*ry; const double ry3 = ry2*ry;
    if(qe[idx] == NULL)
    {
      qe[idx] = new double [16];
      double A[16];
      A[ 0] = this->getxm(this->ym(idx));
      A[ 1] = this->getxm(idx);
      A[ 2] = this->getxm(this->yp(idx));
      A[ 3] = this->getxm(this->yp(this->yp(idx)));

      A[ 4] = this->getym(idx);
      A[ 5] = this->get(idx);
      A[ 6] = this->getyp(idx);
      A[ 7] = this->getyp(this->yp(idx));

      A[ 8] = this->getxp(this->ym(idx));
      A[ 9] = this->getxp(idx);
      A[10] = this->getxp(this->yp(idx));
      A[11] = this->getxp(this->yp(this->yp(idx)));

      A[12] = this->getxp(this->xp(this->ym(idx)));
      A[13] = this->getxp(this->xp(idx));
      A[14] = this->getxp(this->xp(this->yp(idx)));
      A[15] = this->getxp(this->xp(this->yp(this->yp(idx))));

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
    result[0] = ((qe[idx][0]+qe[idx][1]*ry+qe[idx][2]*ry2+qe[idx][3]*ry3)+rx*(qe[idx][4]+qe[idx][5]*ry+qe[idx][6]*ry2+qe[idx][7]*ry3)+rx2*(qe[idx][8]+qe[idx][9]*ry+qe[idx][10]*ry2+qe[idx][11]*ry3)+rx3*(qe[idx][12]+qe[idx][13]*ry+qe[idx][14]*ry2+qe[idx][15]*ry3));
    result[1] = ((qe[idx][4]+qe[idx][5]*ry+qe[idx][6]*ry2+qe[idx][7]*ry3)+2.0f*rx*(qe[idx][8]+qe[idx][9]*ry+qe[idx][10]*ry2+qe[idx][11]*ry3)+3.0f*rx2*(qe[idx][12]+qe[idx][13]*ry+qe[idx][14]*ry2+qe[idx][15]*ry3))*dx;
    result[2] = ((qe[idx][ 1]+qe[idx][ 2]*2.0f*ry+qe[idx][ 3]*3.0f*ry2)
             +rx*(qe[idx][ 5]+qe[idx][ 6]*2.0f*ry+qe[idx][ 7]*3.0f*ry2)
            +rx2*(qe[idx][ 9]+qe[idx][10]*2.0f*ry+qe[idx][11]*3.0f*ry2)
            +rx3*(qe[idx][13]+qe[idx][14]*2.0f*ry+qe[idx][15]*3.0f*ry2))*dy;
    return;
  }
}

void IDArray2D::freeQeIdxAll() const
{ // clear all the interpolation coefficients
  for (auto &elt : qe)
    if (elt != nullptr) {
      delete [] elt;
      elt = nullptr;
    }
}
