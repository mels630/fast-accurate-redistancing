#include "idarray3d.hpp"

IDArray3D::IDArray3D(const int nn, const int _flag) :
  Array3D<double>(nn),
  flag(_flag)
{
  qe.assign(N, nullptr);
}

IDArray3D::IDArray3D(const int mm, const int nn, const int kk, const int _flag) :
  Array3D<double>(mm,nn,kk),
  flag(_flag)
{
  qe.assign(N, nullptr);
}

IDArray3D::IDArray3D(const Array3D<double> &input, const int _flag) :
  Array3D<double>(input,0),
  flag(_flag)
{
  qe.assign(N, nullptr);
}

IDArray3D::~IDArray3D()
{
  for(size_t ii=0; ii<N; ++ii)
    if(qe[ii] != nullptr) {
      delete [] qe[ii];
      qe[ii] = nullptr;
    }
}

double IDArray3D::interpolate(double x, double y, double z)
{
  if((mymin(mymin(x,y),z) < 0.0f) || (mymax(mymax(x,y),z) >= 1.0f))
  {
    x = dinrange2(x);
    y = dinrange2(y);
    z = dinrange2(z);
  }

  if(flag == 2)
  {
    int ii = static_cast<int>(round(y*static_cast<double>(m)));    // [0,m]
    int jj = static_cast<int>(round(x*static_cast<double>(n)));    // [0,n] 
    int kk = static_cast<int>(round(z*static_cast<double>(k)));    // [0,k] 
    int idx = inrange(ii, int(m)) + inrange(jj, int(n))*m + inrange(kk, int(k))*m*n; // nearest grid point 

    const double rx = x*static_cast<double>(n)-static_cast<double>(jj); // in range [-1/2,1/2], not scaled by dx -- cancels with terms in qe matrix   
    const double ry = y*static_cast<double>(m)-static_cast<double>(ii);
    const double rz = z*static_cast<double>(k)-static_cast<double>(kk);
    const double rx2 = rx*rx; 
    const double ry2 = ry*ry;
    const double rz2 = rz*rz;

    if(qe[idx] == NULL)
    { /* need to compute coefficients */
      qe[idx] = new double [27];
      double A[27];

      A[ 0] = this->getzm(this->xm(this->ym(idx)));
      A[ 1] = this->getzm(this->xm(         idx ));
      A[ 2] = this->getzm(this->xm(this->yp(idx)));
      A[ 3] = this->getzm(         this->ym(idx) );
      A[ 4] = this->getzm(                  idx  );
      A[ 5] = this->getzm(         this->yp(idx) );
      A[ 6] = this->getzm(this->xp(this->ym(idx)));
      A[ 7] = this->getzm(this->xp(         idx ));
      A[ 8] = this->getzm(this->xp(this->yp(idx)));
      A[ 9] =          this->getxm(this->ym(idx)) ;
      A[10] =          this->getxm(         idx ) ;
      A[11] =          this->getxm(this->yp(idx)) ;
      A[12] =                   this->getym(idx)  ;
      A[13] =                     this->get(idx)  ;
      A[14] =                   this->getyp(idx)  ;
      A[15] =          this->getxp(this->ym(idx)) ;
      A[16] =          this->getxp(   idx       ) ;
      A[17] =          this->getxp(this->yp(idx)) ;
      A[18] = this->getzp(this->xm(this->ym(idx)));
      A[19] = this->getzp(this->xm(         idx ));
      A[20] = this->getzp(this->xm(this->yp(idx)));
      A[21] = this->getzp(         this->ym(idx) );
      A[22] = this->getzp(                  idx  );
      A[23] = this->getzp(         this->yp(idx) );
      A[24] = this->getzp(this->xp(this->ym(idx)));
      A[25] = this->getzp(this->xp(         idx ));
      A[26] = this->getzp(this->xp(this->yp(idx)));

      qe[idx][ 0] = A[13];

      qe[idx][ 1] = (A[14]-A[12])/2.0;
      qe[idx][ 3] = (A[16]-A[10])/2.0;
      qe[idx][ 9] = (A[22]-A[ 4])/2.0;

      qe[idx][ 2] = (A[14]+A[12])/2.0 - A[13];
      qe[idx][ 6] = (A[16]+A[10])/2.0 - A[13];
      qe[idx][18] = (A[22]+A[ 4])/2.0 - A[13];

      qe[idx][ 4] = (A[17]-A[15]-A[11]+A[ 9])/4.0;
      qe[idx][10] = (A[23]-A[21]-A[ 5]+A[ 3])/4.0;
      qe[idx][12] = (A[25]-A[19]-A[ 7]+A[ 1])/4.0;

      qe[idx][ 7] = (A[11]+A[17]-A[ 9]-A[15])/4.0 + (A[12]-A[14])/2.0;
      qe[idx][15] = (A[19]+A[25]-A[ 1]-A[ 7])/4.0 + (A[ 4]-A[22])/2.0;
      qe[idx][ 5] = (A[17]+A[15]-A[11]-A[ 9])/4.0 + (A[10]-A[16])/2.0;
      qe[idx][11] = (A[23]+A[21]-A[ 5]-A[ 3])/4.0 + (A[ 4]-A[22])/2.0;
      qe[idx][21] = (A[25]+A[ 7]-A[19]-A[ 1])/4.0 + (A[10]-A[16])/2.0;
      qe[idx][19] = (A[23]+A[ 5]-A[21]-A[ 3])/4.0 + (A[12]-A[14])/2.0;
      
      qe[idx][13] = (A[26]-A[ 8]-A[24]+A[ 6]-A[20]+A[ 2]+A[18]-A[ 0])/8.0;

      qe[idx][ 8] = (A[ 9]+A[15]+A[11]+A[17])/4.0 - (A[12]+A[10]+A[16]+A[14])/2.0 + A[13];
      qe[idx][24] = (A[ 1]+A[ 7]+A[19]+A[25])/4.0 - (A[ 4]+A[10]+A[16]+A[22])/2.0 + A[13];
      qe[idx][20] = (A[ 3]+A[ 5]+A[21]+A[23])/4.0 - (A[ 4]+A[12]+A[14]+A[22])/2.0 + A[13];
      
      qe[idx][14] = (A[26]-A[ 8]+A[24]-A[ 6]-A[20]+A[ 2]-A[18]+A[ 0])/8.0 + (A[ 7]-A[25]+A[19]-A[ 1])/4.0;
      qe[idx][16] = (A[20]-A[ 2]+A[26]-A[ 8]-A[18]+A[ 0]-A[24]+A[ 6])/8.0 + (A[ 5]-A[23]+A[21]-A[ 3])/4.0;
      qe[idx][22] = (A[26]-A[24]+A[ 8]-A[ 6]-A[20]+A[18]-A[ 2]+A[ 0])/8.0 + (A[15]-A[17]+A[11]-A[ 9])/4.0;
      
      qe[idx][23] = (A[ 6]-A[ 0]+A[ 8]-A[ 2]+A[24]-A[18]+A[26]-A[20])/8.0 - (A[ 7]-A[ 1]+A[15]-A[ 9]+A[17]-A[11]+A[25]-A[19])/4.0 + (A[16]-A[10])/2.0;
      qe[idx][25] = (A[ 2]-A[ 0]+A[ 8]-A[ 6]+A[20]-A[18]+A[26]-A[24])/8.0 - (A[ 5]-A[ 3]+A[11]-A[ 9]+A[17]-A[15]+A[23]-A[21])/4.0 + (A[14]-A[12])/2.0;
      qe[idx][17] = (A[18]-A[ 0]+A[24]-A[ 6]+A[20]-A[ 2]+A[26]-A[ 8])/8.0 - (A[21]-A[ 3]+A[19]-A[ 1]+A[25]-A[ 7]+A[23]-A[ 5])/4.0 + (A[22]-A[ 4])/2.0;

      qe[idx][26] = (A[ 0]+A[ 6]+A[ 2]+A[ 8]+A[18]+A[24]+A[20]+A[26])/8.0 - (A[ 1]+A[ 3]+A[ 5]+A[ 7]+A[ 9]+A[11]+A[15]+A[17]+A[19]+A[21]+A[23]+A[25])/4.0 + (A[ 4]+A[10]+A[12]+A[14]+A[16]+A[22])/2.0-A[13];
    }
    return( ((qe[idx][ 0]+qe[idx][ 1]*ry+qe[idx][ 2]*ry2) + rx*(qe[idx][ 3]+qe[idx][ 4]*ry+qe[idx][ 5]*ry2) + rx2*(qe[idx][ 6]+qe[idx][ 7]*ry+qe[idx][ 8]*ry2)) +
            rz *((qe[idx][ 9]+qe[idx][10]*ry+qe[idx][11]*ry2) + rx*(qe[idx][12]+qe[idx][13]*ry+qe[idx][14]*ry2) + rx2*(qe[idx][15]+qe[idx][16]*ry+qe[idx][17]*ry2)) +
            rz2*((qe[idx][18]+qe[idx][19]*ry+qe[idx][20]*ry2) + rx*(qe[idx][21]+qe[idx][22]*ry+qe[idx][23]*ry2) + rx2*(qe[idx][24]+qe[idx][25]*ry+qe[idx][26]*ry2)) );
  }
  else if(flag == 1)
  { 
    size_t ii = static_cast<size_t>(floor(y*static_cast<double>(m)));    // [0,m)
    size_t jj = static_cast<size_t>(floor(x*static_cast<double>(n)));    // [0,n) 
    size_t kk = static_cast<size_t>(floor(z*static_cast<double>(k)));    // [0,k)
    size_t idx = inrange(ii,m) + inrange(jj,n)*m + inrange(kk,k)*m*n; // nearest grid point 

    const double rx = x*static_cast<double>(n)-static_cast<double>(jj); // in range [-1/2,1/2], not scaled by dx -- cancels with terms in qe matrix   
    const double ry = y*static_cast<double>(m)-static_cast<double>(ii);
    const double rz = z*static_cast<double>(k)-static_cast<double>(kk);

    if(qe[idx] == NULL)
    { /* need to compute coefficients */
      qe[idx] = new double [8];
      qe[idx][0] = this->get(idx);
      qe[idx][1] = this->getyp(idx);
      qe[idx][2] = this->getxp(idx);
      qe[idx][3] = this->getyp(this->xp(idx));
      qe[idx][4] = this->getzp(idx);
      qe[idx][5] = this->getzp(this->yp(idx));
      qe[idx][6] = this->getzp(this->xp(idx));
      qe[idx][7] = this->getzp(this->xp(this->yp(idx)));
    }
    double xi1 = qe[idx][0]*(1.0-rx) + qe[idx][2]*rx;
    double xi2 = qe[idx][1]*(1.0-rx) + qe[idx][3]*rx;
    double xi3 = qe[idx][4]*(1.0-rx) + qe[idx][6]*rx;
    double xi4 = qe[idx][5]*(1.0-rx) + qe[idx][7]*rx;
    double yi1 = xi1*(1.0-ry) + xi2*ry;
    double yi2 = xi3*(1.0-ry) + xi4*ry;
    return(yi1*(1.0-rz) + yi2*rz);
  }
  else // flag == 3
  {
    cout << "Tricubic interpolation is not yet implemented. Care to have a go? ... Aborting." << endl;
    abort();
  }
}
