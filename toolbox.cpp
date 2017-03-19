#include "toolbox.hpp"

bool fileExist(const string &filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if(infile.good())
  {
    infile.close();
    infile.clear();
    return(true);
  }
  else
    return(false);
}

void readToColon(ifstream &infile)
{
  char thischar = 'x';
  while(thischar != ':')
  {
    if(!infile.eof())
      infile >> thischar;
    else
    {
      cout << "Read to end-of-file. Aborting." << endl;
      abort();
    }
  }
  return;
}

void removeTrailingSlash(string &mystring)
{ // remove null characters AND trailing slash
  for(size_t ii=0; ii<mystring.length(); ++ii)
    while((ii < mystring.length()) && (mystring[ii] == '\0'))
      mystring.erase(ii);
  if(mystring[mystring.length()-1] == '/')
    mystring.erase(mystring.length()-1);
  //cout << " mystring.length = " << mystring.length() << endl;
  return;
}

void myMkdir(const string &dirname)
{
  DIR *d1;
  d1 = opendir(dirname.c_str());
  if(d1 != NULL)
    cout << " Warning, directory " << dirname << " already exists." << endl;
  else
    if(mkdir(dirname.c_str(),S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == -1)
    {
      cout << "Failed to create directory " << dirname << ", aborting..." << endl;
      abort();
    }
  free(d1);
  return;
}

void printn(int n,int newline)
{
  string tmp;
  tmp.assign(n,'*');
  cout << tmp;
  if(newline)
    cout << endl;
  return;
}

void printtextn(string str, int n, int newline)
{
  int midlen = str.length();
  int remlen = n-midlen;
  printn(remlen/2,0);
  cout << str;
  printn(remlen-remlen/2,newline);
  return;
}

void printFileName(string filename)
{
  int linelength = 120;
  string mystr;
  printn(linelength,1);
  mystr = " Running program " + filename + ". ";
  printtextn(mystr,linelength,1);
  printn(linelength,1);
  return;
}

void fillSign(Array2D<double> &u, const double nullval, const double newval)
{
  const double TOL = 0.0000000001;
  size_t curr = 0;
  vector<int> ind;
  ind.reserve(u.getN());

  for(idx_t ii=0; ii<u.getN(); ++ii)
    if(fabs(u.get(ii)-nullval) > TOL)
      ind.push_back(ii);

  // Fill in everywhere that u == nullval with the appropriate signed value
  // Loop over the set values of u until the end of the list ind is reached

  while(curr < ind.size())
  {
    fillSignSubfn(u,nullval,ind,curr,newval,TOL);
    curr++;
  }
}

void fillSignSubfn(Array2D<double> &u, const double nullval, vector<int> &list, const int curr, const double newval, const double TOL)
{
  if(fabs(u.getxp(list[curr]) - nullval) < TOL)
  { 
    u.put(newval * (double)mysign(u.get(list[curr])), u.xp(list[curr]));
    list.push_back(u.xp(list[curr]));
  }
  if(fabs(u.getxm(list[curr]) - nullval) < TOL)
  {  
    u.put(newval * (double)mysign(u.get(list[curr])), u.xm(list[curr]));
    list.push_back(u.xm(list[curr]));
  }
  if(fabs(u.getyp(list[curr]) - nullval) < TOL)
  {  
    u.put(newval * (double)mysign(u.get(list[curr])), u.yp(list[curr]));
    list.push_back(u.yp(list[curr]));
  }
  if(fabs(u.getym(list[curr]) - nullval) < TOL)
  {  
    u.put(newval * (double)mysign(u.get(list[curr])), u.ym(list[curr]));
    list.push_back(u.ym(list[curr]));
  }
}

vector< vector<int> > createLabelList(const Array2D<int> &q)
{ // q is a 2D array of integers, with values in [0,numlabel].
  // Values of 0 are ignored, and a vector of SLists
  // is returned, with lval[ii] containing the indices
  // corresponding to the label (ii+1)

  assert(q.minval() >= 0);

  int numlabel = q.maxval();
  vector< vector<int> > lval;
  lval.resize(numlabel);
  for(idx_t ii=0; ii < q.getN(); ++ii)
    if(q.get(ii) > 0)
      lval[q.get(ii)-1].push_back(ii); // check to see if reserve if
                                       // worthwhile here
  return(lval);
}

Array2D<int> expandList(const vector<int> &ind, const int m, const int n)
{ // expands a list of indices 
  Array2D<int> q(m,n);
  q.fillWithValue(0);
  for(size_t ii=0; ii<ind.size(); ++ii)
    q.put(1,ind[ii]);
  return(q);
}

vector<int> binaryGrow(const vector<int> &idx, const int m, const int n, const int num)
{ // Dilates a given simple list of indices by num cells 
  //  (using octagonal approximation)
  // Inputs:
  //   idx: list of indices to be dilated
  //   (m,n): grid dimensions for dilation
  //   num: number of pixels to grow in each direction
  // Output:
  //   list of indices in dilated version of the grain

  vector<int> lval, gl;
  vector<int> candidate(idx);
  Array2D<int> q(m,n);
  q.fillWithValue(0);

  for(size_t ii=0; ii<candidate.size(); ++ii)
  {
    q.put(1,candidate[ii]);
    gl.push_back(candidate[ii]);
  }
  // find boundary of the given component
  for(int ii=0; ii<num; ++ii)
  {
    vector<int> bn = updateBoundary(q,candidate);
    candidate.clear();
    // loop over all the boundary elements, add those
    // which are not already part of component and mark
    // them for checking at the next step
    for(size_t jj=0; jj<bn.size(); ++jj)
    {
      int curr = bn[jj];
      int nnbor = 4;
      int nbor[8];
      nbor[0] = q.yp(curr);
      nbor[1] = q.ym(curr);
      nbor[2] = q.xp(curr);
      nbor[3] = q.xm(curr);
      if((ii%2)==0) // do 8 neighbors every other step
      {
        nnbor = 8;
        nbor[4] = q.yp(q.xp(curr));
        nbor[5] = q.yp(q.xm(curr));
        nbor[6] = q.ym(q.xp(curr));
        nbor[7] = q.ym(q.xm(curr));
      }
      for(int kk=0; kk < nnbor; ++kk)
        if(q.get(nbor[kk]) == 0)
        {
          q.put(1,nbor[kk]);
          gl.push_back(nbor[kk]);
          candidate.push_back(nbor[kk]);
        }
    }
  }
  return(gl);
}

vector<int> updateBoundary(const Array2D<int> &q, const vector<int> &idx)
{ // returns a simple list to binaryGrowMN. Requires
  // the expanded list (binary function on whole grid)
  // and a simple list of indices to check. 
  vector<int> bn;
  bn.reserve(idx.size());
  for(size_t ii=0; ii<idx.size(); ++ii)
    if(isbdry(q,idx[ii]))
      bn.push_back(idx[ii]);
  return(bn);
}

int isbdry(const Array2D<int> &q, const int idx)
{ // Checks to see if a given pixel is on the boundary of the
  // set to be expanded, i.e. whether one if it's four 
  // neighbors has value 0
  // Assumes that q is truly a binary function
  return((q.fourNborMin(idx)+1)%2);
}

int myLabelDbl(const Array2D<double> &u, const double thres, Array2D<int> &label, const bool sflag)
{ // Labels connected components of u>thres 
  // if sflag == true, label u>0 with + sign and 0>u>thres with - sign (assumes thres < 0) 
  // if sflag == false, label u>thres with + sign 
  // Returns the number of labels used (number of connected components of u > thres). Labeling put in label.
  const int N = u.getN();

  int numlabels = 0;
  
  label.fillWithValue(0);

  for(int ii=0; ii<N; ++ii)
  {
    if((u.get(ii) >= thres) && (label.get(ii) == 0))
    {
      numlabels++;
      vector<int> list;
      list.reserve(N);
      list.push_back(ii);

      assert(list.size() == 1);

      size_t curr = 0;
      label.put(numlabels,list[curr]);
      while(curr < list.size())
      {
        addNborsDbl(u,thres,label,list,list[curr],numlabels,sflag);
        curr++;
      }
    }
  }
  return(numlabels);
}

void addNborsDbl(const Array2D<double> &u, const double thres, Array2D<int> &label, vector<int> &list, const int idx, const int labelnum, const bool sflag)
{
  int checkind[8]; // indices to check
  checkind[0] = label.xp(idx);
  checkind[1] = label.xm(idx);
  checkind[2] = label.yp(idx);
  checkind[3] = label.ym(idx);
  checkind[4] = label.xp(label.yp(idx));
  checkind[5] = label.xp(label.ym(idx));
  checkind[6] = label.xm(label.yp(idx));
  checkind[7] = label.xm(label.ym(idx));
  for(int ii=0; ii<8; ++ii)
  {
    if((u.get(checkind[ii]) >= thres) && (label.get(checkind[ii]) == 0))
    {
      list.push_back(checkind[ii]);
      if(sflag)
        label.put(labelnum*mysign(u.get(checkind[ii])),checkind[ii]);
      else
        label.put(labelnum,checkind[ii]);
    }
  }
}

/*
int myLabelInt(const Array2D<double> &q, Array2D<int> &label, const bool sflag)
{ // Labels components of q of equal absolute value
  // if sflag == 1, label u>0 with + sign and 0>u>thres with - sign (assumes thres < 0) 
  // if sflag == 0, label u>thres with + sign 
  // Returns the number of labels used (number of connected components of u > thres). Labeling put in label.
  const int N = u.getN();

  int numlabels = 0;
  
  label.fillWithValue(0);

  for(int ii=0; ii<N; ++ii)
  {
    if((u.get(ii) >= thres) && (label.get(ii) == 0))
    {
      numlabels++;
      vector<int> list;
      list.reserve(N);
      list.push_back(ii);

      assert(list.size() == 1);

      int curr = 0;
      label.put(numlabels,list[curr]);
      while(curr < list.size())
      {
        addNborsDbl(u,thres,label,list,list[curr],numlabels,sflag);
        curr++;
      }
    }
  }
  return(numlabels);
}

void addNborsDbl(const Array2D<double> &u, const double thres, Array2D<int> &label, vector<int> &list, const int idx, const int labelnum, const bool sflag)
{
  int checkind[8]; // indices to check
  checkind[0] = label.xp(idx);
  checkind[1] = label.xm(idx);
  checkind[2] = label.yp(idx);
  checkind[3] = label.ym(idx);
  checkind[4] = label.xp(label.yp(idx));
  checkind[5] = label.xp(label.ym(idx));
  checkind[6] = label.xm(label.yp(idx));
  checkind[7] = label.xm(label.ym(idx));
  for(int ii=0; ii<8; ++ii)
  {
    if((u.get(checkind[ii]) >= thres) && (label.get(checkind[ii]) == 0))
    {
      list.push_back(checkind[ii]);
      if(sflag)
        label.put(labelnum*mysign(u.get(checkind[ii])),checkind[ii]);
      else
        label.put(labelnum,checkind[ii]);
    }
  }
}
*/

int* intVectorToArray(const vector<int> &intvec)
{
  int* intarr = new int[intvec.size()];
  for(size_t ii=0; ii<intvec.size(); ++ii)
    intarr[ii] = intvec[ii];
  return(intarr);
}

double* doubleVectorToArray(const vector<double> &doublevec)
{
  double* doublearr = new double[doublevec.size()];
  for(size_t ii=0; ii<doublevec.size(); ++ii)
    doublearr[ii] = doublevec[ii];
  return(doublearr);
}

void addUniqueIndex(vector<int> &vec, const int value) 
{ // add value to the vector vec only if it does not already exist in vec.
  for(vector<int>::iterator it=vec.begin(); it != vec.end(); ++it)
    if(*it == value)
      return;
  vec.push_back(value);
}

bool inVector(const vector<int> &vec, const int value)
{
  for(vector<int>::const_iterator it=vec.begin(); it != vec.end(); ++it)
    if(*it == value)
      return(true);
  return(false);
}

double getRadius(const Array2D<double> u)
{ // assumes that u contains a circle with its center at u[n/2+1,n/2+1], estimates radius of circle
  // by looking along half-line x<(n/2+1), y=n/2+1
  double rad;
  const int n = u.getn();
  int yind = n/2-1;
  int xind = 0;
  while((u.get(yind,xind) < 0.0f) && (xind < n/2))
    xind++;
  // u[yind,xind] > 0, u[yind,xind-1] < 0
  if((xind > 0) && (xind < n/2))
  {
    double corr = -u.get(yind,xind-1) / (u.get(yind,xind)-u.get(yind,xind-1));
    rad = (static_cast<double>(n)/2.0f-1.0f-(static_cast<double>(xind)-1.0f)-corr) / static_cast<double>(n);
  }
  else if(xind == 0)
    rad = 1.0f;
  else
    rad = 0.0f;
  return(rad);
}

Array2D<double> makeCircle(const int n, const double r, const double xc, const double yc)
{
  Array2D<double> u(n,n);
  for(int ii=0; ii<n; ++ii)
  {
    double x = static_cast<double>(ii) / static_cast<double>(n);
    for(int jj=0; jj<n; ++jj)
    {
      double y = static_cast<double>(jj) / static_cast<double>(n);
      u.put(r-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)),jj,ii);
    }
  }
  return(u);
}

Array2D<double> makeCircle(const int n, const double r)
{
  const double ctr = static_cast<double>(0.5 - 1.0f/static_cast<double>(n));
  return(makeCircle(n,r,ctr,ctr));
}

Array2D<double> makeBox(const int n)
{ // make a function which is 0.5 inside a square and -0.5 outside.
  Array2D<double> u(n,n);
  for(int ii=0; ii<n; ++ii)
    for(int jj=0; jj<n; ++jj)
      if((ii>=n/4) && (ii<=3*n/4) && (jj>=n/4) && (jj<=3*n/4))
        u.put(0.5f,jj,ii);
      else
        u.put(-0.5f,jj,ii);
  return(u);
}

double l2err(const Array2D<double> u, const Array2D<double> v)
{
  if((u.getn() != v.getn()) || (u.getm() != v.getn()))
    return(-1.0f);
  const double da = 1.0f / static_cast<double>(u.getN());
  double err = 0.0f;
  for(idx_t ii=0; ii<u.getN(); ++ii)
    err += (v.get(ii)-u.get(ii))*(v.get(ii)-u.get(ii));
  return(sqrt(err*da));
}

double normalize(double (&grad)[2])
{ // assume grad = <grad[0],grad[1]>
  double normsq = grad[0]*grad[0]+grad[1]*grad[1];
  if(normsq < 1e-14)
    if(fabs(grad[0]) > fabs(grad[1]))
    {
      grad[0] = mysign(grad[0]);
      grad[1] = 0.0f;
    }
    else
    {
      grad[0] = 0.0f;
      grad[1] = mysign(grad[1]);
    }
  else
  {
    grad[0] /= sqrt(normsq);
    grad[1] /= sqrt(normsq);
  }
  return(sqrt(normsq));
}

double normalize(struct point &vec)
{ 
  double normsq = vec.x*vec.x+vec.y*vec.y;
  if(normsq < 1e-14)
    if(fabs(vec.x) > fabs(vec.y))
    {
      vec.x = mysign(vec.x);
      vec.y = 0.0f;
    }
    else
    {
      vec.x = 0.0f;
      vec.y = mysign(vec.y);
    }
  else
  {
    vec.x /= sqrt(normsq);
    vec.y /= sqrt(normsq);
  }
  return(sqrt(normsq));
}

void ccomb(const double *x1, const double *x2, double *xr, double theta, const double xlen, const double ylen)
{ // forms convex combination xr = theta*x1 + (1-theta)*x2
  double diff[2];
  if(theta < 0.0f || theta > 1.0f || std::isnan(theta))
  {
    cout << "theta = " << theta << endl;
    //abort();
  }
  if(std::isnan(theta))
    theta = 0.5f;

  diff[0] = pdl(x2[0],x1[0],xlen);
  diff[1] = pdl(x2[1],x1[1],ylen);
  xr[0] = x1[0] + (1.0f-theta)*diff[0];
  xr[1] = x1[1] + (1.0f-theta)*diff[1];

  //printf(" theta = %.3e, xr = (%f,%f).\n",theta,xr[0],xr[1]); fflush(stdout);
}

void setStackSize(const int sizeInBytes)
{
  const rlim_t kStackSize = static_cast<rlim_t>(sizeInBytes); 
  struct rlimit rl;
  
  int result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0)
  {
    if (rl.rlim_cur < kStackSize)
    {
      rl.rlim_cur = kStackSize;
      result = setrlimit(RLIMIT_STACK, &rl);
      if (result != 0)
      {
        fprintf(stderr, "setrlimit returned result = %d\n", result);
      }
    }
  }
  else
    fprintf(stderr, "getrlimit returned result = %d\n", result);
}



