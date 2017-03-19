#ifndef _TOOLBOX_HPP_
#define _TOOLBOX_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <climits>
#include "sys/stat.h"
#include "sys/resource.h"
#include "dirent.h"
#include "defs.h"
#include "array2d.hpp"
#include "assert.h"

using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;

// Function definitions 
bool fileExist(const string &filename);     // Return true if the filename is able to be read, false otherwise.
void readToColon(ifstream &infile);         // Read infile until the next ":" is read and removed from the stream.
void removeTrailingSlash(string &mystring); // Removes internal null characters and trailing slashes.
void myMkdir(const string &dirname);        // Makes the directory dirname
void printn(int n,int newline);
void printtextn(string str, int n, int newline);
void printFileName(string filename);

void fillSign(Array2D<double> &u, const double nullval, const double newval);
void fillSignSubfn(Array2D<double> &u, const double nullval, vector<int> &list, const int curr, const double newval, const double TOL);
void addUniqueIndex(vector<int> &vec, const int value);
bool inVector(const vector<int> &vec, const int value);
vector<vector<int> > createLabelList(const Array2D<int> &q);
Array2D<int> expandList(const vector<int> &ind, const int m, const int n);
vector<int> binaryGrow(const vector<int> &idx, const int m, const int n, const int num);
vector<int> updateBoundary(const Array2D<int> &q, const vector<int> &idx);
int isbdry(const Array2D<int> &q, const int idx);
int myLabelDbl(const Array2D<double> &u, const double thres, Array2D<int> &q, const bool sflag);
void addNborsDbl(const Array2D<double> &u, const double thres, Array2D<int> &q, vector<int> &list, const int idx, const int labelnum, const bool sflag);
int* intVectorToArray(const vector<int> &intvec);
double* doubleVectorToArray(const vector<double> &doublevec);
double getRadius(const Array2D<double> u);
Array2D<double> makeCircle(const int n, const double r, const double xc, const double yc);
Array2D<double> makeCircle(const int n, const double r);
Array2D<double> makeBox(const int n);
double l2err(const Array2D<double> u, const Array2D<double> v);
double normalize(double (&grad)[2]);
double normalize(struct point &vec);
inline double dist(const double *x1, const double *x2);
inline double dist(const double *x1, const double *x2, const double lenx, const double leny);
inline double dist(struct point x1, struct point x2);
inline double dotProd(struct point x1, struct point x2);
void ccomb(const double *x1, const double *x2, double *xr, double theta, const double xlen, const double ylen);

void setStackSize(const int sizeInBytes);

// Inline function definitions
inline double dist(const double *x1, const double *x2)
{ return(sqrt(pd(x1[0],x2[0])*pd(x1[0],x2[0])+pd(x1[1],x2[1])*pd(x1[1],x2[1]))); }

inline double dist(const double *x1, const double *x2, const double lenx, const double leny)
{ return(sqrt(pdl(x1[0],x2[0],lenx)*pdl(x1[0],x2[0],lenx)+pdl(x1[1],x2[1],leny)*pdl(x1[1],x2[1],leny))); }

inline double dist(struct point x1, struct point x2)
{ return(sqrt(pd(x1.x,x2.x)*pd(x1.x,x2.x)+pd(x1.y,x2.y)*pd(x1.y,x2.y))); }

inline double dotProd(struct point x1, struct point x2)
{
  return(x1.x*x2.x+x1.y*x2.y);
}

inline int maxInt()
{ return(INT_MAX); }

#endif
