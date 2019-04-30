{\rtf1\ansi\ansicpg1252\cocoartf1348\cocoasubrtf170
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural

\f0\fs24 \cf0 # include <cstdlib>\
# include <iostream>\
# include <iomanip>\
# include <cmath>\
# include <ctime>\
\
//****************************************************************************80\
//  This file is based on the work of John Burkardt. Some functions defined in the \
//  file Functions, have been added to the file, such that the program computes 3 of \
//  the Genz's integrals: Oscillatory, Product Peak and Gaussian.\
//\
//  Reference:\
//  \
//    This code is a collection of several items, including integrand functions, \
//    an early version of ADAPT, a multidimensional quadrature program, and MULTST, \
//    a routine that tests quadrature programs on the test integrands.  \
//\
//  Modified:\
//\
//    26 May 2007\
//\
//  Author:\
//\
//    Original FORTRAN77 version by Alan Genz\
//    C++ version by John Burkardt\
//\
//  Reference:\
//\
//    Alan Genz,\
//    A Package for Testing Multiple Integration Subroutines,\
//    in Numerical Integration:\
//    Recent Developments, Software and Applications,\
//    edited by Patrick Keast, Graeme Fairweather,\
//    D Reidel, 1987, pages 337-340,\
//    LC: QA299.3.N38.\
//\
//****************************************************************************80\
\
using namespace std;\
\
int main ( void );\
void adapt ( int ndim, double a[], double b[], int *minpts, int maxpts, \
  double functn ( int indx, int ndim, double z[], double alpha[], \
  double beta[] ), double rel_tol, int itest, double alpha[], double beta[],\
  int lenwrk, double wrkstr[], double *relerr, double *finest, int *ifail );\
double genz_function ( int indx, int ndim, double z[], double alpha[], \
  double beta[] );\
double genz_integral ( int indx, int ndim, double a[], double b[], \
  double alpha[], double beta[] );\
char *genz_name ( int indx );\
double genz_phi ( double z );\
double genz_random ( int *seed );\
int i4_max ( int i1, int i2 );\
int i4_min ( int i1, int i2 );\
int i4_power ( int i, int j );\
int i4vec_sum ( int n, int a[] );\
void multst ( int nsamp, int tstlim, int tstfns[], int tstmax, double difclt[], \
  double expnts[], int ndiml, int ndims[], char *sbname, \
  void subrtn ( int ndim, double a[], double b[], int *minpts, int maxpts, \
    double functn ( int indx, int ndim, double z[], double alpha[], \
      double beta[] ), \
    double rel_tol, int itest, double alpha[], double beta[], int lenwrk, \
    double wrkstr[], double *errest, double *finest, int *ifail ), \
  double rel_tol, int maxpts );\
double r8_abs ( double x );\
double r8_epsilon ( );\
double r8_max ( double x, double y );\
double r8_min ( double x, double y );\
double r8vec_dot ( int n, double a1[], double a2[] );\
void r8vec_median ( int n, double r[], double rmed[3] );\
double r8vec_product ( int n, double a[] );\
double r8vec_sum ( int n, double a[] );\
void timestamp ( void );\
void tuple_next ( int m1, int m2, int n, int *rank, int x[] );\
\
//****************************************************************************80\
//  Author:\
//\
//    Original FORTRAN77 version by Alan Genz\
//    C++ version by John Burkardt\
//\
//  Parameters:\
//\
//    Input, int INDX, the index of the test function.\
//\
//    Input, int NDIM, the spatial dimension.\
//\
//    Input, double Z[NDIM], the point at which the integrand \
//    is to be evaluated.\
//\
//    Input, double ALPHA[NDIM], BETA[NDIM], parameters \
//    associated with the integrand function.\
//\
//    Output the value of the test function.\
//****************************************************************************80\
\
int main ( void )\
\
\{\
# define NDIML 5\
# define TSTLIM 6\
# define TSTMAX 6\
\
  double difclt[TSTMAX] = \{ 110.0, 600.0, 600.0, 100.0, 150.0, 100.0 \};\
  double expnts[TSTMAX] = \{ 1.5, 2.0, 2.0, 1.0, 2.0, 2.0 \};\
  int i;\
  int maxpts = 10000;\
  int ndims[NDIML] = \{ 2, 3, 4, 6, 8 \};\
  int nsamp = 20;\
  double rel_tol = 1.0E-06;\
  char *sbname = (char*)"ADAPT";\
  int tstfns[TSTLIM] = \{ 1, 2, 3, 4, 5, 6 \};\
\
  timestamp ( );\
  cout << "\\n";\
  cout << "PS3_Data, Exercice 1.\\n";\
  cout << "  This exercice was performed using the work from John Burkardt\\n";\
  cout << "\\n";\
  cout << "  Reference: \\n";\
  cout << "  This code is a collection of several items, including integrand\\n";\
  cout << "  functions, an early version of ADAPT, a multidimensional quadrature\\n";\
  cout << "  program, and MULTST, a routine that tests quadrature programs on the\\n";\
  cout << "  test integrands.\\n";\
  cout << "  Author: \\n";\
  cout << "     Original FORTRAN77 version by Alan Genz\\n";\
  cout << "     C++ version by John Burkardt.\\n"; \
  cout << "\\n";\
  cout << "  Call MULTST, which can test a routine that\\n";\
  cout << "  is designed to estimate multidimensional\\n";\
  cout << "  integrals, by numerical quadrature.\\n";\
  cout << "\\n";\
  cout << "  The routine to be tested here is called ADAPT.\\n";\
  cout << "\\n";\
  cout << "  The test integrands are Genz's standard set functions.\\n";\
  cout << "\\n";\
  cout << "  MULTST, ADAPT and the test integrands were\\n";\
  cout << "  written in FORTRAN77 by Alan Genz.\\n";\
\
  multst ( nsamp, TSTLIM, tstfns, TSTMAX, difclt, \
    expnts, NDIML, ndims, sbname, adapt, rel_tol, maxpts );\
\
  cout << "\\n";\
  cout << "TEST of Genz's functions\\n";\
  cout << "  Normal end of execution\\n";\
  cout << "\\n";\
  timestamp ( );\
\
  return 0;\
# undef NDIML\
# undef TSTLIM\
# undef TSTMAX \
\}}