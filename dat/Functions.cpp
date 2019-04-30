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
# include "main.hpp"\
# include "adapt.hpp"\
\
using namespace std;\
\
// Purpose of Exercice 1 is to evaluate GENZ function: Oscillatory, Product Peak and Gaussian \
// Using reference from John Burkardt and Alan Genz.\
\
//****************************************************************************80\
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
//   Define the value of the GENZ's function. \
//   Use the index of the function to test, the spatial dimension, the point at which the integrand is to be evaluated and the Alpha, Beta parameters associated with integrand function. \
\
double genz_function ( int indx, int ndim, double z[], double alpha[], \
  double beta[] )\
\
// Define the variables \
\{\
  int j;\
  const double pi = 3.14159265358979323844;\
  bool test;\
  double total;\
  double value;\
\
  value = 0.0;\
\
//  GENZ's Oscillatory function\
//  Set up the function value\
   \
  if ( indx == 1 )\
  \{\
    // Define the function f(x)\
    total = 2.0 * pi * beta[0] + r8vec_sum ( ndim, z );\
    value = cos ( total );\
  \}\
\
//  GENZ's Product Peak function\
//  Set up the function value\
  else if ( indx == 2 )\
  \{\
    total = 1.0;\
    for ( j = 0; j < ndim; j++ )\
    \{\
      total = total * (\
        1.0 / pow ( alpha[j], 2) + pow ( z[j] - beta[j], 2 ) ); // Using pow function to return the base of exp(x)\
    \}\
    value = 1.0 / total;\
  \}\
\
//  GENZ's Gaussian function\
//  C++ math library do not accept things like exp ( -700 )!\
//  Let's define a rule to pass the complain on C++ math   \
  else if ( indx == 4 )\
  \{\
    total = 0.0;\
    for ( j = 0; j < ndim; j++ )\
    \{\
      total = total + pow ( alpha[j] * ( z[j] - beta[j] ), 2 );\
    \}\
    total = r8_min ( total, 100.0 );\
    value = exp ( - total );\
  \}\
  return value;\
\}\
        \
//  Let's define the GENZ's integral, using the values set up in the main.hpp and adapt.hpp\
\
double genz_integral ( int indx, int ndim, double a[], double b[], \
     double alpha[], double beta[] )\
\
// Set up the variables that will be used in the functions\
\{\
  double ab;\
  int *ic;\
  int isum;\
  int j;\
  const double pi = 3.14159265358979323844;\
  int rank;\
  double s;\
  double sgndm;\
  double total;\
  double value;\
\
//  Oscillatory function\
  if ( indx == 1 )\
  \{\
    value = 0.0;\
\
//  Generate all sequences of N dimensions of 0's and 1's.\
    rank = 0;\
    ic = new int[ndim];\
\
    for ( ; ; )\
    \{\
      tuple_next ( 0, 1, ndim, &rank, ic );\
\
      if ( rank == 0 )\
      \{\
        break;\
      \}\
\
      total = 2.0 * pi * beta[0];\
      for ( j = 0; j < ndim; j++ )\
      \{\
        if ( ic[j] != 1 )\
        \{\
          total = total + alpha[j];\
        \}\
      \}\
\
      isum = i4vec_sum ( ndim, ic );\
\
      s = 1 + 2 * ( ( isum / 2 ) * 2 - isum );\
\
      if ( ( ndim % 2 ) == 0 )\
      \{\
        value = value + s * cos ( total );\
      \}\
      else\
      \{\
        value = value + s * sin ( total );\
      \}\
    \}\
    delete [] ic;\
\
    if ( 1 < ( ndim % 4 ) )\
    \{\
      value = - value;\
    \}\
  \}\
\
  //  Product Peak.\
  else if ( indx == 2 )\
  \{\
    value = 1.0;\
\
    for ( j = 0; j < ndim; j++ )\
    \{\
      value = value * alpha[j] * ( \
          atan ( ( 1.0 - beta[j] ) * alpha[j] ) \
        + atan (       + beta[j]   * alpha[j] ) );\
    \}\
  \}\
\
//  Gaussian\
  else if ( indx == 4 )\
  \{\
    value = 1.0;\
\
    ab = sqrt ( 2.0 );\
    for ( j = 0; j < ndim; j++ )\
    \{\
      value = value * ( sqrt ( pi ) / alpha[j] ) * \
        (   genz_phi ( ( 1.0 - beta[j] ) * ab * alpha[j] ) \
          - genz_phi (       - beta[j]   * ab * alpha[j] ) );\
    \}\
  \}\
    \
  return value;\
\}\
\
\
// Use the char* to store letters. Define the letters to store as the different names of the GENZ's function\
\
char *genz_name ( int indx ) // Name of the test integrand\
\
\{\
  char *name;\
\
  name = new char[14];\
\
  if ( indx == 1 )\
  \{\
    strcpy ( name, "Oscillatory  " );\
  \}\
  else if ( indx == 2 )\
  \{\
    strcpy ( name, "Product Peak " );\
  \}\
  else if ( indx == 4 )\
  \{\
    strcpy ( name, "Gaussian     " );\
  \}\
  else\
  \{\
    cout << "\\n";\
    cout << "  GENZ_NAME - Fatal error!\\n";\
    cout << "  1 <= INDX <= 6 is required.\\n";\
    exit ( 1 );\
  \}\
  return name;\
\}\
\
// Define the double integral for the GENZ's functions. This was performed using the function GENZ_PHI from John Burkart\
    \
double genz_phi ( double z )\
\
//  Purpose:\
//\
//    GENZ_PHI estimates the normal cumulative density function.\
//\
//  Discussion:\
//\
//    The approximation is accurate to 1.0E-07.\
//\
//    This routine is based upon algorithm 5666 for the error function,\
//    from Hart et al.\
//\
//  Modified:\
//\
//    20 March 2007\
//\
//  Author:\
//\
//    Original FORTRAN77 version by Alan Miller\
//    C++ version by John Burkardt\
//\
//  Reference:\
//\
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,\
//    Charles Mesztenyi, John Rice, Henry Thatcher,\
//    Christoph Witzgall,\
//    Computer Approximations,\
//    Wiley, 1968,\
//    LC: QA297.C64.\
//\
//  Parameters:\
//\
//    Input, double Z, a value which can be regarded as the distance,\
//    in standard deviations, from the mean.\
//\
//    Output, double GENZ_PHI, the integral of the normal PDF from negative\
//    infinity to Z.\
//\
//  Local parameters:\
//\
//    Local, double ROOTPI, despite the name, is actually the \
//    square root of TWO * pi.\
   \
\{\
  double expntl;\
  double p;\
  const double p0 = 220.2068679123761;\
  const double p1 = 221.2135961699311;\
  const double p2 = 112.0792914978709;\
  const double p3 = 33.91286607838300;\
  const double p4 = 6.373962203531650;\
  const double p5 = 0.7003830644436881;\
  const double p6 = 0.03526249659989109;\
  const double q0 = 440.4137358247522;\
  const double q1 = 793.8265125199484;\
  const double q2 = 637.3336333788311;\
  const double q3 = 296.5642487796737;\
  const double q4 = 86.78073220294608;\
  const double q5 = 16.06417757920695;\
  const double q6 = 1.755667163182642;\
  const double q7 = 0.08838834764831844;\
  const double rootpi = 2.506628274631001;\
  double zabs;\
\
  zabs = r8_abs ( z );\
//\
//  12 < |Z|.\
//\
  if ( 12.0 < zabs )\
  \{\
    p = 0.0;\
  \}\
  else\
  \{\
//\
//  |Z| <= 12\
//\
    expntl = exp ( - zabs * zabs / 2.0 );\
//\
//  |Z| < 7\
//\
    if ( zabs < 7.0 )\
    \{\
      p = expntl * (((((( \
                  p6 \
         * zabs + p5 ) \
         * zabs + p4 ) \
         * zabs + p3 ) \
         * zabs + p2 ) \
         * zabs + p1 ) \
         * zabs + p0 ) / ((((((( \
                  q7 \
         * zabs + q6 ) \
         * zabs + q5 ) \
         * zabs + q4 ) \
         * zabs + q3 ) \
         * zabs + q2 ) \
         * zabs + q1 ) \
         * zabs + q0 );\
    \}\
//\
//  CUTOFF <= |Z|\
//\
    else\
    \{\
      p = expntl / ( \
        zabs + 1.0 / (\
        zabs + 2.0 / ( \
        zabs + 3.0 / ( \
        zabs + 4.0 / ( \
        zabs + 0.65 ))))) / rootpi;\
    \}\
  \}\
\
  if ( 0.0 < z )\
  \{\
    p = 1.0 - p;\
  \}\
\
  return p;\
\}\
//****************************************************************************80\
// Using the function genz_random based on Linus Schrage work.\
\
double genz_random ( int *seed )\
\
//****************************************************************************80\
//\
//  Purpose:\
//\
//    GENZ_RANDOM is a portable random number generator\
//\
//  Modified:\
//\
//    27 May 2007\
//\
//  Author:\
//\
//    Original FORTRAN77 version by Linus Schrage\
//\
//  Reference:\
//\
//    Linus Schrage,\
//    A More Portable Fortran Random Number Generator,\
//    ACM Transactions on Mathematical Software,\
//    Volume 5, Number 2, June 1979, pages 132-138.\
//\
//  Parameters:\
//\
//    Input, integer/output, int *SEED, a seed for the random\
//    number generator.\
//\
//    Output, double GENZ_RANDOM, a pseudorandom value.\
//\
\{\
  const int a = 16807;\
  const int b15 = 32768;\
  const int b16 = 65536;\
  int fhi;\
  int k;\
  int leftlo;\
  const int p = 2147483647;\
  double value;\
  int xalo;\
  int xhi;\
\
  xhi = *seed / b16;\
  xalo = ( *seed - xhi * b16 ) * a;\
  leftlo = xalo / b16;\
  fhi = xhi * a + leftlo;\
  k = fhi / b15;\
\
  *seed = ( \
            ( \
              ( xalo - leftlo * b16 ) - p \
            ) \
          + ( fhi - k * b15 ) * b16 \
          ) + k;\
\
  if ( *seed < 0 )\
  \{\
    *seed = *seed + p;\
  \}\
\
  value = ( double ) ( *seed ) / ( double ) ( p );\
\
  return value;\
\}\
//****************************************************************************80\
\
\
// Define the function based on the max between two integers to be compared. The function returns the larger of the two. \
\
int i4_max ( int i1, int i2 )// Integer 1 and Integer 2\
\
\{\
  int value;\
\
  if ( i2 < i1 )\
  \{\
    value = i1;\
  \}\
  else\
  \{\
    value = i2;\
  \}\
  return value;\
\}\
\
// Define the function based on the minimum between two integers to be compared. The function returns the smaller of the two.\
int i4_min ( int i1, int i2 )//Integer 1, Integer 2\
\
\{\
  int value;\
\
  if ( i1 < i2 )\
  \{\
    value = i1;\
  \}\
  else\
  \{\
    value = i2;\
  \}\
  return value;\
\}\
\
// Returns the value of the power of I^J. Careful: the power cannot be negative. \
int i4_power ( int i, int j )// Int i = the base and int j = the power\
\
\{\
  int k;\
  int value;\
\
// Set up the rule for nonnegativity of the power J\
  if ( j < 0 )\
  \{\
    if ( i == 1 )\
    \{\
      value = 1;\
    \}\
    else if ( i == 0 )\
    \{\
      cout << "\\n";\
      cout << "I4_POWER - Fatal error!\\n";\
      cout << "  I^J requested, with I = 0 and J negative.\\n";\
      exit ( 1 );\
    \}\
    else\
    \{\
      value = 0;\
    \}\
  \}\
  else if ( j == 0 )\
  \{\
    if ( i == 0 )\
    \{\
      cout << "\\n";\
      cout << "I4_POWER - Fatal error!\\n";\
      cout << "  I^J requested, with I = 0 and J = 0.\\n";\
      exit ( 1 );\
    \}\
    else\
    \{\
      value = 1;\
    \}\
  \}\
  else if ( j == 1 )\
  \{\
    value = i;\
  \}\
  else\
  \{\
    value = 1;\
    for ( k = 1; k <= j; k++ )\
    \{\
      value = value * i; //I^J\
    \}\
  \}\
  return value;\
\}\
\
// Sum of the entries of a vector.\
int i4vec_sum ( int n, int a[] ) // int n = number of entries in the vector and int[a] is the vector to be summed.\
\
\{\
  int i;\
  int sum;\
\
  sum = 0;\
  for ( i = 0; i < n; i++ )\
  \{\
    sum = sum + a[i];\
  \}\
\
  return sum;\
\}\
\
\
//****************************************************************************80\
// Using MLUST for John Burkardt\
\
void multst ( int nsamp, int tstlim, int tstfns[], int tstmax, double difclt[], \
  double expnts[], int ndiml, int ndims[], char *sbname, \
  void subrtn ( int ndim, double a[], double b[], int *minpts, int maxpts, \
    double functn ( int indx, int ndim, double z[], double alpha[], \
      double beta[] ), \
    double rel_tol, int itest, double alpha[], double beta[], int lenwrk, \
    double wrkstr[], double *errest, double *finest, int *ifail ), \
  double rel_tol, int maxpts )\
\
//****************************************************************************80\
//\
//  Purpose:\
//\
//    MULTST tests a multidimensional integration routine.\
//\
//  Discussion:\
//\
//    The routine uses the test integrand functions, with\
//    the user selecting the particular subset of test integrands,\
//    the set of difficulty factors, and the spatial dimensions.\
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
//  Parameters:\
//\
//    Input, int NSAMP, the number of samples.\
//    1 <= NSAMP.\
//\
//    Input, int TSTLIM, the number of test integrands.\
//\
//    Input, int TSTFNS[TSTLIM], the indices of the test integrands.\
//    Each index is between 1 and 6.\
//\
//    Input, int TSTMAX, the number of difficulty levels to be tried.\
//\
//    Input, double DIFCLT[TSTMAX], difficulty levels.\
//\
//    Input, double EXPNTS[TSTMAX], the difficulty exponents.\
//\
//    Input, int NDIML, the number of sets of variable sizes.\
//\
//    Input, int NDIMS[NDIML], the number of variables for the integrals\
//    in each test.\
//\
//    Input, char *SBNAME, the name of the integration\
//    subroutine to be tested.\
//\
//    Input, external SUBRTN, the integration subroutine to be tested.\
//\
//    Input, double REL_TOL, the relative error tolerance.\
//\
//    Input, int MAXPTS, the maximum number of integrand calls\
//    for all tests.\
//\
\{\
# define MXTSFN 6\
\
  double *a;\
  double *alpha;\
  double *b;\
  double *beta;\
  double callsa[MXTSFN*MXTSFN];\
  double callsb[MXTSFN*MXTSFN];\
  double concof;\
  double dfact;\
  double dfclt;\
  int digits;\
  double errest;\
  double errlog;\
  double ersacb[MXTSFN*MXTSFN];\
  double ersact[MXTSFN*MXTSFN];\
  double ersdsb[MXTSFN*MXTSFN];\
  double ersdsc[MXTSFN*MXTSFN];\
  double ersesb[MXTSFN*MXTSFN];\
  double ersest[MXTSFN*MXTSFN];\
  double ersrel[MXTSFN*MXTSFN];\
  double estlog;\
  double exn;\
  double expons[MXTSFN];\
  double finest;\
  int i;\
  int idfclt[MXTSFN];\
  int ifail;\
  int ifails;\
  int it;\
  int itest;\
  int j;\
  int k;\
  int lenwrk;\
  double medacb[MXTSFN];\
  double medacb_med[3];\
  double *medact;\
  double medact_med[3];\
  double medcla[MXTSFN];\
  double medcla_med[3];\
  double medclb[MXTSFN];\
  double medclb_med[3];\
  double *medcls;\
  double medcls_med[3];\
  double meddsb[MXTSFN];\
  double meddsb_med[3];\
  double *meddsc;\
  double meddsc_med[3];\
  double medesb[MXTSFN];\
  double medesb_med[3];\
  double *medest;\
  double medest_med[3];\
  double medrel;\
  double *medrll;\
  double medrll_med[3];\
  int minpts;\
  int n;\
  char *name;\
  int nconf;\
  int ndim;\
  int ndimv;\
  double qality;\
  double *qallty;\
  double qallty_med[3];\
  double qualty[MXTSFN*MXTSFN];\
  int rcalsa;\
  int rcalsb;\
  double relerr;\
  int rulcls;\
  int seed;\
  double small;\
  double tactrb[MXTSFN];\
  double tactrb_med[3];\
  double tactrs[MXTSFN];\
  double tactrs_med[3];\
  double tcalsa[MXTSFN];\
  double tcalsa_med[3];\
  double tcalsb[MXTSFN];\
  double tcalsb_med[3];\
  double terdsb[MXTSFN];\
  double terdsb_med[3];\
  double terdsc[MXTSFN];\
  double terdsc_med[3];\
  double testrb[MXTSFN];\
  double testrb_med[3];\
  double testrs[MXTSFN];\
  double testrs_med[3];\
  double tqualt[MXTSFN];\
  double tqualt_med[3];\
  double total;\
  double trelib[MXTSFN];\
  double trelib_med[3];\
  double value;\
  double *wrkstr;\
\
  medact = new double[nsamp];\
  medcls = new double[nsamp];\
  meddsc = new double[nsamp];\
  medest = new double[nsamp];\
  medrll = new double[nsamp];\
  qallty = new double[nsamp];\
//\
//  Initialize and compute confidence coefficient.\
//\
  concof = 0.0;\
  nconf = i4_max ( 1, ( 2 * nsamp ) / 5 - 2 );\
\
  for ( i = 1; i <= nconf; i++ )\
  \{\
    concof = 1.0 + ( double ) ( nsamp - nconf + i ) * concof \
      / ( double ) ( nconf - i + 1 );\
  \}\
\
  concof = 1.0 - concof / ( double ) ( i4_power ( 2, nsamp - 1 ) );\
\
  seed = 123456;\
\
  small = r8_epsilon ( );\
\
  for ( i = 0; i < tstlim; i++ )\
  \{\
    idfclt[i] = ( int ) difclt[tstfns[i]-1];\
  \}\
  for ( i = 0; i < tstlim; i++ )\
  \{\
    expons[i] = expnts[tstfns[i]-1];\
  \}\
//\
//  Begin main loop for different numbers of variables.\
//\
  for ( ndimv = 0; ndimv < ndiml; ndimv++ )\
  \{\
    ndim = ndims[ndimv];\
\
    a = new double[ndim];\
    alpha = new double[ndim];\
    b = new double[ndim];\
    beta = new double[ndim];\
\
    if ( ndim <= 15 )\
    \{\
      rulcls = i4_power ( 2, ndim ) + 2 * i4_power ( ndim, 2 ) + 2 * ndim + 1;\
    \}\
    else\
    \{\
      rulcls = ( ndim * ( 14 - ndim * ( 6 - 4 * ndim ) ) ) / 3 + 1;\
    \}\
\
    lenwrk = ( 2 * ndim + 3 ) * ( 1 + maxpts / rulcls ) / 2;\
    wrkstr = new double[lenwrk];\
\
    if ( ( ndimv % 6 ) == 0 )\
    \{\
      cout << "\\n";\
      cout << "  Test results with " << nsamp << " samples per test.\\n";\
      cout << "\\n";\
      cout << "  Difficulty levels";\
      for ( j = 0; j < tstlim; j++ )\
      \{\
        cout << setw(6) << idfclt[j];\
      \}\
      cout << "\\n";\
      cout << "          Exponents";\
      for ( j = 0; j < tstlim; j++ )\
      \{\
        cout << setw(6) << expons[j];\
      \}\
      cout << "\\n";\
\
      digits = ( int ) ( -log10 ( rel_tol ) );\
\
      cout << "\\n";\
      cout << "  Requested digits = " << digits\
           << " Maximum values = " << maxpts << "\\n";\
      cout << "\\n";\
      cout << sbname << " tests, variable results with confidence "\
           << concof << "\\n";\
      cout << "\\n";\
      cout << " Vari-  Integrand     Correct digits   Relia-  Wrong"\
           << "   Integrand   Quality Total\\n";\
      cout << " ables              Estimated   Actual bility Digits"\
           << "    Values             Fails\\n";\
      cout << "\\n";\
    \}\
//\
//  Begin loop for different test integrands.\
//\
    for ( it = 0; it < tstlim; it++ )\
    \{\
      itest = tstfns[it];\
      exn = expnts[itest-1];\
      dfclt = difclt[itest-1];\
\
      for ( j = 0; j < ndim; j++ )\
      \{\
        a[j] = 0.0;\
      \}\
      for ( j = 0; j < ndim; j++ )\
      \{\
        b[j] = 1.0;\
      \}\
      ifails = 0;\
      medrel = 0;\
//\
//  Begin loop for different samples.\
//\
      for ( k = 0; k < nsamp; k++ )\
      \{\
        ifail = 1;\
//\
//  Choose the integrand function parameters at random.\
//\
        for ( n = 0; n < ndim; n++ )\
        \{\
          alpha[n] = genz_random ( &seed );\
          beta[n] = genz_random ( &seed );\
        \}\
//\
//  Modify ALPHA to account for difficulty parameter.\
//\
        total = r8vec_sum ( ndim, alpha );\
        dfact = total * pow ( ndim, exn ) / dfclt;\
        for ( j = 0; j < ndim; j++ )\
        \{\
          alpha[j] = alpha[j] / dfact;\
        \}\
//\
//  For tests 1 and 3, we modify the value of B.\
//\
        if ( itest == 1 || itest == 3 )\
        \{\
          for ( j = 0; j < ndim; j++ )\
          \{\
            b[j] = alpha[j];\
          \}\
        \}\
//\
//  For test 6, we modify the value of BETA.\
//\
        if ( itest == 6 )\
        \{\
          for ( n = 2; n < ndim; n++ )\
          \{\
            beta[n] = 1.0;\
          \}\
        \}\
//\
//  Get the exact value of the integral.\
//\
        value = genz_integral ( itest, ndim, a, b, alpha, beta );//Use the function defined as the Genz integral\
//\
//  Call the integration subroutine.\
//\
        minpts = 4 * i4_power ( 2, ndim );\
\
        subrtn ( ndim, a, b, &minpts, maxpts, genz_function, rel_tol, \
          itest, alpha, beta, lenwrk, wrkstr, &errest, &finest, &ifail );\
\
        relerr = r8_abs ( ( finest - value ) / value );\
        ifails = ifails + i4_min ( ifail, 1 );\
        relerr = r8_max ( r8_min ( 1.0, relerr ), small );\
        errlog = r8_max ( 0.0, -log10 ( relerr ) );\
        errest = r8_max ( r8_min ( 1.0, errest ), small );\
        estlog = r8_max ( 0.0, -log10 ( errest ) );\
        meddsc[k] = r8_max ( 0.0, estlog - errlog );\
        medest[k] = estlog;\
        medact[k] = errlog;\
        medcls[k] = minpts;\
\
        if ( relerr <= errest )\
        \{\
          medrel = medrel + 1;\
        \}\
      \}\
//\
//  End loop for different samples and compute medians.\
//\
      r8vec_median ( nsamp, medest, medest_med );\
      r8vec_median ( nsamp, medact, medact_med );\
      r8vec_median ( nsamp, medcls, medcls_med );\
      r8vec_median ( nsamp, meddsc, meddsc_med );\
\
      medrel = medrel / ( double ) ( nsamp );\
\
      trelib[it] = medrel;\
\
      tactrs[it] = medact_med[1];\
      testrs[it] = medest_med[1];\
      terdsc[it] = meddsc_med[1];\
      tcalsa[it] = medcls_med[1];\
\
      tcalsb[it] = medcls_med[2];\
      tactrb[it] = medact_med[2];\
      testrb[it] = medest_med[2];\
      terdsb[it] = meddsc_med[2];\
\
      ersrel[itest-1+ndimv*MXTSFN] = medrel;\
\
      ersest[itest-1+ndimv*MXTSFN] = medest_med[1];\
      ersact[itest-1+ndimv*MXTSFN] = medact_med[1];\
      ersdsc[itest-1+ndimv*MXTSFN] = meddsc_med[1];\
\
      ersesb[itest-1+ndimv*MXTSFN] = medest_med[2];\
      ersacb[itest-1+ndimv*MXTSFN] = medact_med[2];\
      ersdsb[itest-1+ndimv*MXTSFN] = meddsc_med[2];\
\
      callsa[itest-1+ndimv*MXTSFN] = medcls_med[1];\
\
      callsb[itest-1+ndimv*MXTSFN] = medcls_med[2];\
\
      qality = 0.0;\
\
      if ( medcls_med[0] != 0.0 )\
      \{\
        qality = ( medact_med[0] + 1.0 ) * \
          ( medest_med[0] + 1.0 - meddsc_med[0] ) / log ( medcls_med[0] );\
      \}\
\
      tqualt[it] = qality;\
      qualty[itest-1+ndimv*MXTSFN] = qality;\
      rcalsa = ( int ) medcls_med[1];\
      rcalsb = ( int ) medcls_med[2];\
      name = genz_name ( itest );\
\
      cout << setw(4) << ndim\
           << "  " << setw(14) << name\
           << setprecision(2) << setw(4) << medest_med[1]\
           << setprecision(2) << setw(5) << medest_med[2]\
           << setprecision(2) << setw(5) << medact_med[1]\
           << setprecision(2) << setw(5) << medact_med[2]\
           << setprecision(3) << setw(5) << medrel\
           << setprecision(2) << setw(4) << meddsc_med[1]\
           << setprecision(2) << setw(4) << meddsc_med[2]\
           << setw(7) << rcalsa\
           << setw(8) << rcalsb\
           << setprecision(3) << setw(6) << qality\
           << setw(5) << ifails << "\\n";\
      delete [] name;\
    \}\
//\
//  End loop for different test integrands.\
//\
    r8vec_median ( tstlim, tactrs, tactrs_med );\
    r8vec_median ( tstlim, trelib, trelib_med );\
    r8vec_median ( tstlim, testrs, testrs_med );\
    r8vec_median ( tstlim, terdsc, terdsc_med );\
    r8vec_median ( tstlim, tactrb, tactrb_med );\
    r8vec_median ( tstlim, testrb, testrb_med );\
    r8vec_median ( tstlim, terdsb, terdsb_med );\
    r8vec_median ( tstlim, tqualt, tqualt_med );\
    r8vec_median ( tstlim, tcalsa, tcalsa_med );\
    r8vec_median ( tstlim, tcalsb, tcalsb_med );\
\
    rcalsa = ( int ) tcalsa_med[0];\
    rcalsb = ( int ) tcalsb_med[0];\
\
    cout << setw(4) << ndim\
         << "  Medians       "\
         << setprecision(2) << setw(4) << testrs_med[0]\
         << setprecision(2) << setw(5) << testrb_med[0]\
         << setprecision(2) << setw(5) << testrs_med[0]\
         << setprecision(2) << setw(5) << tactrb_med[0]\
         << setprecision(3) << setw(5) << trelib_med[0]\
         << setprecision(2) << setw(4) << terdsc_med[0]\
         << setprecision(2) << setw(4) << terdsb_med[0]\
         << setw(7) << rcalsa\
         << setw(8) << rcalsb\
         << setprecision(3) << setw(6) << tqualt_med[0] << "\\n";\
\
    cout << "\\n";\
\
    delete [] a;\
    delete [] alpha;\
    delete [] b;\
    delete [] beta;\
    delete [] wrkstr;\
  \}\
//\
//  End loop for different numbers of variables.\
//\
  if ( 1 < ndiml )\
  \{\
    cout << "\\n";\
    cout << "      " << sbname << " Test integrand medians for variables";\
    for ( j = 0; j < ndiml; j++ )\
    \{\
      cout << setw(3) << ndims[j];\
    \}\
    cout << "\\n";\
\
    cout << "\\n";\
    cout << "        Integrand     Correct digits   Relia-  Wrong"\
         << "   Integrand   Quality\\n";\
    cout << "                    Estimated   Actual bility digits"\
         << "     Values\\n";\
    cout << "\\n";\
\
    for ( it = 0; it < tstlim; it++ )\
    \{\
      itest = tstfns[it];\
\
      for ( j = 0; j < ndiml; j++ )\
      \{\
        medact[j] = ersact[itest-1+j*MXTSFN];\
        medest[j] = ersest[itest-1+j*MXTSFN];\
        meddsc[j] = ersdsc[itest-1+j*MXTSFN];\
        medacb[j] = ersacb[itest-1+j*MXTSFN];\
        medesb[j] = ersesb[itest-1+j*MXTSFN];\
        meddsb[j] = ersdsb[itest-1+j*MXTSFN];\
        medrll[j] = ersrel[itest-1+j*MXTSFN];\
        qallty[j] = qualty[itest-1+j*MXTSFN];\
        medcla[j] = callsa[itest-1+j*MXTSFN];\
        medclb[j] = callsb[itest-1+j*MXTSFN];\
      \}\
\
      r8vec_median ( ndiml, medrll, medrll_med );\
      r8vec_median ( ndiml, medact, medact_med );\
      r8vec_median ( ndiml, medest, medest_med );\
      r8vec_median ( ndiml, meddsc, meddsc_med );\
      r8vec_median ( ndiml, medacb, medacb_med );\
      r8vec_median ( ndiml, medesb, medesb_med );\
      r8vec_median ( ndiml, meddsb, meddsb_med );\
      r8vec_median ( ndiml, qallty, qallty_med );\
      r8vec_median ( ndiml, medcla, medcla_med );\
      r8vec_median ( ndiml, medclb, medclb_med );\
\
      rcalsa = ( int ) medcla_med[0];\
      rcalsb = ( int ) medclb_med[0];\
      name = genz_name ( itest );\
\
      cout << "      " << setw(14) << name\
           << setprecision(2) << setw(4) << medest_med[0]\
           << setprecision(2) << setw(5) << medesb_med[0]\
           << setprecision(2) << setw(5) << medact_med[0]\
           << setprecision(2) << setw(5) << medacb_med[0]\
           << setprecision(3) << setw(5) << medrll_med[0]\
           << setprecision(2) << setw(4) << meddsc_med[0]\
           << setprecision(2) << setw(4) << meddsb_med[0]\
           << setw(7) << rcalsa\
           << setw(8) << rcalsb\
           << setprecision(3) << setw(6) << qallty_med[0]\
           << setw(5) << ifails << "\\n";\
\
      delete [] name;\
\
      tactrs[it] = medact_med[0];\
      testrs[it] = medest_med[0];\
      terdsc[it] = meddsc_med[0];\
      tactrb[it] = medacb_med[0];\
      testrb[it] = medesb_med[0];\
      terdsb[it] = meddsb_med[0];\
      tcalsa[it] = medcla_med[0];\
      tcalsb[it] = medclb_med[0];\
      trelib[it] = medrll_med[0];\
      tqualt[it] = qallty_med[0];\
    \}\
\
    r8vec_median ( tstlim, tactrs, tactrs_med );\
    r8vec_median ( tstlim, testrs, testrs_med );\
    r8vec_median ( tstlim, terdsc, terdsc_med );\
    r8vec_median ( tstlim, tactrb, tactrb_med );\
    r8vec_median ( tstlim, testrb, testrb_med );\
    r8vec_median ( tstlim, terdsb, terdsb_med );\
    r8vec_median ( tstlim, trelib, trelib_med );\
    r8vec_median ( tstlim, tqualt, tqualt_med );\
    r8vec_median ( tstlim, tcalsa, tcalsa_med );\
    r8vec_median ( tstlim, tcalsb, tcalsb_med );\
\
    rcalsa = ( int ) tcalsa_med[0];\
    rcalsb = ( int ) tcalsb_med[0];\
\
    cout << "      Global medians"\
         << setprecision(2) << setw(4) << testrs_med[0]\
         << setprecision(2) << setw(5) << testrb_med[0]\
         << setprecision(2) << setw(5) << tactrs_med[0]\
         << setprecision(2) << setw(5) << tactrb_med[0]\
         << setprecision(3) << setw(5) << trelib_med[0]\
         << setprecision(2) << setw(4) << terdsc_med[0]\
         << setprecision(2) << setw(4) << terdsb_med[0]\
         << setw(7) << rcalsa\
         << setw(8) << rcalsb\
         << setprecision(3) << setw(6) << tqualt_med[0]\
         << setw(5) << ifails << "\\n";\
\
    cout << "\\n";\
  \}\
\
  delete [] medact;\
  delete [] medcls;\
  delete [] meddsc;\
  delete [] medest;\
  delete [] medrll;\
  delete [] qallty;\
\
  return;\
# undef MXTSFN\
\}\
//****************************************************************************80\
\
// Define the function, such that it returns the absolute value of X. \
double r8_abs ( double x )\
\
\{\
  double value;\
\
  if ( 0.0 <= x )\
  \{\
    value = x;\
  \} \
  else\
  \{\
    value = -x;\
  \}\
  return value;\
\}\
\
// Define the epsilon function: return the roundoff unit. \
// The function was set up using John Burkardt file.\
\
double r8_epsilon ( )\
\
//****************************************************************************80\
//\
//  Purpose:\
//\
//    R8_EPSILON returns the R8 roundoff unit.\
//\
//  Discussion:\
//\
//    The roundoff unit is a number R which is a power of 2 with the\
//    property that, to the precision of the computer's arithmetic,\
//      1 < 1 + R\
//    but\
//      1 = ( 1 + R / 2 )\
//\
//  Licensing:\
//\
//    This code is distributed under the GNU LGPL license.\
//\
//  Modified:\
//\
//    01 September 2012\
//\
//  Author:\
//\
//    John Burkardt\
//\
//  Parameters:\
//\
//    Output, double R8_EPSILON, the R8 round-off unit.\
//\
\{\
  const double value = 2.220446049250313E-016;\
\
  return value;\
\}\
//****************************************************************************80\
\
// Define the function so that it returns the maximum value between X and Y. \
\
double r8_max ( double x, double y )// The quantities to compare\
\
\{\
  double value;\
\
  if ( y < x )\
  \{\
    value = x;\
  \} \
  else\
  \{\
    value = y;\
  \}\
  return value;\
\}\
\
\
// The function is defined in such as way that it returns the minimum value between X and Y.\
\
double r8_min ( double x, double y )// Values to compare\
\
\{\
  double value;\
\
  if ( y < x )\
  \{\
    value = y;\
  \} \
  else\
  \{\
    value = x;\
  \}\
  return value;\
\}\
\
// This function computes the dot product of a pair of vectors. Where A1[] and A2[] are the two vectors used as value. \
\
double r8vec_dot ( int n, double a1[], double a2[] )// Int n = the number of entries in the vector. \
\
\{\
  int i;\
  double value;\
\
  value = 0.0;\
  for ( i = 0; i < n; i++ )\
  \{\
    value = value + a1[i] * a2[i];\
  \}\
\
  return value;\
\}\
\
//****************************************************************************80\
// Median function: taken form John Burkardt \
\
void r8vec_median ( int n, double r[], double rmed[3] )\
\
//****************************************************************************80\
//\
//  Purpose:\
//\
//    R8VEC_MEDIAN computes the median of an R8VEC.\
//\
//  Modified:\
//\
//    20 March 2007\
//\
//  Author:\
//\
//    Original FORTRAN77 version by Alan Genz\
//    C++ version by John Burkardt\
//\
//  Parameters:\
//\
//    Input, int N, the dimension of the array.\
//\
//    Input, double R[N], the array to be examined.\
//\
//    Output, double RMED[3].  RMED[0] contains the median,\
//    RMED[1] and RMED[2] specify the confidence interval.\
//\
\{\
  int j;\
  int k;\
  int kmax;\
  int nconf;\
  int nd;\
  double rmax;\
\
  for ( j = 0; j < n; j++ )\
  \{\
    kmax = j;\
\
    for ( k = j + 1; k < n; k++ )\
    \{\
      if ( r[kmax] < r[k] )\
      \{\
        kmax = k;\
      \}\
    \}\
    rmax = r[kmax];\
    r[kmax] = r[j];\
    r[j] = rmax;\
  \}\
\
  nd = n / 2;\
\
  if ( ( n % 2 ) == 0 )\
  \{\
    rmed[0] = ( r[nd-1] + r[nd] ) / 2.0;\
  \}\
  else\
  \{\
    rmed[0] = r[nd];\
  \}\
\
  nconf = i4_max ( 1, ( 2 * n ) / 5 - 2 );\
\
  rmed[1] = r[n-nconf];\
  rmed[2] = r[nconf-1];\
\
  return;\
\}\
//****************************************************************************80\
\
double r8vec_product ( int n, double a[] )\
\
//****************************************************************************80\
//\
//  Purpose:\
//\
//    R8VEC_PRODUCT returns the product of the entries of an R8VEC.\
//\
//  Modified:\
//\
//    17 September 2003\
//\
//  Author:\
//\
//    John Burkardt\
//\
//  Parameters:\
//\
//    Input, int N, the number of entries in the vector.\
//\
//    Input, double A[N], the vector.\
//\
//    Output, double R8VEC_PRODUCT, the product of the vector.\
//\
\{\
  int i;\
  double product;\
\
  product = 1.0;\
  for ( i = 0; i < n; i++ )\
  \{\
    product = product * a[i];\
  \}\
\
  return product;\
\}\
//****************************************************************************80\
\
// The function returns the sum of vectors \
double r8vec_sum ( int n, double a[] )//int n = number of entries in the vector and double A[] = the vector\
\
\{\
  int i;\
  double sum;\
\
  sum = 0.0;\
  for ( i = 0; i < n; i++ )\
  \{\
    sum = sum + a[i];\
  \}\
\
  return sum;\
\}\
\
//****************************************************************************80\
//  Using the timestamp function from the work of John Burkardt\
\
void timestamp ( void )\
\
//****************************************************************************80\
//\
//  Purpose:\
//\
//    TIMESTAMP prints the current YMDHMS date as a time stamp.\
//\
//  Example:\
//\
//    May 31 2001 09:45:54 AM\
//\
//  Modified:\
//\
//    24 September 2003\
//\
//  Author:\
//\
//    John Burkardt\
//\
//  Parameters:\
//\
//    None\
//\
\{\
# define TIME_SIZE 40\
\
  static char time_buffer[TIME_SIZE];\
  const struct tm *tm;\
  size_t len;\
  time_t now;\
\
  now = time ( NULL );\
  tm = localtime ( &now );\
\
  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );\
\
  cout << time_buffer << "\\n";\
\
  return;\
# undef TIME_SIZE\
\}\
//****************************************************************************80\
\
void tuple_next ( int m1, int m2, int n, int *rank, int x[] )\
\
//****************************************************************************80\
//\
//  Purpose:\
//\
//    TUPLE_NEXT computes the next element of a tuple space.\
//\
//  Discussion:\
//\
//    The elements are N vectors.  Each entry is constrained to lie\
//    between M1 and M2.  The elements are produced one at a time.\
//    The first element is\
//      (M1,M1,...,M1),\
//    the second element is\
//      (M1,M1,...,M1+1),\
//    and the last element is\
//      (M2,M2,...,M2)\
//    Intermediate elements are produced in lexicographic order.\
//\
//  Example:\
//\
//    N = 2, M1 = 1, M2 = 3\
//\
//    INPUT        OUTPUT\
//    -------      -------\
//    Rank  X      Rank   X\
//    ----  ---    -----  ---\
//    0     * *    1      1 1\
//    1     1 1    2      1 2\
//    2     1 2    3      1 3\
//    3     1 3    4      2 1\
//    4     2 1    5      2 2\
//    5     2 2    6      2 3\
//    6     2 3    7      3 1\
//    7     3 1    8      3 2\
//    8     3 2    9      3 3\
//    9     3 3    0      0 0\
//\
//  Modified:\
//\
//    29 April 2003\
//\
//  Author:\
//\
//    John Burkardt\
//\
//  Parameters:\
//\
//    Input, int M1, M2, the minimum and maximum entries.\
//\
//    Input, int N, the number of components.\
//\
//    Input/output, int *RANK, counts the elements.\
//    On first call, set RANK to 0.  Thereafter, the output value of RANK\
//    will indicate the order of the element returned.  When there are no\
//    more elements, RANK will be returned as 0.\
//\
//    Input/output, int X[N], on input the previous tuple.\
//    On output, the next tuple.\
//\
\{\
  int i;\
  int j;\
\
  if ( m2 < m1 )\
  \{\
    *rank = 0;\
    return;\
  \}\
\
  if ( *rank <= 0 )\
  \{\
    for ( i = 0; i < n; i++ )\
    \{\
      x[i] = m1;\
    \}\
    *rank = 1;\
  \}\
  else\
  \{\
    *rank = *rank + 1;\
    i = n - 1;\
\
    for ( ; ; )\
    \{\
\
      if ( x[i] < m2 )\
      \{\
        x[i] = x[i] + 1;\
        break;\
      \}\
\
      x[i] = m1;\
\
      if ( i == 0 )\
      \{\
        *rank = 0;\
        for ( j = 0; j < n; j++ )\
        \{\
          x[j] = m1;\
        \}\
        break;\
      \}\
      i = i - 1;\
    \}\
  \}\
  return;\
\}}