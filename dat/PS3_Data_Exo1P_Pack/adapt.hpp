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
void adapt ( int ndim, double a[], double b[], int *minpts, int maxpts, \
  double functn ( int indx, int ndim, double z[], double alpha[], \
  double beta[] ), double rel_tol, int itest, double alpha[], double beta[],\
  int lenwrk, double wrkstr[], double *relerr, double *finest, int *ifail )\
\
//****************************************************************************80\
//\
//  Purpose:\
//\
//    ADAPT carries out adaptive multidimensional quadrature.\
//\
//  Modified:\
//\
//    22 March 2007\
//\
//  Author:\
//\
//    Original FORTRAN77 version by Alan Genz\
//\
//  Parameters:\
//\
//    Input, int NDIM, the number of variables.\
//    2 <= NDIM.\
//\
//    Input, double A[NDIM], the lower limits of integration.\
//\
//    Input, double B[NDIM], the upper limits of integration.\
//\
//    Input/output, int *MINPTS, the minimum number of function evaluations\
//    to be allowed,  MINPTS must not exceed MAXPTS.  If MINPTS < 0 then the\
//    routine assumes a previous call has been made with the same integrand\
//    and continues that calculation.\
//\
//    Input, int MAXPTS, the maximum number of function\
//    evaluations allowed, which must be at least RULCLS, where\
//    RULCLS = 2**NDIM + 2 * NDIM**2 + 2 * NDIM + 1, when NDIM <= 15 and\
//    RULCLS = ( NDIM * ( 14 - NDIM * ( 6 - 4 * NDIM ) ) ) / 3 + 1,\
//    when 15 < NDIM.\
//    for NDIM  =  2   3   4   5   6   7   8   9   10   11   12\
//    RULCLS   =  17  33  57  93 149 241 401 693 1245 2313 4409\
//    A suggested starting value for MAXPTS is 100*RULCLS.  If\
//    this is not large enough for the required accuracy, then\
//    MAXPTS and LENWRK should be increased accordingly.\
//\
//    Input, external, double FUNCTN, the user-defined function\
//    to be integrated.  It must have the form\
//      double functn ( int indx, ind ntim, double z[], double alpha[],\
//        double beta[] )\
//    where\
//      INDX is the index of the test function,\
//      NDIM is the spatial dimension,\
//      Z is the evaluation point,\
//      ALPHA is a parameter vector,\
//      BETA is a parameter vector.\
//\
//    Input, double REL_TOL, the user's requested relative accuracy.\
//\
//    Input, int ITEST, the index of the test.\
//\
//    Input, double ALPHA[NDIM], BETA[NDIM], parameters\
//    associated with the integrand function.\
//\
//    Input, int LENWRK, the length of the array WRKSTR.\
//    The routine needs (2*NDIM+3)*(1+MAXPTS/RULCLS)/2 for LENWRK if\
//    MAXPTS function calls are used.\
//\
//    Input/output, double WRKSTR[LENWRK].  This array does not\
//    need to be set or inspected by the user.  However, the output value of\
//    WKRSTR from one call may be needed by the program on a followup call\
//    if the input value of MINPTS < 0, which signals that another calculation\
//    is requested for the same integrand.\
//\
//    Output, double *RELERR, the estimated relative accuracy\
//    of the integral estimate.\
//\
//    Output, double *FINEST, the estimated value of integral.\
//\
//    Output, int *IFAIL\
//    * 0, for normal exit, when estimated relative error RELERR is less\
//    than REL_TOL, and with MAXPTS or less function calls made.\
//    * 1, if MAXPTS was too small for ADAPT to obtain the required relative\
//    error REL_TOL.  In this case ADAPT returns a value of FINEST with\
//    estimated relative error RELERR.\
//    * 2, if LENWRK was too small for MAXPTS function calls.  In\
//    this case ADAPT returns a value of FINEST with estimated error\
//    RELERR using the working storage available, but RELERR is likely to\
//    be greater than REL_TOL.\
//    * 3, if NDIM < 2 or MAXPTS < MINPTS or MAXPTS < RULCLS.\
//\
\{\
  double *center;\
  double df1;\
  double df2;\
  double dif;\
  double difmax;\
  int divaxn;\
  int divaxo;\
  int divflg;\
  double f1;\
  double f2;\
  double f3;\
  double f4;\
  int funcls;\
  int i;\
  int index1;\
  int index2;\
  int j;\
  int k;\
  int l;\
  double lambda2;\
  double lambda4;\
  double lambda5;\
  int m;\
  int n;\
  double ratio;\
  double rgncmp;\
  double rgnerr;\
  int rgnstr = 0;\
  double rgnval;\
  double rgnvol;\
  int rulcls;\
  int sbrgns = 0;\
  int sbtmpp;\
  int subrgn;\
  int subtmp;\
  double sum1;\
  double sum2;\
  double sum3;\
  double sum4;\
  double sum5;\
  double weit1;\
  double weit2;\
  double weit3;\
  double weit4;\
  double weit5;\
  double weitp1;\
  double weitp2;\
  double weitp3;\
  double weitp4;\
  double *width;\
  double *widthl;\
  double *z;\
\
  *ifail = 3;\
  *relerr = 1.0;\
  funcls = 0;\
\
  if ( ndim < 2 )\
  \{\
    *minpts = 0;\
    wrkstr[lenwrk-2] = sbrgns;\
    *relerr = 1.0;\
    *finest = 0.0;\
    *ifail = 3;\
    return;\
  \}\
\
  if ( maxpts < *minpts )\
  \{\
    *minpts = 0;\
    wrkstr[lenwrk-2] = sbrgns;\
    *relerr = 1.0;\
    *finest = 0.0;\
    *ifail = 3;\
    return;\
  \}\
\
  if ( ndim <= 15 )\
  \{\
    rulcls = i4_power ( 2, ndim ) + 2 * ndim * ndim + 2 * ndim + 1;\
  \}\
  else if ( 15 < ndim )\
  \{\
    rulcls = 1 + ( ndim * ( 12 + ( ndim - 1 ) \
      * ( 6 + ( ndim - 2 ) * 4 ) ) ) / 3;\
  \}\
\
  if ( maxpts < rulcls )\
  \{\
    *relerr = 1.0;\
    *finest = 0.0;\
    *ifail = 3;\
    return;\
  \}\
//\
//  Initialization.\
//\
  rgnstr = 2 * ndim + 3;\
  divaxo = 0;\
\
  center = new double[ndim];\
  width = new double[ndim];\
  widthl = new double[ndim];\
  z = new double[ndim];\
//\
//  Basic rule initialization.\
//\
  lambda5 = 9.0 / 19.0;\
\
  if ( ndim <= 15 )\
  \{\
    lambda4 = 9.0 / 10.0;\
    lambda2 = 9.0 / 70.0;\
    weit5 = 1.0 / pow ( 3.0 * lambda5, 3 ) / pow ( 2.0, ndim );\
  \}\
  else\
  \{\
    ratio = ( double ) ( ndim - 2 ) / 9.0;\
\
    lambda4 = ( 1.0 / 5.0 - ratio ) / ( 1.0 / 3.0 - ratio / lambda5 );\
\
    ratio = ( 1.0 - lambda4 / lambda5 ) \
      * ( double ) ( ndim - 1 ) * ratio / 6.0;\
\
    lambda2 = ( 1.0 / 7.0 - lambda4 / 5.0 - ratio ) \
      / ( 1.0 / 5.0 - lambda4 / 3.0 - ratio / lambda5 );\
\
    weit5 = 1.0 / pow ( 6.0 * lambda5, 3 );\
  \}\
\
  weit4 = ( 1.0 / 15.0 - lambda5 / 9.0 ) \
    / ( 4.0 * ( lambda4 - lambda5 ) * lambda4 * lambda4 );\
\
  weit3 = ( 1.0 / 7.0 - ( lambda5 + lambda2 ) / 5.0 \
    + lambda5 * lambda2 / 3.0 ) / ( 2.0 * lambda4 \
    * ( lambda4 - lambda5 ) * ( lambda4 - lambda2 ) ) \
    - 2.0 * ( double ) ( ndim - 1 ) * weit4;\
\
  weit2 = ( 1.0 / 7.0 - ( lambda5 + lambda4 ) / 5.0 \
    + lambda5 * lambda4 / 3.0 ) / ( 2.0 * lambda2 \
    * ( lambda2 - lambda5 ) * ( lambda2 - lambda4 ) );\
\
  if ( ndim <= 15 )\
  \{\
    weit1 = 1.0 - 2.0 * ( double ) ( ndim )\
      * ( weit2 + weit3 + ( double ) ( ndim - 1 ) * weit4 ) \
      - pow ( 2.0, ndim ) * weit5;\
  \}\
  else\
  \{\
    weit1 = 1.0 - 2.0 * ( double ) ndim \
      * ( weit2 + weit3 + ( double ) ( ndim - 1 ) * \
      ( weit4 + 2.0 * ( double ) ( ndim - 2 ) * weit5 / 3.0 ) );\
  \}\
\
  weitp4 = 1.0 / pow ( 6.0 * lambda4, 2 );\
\
  weitp3 = ( 1.0 / 5.0 - lambda2 / 3.0 ) / \
    ( 2.0 * lambda4 * ( lambda4 - lambda2 ) ) \
    - 2.0 * ( double ) ( ndim - 1 ) * weitp4;\
\
  weitp2 = ( 1.0 / 5.0 - lambda4 / 3.0 ) \
    / ( 2.0 * lambda2 * ( lambda2 - lambda4 ) );\
\
  weitp1 = 1.0 - 2.0 * ( double ) ( ndim ) * \
    ( weitp2 + weitp3 + ( double ) ( ndim - 1 ) * weitp4 );\
\
  ratio = lambda2 / lambda4;\
\
  lambda5 = sqrt ( lambda5 );\
  lambda4 = sqrt ( lambda4 );\
  lambda2 = sqrt ( lambda2 );\
//\
//  End basic rule initialization.\
//\
  if ( *minpts < 0 )\
  \{\
    sbrgns = ( int ) wrkstr[lenwrk-2];\
    divflg = 0;\
    subrgn = rgnstr;\
    wrkstr[lenwrk-1] = wrkstr[lenwrk-1] - wrkstr[subrgn-1];\
    *finest = *finest - wrkstr[subrgn-2];\
    divaxo = ( int ) wrkstr[subrgn-3];\
\
    for ( j = 1; j <= ndim; j++ )\
    \{\
      subtmp = subrgn - 2 * ( j + 1 );\
      center[j-1] = wrkstr[subtmp];\
      width[j-1] = wrkstr[subtmp-1];\
    \}\
    width[divaxo-1] = width[divaxo-1] / 2.0;\
    center[divaxo-1] = center[divaxo-1] - width[divaxo-1];\
  \}\
  else\
  \{\
    for ( j = 0; j < ndim; j++ )\
    \{\
      width[j] = ( b[j] - a[j] ) / 2.0;\
    \}\
    for ( j = 0; j < ndim; j++ )\
    \{\
      center[j] = a[j] + width[j];\
    \}\
\
    *finest = 0.0;\
    wrkstr[lenwrk-1] = 0.0;\
    divflg = 1;\
    subrgn = rgnstr;\
    sbrgns = rgnstr;\
  \}\
//\
//  Begin basic rule.\
//\
  for ( ; ; )\
  \{\
    rgnvol = pow ( 2.0, ndim ) * r8vec_product ( ndim, width );\
\
    for ( j = 0; j < ndim; j++ )\
    \{\
      z[j] = center[j];\
    \}\
\
    sum1 = functn ( itest, ndim, z, alpha, beta );\
//\
//  Compute symmetric sums of functn(lambda2,0,0,...,0) and\
//  functn(lambda4,0,0,...,0), and maximum fourth difference.\
//\
    difmax = -1.0;\
    sum2 = 0.0;\
    sum3 = 0.0;\
\
    for ( j = 0; j < ndim; j++ )\
    \{\
      z[j] = center[j] - lambda2 * width[j];\
      f1 = functn ( itest, ndim, z, alpha, beta );\
      z[j] = center[j] + lambda2 * width[j];\
      f2 = functn ( itest, ndim, z, alpha, beta );\
      widthl[j] = lambda4 * width[j];\
      z[j] = center[j] - widthl[j];\
      f3 = functn ( itest, ndim, z, alpha, beta );\
      z[j] = center[j] + widthl[j];\
      f4 = functn ( itest, ndim, z, alpha, beta );\
      sum2 = sum2 + f1 + f2;\
      sum3 = sum3 + f3 + f4;\
      df1 = f1 + f2 - 2.0 * sum1;\
      df2 = f3 + f4 - 2.0 * sum1;\
      dif = r8_abs ( df1 - ratio * df2 );\
\
      if ( difmax < dif )\
      \{\
        difmax = dif;\
        divaxn = j + 1;\
      \}\
      z[j] = center[j];\
    \}\
\
    if ( sum1 == sum1 + difmax / 8.0 )\
    \{\
      divaxn = ( divaxo % ndim ) + 1;\
    \}\
//\
//  Compute symmetric sum of functn(lambda4,lambda4,0,0,...,0).\
//\
    sum4 = 0.0;\
\
    for ( j = 2; j <= ndim; j++ )\
    \{\
      for ( k = j; k <= ndim; k++ )\
      \{\
        for ( l = 1; l <= 2; l++ )\
        \{\
          widthl[j-2] = -widthl[j-2];\
          z[j-2] = center[j-2] + widthl[j-2];\
          for ( m = 1; m <= 2; m++ )\
          \{\
            widthl[k-1] = -widthl[k-1];\
            z[k-1] = center[k-1] + widthl[k-1];\
            sum4 = sum4 + functn ( itest, ndim, z, alpha, beta );\
          \}\
        \}\
        z[k-1] = center[k-1];\
      \}\
      z[j-2] = center[j-2];\
    \}\
//\
//  If NDIM < 16 compute symmetric sum of functn(lambda5,lambda5,...,lambda5).\
//\
    if ( ndim <= 15 )\
    \{\
      sum5 = 0.0;\
\
      for ( j = 0; j < ndim; j++ )\
      \{\
        widthl[j] = -lambda5 * width[j];\
      \}\
      for ( j = 0; j < ndim; j++ )\
      \{\
        z[j] = center[j] + widthl[j];\
      \}\
\
      for ( ; ; )\
      \{\
        sum5 = sum5 + functn ( itest, ndim, z, alpha, beta );\
\
        j = ndim;\
\
        for ( ; ; )\
        \{\
          widthl[j-1] = - widthl[j-1];\
          z[j-1] = center[j-1] + widthl[j-1];\
\
          if ( 0.0 <= widthl[j-1] )\
          \{\
            break;\
          \}\
          j = j - 1;\
\
          if ( j < 1 )\
          \{\
            break;\
          \}\
        \}\
\
        if ( j < 1 )\
        \{\
          break;\
        \}\
      \}\
    \}\
//\
//  If 15 < NDIM, compute symmetric sum of\
//  FUNCTN(lambda5,lambda5,lambda5,0,0,...,0).\
//\
    else\
    \{\
      sum5 = 0.0;\
\
\
      for ( j = 0; j < ndim; j++ )\
      \{\
        widthl[j] = lambda5 * width[j];\
      \}\
\
      for ( i = 3; i <= ndim; i++ )\
      \{\
        for ( j = i; j <= ndim; j++ )\
        \{\
          for ( k = j; k <= ndim; k++ )\
          \{\
            for ( l = 1; l <= 2; l++ )\
            \{\
              widthl[i-3] = -widthl[i-3];\
              z[i-3] = center[i-3] + widthl[i-3];\
              for ( m = 1; m <= 2; m++ )\
              \{\
                widthl[j-2] = -widthl[j-2];\
                z[j-2] = center[j-2] + widthl[j-2];\
                for ( n = 1; n <= 2; n++ )\
                \{\
                  widthl[k-1] = -widthl[k-1];\
                  z[k-1] = center[k-1] + widthl[k-1];\
                  sum5 = sum5 + functn ( itest, ndim, z, alpha, beta );\
                \}\
              \}\
            \}\
            z[k-1] = center[k-1];\
          \}\
          z[j-2] = center[j-2];\
        \}\
        z[i-3] = center[i-3];\
      \}\
    \}\
//\
//  Compute fifth and seventh degree rules and error.\
//\
    rgncmp = rgnvol * ( weitp1 * sum1 \
                      + weitp2 * sum2 \
                      + weitp3 * sum3 \
                      + weitp4 * sum4 );\
\
    rgnval = rgnvol * ( weit1 * sum1 \
                      + weit2 * sum2 \
                      + weit3 * sum3 \
                      + weit4 * sum4 \
                      + weit5 * sum5 );\
\
    rgnerr = r8_abs ( rgnval - rgncmp );\
//\
//  End basic rule.\
//\
    *finest = *finest + rgnval;\
    wrkstr[lenwrk-1] = wrkstr[lenwrk-1] + rgnerr;\
    funcls = funcls + rulcls;\
//\
//  Place results of basic rule into partially ordered list\
//  according to subregion error.\
//\
//  When DIVFLG = 0, start at the top of the list and move down the\
//  list tree to find the correct position for the results from the\
//  first half of the recently divided subregion.\
//\
    if ( divflg != 1 )\
    \{\
      for ( ; ; )\
      \{\
        subtmp = 2 * subrgn;\
        if ( sbrgns < subtmp )\
        \{\
          break;\
        \}\
        if ( subtmp != sbrgns )\
        \{\
          sbtmpp = subtmp + rgnstr;\
          if ( wrkstr[subtmp-1] < wrkstr[sbtmpp-1] )\
          \{\
            subtmp = sbtmpp;\
          \}\
        \}\
        if ( wrkstr[subtmp-1] <= rgnerr )\
        \{\
          break;\
        \}\
        for ( k = 1; k <= rgnstr; k++ )\
        \{\
          wrkstr[subrgn-k] = wrkstr[subtmp-k];\
        \}\
        subrgn = subtmp;\
      \}\
    \}\
//\
//  When DIVFLG = 1 start at bottom right branch and move up list\
//  tree to find correct position for results from second half of\
//  recently divided subregion.\
//\
    else\
    \{\
      for ( ; ; )\
      \{\
        subtmp = ( subrgn / ( 2 * rgnstr ) ) * rgnstr;\
\
        if ( subtmp < rgnstr )\
        \{\
          break;\
        \}\
        if ( rgnerr <= wrkstr[subtmp-1] )\
        \{\
          break;\
        \}\
        for ( k = 1; k <= rgnstr; k++ )\
        \{\
          index1 = subrgn - k + 1;\
          index2 = subtmp - k + 1;\
          wrkstr[index1-1] = wrkstr[index2-1];\
        \}\
        subrgn = subtmp;\
      \}\
    \}\
//\
//  Store results of basic rule in correct position in list.\
//\
    wrkstr[subrgn-1] = rgnerr;\
    wrkstr[subrgn-2] = rgnval;\
    wrkstr[subrgn-3] = divaxn;\
\
    for ( j = 1; j <= ndim; j++ )\
    \{\
      subtmp = subrgn - 2 * ( j + 1 );\
      wrkstr[subtmp] = center[j-1];\
      wrkstr[subtmp-1] = width[j-1];\
    \}\
//\
//  When DIVFLG = 0 prepare for second application of basic rule.\
//\
    if ( divflg != 1 )\
    \{\
      center[divaxo-1] = center[divaxo-1] + 2.0 * width[divaxo-1];\
      sbrgns = sbrgns + rgnstr;\
      subrgn = sbrgns;\
      divflg = 1;\
      continue;\
    \}\
//\
//  End ordering and storage of basic rule results.\
//  Make checks for possible termination of routine.\
//\
    *relerr = 1.0;\
\
    if ( wrkstr[lenwrk-1] <= 0.0 )\
    \{\
      wrkstr[lenwrk-1] = 0.0;\
    \}\
\
    if ( r8_abs ( *finest ) != 0.0 )\
    \{\
      *relerr = wrkstr[lenwrk-1] / r8_abs ( *finest );\
    \}\
\
    if ( 1.0 < *relerr )\
    \{\
      *relerr = 1.0;\
    \}\
\
    if ( lenwrk < sbrgns + rgnstr + 2 )\
    \{\
      *ifail = 2;\
    \}\
\
    if ( maxpts < funcls + 2 * rulcls )\
    \{\
      *ifail = 1;\
    \}\
\
    if ( *relerr < rel_tol && *minpts <= funcls )\
    \{\
      *ifail = 0;\
    \}\
\
    if ( *ifail < 3 )\
    \{\
      *minpts = funcls;\
      wrkstr[lenwrk-2] = sbrgns;\
      break;\
    \}\
//\
//  Prepare to use basic rule on each half of subregion with largest\
//  error.\
//\
    divflg = 0;\
    subrgn = rgnstr;\
    wrkstr[lenwrk-1] = wrkstr[lenwrk-1] - wrkstr[subrgn-1];\
    *finest = *finest - wrkstr[subrgn-2];\
    divaxo = ( int ) wrkstr[subrgn-3];\
\
    for ( j = 1; j <= ndim; j++ )\
    \{\
      subtmp = subrgn - 2 * ( j + 1 );\
      center[j-1] = wrkstr[subtmp];\
      width[j-1] = wrkstr[subtmp-1];\
    \}\
\
    width[divaxo-1] = width[divaxo-1] / 2.0;\
    center[divaxo-1] = center[divaxo-1] - width[divaxo-1];\
  \}\
\
  delete [] center;\
  delete [] width;\
  delete [] widthl;\
  delete [] z;\
\
  return;\
\}\
//****************************************************************************80}