# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "simplex_gm_rule.hpp"

//****************************************************************************80

void comp_next ( int n, int k, int a[], bool *more, int *h, int *t )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//    This routine originally used a STATIC statement to maintain the
//    variables H and T.  I have decided (based on an wasting an
//    entire morning trying to track down a problem) that it is safer
//    to pass these variables as arguments, even though the user should
//    never alter them.  This allows this routine to safely shuffle
//    between several ongoing calculations.
//
//  Example:
//
//    The 28 compositions of 6 into three parts are:
//
//      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
//      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
//      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
//      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
//      0 3 3,  0 2 4,  0 1 5,  0 0 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf,
//    C++ translation by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool *MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them//
//
{
  int i;

  if ( ! ( *more ) )
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < *t )
    {
      *h = 0;
    }

    *h = *h + 1;
    *t = a[*h-1];
    a[*h-1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;

  }

  *more = ( a[k-1] != n );

  return;
}
//****************************************************************************80

void gm_rule_set ( int rule, int dim_num, int point_num, double w[],
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    GM_RULE_SET sets a Grundmann-Moeller rule.
//
//  Discussion:
//
//    This is a revised version of the calculation which seeks to compute
//    the value of the weight in a cautious way that avoids intermediate
//    overflow.  Thanks to John Peterson for pointing out the problem on
//    26 June 2008.
//
//    This rule returns weights and abscissas of a Grundmann-Moeller
//    quadrature rule for the DIM_NUM-dimensional unit simplex.
//
//    The dimension POINT_NUM can be determined by calling GM_RULE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Axel Grundmann, Michael Moeller,
//    Invariant Integration Formulas for the N-Simplex
//    by Combinatorial Methods,
//    SIAM Journal on Numerical Analysis,
//    Volume 15, Number 2, April 1978, pages 282-290.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//    0 <= RULE.
//
//    Input, int DIM_NUM, the spatial dimension.
//    1 <= DIM_NUM.
//
//    Input, int POINT_NUM, the number of points in the rule.
//
//    Output, double W[POINT_NUM], the weights.
//
//    Output, double X[DIM_NUM*POINT_NUM], the abscissas.
//
{
  int *beta;
  int beta_sum;
  int d;
  int dim;
  int h;
  int i;
  int j;
  int j_hi;
  int k;
  bool more;
  int n;
  double one_pm;
  int s;
  int t;
  double weight;

  s = rule;
  d = 2 * s + 1;
  k = 0;
  n = dim_num;
  one_pm = 1.0;

  beta = new int[dim_num+1];

  for ( i = 0; i <= s; i++ )
  {
    weight = ( double ) one_pm;

    j_hi = i4_max ( n, i4_max ( d, d + n - i ) );

    for ( j = 1; j <= j_hi; j++ )
    {
      if ( j <= n )
      {
        weight = weight * ( double ) ( j );
      }
      if ( j <= d )
      {
        weight = weight * ( double ) ( d + n - 2 * i );
      }
      if ( j <= 2 * s )
      {
        weight = weight / 2.0;
      }
      if ( j <= i )
      {
        weight = weight / ( double ) ( j );
      }
      if ( j <= d + n - i )
      {
        weight = weight / ( double ) ( j );
      }
    }

    one_pm = - one_pm;

    beta_sum = s - i;
    more = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      comp_next ( beta_sum, dim_num + 1, beta, &more, &h, &t );

      w[k] = weight;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        x[dim+k*dim_num] = ( double ) ( 2 * beta[dim+1] + 1 )
                         / ( double ) ( d + n - 2 * i );
      }
      k = k + 1;

      if ( !more )
      {
        break;
      }
    }
  }

  delete [] beta;

  return;
}
//****************************************************************************80

void gm_rule_set_old ( int rule, int dim_num, int point_num, double w[],
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    GM_RULE_SET_OLD sets a Grundmann-Moeller rule.  (OBSOLETE VERSION)
//
//  Discussion:
//
//    This version of the computation is no longer used.  The direct
//    application of the formula results in overflows and inaccuracies
//    very quickly.
//
//    This rule returns weights and abscissas of a Grundmann-Moeller
//    quadrature rule for the DIM_NUM-dimensional unit simplex.
//
//    The dimension POINT_NUM can be determined by calling GM_RULE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Axel Grundmann, Michael Moeller,
//    Invariant Integration Formulas for the N-Simplex
//    by Combinatorial Methods,
//    SIAM Journal on Numerical Analysis,
//    Volume 15, Number 2, April 1978, pages 282-290.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//    0 <= RULE.
//
//    Input, int DIM_NUM, the spatial dimension.
//    1 <= DIM_NUM.
//
//    Input, int POINT_NUM, the number of points in the rule.
//
//    Output, double W[POINT_NUM], the weights.
//
//    Output, double X[DIM_NUM*POINT_NUM], the abscissas.
//
{
  int *beta;
  int beta_sum;
  int d;
  int dim;
  int h;
  int i;
  int k;
  bool more;
  int n;
  double one_pm;
  int s;
  int t;
  double weight;

  s = rule;
  d = 2 * s + 1;
  k = 0;
  n = dim_num;
  one_pm = 1.0;

  beta = new int[dim_num+1];

  for ( i = 0; i <= s; i++ )
  {
    weight = r8_factorial ( n ) * one_pm
      * pow ( ( double ) d + n - 2 * i, d )
      / pow ( 2.0, 2 * s )
      / r8_factorial ( i )
      / r8_factorial ( d + n - i );

    one_pm = - one_pm;

    beta_sum = s - i;
    more = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      comp_next ( beta_sum, dim_num + 1, beta, &more, &h, &t );

      w[k] = weight;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        x[dim+k*dim_num] = ( double ) ( 2 * beta[dim+1] + 1 )
                         / ( double ) ( d + n - 2 * i );
      }
      k = k + 1;

      if ( !more )
      {
        break;
      }
    }
  }

  delete [] beta;

  return;
}
//****************************************************************************80

int gm_rule_size ( int rule, int dim_num )

//****************************************************************************80
//
//  Purpose:
//
//    GM_RULE_SIZE determines the size of a Grundmann-Moeller rule.
//
//  Discussion:
//
//    This rule returns the value of POINT_NUM, the number of points associated
//    with a GM rule of given index.
//
//    After calling this rule, the user can use the value of POINT_NUM to
//    allocate space for the weight vector as W(POINT_NUM) and the abscissa
//    vector as X(DIM_NUM,POINT_NUM), and then call GM_RULE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Axel Grundmann, Michael Moeller,
//    Invariant Integration Formulas for the N-Simplex
//    by Combinatorial Methods,
//    SIAM Journal on Numerical Analysis,
//    Volume 15, Number 2, April 1978, pages 282-290.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//    0 <= RULE.
//
//    Input, int DIM_NUM, the spatial dimension.
//    1 <= DIM_NUM.
//
//    Output, int GM_RULE_SIZE, the number of points in the rule.
//
{
  int arg1;
  int point_num;

  arg1 = dim_num + rule + 1;

  point_num = i4_choose ( arg1, rule );

  return point_num;
}
//****************************************************************************80

int i4_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE computes the binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, are the values of N and K.
//
//    Output, int I4_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0;
  }
  else if ( mn == 0 )
  {
    value = 1;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
    }
  }

  return value;
}
//****************************************************************************80

int i4_huge ( void )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" integer.
//
{
  return 2147483647;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of i1 and i2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

double *monomial_value ( int dim_num, int point_num, double x[], int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points at which the
//    monomial is to be evaluated.
//
//    Input, double X[DIM_NUM*POINT_NUM], the point coordinates.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Output, double MONOMIAL_VALUE[POINT_NUM], the value of the monomial.
//
{
  int dim;
  int point;
  double *value;

  value = new double[point_num];

  for ( point = 0; point < point_num; point++ )
  {
    value[point] = 1.0;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0 != expon[dim] )
    {
      for ( point = 0; point < point_num; point++ )
      {
        value[point] = value[point] * pow ( x[dim+point*dim_num], expon[dim] );
      }
    }
  }

  return value;
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

double *r8vec_uniform_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
//
{
  int i;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error//\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge ( );
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double simplex_unit_monomial_int ( int dim_num, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_MONOMIAL_INT integrates a monomial over a simplex.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Output, double SIMPLEX_UNIT_MONOMIAL_INT, the value of the integral
//    of the monomial.
//
{
  int dim;
  int i;
  int k;
  double value;

  value = 1.0;

  k = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    for ( i = 1; i <= expon[dim]; i++ )
    {
      k = k + 1;
      value = value * ( double ) ( i ) / ( double ) ( k );
    }
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    k = k + 1;
    value = value / ( double ) ( k );
  }

  return value;
}
//****************************************************************************80

double simplex_unit_monomial_quadrature ( int dim_num, int expon[],
  int point_num, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_MONOMIAL_QUADRATURE: quadrature of monomials in a unit simplex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Input, int POINT_NUM, the number of points in the rule.
//
//    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
//
//    Input, double W[POINT_NUM], the quadrature weights.
//
//    Output, double SIMPLEX_UNIT_MONOMIAL_QUADRATURE, the quadrature error.
//
{
  double exact = 1.0;
  double quad;
  double quad_error;
  double scale;
  double *value;
  double volume;
//
//  Get the exact value of the integral of the unscaled monomial.
//
  scale = simplex_unit_monomial_int ( dim_num, expon );
//
//  Evaluate the monomial at the quadrature points.
//
  value = monomial_value ( dim_num, point_num, x, expon );
//
//  Compute the weighted sum and divide by the exact value.
//
  volume = simplex_unit_volume ( dim_num );
  quad = volume * r8vec_dot_product ( point_num, w, value ) / scale;

  delete [] value;
//
//  Error:
//
  quad_error = fabs ( quad - exact );

  return quad_error;
}
//****************************************************************************80

double *simplex_unit_sample ( int dim_num, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_SAMPLE returns uniformly random points from a general simplex.
//
//  Discussion:
//
//    The interior of the unit DIM_NUM dimensional simplex is the set of 
//    points X(1:DIM_NUM) such that each X(I) is nonnegative, and 
//    sum(X(1:DIM_NUM)) <= 1.
//
//    This routine is valid for any spatial dimension DIM_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity
//    of Queueing Networks,
//    Krieger, 1992,
//    ISBN: 0894647644,
//    LC: QA298.R79.
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_SIMPLEX01_MAP[DIM_NUM*N], the points.
//
{
  double *e;
  int i;
  int j;
  double total;
  double *x;
//
//  The construction begins by sampling DIM_NUM+1 points from the
//  exponential distribution with parameter 1.
//
  x = new double[dim_num*n];

  for ( j = 0; j < n; j++ )
  {
    e = r8vec_uniform_01 ( dim_num+1, seed );

    for ( i = 0; i <= dim_num; i++ )
    {
      e[i] = -log ( e[i] );
    }

    total = 0.0;
    for ( i = 0; i <= dim_num; i++ )
    {
      total = total + e[i];
    }

    for ( i = 0; i < dim_num; i++ )
    {
      x[i+dim_num*j] = e[i] / total;
    }
    delete [] e;
  }

  return x;
}
//****************************************************************************80

double *simplex_unit_to_general ( int dim_num, int point_num, double t[],
  double ref[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_TO_GENERAL maps the unit simplex to a general simplex.
//
//  Discussion:
//
//    Given that the unit simplex has been mapped to a general simplex
//    with vertices T, compute the images in T, under the same linear
//    mapping, of points whose coordinates in the unit simplex are REF.
//
//    The vertices of the unit simplex are listed as suggested in the
//    following:
//
//      (0,0,0,...,0)
//      (1,0,0,...,0)
//      (0,1,0,...,0)
//      (0,0,1,...,0)
//      (...........)
//      (0,0,0,...,1)
//
//    Thanks to Andrei ("spiritualworlds") for pointing out a mistake in the
//    previous implementation of this routine, 02 March 2008.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points to transform.
//
//    Input, double T[DIM_NUM*(DIM_NUM+1)], the vertices of the
//    general simplex.
//
//    Input, double REF[DIM_NUM*POINT_NUM], points in the
//    reference triangle.
//
//    Output, double SIMPLEX_UNIT_TO_GENERAL[DIM_NUM*POINT_NUM],
//    corresponding points in the physical triangle.
//
{
  int dim;
  double *phy;
  int point;
  int vertex;

  phy = new double[dim_num*point_num];
//
//  The image of each point is initially the image of the origin.
//
//  Insofar as the pre-image differs from the origin in a given vertex
//  direction, add that proportion of the difference between the images
//  of the origin and the vertex.
//
  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      phy[dim+point*dim_num] = t[dim+0*dim_num];

      for ( vertex = 1; vertex < dim_num + 1; vertex++ )
      {
        phy[dim+point*dim_num] = phy[dim+point*dim_num]
        + ( t[dim+vertex*dim_num] - t[dim+0*dim_num] ) * ref[vertex-1+point*dim_num];
      }
    }
  }

  return phy;
}
//****************************************************************************80

double simplex_unit_volume ( int dim_num )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_VOLUME computes the volume of the unit simplex.
//
//  Discussion:
//
//    The formula is simple: volume = 1/DIM_NUM!.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Output, double SIMPLEX_UNIT_VOLUME, the volume of the cone.
//
{
  int i;
  double volume;

  volume = 1.0;
  for ( i = 1; i <= dim_num; i++ )
  {
    volume = volume / ( ( double ) i );
  }

  return volume;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
