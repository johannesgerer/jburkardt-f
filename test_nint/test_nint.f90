function csevl ( x, cs, n )

!*****************************************************************************80
!
!! CSEVL evaluates an N term Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, 1973, pages 254-256.
!
!    Leslie Fox, Ian Parker,
!    Chebyshev Polynomials in Numerical Analysis,
!    Oxford Press, 1968,
!    LC: QA297.F65.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument at which the series is 
!    to be evaluated.
!    X must satisfy -1.0 <= X <= 1.0.
!
!    Input, real ( kind = 8 ) CS(N), the array of N terms of a Chebyshev series.
!    In evaluating CS, only half the first coefficient is summed.
!
!    Input, integer ( kind = 4 ) N, the number of terms in array CS.
!    N must be at least 1, and no more than 1000.
!
!    Output, real ( kind = 8 ) CSEVL, the value of the Chebyshev series.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) cs(n)
  real ( kind = 8 ) csevl
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms N is less than 1.'
    stop
  end if

  if ( 1000 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  The number of terms is more than 1000.'
    stop
  end if

  if ( x < -1.0D+00 .or. 1.0D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  The input argument X is outside the interval [-1,1].'
    stop
  end if

  b1 = 0.0D+00
  b0 = 0.0D+00

  do i = n, 1, -1
    b2 = b1
    b1 = b0
    b0 = 2.0D+00 * x * b1 - b2 + cs(i)
  end do

  csevl = 0.5D+00 * ( b0 - b2 )

  return
end
subroutine get_problem_num ( problem_num )

!*****************************************************************************80
!
!! GET_PROBLEM_NUM returns the number of test integration problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROBLEM_NUM, the number of test problems.
!
  implicit none

  integer ( kind = 4 ) problem_num

  problem_num = 33

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) /  11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) /  30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) /  23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp                                    /   6.0D+00

  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( i4_huge ( ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == i4_huge ( ) ) then
    seed = seed - 1
  end if

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
function inits ( os, nos, eta )

!*****************************************************************************80
!
!! INITS estimates the order of an orthogonal series for a given accuracy.
!
!  Discussion:
!
!    Because this is a function, it is not possible to print out
!    warning error messages.  Therefore, if an error condition is
!    detected, a bogus value of INITS is returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, 1973, pages 254-256.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) OS(NOS), the coefficients in the series.
!
!    Input, integer ( kind = 4 ) NOS, the number of coefficients.  NOS must be
!    at least 1, or an error condition arises.
!
!    Input, real ( kind = 8 ) ETA, the requested accuracy of the series.
!    Ordinarily, ETA will be chosen to be one-tenth machine precision.
!
!    Output, integer ( kind = 4 ) INITS, the order of the series guaranteeing 
!    the given accuracy.  However, on error, INITS will be returned
!    as a negative number.
!
  implicit none

  integer ( kind = 4 ) nos

  real ( kind = 8 ) err
  real ( kind = 8 ) eta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inits
  real ( kind = 8 ) os(nos)
!
!  Fatal error.  Number of coefficients less than 1.
!  But an error message cannot be printed out because this is a function.
!
  if ( nos < 1 ) then
    inits = - huge ( i )
    return
  end if

  err = 0.0D+00

  i = 0

  do ii = 1, nos

    i = nos + 1 - ii
    err = err + abs ( os(i) )

    if ( eta < err ) then
      i = nos + 1 - ii
      exit
    end if

  end do
!
!  Eta may be too small.  Accuracy cannot be guaranteed.
!  But an error message cannot be printed out because this is a function.
!
  if ( i == 0 ) then
    i = - nos
  end if

  inits = i

  return
end
subroutine legendre_compute ( order, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE computes a Gauss-Legendre quadrature rule.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be greater than 0.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2pn
  real ( kind = 8 ) d3pn
  real ( kind = 8 ) d4pn
  real ( kind = 8 ) dp
  real ( kind = 8 ) dpn
  real ( kind = 8 ) e1
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iback
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp1mi
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) nmove
  real ( kind = 8 ) p
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pk
  real ( kind = 8 ) pkm1
  real ( kind = 8 ) pkp1
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x0
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xtemp
  real ( kind = 8 ) weight(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if

  e1 = real ( order * ( order + 1 ), kind = 8 )

  m = ( order + 1 ) / 2

  do i = 1, m

    mp1mi = m + 1 - i

    t = real ( 4 * i - 1, kind = 8 ) * pi &
      / real ( 4 * order + 2, kind = 8 )

    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( order, kind = 8 ) ) &
      / real ( 8 * order * order, kind = 8 ) )

    pkm1 = 1.0D+00
    pk = x0

    do k = 2, order
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do

    d1 = real ( order, kind = 8 ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0 * x0 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 * x0 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 * x0 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) / &
      ( 1.0D+00 - x0 * x0 )

    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )

    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )

    h = h - p / dp

    xtemp = x0 + h

    xtab(mp1mi) = xtemp

    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

    weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp * xtemp ) / ( fx * fx )

  end do

  if ( mod ( order, 2 ) == 1 ) then
    xtab(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = ( order + 1 ) / 2
  ncopy = order - nmove

  do i = 1, nmove
    iback = order + 1 - i
    xtab(iback) = xtab(iback-ncopy)
    weight(iback) = weight(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, order - nmove
    xtab(i) = - xtab(order+1-i)
    weight(i) = weight(order+1-i)
  end do

  return
end
subroutine normal_01_cdf_inv ( cdf, x )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the Normal 01 CDF.
!
!  Discussion:
!
!    Modified to handle case where R = 0 would otherwise require
!    evaluation of LOG(0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2007
!
!  Author:
!
!    Original FORTRAN77 version by JD Beasley, SG Springer,
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    JD Beasley, SG Springer,
!    Algorithm AS 111:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 26, 1977, pages 118-121.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) X, the corresponding argument.
!
  implicit none

  real ( kind = 8 ), parameter :: a0 =   2.50662823884D+00
  real ( kind = 8 ), parameter :: a1 = -18.61500062529D+00
  real ( kind = 8 ), parameter :: a2 =  41.39119773534D+00
  real ( kind = 8 ), parameter :: a3 = -25.44106049637D+00
  real ( kind = 8 ), parameter :: b1 =  -8.47351093090D+00
  real ( kind = 8 ), parameter :: b2 =  23.08336743743D+00
  real ( kind = 8 ), parameter :: b3 = -21.06224101826D+00
  real ( kind = 8 ), parameter :: b4 =   3.13082909833D+00
  real ( kind = 8 ), parameter :: c0 =  -2.78718931138D+00
  real ( kind = 8 ), parameter :: c1 =  -2.29796479134D+00
  real ( kind = 8 ), parameter :: c2 =   4.85014127135D+00
  real ( kind = 8 ), parameter :: c3 =   2.32121276858D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ), parameter :: d1 =   3.54388924762D+00
  real ( kind = 8 ), parameter :: d2 =   1.63706781897D+00
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NORMAL_01_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  q = cdf - 0.5D+00

  q = min ( q, 0.5D+00 )
  q = max ( q, -0.5D+00 )
!
!  0.08 < CDF < 0.92
!
  if ( abs ( q ) <= 0.42D+00 ) then

    r = q * q

    x = q * ( ( ( &
           a3   * r &
         + a2 ) * r &
         + a1 ) * r &
         + a0 ) / ( ( ( ( &
           b4   * r &
         + b3 ) * r &
         + b2 ) * r &
         + b1 ) * r + 1.0D+00 )
!
!  CDF < 0.08 or 0.92 < CDF.
!
  else

    r = min ( cdf, 1.0D+00 - cdf )
    r = max ( r, 0.0D+00 )
    r = min ( r, 1.0D+00 )

    if ( r <= 0.0D+00 ) then

      x = r8_huge ( )

    else

      r = sqrt ( - log ( r ) )

      x = ( ( ( &
             c3   * r &
           + c2 ) * r &
           + c1 ) * r &
           + c0 ) / ( ( &
             d2   * r &
           + d1 ) * r + 1.0D+00 )

    end if

    if ( q < 0.0D+00 ) then
      x = - x
    end if

  end if

  return
end
subroutine p00_box_gl ( problem, dim_num, order_1d, result )

!*****************************************************************************80
!
!! P00_BOX_GL uses Gauss-Legendre quadrature in an N-dimensional box. 
!
!  Discussion:
!
!    The rule is the product of a Gauss-Legendre rule in each spatial dimension.
!
!    Subdivision is NOT used.  This is a pure product rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ORDER_1D, the order of the 1D rule to use.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_1d
  integer ( kind = 4 ), parameter :: point_num = 1

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) indx(dim_num)
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) k
  real ( kind = 8 ) result
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) weight(order_1d)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ) xtab(order_1d)
!
!  Get the integration limits, and the weight adjustment factor.
!
  call p00_lim ( problem, dim_num, a, b )
!
!  Compute the 1D rule.
!
  call legendre_compute ( order_1d, xtab, weight )
!
!  Carry out the product rule.
!
  result = 0.0D+00

  k = 0

  do

    call tuple_next ( 1, order_1d, dim_num, k, indx )

    if ( k == 0  ) then
      exit
    end if

    w = product ( weight(indx(1:dim_num)) )

    x(1:dim_num,1) = 0.5D+00 * ( &
        a(1:dim_num) * ( 1.0D+00 - xtab(indx(1:dim_num)) ) &
      + b(1:dim_num) * ( 1.0D+00 + xtab(indx(1:dim_num)) ) )

    call p00_f ( problem, dim_num, point_num, x, value )

    result = result + w * value(1)

  end do
!
!  Get the volume.
!
  call p00_volume ( problem, dim_num, volume )

  result = result * volume / real ( 2.0D+00**dim_num, kind = 8 )

  return
end
subroutine p00_box_gl05 ( problem, dim_num, sub_num, result )

!*****************************************************************************80
!
!! P00_BOX_GL05 uses Gauss-Legendre quadrature in an N-dimensional box. 
!
!  Discussion:
!
!    The rule is a composite rule.  The region is divided into subregions,
!    over each of which a product rule is used, based on the 5-point 
!    Gauss-Legendre rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SUB_NUM, the number of subdivisions 
!    in each dimension.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), parameter :: norder = 5
  integer ( kind = 4 ), parameter :: point_num = 1

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) a_sub(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ) b_sub(dim_num)
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) indx(dim_num)
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) result
  integer ( kind = 4 ) sub_indx(dim_num)
  integer ( kind = 4 ) sub_num
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ), save, dimension ( norder ) :: weight = (/ &
    0.236926885056189087514264040720D+00, &
    0.478628670499366468041291514836D+00, &
    0.568888888888888888888888888889D+00, &
    0.478628670499366468041291514836D+00, &
    0.236926885056189087514264040720D+00 /)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ), save, dimension ( norder ) :: xtab = (/ &
    - 0.906179845938663992797626878299D+00, &
    - 0.538469310105683091036314420700D+00, &
      0.0D+00, &
      0.538469310105683091036314420700D+00, &
      0.906179845938663992797626878299D+00 /)
!
!  Get the integration limits, and the weight adjustment factor.
!
  call p00_lim ( problem, dim_num, a, b )
!
!  Carry out the product rule.
!
  result = 0.0D+00
!
!  In the outer loop, we pick a sub-box.
!
  j = 0

  do

    call tuple_next ( 1, sub_num, dim_num, j, sub_indx )

    if ( j == 0 ) then
      exit
    end if

    do dim = 1, dim_num

      a_sub(dim) = ( real ( sub_num + 1 - sub_indx(dim), kind = 8 ) * a(dim)   &
                   + real (         - 1 + sub_indx(dim), kind = 8 ) * b(dim) ) &
                   / real ( sub_num,                     kind = 8 )

      b_sub(dim) = ( real ( sub_num     - sub_indx(dim), kind = 8 ) * a(dim)   &
                   + real (               sub_indx(dim), kind = 8 ) * b(dim) ) &
                   / real ( sub_num,                     kind = 8 )

    end do
!
!  In the inner loop, we go through all the points in the
!  cross product of the quadrature rule.
!
    k = 0

    do

      call tuple_next ( 1, norder, dim_num, k, indx )

      if ( k == 0  ) then
        exit
      end if

      w = product ( weight(indx(1:dim_num)) )

      x(1:dim_num,1) = 0.5D+00 * ( &
          a_sub(1:dim_num) * ( 1.0D+00 - xtab(indx(1:dim_num)) ) &
        + b_sub(1:dim_num) * ( 1.0D+00 + xtab(indx(1:dim_num)) ) )

      call p00_f ( problem, dim_num, point_num, x, value )

      result = result + w * value(1)

    end do

  end do
!
!  Get the volume.
!
  call p00_volume ( problem, dim_num, volume )

  result = result * volume / real ( ( 2 * sub_num )**dim_num, kind = 8 )

  return
end
subroutine p00_box_mc ( problem, dim_num, point_num, result )

!*****************************************************************************80
!
!! P00_BOX_MC integrates over an multi-dimensional box using Monte Carlo.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points to use.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) problem
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ), allocatable, dimension ( : ) :: value
  real ( kind = 8 ) volume
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x

  allocate ( value(1:point_num) )
  allocate ( x(1:dim_num,1:point_num) )
!
!  Get the volume.
!
  call p00_volume ( problem, dim_num, volume )
!
!  Get the integration limits.
!
  call p00_lim ( problem, dim_num, a, b )
!
!  Get the random values.
!
  call random_number ( harvest = x(1:dim_num,1:point_num) )
!
!  Map them to the domain.
!
  do dim = 1, dim_num
    x(dim,1:point_num) = ( ( 1.0D+00 - x(dim,1:point_num) ) * a(dim) &
                         + (         + x(dim,1:point_num) ) * b(dim) )
  end do
!
!  Evaluate the function.
!
  call p00_f ( problem, dim_num, point_num, x, value )

  result = sum ( value(1:point_num) ) * volume / real ( point_num, kind = 8 )

  deallocate ( value )
  deallocate ( x )

  return
end
subroutine p00_default ( problem, dim_num )

!*****************************************************************************80
!
!! P00_DEFAULT sets a problem to a default state.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the desired problem.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) problem
  integer ( kind = 4 ) dim_num

  if ( problem == 1 ) then
    call p01_default ( dim_num )
  else if ( problem == 2 ) then
    call p02_default ( dim_num )
  else if ( problem == 3 ) then
    call p03_default ( dim_num )
  else if ( problem == 4 ) then
    call p04_default ( dim_num )
  else if ( problem == 5 ) then
    call p05_default ( dim_num )
  else if ( problem == 6 ) then
    call p06_default ( dim_num )
  else if ( problem == 7 ) then
    call p07_default ( dim_num )
  else if ( problem == 8 ) then
    call p08_default ( dim_num )
  else if ( problem == 9 ) then
    call p09_default ( dim_num )
  else if ( problem == 10 ) then
    call p10_default ( dim_num )
  else if ( problem == 11 ) then
    call p11_default ( dim_num )
  else if ( problem == 12 ) then
    call p12_default ( dim_num )
  else if ( problem == 13 ) then
    call p13_default ( dim_num )
  else if ( problem == 14 ) then
    call p14_default ( dim_num )
  else if ( problem == 15 ) then
    call p15_default ( dim_num )
  else if ( problem == 16 ) then
    call p16_default ( dim_num )
  else if ( problem == 17 ) then
    call p17_default ( dim_num )
  else if ( problem == 18 ) then
    call p18_default ( dim_num )
  else if ( problem == 19 ) then
    call p19_default ( dim_num )
  else if ( problem == 20 ) then
    call p20_default ( dim_num )
  else if ( problem == 21 ) then
    call p21_default ( dim_num )
  else if ( problem == 22 ) then
    call p22_default ( dim_num )
  else if ( problem == 23 ) then
    call p23_default ( dim_num )
  else if ( problem == 24 ) then
    call p24_default ( dim_num )
  else if ( problem == 25 ) then
    call p25_default ( dim_num )
  else if ( problem == 26 ) then
    call p26_default ( dim_num )
  else if ( problem == 27 ) then
    call p27_default ( dim_num )
  else if ( problem == 28 ) then
    call p28_default ( dim_num )
  else if ( problem == 29 ) then
    call p29_default ( dim_num )
  else if ( problem == 30 ) then
    call p30_default ( dim_num )
  else if ( problem == 31 ) then
    call p31_default ( dim_num )
  else if ( problem == 32 ) then
    call p32_default ( dim_num )
  else if ( problem == 33 ) then
    call p33_default ( dim_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_DEFAULT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_exact ( problem, dim_num, exact )

!*****************************************************************************80
!
!! P00_EXACT returns the exact integral for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the desired problem.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) dim_num

  if ( problem == 1 ) then
    call p01_exact ( dim_num, exact )
  else if ( problem == 2 ) then
    call p02_exact ( dim_num, exact )
  else if ( problem == 3 ) then
    call p03_exact ( dim_num, exact )
  else if ( problem == 4 ) then
    call p04_exact ( dim_num, exact )
  else if ( problem == 5 ) then
    call p05_exact ( dim_num, exact )
  else if ( problem == 6 ) then
    call p06_exact ( dim_num, exact )
  else if ( problem == 7 ) then
    call p07_exact ( dim_num, exact )
  else if ( problem == 8 ) then
    call p08_exact ( dim_num, exact )
  else if ( problem == 9 ) then
    call p09_exact ( dim_num, exact )
  else if ( problem == 10 ) then
    call p10_exact ( dim_num, exact )
  else if ( problem == 11 ) then
    call p11_exact ( dim_num, exact )
  else if ( problem == 12 ) then
    call p12_exact ( dim_num, exact )
  else if ( problem == 13 ) then
    call p13_exact ( dim_num, exact )
  else if ( problem == 14 ) then
    call p14_exact ( dim_num, exact )
  else if ( problem == 15 ) then
    call p15_exact ( dim_num, exact )
  else if ( problem == 16 ) then
    call p16_exact ( dim_num, exact )
  else if ( problem == 17 ) then
    call p17_exact ( dim_num, exact )
  else if ( problem == 18 ) then
    call p18_exact ( dim_num, exact )
  else if ( problem == 19 ) then
    call p19_exact ( dim_num, exact )
  else if ( problem == 20 ) then
    call p20_exact ( dim_num, exact )
  else if ( problem == 21 ) then
    call p21_exact ( dim_num, exact )
  else if ( problem == 22 ) then
    call p22_exact ( dim_num, exact )
  else if ( problem == 23 ) then
    call p23_exact ( dim_num, exact )
  else if ( problem == 24 ) then
    call p24_exact ( dim_num, exact )
  else if ( problem == 25 ) then
    call p25_exact ( dim_num, exact )
  else if ( problem == 26 ) then
    call p26_exact ( dim_num, exact )
  else if ( problem == 27 ) then
    call p27_exact ( dim_num, exact )
  else if ( problem == 28 ) then
    call p28_exact ( dim_num, exact )
  else if ( problem == 29 ) then
    call p29_exact ( dim_num, exact )
  else if ( problem == 30 ) then
    call p30_exact ( dim_num, exact )
  else if ( problem == 31 ) then
    call p31_exact ( dim_num, exact )
  else if ( problem == 32 ) then
    call p32_exact ( dim_num, exact )
  else if ( problem == 33 ) then
    call p33_exact ( dim_num, exact )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_EXACT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_f ( problem, dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P00_F evaluates the integrand for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the desired problem.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) problem
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  if ( problem == 1 ) then
    call p01_f ( dim_num, point_num, x, value )
  else if ( problem == 2 ) then
    call p02_f ( dim_num, point_num, x, value )
  else if ( problem == 3 ) then
    call p03_f ( dim_num, point_num, x, value )
  else if ( problem == 4 ) then
    call p04_f ( dim_num, point_num, x, value )
  else if ( problem == 5 ) then
    call p05_f ( dim_num, point_num, x, value )
  else if ( problem == 6 ) then
    call p06_f ( dim_num, point_num, x, value )
  else if ( problem == 7 ) then
    call p07_f ( dim_num, point_num, x, value )
  else if ( problem == 8 ) then
    call p08_f ( dim_num, point_num, x, value )
  else if ( problem == 9 ) then
    call p09_f ( dim_num, point_num, x, value )
  else if ( problem == 10 ) then
    call p10_f ( dim_num, point_num, x, value )
  else if ( problem == 11 ) then
    call p11_f ( dim_num, point_num, x, value )
  else if ( problem == 12 ) then
    call p12_f ( dim_num, point_num, x, value )
  else if ( problem == 13 ) then
    call p13_f ( dim_num, point_num, x, value )
  else if ( problem == 14 ) then
    call p14_f ( dim_num, point_num, x, value )
  else if ( problem == 15 ) then
    call p15_f ( dim_num, point_num, x, value )
  else if ( problem == 16 ) then
    call p16_f ( dim_num, point_num, x, value )
  else if ( problem == 17 ) then
    call p17_f ( dim_num, point_num, x, value )
  else if ( problem == 18 ) then
    call p18_f ( dim_num, point_num, x, value )
  else if ( problem == 19 ) then
    call p19_f ( dim_num, point_num, x, value )
  else if ( problem == 20 ) then
    call p20_f ( dim_num, point_num, x, value )
  else if ( problem == 21 ) then
    call p21_f ( dim_num, point_num, x, value )
  else if ( problem == 22 ) then
    call p22_f ( dim_num, point_num, x, value )
  else if ( problem == 23 ) then
    call p23_f ( dim_num, point_num, x, value )
  else if ( problem == 24 ) then
    call p24_f ( dim_num, point_num, x, value )
  else if ( problem == 25 ) then
    call p25_f ( dim_num, point_num, x, value )
  else if ( problem == 26 ) then
    call p26_f ( dim_num, point_num, x, value )
  else if ( problem == 27 ) then
    call p27_f ( dim_num, point_num, x, value )
  else if ( problem == 28 ) then
    call p28_f ( dim_num, point_num, x, value )
  else if ( problem == 29 ) then
    call p29_f ( dim_num, point_num, x, value )
  else if ( problem == 30 ) then
    call p30_f ( dim_num, point_num, x, value )
  else if ( problem == 31 ) then
    call p31_f ( dim_num, point_num, x, value )
  else if ( problem == 32 ) then
    call p32_f ( dim_num, point_num, x, value )
  else if ( problem == 33 ) then
    call p33_f ( dim_num, point_num, x, value )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_F - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_i4 ( problem, action, name, value )

!*****************************************************************************80
!
!! P00_I4 sets or gets I4 parameters for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the desired problem.
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ) problem
  character ( len = * ) name
  integer ( kind = 4 ) value

  if ( problem == 1 ) then
    call p01_i4 ( action, name, value )
  else if ( problem == 2 ) then
    call p02_i4 ( action, name, value )
  else if ( problem == 3 ) then
    call p03_i4 ( action, name, value )
  else if ( problem == 4 ) then
    call p04_i4 ( action, name, value )
  else if ( problem == 5 ) then
    call p05_i4 ( action, name, value )
  else if ( problem == 6 ) then
    call p06_i4 ( action, name, value )
  else if ( problem == 7 ) then
    call p07_i4 ( action, name, value )
  else if ( problem == 8 ) then
    call p08_i4 ( action, name, value )
  else if ( problem == 9 ) then
    call p09_i4 ( action, name, value )
  else if ( problem == 10 ) then
    call p10_i4 ( action, name, value )
  else if ( problem == 11 ) then
    call p11_i4 ( action, name, value )
  else if ( problem == 12 ) then
    call p12_i4 ( action, name, value )
  else if ( problem == 13 ) then
    call p13_i4 ( action, name, value )
  else if ( problem == 14 ) then
    call p14_i4 ( action, name, value )
  else if ( problem == 15 ) then
    call p15_i4 ( action, name, value )
  else if ( problem == 16 ) then
    call p16_i4 ( action, name, value )
  else if ( problem == 17 ) then
    call p17_i4 ( action, name, value )
  else if ( problem == 18 ) then
    call p18_i4 ( action, name, value )
  else if ( problem == 19 ) then
    call p19_i4 ( action, name, value )
  else if ( problem == 20 ) then
    call p20_i4 ( action, name, value )
  else if ( problem == 21 ) then
    call p21_i4 ( action, name, value )
  else if ( problem == 22 ) then
    call p22_i4 ( action, name, value )
  else if ( problem == 23 ) then
    call p23_i4 ( action, name, value )
  else if ( problem == 24 ) then
    call p24_i4 ( action, name, value )
  else if ( problem == 25 ) then
    call p25_i4 ( action, name, value )
  else if ( problem == 26 ) then
    call p26_i4 ( action, name, value )
  else if ( problem == 27 ) then
    call p27_i4 ( action, name, value )
  else if ( problem == 28 ) then
    call p28_i4 ( action, name, value )
  else if ( problem == 29 ) then
    call p29_i4 ( action, name, value )
  else if ( problem == 30 ) then
    call p30_i4 ( action, name, value )
  else if ( problem == 31 ) then
    call p31_i4 ( action, name, value )
  else if ( problem == 32 ) then
    call p32_i4 ( action, name, value )
  else if ( problem == 33 ) then
    call p33_i4 ( action, name, value )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_I4 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_lim ( problem, dim_num, a, b )

!*****************************************************************************80
!
!! P00_LIM returns the integration limits for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the test problem.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper 
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_lim ( dim_num, a, b )
  else if ( problem == 2 ) then
    call p02_lim ( dim_num, a, b )
  else if ( problem == 3 ) then
    call p03_lim ( dim_num, a, b )
  else if ( problem == 4 ) then
    call p04_lim ( dim_num, a, b )
  else if ( problem == 5 ) then
    call p05_lim ( dim_num, a, b )
  else if ( problem == 6 ) then
    call p06_lim ( dim_num, a, b )
  else if ( problem == 7 ) then
    call p07_lim ( dim_num, a, b )
  else if ( problem == 8 ) then
    call p08_lim ( dim_num, a, b )
  else if ( problem == 9 ) then
    call p09_lim ( dim_num, a, b )
  else if ( problem == 10 ) then
    call p10_lim ( dim_num, a, b )
  else if ( problem == 11 ) then
    call p11_lim ( dim_num, a, b )
  else if ( problem == 12 ) then
    call p12_lim ( dim_num, a, b )
  else if ( problem == 13 ) then
    call p13_lim ( dim_num, a, b )
  else if ( problem == 14 ) then
    call p14_lim ( dim_num, a, b )
  else if ( problem == 15 ) then
    call p15_lim ( dim_num, a, b )
  else if ( problem == 16 ) then
    call p16_lim ( dim_num, a, b )
  else if ( problem == 17 ) then
    call p17_lim ( dim_num, a, b )
  else if ( problem == 18 ) then
    call p18_lim ( dim_num, a, b )
  else if ( problem == 19 ) then
    call p19_lim ( dim_num, a, b )
  else if ( problem == 20 ) then
    call p20_lim ( dim_num, a, b )
  else if ( problem == 21 ) then
    call p21_lim ( dim_num, a, b )
  else if ( problem == 22 ) then
    call p22_lim ( dim_num, a, b )
  else if ( problem == 23 ) then
    call p23_lim ( dim_num, a, b )
  else if ( problem == 24 ) then
    call p24_lim ( dim_num, a, b )
  else if ( problem == 25 ) then
    call p25_lim ( dim_num, a, b )
  else if ( problem == 26 ) then
    call p26_lim ( dim_num, a, b )
  else if ( problem == 27 ) then
    call p27_lim ( dim_num, a, b )
  else if ( problem == 28 ) then
    call p28_lim ( dim_num, a, b )
  else if ( problem == 29 ) then
    call p29_lim ( dim_num, a, b )
  else if ( problem == 30 ) then
    call p30_lim ( dim_num, a, b )
  else if ( problem == 31 ) then
    call p31_lim ( dim_num, a, b )
  else if ( problem == 32 ) then
    call p32_lim ( dim_num, a, b )
  else if ( problem == 33 ) then
    call p33_lim ( dim_num, a, b )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_LIM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_name ( problem, name )

!*****************************************************************************80
!
!! P00_NAME returns the name of the problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the index of the test problem.
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  integer ( kind = 4 ) problem
  character ( len = * ) name

  if ( problem == 1 ) then
    call p01_name ( name )
  else if ( problem == 2 ) then
    call p02_name ( name )
  else if ( problem == 3 ) then
    call p03_name ( name )
  else if ( problem == 4 ) then
    call p04_name ( name )
  else if ( problem == 5 ) then
    call p05_name ( name )
  else if ( problem == 6 ) then
    call p06_name ( name )
  else if ( problem == 7 ) then
    call p07_name ( name )
  else if ( problem == 8 ) then
    call p08_name ( name )
  else if ( problem == 9 ) then
    call p09_name ( name )
  else if ( problem == 10 ) then
    call p10_name ( name )
  else if ( problem == 11 ) then
    call p11_name ( name )
  else if ( problem == 12 ) then
    call p12_name ( name )
  else if ( problem == 13 ) then
    call p13_name ( name )
  else if ( problem == 14 ) then
    call p14_name ( name )
  else if ( problem == 15 ) then
    call p15_name ( name )
  else if ( problem == 16 ) then
    call p16_name ( name )
  else if ( problem == 17 ) then
    call p17_name ( name )
  else if ( problem == 18 ) then
    call p18_name ( name )
  else if ( problem == 19 ) then
    call p19_name ( name )
  else if ( problem == 20 ) then
    call p20_name ( name )
  else if ( problem == 21 ) then
    call p21_name ( name )
  else if ( problem == 22 ) then
    call p22_name ( name )
  else if ( problem == 23 ) then
    call p23_name ( name )
  else if ( problem == 24 ) then
    call p24_name ( name )
  else if ( problem == 25 ) then
    call p25_name ( name )
  else if ( problem == 26 ) then
    call p26_name ( name )
  else if ( problem == 27 ) then
    call p27_name ( name )
  else if ( problem == 28 ) then
    call p28_name ( name )
  else if ( problem == 29 ) then
    call p29_name ( name )
  else if ( problem == 30 ) then
    call p30_name ( name )
  else if ( problem == 31 ) then
    call p31_name ( name )
  else if ( problem == 32 ) then
    call p32_name ( name )
  else if ( problem == 33 ) then
    call p33_name ( name )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_NAME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_r8vec ( problem, action, name, dim_num, value )

!*****************************************************************************80
!
!! P00_R8VEC sets or gets R8VEC parameters for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the number of the test problem.
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the  object to a default value.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the object should be set to VALUE.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the variable.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the variable.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  integer ( kind = 4 ) problem
  character ( len = * ) name
  real ( kind = 8 ) value(dim_num)

  if ( problem == 1 ) then

  else if ( problem == 2 ) then

  else if ( problem == 3 ) then

  else if ( problem == 4 ) then

  else if ( problem == 5 ) then

  else if ( problem == 6 ) then

  else if ( problem == 7 ) then

  else if ( problem == 8 ) then

  else if ( problem == 9 ) then
    call p09_r8vec ( action, name, dim_num, value )
  else if ( problem == 10 ) then

  else if ( problem == 11 ) then

  else if ( problem == 12 ) then

  else if ( problem == 13 ) then

  else if ( problem == 14 ) then

  else if ( problem == 15 ) then

  else if ( problem == 16 ) then
    call p16_r8vec ( action, name, dim_num, value )
  else if ( problem == 17 ) then
    call p17_r8vec ( action, name, dim_num, value )
  else if ( problem == 18 ) then
    call p18_r8vec ( action, name, dim_num, value )
  else if ( problem == 19 ) then
    call p19_r8vec ( action, name, dim_num, value )
  else if ( problem == 20 ) then

  else if ( problem == 21 ) then

  else if ( problem == 22 ) then

  else if ( problem == 23 ) then

  else if ( problem == 24 ) then
    call p24_r8vec ( action, name, dim_num, value )
  else if ( problem == 25 ) then

  else if ( problem == 26 ) then
    call p26_r8vec ( action, name, dim_num, value )
  else if ( problem == 27 ) then
    call p27_r8vec ( action, name, dim_num, value )
  else if ( problem == 28 ) then
    call p28_r8vec ( action, name, dim_num, value )
  else if ( problem == 29 ) then
    call p29_r8vec ( action, name, dim_num, value )
  else if ( problem == 30 ) then
    call p30_r8vec ( action, name, dim_num, value )
  else if ( problem == 31 ) then
    call p31_r8vec ( action, name, dim_num, value )
  else if ( problem == 32 ) then
    call p32_r8vec ( action, name, dim_num, value )
  else if ( problem == 33 ) then

  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_R8VEC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_region ( problem, region )

!*****************************************************************************80
!
!! P00_REGION returns the name of the integration region for any problem.
!
!  Discussion:
!
!    I thought I was going to use this idea a lot, but most of my test
!    regions are boxes.
!
!    BALL
!      the interior of a 2D circle,
!      the interior of a 3D sphere,
!      the interior of an ND sphere.
!
!    BOX
!      a 1D finite line segment,
!      a 2D finite rectangle,
!      a 3D box,
!      an ND box.
!
!    SIMPLEX 
!      a 2D triangle,
!      a 3D tetrahedron,
!      an ND simplex.
!      The "unit simplex" in ND is the set of nonnegative points X 
!      such that sum ( X(1:N) ) <= 1.
!
!    SPACE
!      a 1D infinite line,
!      a 2D infinite place,
!      a 3D space,
!      an ND space.
!
!    SPHERE
!      the circumference of a 2D circle,
!      the surface of a 3D sphere,
!      the surface of an ND sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the number of the test problem.
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  integer ( kind = 4 ) problem
  character ( len = * ) region

  if ( problem == 1 ) then
    call p01_region ( region )
  else if ( problem == 2 ) then
    call p02_region ( region )
  else if ( problem == 3 ) then
    call p03_region ( region )
  else if ( problem == 4 ) then
    call p04_region ( region )
  else if ( problem == 5 ) then
    call p05_region ( region )
  else if ( problem == 6 ) then
    call p06_region ( region )
  else if ( problem == 7 ) then
    call p07_region ( region )
  else if ( problem == 8 ) then
    call p08_region ( region )
  else if ( problem == 9 ) then
    call p09_region ( region )
  else if ( problem == 10 ) then
    call p10_region ( region )
  else if ( problem == 11 ) then
    call p11_region ( region )
  else if ( problem == 12 ) then
    call p12_region ( region )
  else if ( problem == 13 ) then
    call p13_region ( region )
  else if ( problem == 14 ) then
    call p14_region ( region )
  else if ( problem == 15 ) then
    call p15_region ( region )
  else if ( problem == 16 ) then
    call p16_region ( region )
  else if ( problem == 17 ) then
    call p17_region ( region )
  else if ( problem == 18 ) then
    call p18_region ( region )
  else if ( problem == 19 ) then
    call p19_region ( region )
  else if ( problem == 20 ) then
    call p20_region ( region )
  else if ( problem == 21 ) then
    call p21_region ( region )
  else if ( problem == 22 ) then
    call p22_region ( region )
  else if ( problem == 23 ) then
    call p23_region ( region )
  else if ( problem == 24 ) then
    call p24_region ( region )
  else if ( problem == 25 ) then
    call p25_region ( region )
  else if ( problem == 26 ) then
    call p26_region ( region )
  else if ( problem == 27 ) then
    call p27_region ( region )
  else if ( problem == 28 ) then
    call p28_region ( region )
  else if ( problem == 29 ) then
    call p29_region ( region )
  else if ( problem == 30 ) then
    call p30_region ( region )
  else if ( problem == 31 ) then
    call p31_region ( region )
  else if ( problem == 32 ) then
    call p32_region ( region )
  else if ( problem == 33 ) then
    call p33_region ( region )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_REGION - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_remap01 ( problem, dim_num, point_num, x01, xab )

!*****************************************************************************80
!
!! P00_REMAP01 remaps points in [0,1]^DIM_NUM into [A(1:DIM_NUM),B(1:DIM_NUM)].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X01(DIM_NUM,POINT_NUM), the points, in [0,1], 
!    to be transformed.
!
!    Output, real ( kind = 8 ) XAB(DIM_NUM,POINT_NUM), the transformed points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x01(dim_num,point_num)
  real ( kind = 8 ) xab(dim_num,point_num)

  call p00_lim ( problem, dim_num, a, b )

  do dim = 1, dim_num
    xab(dim,1:point_num) =  a(dim) + ( b(dim) - a(dim) ) * x01(dim,1:point_num)
  end do

  return
end
subroutine p00_title ( problem )

!*****************************************************************************80
!
!! P00_TITLE prints a title for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the number of the test problem.
!
  implicit none

  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_title ( )
  else if ( problem == 2 ) then
    call p02_title ( )
  else if ( problem == 3 ) then
    call p03_title ( )
  else if ( problem == 4 ) then
    call p04_title ( )
  else if ( problem == 5 ) then
    call p05_title ( )
  else if ( problem == 6 ) then
    call p06_title ( )
  else if ( problem == 7 ) then
    call p07_title ( )
  else if ( problem == 8 ) then
    call p08_title ( )
  else if ( problem == 9 ) then
    call p09_title ( )
  else if ( problem == 10 ) then
    call p10_title ( )
  else if ( problem == 11 ) then
    call p11_title ( )
  else if ( problem == 12 ) then
    call p12_title ( )
  else if ( problem == 13 ) then
    call p13_title ( )
  else if ( problem == 14 ) then
    call p14_title ( )
  else if ( problem == 15 ) then
    call p15_title ( )
  else if ( problem == 16 ) then
    call p16_title ( )
  else if ( problem == 17 ) then
    call p17_title ( )
  else if ( problem == 18 ) then
    call p18_title ( )
  else if ( problem == 19 ) then
    call p19_title ( )
  else if ( problem == 20 ) then
    call p20_title ( )
  else if ( problem == 21 ) then
    call p21_title ( )
  else if ( problem == 22 ) then
    call p22_title ( )
  else if ( problem == 23 ) then
    call p23_title ( )
  else if ( problem == 24 ) then
    call p24_title ( )
  else if ( problem == 25 ) then
    call p25_title ( )
  else if ( problem == 26 ) then
    call p26_title ( )
  else if ( problem == 27 ) then
    call p27_title ( )
  else if ( problem == 28 ) then
    call p28_title ( )
  else if ( problem == 29 ) then
    call p29_title ( )
  else if ( problem == 30 ) then
    call p30_title ( )
  else if ( problem == 31 ) then
    call p31_title ( )
  else if ( problem == 32 ) then
    call p32_title ( )
  else if ( problem == 33 ) then
    call p33_title ( )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_volume ( problem, dim_num, volume )

!*****************************************************************************80
!
!! P00_VOLUME returns the volume of the integration region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the number of the test problem.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the integration region.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) problem
  real ( kind = 8 ) r
  real ( kind = 8 ) volume

  if ( 1 <= problem .and. problem <= 20 ) then

    call p00_lim ( problem, dim_num, a, b )

    volume = product ( b(1:dim_num) - a(1:dim_num) )

  else if ( problem == 21 ) then

    call sphere_unit_area_nd ( dim_num, volume )

  else if ( problem == 22 ) then

    call p22_r8 ( 'G', 'R', r )

    call sphere_volume_nd ( dim_num, r, volume )

  else if ( problem == 23 ) then

    call simplex_unit_volume_nd ( dim_num, volume )

  else if ( 24 <= problem .and. problem <= 32 ) then

    call p00_lim ( problem, dim_num, a, b )

    volume = product ( b(1:dim_num) - a(1:dim_num) )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_VOLUME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p01_default ( dim_num )

!*****************************************************************************80
!
!! P01_DEFAULT sets default values for problem 01.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) i4

  call p01_i4 ( 'D', '*', i4 )

  return
end
subroutine p01_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P01_EXACT returns the exact integral for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num

  exact = real ( ( dim_num ) * ( 3 * dim_num + 1 ), kind = 8 ) / 12.0D+00

  return
end
subroutine p01_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P01_F evaluates the integrand for problem 01.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    f(x) = ( sum ( x(1:dim_num) ) )^2
!
!  Exact Integral:
!
!    dim_num * ( 3 * dim_num + 1 ) / 12
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = ( sum ( x(1:dim_num,point) ) )**2
  end do

  call p01_i4 ( 'I', '#', point_num )

  return
end
subroutine p01_i4 ( action, name, value )

!*****************************************************************************80
!
!! P01_I4 sets or gets I4 parameters for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P01_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P01_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P01_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P01_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p01_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P01_LIM returns the integration limits for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  return
end
subroutine p01_name ( name )

!*****************************************************************************80
!
!! P01_NAME returns the name of problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'SquareSum'

  return
end
subroutine p01_region ( region )

!*****************************************************************************80
!
!! P01_REGION returns the name of the integration region for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p01_title ( )

!*****************************************************************************80
!
!! P01_TITLE prints a title for problem 01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 01'
  write ( *, '(a)' ) '  Name:       SquareSum'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = ( sum ( X(i) ) )^2'

  return
end
subroutine p02_default ( dim_num )

!*****************************************************************************80
!
!! P02_DEFAULT sets default values for problem 02.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p02_i4 ( 'D', '*', i4 )

  return
end
subroutine p02_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P02_EXACT returns the exact integral for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) exact

  exact = real ( dim_num * ( 5 * dim_num - 2 ), kind = 8 ) / 15.0D+00

  return
end
subroutine p02_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P02_F evaluates the integrand for problem 02.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    ( sum ( 2 * x(1:dim_num) - 1 ) )^4
!
!  Exact Integral:
!
!    DIM_NUM * (5*DIM_NUM-2) / 15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = ( sum ( 2.0D+00 * x(1:dim_num,point) - 1.0D+00 ) )**4
  end do

  call p02_i4 ( 'I', '#', point_num )

  return
end
subroutine p02_i4 ( action, name, value )

!*****************************************************************************80
!
!! P02_I4 sets or gets I4 parameters for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P02_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P02_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P02_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P02_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p02_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P02_LIM returns the integration limits for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the limits
!    of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p02_name ( name )

!*****************************************************************************80
!
!! P02_NAME returns the name of problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'QuadSum'

  return
end
subroutine p02_region ( region )

!*****************************************************************************80
!
!! P02_REGION returns the name of the integration region for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p02_title ( )

!*****************************************************************************80
!
!! P02_TITLE prints a title for problem 02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 02'
  write ( *, '(a)' ) '  Name:       QuadSum'
  write ( *, '(a)' ) '              Davis, Rabinowitz, page 370, #1.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = ( sum ( 2 * X(i) - 1 ) )^4'

  return
end
subroutine p03_default ( dim_num )

!*****************************************************************************80
!
!! P03_DEFAULT sets default values for problem 03.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p03_i4 ( 'D', '*', i4 )

  return
end
subroutine p03_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P03_EXACT returns the exact integral for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num

  exact = 0.0D+00

  return
end
subroutine p03_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P03_F evaluates the integrand for problem 03.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    ( sum ( 2 * x(1:dim_num) - 1 ) )^5
!
!  Exact Integral:
!
!    0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = ( sum ( 2.0D+00 * x(1:dim_num,point) - 1.0D+00 ) )**5
  end do

  call p03_i4 ( 'I', '#', point_num )

  return
end
subroutine p03_i4 ( action, name, value )

!*****************************************************************************80
!
!! P03_I4 sets or gets I4 parameters for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P03_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P03_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P03_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P03_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p03_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P03_LIM returns the integration limits for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p03_name ( name )

!*****************************************************************************80
!
!! P03_NAME returns the name of problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'QuintSum'

  return
end
subroutine p03_region ( region )

!*****************************************************************************80
!
!! P03_REGION returns the name of the integration region for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p03_title ( )

!*****************************************************************************80
!
!! P03_TITLE prints a title for problem 03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 03'
  write ( *, '(a)' ) '  Name:       QuintSum'
  write ( *, '(a)' ) '              Davis, Rabinowitz, page 370, #3.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = ( sum ( X(i) ) )^5'

  return
end
subroutine p04_default ( dim_num )

!*****************************************************************************80
!
!! P04_DEFAULT sets default values for problem 04.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p04_i4 ( 'D', '*', i4 )

  return
end
subroutine p04_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P04_EXACT returns the exact integral for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) exact

  exact = real ( dim_num * ( 7 * ( dim_num - 1 ) * ( 5 * dim_num - 1 ) + 9 ), &
    kind = 8 ) / 63.0D+00

  return
end
subroutine p04_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P04_F evaluates the integrand for problem 04.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    ( sum ( 2 * x(1:dim_num) - 1 ) )^6
!
!  Exact Integral:
!
!    DIM_NUM * ( 7 * (DIM_NUM-1) * (5*DIM_NUM-1) + 9 ) / 63
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = ( sum ( 2.0D+00 * x(1:dim_num,point) - 1.0D+00 ) )**6
  end do

  call p04_i4 ( 'I', '#', point_num )

  return
end
subroutine p04_i4 ( action, name, value )

!*****************************************************************************80
!
!! P04_I4 sets or gets I4 parameters for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P04_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P04_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P04_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P04_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p04_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P04_LIM returns the integration limits for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM),B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  return
end
subroutine p04_name ( name )

!*****************************************************************************80
!
!! P04_NAME returns the name of problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'HexSum'

  return
end
subroutine p04_region ( region )

!*****************************************************************************80
!
!! P04_REGION returns the name of the integration region for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p04_title ( )

!*****************************************************************************80
!
!! P04_TITLE prints a title for problem 04.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 04'
  write ( *, '(a)' ) '  Name:       HexSum'
  write ( *, '(a)' ) '              Davis, Rabinowitz, page 370, #2.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = ( sum ( 2 * X(i) - 1 ) )^6'

  return
end
subroutine p05_default ( dim_num )

!*****************************************************************************80
!
!! P05_DEFAULT sets default values for problem 05.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p05_i4 ( 'D', '*', i4 )

  return
end
subroutine p05_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P05_EXACT returns the exact integral for problem 05.
!
!  Discussion:
!
!    The exact value is given only for DIM_NUM = 1, 2, 3, 4 or 5.
!    For other cases, the value HUGE is returned instead.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral,
!    or HUGE if the exact value is not known.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num

  if ( dim_num == 1 ) then

    exact = log ( 3.0D+00 )

  else if ( dim_num == 2 ) then

    exact = 5.0D+00 * log ( 5.0D+00 ) - 6.0D+00 * log ( 3.0D+00 )

  else if ( dim_num == 3 ) then

    exact = 0.5D+00 * ( 49.0D+00 * log ( 7.0D+00 ) &
      - 75.0D+00 * log ( 5.0D+00 ) + 27.0D+00 * log ( 3.0D+00 ) )

  else if ( dim_num == 4 ) then

    exact = 225.0D+00 * log ( 3.0D+00 ) + 125.0D+00 * log ( 5.0D+00 ) &
      - 686.0D+00 * log ( 7.0D+00 ) / 3.0D+00

  else if ( dim_num == 5 ) then

    exact = ( &
      - 65205.0D+00 * log ( 3.0D+00 ) &
      - 6250.0D+00 * log ( 5.0D+00 ) &
      + 24010.0D+00 * log ( 7.0D+00 ) &
      + 14641.0D+00 * log ( 11.0D+00 ) ) / 24.0D+00

  else

    exact = huge ( exact )

  end if

  return
end
subroutine p05_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P05_F evaluates the integrand for problem 05.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    2^DIM_NUM / ( 1 + sum ( 2 * x(1:dim_num) ) )
!
!  Exact Integral:
!
!    For DIM_NUM = 1:
!
!      ln ( 3 )
!
!    For DIM_NUM = 2:
!
!      ln ( 3125 / 729 )
!
!    For DIM_NUM = 3:
!
!      0.5 * ( 49 * ln ( 7 ) - 75 * ln ( 5 ) + 27 * ln ( 3 ) )
! 
!    For DIM_NUM = 4:
!
!      225 * ln ( 3 ) + 125 * ln ( 5 ) - 686 * ln ( 7 ) / 3
!
!    For DIM_NUM = 5:
!
!      ( -65205 * ln ( 3 ) - 6250 * ln ( 5 ) + 24010 * ln ( 7 ) 
!      + 14641 * ln ( 11 ) ) / 24
!
!  Approximate Integral:
!
!    DIM_NUM  VALUE
!
!       1  1.098612289
!       2  1.455514830
!       3  2.152142833
!       4  3.402716587
!       5  5.620255523
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = 2.0D+00**dim_num &
      / ( 1.0D+00 + sum ( 2.0D+00 * x(1:dim_num,point) ) )
  end do

  call p05_i4 ( 'I', '#', point_num )

  return
end
subroutine p05_i4 ( action, name, value )

!*****************************************************************************80
!
!! P05_I4 sets or gets I4 parameters for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P05_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P05_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P05_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P05_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p05_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P05_LIM returns the integration limits for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p05_name ( name )

!*****************************************************************************80
!
!! P05_NAME returns the name of problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'ST04'

  return
end
subroutine p05_region ( region )

!*****************************************************************************80
!
!! P05_REGION returns the name of the integration region for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p05_title ( )

!*****************************************************************************80
!
!! P05_TITLE prints a title for problem 05.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 05'
  write ( *, '(a)' ) '  Name:       ST04'
  write ( *, '(a)' ) '              Stroud #4, page 26.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = 1 / ( 1 + sum ( 2 * X(i) ) )'

  return
end
subroutine p06_default ( dim_num )

!*****************************************************************************80
!
!! P06_DEFAULT sets default values for problem 06.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p06_i4 ( 'D', '*', i4 )

  return
end
subroutine p06_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P06_EXACT returns the exact integral for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num

  exact = 1.0D+00

  return
end
subroutine p06_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P06_F evaluates the integrand for problem 06.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    product ( 2 * abs ( 2 * x(1:dim_num) - 1 ) )
!
!  Exact Integral:
!
!    1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Implementation and Tests of Low-Discrepancy Sequences,
!    ACM Transactions on Modeling and Computer Simulation,
!    Volume 2, Number 3, pages 195-213, 1992.
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = product ( 2.0D+00 &
      * abs ( 2.0D+00 * x(1:dim_num,point) - 1.0D+00 ) )
  end do

  call p06_i4 ( 'I', '#', point_num )

  return
end
subroutine p06_i4 ( action, name, value )

!*****************************************************************************80
!
!! P06_I4 sets or gets I4 parameters for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P06_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P06_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P06_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p06_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P06_LIM returns the integration limits for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p06_name ( name )

!*****************************************************************************80
!
!! P06_NAME returns the name of problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'DR4061'

  return
end
subroutine p06_region ( region )

!*****************************************************************************80
!
!! P06_REGION returns the name of the integration region for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p06_title ( )

!*****************************************************************************80
!
!! P06_TITLE prints a title for problem 06.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 07'
  write ( *, '(a)' ) '  Name:       DR4061'
  write ( *, '(a)' ) '              Davis, Rabinowitz, page 406, #1.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = product ( abs ( 4 * X(i) - 2 ) )'

  return
end
subroutine p07_default ( dim_num )

!*****************************************************************************80
!
!! P07_DEFAULT sets default values for problem 07.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p07_i4 ( 'D', '*', i4 )

  return
end
subroutine p07_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P07_EXACT returns the exact integral for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num

  exact = 1.0D+00

  return
end
subroutine p07_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P07_F evaluates the integrand for problem 07.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    product ( pi / 2 ) * sin ( pi * x(1:dim_num) )
!
!  Exact Integral:
!
!    1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = product ( 0.5D+00 * pi * sin ( pi * x(1:dim_num,point) ) )
  end do

  call p07_i4 ( 'I', '#', point_num )

  return
end
subroutine p07_i4 ( action, name, value )

!*****************************************************************************80
!
!! P07_I4 sets or gets I4 parameters for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P07_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P07_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P07_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P07_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p07_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P07_LIM returns the integration limits for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper 
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p07_name ( name )

!*****************************************************************************80
!
!! P07_NAME returns the name of problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'DR4062'

  return
end
subroutine p07_region ( region )

!*****************************************************************************80
!
!! P07_REGION returns the name of the integration region for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p07_title ( )

!*****************************************************************************80
!
!! P07_TITLE prints a title for problem 07.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 07'
  write ( *, '(a)' ) '  Name:       DR4062'
  write ( *, '(a)' ) '              Davis, Rabinowitz, page 406, #2.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = prod ( pi * sin ( pi * X(i) ) / 2 )'

  return
end
subroutine p08_default ( dim_num )

!*****************************************************************************80
!
!! P08_DEFAULT sets default values for problem 08.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p08_i4 ( 'D', '*', i4 )

  return
end
subroutine p08_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P08_EXACT returns the exact integral for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = 0.5D+00 - sqrt ( 2.0D+00**( 3 * dim_num ) ) &
    * cos ( 0.25D+00 * pi * real ( dim_num, kind = 8 ) ) &
    / ( 2.0D+00 * pi**dim_num )

  return
end
subroutine p08_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P08_F evaluates the integrand for problem 08.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    ( sin ( (pi/4) * sum ( x(1:dim_num) ) ) )^2
!
!  Exact Integral:
!
!    1/2 - sqrt ( 2^(3*DIM_NUM) ) * cos ( DIM_NUM * pi / 4 ) ) 
!      / ( 2 * pi^DIM_NUM )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Crandall,
!    Projects in Scientific Computing,
!    Springer, 2000,
!    ISBN: 0387950095,
!    LC: Q183.9.C733.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = ( sin ( pi * sum ( x(1:dim_num,point) ) / 4.0D+00 ) )**2
  end do

  call p08_i4 ( 'I', '#', point_num )

  return
end
subroutine p08_i4 ( action, name, value )

!*****************************************************************************80
!
!! P08_I4 sets or gets I4 parameters for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P08_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P08_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P08_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P08_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p08_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P08_LIM returns the integration limits for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper 
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p08_name ( name )

!*****************************************************************************80
!
!! P08_NAME returns the name of problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'RC01'

  return
end
subroutine p08_region ( region )

!*****************************************************************************80
!
!! P08_REGION returns the name of the integration region for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p08_title ( )

!*****************************************************************************80
!
!! P08_TITLE prints a title for problem 08.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 08'
  write ( *, '(a)' ) '  Name:       RC01'
  write ( *, '(a)' ) '              Crandall, page 49, #1'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = sin^2 ( pi/4 * sum ( X(i) ) )'

  return
end
subroutine p09_default ( dim_num )

!*****************************************************************************80
!
!! P09_DEFAULT sets default values for problem 09.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p09_i4 ( 'D', '*', i4 )
  call p09_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p09_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P09_EXACT returns the exact integral for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) c(dim_num)
  real ( kind = 8 ) exact

  call p09_r8vec ( 'G', 'C', dim_num, c )

  exact = product ( ( exp ( c(1:dim_num) ) - 1.0D+00 ) / c(1:dim_num) )

  return
end
subroutine p09_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P09_F evaluates the integrand for problem 09.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    The integral depends on a parameter vector C(1:N).
!
!    The reference suggests choosing C at random in [0,1]
!    and then multiplying by the normalizing factor (60/N).
!
!    To get or set C, call P09_R8VEC.
!
!    The default value of C(1:N) is 1/N.
!
!  Integrand:
!
!    exp ( sum ( c(1:dim_num) * x(1:dim_num) ) )
!
!  Exact Integral:
!
!    product ( ( exp ( c(1:n) - 1 ) / c(1:n) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Patterson,
!    [Integral #7]
!    On the Construction of a Practical Ermakov-Zolotukhin 
!    Multiple Integrator,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D. Reidel, 1987, pages 269-290,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p09_r8vec ( 'G', 'C', dim_num, c )

  do point = 1, point_num
    value(point) = exp ( dot_product ( c(1:dim_num), x(1:dim_num,point) ) )
  end do

  call p09_i4 ( 'I', '#', point_num )

  return
end
subroutine p09_i4 ( action, name, value )

!*****************************************************************************80
!
!! P09_I4 sets or gets I4 parameters for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P09_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P09_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P09_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p09_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P09_LIM returns the integration limits for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper 
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p09_name ( name )

!*****************************************************************************80
!
!! P09_NAME returns the name of problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Patterson #7'

  return
end
subroutine p09_region ( region )

!*****************************************************************************80
!
!! P09_REGION returns the name of the integration region for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p09_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P09_R8VEC sets or gets R8VEC parameters for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the coefficient vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: c
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( c ) ) then
      deallocate ( c )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( c(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c(1:dim_num) = 1.0D+00 / real ( dim_num, kind = 8 )
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P09_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c(1:dim_num) )
      c(1:dim_num) = c(1:dim_num) * 60.0D+00 / real ( dim_num, kind = 8 )
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P09_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P09_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p09_title ( )

!*****************************************************************************80
!
!! P09_TITLE prints a title for problem 09.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 09'
  write ( *, '(a)' ) '  Name:       Patterson #7'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = exp ( sum ( C(i) * X(i) ) )'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              C(1:DIM_NUM) defaults to 1/DIM_NUM.'

  return
end
subroutine p10_default ( dim_num )

!*****************************************************************************80
!
!! P10_DEFAULT sets default values for problem 10.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p10_i4 ( 'D', '*', i4 )

  return
end
subroutine p10_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P10_EXACT returns the exact integral for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num

  exact = real ( dim_num, kind = 8 ) / 4.0D+00

  return
end
subroutine p10_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P10_F evaluates the integrand for problem 10.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    sum ( abs ( x(1:dim_num) - 0.5 ) )
!
!  Exact Integral:
!
!    DIM_NUM / 4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Thomas Patterson,
!    [Integral #4],
!    On the Construction of a Practical Ermakov-Zolotukhin 
!    Multiple Integrator,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast, Graeme Fairweather,
!    D. Reidel, 1987, pages 269-290,
!    LC: QA299.3.N38.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = sum ( abs ( x(1:dim_num,point) - 0.5D+00 ) )
  end do

  call p10_i4 ( 'I', '#', point_num )

  return
end
subroutine p10_i4 ( action, name, value )

!*****************************************************************************80
!
!! P10_I4 sets or gets I4 parameters for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P10_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P10_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P10_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P10_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p10_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P10_LIM returns the integration limits for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p10_name ( name )

!*****************************************************************************80
!
!! P10_NAME returns the name of problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Patterson #4'

  return
end
subroutine p10_region ( region )

!*****************************************************************************80
!
!! P10_REGION returns the name of the integration region for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p10_title ( )

!*****************************************************************************80
!
!! P10_TITLE prints a title for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 10'
  write ( *, '(a)' ) '  Name:       Patterson #4'
  write ( *, '(a)' ) '              Stroud, page ?'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = sum ( abs ( X(i) - 0.5 ) )'

  return
end
subroutine p11_default ( dim_num )

!*****************************************************************************80
!
!! P11_DEFAULT sets default values for problem 11.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p11_i4 ( 'D', '*', i4 )

  return
end
subroutine p11_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P11_EXACT returns the exact integral for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num

  exact = ( exp ( 1.0D+00 ) - 1.0D+00 )**dim_num

  return
end
subroutine p11_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P11_F evaluates the integrand for problem 11.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    exp ( sum ( abs ( 2 * X(1:DIM_NUM) - 1 ) ) )
!
!  Exact Integral:
!
!    ( E - 1.0 )^DIM_NUM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Patterson,
!    [Integral #2],
!    On the Construction of a Practical Ermakov-Zolotukhin 
!    Multiple Integrator,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D. Reidel, 1987, pages 269-290,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = exp ( &
      sum ( abs ( 2.0D+00 * x(1:dim_num,point) - 1.0D+00 ) ) )
  end do

  call p11_i4 ( 'I', '#', point_num )

  return
end
subroutine p11_i4 ( action, name, value )

!*****************************************************************************80
!
!! P11_I4 sets or gets I4 parameters for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P11_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P11_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P11_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P11_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p11_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P11_LIM returns the integration limits for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  return
end
subroutine p11_name ( name )

!*****************************************************************************80
!
!! P11_NAME returns the name of problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Patterson #2, exp(sum(abs(X)))'

  return
end
subroutine p11_region ( region )

!*****************************************************************************80
!
!! P11_REGION returns the name of the integration region for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p11_title ( )

!*****************************************************************************80
!
!! P11_TITLE prints a title for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 11'
  write ( *, '(a)' ) '  Name:       Patterson #2, exp(sum(abs(X)))'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = exp ( sum ( abs ( X(i) )))'

  return
end
subroutine p12_default ( dim_num )

!*****************************************************************************80
!
!! P12_DEFAULT sets default values for problem 12.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p12_i4 ( 'D', '*', i4 )

  return
end
subroutine p12_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P12_EXACT returns the exact integral for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) exact

  exact = 1.0D+00
  do dim = 1, dim_num
    exact = exact * sin ( real ( dim, kind = 8 ) )
  end do

  return
end
subroutine p12_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P12_F evaluates the integrand for problem 12.
!
!  Discussion:
!
!    The highly oscillatory nature of the integrand makes this
!    a difficult and perhaps even dubious test.
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    product ( 1 <= i <= dim_num ) ( i * cos ( i * x(i) ) )
!
!  Exact Integral:
!
!    product ( 1 <= I <= DIM_NUM ) sin ( i )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Implementation and Tests of Low-Discrepancy Sequences,
!    ACM Transactions on Modeling and Computer Simulation,
!    Volume 2, Number 3, July 1992, pages 195-213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  value(1:point_num) = 1.0D+00

  do point = 1, point_num
    do dim = 1, dim_num
      value(point) = value(point) * real ( dim, kind = 8 ) &
        * cos ( real ( dim, kind = 8 ) * x(dim,point) )
    end do
  end do

  call p12_i4 ( 'I', '#', point_num )

  return
end
subroutine p12_i4 ( action, name, value )

!*****************************************************************************80
!
!! P12_I4 sets or gets I4 parameters for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P12_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P12_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P12_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p12_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P12_LIM returns the integration limits for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p12_name ( name )

!*****************************************************************************80
!
!! P12_NAME returns the name of problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'BFN02'

  return
end
subroutine p12_region ( region )

!*****************************************************************************80
!
!! P12_REGION returns the name of the integration region for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p12_title ( )

!*****************************************************************************80
!
!! P12_TITLE prints a title for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 12'
  write ( *, '(a)' ) '  Name:       BFN02'
  write ( *, '(a)' ) '              Bratley, Fox, Niederreiter, #2'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = product ( i * cos ( X(i) ) )'

  return
end
subroutine p13_default ( dim_num )

!*****************************************************************************80
!
!! P13_DEFAULT sets default values for problem 13.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p13_i4 ( 'D', '*', i4 )

  return
end
subroutine p13_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P13_EXACT returns the exact integral for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num

  exact = 0.0D+00

  return
end
subroutine p13_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P13_F evaluates the integrand for problem 13.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    product ( 1 <= i <= dim_num ) t(n(i))(2*x(i)-1)
!
!    where T(N)(X) is the Chebyshev polynomial of order N,
!    and N(I) = mod ( i, 4 ) + 1.
!
!  Exact Integral:
!
!    0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Implementation and Tests of Low-Discrepancy Sequences,
!    ACM Transactions on Modeling and Computer Simulation,
!    Volume 2, Number 3, July 1992, pages 195-213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  real ( kind = 8 ) factor
  integer ( kind = 4 ) k
  integer ( kind = 4 ) point
  real ( kind = 8 ) t
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  value(1:point_num) = 1.0D+00

  do point = 1, point_num

    do dim = 1, dim_num

      t = 2.0D+00 * x(dim,point) - 1.0D+00
      k = mod ( dim, 4 )

      if ( k == 1 ) then
        factor = t
      else if ( k == 2 ) then
        factor = 2.0D+00 * t - 1.0D+00
      else if ( k == 3 ) then
        factor = ( 4.0D+00 * t - 3.0D+00 ) * t
      else if ( k == 4 ) then
        factor = ( 8.0D+00 * t - 8.0D+00 * t + 1.0D+00 )
      end if

      value(point) = value(point) * factor

    end do

  end do

  call p13_i4 ( 'I', '#', point_num )

  return
end
subroutine p13_i4 ( action, name, value )

!*****************************************************************************80
!
!! P13_I4 sets or gets I4 parameters for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P13_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P13_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P13_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P13_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p13_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P13_LIM returns the integration limits for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p13_name ( name )

!*****************************************************************************80
!
!! P13_NAME returns the name of problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'BFN03'

  return
end
subroutine p13_region ( region )

!*****************************************************************************80
!
!! P13_REGION returns the name of the integration region for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p13_title ( )

!*****************************************************************************80
!
!! P13_TITLE prints a title for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 13'
  write ( *, '(a)' ) '  Name:       BFN03'
  write ( *, '(a)' ) '              Bratley, Fox, Niederreiter, #3'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = product ( low order Chebyshevs )'

  return
end
subroutine p14_default ( dim_num )

!*****************************************************************************80
!
!! P14_DEFAULT sets default values for problem 14.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4

  call p14_i4 ( 'D', '*', i4 )

  return
end
subroutine p14_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P14_EXACT returns the exact integral for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num

  exact = ( - 1.0D+00 / 3.0D+00 ) &
    * ( 1.0D+00 - ( -1.0D+00 / 2.0D+00 )**dim_num )

  return
end
subroutine p14_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P14_F evaluates the integrand for problem 14.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    sum ( 1 <= i <= dim_num ) (-1)^i * product ( 1 <= j <= i ) x(j)
!
!  Exact Integral:
!
!    -1/3 ( 1 - (-1/2)^DIM_NUM )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Implementation and Tests of Low-Discrepancy Sequences,
!    ACM Transactions on Modeling and Computer Simulation,
!    Volume 2, Number 3, July 1992, pages 195-213.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  real ( kind = 8 ) factor
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  value(1:point_num) = 0.0D+00

  do point = 1, point_num

    factor = 1.0D+00

    do dim = 1, dim_num

      factor = - factor * x(dim,point)

      value(point) = value(point) + factor

    end do

  end do

  call p14_i4 ( 'I', '#', point_num )

  return
end
subroutine p14_i4 ( action, name, value )

!*****************************************************************************80
!
!! P14_I4 sets or gets I4 parameters for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P14_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P14_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P14_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P14_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p14_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P14_LIM returns the integration limits for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p14_name ( name )

!*****************************************************************************80
!
!! P14_NAME returns the name of problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'BFN04'

  return
end
subroutine p14_region ( region )

!*****************************************************************************80
!
!! P14_REGION returns the name of the integration region for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p14_title ( )

!*****************************************************************************80
!
!! P14_TITLE prints a title for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 14'
  write ( *, '(a)' ) '  Name:       BFN04'
  write ( *, '(a)' ) '              Bratley, Fox, Niederreiter, #4'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = sum ( -1^I * product(X(1:I)) )'

  return
end
subroutine p15_default ( dim_num )

!*****************************************************************************80
!
!! P15_DEFAULT sets default values for problem 15.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) n

  call p15_i4 ( 'D', '*', i4 )

  n = ( 3 * dim_num ) / 2
  call p15_i4 ( 'S', 'N', n )

  return
end
subroutine p15_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P15_EXACT returns the exact integral for problem 15.
!
!  Discussion:
!
!    Thanks to Jeffrey Sax of Extreme Optimization for suggesting a revision 
!    of this routine to improve its clarity, 13 August 2008.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  real ( kind = 8 ) exact
  integer ( kind = 4 ) fraction
  integer ( kind = 4 ) n
  integer ( kind = 4 ) remainder

  call p15_i4 ( 'G', 'N', n )

  fraction = n / dim_num
  remainder = n - fraction * dim_num

  exact = 1.0D+00
  do dim = 1, dim_num
    if ( dim <= remainder ) then
      exact = exact / real ( fraction + 2, kind = 8 )
    else
      exact = exact / real ( fraction + 1, kind = 8 )
    end if
  end do

  return
end
subroutine p15_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P15_F evaluates the integrand for problem 15.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    N determines the order of the product.  It defaults to 1.
!    You can modify N by calling P15_I4.
!
!  Integrand:
!
!    f(x) = product ( 1 <= I <= N ) X(MOD(I-1,DIM_NUM)+1)
!
!  Exact integral:
!
!    product ( 1 / exponent(1:DIM_NUM) )
!
!    where, if I <= N - DIM_NUM * ( N/DIM_NUM),
!    
!      exponent ( I ) = ( N / dim_num ) + 2 
!
!    else
!
!      exponent ( I ) = ( N / dim_num ) + 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) fraction
  integer ( kind = 4 ) n
  integer ( kind = 4 ) point
  integer ( kind = 4 ) remainder
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p15_i4 ( 'G', 'N', n )

  fraction = n / dim_num
  remainder = n - fraction * dim_num

  value(1:point_num) = 1.0D+00

  do point = 1, point_num

    if ( any ( x(1:dim_num,point) == 0.0D+00 ) ) then

      value(point) = 0.0D+00

    else 

      do dim = 1, dim_num
        if ( dim <= remainder ) then
          value(point) = value(point) * x(dim,point) ** ( fraction + 1 )
        else if ( fraction /= 0 ) then
          value(point) = value(point) * x(dim,point) ** fraction
        end if
      end do

    end if

  end do

  call p15_i4 ( 'I', '#', point_num )

  return
end
subroutine p15_i4 ( action, name, value )

!*****************************************************************************80
!
!! P15_I4 sets or gets I4 parameters for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!    'N' is the number of factors to multiply in the integrand.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ), save :: n = 1
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

    if ( name == 'n' .or. name == 'N' .or. name == '*' ) then
      n = 1
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else if ( name == 'N' .or. name == 'n' ) then
      value = n
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P15_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else if ( name == 'N' .or. name == 'n' ) then
      n = n + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P15_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else if ( name == 'N' .or. name == 'n' ) then
      n = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P15_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P15_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p15_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P15_LIM returns the integration limits for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p15_name ( name )

!*****************************************************************************80
!
!! P15_NAME returns the name of problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Partial product ( X(1:N) )'

  return
end
subroutine p15_region ( region )

!*****************************************************************************80
!
!! P15_REGION returns the name of the integration region for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p15_title ( )

!*****************************************************************************80
!
!! P15_TITLE prints a title for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 15'
  write ( *, '(a)' ) '  Name:       Partial product ( X(1:N) )'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = product ( X(1:N) )'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              N, defaults to 1'

  return
end
subroutine p16_default ( dim_num )

!*****************************************************************************80
!
!! P16_DEFAULT sets default values for problem 16.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p16_i4 ( 'D', '*', i4 )
  call p16_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p16_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P16_EXACT returns the exact integral for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) exact
  real ( kind = 8 ) integral
  real ( kind = 8 ) volume
  real ( kind = 8 ) z(dim_num)
!
!  Get the limits of integration.
!
  call p16_lim ( dim_num, a, b )
!
!  Get the location of Z.
!
  call p16_r8vec ( 'G', 'Z', dim_num, z )
!
!  The value of the DIM_NUM dimensional integral can be broken down
!  into the weighted sum of 1 dimensional integrals.
!
  exact = 0.0D+00

  volume = product ( b(1:dim_num) - a(1:dim_num) )

  do dim = 1, dim_num
!
!  Z < A < B
!
    if ( z(dim) < a(dim) ) then

      integral = &
        0.5D+00 * ( b(dim) - a(dim) ) * ( b(dim) + a(dim) - 2.0D+00 * z(dim) )
!
!  A < Z < B
!
    else if ( z(dim) < b(dim) ) then

      integral = &
        0.5D+00 * ( ( b(dim) - z(dim) )**2 + ( z(dim) - a(dim) )**2 )
!
!  A < B < Z
!
    else

      integral = &
        0.5D+00 * ( b(dim) - a(dim) ) * ( 2.0D+00 * z(dim) - a(dim) - b(dim) )

    end if

    exact = exact + volume * integral / ( b(dim) - a(dim) )

  end do

  return
end
subroutine p16_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P16_F evaluates the integrand for problem 16.
!
!  Discussion:
!
!    The integrand can be regarded as the L1 norm of X - Z.
!
!    It would be nice to allow the use to specify several
!    base points Z, to make the function more jagged more places!
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    The integrand can be regarded as the L1 norm of X - Z.
!
!    There is a basis point Z associated with the integrand.
!    Z(1:DIM_NUM) defaults to ( 0.5, 0.5, ..., 0.5 ).
!    The user can set, get, or randomize this value by calling
!    P16_R8VEC.
!
!  Integrand:
!
!    sum ( abs ( x(1:dim_num) - z(1:dim_num) ) )
!
!  Exact Integral:
!
!    The integral is separable into
!
!       Int ( A(1) <= X(1) <= B(1) ) abs ( X(1) - Z(1) ) 
!     * Product ( B(1:N)-A(1:N), skip index 1 )
!     + Int ( A(2) <= X(2) <= B(2) ) abs ( X(2) - Z(2) )
!     * Product ( B(1:N)-A(1:N), skip index 2 )
!     ...
!     + Int ( A(N) <= X(N) <= B(N) ) abs ( X(N) - Z(N) )
!     * Product ( B(1:N)-A(1:N), skip index N )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ) z(dim_num)

  call p16_r8vec ( 'G', 'Z', dim_num, z )

  do point = 1, point_num
    value(point) = sum ( abs ( x(1:dim_num,point) - z(1:dim_num) ) )
  end do

  call p16_i4 ( 'I', '#', point_num )

  return
end
subroutine p16_i4 ( action, name, value )

!*****************************************************************************80
!
!! P16_I4 sets or gets I4 parameters for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P16_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P16_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P16_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P16_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p16_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P16_LIM returns the integration limits for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p16_name ( name )

!*****************************************************************************80
!
!! P16_NAME returns the name of problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'L1(X-Z)'

  return
end
subroutine p16_region ( region )

!*****************************************************************************80
!
!! P16_REGION returns the name of the integration region for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p16_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P16_R8VEC sets or gets R8VEC parameters for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!      be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'Z' is the base vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: z

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( z ) ) then
      deallocate ( z )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( z(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'Z' .or. name == 'z' .or. name == '*' ) then
      z(1:dim_num) = 0.5D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then
      
    if ( name == 'Z' .or. name == 'z' ) then
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P16_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'Z' .or. name == 'z' ) then
      call random_number ( harvest = z(1:dim_num) )
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P16_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'Z' .or. name == 'z' ) then
      z(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P16_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P16_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p16_title ( )

!*****************************************************************************80
!
!! P16_TITLE prints a title for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 16'
  write ( *, '(a)' ) '  Name:       L1(X-Z)'
  write ( *, '(a)' ) '              Lipschitz continuous.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = sum ( | X(i) - Z(i) | )'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              Z(1:DIM_NUM) defaults to (0.5,0.5,...)'

  return
end
subroutine p17_default ( dim_num )

!*****************************************************************************80
!
!! P17_DEFAULT sets default values for problem 17.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p17_i4 ( 'D', '*', i4 )
  call p17_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p17_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P17_EXACT returns the exact integral for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) exact
  real ( kind = 8 ) integral
  real ( kind = 8 ) volume
  real ( kind = 8 ) z(dim_num)
!
!  Get the limits of integration.
!
  call p17_lim ( dim_num, a, b )
!
!  Get the location of Z.
!
  call p17_r8vec ( 'G', 'Z', dim_num, z )
!
!  The value of the DIM_NUM dimensional integral can be broken
!  into the weighted sum of integrals over each dimension.
!
  exact = 0.0D+00

  volume = product ( b(1:dim_num) - a(1:dim_num) )

  do dim = 1, dim_num
!
!  Z < A < B
!
    if ( z(dim) < a(dim) ) then

      integral = ( ( b(dim) - z(dim) )**3 &
                 - ( a(dim) - z(dim) )**3 ) / 3.0D+00
!
!  A < Z < B
!
    else if ( z(dim) < b(dim) ) then

      integral = ( ( b(dim) - z(dim) )**3 &
                 + ( z(dim) - a(dim) )**3 ) / 3.0D+00
!
!  A < B < Z
!
    else

      integral = ( ( z(dim) - b(dim) )**3 &
                 - ( z(dim) - a(dim) )**3 ) / 3.0D+00

    end if

    exact = exact + volume * integral / ( b(dim) - a(dim) )

  end do

  return
end
subroutine p17_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P17_F evaluates the integrand for problem 17.
!
!  Discussion:
!
!    This integrand can be regarded as the square of the L2
!    norm of X - Z.
!
!    This integrand has the advantage of symmetry under rotation
!    about Z.  Thus, it is possible to test whether a particular
!    quadrature rule has an occasional advantage because it
!    "lines up" with the X and Y coordinate axes and hence can
!    integrate some separable integrals very well.
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    There is a basis point Z associated with the integrand.
!    Z(1:DIM_NUM) defaults to ( 0.5, 0.5, ..., 0.5 ).
!    The user can set, get, or randomize this value by calling
!    P17_R8VEC.
!
!  Integrand:
!
!    sum ( ( x(1:dim_num) - z(1:dim_num) )^2 )
!
!  Exact Integral:
!
!    The integral may be broken into the sum of weighted 
!    one dimensional integrals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ) z(dim_num)

  call p17_r8vec ( 'G', 'Z', dim_num, z )

  do point = 1, point_num
    value(point) = sum ( ( x(1:dim_num,point) - z(1:dim_num) )**2 )
  end do

  call p17_i4 ( 'I', '#', point_num )

  return
end
subroutine p17_i4 ( action, name, value )

!*****************************************************************************80
!
!! P17_I4 sets or gets I4 parameters for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P17_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P17_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P17_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P17_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p17_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P17_LIM returns the integration limits for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p17_name ( name )

!*****************************************************************************80
!
!! P17_NAME returns the name of problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'L2(X-Z)^2'

  return
end
subroutine p17_region ( region )

!*****************************************************************************80
!
!! P17_REGION returns the name of the integration region for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p17_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P17_R8VEC sets or gets R8VEC parameters for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!      be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'Z' is the base vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: z

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( z ) ) then
      deallocate ( z )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( z(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'Z' .or. name == 'z' .or. name == '*' ) then
      z(1:dim_num) = 0.5D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then
      
    if ( name == 'Z' .or. name == 'z' ) then
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P17_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'Z' .or. name == 'z' ) then
      call random_number ( harvest = z(1:dim_num) )
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P17_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'Z' .or. name == 'z' ) then
      z(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P17_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P17_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p17_title ( )

!*****************************************************************************80
!
!! P17_TITLE prints a title for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 17'
  write ( *, '(a)' ) '  Name:       L2(X-Z)^2'
  write ( *, '(a)' ) '              Zero at point Z.'
  write ( *, '(a)' ) '              Radially symmetric.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = sum ( ( X(i) - Z(i) )^2 )'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              Z(1:DIM_NUM) defaults to (0.5,0.5,...)'

  return
end
subroutine p18_default ( dim_num )

!*****************************************************************************80
!
!! P18_DEFAULT sets default values for problem 18.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8
  real ( kind = 8 ) r8vec(1)

  call p18_i4 ( 'D', '*', i4 )
  call p18_r8 ( 'D', '*', r8 )
  call p18_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p18_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P18_EXACT returns the exact integral for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) r

  call p18_r8 ( 'G', 'R', r )

  call sphere_volume_nd ( dim_num, r, exact )

  return
end
subroutine p18_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P18_F evaluates the integrand for problem 18.
!
!  Discussion:
!
!    This is the characteristic function of the interior of the 
!    N sphere of radius R and center Z, to be integrated within the
!    unit hypercube [0,1]^N.  If the user picks a combination of R
!    and Z that causes the volume of the sphere to lie at least
!    partially outside the unit hypercube, the formula for the
!    exact integral will no longer be correct.
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    R defaults to 0.50.  
!    You can change R by calling P18_R8.
!
!    Z(1:DIM_NUM) defaults to (0.5,0.5,...0.5).  
!    You can change Z by calling P18_R8VEC.
!
!  Integrand:
!
!    f(x) = 1 if X(1:DIM_NUM) is less than R from Z(1:DIM_NUM),
!           0 otherwise.
!
!  Exact Integral:
!
!    The volume of the DIM_NUM sphere of radius R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) d
  integer ( kind = 4 ) point
  real ( kind = 8 ) r
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ) z(dim_num)

  call p18_r8vec ( 'G', 'Z', dim_num, z )

  call p18_r8 ( 'G', 'R', r )

  do point = 1, point_num

    d = sqrt ( sum ( ( x(1:dim_num,point) - z(1:dim_num) )**2 ) )

    if ( d <= r ) then
      value(point) = 1.0D+00
    else
      value(point) = 0.0D+00
    end if

  end do

  call p18_i4 ( 'I', '#', point_num )

  return
end
subroutine p18_i4 ( action, name, value )

!*****************************************************************************80
!
!! P18_I4 sets or gets I4 parameters for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P18_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p18_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P18_LIM returns the integration limits for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00
 
  return
end
subroutine p18_name ( name )

!*****************************************************************************80
!
!! P18_NAME returns the name of problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Disk'

  return
end
subroutine p18_r8 ( action, name, value )

!*****************************************************************************80
!
!! P18_R8 sets or gets R8 parameters for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' sets the parameter to its default value;
!    'G' gets a parameter.
!    'R' sets the parameter to a random value.
!    'S' sets a parameter,
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    'R' is the radius.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G' or 'R', then VALUE is an output quantity, 
!      and is the current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  real ( kind = 8 ), save :: r = 0.50D+00
  real ( kind = 8 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'R' .or. name == 'r' .or. name == '*' ) then
      r = 0.5D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'R' .or. name == 'r' ) then
      value = r
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'R' .or. name == 'r' ) then
      call random_number ( harvest = r )
      value = r
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'R' .or. name == 'r' ) then
      r = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P18_R8 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p18_region ( region )

!*****************************************************************************80
!
!! P18_REGION returns the name of the integration region for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p18_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P18_R8VEC sets or gets R8VEC parameters for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!      be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'Z' is the base vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: z

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( z ) ) then
      deallocate ( z )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( z(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'Z' .or. name == 'z' .or. name == '*' ) then
      z(1:dim_num) = 0.5D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'Z' .or. name == 'z' ) then
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'Z' .or. name == 'z' ) then
      call random_number ( harvest = z(1:dim_num) )
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'Z' .or. name == 'z' ) then
      z(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P18_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p18_title ( )

!*****************************************************************************80
!
!! P18_TITLE prints a title for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 18'
  write ( *, '(a)' ) '  Name:       Disk'
  write ( *, '(a)' ) '              Disk of radius R centered at Z.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = sphere interior characteristic'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              R, defaults to 0.5'
  write ( *, '(a)' ) '              Z(1:DIM_NUM) defaults to (0.5,0.5,...0.5)'

  return
end
subroutine p19_default ( dim_num )

!*****************************************************************************80
!
!! P19_DEFAULT sets default values for problem 19.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p19_i4 ( 'D', '*', i4 )
  call p19_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p19_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P19_EXACT returns the exact integral for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  real ( kind = 8 ) exact
  real ( kind = 8 ) z(dim_num)

  call p19_r8vec ( 'G', 'Z', dim_num, z )

  exact = 1.0D+00

  do dim = 1, dim_num 

    if ( z(dim) <= 0.0D+00 ) then
      exact = exact * ( 2.0D+00 / 3.0D+00 ) * &
        ( sqrt ( ( 1.0D+00 - z(dim) )**3 ) - sqrt ( - ( z(dim)**3 ) ) )
    else if ( z(dim) < 1.0D+00 ) then
      exact = exact * ( 2.0D+00 / 3.0D+00 ) * &
        ( sqrt ( z(dim)**3 ) + sqrt ( ( 1.0D+00 - z(dim) )**3 ) )
    else if ( 1.0D+00 <= z(dim) ) then
      exact = exact * ( 2.0D+00 / 3.0D+00 ) * &
        ( sqrt ( z(dim)**3 ) - sqrt ( ( z(dim) - 1.0D+00 )**3 ) )
    end if

  end do

  return
end
subroutine p19_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P19_F evaluates the integrand for problem 19.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    Z defaults to (1/3,1/3,...,1/3).  
!    You can reset Z by calling P19_R8VEC.
!
!  Integrand:
!
!    f(x) = product ( sqrt ( abs ( x(1:dim_num) - z(1:dim_num) ) ) )
!
!  Exact Integral:
!
!    With Z as given, 
!
!      (2/3)^DIM_NUM * ( (2/3)^(3/2) + (1/3)^(3/2) )^DIM_NUM
!
!    or approximately 0.49^DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arnold Krommer, Christoph Ueberhuber,
!    Numerical Integration on Advanced Systems,
!    Springer, 1994,
!    ISBN: 3540584102,
!    LC: QA299.3.K76.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ) z(dim_num)

  call p19_r8vec ( 'G', 'Z', dim_num, z )

  do point = 1, point_num
    value(point) = product ( &
      sqrt ( abs ( x(1:dim_num,point) - z(1:dim_num) ) ) )
  end do

  call p19_i4 ( 'I', '#', point_num )

  return
end
subroutine p19_i4 ( action, name, value )

!*****************************************************************************80
!
!! P19_I4 sets or gets I4 parameters for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P19_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P19_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P19_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P19_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p19_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P19_LIM returns the integration limits for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  return
end
subroutine p19_name ( name )

!*****************************************************************************80
!
!! P19_NAME returns the name of problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Sqrt-Prod'

  return
end
subroutine p19_region ( region )

!*****************************************************************************80
!
!! P19_REGION returns the name of the integration region for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p19_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P19_R8VEC sets or gets R8VEC parameters for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!      be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'Z' is the base vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: z

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( z ) ) then
      deallocate ( z )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( z(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'Z' .or. name == 'z' .or. name == '*' ) then
      z(1:dim_num) = 1.0D+00 / 3.0D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'Z' .or. name == 'z' ) then
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P19_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then
      
    if ( name == 'Z' .or. name == 'z' ) then
      call random_number ( harvest = z(1:dim_num) )
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P19_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'Z' .or. name == 'z' ) then
      z(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P19_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P19_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p19_title ( )

!*****************************************************************************80
!
!! P19_TITLE prints a title for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 19'
  write ( *, '(a)' ) '  Name:       Sqrt-Prod'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  ' // &
    'F(X) = prod ( sqrt ( | X(i) - Z(i) | ) )'
  write ( *, '(a)' ) '  Parameters: '
  write ( *, '(a)' ) '              Z(1:DIM_NUM) defaults to (1/3,1/3,...,1/3)'

  return
end
subroutine p20_default ( dim_num )

!*****************************************************************************80
!
!! P20_DEFAULT sets default values for problem 20.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8

  call p20_i4 ( 'D', '*', i4 )
  call p20_r8 ( 'D', '*', r8 )

  return
end
subroutine p20_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P20_EXACT returns the exact integral for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) aval
  real ( kind = 8 ) bval
  integer ( kind = 4 ) dim
  real ( kind = 8 ) exact
  real ( kind = 8 ) exponent
  real ( kind = 8 ) minus_one
  real ( kind = 8 ) p
  real ( kind = 8 ) r8_choose

  call p20_r8 ( 'G', 'A', aval )
  call p20_r8 ( 'G', 'B', bval )
  call p20_r8 ( 'G', 'P', p )

  exact = 0.0D+00
  exponent = real ( dim_num + p, kind = 8 )

  minus_one = -1.0D+00
  do dim = 0, dim_num
    minus_one = - minus_one
    exact = exact + minus_one * r8_choose ( dim_num, dim ) &
      * ( real ( dim_num - dim, kind = 8 ) * bval + real ( dim, kind = 8 ) &
      * aval )**exponent
  end do

  do dim = 1, dim_num
    exact = exact / ( p + real ( dim, kind = 8 ) )
  end do

  return
end
subroutine p20_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P20_F evaluates the integrand for problem 20.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    A <= X(1:DIM_NUM) <= B
!
!  Integrand:
!
!    f(x) = ( sum ( x(1:dim_num) ) )^p
!
!    P is greater than -dim_num, and is not a negative integer.
!
!  Exact Integral:
!
!    sum ( 0 <= i <= dim_num ) (-1)^i * choose(dim_num,i) 
!      * ((dim_num-i)*b+i*a)^(dim_num+p)
!      / ( (t+1) * (t+2) * ... * (t+dim_num) )
!
!  Parameters:
!
!    A defaults to 0.0.  You can change A by calling P20_R8.
!
!    B defaults to 1.0.  You can change B by calling P20_R8.
!
!    P defaults to 2.0.  You can change P by calling P20_R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) p
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p20_r8 ( 'G', 'P', p )

  do point = 1, point_num
    value(point) = ( sum ( x(1:dim_num,point) ) )**p
  end do

  call p20_i4 ( 'I', '#', point_num )

  return
end
subroutine p20_i4 ( action, name, value )

!*****************************************************************************80
!
!! P20_I4 sets or gets I4 parameters for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P20_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p20_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P20_LIM returns the integration limits for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) aval
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ) bval

  call p20_r8 ( 'G', 'A', aval )
  call p20_r8 ( 'G', 'B', bval )

  a(1:dim_num) = aval
  b(1:dim_num) = bval

  return
end
subroutine p20_name ( name )

!*****************************************************************************80
!
!! P20_NAME returns the name of problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Sum^P'

  return
end
subroutine p20_r8 ( action, name, value )

!*****************************************************************************80
!
!! P20_R8 sets or gets R8 parameters for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' sets the parameter to its default value;
!    'G' gets a parameter.
!    'R' sets the parameter to a random value.
!    'S' sets a parameter,
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'A' is the value of all lower integral bounds.
!    'B' is the value of all upper integral bounds.
!    'P' is the value of the exponent in the integrand.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G' or 'R', then VALUE is an output quantity, 
!      and is the current value of NAME.
!
  implicit none

  real ( kind = 8 ), save :: a = 0.0D+00
  character ( len = * ) action
  real ( kind = 8 ), save :: b = 1.0D+00
  character ( len = * ) name
  real ( kind = 8 ) value
  real ( kind = 8 ), save :: p = 2.0D+00

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'A' .or. name == 'a' .or. name == '*' ) then
      a = 0.0D+00
    end if

    if ( name == 'B' .or. name == 'b' .or. name == '*' ) then
      b = 1.0D+00
    end if

    if ( name == 'P' .or. name == 'p' .or. name == '*' ) then
      p = 2.0D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'A' .or. name == 'a' ) then
      value = a
    else if ( name == 'B' .or. name == 'b' ) then
      value = b
    else if ( name == 'P' .or. name == 'p' ) then
      value = p
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'A' .or. name == 'a' ) then
      call random_number ( harvest = a )
      value = a
    else if ( name == 'B' .or. name == 'b' ) then
      call random_number ( harvest = b )
      value = b
    else if ( name == 'P' .or. name == 'p' ) then
      call random_number ( harvest = p )
      value = p
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'A' .or. name == 'a' ) then
      a = value
    else if ( name == 'B' .or. name == 'b' ) then
      b = value
    else if ( name == 'P' .or. name == 'p' ) then
      p = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P20_R8 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p20_region ( region )

!*****************************************************************************80
!
!! P20_REGION returns the name of the integration region for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p20_title ( )

!*****************************************************************************80
!
!! P20_TITLE prints a title for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 20'
  write ( *, '(a)' ) '  Name:       Sum^P'
  write ( *, '(a)' ) '  Region:     A <= X(i) <= B'
  write ( *, '(a)' ) '  Integrand:  F(X) = ( sum ( X(i) ) )^p'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              A, defaults to 0.0.'
  write ( *, '(a)' ) '              B, defaults to 1.0.'
  write ( *, '(a)' ) '              P, defaults to 2.0.'

  return
end
subroutine p21_default ( dim_num )

!*****************************************************************************80
!
!! P21_DEFAULT sets default values for problem 21.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4vec(1)
  real ( kind = 8 ) r8

  call p21_i4 ( 'D', '*', i4 )
  call p21_r8 ( 'D', '*', r8 )
  call p21_i4vec ( 'D', '*', dim_num, i4vec )

  return
end
subroutine p21_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P21_EXACT returns the exact integral for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) arg
  real ( kind = 8 ) c
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) e(dim_num)
  real ( kind = 8 ) exact
  real ( kind = 8 ) r8_gamma

  call p21_r8 ( 'G', 'C', c )
  call p21_i4vec ( 'G', 'E', dim_num, e )

  if ( any ( mod ( e(1:dim_num), 2 ) == 1 ) ) then
    exact = 0.0D+00
    return
  end if

  exact = 2.0D+00 * c
  do dim = 1, dim_num
    arg = real ( e(dim) + 1, kind = 8 ) / 2.0D+00
    exact = exact * r8_gamma ( arg )
  end do

  arg = real ( sum ( e(1:dim_num) ) + dim_num, kind = 8 ) / 2.0D+00

  exact = exact / r8_gamma ( arg )

  return
end
subroutine p21_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P21_F evaluates the integrand for problem 21.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    The (surface of the) unit sphere
!
!  Integral Parameters:
!
!    C defaults to 1.  
!    Call P21_R8 to get or set this value.
!
!    E(1:DIM_NUM) defaults to (/ 2, 2, ..., 2 /).  
!    Call P21_I4VEC to get or set this value.
!
!  Integrand:
!
!    F(X) = C * X1^E1 * X2^E2 * ... * Xn^En
!
!    C is real, all exponents E are nonnegative integers.
!
!  Exact Integral:
!
!    0, if any exponent is odd.
!    2 * C * Gamma((E1+1)/2) * Gamma((E2+1)/2) * ... * Gamma((En+1)/2) 
!      / Gamma( (E1+E2+...+En+N)/2 ), otherwise.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Gerald Folland,
!    How to Integrate a Polynomial Over a Sphere,
!    American Mathematical Monthly, 
!    May 2001, pages 446-448.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c
  integer ( kind = 4 ) dim
  real ( kind = 8 ) e(dim_num)
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p21_r8 ( 'G', 'C', c )
  call p21_i4vec ( 'G', 'E', dim_num, e )

  value(1:point_num) = c

  do dim = 1, dim_num
    value(1:point_num) = value(1:point_num) * x(dim,1:point_num)**e(dim)
  end do

  call p21_i4 ( 'I', '#', point_num )

  return
end
subroutine p21_i4 ( action, name, value )

!*****************************************************************************80
!
!! P21_I4 sets or gets I4 parameters for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P21_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P21_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P21_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P21_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p21_i4vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P21_I4VEC sets or gets I4VEC parameters for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the object to a default value.
!    'G' means the current value of the object should be returned.
!    'S' means the input values of the object and its dimension
!    should be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'E' is the exponent vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, integer ( kind = 4 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  integer ( kind = 4 ), allocatable, save, dimension ( : ) :: e
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  integer ( kind = 4 ) value(dim_num)

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( e ) ) then
      deallocate ( e )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( e(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'E' .or. name == 'e' .or. name == '*' ) then
      e(1:dim_num) = 2
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'E' .or. name == 'e' ) then
      value(1:dim_num) = e(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P21_I4VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'E' .or. name == 'e' ) then
      e(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P21_I4VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P21_I4VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p21_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P21_LIM returns the integration limits for problem 21.
!
!  Discussion:
!
!    Because the integration region is the surface of the unit sphere,
!    the integration limits simply specify the limits of a box
!    containing the integration region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) = -1.0D+00
  b(1:dim_num) = +1.0D+00

  return
end
subroutine p21_name ( name )

!*****************************************************************************80
!
!! P21_NAME returns the name of problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'SphereMonomial'

  return
end
subroutine p21_r8 ( action, name, value )

!*****************************************************************************80
!
!! P21_R8 sets or gets R8 parameters for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' sets the parameter to its default value;
!    'G' gets a parameter.
!    'R' sets the parameter to a random value.
!    'S' sets a parameter,
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the coefficient.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G' or 'R', then VALUE is an output quantity, 
!      and is the current value of NAME.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: c = 1.0D+00
  character ( len = * ) name
  real ( kind = 8 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c = 1.0D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value = c
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P21_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c )
      value = c
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P21_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P21_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P21_R8 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p21_region ( region )

!*****************************************************************************80
!
!! P21_REGION returns the name of the integration region for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'sphere'

  return
end
subroutine p21_title ( )

!*****************************************************************************80
!
!! P21_TITLE prints a title for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 21'
  write ( *, '(a)' ) '  Name:       SphereMonomial'
  write ( *, '(a)' ) '  Region:     Sphere surface, radius 1, center 0'
  write ( *, '(a)' ) '  Integrand:  F(X) = C * product ( X(i)^E(i) )'
  write ( *, '(a)' ) '  Parameters:' 
  write ( *, '(a)' ) '              C, defaults to 1.0.'
  write ( *, '(a)' ) '              E(1:DIM_NUM) defaults to 2.'

  return
end
subroutine p22_default ( dim_num )

!*****************************************************************************80
!
!! P22_DEFAULT sets default values for problem 22.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4vec(1)
  real ( kind = 8 ) r8

  call p22_i4 ( 'D', '*', i4 )
  call p22_r8 ( 'D', '*', r8 )
  call p22_i4vec ( 'D', '*', dim_num, i4vec )

  return
end
subroutine p22_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P22_EXACT returns the exact integral for problem 22.
!
!  Discussion:
!
!    Thanks to Jeffrey Sax of Extreme Optimization for pointing out a mistake
!    in a previous version of this routine, 28 May 2008.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) arg
  real ( kind = 8 ) c
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) e(dim_num)
  integer ( kind = 4 ) e_sum
  real ( kind = 8 ) exact
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_gamma

  call p22_r8 ( 'G', 'C', c )
  call p22_r8 ( 'G', 'R', r )
  call p22_i4vec ( 'G', 'E', dim_num, e )

  if ( any ( mod ( e(1:dim_num), 2 ) == 1 ) ) then
    exact = 0.0D+00
    return
  end if

  e_sum = sum ( e(1:dim_num) ) + dim_num

  exact = 2.0D+00 * c
  do dim = 1, dim_num
    arg = real ( e(dim) + 1, kind = 8 ) / 2.0D+00
    exact = exact * r8_gamma ( arg )
  end do

  arg = real ( e_sum, kind = 8 ) / 2.0D+00

  exact = exact / r8_gamma ( arg )

  exact = exact * r**( e_sum ) / real ( e_sum, kind = 8 )

  return
end
subroutine p22_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P22_F evaluates the integrand for problem 22.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    The interior of a sphere of radius R centered at the origin.
!
!  Integral Parameters:
!
!    C defaults to 1.0.  
!    Call P22_R8 to get or set this value.
!
!    R defaults to 1.0.  
!    Call P22_R8 to get or set this value.
!
!    E(1:DIM_NUM) defaults to (/ 2, 2, ..., 2 /).  
!    Call P22_I4VEC to get or set this value.
!
!  Integrand:
!
!    F(X) = C * X1^E1 * X2^E2 * ... * Xn^En
!
!    C is real, all exponents E are nonnegative integers.
!
!  Exact Integral:
!
!    0, if any exponent is odd.
!    2 * C * R^(E1+E2+...+EN+N) 
!      * Gamma((E1+1)/2) * Gamma((E2+1)/2) * ... * Gamma((En+1)/2) 
!      / ( Gamma( (E1+E2+...+En+N)/2 ) * ( E1+E2+...+EN+N) ), otherwise.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Gerald Folland,
!    How to Integrate a Polynomial Over a Sphere,
!    American Mathematical Monthly, 
!    May 2001, pages 446-448.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c
  integer ( kind = 4 ) dim
  real ( kind = 8 ) e(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p22_r8 ( 'G', 'C', c )
  call p22_i4vec ( 'G', 'E', dim_num, e )

  value(1:point_num) = c

  do point = 1, point_num
    do dim = 1, dim_num
      value(point) = value(point) * x(dim,point)**e(dim)
    end do
  end do

  call p22_i4 ( 'I', '#', point_num )

  return
end
subroutine p22_i4 ( action, name, value )

!*****************************************************************************80
!
!! P22_I4 sets or gets I4 parameters for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P22_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P22_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P22_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P22_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p22_i4vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P22_I4VEC sets or gets I4VEC parameters for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the object to a default value.
!    'G' means the current value of the object should be returned.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'E' is the exponent vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, integer ( kind = 4 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  integer ( kind = 4 ), allocatable, save, dimension ( : ) :: e
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  integer ( kind = 4 ) value(dim_num)

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( e ) ) then
      deallocate ( e )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( e(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'E' .or. name == 'e' .or. name == '*' ) then
      e(1:dim_num) = 2
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'E' .or. name == 'e' ) then
      value(1:dim_num) = e(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P22_I4VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'E' .or. name == 'e' ) then
      e(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P22_I4VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P22_I4VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p22_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P22_LIM returns the integration limits for problem 22.
!
!  Discussion:
!
!    Because the integration region is the interior of a sphere
!    of radius R centered at the origin, the integration limits simply 
!    specify the limits of a box containing the integration region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ) r

  call p22_r8 ( 'G', 'R', r )

  a(1:dim_num) = -r
  b(1:dim_num) = +r

  return
end
subroutine p22_name ( name )

!*****************************************************************************80
!
!! P22_NAME returns the name of problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'BallMonomial'

  return
end
subroutine p22_r8 ( action, name, value )

!*****************************************************************************80
!
!! P22_R8 sets or gets R8 parameters for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' sets the parameter to its default value;
!    'G' gets a parameter.
!    'R' sets the parameter to a random value.
!    'S' sets a parameter,
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the coefficient.
!    'R' is the radius of the sphere.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'D' or 'G' or 'R', then VALUE is an output quantity, 
!      and is the current value of NAME.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: c = 1.0D+00
  character ( len = * ) name
  real ( kind = 8 ), save :: r = 1.0D+00
  real ( kind = 8 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c = 1.0D+00
    end if

    if ( name == 'R' .or. name == 'r' .or. name == '*' ) then
      r = 1.0D+00
     end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value = c
    else if ( name == 'R' .or. name == 'r' ) then
      value = r
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P22_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c )
      value = c
    else if ( name == 'R' .or. name == 'r' ) then
      call random_number ( harvest = r )
      value = r
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P22_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c = value
    else if ( name == 'R' .or. name == 'r' ) then
      r = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P22_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P22_R8 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p22_region ( region )

!*****************************************************************************80
!
!! P22_REGION returns the name of the integration region for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'ball'

  return
end
subroutine p22_title ( )

!*****************************************************************************80
!
!! P22_TITLE prints a title for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 22'
  write ( *, '(a)' ) '  Name:       BallMonomial'
  write ( *, '(a)' ) '  Region:     Sphere interior, radius R, center 0'
  write ( *, '(a)' ) '  Integrand:  F(X) = C * product ( X(i)^E(i) )'
  write ( *, '(a)' ) '  Parameters:' 
  write ( *, '(a)' ) '              C, defaults to 1.0.'
  write ( *, '(a)' ) '              R, defaults to 1.0.'
  write ( *, '(a)' ) '              E(1:DIM_NUM) defaults to 2.'

  return
end
subroutine p23_default ( dim_num )

!*****************************************************************************80
!
!! P23_DEFAULT sets default values for problem 23.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4vec
  real ( kind = 8 ) r8

  call p23_i4 ( 'D', '*', i4 )
  call p23_r8 ( 'D', '*', r8 )
  call p23_i4vec ( 'D', '*', dim_num, i4vec )

  return
end
subroutine p23_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P23_EXACT returns the exact integral for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) c
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) e(dim_num)
  real ( kind = 8 ) exact
  real ( kind = 8 ) r8_gamma

  call p23_r8 ( 'G', 'C', c )
  call p23_i4vec ( 'G', 'E', dim_num, e )

  exact = c
  do dim = 1, dim_num
    exact = exact * r8_gamma ( real ( e(dim) + 1, kind = 8 ) )
  end do

  exact = exact / r8_gamma ( real ( sum ( e(1:dim_num) ) + 1, kind = 8 ) )

  return
end
subroutine p23_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P23_F evaluates the integrand for problem 23.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    The interior of the unit simplex, for which all X's are nonnegative,
!    and sum ( X(1:N) ) <= 1.
!
!  Integral Parameters:
!
!    C defaults to 1.0.  
!    Call P23_R8 to get or set this value.
!
!    E(1:DIM_NUM) defaults to (/ 2, 2, ..., 2 /).  
!    Call P23_I4VEC to get or set this value.
!
!  Integrand:
!
!    F(X) = C * X1^E1 * X2^E2 * ... * Xn^En
!
!    C is real, all exponents E are nonnegative integers.
!
!  Exact Integral:
!
!    C * Gamma(E1+1) * Gamma(E2+1) * ... * Gamma(En+1) / Gamma(E1+E2+...+En+1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c
  integer ( kind = 4 ) dim
  real ( kind = 8 ) e(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p23_r8 ( 'G', 'C', c )
  call p23_i4vec ( 'G', 'E', dim_num, e )

  value(1:point_num) = c

  do point = 1, point_num
    do dim = 1, dim_num
      value(point) = value(point) * x(dim,point)**e(dim)
    end do
  end do

  call p23_i4 ( 'I', '#', point_num )

  return
end
subroutine p23_i4 ( action, name, value )

!*****************************************************************************80
!
!! P23_I4 sets or gets I4 parameters for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P23_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P23_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P23_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P23_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p23_i4vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P23_I4VEC sets or gets I4VEC parameters for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the object to a default value.
!    'G' means the current value of the object should be returned.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'E' is the exponent vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, integer ( kind = 4 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  integer ( kind = 4 ), allocatable, save, dimension ( : ) :: e
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  integer ( kind = 4 ) value(dim_num)

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( e ) ) then
      deallocate ( e )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( e(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'E' .or. name == 'e' .or. name == '*' ) then
      e(1:dim_num) = 2
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'E' .or. name == 'e' ) then
      value(1:dim_num) = e(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P23_I4VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'E' .or. name == 'e' ) then
      e(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P23_I4VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P23_I4VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p23_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P23_LIM returns the integration limits for problem 23.
!
!  Discussion:
!
!    Because the integration region is the interior of the unit simplex,
!    the integration limits simply specify the limits of a box containing 
!    the integration region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  return
end
subroutine p23_name ( name )

!*****************************************************************************80
!
!! P23_NAME returns the name of problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'SimplexMonomial'

  return
end
subroutine p23_r8 ( action, name, value )

!*****************************************************************************80
!
!! P23_R8 sets or gets R8 parameters for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' sets the parameter to its default value;
!    'G' gets a parameter.
!    'R' sets the parameter to a random value.
!    'S' sets a parameter,
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the coefficient.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G' or 'R', then VALUE is an output quantity, 
!      and is the current value of NAME.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: c = 1.0D+00
  character ( len = * ) name
  real ( kind = 8 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c = 1.0D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value = c
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P23_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c )
      value = c
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P23_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P23_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P23_R8 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p23_region ( region )

!*****************************************************************************80
!
!! P23_REGION returns the name of the integration region for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'simplex'

  return
end
subroutine p23_title ( )

!*****************************************************************************80
!
!! P23_TITLE prints a title for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 23'
  write ( *, '(a)' ) '  Name:       SimplexMonomial'
  write ( *, '(a)' ) '  Region:     Interior of unit simplex'
  write ( *, '(a)' ) '  Integrand:  F(X) = C * product ( X(i)^E(i) )'
  write ( *, '(a)' ) '  Parameters:' 
  write ( *, '(a)' ) '              C, defaults to 1.0.'
  write ( *, '(a)' ) '              E(1:DIM_NUM) defaults to 2.'

  return
end
subroutine p24_default ( dim_num )

!*****************************************************************************80
!
!! P24_DEFAULT sets default values for problem 24.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p24_i4 ( 'D', '*', i4 )
  call p24_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p24_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P24_EXACT returns the exact integral for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) exact

  exact = 1.0D+00

  return
end
subroutine p24_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P24_F evaluates the integrand for problem 24.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!    Note that as the dimension increases, the product can grow very large.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    C(1:DIM_NUM) defaults to 0.0.
!
!    A typical, more difficult problem, has
!      C(I) = I^(1/3)
! 
!    Call P24_R8VEC to get or set C.
!
!  Integrand:
!
!    F(X) = product (   ( abs ( 4 * X(1:DIM_NUM) - 2 ) + C(1:DIM_NUM) ) 
!                     / ( 1 + C(1:DIM_NUM) ) 
!                   )
!
!  Exact Integral:
!
!    1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Joe, Frances Kuo,
!    Remark on Algorithm 659:
!    Implementing Sobol's Quasirandom Seqence Generator,
!    ACM Transactions on Mathematical Software,
!    Volume 29, Number 1, March 2003, pages 49-57.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p24_r8vec ( 'G', 'C', dim_num, c )

  do point = 1, point_num

    value(point) = product &
    ( &
        ( abs ( 4.0D+00 * x(1:dim_num,point) - 2.0D+00 ) + c(1:dim_num) ) &
      / ( 1.0D+00 + c(1:dim_num) ) &
    )
  end do

  call p24_i4 ( 'I', '#', point_num )

  return
end
subroutine p24_i4 ( action, name, value )

!*****************************************************************************80
!
!! P24_I4 sets or gets I4 parameters for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P24_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P24_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P24_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P24_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p24_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P24_LIM returns the integration limits for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) = 0.0D+00
  b(1:dim_num) = 1.0D+00

  return
end
subroutine p24_name ( name )

!*****************************************************************************80
!
!! P24_NAME returns the name of problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = '(|4X-2|+c)/(1+c)'

  return
end
subroutine p24_region ( region )

!*****************************************************************************80
!
!! P24_REGION returns the name of the integration region for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p24_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P24_R8VEC sets or gets R8VEC parameters for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the coefficient vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: c
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( c ) ) then
      deallocate ( c )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( c(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c(1:dim_num) = 0.0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P24_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c(1:dim_num) )
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P24_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P24_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P24_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p24_title ( )

!*****************************************************************************80
!
!! P24_TITLE prints a title for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 24'
  write ( *, '(a)' ) '  Name:       (|4X-2|+C)/(1+C)'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  ' // &
    'F(X) = product ( ( |4*X(i)-2| + C(i) ) / ( 1 + C(i) ) )'
  write ( *, '(a)' ) '  Parameters:' 
  write ( *, '(a)' ) '              C(1:DIM_NUM) defaults to 0.0'

  return
end
subroutine p25_default ( dim_num )

!*****************************************************************************80
!
!! P25_DEFAULT sets default values for problem 25.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8

  call p25_i4 ( 'D', '*', i4 )
  call p25_r8 ( 'D', '*', r8 )

  return
end
subroutine p25_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P25_EXACT returns the exact integral for problem 25.
!
!  Discussion:
!
!    The formula in the reference seems to yield a result
!    that is too small by 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) roundoff
  real ( kind = 8 ) term
  real ( kind = 8 ) term2

  call p25_r8 ( 'G', 'C', c )

  roundoff = epsilon ( roundoff )

  exact = 1.0D+00

  term = 1.0D+00
  i = 0

  do 

    i = i + 1
    term = term * c / real ( i, kind = 8 )

    term2 = term / real ( ( i + 1 )**dim_num, kind = 8 )

    if ( abs ( term2 ) <= roundoff * ( 1.0D+00 + abs ( exact ) ) ) then
      exit
    end if

    exact = exact + term2

  end do

  return
end
subroutine p25_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P25_F evaluates the integrand for problem 25.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integrand:
!
!    exp ( C * product ( X(1:N) ) )
!
!  Parameters:
!
!    C defaults to 0.3, and can be changed by calling P25_R8.
!
!  Exact Integral:
!
!    sum ( 1 <= i <= +oo ) C^i / ( ( i + 1 )^N * i! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Patterson,
!    [Integral #3],
!    On the Construction of a Practical Ermakov-Zolotukhin 
!    Multiple Integrator,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D. Reidel, 1987, pages 269-290,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p25_r8 ( 'G', 'C', c )

  do point = 1, point_num
    value(point) = exp ( c * product ( x(1:dim_num,point) ) )
  end do

  call p25_i4 ( 'I', '#', point_num )

  return
end
subroutine p25_i4 ( action, name, value )

!*****************************************************************************80
!
!! P25_I4 sets or gets I4 parameters for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P25_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P25_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P25_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P25_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p25_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P25_LIM returns the integration limits for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p25_name ( name )

!*****************************************************************************80
!
!! P25_NAME returns the name of problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Patterson #3, exp(c*X)'

  return
end
subroutine p25_r8 ( action, name, value )

!*****************************************************************************80
!
!! P25_R8 sets or gets R8 parameters for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' sets the parameter to its default value;
!    'G' gets a parameter.
!    'R' sets the parameter to a random value.
!    'S' sets a parameter,
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    'C' is the coefficient.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G' or 'R', then VALUE is an output quantity, 
!      and is the current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  real ( kind = 8 ), save :: c = 0.3D+00
  real ( kind = 8 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c = 0.3D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value = c
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P25_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c )
      value = c
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P25_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P25_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P25_R8 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p25_region ( region )

!*****************************************************************************80
!
!! P25_REGION returns the name of the integration region for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p25_title ( )

!*****************************************************************************80
!
!! P25_TITLE prints a title for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 25'
  write ( *, '(a)' ) '  Name:       Patterson #3, exp(c*X)'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = exp ( C * product ( X(i) ) )'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              C, defaults to 0.3.'

  return
end
subroutine p26_default ( dim_num )

!*****************************************************************************80
!
!! P26_DEFAULT sets default values for problem 26.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p26_i4 ( 'D', '*', i4 )
  call p26_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p26_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P26_EXACT returns the exact integral for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) c(dim_num)
  real ( kind = 8 ) exact

  call p26_r8vec ( 'G', 'C', dim_num, c )

  exact = product ( 1.0D+00 - exp ( -c(1:dim_num) ) )

  return
end
subroutine p26_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P26_F evaluates the integrand for problem 26.
!
!  Discussion:
!
!    The integrand is similar to that for the Patterson integral #7,
!    except for a normalization of the constants, and a (random) constant
!    factor in the integrand.
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    The integral depends on a parameter vector C(1:N).
!
!    The reference suggests choosing C at random in [0,1].
!    C(1:N) defaults to 1/N.
!
!    To get or set C, call P26_R8VEC.
!
!  Integrand:
!
!    product ( c(1:dim_num) * ( exp ( - c(1:dim_num) * x(1:dim_num) ) ) )
!    = product ( c(1:dim_num) ) * exp ( - sum ( c(1:dim_num) * x(1:dim_num) ) )
!
!  Exact Integral:
!
!    product ( 1 - exp ( c(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Patterson,
!    [Integral #1],
!    On the Construction of a Practical Ermakov-Zolotukhin 
!    Multiple Integrator,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D. Reidel, 1987, pages 269-290,
!    LC: QA299.3.N38
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p26_r8vec ( 'G', 'C', dim_num, c )

  do point = 1, point_num
    value(point) = product ( c(1:dim_num) ) &
      * exp ( - dot_product ( c(1:dim_num), x(1:dim_num,point) ) )
  end do

  call p26_i4 ( 'I', '#', point_num )

  return
end
subroutine p26_i4 ( action, name, value )

!*****************************************************************************80
!
!! P26_I4 sets or gets I4 parameters for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P26_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P26_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P26_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P26_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p26_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P26_LIM returns the integration limits for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p26_name ( name )

!*****************************************************************************80
!
!! P26_NAME returns the name of problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Patterson #1'

  return
end
subroutine p26_region ( region )

!*****************************************************************************80
!
!! P26_REGION returns the name of the integration region for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p26_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P26_R8VEC sets or gets R8VEC parameters for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the coefficient vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: c
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( c ) ) then
      deallocate ( c )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( c(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c(1:dim_num) = 1.0D+00 / real ( dim_num, kind = 8 )
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P26_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c(1:dim_num) )
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P26_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P26_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P26_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p26_title ( )

!*****************************************************************************80
!
!! P26_TITLE prints a title for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 26'
  write ( *, '(a)' ) '  Name:       Patterson #1'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) &
    '  Integrand:  F(X) = product ( C(i) * exp ( - C(i) * X(i) ) )'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              C(1:DIM_NUM) defaults to 1/DIM_NUM.'

  return
end
subroutine p27_default ( dim_num )

!*****************************************************************************80
!
!! P27_DEFAULT sets default values for problem 27.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8
  real ( kind = 8 ) r8vec(1)

  call p27_i4 ( 'D', '*', i4 )
  call p27_r8 ( 'D', '*', r8 )
  call p27_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p27_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P27_EXACT returns the exact integral for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) c(dim_num)
  real ( kind = 8 ) exact
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  call p27_r8 ( 'G', 'R', r )
  call p27_r8vec ( 'G', 'C', dim_num, c )

  exact = 2.0D+00**dim_num * &
    cos ( 2.0D+00 * pi * r + 0.5D+00 * sum ( c(1:dim_num) ) ) * &
    product ( sin ( 0.5D+00 * c(1:dim_num) ) / c(1:dim_num ) )

  return
end
subroutine p27_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P27_F evaluates the integrand for problem 27.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    The integral depends on a parameter R and vector C(1:N).
!
!    R defaults to 0.3.
!
!    The reference suggests choosing C at random in [0,1]
!    and then multiplying by the normalizing factor (25/N).
!    C(1:N) defaults to 1/N.
!
!    To get or set R, call P27_R8.
!    To get or set C, call P27_R8VEC.
!
!  Integrand:
!
!    cos ( 2 * pi * R + sum ( c(1:n) * x(1:n) ) )
!
!  Exact Integral:
!
!    2^N * cos ( 2 * pi * R + 0.5 * sum ( c(1:n) ) )
!      * product ( sin ( 0.5 * c(1:n) ) / c(1:n) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz,
!    [Integral #1]
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!    Thomas Patterson,
!    [Integral #5],
!    On the Construction of a Practical Ermakov-Zolotukhin 
!    Multiple Integrator,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D. Reidel, 1987, pages 269-290,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) arg
  real ( kind = 8 ) c(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) point
  real ( kind = 8 ) r
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p27_r8 ( 'G', 'R', r )
  call p27_r8vec ( 'G', 'C', dim_num, c )

  do point = 1, point_num
    arg = 2.0D+00 * pi * r + dot_product ( c(1:dim_num), x(1:dim_num,point) )
    value(point) = cos ( arg )
  end do

  call p27_i4 ( 'I', '#', point_num )

  return
end
subroutine p27_i4 ( action, name, value )

!*****************************************************************************80
!
!! P27_I4 sets or gets I4 parameters for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P27_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P27_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P27_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P27_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p27_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P27_LIM returns the integration limits for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p27_name ( name )

!*****************************************************************************80
!
!! P27_NAME returns the name of problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Genz #1 / Patterson #5, Oscillatory'

  return
end
subroutine p27_r8 ( action, name, value )

!*****************************************************************************80
!
!! P27_R8 sets or gets R8 parameters for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' sets the parameter to its default value;
!    'G' gets a parameter.
!    'R' sets the parameter to a random value.
!    'S' sets a parameter,
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    'R' is a factor in the integrand.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G' or 'R', then VALUE is an output quantity, 
!      and is the current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  real ( kind = 8 ), save :: r = 0.30D+00
  real ( kind = 8 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'R' .or. name == 'r' .or. name == '*' ) then
      r = 0.3D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'R' .or. name == 'r' ) then
      value = r
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P27_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'R' .or. name == 'r' ) then
      call random_number ( harvest = r )
      value = r
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P27_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'R' .or. name == 'r' ) then
      r = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P27_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P27_R8 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p27_region ( region )

!*****************************************************************************80
!
!! P27_REGION returns the name of the integration region for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p27_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P27_R8VEC sets or gets R8VEC parameters for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the coefficient vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: c
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( c ) ) then
      deallocate ( c )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( c(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c(1:dim_num) = 1.0D+00 / real ( dim_num, kind = 8 )
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then
      
    if ( name == 'C' .or. name == 'c' ) then
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P27_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c(1:dim_num) )
      c(1:dim_num) = c(1:dim_num) * 25.0D+00 / real ( dim_num, kind = 8 )
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P27_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P27_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P27_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p27_title ( )

!*****************************************************************************80
!
!! P27_TITLE prints a title for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 27'
  write ( *, '(a)' ) '  Name:       Genz #1 / Patterson #5, Oscillatory'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = cos ( 2 * pi * R ' // &
                        '+ sum ( C(i) * X(i) ) )'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              R, defaults to 0.3'
  write ( *, '(a)' ) '              C(1:DIM_NUM) defaults to 1/DIM_NUM'

  return
end
subroutine p28_default ( dim_num )

!*****************************************************************************80
!
!! P28_DEFAULT sets default values for problem 28.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p28_i4 ( 'D', '*', i4 )
  call p28_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p28_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P28_EXACT returns the exact integral for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) c(dim_num)
  real ( kind = 8 ) exact
  real ( kind = 8 ) z(dim_num)

  call p28_r8vec ( 'G', 'C', dim_num, c )
  call p28_r8vec ( 'G', 'Z', dim_num, z )

  exact = product ( &
                    (   atan ( ( 1.0D+00 - z(1:dim_num) ) / c(1:dim_num) ) &
                      + atan (             z(1:dim_num)   / c(1:dim_num) ) &
                    ) &
                    / c(1:dim_num) &
                  )

  return
end
subroutine p28_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P28_F evaluates the integrand for problem 28.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    The integral depends on vectors C(1:N) and Z(1:N).
!    To get or set C or Z, call P28_R8VEC.
!
!    The reference suggests choosing C by initializing
!    it to random values in [0,1], and then normalizing so that
!
!      sum ( 1/C(1:N)^2 ) = 170 / N^(7/2)
!
!    C(1:N) used to default to N^(9/4) / sqrt(170)
!    but this is INSUPPORTABLE for large dimension N.
!
!    So now we're setting C(1:N) to default to 1.0
!
!    The reference suggests choosing Z at random in [0,1].
!
!    Z(1:N) defaults to 0.5.
!
!  Integrand:
!
!    1 / product ( C(1:DIM_NUM)^2 + ( X(1:DIM_NUM) - Z(1:DIM_NUM) )^2 )
!
!  Exact Integral:
!
!    product ( (   arctan ( ( 1 - Z(1:DIM_NUM) ) / C(1:DIM_NUM) )
!                + arctan (       Z(1:DIM_NUM)   / C(1:DIM_NUM) ) 
!              ) / C(1:DIM_NUM)
!            )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz,
!    [Integral #2]
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!    Thomas Patterson,
!    [Integral #6],
!    On the Construction of a Practical Ermakov-Zolotukhin 
!    Multiple Integrator,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D. Reidel, 1987, pages 269-290,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ) z(dim_num)

  call p28_r8vec ( 'G', 'C', dim_num, c )
  call p28_r8vec ( 'G', 'Z', dim_num, z )

  do point = 1, point_num
    value(point) = 1.0D+00 &
      / product ( c(1:dim_num)**2 + ( x(1:dim_num,point) - z(1:dim_num) )**2 )
  end do

  call p28_i4 ( 'I', '#', point_num )

  return
end
subroutine p28_i4 ( action, name, value )

!*****************************************************************************80
!
!! P28_I4 sets or gets I4 parameters for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P28_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P28_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P28_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P28_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p28_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P28_LIM returns the integration limits for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p28_name ( name )

!*****************************************************************************80
!
!! P28_NAME returns the name of problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Genz #2 / Patterson #6, Product Peak'

  return
end
subroutine p28_region ( region )

!*****************************************************************************80
!
!! P28_REGION returns the name of the integration region for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p28_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P28_R8VEC sets or gets R8VEC parameters for problem 28.
!
!  Discussion:
!
!    For DIM_NUM = 16, computing sqrt(sqrt(dim_num^9)) was overflowing, so
!    we replaced it with ( real ( dim_num ) )^2.25
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the first vector.
!    'Z' is the base vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: c
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) s
  real ( kind = 8 ) value(dim_num)
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: z

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( c ) ) then
      deallocate ( c )
    end if
    if ( allocated ( z ) ) then
      deallocate ( z )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( c(1:dim_num) )
    allocate ( z(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c(1:dim_num) = 1.0D+00
    end if

    if ( name == 'Z' .or. name == 'z' .or. name == '*' ) then
      z(1:dim_num) = 0.5D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value(1:dim_num) = c(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P28_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c(1:dim_num) )
      s = sqrt ( &
        sqrt ( real ( dim_num**7, kind = 8 ) ) &
        * sum ( 1.0D+00 / c(1:dim_num)**2 ) / 170.0D+00 )
      c(1:dim_num) = s * c(1:dim_num)
      value(1:dim_num) = c(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      call random_number ( harvest = z(1:dim_num) )
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P28_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c(1:dim_num) = value(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      z(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P28_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P28_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p28_title ( )

!*****************************************************************************80
!
!! P28_TITLE prints a title for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 28'
  write ( *, '(a)' ) '  Name:       Genz #2 / Patterson #6, Product Peak'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  ' // &
    'F(X) = 1 / product ( C(i)^2 + ( X(i) - Z(i) )^2 )'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) &
    '              C(1:DIM_NUM) defaults to DIM_NUM^(9/4)/sqrt(170)'
  write ( *, '(a)' ) '              Z(1:DIM_NUM) defaults to 0.5.'

  return
end
subroutine p29_default ( dim_num )

!*****************************************************************************80
!
!! P29_DEFAULT sets default values for problem 29.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8
  real ( kind = 8 ) r8vec(1)

  call p29_i4 ( 'D', '*', i4 )
  call p29_r8 ( 'D', '*', r8 )
  call p29_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p29_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P29_EXACT returns the exact integral for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a
  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) exact
  integer ( kind = 4 ) ivec(dim_num)
  real ( kind = 8 ) r
  integer ( kind = 4 ) rank
  real ( kind = 8 ) total

  call p29_r8 ( 'G', 'R', r )
  call p29_r8vec ( 'G', 'C', dim_num, c )
!
!  Here, we need to generate all possible DIM_NUM tuples with
!  values of 0 or 1.
! 
  total = 0.0D+00
  rank = 0

  do

    call tuple_next ( 0, 1, dim_num, rank, ivec )

    if ( rank == 0 ) then
      exit
    end if

    total = total + (-1.0D+00)**sum(ivec(1:dim_num)) &
      / ( 1.0D+00 &
      + dot_product ( real ( ivec(1:dim_num), kind = 8 ), c(1:dim_num) ) )**r

  end do

  a = 1.0D+00
  do dim = 0, dim_num-1
    a = a * ( r + real ( dim, kind = 8 ) )
  end do

  exact = total / ( a * product ( c(1:dim_num) ) )

  return
end
subroutine p29_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P29_F evaluates the integrand for problem 29.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral parameters:
!
!    The integral depends on a parameter R,
!    and a vector C(1:N).
!
!    The Genz reference uses R = 1.  The Patterson reference
!    suggests using R = 0.3.  Note that R should NOT equal
!    any integer value between 1-DIM_NUM and 0, otherwise the
!    formula breaks down.  R should normally be strictly positive.
!
!    The Patterson reference suggests choosing C at random in [0,1]
!    and then multiplying by the normalizing factor (80/N^2).
!    This is what you will get if you "RANDOMIZE" C.
!
!    C defaults to 1/DIM_NUM.
!
!    To get or set R, call P29_R8.
!    To get or set C, call P29_R8VEC.
!
!  Integrand:
!
!    1 / ( 1 + sum ( c(1:dim_num) * x(1:dim_num) ) )^(r+dim_num)
!
!  Exact Integral:
!
!    (1/A) * ( 1 / product ( c(1:dim_num) ) ) *
!    sum(0<=I(1)<=1) sum (0<=I(2)<=1) ... sum(0<=I(dim_num)<=1)
!    (-1)^sum(I(1:dim_num)) / ( 1 + sum ( i(1:dim_num)*c(1:dim_num) ) )^r
!
!    with A = r * ( r + 1 ) * ... * ( r + dim_num - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz,
!    [Integral #3]
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!    Thomas Patterson,
!    [Integral #8],
!    On the Construction of a Practical Ermakov-Zolotukhin 
!    Multiple Integrator,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D. Reidel, 1987, pages 269-290,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) r
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  call p29_r8 ( 'G', 'R', r )
  call p29_r8vec ( 'G', 'C', dim_num, c )

  do point = 1, point_num
    value(point) = 1.0D+00 / ( 1.0D+00 &
      + dot_product ( c(1:dim_num), x(1:dim_num,point) ) )**( r + dim_num )
  end do

  call p29_i4 ( 'I', '#', point_num )

  return
end
subroutine p29_i4 ( action, name, value )

!*****************************************************************************80
!
!! P29_I4 sets or gets I4 parameters for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P29_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P29_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P29_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P29_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p29_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P29_LIM returns the integration limits for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p29_name ( name )

!*****************************************************************************80
!
!! P29_NAME returns the name of problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Genz #3 / Patterson #8, Corner Peak'

  return
end
subroutine p29_r8 ( action, name, value )

!*****************************************************************************80
!
!! P29_R8 sets or gets R8 parameters for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' sets the parameter to its default value;
!    'G' gets a parameter.
!    'R' sets the parameter to a random value.
!    'S' sets a parameter,
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    'R' is a factor in the integrand.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G' or 'R', then VALUE is an output quantity, 
!      and is the current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  real ( kind = 8 ), save :: r = 0.30D+00
  real ( kind = 8 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'R' .or. name == 'r' .or. name == '*' ) then
      r = 0.3D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'R' .or. name == 'r' ) then
      value = r
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P29_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'R' .or. name == 'r' ) then
      call random_number ( harvest = r )
      value = r
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P29_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'R' .or. name == 'r' ) then
      r = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P29_R8 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P29_R8 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p29_region ( region )

!*****************************************************************************80
!
!! P29_REGION returns the name of the integration region for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p29_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P29_R8VEC sets or gets R8VEC parameters for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the coefficient vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: c
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( c ) ) then
      deallocate ( c )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( c(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c(1:dim_num) = 1.0D+00 / real ( dim_num, kind = 8 )
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P29_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then
      
    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c(1:dim_num) )
      c(1:dim_num) = c(1:dim_num) * 80.0D+00 / real ( dim_num**2, kind = 8 )
      value(1:dim_num) = c(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P29_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P29_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P29_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p29_title ( )

!*****************************************************************************80
!
!! P29_TITLE prints a title for problem 29.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 29'
  write ( *, '(a)' ) '  Name:       Genz #3 / Patterson #8, Corner Peak'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  F(X) = 1 / ( 1 + sum( C(i) * X(i) ) )^R'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              R, defaults to 0.3'
  write ( *, '(a)' ) '              C(1:DIM_NUM) defaults to 1/DIM_NUM.'

  return
end
subroutine p30_default ( dim_num )

!*****************************************************************************80
!
!! P30_DEFAULT sets default values for problem 30.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p30_i4 ( 'D', '*', i4 )
  call p30_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p30_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P30_EXACT returns the exact integral for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_error
  real ( kind = 8 ) z(dim_num)

  call p30_r8vec ( 'G', 'C', dim_num, c )
  call p30_r8vec ( 'G', 'Z', dim_num, z )

  exact = 1.0D+00

  do dim = 1, dim_num

    exact = exact * &
      sqrt ( pi ) &
      * ( r8_error ( c(dim) * ( 1.0D+00 - z(dim) ) ) &
        + r8_error ( c(dim) *             z(dim) ) ) &
      / ( 2.0D+00 * c(dim) )

  end do

  return
end
subroutine p30_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P30_F evaluates the integrand for problem 30.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    The integral depends on vectors C(1:N) and Z(1:N).
!
!    The reference suggests choosing C at random in [0,1]
!    and then multiplying by the normalizing factor sqrt(140/N^(3/2)).
!
!    C(1:N) defaults to 1/N.
!    Z(1:N) defaults to 0.5.
!
!    To get or set C or Z, call P30_R8VEC.
!
!  Integrand:
!
!    exp ( - sum ( c(1:n)^2 * ( x(1:n) - z(1:n) )^2 )
!
!  Exact Integral:
!
!    product
!    ( sqrt ( pi )
!      * (   erf ( c(1:n) * ( 1 - z(1:n) ) ) 
!          + erf ( c(1:n) *       z(1:n)   ) ) 
!      / ( 2 * c(1:n) )
!    )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz,
!    [Integral #4]
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!    Thomas Patterson,
!    [Integral #9],
!    On the Construction of a Practical Ermakov-Zolotukhin 
!    Multiple Integrator,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D. Reidel, 1987, pages 269-290,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ) z(dim_num)

  call p30_r8vec ( 'G', 'C', dim_num, c )
  call p30_r8vec ( 'G', 'Z', dim_num, z )

  do point = 1, point_num
    value(point) = exp ( - sum ( c(1:dim_num)**2 &
      * ( x(1:dim_num,point) - z(1:dim_num) )**2 ) )
  end do

  call p30_i4 ( 'I', '#', point_num )

  return
end
subroutine p30_i4 ( action, name, value )

!*****************************************************************************80
!
!! P30_I4 sets or gets I4 parameters for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P30_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P30_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P30_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P30_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p30_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P30_LIM returns the integration limits for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p30_name ( name )

!*****************************************************************************80
!
!! P30_NAME returns the name of problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Genz #4 / Patterson #9, Gaussian'

  return
end
subroutine p30_region ( region )

!*****************************************************************************80
!
!! P30_REGION returns the name of the integration region for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p30_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P30_R8VEC sets or gets R8VEC parameters for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the first vector.
!    'Z' is the base vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: c
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: z

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( c ) ) then
      deallocate ( c )
    end if
    if ( allocated ( z ) ) then
      deallocate ( z )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( c(1:dim_num) )
    allocate ( z(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c(1:dim_num) = 1.0D+00 / real ( dim_num, kind = 8 )
    end if

    if ( name == 'Z' .or. name == 'z' .or. name == '*' ) then
      z(1:dim_num) = 0.5D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value(1:dim_num) = c(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P30_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c(1:dim_num) )
      c(1:dim_num) = c(1:dim_num) * sqrt ( 140.0D+00 &
        / sqrt ( real ( dim_num**3, kind = 8 ) ) ) &
        / sqrt ( 170.0D+00 )
      value(1:dim_num) = c(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      call random_number ( harvest = z(1:dim_num) )
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P30_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c(1:dim_num) = value(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      z(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P30_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P30_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p30_title ( )

!*****************************************************************************80
!
!! P30_TITLE prints a title for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 30'
  write ( *, '(a)' ) '  Name:       Genz #4 / Patterson #9, Gaussian'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) &
    '  Integrand:  F(X) = exp ( - sum ( C(i)^2 * ( X(i) - Z(i) )^2 )'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              C(1:DIM_NUM) defaults to 1/DIM_NUM.'
  write ( *, '(a)' ) '              Z(1:DIM_NUM) defaults to 0.5.'

  return
end
subroutine p31_default ( dim_num )

!*****************************************************************************80
!
!! P31_DEFAULT sets default values for problem 31.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p31_i4 ( 'D', '*', i4 )
  call p31_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p31_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P31_EXACT returns the exact integral for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kenneth Hanson,
!    Quasi-Monte Carlo: halftoning in high dimensions?
!    in Computatinal Imaging,
!    Edited by CA Bouman and RL Stevenson,
!    Proceedings SPIE,
!    Volume 5016, 2003, pages 161-172.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) exact
  real ( kind = 8 ) z(dim_num)
!
!  Get the limits of integration.
!
  call p31_lim ( dim_num, a, b )
!
!  Get the coefficient vector C.
!
  call p31_r8vec ( 'G', 'C', dim_num, c )
!
!  Get the location of Z.
!
  call p31_r8vec ( 'G', 'Z', dim_num, z )
!
!  The value of the DIM_NUM dimensional integral is separable
!  into the product of integrals over each dimension.
!
!  Each of these 1 dimensional integrals, in turn, is
!  easily computed, depending on where Z(I) lies with
!  respect to the limits of integration A(I) and B(I).
!
  exact = 1.0D+00

  do dim = 1, dim_num
!
!  Z < A < B
!
!  | X - Z | = X - Z from A to B.
!
    if ( z(dim) < a(dim) ) then

      exact = exact * &
      ( exp ( - c(dim) * ( a(dim) - z(dim) ) ) &
      - exp ( - c(dim) * ( b(dim) - z(dim) ) ) ) / c(dim)
!
!  A < Z < B
!
!  | X - Z | = Z - X from B to Z, 
!            = X - Z from      Z to A.
!
    else if ( z(dim) < b(dim) ) then

      exact = exact * ( 2.0D+00 &
          - exp ( - c(dim) * ( z(dim) - a(dim) ) ) &
          - exp ( - c(dim) * ( b(dim) - z(dim) ) ) ) / c(dim)
!
!  A < B < Z
!
!  | X - Z | = Z - X from A to B.
!
    else

      exact = exact * &
      ( exp ( - c(dim) * ( z(dim) - b(dim) ) ) &
      - exp ( - c(dim) * ( z(dim) - a(dim) ) ) ) / c(dim)

    end if

  end do

  return
end
subroutine p31_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P31_F evaluates the integrand for problem 31.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    There is a basis point Z associated with the integrand.
!    Z(1:DIM_NUM) defaults to ( 0.5, 0.5, ..., 0.5 ).
!    The user can set, get, or randomize this value by calling
!    P31_R8VEC.
!
!    The coefficient vector C (whose entries are usually positive)
!    controls the steepness and circularity of the pseudo-Gaussian.
!    C(1:DIM_NUM) defaults to 2.0.
!    The user can set, get, or randomize this value by calling
!    P31_R8VEC.
!
!  Integrand:
!
!    exp ( - sum ( c(1:dim_num) * abs ( x(1:dim_num) - z(1:dim_num) ) ) )
!
!  Exact Integral:
!
!    The integral is separable into
!
!      Int ( A(1) <= X(1) <= B(1) ) exp ( - C(1) * abs ( X(1) - Z(1) ) ) 
!        * Int ( A(2) <= X(2) <= B(2) ) exp ( - C(2) * abs ( X(2) - Z(2) ) )
!          * ...
!
!    Hence, the exact integral is computed as the product of
!    one dimensional integrals.  Each of these is easily computed
!    once the location of Z(I) with respect to A(I) and B(I) is
!    determined.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz,
!    [Integral #5]
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D. Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!    Kenneth Hanson,
!    Quasi-Monte Carlo: halftoning in high dimensions?
!    in Computatinal Imaging,
!    Edited by CA Bouman and RL Stevenson,
!    Proceedings SPIE,
!    Volume 5016, 2003, pages 161-172.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ) z(dim_num)

  call p31_r8vec ( 'G', 'C', dim_num, c )
  call p31_r8vec ( 'G', 'Z', dim_num, z )

  do point = 1, point_num
    value(point) = exp ( - sum ( c(1:dim_num) &
      * abs ( x(1:dim_num,point) - z(1:dim_num) ) ) )
  end do

  call p31_i4 ( 'I', '#', point_num )

  return
end
subroutine p31_i4 ( action, name, value )

!*****************************************************************************80
!
!! P31_I4 sets or gets I4 parameters for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P31_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P31_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P31_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P31_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p31_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P31_LIM returns the integration limits for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p31_name ( name )

!*****************************************************************************80
!
!! P31_NAME returns the name of problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Genz #5, Continuous'

  return
end
subroutine p31_region ( region )

!*****************************************************************************80
!
!! P31_REGION returns the name of the integration region for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p31_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P31_R8VEC sets or gets R8VEC parameters for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!      be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the coefficient vector.
!    'Z' is the base vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: c
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: z

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( c ) ) then
      deallocate ( c )
    end if
    if ( allocated ( z ) ) then
      deallocate ( z )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( c(1:dim_num) )
    allocate ( z(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c(1:dim_num) = 2.0D+00
    end if

    if ( name == 'Z' .or. name == 'z' .or. name == '*' ) then
      z(1:dim_num) = 0.5D+00
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value(1:dim_num) = c(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P31_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c(1:dim_num) )
      c(1:dim_num) = 4.0D+00 * c(1:dim_num)
      value(1:dim_num) = c(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      call random_number ( harvest = z(1:dim_num) )
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P31_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c(1:dim_num) = value(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      z(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P31_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P31_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p31_title ( )

!*****************************************************************************80
!
!! P31_TITLE prints a title for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 31'
  write ( *, '(a)' ) '  Name:       Genz #5, Continuous'
  write ( *, '(a)' ) '              Nondifferentiable peak at point Z.'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) &
    '  Integrand:  F(X) = exp ( -sum ( C(i) * | X(i) - Z(i) | ))'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              C(1:DIM_NUM) defaults to 2.0;'
  write ( *, '(a)' ) '              Z(1:DIM_NUM) defaults to 0.5;'

  return
end
subroutine p32_default ( dim_num )

!*****************************************************************************80
!
!! P32_DEFAULT sets default values for problem 32.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4
  real ( kind = 8 ) r8vec(1)

  call p32_i4 ( 'D', '*', i4 )
  call p32_r8vec ( 'D', '*', dim_num, r8vec )

  return
end
subroutine p32_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P32_EXACT returns the exact integral for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) exact
  real ( kind = 8 ) z(dim_num)

  call p32_r8vec ( 'G', 'C', dim_num, c )
  call p32_r8vec ( 'G', 'Z', dim_num, z )
  call p32_lim ( dim_num, a, b )

  exact = 1.0D+00

  do dim = 1, dim_num

    if ( z(dim) <= a(dim) ) then

      exact = exact * 0.0D+00

    else if ( z(dim) <= b(dim) ) then

      if ( c(dim) == 0.0D+00 ) then
        exact = exact * ( z(dim) - a(dim) )
      else
        exact = exact &
          * ( exp ( c(dim) * z(dim) ) - exp ( c(dim) * a(dim) ) ) / c(dim)
      end if

    else

      if ( c(dim) == 0.0D+00 ) then
        exact = exact * ( b(dim) - a(dim) )
      else
        exact = exact &
          * ( exp ( c(dim) * z(dim) ) - exp ( c(dim) * a(dim) ) ) / c(dim)
      end if

    end if

  end do

  return
end
subroutine p32_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P32_F evaluates the integrand for problem 32.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    0 <= X(1:DIM_NUM) <= 1
!
!  Integral Parameters:
!
!    The integral depends on vectors C(1:N) and Z(1:N).
!
!    The reference suggests choosing C at random in [0,1]
!    and then multiplying by the normalizing factor sqrt(140/N**(3/2)).
!
!    The default value of C(1:N) is (1/2)^(1/N).
!
!    The default value of Z(1:N) is (1/2)^(1/N).
!
!  Integrand:
!
!    exp ( c(1:n)*x(1:n) ) if all x(1:n) <= z(1:n)
!    0                        otherwise
!
!  Exact Integral:
!
!    product ( g(1:n)(x,z,a,b,c) )
!
!    where g(i)(x,z,a,b,c) =
!
!      0                                         if z(i) <= a(i)
!      ( e^(c(i)*z(i) ) - e^(c(i)*a(i)) ) / c(i) if a(i) <= z(i) <= b(i)
!      ( e^(c(i)*b(i) ) - e^(c(i)*a(i)) ) / c(i) if b(i) <= z(i)
!      
!    with obvious modifications when c(i) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz,
!    [Integral #6]
!    A Package for Testing Multiple Integration Subroutines,
!    in Numerical Integration: Recent Developments, Software
!    and Applications,
!    edited by Patrick Keast and Graeme Fairweather,
!    D Reidel, 1987, pages 337-340,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) c(dim_num)
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
  real ( kind = 8 ) z(dim_num)

  call p32_r8vec ( 'G', 'C', dim_num, c )
  call p32_r8vec ( 'G', 'Z', dim_num, z )

  do point = 1, point_num

    if ( all ( x(1:dim_num,point) <= z(1:dim_num) ) ) then
      value(point) = exp ( dot_product ( c(1:dim_num), x(1:dim_num,point) ) )
    else
      value(point) = 0.0D+00
    end if
  end do

  call p32_i4 ( 'I', '#', point_num )

  return
end
subroutine p32_i4 ( action, name, value )

!*****************************************************************************80
!
!! P32_I4 sets or gets I4 parameters for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P32_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P32_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P32_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P32_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p32_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P32_LIM returns the integration limits for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) =  0.0D+00
  b(1:dim_num) =  1.0D+00

  return
end
subroutine p32_name ( name )

!*****************************************************************************80
!
!! P32_NAME returns the name of problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Genz #6, Discontinuous'

  return
end
subroutine p32_region ( region )

!*****************************************************************************80
!
!! P32_REGION returns the name of the integration region for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'box'

  return
end
subroutine p32_r8vec ( action, name, dim_num, value )

!*****************************************************************************80
!
!! P32_R8VEC sets or gets R8VEC parameters for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the action.
!    'D' sets the internal value of the object to a default value.
!    If NAME = '*', then all variables are defaulted.
!    'G' means the current value of the object should be returned.
!    'R' means randomize the object and return it.
!    'S' means the input values of the object and its dimension should
!    be stored.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'C' is the first vector.
!    'Z' is the base vector.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the object.
!
!    Input/output, real ( kind = 8 ) VALUE(DIM_NUM), the value of the object.
!
  implicit none

  integer ( kind = 4 ) dim_num

  character ( len = * ) action
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: c
  character ( len = * ) name
  integer ( kind = 4 ), save :: dim_num_save = 0
  real ( kind = 8 ) value(dim_num)
  real ( kind = 8 ), allocatable, save, dimension ( : ) :: z

  if ( dim_num_save /= dim_num ) then
    dim_num_save = 0
    if ( allocated ( c ) ) then
      deallocate ( c )
    end if
    if ( allocated ( z ) ) then
      deallocate ( z )
    end if
  end if

  if ( dim_num_save == 0 ) then
    dim_num_save = dim_num
    allocate ( c(1:dim_num) )
    allocate ( z(1:dim_num) )
  end if

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == 'C' .or. name == 'c' .or. name == '*' ) then
      c(1:dim_num) = 0.5D+00**( 1.0D+00 / real ( dim_num, kind = 8 ) )
    end if

    if ( name == 'Z' .or. name == 'z' .or. name == '*' ) then
      z(1:dim_num) = 0.5D+00**( 1.0D+00 / real ( dim_num, kind = 8 ) )
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value(1:dim_num) = c(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P32_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'R' .or. action(1:1) == 'r' ) then

    if ( name == 'C' .or. name == 'c' ) then
      call random_number ( harvest = c(1:dim_num) )
      value(1:dim_num) = c(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      call random_number ( harvest = z(1:dim_num) )
      value(1:dim_num) = z(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P32_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c(1:dim_num) = value(1:dim_num)
    else if ( name == 'Z' .or. name == 'z' ) then
      z(1:dim_num) = value(1:dim_num)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P32_R8VEC - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P32_R8VEC - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p32_title ( )

!*****************************************************************************80
!
!! P32_TITLE prints a title for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 32'
  write ( *, '(a)' ) '  Name:       Genz #6, Discontinuous'
  write ( *, '(a)' ) '  Region:     0 <= X(i) <= 1'
  write ( *, '(a)' ) '  Integrand:  ' // &
    'F(X) = exp ( C(i) * X(i) ) if X <= Z, 0 otherwise.'
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) '              C(1:DIM_NUM) defaults to 1/DIM_NUM.'
  write ( *, '(a)' ) '              Z(1:DIM_NUM) defaults to 0.5.'

  return
end
subroutine p33_default ( dim_num )

!*****************************************************************************80
!
!! P33_DEFAULT sets default values for problem 33.
!
!  Discussion:
!
!    If a problem uses vector parameters, then the spatial dimension
!    DIM_NUM is needed as input, so that the vector parameters can be
!    properly dimensioned.
!
!    Every problem keeps a count of the number of function calls.  Calling
!    this routine causes that count to be reset to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the problem.
!
  implicit none

  integer ( kind = 4 ) dim_num

  return
end
subroutine p33_exact ( dim_num, exact )

!*****************************************************************************80
!
!! P33_EXACT returns the exact integral for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_gamma

  exact = real ( dim_num, kind = 8 ) / 2.0D+00 * sqrt ( pi**dim_num ) &
    / r8_gamma ( real ( dim_num + 4 ) / 2.0D+00 )

  return
end
subroutine p33_f ( dim_num, point_num, x, value )

!*****************************************************************************80
!
!! P33_F evaluates the integrand for problem 33.
!
!  Discussion:
!
!    The spatial dimension DIM_NUM is arbitrary.
!
!  Region:
!
!    The interior of a sphere of radius 1 centered at the origin.
!
!  Integrand:
!
!    F(X) = sum ( 1 <= I <= N ) ( X(I)^2 )
!
!  Exact Integral:
!
!    N/2 * Pi^(N/2) / Gamma ( ( N + 4 ) / 2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the function values.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  do point = 1, point_num
    value(point) = sum ( x(1:dim_num,point)**2 )
  end do

  return
end
subroutine p33_i4 ( action, name, value )

!*****************************************************************************80
!
!! P33_I4 sets or gets I4 parameters for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'D' to set a parameter to its default value.
!    'S' to set a parameter,
!    'G' to get a parameter,
!    'I' to increment a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    '#' is the number of calls to the integrand routine.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    * If ACTION = 'I', then VALUE is an input quantity, and is the
!      new value to be added to NAME.
!    * If ACTION = 'S', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'G', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  character ( len = * ) name
  integer ( kind = 4 ), save :: calls = 0
  integer ( kind = 4 ) value

  if ( action(1:1) == 'D' .or. action(1:1) == 'd' ) then

    if ( name == '#' .or. name == '*' ) then
      calls = 0
    end if

  else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == '#' ) then
      value = calls
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P33_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

    if ( name == '#' ) then
      calls = calls + value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P33_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == '#' ) then
      calls = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P33_I4 - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P33_I4 - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p33_lim ( dim_num, a, b )

!*****************************************************************************80
!
!! P33_LIM returns the integration limits for problem 33.
!
!  Discussion:
!
!    Because the integration region is the interior of a sphere
!    of radius 1 centered at the origin, the integration limits simply 
!    specify the limits of a box containing the integration region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the argument.
!
!    Output, real ( kind = 8 ) A(DIM_NUM), B(DIM_NUM), the lower and upper
!    limits of integration.
!    Note that if A = -HUGE(A), the lower limit is
!    actually negative infinity, and if B = HUGE(B), the upper limit
!    is actually infinity.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) b(dim_num)

  a(1:dim_num) = -1.0D+00
  b(1:dim_num) = +1.0D+00

  return
end
subroutine p33_name ( name )

!*****************************************************************************80
!
!! P33_NAME returns the name of problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME, the name of the problem.
!
  implicit none

  character ( len = * ) name

  name = 'Ball R^2'

  return
end
subroutine p33_region ( region )

!*****************************************************************************80
!
!! P33_REGION returns the name of the integration region for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) REGION, the name of the integration region.
!
  implicit none

  character ( len = * ) region

  region = 'ball'

  return
end
subroutine p33_title ( )

!*****************************************************************************80
!
!! P33_TITLE prints a title for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Problem 33'
  write ( *, '(a)' ) '  Name:       Ball R^2'
  write ( *, '(a)' ) '  Region:     Sphere interior, radius 1, center 0'
  write ( *, '(a)' ) '  Integrand:  F(X) = sum ( X(1:N)^2 )'

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the combinatorial coefficient C(N,K).
!
!  Discussion:
!
!    C(N,K) is the number of distinct combinations of K objects
!    chosen from a set of N distinct objects.  A combination is
!    like a set, in that order does not matter.
!
!    The formula is:
!
!      C(N,K) = N! / ( (N-K)! * K! )
!
!    Real arithmetic is used, and C(N,K) is computed directly, via
!    Gamma functions, rather than recursively.
!
!    For example, the number of combinations of 2 things chosen from 
!    5 is 10.  Our formula is
!
!      C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
!
!    The actual combinations may be represented as:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3),
!      (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value of N.
!
!    Input, integer ( kind = 4 ) K, the value of K.
!
!    Output, real ( kind = 8 ) R8_CHOOSE, the value of C(N,K)
!
  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) fack
  real ( kind = 8 ) facn
  real ( kind = 8 ) facnmk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) value

  if ( n < 0 ) then

    value = 0.0D+00

  else if ( k == 0 ) then

    value = 1.0D+00

  else if ( k == 1 ) then

    value = real ( n, kind = 8 )

  else if ( 1 < k .and. k < n - 1 ) then

    arg = real ( n + 1, kind = 8 )
    facn = r8_gamma_log ( arg )

    arg = real ( k + 1, kind = 8 )
    fack = r8_gamma_log ( arg )

    arg = real ( n - k + 1, kind = 8 )
    facnmk = r8_gamma_log ( arg )

    value = real ( nint ( exp ( facn - fack - facnmk ) ), kind = 8 )

  else if ( k == n-1 ) then

    value = real ( n, kind = 8 )

  else if ( k == n ) then

    value = 1.0D+00

  else

    value = 0.0D+00

  end if

  r8_choose = value

  return
end
function r8_error ( x )

!*****************************************************************************80
!
!! R8_ERROR computes the error function.
!
!  Discussion:
!
!    This function was renamed "R8_ERROR" from "ERF", to avoid a conflict
!    with the name of a corresponding routine often, but not always,
!    supplied as part of the math support library.
!
!    The definition of the error function is:
!
!      ERF(X) = ( 2 / SQRT ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( -T**2 ) dT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) R8_ERROR, the value of the error function at X.
!
  implicit none

  real ( kind = 8 ) csevl
  real ( kind = 8 ) r8_error
  real ( kind = 8 ) r8_errorc
  real ( kind = 8 ), parameter, dimension ( 13 ) :: erfcs = (/ &
    -0.049046121234691808D+00, -0.14226120510371364D+00, &
     0.010035582187599796D+00, -0.000576876469976748D+00, &
     0.000027419931252196D+00, -0.000001104317550734D+00, &
     0.000000038488755420D+00, -0.000000001180858253D+00, &
     0.000000000032334215D+00, -0.000000000000799101D+00, &
     0.000000000000017990D+00, -0.000000000000000371D+00, &
     0.000000000000000007D+00 /)
  integer ( kind = 4 ) inits
  integer ( kind = 4 ), save :: nterf = 0
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ), parameter :: sqrtpi = 1.7724538509055160D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xbig = 0.0D+00
  real ( kind = 8 ) y
!
!  Initialize the Chebyshev series.
!
  if ( nterf == 0 ) then
    nterf = inits ( erfcs, 13, 0.1D+00 * epsilon ( erfcs ) )
    xbig = sqrt ( - log ( sqrtpi * epsilon ( xbig ) ) )
    sqeps = sqrt ( 2.0D+00 * epsilon ( sqeps ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    r8_error = 2.0D+00 * x / sqrtpi
  else if ( y <= 1.0D+00 ) then
    r8_error = x * ( 1.0D+00 + csevl ( 2.0D+00 * x**2 - 1.0D+00, erfcs, nterf ) )
  else if ( y <= xbig ) then
    r8_error = sign ( 1.0D+00 - r8_errorc ( y ), x )
  else
    r8_error = sign ( 1.0D+00, x )
  end if

  return
end
function r8_errorc ( x )

!*****************************************************************************80
!
!! R8_ERRORC computes the complementary error function.
!
!  Discussion:
!
!    This function was renamed "R8_ERRORC" from "ERFC", to avoid a conflict
!    with the name of a corresponding routine often, but not always,
!    supplied as part of the math support library.
!
!    The definition of the complementary error function is:
!
!      ERFC(X) = 1 - ERF(X)
!
!    where ERF(X) is the error function.
!
!  Modified:
!
!    26 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) R8_ERRORC, the value of the complementary 
!    error function at X.
!
  implicit none

  real ( kind = 8 ) csevl
  real ( kind = 8 ), parameter, dimension ( 13 ) :: erfcs = (/ &
   -0.049046121234691808D+00, &
   -0.14226120510371364D+00,  &
    0.010035582187599796D+00, &
   -0.000576876469976748D+00, &
    0.000027419931252196D+00, &
   -0.000001104317550734D+00, &
    0.000000038488755420D+00, &
   -0.000000001180858253D+00, &
    0.000000000032334215D+00, &
   -0.000000000000799101D+00, &
    0.000000000000017990D+00, &
   -0.000000000000000371D+00, &
    0.000000000000000007D+00 /)
  real ( kind = 8 ), parameter, dimension ( 24 ) :: erfccs = (/ &
    0.0715179310202925D+00, &
   -0.026532434337606719D+00, &
    0.001711153977920853D+00, &
   -0.000163751663458512D+00, &
    0.000019871293500549D+00, &
   -0.000002843712412769D+00, &
    0.000000460616130901D+00, &
   -0.000000082277530261D+00, &
    0.000000015921418724D+00, &
   -0.000000003295071356D+00, &
    0.000000000722343973D+00, &
   -0.000000000166485584D+00, &
    0.000000000040103931D+00, &
   -0.000000000010048164D+00, &
    0.000000000002608272D+00, &
   -0.000000000000699105D+00, &
    0.000000000000192946D+00, &
   -0.000000000000054704D+00, &
    0.000000000000015901D+00, &
   -0.000000000000004729D+00, &
    0.000000000000001432D+00, &
   -0.000000000000000439D+00, &
    0.000000000000000138D+00, &
   -0.000000000000000048D+00 /)
  real ( kind = 8 ), parameter, dimension ( 23 ) :: erc2cs = (/ &
   -0.069601346602309501D+00, &
   -0.041101339362620893D+00, &
    0.003914495866689626D+00, &
   -0.000490639565054897D+00, &
    0.000071574790013770D+00, &
   -0.000011530716341312D+00, &
    0.000001994670590201D+00, &
   -0.000000364266647159D+00, &
    0.000000069443726100D+00, &
   -0.000000013712209021D+00, &
    0.000000002788389661D+00, &
   -0.000000000581416472D+00, &
    0.000000000123892049D+00, &
   -0.000000000026906391D+00, &
    0.000000000005942614D+00, &
   -0.000000000001332386D+00, &
    0.000000000000302804D+00, &
   -0.000000000000069666D+00, &
    0.000000000000016208D+00, &
   -0.000000000000003809D+00, &
    0.000000000000000904D+00, &
   -0.000000000000000216D+00, &
    0.000000000000000052D+00 /)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) inits
  integer ( kind = 4 ), save :: nterc2 = 0
  integer ( kind = 4 ), save :: nterf = 0
  integer ( kind = 4 ), save :: nterfc = 0
  real ( kind = 8 ) r8_errorc
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ), parameter :: sqrtpi = 1.7724538509055160D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xmax = 0.0D+00
  real ( kind = 8 ), save :: xsml = 0.0D+00
  real ( kind = 8 ) y

  if ( nterf == 0 ) then

    eta = 0.1D+00 * epsilon ( eta )
    nterf = inits ( erfcs, 13, eta )
    nterfc = inits ( erfccs, 24, eta )
    nterc2 = inits ( erc2cs, 23, eta )

    xsml = -sqrt ( - log ( sqrtpi * epsilon ( xsml ) ) )
    xmax = sqrt ( - log ( sqrtpi * tiny ( xmax ) ) )
    xmax = xmax - 0.5D+00 * log ( xmax ) / xmax - 0.01D+00
    sqeps = sqrt ( 2.0D+00 * epsilon ( sqeps ) )

  end if

  if ( x <= xsml ) then
    r8_errorc = 2.0D+00
    return
  end if
!
!  X so big that ERFC will underflow.
!
  if ( xmax < x ) then
    r8_errorc = 0.0D+00
    return
  end if

  y = abs ( x )
!
!  erfc(x) = 1.0D+00 - erf(x) for -1 <= x <= 1.
!
  if ( y <= 1.0D+00 ) then

    if ( y < sqeps ) then
      r8_errorc = 1.0D+00 - 2.0D+00 * x / sqrtpi
    else if ( sqeps <= y ) then
      r8_errorc = 1.0D+00 - x * ( 1.0D+00 + &
        csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf ) )
    end if

    return

  end if
!
!  For 1 < |x| <= xmax, erfc(x) = 1.0D+00 - erf(x)
!
  y = y * y

  if ( y <= 4.0D+00 ) then
    r8_errorc = exp ( -y ) / abs ( x ) * ( 0.5D+00 &
      + csevl ( ( 8.0D+00 / y - 5.0D+00 ) / 3.0D+00, erc2cs, nterc2 ) )
  else
    r8_errorc = exp ( -y ) / abs ( x ) * ( 0.5D+00 &
      + csevl ( 8.0D+00 / y - 1.0D+00, erfccs, nterfc ) )
  end if

  if ( x < 0.0D+00 ) then
    r8_errorc = 2.0D+00 - r8_errorc
  end if

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Discussion:
!
!    The program uses rational functions that theoretically approximate
!    log ( GAMMA(X) ) to at least 18 significant decimal digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Mathematics of Computation,
!    Volume 21, 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    Hart, Ward Cheney, Charles Lawson, Maehly, Charles Mesztenyi, 
!    John Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function. 
!    X must be positive.
!
!    Output, real ( kind = 8 ) R8_GAMMA_LOG, the logarithm of the Gamma 
!    function of X.  If X <= 0.0, or if overflow would occur, the
!    program returns the value HUGE().
!
!  Machine-dependent constants:
!
!    BETA   - radix for the floating-point representation.
!
!    MAXEXP - the smallest positive power of BETA that overflows.
!
!    XBIG   - largest argument for which LN(GAMMA(X)) is representable
!             in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!    XINF   - largest machine representable floating-point number;
!             approximately BETA**MAXEXP.
!
!    FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!    Approximate values for some important machines are:
!
!                              BETA      MAXEXP         XBIG
!
!    CRAY-1        (S.P.)        2        8191       9.62D+2461
!    Cyber 180/855
!      under NOS   (S.P.)        2        1070       1.72D+319
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)        2         128       4.08D+36
!    IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!    IBM 3033      (D.P.)       16          63       4.29D+73
!    VAX D-Format  (D.P.)        2         127       2.05D+36
!    VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                            FRTBIG
!
!    CRAY-1        (S.P.)   3.13D+615
!    Cyber 180/855
!      under NOS   (S.P.)   6.44D+79
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)   1.42D+9
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)   2.25D+76
!    IBM 3033      (D.P.)   2.56D+18
!    VAX D-Format  (D.P.)   1.20D+9
!    VAX G-Format  (D.P.)   1.89D+76
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ), parameter :: d1 = - 5.772156649015328605195174D-01
  real ( kind = 8 ), parameter :: d2 =   4.227843350984671393993777D-01
  real ( kind = 8 ), parameter :: d4 =   1.791759469228055000094023D+00
  real ( kind = 8 ) eps
  real ( kind = 8 ), parameter :: frtbig = 1.42D+09
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 4.08D+36
  real ( kind = 8 ) xden
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0D+00 .or. xbig < x ) then
    r8_gamma_log = huge ( r8_gamma_log )
    return
  end if

  eps = epsilon ( eps )

  if ( x <= eps ) then

    res = - log ( x )

  else if ( x <= 1.5D+00 ) then

    if ( x < pnt68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0D+00
      xm1 = ( x - 0.5D+00 ) - 0.5D+00
    end if

    if ( x <= 0.5D+00 .or. pnt68 <= x ) then

      xden = 1.0D+00
      xnum = 0.0D+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5D+00 ) - 0.5D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0D+00 ) then

    xm2 = x - 2.0D+00
    xden = 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0D+00 ) then

    xm4 = x - 4.0D+00
    xden = -1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0D+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5D+00 * corr
    res = res + x * ( corr - 1.0D+00 )

  end if

  r8_gamma_log = res

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    If all the entries are integers, the data if printed
!    in integer format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
    do i = 1, n
      write ( *, '(i8,i8)' ) i, int ( a(i) )
    end do
  else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
    do i = 1, n
      write ( *, '(i8,f14.6)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i8,g14.6)' ) i, a(i)
    end do
  end if

  return
end
subroutine random_initialize ( seed_input )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee.
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator
!    a bit, this routine goes through the tedious process of getting the
!    size of the random number seed, making up values based on the current
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed.
!    However, if the input value is 0, the routine will come up with
!    its own "suggestion", based on the system clock.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_input
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 8 ) t
  integer ( kind = 4 ), parameter :: warm_up = 100

  seed = seed_input
!
!  Initialize the random seed routine.
!
  call random_seed ( )
!
!  Determine the size of the random number seed vector.
!
  call random_seed ( size = seed_size )
!
!  Allocate a vector of the right size to be used as a random seed.
!
  allocate ( seed_vector(seed_size) )
!
!  If the user supplied a SEED value, use that.
!
!  Otherwise, use the system clock value to make up a value that is
!  likely to change based on when this routine is called.
!
  if ( seed /= 0 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
        seed
    end if

  end if
!
!  Set the seed vector.  We don't know the significance of the
!  individual entries of the internal seed vector, so we'll just set
!  all entries to SEED.
!
  seed_vector(1:seed_size) = seed
!
!  Now call RANDOM_SEED, and tell it to use this seed vector.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times just to "warm it up".
!
  do i = 1, warm_up
    call random_number ( harvest = t )
  end do

  return
end
subroutine simplex_unit_volume_nd ( dim_num, volume )

!*****************************************************************************80
!
!! SIMPLEX_UNIT_VOLUME_ND computes the volume of the unit simplex in ND.
!
!  Discussion:
!
!    The formula is simple: volume = 1/N!.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the cone.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) volume

  volume = 1.0D+00
  do dim = 1, dim_num
    volume = volume / real ( dim, kind = 8 )
  end do

  return
end
subroutine sphere_unit_area_nd ( dim_num, area )

!*****************************************************************************80
!
!! SPHERE_UNIT_AREA_ND computes the surface area of a unit sphere in ND.
!
!  Discussion:
!
!    The unit sphere in ND satisfies:
!
!      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
!
!    Results for the first few values of N are:
!
!    DIM_NUM   Area
!
!     2    2        * PI
!     3    4        * PI
!     4  ( 2 /   1) * PI**2
!     5  ( 8 /   3) * PI**2
!     6  ( 1 /   1) * PI**3
!     7  (16 /  15) * PI**3
!     8  ( 1 /   3) * PI**4
!     9  (32 / 105) * PI**4
!    10  ( 1 /  12) * PI**5
!
!    For the unit sphere, Area(DIM_NUM) = DIM_NUM * Volume(DIM_NUM)
!
!    Sphere_Unit_Area ( DIM_NUM ) = 2 * PI**(DIM_NUM/2) / Gamma ( DIM_NUM / 2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Output, real ( kind = 8 ) AREA, the area of the sphere.
!
  implicit none

  real ( kind = 8 ) area
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  if ( mod ( dim_num, 2 ) == 0 ) then
    m = dim_num / 2
    area = 2.0D+00 * ( pi )**m
    do i = 1, m-1
      area = area / real ( i, kind = 8 )
    end do
  else
    m = ( dim_num - 1 ) / 2
    area = ( pi )**m * 2.0D+00**dim_num
    do i = m+1, 2*m
      area = area / real ( i,  kind = 8 )
    end do
  end if

  return
end
subroutine sphere_unit_volume_nd ( dim_num, volume )

!*****************************************************************************80
!
!! SPHERE_UNIT_VOLUME_ND: volume of a unit sphere in multi-dimensions.
!
!  Discussion:
!
!    DIM_NUM  Volume
!
!    2             PI
!    3  (4/3)    * PI
!    4  (1/2)    * PI^2
!    5  (8/15)   * PI^2
!    6  (1/6)    * PI^3
!    7  (16/105) * PI^3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the sphere.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) volume

  volume = 1.0D+00

  if ( mod ( dim_num, 2 ) == 0 ) then
    m = dim_num / 2
    do i = 1, m
      volume = volume * pi / real ( i, kind = 8 )
    end do
  else
    m = ( dim_num - 1 ) / 2
    do i = 1, m
      volume = volume * pi * 2.0D+00
    end do
    do i = m+1, 2*m+1
      volume = volume * 2.0D+00 / real ( i, kind = 8 )
    end do
  end if

  return
end
subroutine sphere_volume_nd ( dim_num, r, volume )

!*****************************************************************************80
!
!! SPHERE_VOLUME_ND computes the volume of a sphere in multi-dimensions.
!
!  Discussion:
!
!    A sphere in DIM_NUM dimensions satisfies the equation:
!
!      sum ( ( X(1:DIM_NUM) - XC(1:DIM_NUM) )^2 ) = R^2
!
!    where R is the radius, and XC(1:DIM_NUM) is the center.
!
!    Some sample values are
!
!      DIM_NUM  Volume
!
!      2             PI   * R^2
!      3  (4/3)    * PI   * R^3
!      4  (1/2)    * PI^2 * R^4
!      5  (8/15)   * PI^2 * R^5
!      6  (1/6)    * PI^3 * R^6
!      7  (16/105) * PI^3 * R^7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the sphere.
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) r
  real ( kind = 8 ) volume

  call sphere_unit_volume_nd ( dim_num, volume )

  volume = volume * r**dim_num

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tuple_next ( m1, m2, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT computes the next element of a tuple space.
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between M1 and M2.  The elements are produced one at a time.
!    The first element is
!      (M1,M1,...,M1),
!    the second element is
!      (M1,M1,...,M1+1),
!    and the last element is
!      (M2,M2,...,M2)
!    Intermediate elements are produced in lexicographic order.
!
!    For example, if our scalar input was
!
!      N = 2,
!      M1 = 1,
!      M2 = 3
!
!    then the sequence of input and output vectors would be:
!
!      INPUT        OUTPUT
!      -------      -------
!      Rank  X      Rank   X
!      ----  ---    -----  ---
!      0     * *    1      1 1
!      1     1 1    2      1 2
!      2     1 2    3      1 3
!      3     1 3    4      2 1
!      4     2 1    5      2 2
!      5     2 2    6      2 3
!      6     2 3    7      3 1
!      7     3 1    8      3 2
!      8     3 2    9      3 3
!      9     3 3    0      0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M1, M2, the minimum and maximum entries.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input/output, integer ( kind = 4 ) RANK, counts the elements.
!    On first call, set RANK to 0.  Thereafter, the output value of RANK
!    will indicate the order of the element returned.  When there are no
!    more elements, RANK will be returned as 0.
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( m2 < m1 ) then
    rank = 0
    return
  end if

  if ( rank <= 0 ) then

    x(1:n) = m1
    rank = 1

  else

    rank = rank + 1
    i = n

    do

      if ( x(i) < m2 ) then
        x(i) = x(i) + 1
        exit
      end if

      x(i) = m1

      if ( i == 1 ) then
        rank = 0
        x(1:n) = m1
        exit
      end if

      i = i - 1

    end do

  end if

  return
end
