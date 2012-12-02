subroutine bashforth_set ( n, x, w )

!*****************************************************************************80
!
!! BASHFORTH_SET sets abscissas and weights for Adams-Bashforth quadrature.
!
!  Discussion:
!
!    Adams-Bashforth quadrature formulas are normally used in solving
!    ordinary differential equations, and are not really suitable for
!    general quadrature computations.  However, an Adams-Bashforth formula
!    is equivalent to approximating the integral of F(Y(X)) between X(M)
!    and X(M+1), using an explicit formula that relies only on known values
!    of F(Y(X)) at X(M-N+1) through X(M).  For this reason, the formulas
!    have been included here.
!
!    Suppose the unknown function is denoted by Y(X), with derivative
!    F(Y(X)), and that approximate values of the function are known at a
!    series of X values, which we write as X(1), X(2), ..., X(M).  We write
!    the value Y(X(1)) as Y(1) and so on.
!
!    Then the solution of the ODE Y'=F(X,Y) at the next point X(M+1) is
!    computed by:
!
!      Y(M+1) = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dx
!             = Y(M) + H * Sum ( 1 <= I <= N ) W(I) * F(Y(M+1-I)) approximately.
!
!    In the documentation that follows, we replace F(Y(X)) by F(X).
!
!    The integral:
!
!      integral ( 0 <= x <= 1 ) f(x) dx.
!
!    The quadrature formula:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The Adams-Bashforth formulas require equally spaced data.
!
!    Here is how the formula is applied in the case with non-unit spacing:
!
!      Integral ( A <= X <= A+H ) F(X) dx =
!      H * Sum ( 1 <= I <= N ) W(I) * F ( A - (I-1)*H ),
!      approximately.
!
!    The reference lists the second coefficient of the order 8 Adams-Bashforth
!    formula as
!      w(2) =  -1162169.0D+00 / 120960.0D+00
!    but this should be
!      w(2) =  -1152169.0D+00 / 120960.0D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Leon Lapidus, John Seinfeld,
!    Numerical Solution of Ordinary Differential Equations,
!    Academic Press, 1971,
!    ISBN: 0124366503,
!    LC: QA3.M32.v74.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N should be between 1 and 10, 12, 14, 16, 18 or 20.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!    W(1) is the weight at X(1) = 0,
!    W(2) the weight at X(2) = -1,
!    and so on.  The weights are rational, and should sum to 1.  Some
!    weights may be negative.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    w(1) = 1.0D+00

  else if ( n == 2 ) then

    d = 2.0D+00

    w(1) =   3.0D+00 / d
    w(2) = - 1.0D+00 / d

  else if ( n == 3 ) then

    d = 12.0D+00

    w(1) =   23.0D+00 / d
    w(2) = - 16.0D+00 / d
    w(3) =    5.0D+00 / d

  else if ( n == 4 ) then

    d = 24.0D+00

    w(1) =   55.0D+00 / d
    w(2) = - 59.0D+00 / d
    w(3) =   37.0D+00 / d
    w(4) =  - 9.0D+00 / d

  else if ( n == 5 ) then

    d = 720.0D+00

    w(1) =   1901.0D+00 / d
    w(2) = - 2774.0D+00 / d
    w(3) =   2616.0D+00 / d
    w(4) = - 1274.0D+00 / d
    w(5) =    251.0D+00 / d

  else if ( n == 6 ) then

    d = 1440.0D+00

    w(1) =   4277.0D+00 / d
    w(2) = - 7923.0D+00 / d
    w(3) =   9982.0D+00 / d
    w(4) = - 7298.0D+00 / d
    w(5) =   2877.0D+00 / d
    w(6) =  - 475.0D+00 / d

  else if ( n == 7 ) then

    d = 60480.0D+00

    w(1) =    198721.0D+00 / d
    w(2) =  - 447288.0D+00 / d
    w(3) =    705549.0D+00 / d
    w(4) =  - 688256.0D+00 / d
    w(5) =    407139.0D+00 / d
    w(6) =  - 134472.0D+00 / d
    w(7) =     19087.0D+00 / d

  else if ( n == 8 ) then

    d = 120960.0D+00

    w(1) =     434241.0D+00 / d
    w(2) =  - 1152169.0D+00 / d
    w(3) =    2183877.0D+00 / d
    w(4) =  - 2664477.0D+00 / d
    w(5) =    2102243.0D+00 / d
    w(6) =  - 1041723.0D+00 / d
    w(7) =     295767.0D+00 / d
    w(8) =    - 36799.0D+00 / d

 else if ( n == 9 ) then

    d = 3628800.0D+00

    w(1) =   14097247.0D+00 / d
    w(2) =  -43125206.0D+00 / d
    w(3) =   95476786.0D+00 / d
    w(4) = -139855262.0D+00 / d
    w(5) =  137968480.0D+00 / d
    w(6) =  -91172642.0D+00 / d
    w(7) =   38833486.0D+00 / d
    w(8) =   -9664106.0D+00 / d
    w(9) =    1070017.0D+00 / d

  else if ( n == 10 ) then

    d = 7257600.0D+00

    w( 1) =   30277247.0D+00 / d
    w( 2) = -104995189.0D+00 / d
    w( 3) =  265932680.0D+00 / d
    w( 4) = -454661776.0D+00 / d
    w( 5) =  538363838.0D+00 / d
    w( 6) = -444772162.0D+00 / d
    w( 7) =  252618224.0D+00 / d
    w( 8) =  -94307320.0D+00 / d
    w( 9) =   20884811.0D+00 / d
    w(10) =   -2082753.0D+00 / d

  else if ( n == 12 ) then

    d = 958003200.0D+00

    w( 1) =    4527766399.0D+00 / d
    w( 2) =  -19433810163.0D+00 / d
    w( 3) =   61633227185.0D+00 / d
    w( 4) = -135579356757.0D+00 / d
    w( 5) =  214139355366.0D+00 / d
    w( 6) = -247741639374.0D+00 / d
    w( 7) =  211103573298.0D+00 / d
    w( 8) = -131365867290.0D+00 / d
    w( 9) =   58189107627.0D+00 / d
    w(10) =  -17410248271.0D+00 / d
    w(11) =    3158642445.0D+00 / d
    w(12) =    -262747265.0D+00 / d

  else if ( n == 14 ) then

    d = 5230697472000.0D+00

    w( 1) =    27511554976875.0D+00 / d
    w( 2) =  -140970750679621.0D+00 / d
    w( 3) =   537247052515662.0D+00 / d
    w( 4) = -1445313351681906.0D+00 / d
    w( 5) =  2854429571790805.0D+00 / d
    w( 6) = -4246767353305755.0D+00 / d
    w( 7) =  4825671323488452.0D+00 / d
    w( 8) = -4204551925534524.0D+00 / d
    w( 9) =  2793869602879077.0D+00 / d
    w(10) = -1393306307155755.0D+00 / d
    w(11) =   505586141196430.0D+00 / d
    w(12) =  -126174972681906.0D+00 / d
    w(13) =    19382853593787.0D+00 / d
    w(14) =    -1382741929621.0D+00 / d

  else if ( n == 16 ) then

    d = 62768369664000.0D+00

    w( 1) =     362555126427073.0D+00 / d
    w( 2) =   -2161567671248849.0D+00 / d
    w( 3) =    9622096909515337.0D+00 / d
    w( 4) =  -30607373860520569.0D+00 / d
    w( 5) =   72558117072259733.0D+00 / d
    w( 6) = -131963191940828581.0D+00 / d
    w( 7) =  187463140112902893.0D+00 / d
    w( 8) = -210020588912321949.0D+00 / d
    w( 9) =  186087544263596643.0D+00 / d
    w(10) = -129930094104237331.0D+00 / d
    w(11) =   70724351582843483.0D+00 / d
    w(12) =  -29417910911251819.0D+00 / d
    w(13) =    9038571752734087.0D+00 / d
    w(14) =   -1934443196892599.0D+00 / d
    w(15) =     257650275915823.0D+00 / d
    w(16) =     -16088129229375.0D+00 / d

  else if ( n == 18 ) then

    d = 64023737057280000.0D+00

    w( 1) =     401972381695456831.0D+00 / d
    w( 2) =   -2735437642844079789.0D+00 / d
    w( 3) =   13930159965811142228.0D+00 / d
    w( 4) =  -51150187791975812900.0D+00 / d
    w( 5) =  141500575026572531760.0D+00 / d
    w( 6) = -304188128232928718008.0D+00 / d
    w( 7) =  518600355541383671092.0D+00 / d
    w( 8) = -710171024091234303204.0D+00 / d
    w( 9) =  786600875277595877750.0D+00 / d
    w(10) = -706174326992944287370.0D+00 / d
    w(11) =  512538584122114046748.0D+00 / d
    w(12) = -298477260353977522892.0D+00 / d
    w(13) =  137563142659866897224.0D+00 / d
    w(14) =  -49070094880794267600.0D+00 / d
    w(15) =   13071639236569712860.0D+00 / d
    w(16) =   -2448689255584545196.0D+00 / d
    w(17) =     287848942064256339.0D+00 / d
    w(18) =     -15980174332775873.0D+00 / d

  else if ( n == 20 ) then

    d = 102181884343418880000.0D+00

    w( 1) =      691668239157222107697.0D+00 / d
    w( 2) =    -5292843584961252933125.0D+00 / d
    w( 3) =    30349492858024727686755.0D+00 / d
    w( 4) =  -126346544855927856134295.0D+00 / d
    w( 5) =   399537307669842150996468.0D+00 / d
    w( 6) =  -991168450545135070835076.0D+00 / d
    w( 7) =  1971629028083798845750380.0D+00 / d
    w( 8) = -3191065388846318679544380.0D+00 / d
    w( 9) =  4241614331208149947151790.0D+00 / d
    w(10) = -4654326468801478894406214.0D+00 / d
    w(11) =  4222756879776354065593786.0D+00 / d
    w(12) = -3161821089800186539248210.0D+00 / d
    w(13) =  1943018818982002395655620.0D+00 / d
    w(14) =  -970350191086531368649620.0D+00 / d
    w(15) =   387739787034699092364924.0D+00 / d
    w(16) =  -121059601023985433003532.0D+00 / d
    w(17) =    28462032496476316665705.0D+00 / d
    w(18) =    -4740335757093710713245.0D+00 / d
    w(19) =      498669220956647866875.0D+00 / d
    w(20) =      -24919383499187492303.0D+00 / d

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASHFORTH_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 through 10, 12, 14, 16, 18 or 20.'
    stop

  end if

  do i = 1, n
    x(i) = real ( 1 - i, kind = 8 )
  end do

  return
end
subroutine bdf_set ( n, alpha, beta, gamma )

!*****************************************************************************80
!
!! BDF_SET sets weights for ODE backward differentiation.
!
!  Discussion:
!
!    GAMMA * Y(N+1) = Sum ( 1 <= I <= N ) ALPHA(I) * Y(N+1-I)
!                     + dx * BETA * Y'(X(N+1),Y(N+1))
!
!    This is equivalent to the backward differentiation corrector formulas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N should be between 1 and 6.
!
!    Output, real ( kind = 8 ) ALPHA(N), BETA, GAMMA, the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha(n)
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma

  if ( n == 1 ) then
    beta =     1.0D+00
    gamma =    1.0D+00
    alpha(1) = 1.0D+00
  else if ( n == 2 ) then
    beta =       2.0D+00
    gamma =      3.0D+00
    alpha(1) =   4.0D+00
    alpha(2) = - 1.0D+00
  else if ( n == 3 ) then
    beta =       6.0D+00
    gamma =     11.0D+00
    alpha(1) =  18.0D+00
    alpha(2) = - 9.0D+00
    alpha(3) =   2.0D+00
  else if ( n == 4 ) then
    beta =       12.0D+00
    gamma =      25.0D+00
    alpha(1) =   48.0D+00
    alpha(2) = - 36.0D+00
    alpha(3) =   16.0D+00
    alpha(4) =  - 3.0D+00
  else if ( n == 5 ) then
    beta =        60.0D+00
    gamma =      137.0D+00
    alpha(1) =   300.0D+00
    alpha(2) = - 300.0D+00
    alpha(3) =   200.0D+00
    alpha(4) =  - 75.0D+00
    alpha(5) =    12.0D+00
  else if ( n == 6 ) then
    beta =        60.0D+00
    gamma =      147.0D+00
    alpha(1) =   360.0D+00
    alpha(2) = - 450.0D+00
    alpha(3) =   400.0D+00
    alpha(4) = - 225.0D+00
    alpha(5) =    72.0D+00
    alpha(6) =  - 10.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BDF_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal order N = ', n
    stop
  end if

  return
end
subroutine bdf_sum ( func, n, x, w, result )

!*****************************************************************************80
!
!! BDF_SUM: an explicit backward difference quadrature rule for [0,1].
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * bdf^(i-1) f ( 0 )
!
!    The integral from 0 to 1 is approximated using data at X = 0,
!    -1, -2, ..., -N+1.  This is a form of extrapolation, and
!    the approximation can become poor as N increases.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the function which evaluates
!    the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Input, real ( kind = 8 ) W(N), the weights.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) diftab(n)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) result
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  do i = 1, n
    diftab(i) = func ( x(i) )
  end do

  do i = 2, n
    diftab(n:i:-1) = ( diftab(n-1:i-1:-1) - diftab(n:i:-1) )
  end do

  result = dot_product ( w(1:n), diftab(1:n) )

  return
end
subroutine bdfc_set ( n, x, w )

!*****************************************************************************80
!
!! BDFC_SET sets weights for backward differentiation corrector quadrature.
!
!  Discussion:
!
!    A backward differentiation corrector formula is defined for a set
!    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
!    that the values of the function to be integrated are known at the
!    abscissas, the formula is written in terms of the function value at
!    X(1), and the backward differences at X(1) that approximate the
!    derivatives there.
!
!    The integral:
!
!      integral ( 0 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * bd^(i-1) f ( 1 ).
!
!    Here, "BD^(I-1) F ( 1 )" denotes the (I-1)st backward difference
!    of F at X = 1, using a spacing of 1.  In particular,
!
!    BD^0 F(1) = F(1)
!    BD^1 F(1) = F(1) - F(0)
!    BD^2 F(1) = F(1) - 2 * F(0) + F(-1 )
!
!    The relationship between a backward difference corrector and the
!    corresponding Adams-Moulton formula may be illustrated for the
!    BDF corrector of order 4:
!
!      BD^0 F(1) - 1/2 * BD^1 F(1) - 1/12 * BD^2 F(1) - 1/24 * BDF^3 F(1)
!      =            F(1)
!        -  1/2 * ( F(1) -         F(0) )
!        - 1/12 * ( F(1) - 2     * F(0) +        F(-1) )
!        - 1/24 * ( F(1) - 3     * F(0) + 3    * F(-1) -        F(-2) )
!      =   9/24 *   F(1) + 19/24 * F(0) - 5/24 * F(-1) + 1/24 * F(-2)
!
!    which is the Adams-Moulton formula of order 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Simeon Fatunla,
!    Numerical Methods for Initial Value Problems in Ordinary
!    Differential Equations,
!    Academic Press, 1988,
!    ISBN: 0122499301,
!    LC: QA372.F35.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N should be between 1 and 19.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 19

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) w_save(n_max)
  real ( kind = 8 ) x(n)

  w_save(1) =                 1.0D+00
  w_save(2) =               - 1.0D+00 /               2.0D+00
  w_save(3) =               - 1.0D+00 /              12.0D+00
  w_save(4) =               - 1.0D+00 /              24.0D+00
  w_save(5) =              - 19.0D+00 /             720.0D+00
  w_save(6) =               - 3.0D+00 /             160.0D+00
  w_save(7) =             - 863.0D+00 /           60480.0D+00
  w_save(8) =             - 275.0D+00 /           24792.0D+00
  w_save(9) =           - 33953.0D+00 /         3628800.0D+00
  w_save(10) =           - 8183.0D+00 /         1036800.0D+00
  w_save(11) =        - 3250433.0D+00 /       479001600.0D+00
  w_save(12) =           - 4671.0D+00 /          788480.0D+00
  w_save(13) =    - 13695779093.0D+00 /   2615348736000.0D+00
  w_save(14) =     - 2224234463.0D+00 /    475517952000.0D+00
  w_save(15) =   - 132282840127.0D+00 /  31384184832000.0D+00
  w_save(16) =     - 2639651053.0D+00 /    689762304000.0D+00
  w_save(17) =  111956703448001.0D+00 /   3201186852864.0D+00
  w_save(18) =         50188465.0D+00 /     15613165568.0D+00
  w_save(19) = 2334028946344463.0D+00 / 786014494949376.0D+00

  if ( n_max < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BDFC_SET - Fatal error!'
    write ( *, '(a)' ) '  Input order N = ', n
    write ( *, '(a)' ) '  exceeds maximum order N_MAX = ', n_max
    stop
  end if

  w(1:n) = w_save(1:n)

  do i = 1, n
    x(i) = real ( 2 - i, kind = 8 )
  end do

  return
end
subroutine bdfp_set ( n, x, w )

!*****************************************************************************80
!
!! BDFP_SET sets weights for backward differentiation predictor quadrature.
!
!  Discussion:
!
!    A backward differentiation predictor formula is defined for a set
!    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
!    that the values of the function to be integrated are known at the
!    abscissas, the formula is written in terms of the function value at
!    X(2), and the backward differences at X(2) that approximate the
!    derivatives there.  A backward differentiation predictor formula
!    is equivalent to an Adams-Bashforth formula of the same order.
!
!    The integral:
!
!      integral ( 0 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * bd^(i-1) f ( 0 ),
!
!    Here, "BD^(I-1) F ( 0 )" denotes the (I-1)st backward difference
!    of F at X = 0, using a spacing of 1.  In particular,
!
!    BD^0 F(0) = F(0)
!    BD^1 F(0) = F(0) - F(-1)
!    BD^2 F(0) = F(0) - 2 * F(-1) + F(-2 )
!
!    The relationship between a backward difference predictor and the
!    corresponding Adams-Bashforth formula may be illustrated for the
!    BDF predictor of order 3:
!
!      BD^0 F(0) + 0.5 * BD^1 F(0) + 5/12 * BD^2 F(0)
!      =            F(0)
!        + 1/2  * ( F(0) -         F(1) )
!        + 5/12 * ( F(0) - 2     * F(-1) +      F(-2) )
!      =  23/12 *   F(0) - 16/12 * F(-1) + 5/12 F(-2)
!
!    which is the Adams-Bashforth formula of order 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Simeon Fatunla,
!    Numerical Methods for Initial Value Problems in Ordinary
!    Differential Equations,
!    Academic Press, 1988,
!    ISBN: 0122499301,
!    LC: QA372.F35.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N should be between 1 and 19.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 19

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) w_save(n_max)
  real ( kind = 8 ) x(n)

  w_save(1) =                       1.0D+00
  w_save(2) =                       1.0D+00 /                2.0D+00
  w_save(3) =                       5.0D+00 /               12.0D+00
  w_save(4) =                       3.0D+00 /                8.0D+00
  w_save(5) =                     251.0D+00 /              720.0D+00
  w_save(6) =                      95.0D+00 /              288.0D+00
  w_save(7) =                   19087.0D+00 /            60480.0D+00
  w_save(8) =                    5257.0D+00 /            17280.0D+00
  w_save(9) =                 1070017.0D+00 /          3628800.0D+00
  w_save(10) =                  25713.0D+00 /            89600.0D+00
  w_save(11) =               26842253.0D+00 /         95800320.0D+00
  w_save(12) =                4777223.0D+00 /         17418240.0D+00
  w_save(13) =           703604254357.0D+00 /    2615348736000.0D+00
  w_save(14) =           106364763817.0D+00 /     402361344000.0D+00
  w_save(15) =          1166309819657.0D+00 /    4483454976000.0D+00
  w_save(16) =               25221445.0D+00 /         98402304.0D+00
  w_save(17) =       8092989203533249.0D+00 /    3201186852864.0D+00
  w_save(18) =         85455477715379.0D+00 /      34237292544.0D+00
  w_save(19) =   12600467236042756559.0D+00 / 5109094217170944.0D+00

  if ( n_max < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BDFP_SET - Fatal error!'
    write ( *, '(a)' ) '  Input order N = ', n
    write ( *, '(a)' ) '  exceeds maximum order N_MAX = ', n_max
    stop
  end if

  w(1:n) = w_save(1:n)

  do i = 1, n
    x(i) = real ( 1 - i, kind = 8 )
  end do

  return
end
subroutine cheb_set ( n, x, w )

!*****************************************************************************80
!
!! CHEB_SET sets abscissas and weights for Chebyshev quadrature.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The Chebyshev rule is distinguished by using equal weights.
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Hermann Engels,
!    Numerical Quadrature and Cubature,
!    Academic Press, 1980,
!    ISBN: 012238850X,
!    LC: QA299.3E5.
!
!    Zdenek Kopal,
!    Numerical Analysis,
!    John Wiley, 1955,
!    LC: QA297.K6.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N may only have the values 1, 2, 3, 4, 5, 6, 7 or 9.
!    There are NO other Chebyshev rules with real abscissas.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = 0.0D+00

  else if ( n == 2 ) then

    x(1) = - 1.0D+00 / sqrt ( 3.0D+00 )
    x(2) =   1.0D+00 / sqrt ( 3.0D+00 )

  else if ( n == 3 ) then

    x(1) = - 1.0D+00 / sqrt ( 2.0D+00 )
    x(2) =   0.0D+00
    x(3) =   1.0D+00 / sqrt ( 2.0D+00 )

  else if ( n == 4 ) then

    x(1) =   - sqrt ( ( 1.0D+00 + 2.0D+00 / sqrt ( 5.0D+00 ) ) / 3.0D+00 )
    x(2) =   - sqrt ( ( 1.0D+00 - 2.0D+00 / sqrt ( 5.0D+00 ) ) / 3.0D+00 )
    x(3) =     sqrt ( ( 1.0D+00 - 2.0D+00 / sqrt ( 5.0D+00 ) ) / 3.0D+00 )
    x(4) =     sqrt ( ( 1.0D+00 + 2.0D+00 / sqrt ( 5.0D+00 ) ) / 3.0D+00 )

  else if ( n == 5 ) then

    x(1) = - sqrt ( ( 5.0D+00 + sqrt ( 11.0D+00 ) ) / 12.0D+00 )
    x(2) = - sqrt ( ( 5.0D+00 - sqrt ( 11.0D+00 ) ) / 12.0D+00 )
    x(3) =   0.0D+00
    x(4) =   sqrt ( ( 5.0D+00 - sqrt ( 11.0D+00 ) ) / 12.0D+00 )
    x(5) =   sqrt ( ( 5.0D+00 + sqrt ( 11.0D+00 ) ) / 12.0D+00 )

  else if ( n == 6 ) then

    x(1) = - 0.866246818107820591383598D+00
    x(2) = - 0.422518653761111529118546D+00
    x(3) = - 0.266635401516704720331534D+00
    x(4) =   0.266635401516704720331534D+00
    x(5) =   0.422518653761111529118546D+00
    x(6) =   0.866246818107820591383598D+00

  else if ( n == 7 ) then

    x(1) = - 0.883861700758049035704224D+00
    x(2) = - 0.529656775285156811385048D+00
    x(3) = - 0.323911810519907637519673D+00
    x(4) =   0.0D+00
    x(5) =   0.323911810519907637519673D+00
    x(6) =   0.529656775285156811385048D+00
    x(7) =   0.883861700758049035704224D+00

  else if ( n == 9 ) then

    x(1) = - 0.911589307728434473664949D+00
    x(2) = - 0.601018655380238071428128D+00
    x(3) = - 0.528761783057879993260181D+00
    x(4) = - 0.167906184214803943068031D+00
    x(5) =   0.0D+00
    x(6) =   0.167906184214803943068031D+00
    x(7) =   0.528761783057879993260181D+00
    x(8) =   0.601018655380238071428128D+00
    x(9) =   0.911589307728434473664949D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEB_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 through 7, and 9.'
    stop

  end if

  w(1:n) = 2.0D+00 / real ( n, kind = 8 )

  return
end
subroutine chebyshev1_compute ( n, x, w )

!*****************************************************************************80
!
!! CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) / sqrt ( 1 - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2008
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
!    Input, integer ( kind = 4 ) N, the order.
!    0 < N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEBYSHEV1_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop
  end if

  w(1:n) = pi / real ( n, kind = 8 )

  do i = 1, n
    x(i) = cos ( pi * real ( 2 * n + 1 - 2 * i, kind = 8 ) &
                    / real ( 2 * n, kind = 8 ) )
  end do

  return
end
subroutine chebyshev1_integral ( expon, exact )

!*****************************************************************************80
!
!! CHEBYSHEV1_INTEGRAL evaluates a monomial Chebyshev type 1 integral.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= +1 ) x^n / sqrt ( 1 - x * x ) dx
!
!    This routine is given the value of the exponent, and returns the
!    exact value of the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Output, real ( kind = 8 ) EXACT, the value of the exact integral.
!
  implicit none

  real ( kind = 8 ) bot
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) top
!
!  Get the exact value of the integral.
!
  if ( mod ( expon, 2 ) == 0 ) then

    top = 1
    bot = 1
    do i = 2, expon, 2
      top = top * ( i - 1 )
      bot = bot *   i
    end do
    
    exact = pi * real ( top, kind = 8 ) / real ( bot, kind = 8 )

  else

    exact = 0.0D+00
    
  end if

  return
end
subroutine chebyshev2_compute ( n, x, w )

!*****************************************************************************80
!
!! CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) sqrt ( 1 - x * x )  dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2008
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
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEBYSHEV2_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop
  end if

  do i = 1, n
    angle = pi * real ( n + 1 - i, kind = 8 ) / real ( n + 1, kind = 8 )
    w(i) = pi / real ( n + 1, kind = 8 ) * ( sin ( angle ) )**2
    x(i) = cos ( angle )
  end do

  return
end
subroutine chebyshev2_integral ( expon, exact )

!*****************************************************************************80
!
!! CHEBYSHEV2_INTEGRAL evaluates a monomial Chebyshev type 2 integral.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= +1 ) x^n * sqrt ( 1 - x * x ) dx
!
!    This routine is given the value of the exponent, and returns the
!    exact value of the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Output, real ( kind = 8 ) EXACT, the value of the exact integral.
!
  implicit none

  real ( kind = 8 ) bot
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) top
!
!  Get the exact value of the integral.
!
  if ( mod ( expon, 2 ) == 0 ) then

    top = 1
    bot = 1
    do i = 2, expon, 2
      top = top * ( i - 1 )
      bot = bot *   i
    end do

    bot = bot * real ( expon + 2, kind = 8 )

    exact = pi * real ( top, kind = 8 ) / real ( bot, kind = 8 )

  else

    exact = 0.0D+00
    
  end if

  return
end
subroutine chebyshev3_compute ( n, x, w )

!*****************************************************************************80
!
!! CHEBYSHEV3_COMPUTE computs a Gauss-Chebyshev type 3 quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) / sqrt ( 1 - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The N = 1 rule is exceptional.  It consists of a single
!    point at 0, with weight PI.
!
!    For rules with N = 2 or greater, the following remarks apply:
!
!    If N points are used, then Gauss-Chebyshev quadrature
!    will compute the integral exactly, whenever F(X) is a polynomial
!    of degree 2*N-3 or less.
!
!    The abscissas include -1 and 1.
!
!    The first and last weights are 0.5 * PI / ( N - 1),
!    and all other weights are PI / ( N - 1 ).
!
!    If the order is doubled, the abscissas of the new rule include
!    all the points of the old rule.  This fact can be used to
!    efficiently implement error estimation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEBYSHEV3_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 1.'
    write ( *, '(a,i8)' ) '  The input value was N = ', n
    stop
  end if
!
!  Take care of the special case N = 1.
!
  if ( n == 1 ) then
    x(1) = 0.0D+00
    w(1) = pi
    return
  end if

  do i = 1, n

    angle = real ( n - i, kind = 8 ) * pi &
          / real ( n - 1, kind = 8 )
    x(i) = cos ( angle )

  end do

  w(1)     = pi / real ( 2 * ( n - 1 ), kind = 8 )
  w(2:n-1) = pi / real (       n - 1,   kind = 8 )
  w(n)     = pi / real ( 2 * ( n - 1 ), kind = 8 )

  return
end
subroutine clenshaw_curtis_compute ( n, x, w )

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
!
!  Discussion:
!
!    This method uses a direct approach.  The paper by Waldvogel
!    exhibits a more efficient approach using Fourier transforms.
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The abscissas for the rule of order N can be regarded
!    as the cosines of equally spaced angles between 180 and 0 degrees:
!
!      X(I) = cos ( ( N - I ) * PI / ( N - 1 ) )
!
!    except for the basic case N = 1, when
!
!      X(1) = 0.
!
!    A Clenshaw-Curtis rule that uses N points will integrate
!    exactly all polynomials of degrees 0 through N-1.  If N
!    is odd, then by symmetry the polynomial of degree N will
!    also be integrated exactly.
!
!    If the value of N is increased in a sensible way, then
!    the new set of abscissas will include the old ones.  One such
!    sequence would be N(K) = 2*K+1 for K = 0, 1, 2, ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLENSHAW_CURTIS_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( n == 1 ) then
    x(1) = 0.0D+00
    w(1) = 2.0D+00
    return
  end if

  do i = 1, n
    theta(i) = real ( n - i, kind = 8 ) * pi &
             / real ( n - 1, kind = 8 )
  end do

  x(1:n) = cos ( theta(1:n) )

  do i = 1, n

    w(i) = 1.0D+00

    do j = 1, ( n - 1 ) / 2

      if ( 2 * j == ( n - 1 ) ) then
        b = 1.0D+00
      else
        b = 2.0D+00
      end if

      w(i) = w(i) - b * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta(i) ) &
           / real ( 4 * j * j - 1, kind = 8 )

    end do

  end do

  w(1)     =           w(1)     / real ( n - 1, kind = 8 )
  w(2:n-1) = 2.0D+00 * w(2:n-1) / real ( n - 1, kind = 8 )
  w(n)     =           w(n)     / real ( n - 1, kind = 8 )

  return
end
subroutine clenshaw_curtis_set ( n, x, w )

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_SET sets a Clenshaw-Curtis quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The abscissas for the rule of order N can be regarded
!    as the cosines of equally spaced angles between 180 and 0 degrees:
!
!      X(I) = cos ( ( I - 1 ) * PI / ( N - 1 ) )
!
!    except for the basic case N = 1, when
!
!      X(1) = 0.
!
!    A Clenshaw-Curtis rule that uses N points will integrate
!    exactly all polynomials of degrees 0 through N-1.  If N
!    is odd, then by symmetry the polynomial of degree N will
!    also be integrated exactly.
!
!    If the value of N is increased in a sensible way, then
!    the new set of abscissas will include the old ones.  One such
!    sequence would be N(K) = 2*K+1 for K = 0, 1, 2, ...
!    Thus, in the table below, the abscissas for order 9 include
!    those for order 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 17, 33, 65 or 129.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =  0.00000000000000000000D+00
    w(1) =  2.00000000000000000000D+00

  else if ( n == 2 ) then

    x(1) = -1.00000000000000000000D+00
    x(2) =  1.00000000000000000000D+00

    w(1) =  1.00000000000000000000D+00
    w(2) =  1.00000000000000000000D+00

  else if ( n == 3 ) then

    x(1) = -1.00000000000000000000D+00
    x(2) =  0.00000000000000000000D+00
    x(3) =  1.00000000000000000000D+00

    w(1) =  0.33333333333333333333D+00
    w(2) =  1.33333333333333333333D+00
    w(3) =  0.33333333333333333333D+00

  else if ( n == 4 ) then

    x(1) = -1.00000000000000000000D+00
    x(2) = -0.50000000000000000000D+00
    x(3) =  0.50000000000000000000D+00
    x(4) =  1.00000000000000000000D+00

    w(1) =  0.11111111111111111111D+00
    w(2) =  0.88888888888888888889D+00
    w(3) =  0.88888888888888888889D+00
    w(4) =  0.11111111111111111111D+00

  else if ( n == 5 ) then

    x(1) = -1.00000000000000000000D+00
    x(2) = -0.70710678118654752440D+00
    x(3) =  0.00000000000000000000D+00
    x(4) =  0.70710678118654752440D+00
    x(5) =  1.00000000000000000000D+00

    w(1) =  0.06666666666666666667D+00
    w(2) =  0.53333333333333333333D+00
    w(3) =  0.80000000000000000000D+00
    w(4) =  0.53333333333333333333D+00
    w(5) =  0.06666666666666666667D+00

  else if ( n == 6 ) then

    x(1) = -1.00000000000000000000D+00
    x(2) = -0.80901699437494742410D+00
    x(3) = -0.30901699437494742410D+00
    x(4) =  0.30901699437494742410D+00
    x(5) =  0.80901699437493732410D+00
    x(6) =  1.00000000000000000000D+00

    w(1) =  0.04000000000000000000D+00
    w(2) =  0.36074304120001121619D+00
    w(3) =  0.59925695879998878381D+00
    w(4) =  0.59925695879998878381D+00
    w(5) =  0.36074304120001121619D+00
    w(6) =  0.04000000000000000000D+00

  else if ( n == 7 ) then

    x(1) = -1.00000000000000000000D+00
    x(2) = -0.86602540378443864676D+00
    x(3) = -0.50000000000000000000D+00
    x(4) =  0.00000000000000000000D+00
    x(5) =  0.50000000000000000000D+00
    x(6) =  0.86602540378443864676D+00
    x(7) =  1.00000000000000000000D+00

    w(1) =  0.02857142857142857143D+00
    w(2) =  0.25396825396825396825D+00
    w(3) =  0.45714285714285714286D+00
    w(4) =  0.52063492063492063492D+00
    w(5) =  0.45714285714285714286D+00
    w(6) =  0.25396825396825396825D+00
    w(7) =  0.02857142857142857143D+00

  else if ( n == 8 ) then

    x(1) = -1.00000000000000000000D+00
    x(2) = -0.90096886790241912624D+00
    x(3) = -0.62348980185873353053D+00
    x(4) = -0.22252093395631440429D+00
    x(5) =  0.22252093395631440429D+00
    x(6) =  0.62348980185873353053D+00
    x(7) =  0.90096886790241910624D+00
    x(8) =  1.00000000000000000000D+00

    w(1) =  0.02040816326530612245D+00
    w(2) =  0.19014100721820835178D+00
    w(3) =  0.35224242371815911533D+00
    w(4) =  0.43720840579832641044D+00
    w(5) =  0.43720840579832641044D+00
    w(6) =  0.35224242371815911533D+00
    w(7) =  0.19014100721820835178D+00
    w(8) =  0.02040816326530612245D+00

  else if ( n == 9 ) then

    x(1) = -1.00000000000000000000D+00
    x(2) = -0.92387953251128675613D+00
    x(3) = -0.70710678118654752440D+00
    x(4) = -0.38268343236508977173D+00
    x(5) =  0.00000000000000000000D+00
    x(6) =  0.38268343236508977173D+00
    x(7) =  0.70710678118654752440D+00
    x(8) =  0.92387953251128675613D+00
    x(9) =  1.00000000000000000000D+00

    w(1) =  0.01587301587301587302D+00
    w(2) =  0.14621864921601815501D+00
    w(3) =  0.27936507936507936508D+00
    w(4) =  0.36171785872048978150D+00
    w(5) =  0.39365079365079365079D+00
    w(6) =  0.36171785872048978150D+00
    w(7) =  0.27936507936507936508D+00
    w(8) =  0.14621864921601815501D+00
    w(9) =  0.01587301587301587302D+00

  else if ( n == 10 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.93969262078590838405D+00
    x(3)  = -0.76604444311897903520D+00
    x(4)  = -0.50000000000000000000D+00
    x(5)  = -0.17364817766693034885D+00
    x(6)  =  0.17364817766693034885D+00
    x(7)  =  0.50000000000000000000D+00
    x(8)  =  0.76604444311897903520D+00
    x(9)  =  0.93969262078590838405D+00
    x(10) =  1.00000000000000000000D+00

    w(1)  =  0.01234567901234567901D+00
    w(2)  =  0.11656745657203712296D+00
    w(3)  =  0.22528432333810440813D+00
    w(4)  =  0.30194003527336860670D+00
    w(5)  =  0.34386250580414418320D+00
    w(6)  =  0.34386250580414418320D+00
    w(7)  =  0.30194003527336860670D+00
    w(8)  =  0.22528432333810440813D+00
    w(9)  =  0.11656745657203712296D+00
    w(10) =  0.01234567901234567901D+00

  else if ( n == 11 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.95105651629515357212D+00
    x(3)  = -0.80901699437494742410D+00
    x(4)  = -0.58778525229247312917D+00
    x(5)  = -0.30901699437494742410D+00
    x(6)  =  0.00000000000000000000D+00
    x(7)  =  0.30901699437494742410D+00
    x(8)  =  0.58778525229247312917D+00
    x(9)  =  0.80901699437494742410D+00
    x(10) =  0.95105651629515357212D+00
    x(11) =  1.00000000000000000000D+00

    w(1)  =  0.01010101010101010101D+00
    w(2)  =  0.09457905488370156116D+00
    w(3)  =  0.18563521442424776529D+00
    w(4)  =  0.25358833328368660623D+00
    w(5)  =  0.29921327042423708320D+00
    w(6)  =  0.31376623376623376623D+00
    w(7)  =  0.29921327042423708320D+00
    w(8)  =  0.25358833328368660623D+00
    w(9)  =  0.18563521442424776529D+00
    w(10) =  0.09457905488370156116D+00
    w(11) =  0.01010101010101010101D+00

  else if ( n == 12 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.95949297361449738989D+00
    x(3)  = -0.84125353283118116886D+00
    x(4)  = -0.65486073394528506406D+00
    x(5)  = -0.41541501300188642553D+00
    x(6)  = -0.14231483827328514044D+00
    x(7)  =  0.14231483827328514044D+00
    x(8)  =  0.41541501300188642553D+00
    x(9)  =  0.65486073394528506406D+00
    x(10) =  0.84125353283118116886D+00
    x(11) =  0.95949297361449738989D+00
    x(12) =  1.00000000000000000000D+00

    w(1)  =  0.00826446280991735537D+00
    w(2)  =  0.07856015374620000543D+00
    w(3)  =  0.15504045508256136552D+00
    w(4)  =  0.21556254600086858099D+00
    w(5)  =  0.25991734106691617602D+00
    w(6)  =  0.28265504129353651666D+00
    w(7)  =  0.28265504129353651666D+00
    w(8)  =  0.25991734106691617602D+00
    w(9)  =  0.21556254600086858099D+00
    w(10) =  0.15504045508256136552D+00
    w(11) =  0.07856015374620000543D+00
    w(12) =  0.00826446280991735537D+00

  else if ( n == 13 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.96592582628906828675D+00
    x(3)  = -0.86602540378443864676D+00
    x(4)  = -0.70710678118654752440D+00
    x(5)  = -0.50000000000000000000D+00
    x(6)  = -0.25881904510252076235D+00
    x(7)  =  0.00000000000000000000D+00
    x(8)  =  0.25881904510252076235D+00
    x(9)  =  0.50000000000000000000D+00
    x(10) =  0.70710678118654752440D+00
    x(11) =  0.86602540378443864676D+00
    x(12) =  0.96592582628906828675D+00
    x(13) =  1.00000000000000000000D+00

    w(1)  =  0.00699300699300699301D+00
    w(2)  =  0.06605742495207439452D+00
    w(3)  =  0.13154253154253154253D+00
    w(4)  =  0.18476338476338476338D+00
    w(5)  =  0.22697302697302697303D+00
    w(6)  =  0.25267569378104433860D+00
    w(7)  =  0.26198986198986198986D+00
    w(8)  =  0.25267569378104433860D+00
    w(9)  =  0.22697302697302697303D+00
    w(10) =  0.18476338476338476338D+00
    w(11) =  0.13154253154253154253D+00
    w(12) =  0.06605742495207439452D+00
    w(13) =  0.00699300699300699301D+00

  else if ( n == 14 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.97094181742605202716D+00
    x(3)  = -0.88545602565320989590D+00
    x(4)  = -0.74851074817110109863D+00
    x(5)  = -0.56806474673115580251D+00
    x(6)  = -0.35460488704253562597D+00
    x(7)  = -0.12053668025532305335D+00
    x(8)  =  0.12053668025532305335D+00
    x(9)  =  0.35460488704253562597D+00
    x(10) =  0.56806474673115580251D+00
    x(11) =  0.74851074817110109863D+00
    x(12) =  0.88545602565320989590D+00
    x(13) =  0.97094181742605202716D+00
    x(14) =  1.00000000000000000000D+00

    w(1)  =  0.00591715976331360947D+00
    w(2)  =  0.05646531376341444627D+00
    w(3)  =  0.11276867248985655881D+00
    w(4)  =  0.16003802611671868523D+00
    w(5)  =  0.19899241036578321848D+00
    w(6)  =  0.22590304977856444935D+00
    w(7)  =  0.23991536772234903239D+00
    w(8)  =  0.23991536772234903239D+00
    w(9)  =  0.22590304977856444935D+00
    w(10) =  0.19899241036578321848D+00
    w(11) =  0.16003802611671868523D+00
    w(12) =  0.11276867248985655881D+00
    w(13) =  0.05646531376341444627D+00
    w(14) =  0.00591715976331360947D+00

  else if ( n == 15 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.97492791218182360702D+00
    x(3)  = -0.90096886790241912624D+00
    x(4)  = -0.78183148246802980871D+00
    x(5)  = -0.62348980185873353053D+00
    x(6)  = -0.43388373911755812048D+00
    x(7)  = -0.22252093395631440429D+00
    x(8)  =  0.00000000000000000000D+00
    x(9)  =  0.22252093395631440429D+00
    x(10) =  0.43388373911755812048D+00
    x(11) =  0.62348980185873353053D+00
    x(12) =  0.78183148246802980871D+00
    x(13) =  0.90096886790241912624D+00
    x(14) =  0.97492791218182360702D+00
    x(15) =  1.00000000000000000000D+00

    w(1)  =  0.00512820512820512821D+00
    w(2)  =  0.04869938729508823855D+00
    w(3)  =  0.09782039167605215913D+00
    w(4)  =  0.13966507849560431803D+00
    w(5)  =  0.17560578900106674677D+00
    w(6)  =  0.20205146748238357364D+00
    w(7)  =  0.21888151163057340180D+00
    w(8)  =  0.22429633858205286777D+00
    w(9)  =  0.21888151163057340180D+00
    w(10) =  0.20205146748238357364D+00
    w(11) =  0.17560578900106674677D+00
    w(12) =  0.13966507849560431803D+00
    w(13) =  0.09782039167605215913D+00
    w(14) =  0.04869938729508823855D+00
    w(15) =  0.00512820512820512821D+00

  else if ( n == 16 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.97814760073380563793D+00
    x(3)  = -0.91354545764260089550D+00
    x(4)  = -0.80901699437494742410D+00
    x(5)  = -0.66913060635885821383D+00
    x(6)  = -0.50000000000000000000D+00
    x(7)  = -0.30901699437494742410D+00
    x(8)  = -0.10452846326765347140D+00
    x(9)  =  0.10452846326765347140D+00
    x(10) =  0.30901699437494742410D+00
    x(11) =  0.50000000000000000000D+00
    x(12) =  0.66913060635885821383D+00
    x(13) =  0.80901699437494742410D+00
    x(14) =  0.91354545764260089550D+00
    x(15) =  0.97814760073380563793D+00
    x(16) =  1.00000000000000000000D+00

    w(1)  =  0.00444444444444444444D+00
    w(2)  =  0.04251476624752508988D+00
    w(3)  =  0.08553884025933288291D+00
    w(4)  =  0.12294010082849361533D+00
    w(5)  =  0.15573317603967369176D+00
    w(6)  =  0.18132978132978132978D+00
    w(7)  =  0.19921478132638853955D+00
    w(8)  =  0.20828410952436040635D+00
    w(9)  =  0.20828410952436040635D+00
    w(10) =  0.19921478132638853955D+00
    w(11) =  0.18132978132978132978D+00
    w(12) =  0.15573317603967369176D+00
    w(13) =  0.12294010082849361533D+00
    w(14) =  0.08553884025933288291D+00
    w(15) =  0.04251476624752508988D+00
    w(16) =  0.00444444444444444444D+00

  else if ( n == 17 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.98078528040323044913D+00
    x(3)  = -0.92387953251128675613D+00
    x(4)  = -0.83146961230254523708D+00
    x(5)  = -0.70710678118654752440D+00
    x(6)  = -0.55557023301960222474D+00
    x(7)  = -0.38268343236508977173D+00
    x(8)  = -0.19509032201612826785D+00
    x(9)  =  0.00000000000000000000D+00
    x(10) =  0.19509032201612826785D+00
    x(11) =  0.38268343236508977173D+00
    x(12) =  0.55557023301960222474D+00
    x(13) =  0.70710678118654752440D+00
    x(14) =  0.83146961230254523708D+00
    x(15) =  0.92387953251128675613D+00
    x(16) =  0.98078528040323044913D+00
    x(17) =  1.00000000000000000000D+00

    w(1)  =  0.00392156862745098039D+00
    w(2)  =  0.03736870283720561032D+00
    w(3)  =  0.07548233154315183441D+00
    w(4)  =  0.10890555258189093044D+00
    w(5)  =  0.13895646836823307412D+00
    w(6)  =  0.16317266428170330256D+00
    w(7)  =  0.18147378423649335700D+00
    w(8)  =  0.19251386461292564687D+00
    w(9)  =  0.19641012582189052777D+00
    w(10) =  0.19251386461292564687D+00
    w(11) =  0.18147378423649335700D+00
    w(12) =  0.16317266428170330256D+00
    w(13) =  0.13895646836823307412D+00
    w(14) =  0.10890555258189093044D+00
    w(15) =  0.07548233154315183441D+00
    w(16) =  0.03736870283720561032D+00
    w(17) =  0.00392156862745098039D+00

  else if ( n == 33 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.99518472667219688624D+00
    x(3)  = -0.98078528040323044913D+00
    x(4)  = -0.95694033573220886494D+00
    x(5)  = -0.92387953251128675613D+00
    x(6)  = -0.88192126434835502971D+00
    x(7)  = -0.83146961230254523708D+00
    x(8)  = -0.77301045336273696081D+00
    x(9)  = -0.70710678118654752440D+00
    x(10) = -0.63439328416364549822D+00
    x(11) = -0.55557023301960222474D+00
    x(12) = -0.47139673682599764856D+00
    x(13) = -0.38268343236508977173D+00
    x(14) = -0.29028467725446236764D+00
    x(15) = -0.19509032201612826785D+00
    x(16) = -0.098017140329560601994D+00
    x(17) =  0.000000000000000000000D+00
    x(18) =  0.098017140329560601994D+00
    x(19) =  0.19509032201612826785D+00
    x(20) =  0.29028467725446236764D+00
    x(21) =  0.38268343236508977173D+00
    x(22) =  0.47139673682599764856D+00
    x(23) =  0.55557023301960222474D+00
    x(24) =  0.63439328416364549822D+00
    x(25) =  0.70710678118654752440D+00
    x(26) =  0.77301045336273696081D+00
    x(27) =  0.83146961230254523708D+00
    x(28) =  0.88192126434835502971D+00
    x(29) =  0.92387953251128675613D+00
    x(30) =  0.95694033573220886494D+00
    x(31) =  0.98078528040323044913D+00
    x(32) =  0.99518472667219688624D+00
    x(33) =  1.00000000000000000000D+00

    w(1)  =  0.00097751710654936461D+00
    w(2)  =  0.00939319796295501470D+00
    w(3)  =  0.01923424513268114918D+00
    w(4)  =  0.02845791667723369009D+00
    w(5)  =  0.03759434191404720602D+00
    w(6)  =  0.04626276283775174949D+00
    w(7)  =  0.05455501630398031044D+00
    w(8)  =  0.06227210954529400455D+00
    w(9)  =  0.06942757563043545090D+00
    w(10) =  0.07588380044138847048D+00
    w(11) =  0.08163481765493851023D+00
    w(12) =  0.08657753844182743544D+00
    w(13) =  0.09070611286772099874D+00
    w(14) =  0.09394324443876873573D+00
    w(15) =  0.09629232594548817919D+00
    w(16) =  0.09769818820805558182D+00
    w(17) =  0.09817857778176829677D+00
    w(18) =  0.09769818820805558182D+00
    w(19) =  0.09629232594548817919D+00
    w(20) =  0.09394324443876873573D+00
    w(21) =  0.09070611286772099874D+00
    w(22) =  0.08657753844182743544D+00
    w(23) =  0.08163481765493851023D+00
    w(24) =  0.07588380044138847048D+00
    w(25) =  0.06942757563043545090D+00
    w(26) =  0.06227210954529400455D+00
    w(27) =  0.05455501630398031044D+00
    w(28) =  0.04626276283775174949D+00
    w(29) =  0.03759434191404720602D+00
    w(30) =  0.02845791667723369009D+00
    w(31) =  0.01923424513268114918D+00
    w(32) =  0.00939319796295501470D+00
    w(33) =  0.00097751710654936461D+00

  else if ( n == 65 ) then

    x(1)  = -1.00000000000000000000D+00
    x(2)  = -0.99879545620517239271D+00
    x(3)  = -0.99518472667219688624D+00
    x(4)  = -0.98917650996478097345D+00
    x(5)  = -0.98078528040323044913D+00
    x(6)  = -0.97003125319454399260D+00
    x(7)  = -0.95694033573220886494D+00
    x(8)  = -0.94154406518302077841D+00
    x(9)  = -0.92387953251128675613D+00
    x(10) = -0.90398929312344333159D+00
    x(11) = -0.88192126434835502971D+00
    x(12) = -0.85772861000027206990D+00
    x(13) = -0.83146961230254523708D+00
    x(14) = -0.80320753148064490981D+00
    x(15) = -0.77301045336273696081D+00
    x(16) = -0.74095112535495909118D+00
    x(17) = -0.70710678118654752440D+00
    x(18) = -0.67155895484701840063D+00
    x(19) = -0.63439328416364549822D+00
    x(20) = -0.59569930449243334347D+00
    x(21) = -0.55557023301960222474D+00
    x(22) = -0.51410274419322172659D+00
    x(23) = -0.47139673682599764856D+00
    x(24) = -0.42755509343028209432D+00
    x(25) = -0.38268343236508977173D+00
    x(26) = -0.33688985339222005069D+00
    x(27) = -0.29028467725446236764D+00
    x(28) = -0.24298017990326388995D+00
    x(29) = -0.19509032201612826785D+00
    x(30) = -0.14673047445536175166D+00
    x(31) = -0.098017140329560601994D+00
    x(32) = -0.049067674327418014255D+00
    x(33) =  0.000000000000000000000D+00
    x(34) =  0.049067674327418014255D+00
    x(35) =  0.098017140329560601994D+00
    x(36) =  0.14673047445536175166D+00
    x(37) =  0.19509032201612826785D+00
    x(38) =  0.24298017990326388995D+00
    x(39) =  0.29028467725446236764D+00
    x(40) =  0.33688985339222005069D+00
    x(41) =  0.38268343236508977173D+00
    x(42) =  0.42755509343028209432D+00
    x(43) =  0.47139673682599764856D+00
    x(44) =  0.51410274419322172659D+00
    x(45) =  0.55557023301960222474D+00
    x(46) =  0.59569930449243334347D+00
    x(47) =  0.63439328416364549822D+00
    x(48) =  0.67155895484701840063D+00
    x(49) =  0.70710678118654752440D+00
    x(50) =  0.74095112535495909118D+00
    x(51) =  0.77301045336273696081D+00
    x(52) =  0.80320753148064490981D+00
    x(53) =  0.83146961230254523708D+00
    x(54) =  0.85772861000027206990D+00
    x(55) =  0.88192126434835502971D+00
    x(56) =  0.90398929312344333159D+00
    x(57) =  0.92387953251128675613D+00
    x(58) =  0.94154406518302077841D+00
    x(59) =  0.95694033573220886494D+00
    x(60) =  0.97003125319454399260D+00
    x(61) =  0.98078528040323044913D+00
    x(62) =  0.98917650996478097345D+00
    x(63) =  0.99518472667219688624D+00
    x(64) =  0.99879545620517239271D+00
    x(65) =  1.00000000000000000000D+00

    w(1)  =  0.00024420024420024420D+00
    w(2)  =  0.00235149067531170332D+00
    w(3)  =  0.00483146544879091264D+00
    w(4)  =  0.00719269316173611402D+00
    w(5)  =  0.00958233879528379039D+00
    w(6)  =  0.01192339471421277160D+00
    w(7)  =  0.01425206043235199679D+00
    w(8)  =  0.01653498765728958965D+00
    w(9)  =  0.01878652974179578354D+00
    w(10) =  0.02098627442973743378D+00
    w(11) =  0.02314069493435819848D+00
    w(12) =  0.02523506498175476590D+00
    w(13) =  0.02727225714146838686D+00
    w(14) =  0.02924065319746833770D+00
    w(15) =  0.03114129710406762447D+00
    w(16) =  0.03296454656997632997D+00
    w(17) =  0.03471049818092511427D+00
    w(18) =  0.03637092028663918309D+00
    w(19) =  0.03794545992128481711D+00
    w(20) =  0.03942698871295609976D+00
    w(21) =  0.04081501340035783384D+00
    w(22) =  0.04210333111141810203D+00
    w(23) =  0.04329151496169082935D+00
    w(24) =  0.04437417923925731580D+00
    w(25) =  0.04535110955166067221D+00
    w(26) =  0.04621766751092557684D+00
    w(27) =  0.04697395904661414870D+00
    w(28) =  0.04761604458525019296D+00
    w(29) =  0.04814443257251220341D+00
    w(30) =  0.04855584485714105274D+00
    w(31) =  0.04885125664306609371D+00
    w(32) =  0.04902801843102555294D+00
    w(33) =  0.04908762351494245585D+00
    w(34) =  0.04902801843102555294D+00
    w(35) =  0.04885125664306609371D+00
    w(36) =  0.04855584485714105274D+00
    w(37) =  0.04814443257251220341D+00
    w(38) =  0.04761604458525019296D+00
    w(39) =  0.04697395904661414870D+00
    w(40) =  0.04621766751092557684D+00
    w(41) =  0.04535110955166067221D+00
    w(42) =  0.04437417923925731580D+00
    w(43) =  0.04329151496169082935D+00
    w(44) =  0.04210333111141810203D+00
    w(45) =  0.04081501340035783384D+00
    w(46) =  0.03942698871295609976D+00
    w(47) =  0.03794545992128481711D+00
    w(48) =  0.03637092028663918309D+00
    w(49) =  0.03471049818092511427D+00
    w(50) =  0.03296454656997632997D+00
    w(51) =  0.03114129710406762447D+00
    w(52) =  0.02924065319746833770D+00
    w(53) =  0.02727225714146838686D+00
    w(54) =  0.02523506498175476590D+00
    w(55) =  0.02314069493435819848D+00
    w(56) =  0.02098627442973743378D+00
    w(57) =  0.01878652974179578354D+00
    w(58) =  0.01653498765728958965D+00
    w(59) =  0.01425206043235199679D+00
    w(60) =  0.01192339471421277160D+00
    w(61) =  0.00958233879528379039D+00
    w(62) =  0.00719269316173611402D+00
    w(63) =  0.00483146544879091264D+00
    w(64) =  0.00235149067531170332D+00
    w(65) =  0.00024420024420024420D+00

  else if ( n == 129 ) then

    x(1)   = -1.00000000000000000000D+00
    x(2)   = -0.99969881869620422012D+00
    x(3)   = -0.99879545620517239271D+00
    x(4)   = -0.99729045667869021614D+00
    x(5)   = -0.99518472667219688624D+00
    x(6)   = -0.99247953459870999816D+00
    x(7)   = -0.98917650996478097345D+00
    x(8)   = -0.98527764238894124477D+00
    x(9)   = -0.98078528040323044913D+00
    x(10)  = -0.97570213003852854446D+00
    x(11)  = -0.97003125319454399260D+00
    x(12)  = -0.96377606579543986669D+00
    x(13)  = -0.95694033573220886494D+00
    x(14)  = -0.94952818059303666720D+00
    x(15)  = -0.94154406518302077841D+00
    x(16)  = -0.93299279883473888771D+00
    x(17)  = -0.92387953251128675613D+00
    x(18)  = -0.91420975570353065464D+00
    x(19)  = -0.90398929312344333159D+00
    x(20)  = -0.89322430119551532034D+00
    x(21)  = -0.88192126434835502971D+00
    x(22)  = -0.87008699110871141865D+00
    x(23)  = -0.85772861000027206990D+00
    x(24)  = -0.84485356524970707326D+00
    x(25)  = -0.83146961230254523708D+00
    x(26)  = -0.81758481315158369650D+00
    x(27)  = -0.80320753148064490981D+00
    x(28)  = -0.78834642762660626201D+00
    x(29)  = -0.77301045336273696081D+00
    x(30)  = -0.75720884650648454758D+00
    x(31)  = -0.74095112535495909118D+00
    x(32)  = -0.72424708295146692094D+00
    x(33)  = -0.70710678118654752440D+00
    x(34)  = -0.68954054473706692462D+00
    x(35)  = -0.67155895484701840063D+00
    x(36)  = -0.65317284295377676408D+00
    x(37)  = -0.63439328416364549822D+00
    x(38)  = -0.61523159058062684548D+00
    x(39)  = -0.59569930449243334347D+00
    x(40)  = -0.57580819141784530075D+00
    x(41)  = -0.55557023301960222474D+00
    x(42)  = -0.53499761988709721066D+00
    x(43)  = -0.51410274419322172659D+00
    x(44)  = -0.49289819222978403687D+00
    x(45)  = -0.47139673682599764856D+00
    x(46)  = -0.44961132965460660005D+00
    x(47)  = -0.42755509343028209432D+00
    x(48)  = -0.40524131400498987091D+00
    x(49)  = -0.38268343236508977173D+00
    x(50)  = -0.35989503653498814878D+00
    x(51)  = -0.33688985339222005069D+00
    x(52)  = -0.31368174039889147666D+00
    x(53)  = -0.29028467725446236764D+00
    x(54)  = -0.26671275747489838633D+00
    x(55)  = -0.24298017990326388995D+00
    x(56)  = -0.21910124015686979723D+00
    x(57)  = -0.19509032201612826785D+00
    x(58)  = -0.17096188876030122636D+00
    x(59)  = -0.14673047445536175166D+00
    x(60)  = -0.12241067519921619850D+00
    x(61)  = -0.098017140329560601994D+00
    x(62)  = -0.073564563599667423529D+00
    x(63)  = -0.049067674327418014255D+00
    x(64)  = -0.024541228522912288032D+00
    x(65)  =  0.00000000000000000000D+00
    x(66)  =  0.024541228522912288032D+00
    x(67)  =  0.049067674327418014255D+00
    x(68)  =  0.073564563599667423529D+00
    x(69)  =  0.098017140329560601994D+00
    x(70)  =  0.12241067519921619850D+00
    x(71)  =  0.14673047445536175166D+00
    x(72)  =  0.17096188876030122636D+00
    x(73)  =  0.19509032201612826785D+00
    x(74)  =  0.21910124015686979723D+00
    x(75)  =  0.24298017990326388995D+00
    x(76)  =  0.26671275747489838633D+00
    x(77)  =  0.29028467725446236764D+00
    x(78)  =  0.31368174039889147666D+00
    x(79)  =  0.33688985339222005069D+00
    x(80)  =  0.35989503653498814878D+00
    x(81)  =  0.38268343236508977173D+00
    x(82)  =  0.40524131400498987091D+00
    x(83)  =  0.42755509343028209432D+00
    x(84)  =  0.44961132965460660005D+00
    x(85)  =  0.47139673682599764856D+00
    x(86)  =  0.49289819222978403687D+00
    x(87)  =  0.51410274419322172659D+00
    x(88)  =  0.53499761988709721066D+00
    x(89)  =  0.55557023301960222474D+00
    x(90)  =  0.57580819141784530075D+00
    x(91)  =  0.59569930449243334347D+00
    x(92)  =  0.61523159058062684548D+00
    x(93)  =  0.63439328416364549822D+00
    x(94)  =  0.65317284295377676408D+00
    x(95)  =  0.67155895484701840063D+00
    x(96)  =  0.68954054473706692462D+00
    x(97)  =  0.70710678118654752440D+00
    x(98)  =  0.72424708295146692094D+00
    x(99)  =  0.74095112535495909118D+00
    x(100) =  0.75720884650648454758D+00
    x(101) =  0.77301045336273696081D+00
    x(102) =  0.78834642762660626201D+00
    x(103) =  0.80320753148064490981D+00
    x(104) =  0.81758481315158369650D+00
    x(105) =  0.83146961230254523708D+00
    x(106) =  0.84485356524970707326D+00
    x(107) =  0.85772861000027206990D+00
    x(108) =  0.87008699110871141865D+00
    x(109) =  0.88192126434835502971D+00
    x(110) =  0.89322430119551532034D+00
    x(111) =  0.90398929312344333159D+00
    x(112) =  0.91420975570353065464D+00
    x(113) =  0.92387953251128675613D+00
    x(114) =  0.93299279883473888771D+00
    x(115) =  0.94154406518302077841D+00
    x(116) =  0.94952818059303666720D+00
    x(117) =  0.95694033573220886494D+00
    x(118) =  0.96377606579543986669D+00
    x(119) =  0.97003125319454399260D+00
    x(120) =  0.97570213003852854446D+00
    x(121) =  0.98078528040323044913D+00
    x(122) =  0.98527764238894124477D+00
    x(123) =  0.98917650996478097345D+00
    x(124) =  0.99247953459870999816D+00
    x(125) =  0.99518472667219688624D+00
    x(126) =  0.99729045667869021614D+00
    x(127) =  0.99879545620517239271D+00
    x(128) =  0.99969881869620422012D+00
    x(129) =  1.00000000000000000000D+00

    w(1)   =  0.00006103888176768602D+00
    w(2)   =  0.00058807215382869754D+00
    w(3)   =  0.00120930061875273991D+00
    w(4)   =  0.00180308126695362360D+00
    w(5)   =  0.00240715327877140915D+00
    w(6)   =  0.00300345869904497128D+00
    w(7)   =  0.00360197835812614147D+00
    w(8)   =  0.00419553798718534675D+00
    w(9)   =  0.00478862143341336763D+00
    w(10)  =  0.00537724746840184621D+00
    w(11)  =  0.00596388034730799521D+00
    w(12)  =  0.00654590843862298928D+00
    w(13)  =  0.00712483332325489785D+00
    w(14)  =  0.00769875778896082811D+00
    w(15)  =  0.00826865154203087108D+00
    w(16)  =  0.00883303867470133581D+00
    w(17)  =  0.00939256583934814871D+00
    w(18)  =  0.00994602784923457905D+00
    w(19)  =  0.01049386202576892125D+00
    w(20)  =  0.01103504877427254184D+00
    w(21)  =  0.01156988348290849967D+00
    w(22)  =  0.01209748052807164113D+00
    w(23)  =  0.01261803597977743271D+00
    w(24)  =  0.01313076516693974630D+00
    w(25)  =  0.01363579321293772047D+00
    w(26)  =  0.01413241437853094133D+00
    w(27)  =  0.01462070254634350205D+00
    w(28)  =  0.01510001572479266783D+00
    w(29)  =  0.01557039073899425960D+00
    w(30)  =  0.01603123858745057916D+00
    w(31)  =  0.01648256956220377909D+00
    w(32)  =  0.01692383985846499368D+00
    w(33)  =  0.01735504125411394958D+00
    w(34)  =  0.01777566938875279997D+00
    w(35)  =  0.01818570377926339481D+00
    w(36)  =  0.01858467519566908661D+00
    w(37)  =  0.01897255587067948426D+00
    w(38)  =  0.01934890842392451844D+00
    w(39)  =  0.01971370183700155725D+00
    w(40)  =  0.02006652805198357604D+00
    w(41)  =  0.02040735612003867863D+00
    w(42)  =  0.02073580533490147816D+00
    w(43)  =  0.02105184759002011131D+00
    w(44)  =  0.02135512797425970725D+00
    w(45)  =  0.02164562356712882440D+00
    w(46)  =  0.02192300400598756892D+00
    w(47)  =  0.02218725355897195088D+00
    w(48)  =  0.02243806539722630184D+00
    w(49)  =  0.02267543270456671718D+00
    w(50)  =  0.02289907134390605882D+00
    w(51)  =  0.02310898491627407168D+00
    w(52)  =  0.02330491126131143273D+00
    w(53)  =  0.02348686571193163505D+00
    w(54)  =  0.02365460746057766523D+00
    w(55)  =  0.02380816473024258975D+00
    w(56)  =  0.02394731750476901502D+00
    w(57)  =  0.02407210792327850000D+00
    w(58)  =  0.02418233623893147567D+00
    w(59)  =  0.02427805942075745923D+00
    w(60)  =  0.02435909748927643184D+00
    w(61)  =  0.02442552306156708690D+00
    w(62)  =  0.02447717542743444284D+00
    w(63)  =  0.02451414358881568292D+00
    w(64)  =  0.02453628559651495473D+00
    w(65)  =  0.02454370750551418263D+00
    w(66)  =  0.02453628559651495473D+00
    w(67)  =  0.02451414358881568292D+00
    w(68)  =  0.02447717542743444284D+00
    w(69)  =  0.02442552306156708690D+00
    w(70)  =  0.02435909748927643184D+00
    w(71)  =  0.02427805942075745923D+00
    w(72)  =  0.02418233623893147567D+00
    w(73)  =  0.02407210792327850000D+00
    w(74)  =  0.02394731750476901502D+00
    w(75)  =  0.02380816473024258975D+00
    w(76)  =  0.02365460746057766523D+00
    w(77)  =  0.02348686571193163505D+00
    w(78)  =  0.02330491126131143273D+00
    w(79)  =  0.02310898491627407168D+00
    w(80)  =  0.02289907134390605882D+00
    w(81)  =  0.02267543270456671718D+00
    w(82)  =  0.02243806539722630184D+00
    w(83)  =  0.02218725355897195088D+00
    w(84)  =  0.02192300400598756892D+00
    w(85)  =  0.02164562356712882440D+00
    w(86)  =  0.02135512797425970725D+00
    w(87)  =  0.02105184759002011131D+00
    w(88)  =  0.02073580533490147816D+00
    w(89)  =  0.02040735612003867863D+00
    w(90)  =  0.02006652805198357604D+00
    w(91)  =  0.01971370183700155725D+00
    w(92)  =  0.01934890842392451844D+00
    w(93)  =  0.01897255587067948426D+00
    w(94)  =  0.01858467519566908661D+00
    w(95)  =  0.01818570377926339481D+00
    w(96)  =  0.01777566938875279997D+00
    w(97)  =  0.01735504125411394958D+00
    w(98)  =  0.01692383985846499368D+00
    w(99)  =  0.01648256956220377909D+00
    w(100) =  0.01603123858745057916D+00
    w(101) =  0.01557039073899425960D+00
    w(102) =  0.01510001572479266783D+00
    w(103) =  0.01462070254634350205D+00
    w(104) =  0.01413241437853094133D+00
    w(105) =  0.01363579321293772047D+00
    w(106) =  0.01313076516693974630D+00
    w(107) =  0.01261803597977743271D+00
    w(108) =  0.01209748052807164113D+00
    w(109) =  0.01156988348290849967D+00
    w(110) =  0.01103504877427254184D+00
    w(111) =  0.01049386202576892125D+00
    w(112) =  0.00994602784923457905D+00
    w(113) =  0.00939256583934814871D+00
    w(114) =  0.00883303867470133581D+00
    w(115) =  0.00826865154203087108D+00
    w(116) =  0.00769875778896082811D+00
    w(117) =  0.00712483332325489785D+00
    w(118) =  0.00654590843862298928D+00
    w(119) =  0.00596388034730799521D+00
    w(120) =  0.00537724746840184621D+00
    w(121) =  0.00478862143341336763D+00
    w(122) =  0.00419553798718534675D+00
    w(123) =  0.00360197835812614147D+00
    w(124) =  0.00300345869904497128D+00
    w(125) =  0.00240715327877140915D+00
    w(126) =  0.00180308126695362360D+00
    w(127) =  0.00120930061875273991D+00
    w(128) =  0.00058807215382869754D+00
    w(129) =  0.00006103888176768602D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLENSHAW_CURTIS_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 to 17, 33, 65 or 129'
    stop

  end if

  return
end
subroutine fejer1_compute ( n, x, w )

!*****************************************************************************80
!
!! FEJER1_COMPUTE computes a Fejer type 1 quadrature rule.
!
!  Discussion:
!
!    This method uses a direct approach.  The paper by Waldvogel
!    exhibits a more efficient approach using Fourier transforms.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2006
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
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEJER1_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( n == 1 ) then
    x(1) = 0.0D+00
    w(1) = 2.0D+00
    return
  end if

  do i = 1, n
    theta(i) = real ( 2 * ( n + 1 - i ) - 1, kind = 8 ) * pi &
             / real ( 2 * n,     kind = 8 )
  end do

  x(1:n) = cos ( theta(1:n) )

  do i = 1, n
    w(i) = 1.0D+00
    do j = 1, ( n / 2 )
      w(i) = w(i) - 2.0D+00 &
        * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta(i) ) &
        / real ( 4 * j * j - 1, kind = 8 )
    end do
  end do

  w(1:n) = 2.0D+00 * w(1:n) / real ( n, kind = 8 )

  return
end
subroutine fejer1_set ( n, x, w )

!*****************************************************************************80
!
!! FEJER1_SET sets abscissas and weights for Fejer type 1 quadrature.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2007
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
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N should be between 1 and 10.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1)   =  0.000000000000000D+00
    w(1) =  2.000000000000000D+00

  else if ( n == 2 ) then

    x(1) =   -0.7071067811865475D+00
    x(2) =    0.7071067811865475D+00

    w(1) =  1.000000000000000D+00
    w(2) =  1.000000000000000D+00

  else if ( n == 3 ) then

    x( 1) =  -0.8660254037844387D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   0.8660254037844387D+00

    w(1) =  0.4444444444444444D+00
    w(2) =  1.1111111111111111D+00
    w(3) =  0.4444444444444444D+00

  else if ( n == 4 ) then

    x( 1) =  -0.9238795325112867D+00
    x( 2) =  -0.3826834323650897D+00
    x( 3) =   0.3826834323650898D+00
    x( 4) =   0.9238795325112867D+00

    w( 1) = 0.2642977396044841D+00
    w( 2) = 0.7357022603955158D+00
    w( 3) = 0.7357022603955158D+00
    w( 4) = 0.2642977396044841D+00

  else if ( n == 5 ) then

    x( 1) =  -0.9510565162951535D+00
    x( 2) =  -0.5877852522924730D+00
    x( 3) =   0.0000000000000000D+00
    x( 4) =   0.5877852522924731D+00
    x( 5) =   0.9510565162951535D+00

    w( 1) = 0.1677812284666835D+00
    w( 2) = 0.5255521048666498D+00
    w( 3) = 0.6133333333333333D+00
    w( 4) = 0.5255521048666498D+00
    w( 5) = 0.1677812284666835D+00

  else if ( n == 6 ) then

    x( 1) =  -0.9659258262890682D+00
    x( 2) =  -0.7071067811865475D+00
    x( 3) =  -0.2588190451025206D+00
    x( 4) =   0.2588190451025207D+00
    x( 5) =   0.7071067811865476D+00
    x( 6) =   0.9659258262890683D+00

    w( 1) = 0.1186610213812358D+00
    w( 2) = 0.3777777777777778D+00
    w( 3) = 0.5035612008409863D+00
    w( 4) = 0.5035612008409863D+00
    w( 5) = 0.3777777777777778D+00
    w( 6) = 0.1186610213812358D+00

  else if ( n == 7 ) then

    x( 1) =  -0.9749279121818237D+00
    x( 2) =  -0.7818314824680295D+00
    x( 3) =  -0.4338837391175581D+00
    x( 4) =   0.0000000000000000D+00
    x( 5) =   0.4338837391175582D+00
    x( 6) =   0.7818314824680298D+00
    x( 7) =   0.9749279121818236D+00

    w( 1) = 0.08671618072672234D+00
    w( 2) = 0.2878313947886919D+00
    w( 3) = 0.3982415401308441D+00
    w( 4) = 0.4544217687074830D+00
    w( 5) = 0.3982415401308441D+00
    w( 6) = 0.2878313947886919D+00
    w( 7) = 0.08671618072672234D+00

  else if ( n == 8 ) then

    x( 1) =  -0.9807852804032304D+00
    x( 2) =  -0.8314696123025453D+00
    x( 3) =  -0.5555702330196020D+00
    x( 4) =  -0.1950903220161282D+00
    x( 5) =   0.1950903220161283D+00
    x( 6) =   0.5555702330196023D+00
    x( 7) =   0.8314696123025452D+00
    x( 8) =   0.9807852804032304D+00

    w( 1) = 0.06698294569858981D+00
    w( 2) = 0.2229879330145788D+00
    w( 3) = 0.3241525190645244D+00
    w( 4) = 0.3858766022223071D+00
    w( 5) = 0.3858766022223071D+00
    w( 6) = 0.3241525190645244D+00
    w( 7) = 0.2229879330145788D+00
    w( 8) = 0.06698294569858981D+00

  else if ( n == 9 ) then

    x( 1) =  -0.9848077530122080D+00
    x( 2) =  -0.8660254037844385D+00
    x( 3) =  -0.6427876096865394D+00
    x( 4) =  -0.3420201433256685D+00
    x( 5) =   0.0000000000000000D+00
    x( 6) =   0.3420201433256688D+00
    x( 7) =   0.6427876096865394D+00
    x( 8) =   0.8660254037844387D+00
    x( 9) =   0.9848077530122080D+00

    w( 1) = 0.05273664990990676D+00
    w( 2) = 0.1791887125220458D+00
    w( 3) = 0.2640372225410044D+00
    w( 4) = 0.3308451751681364D+00
    w( 5) = 0.3463844797178130D+00
    w( 6) = 0.3308451751681364D+00
    w( 7) = 0.2640372225410044D+00
    w( 8) = 0.1791887125220458D+00
    w( 9) = 0.05273664990990676D+00

  else if ( n == 10 ) then

    x( 1) =  -0.9876883405951377D+00
    x( 2) =  -0.8910065241883678D+00
    x( 3) =  -0.7071067811865475D+00
    x( 4) =  -0.4539904997395467D+00
    x( 5) =  -0.1564344650402306D+00
    x( 6) =   0.1564344650402309D+00
    x( 7) =   0.4539904997395468D+00
    x( 8) =   0.7071067811865476D+00
    x( 9) =   0.8910065241883679D+00
    x(10) =   0.9876883405951378D+00

    w( 1) = 0.04293911957413078D+00
    w( 2) = 0.1458749193773909D+00
    w( 3) = 0.2203174603174603D+00
    w( 4) = 0.2808792186638755D+00
    w( 5) = 0.3099892820671425D+00
    w( 6) = 0.3099892820671425D+00
    w( 7) = 0.2808792186638755D+00
    w( 8) = 0.2203174603174603D+00
    w( 9) = 0.1458749193773909D+00
    w(10) = 0.04293911957413078D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEJER1_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 through 10.'
    stop

  end if

  return
end
subroutine fejer2_compute ( n, x, w )

!*****************************************************************************80
!
!! FEJER2_COMPUTE computes a Fejer type 2 quadrature rule.
!
!  Discussion:
!
!    This method uses a direct approach.  The paper by Waldvogel
!    exhibits a more efficient approach using Fourier transforms.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2007
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
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEJER2_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( n == 1 ) then
    x(1) = 0.0D+00
    w(1) = 2.0D+00
    return
  else if ( n == 2 ) then
    x(1) = -0.5D+00
    x(2) =  0.5D+00
    w(1:2) = 1.0D+00
    return
  end if

  do i = 1, n
    theta(i) = real ( n + 1 - i, kind = 8 ) * pi &
             / real ( n + 1, kind = 8 )
  end do

  x(1:n) = cos ( theta(1:n) )

  do i = 1, n

    w(i) = 1.0D+00

    do j = 1, ( ( n - 1 ) / 2 )
      w(i) = w(i) - 2.0D+00 &
        * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta(i) ) &
        / real ( 4 * j * j - 1, kind = 8 )
    end do

    if ( 2 < n ) then
      p = 2.0D+00 * real ( ( ( n + 1 ) / 2 ), kind = 8 ) - 1.0D+00
      w(i) = w(i) - cos ( ( p + 1.0D+00 ) * theta(i) ) / p
    end if

  end do

  w(1:n) = 2.0D+00 * w(1:n) / real ( n + 1, kind = 8 )

  return
end
subroutine fejer2_set ( n, x, w )

!*****************************************************************************80
!
!! FEJER2_SET sets abscissas and weights for Fejer type 2 quadrature.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2007
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
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N should be between 1 and 10.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1)   =  0.000000000000000D+00
    w(1) =  2.000000000000000D+00

  else if ( n == 2 ) then

    x(1) =   -0.5000000000000000D+00
    x(2) =    0.5000000000000000D+00

    w(1) =  1.0000000000000000D+00
    w(2) =  1.0000000000000000D+00

  else if ( n == 3 ) then

    x( 1) =  -0.7071067811865476D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   0.7071067811865476D+00

    w(1) =  0.6666666666666666D+00
    w(2) =  0.6666666666666666D+00
    w(3) =  0.6666666666666666D+00

  else if ( n == 4 ) then

    x( 1) =  -0.8090169943749475D+00
    x( 2) =  -0.3090169943749475D+00
    x( 3) =   0.3090169943749475D+00
    x( 4) =   0.8090169943749475D+00

    w( 1) = 0.4254644007500070D+00
    w( 2) = 0.5745355992499930D+00
    w( 3) = 0.5745355992499930D+00
    w( 4) = 0.4254644007500070D+00

  else if ( n == 5 ) then

    x( 1) =  -0.8660254037844387D+00
    x( 2) =  -0.5000000000000000D+00
    x( 3) =   0.0000000000000000D+00
    x( 4) =   0.5000000000000000D+00
    x( 5) =   0.8660254037844387D+00

    w( 1) = 0.3111111111111111D+00
    w( 2) = 0.4000000000000000D+00
    w( 3) = 0.5777777777777777D+00
    w( 4) = 0.4000000000000000D+00
    w( 5) = 0.3111111111111111D+00

  else if ( n == 6 ) then

    x( 1) =  -0.9009688679024191D+00
    x( 2) =  -0.6234898018587336D+00
    x( 3) =  -0.2225209339563144D+00
    x( 4) =   0.2225209339563144D+00
    x( 5) =   0.6234898018587336D+00
    x( 6) =   0.9009688679024191D+00

    w( 1) = 0.2269152467244296D+00
    w( 2) = 0.3267938603769863D+00
    w( 3) = 0.4462908928985841D+00
    w( 4) = 0.4462908928985841D+00
    w( 5) = 0.3267938603769863D+00
    w( 6) = 0.2269152467244296D+00

  else if ( n == 7 ) then

    x( 1) =  -0.9238795325112867D+00
    x( 2) =  -0.7071067811865476D+00
    x( 3) =  -0.3826834323650898D+00
    x( 4) =   0.0000000000000000D+00
    x( 5) =   0.3826834323650898D+00
    x( 6) =   0.7071067811865476D+00
    x( 7) =   0.9238795325112867D+00

    w( 1) = 0.1779646809620499D+00
    w( 2) = 0.2476190476190476D+00
    w( 3) = 0.3934638904665215D+00
    w( 4) = 0.3619047619047619D+00
    w( 5) = 0.3934638904665215D+00
    w( 6) = 0.2476190476190476D+00
    w( 7) = 0.1779646809620499D+00

  else if ( n == 8 ) then

    x( 1) =  -0.9396926207859084D+00
    x( 2) =  -0.7660444431189780D+00
    x( 3) =  -0.5000000000000000D+00
    x( 4) =  -0.1736481776669304D+00
    x( 5) =   0.1736481776669304D+00
    x( 6) =   0.5000000000000000D+00
    x( 7) =   0.7660444431189780D+00
    x( 8) =   0.9396926207859084D+00

    w( 1) = 0.1397697435050225D+00
    w( 2) = 0.2063696457302284D+00
    w( 3) = 0.3142857142857143D+00
    w( 4) = 0.3395748964790348D+00
    w( 5) = 0.3395748964790348D+00
    w( 6) = 0.3142857142857143D+00
    w( 7) = 0.2063696457302284D+00
    w( 8) = 0.1397697435050225D+00

  else if ( n == 9 ) then

    x( 1) =  -0.9510565162951535D+00
    x( 2) =  -0.8090169943749475D+00
    x( 3) =  -0.5877852522924731D+00
    x( 4) =  -0.3090169943749475D+00
    x( 5) =   0.0000000000000000D+00
    x( 6) =   0.3090169943749475D+00
    x( 7) =   0.5877852522924731D+00
    x( 8) =   0.8090169943749475D+00
    x( 9) =   0.9510565162951535D+00

    w( 1) = 0.1147810750857217D+00
    w( 2) = 0.1654331942222276D+00
    w( 3) = 0.2737903534857068D+00
    w( 4) = 0.2790112502222170D+00
    w( 5) = 0.3339682539682539D+00
    w( 6) = 0.2790112502222170D+00
    w( 7) = 0.2737903534857068D+00
    w( 8) = 0.1654331942222276D+00
    w( 9) = 0.1147810750857217D+00

  else if ( n == 10 ) then

    x( 1) =  -0.9594929736144974D+00
    x( 2) =  -0.8412535328311812D+00
    x( 3) =  -0.6548607339452851D+00
    x( 4) =  -0.4154150130018864D+00
    x( 5) =  -0.1423148382732851D+00
    x( 6) =   0.1423148382732851D+00
    x( 7) =   0.4154150130018864D+00
    x( 8) =   0.6548607339452851D+00
    x( 9) =   0.8412535328311812D+00
    x(10) =   0.9594929736144974D+00

    w( 1) = 0.09441954173982806D+00
    w( 2) = 0.1411354380109716D+00
    w( 3) = 0.2263866903636005D+00
    w( 4) = 0.2530509772156453D+00
    w( 5) = 0.2850073526699544D+00
    w( 6) = 0.2850073526699544D+00
    w( 7) = 0.2530509772156453D+00
    w( 8) = 0.2263866903636005D+00
    w( 9) = 0.1411354380109716D+00
    w(10) = 0.09441954173982806D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEJER2_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 through 10.'
    stop

  end if

  return
end
subroutine gegenbauer_compute ( n, alpha, x, w )

!*****************************************************************************80
!
!! GEGENBAUER_COMPUTE computes a Gauss-Gegenbauer quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) (1-x*x)^alpha * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    Thanks to Janiki Raman for pointing out a problem in an earlier
!    version of the code that occurred when ALPHA was -0.5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X*X) in the weight.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) an
  real ( kind = 8 ) bn
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) cc
  real ( kind = 8 ) delta
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval
!
!  Check N.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEGENBAUER_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  1 <= N is required.'
    stop
  end if
!
!  Check ALPHA.
!
  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEGENBAUER_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  -1.0 < ALPHA is required.'
    stop
  end if
!
!  Set the recursion coefficients.
!
  c(1) = 0.0D+00

  if ( 2 <= n ) then
    c(2) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
  end if

  do i = 3, n

    c(i) = real ( i - 1, kind = 8 ) &
          * ( alpha + alpha + real ( i - 1, kind = 8 ) ) / &
          ( ( alpha + alpha + real ( 2 * i - 1, kind = 8 ) ) &
          * ( alpha + alpha + real ( 2 * i - 3, kind = 8 ) ) )

  end do

  delta = r8_gamma ( alpha         + 1.0D+00 ) &
        * r8_gamma (         alpha + 1.0D+00 ) &
        / r8_gamma ( alpha + alpha + 2.0D+00 )

  cc = delta * 2.0D+00**( 2.0D+00 * alpha + 1.0D+00 ) * product ( c(2:n) )

  do i = 1, n

    if ( i == 1 ) then

      an = alpha / real ( n, kind = 8 )

      r1 = ( 1.0D+00 + alpha ) &
        * ( 2.78D+00 / ( 4.0D+00 + real ( n * n, kind = 8 ) ) &
        + 0.768D+00 * an / real ( n, kind = 8 ) )

      r2 = 1.0D+00 + 2.44D+00 * an + 1.282D+00 * an * an

      xval = ( r2 - r1 ) / r2

    else if ( i == 2 ) then

      r1 = ( 4.1D+00 + alpha ) / &
        ( ( 1.0D+00 + alpha ) * ( 1.0D+00 + 0.156D+00 * alpha ) )

      r2 = 1.0D+00 + 0.06D+00 * ( real ( n, kind = 8 ) - 8.0D+00 ) * &
        ( 1.0D+00 + 0.12D+00 * alpha ) / real ( n, kind = 8 )

      r3 = 1.0D+00 + 0.012D+00 * alpha * &
        ( 1.0D+00 + 0.25D+00 * abs ( alpha ) ) / real ( n, kind = 8 )

      xval = xval - r1 * r2 * r3 * ( 1.0D+00 - xval )

    else if ( i == 3 ) then

      r1 = ( 1.67D+00 + 0.28D+00 * alpha ) / ( 1.0D+00 + 0.37D+00 * alpha )

      r2 = 1.0D+00 + 0.22D+00 * ( real ( n, kind = 8 ) - 8.0D+00 ) &
        / real ( n, kind = 8 )

      r3 = 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( n * n, kind = 8 ) )

      xval = xval - r1 * r2 * r3 * ( x(1) - xval )

    else if ( i < n - 1 ) then

      xval = 3.0D+00 * x(i-1) - 3.0D+00 * x(i-2) + x(i-3)

    else if ( i == n - 1 ) then

      r1 = ( 1.0D+00 + 0.235D+00 * alpha ) / ( 0.766D+00 + 0.119D+00 * alpha )

      r2 = 1.0D+00 / ( 1.0D+00 + 0.639D+00 &
        * ( real ( n, kind = 8 ) - 4.0D+00 ) &
        / ( 1.0D+00 + 0.71D+00 * ( real ( n, kind = 8 ) - 4.0D+00 ) ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 20.0D+00 * alpha / ( ( 7.5D+00 + alpha ) * &
        real ( n * n, kind = 8 ) ) )

      xval = xval + r1 * r2 * r3 * ( xval - x(i-2) )

    else if ( i == n ) then

      r1 = ( 1.0D+00 + 0.37D+00 * alpha ) / ( 1.67D+00 + 0.28D+00 * alpha )

      r2 = 1.0D+00 / &
        ( 1.0D+00 + 0.22D+00 * ( real ( n, kind = 8 ) - 8.0D+00 ) &
        / real ( n, kind = 8 ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( n * n, kind = 8 ) ) )

      xval = xval + r1 * r2 * r3 * ( xval - x(i-2) )

    end if

    call gegenbauer_root ( xval, n, alpha, dp2, p1, c )

    x(i) = xval
    w(i) = cc / ( dp2 * p1 )

  end do
!
!  Reverse the data.
!
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

  return
end
subroutine gegenbauer_integral ( expon, alpha, value )

!*****************************************************************************80
!
!! GEGENBAUER_INTEGRAL: integral of a monomial with Gegenbauer weight.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= +1 ) x^expon (1-x*x)^alpha dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X*X) in the weight.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) c
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) value
  real ( kind = 8 ) value1

  if ( mod ( expon, 2 ) == 1 ) then
    value = 0.0D+00
    return
  end if

  c = real ( expon, kind = 8 )

  arg1 = - alpha
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + alpha + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

  value = 2.0D+00 * r8_gamma ( 1.0D+00 + c ) * r8_gamma ( 1.0D+00 + alpha ) &
    * value1 / r8_gamma ( 2.0D+00 + alpha + c )

  return
end
subroutine gegenbauer_recur ( p2, dp2, p1, x, n, alpha, c )

!*****************************************************************************80
!
!! GEGENBAUER_RECUR finds the value and derivative of a Gegenbauer polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of J(N)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of J'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(N-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X*X) in the weight.
!    -1.0 < ALPHA.
!
!    Input, real ( kind = 8 ) C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, n

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = x * p1 - c(i) * p0
    dp2 = x * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine gegenbauer_root ( x, n, alpha, dp2, p1, c )

!*****************************************************************************80
!
!! GEGENBAUER_ROOT improves an approximate root of a Gegenbauer polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X^2) in the weight.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) DP2, the value of J'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(N-1)(X).
!
!    Input, real ( kind = 8 ) C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = epsilon ( eps )

  do step = 1, step_max

    call gegenbauer_recur ( p2, dp2, p1, x, n, alpha, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine gen_hermite_dr_compute ( n, alpha, x, w )

!*****************************************************************************80
!
!! GEN_HERMITE_DR_COMPUTE computes a generalized Gauss-Hermite rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo < x < +oo ) |x|^alpha exp(-x^2) f(x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2011
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
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha_laguerre
  real ( kind = 8 ) arg
  integer ( kind = 4 ) n_laguerre
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ), allocatable, dimension ( : ) :: w_laguerre
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), allocatable, dimension ( : ) :: x_laguerre
!
!  Generate the related generalized Laguerre rule.
!
  if ( n == 1 ) then
    arg = ( alpha + 1.0D+00 ) / 2.0D+00
    x(1) = 0.0D+00
    w(1) = r8_gamma ( arg )
    return
  end if

  if ( mod ( n, 2 ) == 0 ) then
    n_laguerre = n / 2
    alpha_laguerre = ( alpha - 1.0D+00 ) / 2.0D+00
  else if ( mod ( n, 2 ) == 1 ) then
    n_laguerre = ( n - 1 ) / 2
    alpha_laguerre = ( alpha + 1.0D+00 ) / 2.0D+00
  end if

  allocate ( w_laguerre(n_laguerre) )
  allocate ( x_laguerre(n_laguerre) )

  call gen_laguerre_ss_compute ( n_laguerre, alpha_laguerre, x_laguerre, &
    w_laguerre )

  if ( mod ( n, 2 ) == 0 ) then

    x(1:n_laguerre) = - sqrt ( x_laguerre(n_laguerre:1:-1) )
    x(n_laguerre+1:n_laguerre+n_laguerre) &
      = sqrt ( x_laguerre(1:n_laguerre) )
    
    w(1:n_laguerre) = 0.5D+00 * w_laguerre(n_laguerre:1:-1)
    w(n_laguerre+1:n_laguerre+n_laguerre) &
      = 0.5D+00 * w_laguerre(1:n_laguerre)
    
  else if ( mod ( n, 2 ) == 1 ) then

    x(1:n_laguerre) = - sqrt ( x_laguerre(n_laguerre:1:-1) )
    x(n_laguerre+1) = 0.0D+00
    x(n_laguerre+2:n_laguerre+n_laguerre+1) &
      = sqrt ( x_laguerre(1:n_laguerre) )

    w(1:n_laguerre) = 0.5D+00 * w_laguerre(n_laguerre:1:-1) &
                                  / x_laguerre(n_laguerre:1:-1)

    arg = ( alpha + 1.0D+00 ) / 2.0D+00
    w(n_laguerre+1) &
      = r8_gamma ( arg ) &
      - sum ( w_laguerre(1:n_laguerre) / x_laguerre(1:n_laguerre) )

    w(n_laguerre+2:n_laguerre+n_laguerre+1) &
      = 0.5D+00 * w_laguerre(1:n_laguerre) &
                / x_laguerre(1:n_laguerre)

  end if

  deallocate ( w_laguerre )
  deallocate ( x_laguerre )

  return
end
subroutine gen_hermite_ek_compute ( n, alpha, x, w )

!*****************************************************************************80
!
!! GEN_HERMITE_EK_COMPUTE computes a generalized Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The code uses an algorithm by Elhay and Kautsky.
!
!    The integral:
!
!      integral ( -oo < x < +oo ) |x|^alpha exp(-x^2) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of abscissas.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    if ( mod ( i, 2 ) == 1 ) then
      bj(i) = ( i_r8 + alpha ) / 2.0D+00
    else
      bj(i) = i_r8 / 2.0D+00
    end if
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0D+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine gen_hermite_integral ( expon, alpha, value )

!*****************************************************************************80
!
!! GEN_HERMITE_INTEGRAL evaluates a monomial generalized Hermite integral.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo < x < +oo ) x^n |x|^alpha exp(-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent of the monomial.
!    0 <= EXPON.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of |X| in the
!    weight function.  -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) value

  if ( mod ( expon, 2 ) == 1 ) then

    value = 0.0D+00

  else

    a = alpha + real ( expon, kind = 8 )

    if ( a <= -1.0D+00 ) then

      value = - huge ( value )

    else

      value = r8_gamma ( ( a + 1.0D+00 ) / 2.0D+00 )

    end if

  end if

  return
end
subroutine gen_laguerre_ek_compute ( n, alpha, x, w )

!*****************************************************************************80
!
!! GEN_LAGUERRE_EK_COMPUTE: generalized Gauss-Laguerre quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) exp ( - x ) * x^alpha * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^alpha * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * exp ( x(i) ) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = gamma ( alpha + 1.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    bj(i) = i_r8 * ( i_r8 + alpha )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  do i = 1, n
    i_r8 = real ( i, kind = 8 )
    x(i) = 2.0D+00 * i_r8 - 1.0D+00 + alpha
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine gen_laguerre_integral ( expon, alpha, exact )

!*****************************************************************************80
!
!! GEN_LAGUERRE_INTEGRAL evaluates a monomial genearlized Laguerre integral.
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^n * x^alpha * exp ( -x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma

  arg = alpha + real ( expon + 1, kind = 8 )

  exact = r8_gamma ( arg )

  return
end
subroutine gen_laguerre_ss_compute ( n, alpha, x, w )

!*****************************************************************************80
!
!! GEN_LAGUERRE_SS_COMPUTE: generalized Gauss-Laguerre quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) exp ( - x ) * x^alpha * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^alpha * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * exp ( x(i) ) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be at least 1.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) ratio
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval
!
!  Set the recursion coefficients.
!
  do i = 1, n
    b(i) = ( alpha + real ( 2 * i - 1, kind = 8 ) )
  end do

  do i = 1, n
    c(i) = real ( i - 1, kind = 8 ) * ( alpha + real ( i - 1, kind = 8 ) )
  end do

  cc = r8_gamma ( alpha + 1.0D+00 ) * product ( c(2:n) )

  do i = 1, n
!
!  Compute an estimate for the root.
!
    if ( i == 1 ) then

      xval = ( 1.0D+00 + alpha ) * ( 3.0D+00 + 0.92 * alpha ) / &
        ( 1.0D+00 + 2.4D+00 * real ( n, kind = 8 ) + 1.8D+00 * alpha )

    else if ( i == 2 ) then

      xval = xval + ( 15.0D+00 + 6.25D+00 * alpha ) / &
        ( 1.0D+00 + 0.9D+00 * alpha + 2.5D+00 * real ( n, kind = 8 ) )

    else

      r1 = ( 1.0D+00 + 2.55D+00 * real ( i - 2, kind = 8 ) ) &
        / ( 1.9D+00 * real ( i - 2, kind = 8 ) )

      r2 = 1.26D+00 * real ( i - 2, kind = 8 ) * alpha / &
        ( 1.0D+00 + 3.5D+00 * real ( i - 2, kind = 8 ) )

      ratio = ( r1 + r2 ) / ( 1.0D+00 + 0.3D+00 * alpha )

      xval = xval + ratio * ( xval - x(i-2) )

    end if
!
!  Use iteration to find the root.
!
    call gen_laguerre_ss_root ( xval, n, alpha, dp2, p1, b, c )
!
!  Set the abscissa and weight.
!
    x(i) = xval
    w(i) = ( cc / dp2 ) / p1

  end do

  return
end
subroutine gen_laguerre_ss_recur ( p2, dp2, p1, x, n, alpha, b, c )

!*****************************************************************************80
!
!! GEN_LAGUERRE_SS_RECUR evaluates a generalized Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of L(N)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of L'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(N-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor in the
!    integrand.
!
!    Input, real ( kind = 8 ) B(N), C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x - alpha - 1.0D+00
  dp2 = 1.0D+00

  do i = 2, n

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine gen_laguerre_ss_root ( x, n, alpha, dp2, p1, b, c )

!*****************************************************************************80
!
!! GEN_LAGUERRE_SS_ROOT seeks roots of a generalized Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!
!    Output, real ( kind = 8 ) DP2, the value of L'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(N-1)(X).
!
!    Input, real ( kind = 8 ) B(N), C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = epsilon ( eps )

  do step = 1, step_max

    call gen_laguerre_ss_recur ( p2, dp2, p1, x, n, alpha, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine hermite_ek_compute ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_EK_COMPUTE computes a Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The code uses an algorithm by Elhay and Kautsky.
!
!    The abscissas are the zeros of the N-th order Hermite polynomial.
!
!    The integral:
!
!      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2012
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of abscissas.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_gamma ( 1.0D+00 / 2.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = 8 ) / 2.0D+00
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0D+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )
!
!  If N is odd, force the center X to be exactly 0.
!
  if ( mod ( n, 2 ) == 1 ) then
    x((n+1)/2) = 0.0D+00
  end if

  w(1:n) = w(1:n)**2

  return
end
subroutine hermite_gk16_set ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_GK16_SET sets a Hermite Genz-Keister 16 rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+16, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and 35.
!
!    The precisions of these rules are P = 1, 5, 15, 29, and 51.
!
!    Consider, however, the special cases where a rule of precision
!    at least 7, 17, 31 or 33 is desired.  Ordinarily, this would
!    suggest using the nested rule of order 9, 19, 51 or 51 respectively.
!    In these cases, however, the order of the rule that is used exceeds
!    the precision requested.  Hence, it is possible simply to select
!    a subset of the points in the higher precision rule and get a
!    rule of lower order and the desired precision.  This accounts for
!    the four extra rules in this family.
!
!    The entire list of rules is therefore:
!
!    L   P   N
!
!    0   1   1  <-- Full rule 0
!    1   5   3  <-- Full rule 1
!    2   7   7  <-- Partial rule
!    3  15   9  <-- Full rule 2
!    4  17  17  <-- Partial rule
!    5  29  19  <-- Full rule 3
!    6  31  31  <-- Partial rule
!    7  33  33  <-- Partial rule
!    8  51  35  <-- Full rule 4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be 1, 3, 7, 9, 17, 19, 31, 33 or 35.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x( 1) =   0.0000000000000000D+00

    w( 1) =   1.7724538509055159D+00

  else if ( n == 3 ) then

    x( 1) =  -1.2247448713915889D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   1.2247448713915889D+00

    w( 1) =   2.9540897515091930D-01
    w( 2) =   1.1816359006036772D+00
    w( 3) =   2.9540897515091930D-01

  else if ( n == 7 ) then

    x( 1) =  -2.9592107790638380D+00
    x( 2) =  -1.2247448713915889D+00
    x( 3) =  -5.2403354748695763D-01
    x( 4) =   0.0000000000000000D+00
    x( 5) =   5.2403354748695763D-01
    x( 6) =   1.2247448713915889D+00
    x( 7) =   2.9592107790638380D+00

    w( 1) =   1.2330680655153448D-03
    w( 2) =   2.4557928535031393D-01
    w( 3) =   2.3286251787386100D-01
    w( 4) =   8.1310410832613500D-01
    w( 5) =   2.3286251787386100D-01
    w( 6) =   2.4557928535031393D-01
    w( 7) =   1.2330680655153448D-03

  else if ( n == 9 ) then

    x( 1) =  -2.9592107790638380D+00
    x( 2) =  -2.0232301911005157D+00
    x( 3) =  -1.2247448713915889D+00
    x( 4) =  -5.2403354748695763D-01
    x( 5) =   0.0000000000000000D+00
    x( 6) =   5.2403354748695763D-01
    x( 7) =   1.2247448713915889D+00
    x( 8) =   2.0232301911005157D+00
    x( 9) =   2.9592107790638380D+00

    w( 1) =   1.6708826306882348D-04
    w( 2) =   1.4173117873979098D-02
    w( 3) =   1.6811892894767771D-01
    w( 4) =   4.7869428549114124D-01
    w( 5) =   4.5014700975378197D-01
    w( 6) =   4.7869428549114124D-01
    w( 7) =   1.6811892894767771D-01
    w( 8) =   1.4173117873979098D-02
    w( 9) =   1.6708826306882348D-04

  else if ( n == 17 ) then

    x( 1) =  -4.4995993983103881D+00
    x( 2) =  -3.6677742159463378D+00
    x( 3) =  -2.9592107790638380D+00
    x( 4) =  -2.0232301911005157D+00
    x( 5) =  -1.8357079751751868D+00
    x( 6) =  -1.2247448713915889D+00
    x( 7) =  -8.7004089535290285D-01
    x( 8) =  -5.2403354748695763D-01
    x( 9) =   0.0000000000000000D+00
    x(10) =   5.2403354748695763D-01
    x(11) =   8.7004089535290285D-01
    x(12) =   1.2247448713915889D+00
    x(13) =   1.8357079751751868D+00
    x(14) =   2.0232301911005157D+00
    x(15) =   2.9592107790638380D+00
    x(16) =   3.6677742159463378D+00
    x(17) =   4.4995993983103881D+00

    w( 1) =   3.7463469943051758D-08
    w( 2) =  -1.4542843387069391D-06
    w( 3) =   1.8723818949278350D-04
    w( 4) =   1.2466519132805918D-02
    w( 5) =   3.4840719346803800D-03
    w( 6) =   1.5718298376652240D-01
    w( 7) =   2.5155825701712934D-02
    w( 8) =   4.5119803602358544D-01
    w( 9) =   4.7310733504965385D-01
    w(10) =   4.5119803602358544D-01
    w(11) =   2.5155825701712934D-02
    w(12) =   1.5718298376652240D-01
    w(13) =   3.4840719346803800D-03
    w(14) =   1.2466519132805918D-02
    w(15) =   1.8723818949278350D-04
    w(16) =  -1.4542843387069391D-06
    w(17) =   3.7463469943051758D-08

  else if ( n == 19 ) then

    x( 1) =  -4.4995993983103881D+00
    x( 2) =  -3.6677742159463378D+00
    x( 3) =  -2.9592107790638380D+00
    x( 4) =  -2.2665132620567876D+00
    x( 5) =  -2.0232301911005157D+00
    x( 6) =  -1.8357079751751868D+00
    x( 7) =  -1.2247448713915889D+00
    x( 8) =  -8.7004089535290285D-01
    x( 9) =  -5.2403354748695763D-01
    x(10) =   0.0000000000000000D+00
    x(11) =   5.2403354748695763D-01
    x(12) =   8.7004089535290285D-01
    x(13) =   1.2247448713915889D+00
    x(14) =   1.8357079751751868D+00
    x(15) =   2.0232301911005157D+00
    x(16) =   2.2665132620567876D+00
    x(17) =   2.9592107790638380D+00
    x(18) =   3.6677742159463378D+00
    x(19) =   4.4995993983103881D+00

    w( 1) =   1.5295717705322357D-09
    w( 2) =   1.0802767206624762D-06
    w( 3) =   1.0656589772852267D-04
    w( 4) =   5.1133174390883855D-03
    w( 5) = - 1.1232438489069229D-02
    w( 6) =   3.2055243099445879D-02
    w( 7) =   1.1360729895748269D-01
    w( 8) =   1.0838861955003017D-01
    w( 9) =   3.6924643368920851D-01
    w(10) =   5.3788160700510168D-01
    w(11) =   3.6924643368920851D-01
    w(12) =   1.0838861955003017D-01
    w(13) =   1.1360729895748269D-01
    w(14) =   3.2055243099445879D-02
    w(15) = - 1.1232438489069229D-02
    w(16) =   5.1133174390883855D-03
    w(17) =   1.0656589772852267D-04
    w(18) =   1.0802767206624762D-06
    w(19) =   1.5295717705322357D-09

  else if ( n == 31 ) then

    x( 1) =  -6.3759392709822356D+00
    x( 2) =  -5.6432578578857449D+00
    x( 3) =  -5.0360899444730940D+00
    x( 4) =  -4.4995993983103881D+00
    x( 5) =  -3.6677742159463378D+00
    x( 6) =  -2.9592107790638380D+00
    x( 7) =  -2.5705583765842968D+00
    x( 8) =  -2.2665132620567876D+00
    x( 9) =  -2.0232301911005157D+00
    x(10) =  -1.8357079751751868D+00
    x(11) =  -1.5794121348467671D+00
    x(12) =  -1.2247448713915889D+00
    x(13) =  -8.7004089535290285D-01
    x(14) =  -5.2403354748695763D-01
    x(15) =  -1.7606414208200893D-01
    x(16) =   0.0000000000000000D+00
    x(17) =   1.7606414208200893D-01
    x(18) =   5.2403354748695763D-01
    x(19) =   8.7004089535290285D-01
    x(20) =   1.2247448713915889D+00
    x(21) =   1.5794121348467671D+00
    x(22) =   1.8357079751751868D+00
    x(23) =   2.0232301911005157D+00
    x(24) =   2.2665132620567876D+00
    x(25) =   2.5705583765842968D+00
    x(26) =   2.9592107790638380D+00
    x(27) =   3.6677742159463378D+00
    x(28) =   4.4995993983103881D+00
    x(29) =   5.0360899444730940D+00
    x(30) =   5.6432578578857449D+00
    x(31) =   6.3759392709822356D+00

    w( 1) =   2.2365645607044459D-15
    w( 2) =  -2.6304696458548942D-13
    w( 3) =   9.0675288231679823D-12
    w( 4) =   1.4055252024722478D-09
    w( 5) =   1.0889219692128120D-06
    w( 6) =   1.0541662394746661D-04
    w( 7) =   2.6665159778939428D-05
    w( 8) =   4.8385208205502612D-03
    w( 9) =  -9.8566270434610019D-03
    w(10) =   2.9409427580350787D-02
    w(11) =   3.1210210352682834D-03
    w(12) =   1.0939325071860877D-01
    w(13) =   1.1594930984853116D-01
    w(14) =   3.5393889029580544D-01
    w(15) =   4.9855761893293160D-02
    w(16) =   4.5888839636756751D-01
    w(17) =   4.9855761893293160D-02
    w(18) =   3.5393889029580544D-01
    w(19) =   1.1594930984853116D-01
    w(20) =   1.0939325071860877D-01
    w(21) =   3.1210210352682834D-03
    w(22) =   2.9409427580350787D-02
    w(23) =  -9.8566270434610019D-03
    w(24) =   4.8385208205502612D-03
    w(25) =   2.6665159778939428D-05
    w(26) =   1.0541662394746661D-04
    w(27) =   1.0889219692128120D-06
    w(28) =   1.4055252024722478D-09
    w(29) =   9.0675288231679823D-12
    w(30) =  -2.6304696458548942D-13
    w(31) =   2.2365645607044459D-15

  else if ( n == 33 ) then

    x( 1) =  -6.3759392709822356D+00
    x( 2) =  -5.6432578578857449D+00
    x( 3) =  -5.0360899444730940D+00
    x( 4) =  -4.4995993983103881D+00
    x( 5) =  -4.0292201405043713D+00
    x( 6) =  -3.6677742159463378D+00
    x( 7) =  -2.9592107790638380D+00
    x( 8) =  -2.5705583765842968D+00
    x( 9) =  -2.2665132620567876D+00
    x(10) =  -2.0232301911005157D+00
    x(11) =  -1.8357079751751868D+00
    x(12) =  -1.5794121348467671D+00
    x(13) =  -1.2247448713915889D+00
    x(14) =  -8.7004089535290285D-01
    x(15) =  -5.2403354748695763D-01
    x(16) =  -1.7606414208200893D-01
    x(17) =   0.0000000000000000D+00
    x(18) =   1.7606414208200893D-01
    x(19) =   5.2403354748695763D-01
    x(20) =   8.7004089535290285D-01
    x(21) =   1.2247448713915889D+00
    x(22) =   1.5794121348467671D+00
    x(23) =   1.8357079751751868D+00
    x(24) =   2.0232301911005157D+00
    x(25) =   2.2665132620567876D+00
    x(26) =   2.5705583765842968D+00
    x(27) =   2.9592107790638380D+00
    x(28) =   3.6677742159463378D+00
    x(29) =   4.0292201405043713D+00
    x(30) =   4.4995993983103881D+00
    x(31) =   5.0360899444730940D+00
    x(32) =   5.6432578578857449D+00
    x(33) =   6.3759392709822356D+00

    w( 1) =  -1.7602932805372496D-15
    w( 2) =   4.7219278666417693D-13
    w( 3) =  -3.4281570530349562D-11
    w( 4) =   2.7547825138935901D-09
    w( 5) =  -2.3903343382803510D-08
    w( 6) =   1.2245220967158438D-06
    w( 7) =   9.8710009197409173D-05
    w( 8) =   1.4753204901862772D-04
    w( 9) =   3.7580026604304793D-03
    w(10) =  -4.9118576123877555D-03
    w(11) =   2.0435058359107205D-02
    w(12) =   1.3032872699027960D-02
    w(13) =   9.6913444944583621D-02
    w(14) =   1.3726521191567551D-01
    w(15) =   3.1208656194697448D-01
    w(16) =   1.8411696047725790D-01
    w(17) =   2.4656644932829619D-01
    w(18) =   1.8411696047725790D-01
    w(19) =   3.1208656194697448D-01
    w(20) =   1.3726521191567551D-01
    w(21) =   9.6913444944583621D-02
    w(22) =   1.3032872699027960D-02
    w(23) =   2.0435058359107205D-02
    w(24) =  -4.9118576123877555D-03
    w(25) =   3.7580026604304793D-03
    w(26) =   1.4753204901862772D-04
    w(27) =   9.8710009197409173D-05
    w(28) =   1.2245220967158438D-06
    w(29) =  -2.3903343382803510D-08
    w(30) =   2.7547825138935901D-09
    w(31) =  -3.4281570530349562D-11
    w(32) =   4.7219278666417693D-13
    w(33) =  -1.7602932805372496D-15

  else if ( n == 35 ) then

    x( 1) =  -6.3759392709822356D+00
    x( 2) =  -5.6432578578857449D+00
    x( 3) =  -5.0360899444730940D+00
    x( 4) =  -4.4995993983103881D+00
    x( 5) =  -4.0292201405043713D+00
    x( 6) =  -3.6677742159463378D+00
    x( 7) =  -3.3491639537131945D+00
    x( 8) =  -2.9592107790638380D+00
    x( 9) =  -2.5705583765842968D+00
    x(10) =  -2.2665132620567876D+00
    x(11) =  -2.0232301911005157D+00
    x(12) =  -1.8357079751751868D+00
    x(13) =  -1.5794121348467671D+00
    x(14) =  -1.2247448713915889D+00
    x(15) =  -8.7004089535290285D-01
    x(16) =  -5.2403354748695763D-01
    x(17) =  -1.7606414208200893D-01
    x(18) =   0.0000000000000000D+00
    x(19) =   1.7606414208200893D-01
    x(20) =   5.2403354748695763D-01
    x(21) =   8.7004089535290285D-01
    x(22) =   1.2247448713915889D+00
    x(23) =   1.5794121348467671D+00
    x(24) =   1.8357079751751868D+00
    x(25) =   2.0232301911005157D+00
    x(26) =   2.2665132620567876D+00
    x(27) =   2.5705583765842968D+00
    x(28) =   2.9592107790638380D+00
    x(29) =   3.3491639537131945D+00
    x(30) =   3.6677742159463378D+00
    x(31) =   4.0292201405043713D+00
    x(32) =   4.4995993983103881D+00
    x(33) =   5.0360899444730940D+00
    x(34) =   5.6432578578857449D+00
    x(35) =   6.3759392709822356D+00

    w( 1) =   1.8684014894510604D-18
    w( 2) =   9.6599466278563243D-15
    w( 3) =   5.4896836948499462D-12
    w( 4) =   8.1553721816916897D-10
    w( 5) =   3.7920222392319532D-08
    w( 6) =   4.3737818040926989D-07
    w( 7) =   4.8462799737020461D-06
    w( 8) =   6.3328620805617891D-05
    w( 9) =   4.8785399304443770D-04
    w(10) =   1.4515580425155904D-03
    w(11) =   4.0967527720344047D-03
    w(12) =   5.5928828911469180D-03
    w(13) =   2.7780508908535097D-02
    w(14) =   8.0245518147390893D-02
    w(15) =   1.6371221555735804D-01
    w(16) =   2.6244871488784277D-01
    w(17) =   3.3988595585585218D-01
    w(18) =   9.1262675363737921D-04
    w(19) =   3.3988595585585218D-01
    w(20) =   2.6244871488784277D-01
    w(21) =   1.6371221555735804D-01
    w(22) =   8.0245518147390893D-02
    w(23) =   2.7780508908535097D-02
    w(24) =   5.5928828911469180D-03
    w(25) =   4.0967527720344047D-03
    w(26) =   1.4515580425155904D-03
    w(27) =   4.8785399304443770D-04
    w(28) =   6.3328620805617891D-05
    w(29) =   4.8462799737020461D-06
    w(30) =   4.3737818040926989D-07
    w(31) =   3.7920222392319532D-08
    w(32) =   8.1553721816916897D-10
    w(33) =   5.4896836948499462D-12
    w(34) =   9.6599466278563243D-15
    w(35) =   1.8684014894510604D-18

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_GK16_SET - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 7, 9, 17, 19, 31, 33 or 35.'
    stop

  end if

  return
end
subroutine hermite_gk18_set ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_GK18_SET sets a Hermite Genz-Keister 18 rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+18, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and 37.
!
!    The precisions of these rules are P = 1, 5, 15, 29, and 55.
!
!    Some of the data in this function was kindly supplied directly by
!    Alan Genz on 24 April 2011.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be 1, 3, 9, 19, or 37.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x( 1) =   0.0000000000000000D+00

    w( 1) =   1.7724538509055159D+00

  else if ( n == 3 ) then

    x( 1) =  -1.2247448713915889D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   1.2247448713915889D+00

    w( 1) =   2.9540897515091930D-01
    w( 2) =   1.1816359006036772D+00
    w( 3) =   2.9540897515091930D-01

  else if ( n == 9 ) then

    x( 1) =  -2.9592107790638380D+00
    x( 2) =  -2.0232301911005157D+00
    x( 3) =  -1.2247448713915889D+00
    x( 4) =  -5.2403354748695763D-01
    x( 5) =   0.0000000000000000D+00
    x( 6) =   5.2403354748695763D-01
    x( 7) =   1.2247448713915889D+00
    x( 8) =   2.0232301911005157D+00
    x( 9) =   2.9592107790638380D+00

    w( 1) =   1.6708826306882348D-04
    w( 2) =   1.4173117873979098D-02
    w( 3) =   1.6811892894767771D-01
    w( 4) =   4.7869428549114124D-01
    w( 5) =   4.5014700975378197D-01
    w( 6) =   4.7869428549114124D-01
    w( 7) =   1.6811892894767771D-01
    w( 8) =   1.4173117873979098D-02
    w( 9) =   1.6708826306882348D-04

  else if ( n == 19 ) then

    x( 1) =  -4.4995993983103881D+00
    x( 2) =  -3.6677742159463378D+00
    x( 3) =  -2.9592107790638380D+00
    x( 4) =  -2.2665132620567876D+00
    x( 5) =  -2.0232301911005157D+00
    x( 6) =  -1.8357079751751868D+00
    x( 7) =  -1.2247448713915889D+00
    x( 8) =  -8.7004089535290285D-01
    x( 9) =  -5.2403354748695763D-01
    x(10) =   0.0000000000000000D+00
    x(11) =   5.2403354748695763D-01
    x(12) =   8.7004089535290285D-01
    x(13) =   1.2247448713915889D+00
    x(14) =   1.8357079751751868D+00
    x(15) =   2.0232301911005157D+00
    x(16) =   2.2665132620567876D+00
    x(17) =   2.9592107790638380D+00
    x(18) =   3.6677742159463378D+00
    x(19) =   4.4995993983103881D+00

    w( 1) =   1.5295717705322357D-09
    w( 2) =   1.0802767206624762D-06
    w( 3) =   1.0656589772852267D-04
    w( 4) =   5.1133174390883855D-03
    w( 5) =  -1.1232438489069229D-02
    w( 6) =   3.2055243099445879D-02
    w( 7) =   1.1360729895748269D-01
    w( 8) =   1.0838861955003017D-01
    w( 9) =   3.6924643368920851D-01
    w(10) =   5.3788160700510168D-01
    w(11) =   3.6924643368920851D-01
    w(12) =   1.0838861955003017D-01
    w(13) =   1.1360729895748269D-01
    w(14) =   3.2055243099445879D-02
    w(15) =  -1.1232438489069229D-02
    w(16) =   5.1133174390883855D-03
    w(17) =   1.0656589772852267D-04
    w(18) =   1.0802767206624762D-06
    w(19) =   1.5295717705322357D-09

  else if ( n == 37 ) then

    x( 1) =  -6.853200069757519D+00
    x( 2) =  -6.124527854622158D+00
    x( 3) =  -5.521865209868350D+00
    x( 4) =  -4.986551454150765D+00
    x( 5) =  -4.499599398310388D+00
    x( 6) =  -4.057956316089741D+00
    x( 7) =  -3.667774215946338D+00
    x( 8) =  -3.315584617593290D+00
    x( 9) =  -2.959210779063838D+00
    x(10) =  -2.597288631188366D+00
    x(11) =  -2.266513262056788D+00
    x(12) =  -2.023230191100516D+00
    x(13) =  -1.835707975175187D+00
    x(14) =  -1.561553427651873D+00
    x(15) =  -1.224744871391589D+00
    x(16) =  -0.870040895352903D+00
    x(17) =  -0.524033547486958D+00
    x(18) =  -0.214618180588171D+00
    x(19) =   0.000000000000000D+00
    x(20) =   0.214618180588171D+00
    x(21) =   0.524033547486958D+00
    x(22) =   0.870040895352903D+00
    x(23) =   1.224744871391589D+00
    x(24) =   1.561553427651873D+00
    x(25) =   1.835707975175187D+00
    x(26) =   2.023230191100516D+00
    x(27) =   2.266513262056788D+00
    x(28) =   2.597288631188366D+00
    x(29) =   2.959210779063838D+00
    x(30) =   3.315584617593290D+00
    x(31) =   3.667774215946338D+00
    x(32) =   4.057956316089741D+00
    x(33) =   4.499599398310388D+00
    x(34) =   4.986551454150765D+00
    x(35) =   5.521865209868350D+00
    x(36) =   6.124527854622158D+00
    x(37) =   6.853200069757519D+00

    w( 1) = 0.337304188079177058D-20
    w( 2) = 0.332834739632930463D-16
    w( 3) = 0.323016866782871498D-13
    w( 4) = 0.809333688669950037D-11
    w( 5) = 0.748907559239519284D-09
    w( 6) = 0.294146671497083432D-07
    w( 7) = 0.524482423744884136D-06
    w( 8) = 0.586639457073896277D-05
    w( 9) = 0.571885531470621903D-04
    w(10) = 0.41642095727577091D-03
    w(11) = 0.174733389581099482D-02
    w(12) = 0.313373786000304381D-02
    w(13) = 0.768092665770660459D-02
    w(14) = 0.274962713372148476D-01
    w(15) = 0.783630990508037449D-01
    w(16) = 0.16611584261479281D+00
    w(17) = 0.253636910481387185D+00
    w(18) = 0.261712932511430884D+00
    w(19) = 0.171719680968980257D+00
    w(20) = 0.261712932511430884D+00
    w(21) = 0.253636910481387185D+00
    w(22) = 0.16611584261479281D+00
    w(23) = 0.783630990508037449D-01
    w(24) = 0.274962713372148476D-01
    w(25) = 0.768092665770660459D-02
    w(26) = 0.313373786000304381D-02
    w(27) = 0.174733389581099482D-02
    w(28) = 0.41642095727577091D-03
    w(29) = 0.571885531470621903D-04
    w(30) = 0.586639457073896277D-05
    w(31) = 0.524482423744884136D-06
    w(32) = 0.294146671497083432D-07
    w(33) = 0.748907559239519284D-09
    w(34) = 0.809333688669950037D-11
    w(35) = 0.323016866782871498D-13
    w(36) = 0.332834739632930463D-16
    w(37) = 0.337304188079177058D-20

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_GK18_SET - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 9, 19, or 37.'
    stop

  end if

  return
end
subroutine hermite_gk22_set ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_GK22_SET sets a Hermite Genz-Keister 22 rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+22, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and 41.
!
!    The precisions of these rules are P = 1, 5, 15, 29, and 63.
!
!    Some of the data in this function was kindly supplied directly by
!    Alan Genz on 24 April 2011.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 9, 19, and 41.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x( 1) =   0.0000000000000000D+00

    w( 1) =   1.7724538509055159D+00

  else if ( n == 3 ) then

    x( 1) =  -1.2247448713915889D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   1.2247448713915889D+00

    w( 1) =   2.9540897515091930D-01
    w( 2) =   1.1816359006036772D+00
    w( 3) =   2.9540897515091930D-01

  else if ( n == 9 ) then

    x( 1) =  -2.9592107790638380D+00
    x( 2) =  -2.0232301911005157D+00
    x( 3) =  -1.2247448713915889D+00
    x( 4) =  -0.52403354748695763D+00
    x( 5) =   0.0000000000000000D+00
    x( 6) =   0.52403354748695763D+00
    x( 7) =   1.2247448713915889D+00
    x( 8) =   2.0232301911005157D+00
    x( 9) =   2.9592107790638380D+00

    w( 1) =   1.6708826306882348D-04
    w( 2) =   1.4173117873979098D-02
    w( 3) =   1.6811892894767771D-01
    w( 4) =   4.7869428549114124D-01
    w( 5) =   4.5014700975378197D-01
    w( 6) =   4.7869428549114124D-01
    w( 7) =   1.6811892894767771D-01
    w( 8) =   1.4173117873979098D-02
    w( 9) =   1.6708826306882348D-04

  else if ( n == 19 ) then

    x( 1) =  -4.4995993983103881D+00
    x( 2) =  -3.6677742159463378D+00
    x( 3) =  -2.9592107790638380D+00
    x( 4) =  -2.2665132620567876D+00
    x( 5) =  -2.0232301911005157D+00
    x( 6) =  -1.8357079751751868D+00
    x( 7) =  -1.2247448713915889D+00
    x( 8) =  -0.87004089535290285D+00
    x( 9) =  -0.52403354748695763D+00
    x(10) =   0.0000000000000000D+00
    x(11) =   0.52403354748695763D+00
    x(12) =   0.87004089535290285D+00
    x(13) =   1.2247448713915889D+00
    x(14) =   1.8357079751751868D+00
    x(15) =   2.0232301911005157D+00
    x(16) =   2.2665132620567876D+00
    x(17) =   2.9592107790638380D+00
    x(18) =   3.6677742159463378D+00
    x(19) =   4.4995993983103881D+00

    w( 1) =   1.5295717705322357D-09
    w( 2) =   1.0802767206624762D-06
    w( 3) =   1.0656589772852267D-04
    w( 4) =   5.1133174390883855D-03
    w( 5) =  -1.1232438489069229D-02
    w( 6) =   3.2055243099445879D-02
    w( 7) =   1.1360729895748269D-01
    w( 8) =   1.0838861955003017D-01
    w( 9) =   3.6924643368920851D-01
    w(10) =   5.3788160700510168D-01
    w(11) =   3.6924643368920851D-01
    w(12) =   1.0838861955003017D-01
    w(13) =   1.1360729895748269D-01
    w(14) =   3.2055243099445879D-02
    w(15) =  -1.1232438489069229D-02
    w(16) =   5.1133174390883855D-03
    w(17) =   1.0656589772852267D-04
    w(18) =   1.0802767206624762D-06
    w(19) =   1.5295717705322357D-09

  else if ( n == 41 ) then

    x( 1) =  -7.251792998192644D+00
    x( 2) =  -6.547083258397540D+00
    x( 3) =  -5.961461043404500D+00
    x( 4) =  -5.437443360177798D+00
    x( 5) =  -4.953574342912980D+00
    x( 6) =  -4.4995993983103881D+00
    x( 7) =  -4.070919267883068D+00
    x( 8) =  -3.6677742159463378D+00
    x( 9) =  -3.296114596212218D+00
    x(10) =  -2.9592107790638380D+00
    x(11) =  -2.630415236459871D+00
    x(12) =  -2.2665132620567876D+00
    x(13) =  -2.043834754429505D+00
    x(14) =  -2.0232301911005157D+00
    x(15) =  -1.8357079751751868D+00
    x(16) =  -1.585873011819188D+00
    x(17) =  -1.2247448713915889D+00
    x(18) =  -0.87004089535290285D+00
    x(19) =  -0.52403354748695763D+00
    x(20) =  -0.195324784415805D+00
    x(21) =   0.0000000000000000D+00
    x(22) =   0.195324784415805D+00
    x(23) =   0.52403354748695763D+00
    x(24) =   0.87004089535290285D+00
    x(25) =   1.2247448713915889D+00
    x(26) =   1.585873011819188D+00
    x(27) =   1.8357079751751868D+00
    x(28) =   2.0232301911005157D+00
    x(29) =   2.043834754429505D+00
    x(30) =   2.2665132620567876D+00
    x(31) =   2.630415236459871D+00
    x(32) =   2.9592107790638380D+00
    x(33) =   3.296114596212218D+00
    x(34) =   3.6677742159463378D+00
    x(35) =   4.070919267883068D+00
    x(36) =   4.4995993983103881D+00
    x(37) =   4.953574342912980D+00
    x(38) =   5.437443360177798D+00
    x(39) =   5.961461043404500D+00
    x(40) =   6.547083258397540D+00
    x(41) =   7.251792998192644D+00

    w( 1) =   0.117725656974405367D-22
    w( 2) =   0.152506745534300636D-18
    w( 3) =   0.202183949965101288D-15
    w( 4) =   0.724614869051195508D-13
    w( 5) =   0.103121966469463034D-10
    w( 6) =   0.710371395169350952D-09
    w( 7) =   0.264376044449260516D-07
    w( 8) =   0.558982787078644997D-06
    w( 9) =   0.675628907134744976D-05
    w(10) =   0.512198007019776873D-04
    w(11) =   0.335013114947200879D-03
    w(12) =   0.249379691096933139D-02
    w(13) = - 0.25616995850607458D-01
    w(14) =   0.317007878644325588D-01
    w(15) =   0.125041498584003435D-02
    w(16) =   0.293244560924894295D-01
    w(17) =   0.799536390803302298D-01
    w(18) =   0.164543666806555251D+00
    w(19) =   0.258718519718241095D+00
    w(20) =   0.293588795735908566D+00
    w(21) =   0.997525375254611951D-01
    w(22) =   0.293588795735908566D+00
    w(23) =   0.258718519718241095D+00
    w(24) =   0.164543666806555251D+00
    w(25) =   0.799536390803302298D-01
    w(26) =   0.293244560924894295D-01
    w(27) =   0.125041498584003435D-02
    w(28) =   0.317007878644325588D-01
    w(29) = - 0.25616995850607458D-01
    w(30) =   0.249379691096933139D-02
    w(31) =   0.335013114947200879D-03
    w(32) =   0.512198007019776873D-04
    w(33) =   0.675628907134744976D-05
    w(34) =   0.558982787078644997D-06
    w(35) =   0.264376044449260516D-07
    w(36) =   0.710371395169350952D-09
    w(37) =   0.103121966469463034D-10
    w(38) =   0.724614869051195508D-13
    w(39) =   0.202183949965101288D-15
    w(40) =   0.152506745534300636D-18
    w(41) =   0.117725656974405367D-22

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_GK22_SET - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 9, 19, or 41.'
    stop

  end if

  return
end
subroutine hermite_gk24_set ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_GK24_SET sets a Hermite Genz-Keister 24 rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A nested family of rules for the Hermite integration problem
!    was produced by Genz and Keister.  The structure of the nested
!    family was denoted by 1+2+6+10+24, that is, it comprised rules
!    of successive orders O = 1, 3, 9, 19, and 43.
!
!    The precisions of these rules are P = 1, 5, 15, 29, and 67.
!
!    Some of the data in this function was kindly supplied directly by
!    Alan Genz on 24 April 2011.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    Legal values are 1, 3, 9, 19, and 43.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x( 1) =   0.0000000000000000D+00

    w( 1) =   1.7724538509055159D+00

  else if ( n == 3 ) then

    x( 1) =  -1.2247448713915889D+00
    x( 2) =   0.0000000000000000D+00
    x( 3) =   1.2247448713915889D+00

    w( 1) =   2.9540897515091930D-01
    w( 2) =   1.1816359006036772D+00
    w( 3) =   2.9540897515091930D-01

  else if ( n == 9 ) then

    x( 1) =  -2.9592107790638380D+00
    x( 2) =  -2.0232301911005157D+00
    x( 3) =  -1.2247448713915889D+00
    x( 4) =  -0.52403354748695763D+00
    x( 5) =   0.0000000000000000D+00
    x( 6) =   0.52403354748695763D+00
    x( 7) =   1.2247448713915889D+00
    x( 8) =   2.0232301911005157D+00
    x( 9) =   2.9592107790638380D+00

    w( 1) =   1.6708826306882348D-04
    w( 2) =   1.4173117873979098D-02
    w( 3) =   1.6811892894767771D-01
    w( 4) =   4.7869428549114124D-01
    w( 5) =   4.5014700975378197D-01
    w( 6) =   4.7869428549114124D-01
    w( 7) =   1.6811892894767771D-01
    w( 8) =   1.4173117873979098D-02
    w( 9) =   1.6708826306882348D-04

  else if ( n == 19 ) then

    x( 1) =  -4.4995993983103881D+00
    x( 2) =  -3.6677742159463378D+00
    x( 3) =  -2.9592107790638380D+00
    x( 4) =  -2.2665132620567876D+00
    x( 5) =  -2.0232301911005157D+00
    x( 6) =  -1.8357079751751868D+00
    x( 7) =  -1.2247448713915889D+00
    x( 8) =  -0.87004089535290285D+00
    x( 9) =  -0.52403354748695763D+00
    x(10) =   0.0000000000000000D+00
    x(11) =   0.52403354748695763D+00
    x(12) =   0.87004089535290285D+00
    x(13) =   1.2247448713915889D+00
    x(14) =   1.8357079751751868D+00
    x(15) =   2.0232301911005157D+00
    x(16) =   2.2665132620567876D+00
    x(17) =   2.9592107790638380D+00
    x(18) =   3.6677742159463378D+00
    x(19) =   4.4995993983103881D+00

    w( 1) =   1.5295717705322357D-09
    w( 2) =   1.0802767206624762D-06
    w( 3) =   1.0656589772852267D-04
    w( 4) =   5.1133174390883855D-03
    w( 5) =  -1.1232438489069229D-02
    w( 6) =   3.2055243099445879D-02
    w( 7) =   1.1360729895748269D-01
    w( 8) =   1.0838861955003017D-01
    w( 9) =   3.6924643368920851D-01
    w(10) =   5.3788160700510168D-01
    w(11) =   3.6924643368920851D-01
    w(12) =   1.0838861955003017D-01
    w(13) =   1.1360729895748269D-01
    w(14) =   3.2055243099445879D-02
    w(15) =  -1.1232438489069229D-02
    w(16) =   5.1133174390883855D-03
    w(17) =   1.0656589772852267D-04
    w(18) =   1.0802767206624762D-06
    w(19) =   1.5295717705322357D-09

  else if ( n == 43 ) then

    x( 1) = -10.167574994881873D+00
    x( 2) =  -7.231746029072501D+00
    x( 3) =  -6.535398426382995D+00
    x( 4) =  -5.954781975039809D+00
    x( 5) =  -5.434053000365068D+00
    x( 6) =  -4.952329763008589D+00
    x( 7) =  -4.4995993983103881D+00
    x( 8) =  -4.071335874253583D+00
    x( 9) =  -3.6677742159463378D+00
    x(10) =  -3.295265921534226D+00
    x(11) =  -2.9592107790638380D+00
    x(12) =  -2.633356763661946D+00
    x(13) =  -2.2665132620567876D+00
    x(14) =  -2.089340389294661D+00
    x(15) =  -2.0232301911005157D+00
    x(16) =  -1.8357079751751868D+00
    x(17) =  -1.583643465293944D+00
    x(18) =  -1.2247448713915889D+00
    x(19) =  -0.87004089535290285D+00
    x(20) =  -0.52403354748695763D+00
    x(21) =  -0.196029453662011D+00
    x(22) =   0.0000000000000000D+00
    x(23) =   0.196029453662011D+00
    x(24) =   0.52403354748695763D+00
    x(25) =   0.87004089535290285D+00
    x(26) =   1.2247448713915889D+00
    x(27) =   1.583643465293944D+00
    x(28) =   1.8357079751751868D+00
    x(29) =   2.0232301911005157D+00
    x(30) =   2.089340389294661D+00
    x(31) =   2.2665132620567876D+00
    x(32) =   2.633356763661946D+00
    x(33) =   2.9592107790638380D+00
    x(34) =   3.295265921534226D+00
    x(35) =   3.6677742159463378D+00
    x(36) =   4.071335874253583D+00
    x(37) =   4.4995993983103881D+00
    x(38) =   4.952329763008589D+00
    x(39) =   5.434053000365068D+00
    x(40) =   5.954781975039809D+00
    x(41) =   6.535398426382995D+00
    x(42) =   7.231746029072501D+00
    x(43) =  10.167574994881873D+00

    w( 1) =   0.968100020641528185D-37
    w( 2) =   0.15516931262860431D-22
    w( 3) =   0.175937309107750992D-18
    w( 4) =   0.217337608710893738D-15
    w( 5) =   0.747837010380540069D-13
    w( 6) =   0.104028132097205732D-10
    w( 7) =   0.70903573389336778D-09
    w( 8) =   0.263481722999966618D-07
    w( 9) =   0.560127964848432175D-06
    w(10) =   0.680410934802210232D-05
    w(11) =   0.508343873102544037D-04
    w(12) =   0.32753080006610181D-03
    w(13) =   0.267479828788552937D-02
    w(14) = - 0.687704270963253854D-02
    w(15) =   0.119383201790913588D-01
    w(16) =   0.248083722871002796D-02
    w(17) =   0.29000335749726387D-01
    w(18) =   0.798689557875757008D-01
    w(19) =   0.164609842422580606D+00
    w(20) =   0.258535954731607738D+00
    w(21) =   0.292243810406117141D+00
    w(22) =   0.102730713753441829D+00
    w(23) =   0.292243810406117141D+00
    w(24) =   0.258535954731607738D+00
    w(25) =   0.164609842422580606D+00
    w(26) =   0.798689557875757008D-01
    w(27) =   0.29000335749726387D-01
    w(28) =   0.248083722871002796D-02
    w(29) =   0.119383201790913588D-01
    w(30) = - 0.687704270963253854D-02
    w(31) =   0.267479828788552937D-02
    w(32) =   0.32753080006610181D-03
    w(33) =   0.508343873102544037D-04
    w(34) =   0.680410934802210232D-05
    w(35) =   0.560127964848432175D-06
    w(36) =   0.263481722999966618D-07
    w(37) =   0.70903573389336778D-09
    w(38) =   0.104028132097205732D-10
    w(39) =   0.747837010380540069D-13
    w(40) =   0.217337608710893738D-15
    w(41) =   0.175937309107750992D-18
    w(42) =   0.15516931262860431D-22
    w(43) =   0.968100020641528185D-37
 
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_GK24_SET - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 9, 19, or 43.'
    stop

  end if

  return
end
subroutine hermite_integral ( n, value )

!*****************************************************************************80
!
!! HERMITE_INTEGRAL evaluates a monomial Hermite integral.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo < x < +oo ) x^n exp(-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the integral.
!    0 <= N.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_factorial2
  real ( kind = 8 ) value

  if ( n < 0 ) then

    value = - huge ( value )

  else if ( mod ( n, 2 ) == 1 ) then

    value = 0.0D+00

  else

    value = r8_factorial2 ( n - 1 ) * sqrt ( pi ) / 2.0D+00**( n / 2 )

  end if

  return
end
subroutine hermite_set ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_SET sets abscissas and weights for Hermite quadrature.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    Mathematica can numerically estimate the abscissas of the rule
!    of order N to P digits by the command:
!
!      NSolve [ HermiteH [ n, x ] == 0, x, p ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    Dover, 2006,
!    ISBN: 0486445798,
!    LC: QA311.K713.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 20, 31/32/33, 63/64/65, 127/128/129.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = 0.0D+00

    w(1) = 1.77245385090551602729816748334D+00

  else if ( n == 2 ) then

    x(1) = - 0.707106781186547524400844362105D+00
    x(2) =   0.707106781186547524400844362105D+00

    w(1) = 0.886226925452758013649083741671D+00
    w(2) = 0.886226925452758013649083741671D+00

  else if ( n == 3 ) then

    x(1) = - 0.122474487139158904909864203735D+01
    x(2) =   0.0D+00
    x(3) =   0.122474487139158904909864203735D+01

    w(1) = 0.295408975150919337883027913890D+00
    w(2) = 0.118163590060367735153211165556D+01
    w(3) = 0.295408975150919337883027913890D+00

  else if ( n == 4 ) then

    x(1) = - 0.165068012388578455588334111112D+01
    x(2) = - 0.524647623275290317884060253835D+00
    x(3) =   0.524647623275290317884060253835D+00
    x(4) =   0.165068012388578455588334111112D+01

    w(1) = 0.813128354472451771430345571899D-01
    w(2) = 0.804914090005512836506049184481D+00
    w(3) = 0.804914090005512836506049184481D+00
    w(4) = 0.813128354472451771430345571899D-01

  else if ( n == 5 ) then

    x(1) = - 0.202018287045608563292872408814D+01
    x(2) = - 0.958572464613818507112770593893D+00
    x(3) =   0.0D+00
    x(4) =   0.958572464613818507112770593893D+00
    x(5) =   0.202018287045608563292872408814D+01

    w(1) = 0.199532420590459132077434585942D-01
    w(2) = 0.393619323152241159828495620852D+00
    w(3) = 0.945308720482941881225689324449D+00
    w(4) = 0.393619323152241159828495620852D+00
    w(5) = 0.199532420590459132077434585942D-01

  else if ( n == 6 ) then

    x(1) = - 0.235060497367449222283392198706D+01
    x(2) = - 0.133584907401369694971489528297D+01
    x(3) = - 0.436077411927616508679215948251D+00
    x(4) =   0.436077411927616508679215948251D+00
    x(5) =   0.133584907401369694971489528297D+01
    x(6) =   0.235060497367449222283392198706D+01

    w(1) = 0.453000990550884564085747256463D-02
    w(2) = 0.157067320322856643916311563508D+00
    w(3) = 0.724629595224392524091914705598D+00
    w(4) = 0.724629595224392524091914705598D+00
    w(5) = 0.157067320322856643916311563508D+00
    w(6) = 0.453000990550884564085747256463D-02

  else if ( n == 7 ) then

    x(1) = - 0.265196135683523349244708200652D+01
    x(2) = - 0.167355162876747144503180139830D+01
    x(3) = - 0.816287882858964663038710959027D+00
    x(4) =   0.0D+00
    x(5) =   0.816287882858964663038710959027D+00
    x(6) =   0.167355162876747144503180139830D+01
    x(7) =   0.265196135683523349244708200652D+01

    w(1) = 0.971781245099519154149424255939D-03
    w(2) = 0.545155828191270305921785688417D-01
    w(3) = 0.425607252610127800520317466666D+00
    w(4) = 0.810264617556807326764876563813D+00
    w(5) = 0.425607252610127800520317466666D+00
    w(6) = 0.545155828191270305921785688417D-01
    w(7) = 0.971781245099519154149424255939D-03

  else if ( n == 8 ) then

    x(1) = - 0.293063742025724401922350270524D+01
    x(2) = - 0.198165675669584292585463063977D+01
    x(3) = - 0.115719371244678019472076577906D+01
    x(4) = - 0.381186990207322116854718885584D+00
    x(5) =   0.381186990207322116854718885584D+00
    x(6) =   0.115719371244678019472076577906D+01
    x(7) =   0.198165675669584292585463063977D+01
    x(8) =   0.293063742025724401922350270524D+01

    w(1) = 0.199604072211367619206090452544D-03
    w(2) = 0.170779830074134754562030564364D-01
    w(3) = 0.207802325814891879543258620286D+00
    w(4) = 0.661147012558241291030415974496D+00
    w(5) = 0.661147012558241291030415974496D+00
    w(6) = 0.207802325814891879543258620286D+00
    w(7) = 0.170779830074134754562030564364D-01
    w(8) = 0.199604072211367619206090452544D-03

  else if ( n == 9 ) then

    x(1) = - 0.319099320178152760723004779538D+01
    x(2) = - 0.226658058453184311180209693284D+01
    x(3) = - 0.146855328921666793166701573925D+01
    x(4) = - 0.723551018752837573322639864579D+00
    x(5) =   0.0D+00
    x(6) =   0.723551018752837573322639864579D+00
    x(7) =   0.146855328921666793166701573925D+01
    x(8) =   0.226658058453184311180209693284D+01
    x(9) =   0.319099320178152760723004779538D+01

    w(1) = 0.396069772632643819045862946425D-04
    w(2) = 0.494362427553694721722456597763D-02
    w(3) = 0.884745273943765732879751147476D-01
    w(4) = 0.432651559002555750199812112956D+00
    w(5) = 0.720235215606050957124334723389D+00
    w(6) = 0.432651559002555750199812112956D+00
    w(7) = 0.884745273943765732879751147476D-01
    w(8) = 0.494362427553694721722456597763D-02
    w(9) = 0.396069772632643819045862946425D-04

  else if ( n == 10 ) then

    x(1) =  - 0.343615911883773760332672549432D+01
    x(2) =  - 0.253273167423278979640896079775D+01
    x(3) =  - 0.175668364929988177345140122011D+01
    x(4) =  - 0.103661082978951365417749191676D+01
    x(5) =  - 0.342901327223704608789165025557D+00
    x(6) =    0.342901327223704608789165025557D+00
    x(7) =    0.103661082978951365417749191676D+01
    x(8) =    0.175668364929988177345140122011D+01
    x(9) =    0.253273167423278979640896079775D+01
    x(10) =   0.343615911883773760332672549432D+01

    w(1) =  0.764043285523262062915936785960D-05
    w(2) =  0.134364574678123269220156558585D-02
    w(3) =  0.338743944554810631361647312776D-01
    w(4) =  0.240138611082314686416523295006D+00
    w(5) =  0.610862633735325798783564990433D+00
    w(6) =  0.610862633735325798783564990433D+00
    w(7) =  0.240138611082314686416523295006D+00
    w(8) =  0.338743944554810631361647312776D-01
    w(9) =  0.134364574678123269220156558585D-02
    w(10) = 0.764043285523262062915936785960D-05

  else if ( n == 11 ) then

    x(1) =  - 0.366847084655958251845837146485D+01
    x(2) =  - 0.278329009978165177083671870152D+01
    x(3) =  - 0.202594801582575533516591283121D+01
    x(4) =  - 0.132655708449493285594973473558D+01
    x(5) =  - 0.656809566882099765024611575383D+00
    x(6) =    0.0D+00
    x(7) =    0.656809566882099765024611575383D+00
    x(8) =    0.132655708449493285594973473558D+01
    x(9) =    0.202594801582575533516591283121D+01
    x(10) =   0.278329009978165177083671870152D+01
    x(11) =   0.366847084655958251845837146485D+01

    w(1) =  0.143956039371425822033088366032D-05
    w(2) =  0.346819466323345510643413772940D-03
    w(3) =  0.119113954449115324503874202916D-01
    w(4) =  0.117227875167708503381788649308D+00
    w(5) =  0.429359752356125028446073598601D+00
    w(6) =  0.654759286914591779203940657627D+00
    w(7) =  0.429359752356125028446073598601D+00
    w(8) =  0.117227875167708503381788649308D+00
    w(9) =  0.119113954449115324503874202916D-01
    w(10) = 0.346819466323345510643413772940D-03
    w(11) = 0.143956039371425822033088366032D-05

  else if ( n == 12 ) then

    x(1) =  - 0.388972489786978191927164274724D+01
    x(2) =  - 0.302063702512088977171067937518D+01
    x(3) =  - 0.227950708050105990018772856942D+01
    x(4) =  - 0.159768263515260479670966277090D+01
    x(5) =  - 0.947788391240163743704578131060D+00
    x(6) =  - 0.314240376254359111276611634095D+00
    x(7) =    0.314240376254359111276611634095D+00
    x(8) =    0.947788391240163743704578131060D+00
    x(9) =    0.159768263515260479670966277090D+01
    x(10) =   0.227950708050105990018772856942D+01
    x(11) =   0.302063702512088977171067937518D+01
    x(12) =   0.388972489786978191927164274724D+01

    w(1) =  0.265855168435630160602311400877D-06
    w(2) =  0.857368704358785865456906323153D-04
    w(3) =  0.390539058462906185999438432620D-02
    w(4) =  0.516079856158839299918734423606D-01
    w(5) =  0.260492310264161129233396139765D+00
    w(6) =  0.570135236262479578347113482275D+00
    w(7) =  0.570135236262479578347113482275D+00
    w(8) =  0.260492310264161129233396139765D+00
    w(9) =  0.516079856158839299918734423606D-01
    w(10) = 0.390539058462906185999438432620D-02
    w(11) = 0.857368704358785865456906323153D-04
    w(12) = 0.265855168435630160602311400877D-06

  else if ( n == 13 ) then

    x(1) =  - 0.410133759617863964117891508007D+01
    x(2) =  - 0.324660897837240998812205115236D+01
    x(3) =  - 0.251973568567823788343040913628D+01
    x(4) =  - 0.185310765160151214200350644316D+01
    x(5) =  - 0.122005503659074842622205526637D+01
    x(6) =  - 0.605763879171060113080537108602D+00
    x(7) =    0.0D+00
    x(8) =    0.605763879171060113080537108602D+00
    x(9) =    0.122005503659074842622205526637D+01
    x(10) =   0.185310765160151214200350644316D+01
    x(11) =   0.251973568567823788343040913628D+01
    x(12) =   0.324660897837240998812205115236D+01
    x(13) =   0.410133759617863964117891508007D+01

    w(1) =  0.482573185007313108834997332342D-07
    w(2) =  0.204303604027070731248669432937D-04
    w(3) =  0.120745999271938594730924899224D-02
    w(4) =  0.208627752961699392166033805050D-01
    w(5) =  0.140323320687023437762792268873D+00
    w(6) =  0.421616296898543221746893558568D+00
    w(7) =  0.604393187921161642342099068579D+00
    w(8) =  0.421616296898543221746893558568D+00
    w(9) =  0.140323320687023437762792268873D+00
    w(10) = 0.208627752961699392166033805050D-01
    w(11) = 0.120745999271938594730924899224D-02
    w(12) = 0.204303604027070731248669432937D-04
    w(13) = 0.482573185007313108834997332342D-07

  else if ( n == 14 ) then

    x(1) =  - 0.430444857047363181262129810037D+01
    x(2) =  - 0.346265693360227055020891736115D+01
    x(3) =  - 0.274847072498540256862499852415D+01
    x(4) =  - 0.209518325850771681573497272630D+01
    x(5) =  - 0.147668273114114087058350654421D+01
    x(6) =  - 0.878713787329399416114679311861D+00
    x(7) =  - 0.291745510672562078446113075799D+00
    x(8) =    0.291745510672562078446113075799D+00
    x(9) =    0.878713787329399416114679311861D+00
    x(10) =   0.147668273114114087058350654421D+01
    x(11) =   0.209518325850771681573497272630D+01
    x(12) =   0.274847072498540256862499852415D+01
    x(13) =   0.346265693360227055020891736115D+01
    x(14) =   0.430444857047363181262129810037D+01

    w(1) =  0.862859116812515794532041783429D-08
    w(2) =  0.471648435501891674887688950105D-05
    w(3) =  0.355092613551923610483661076691D-03
    w(4) =  0.785005472645794431048644334608D-02
    w(5) =  0.685055342234652055387163312367D-01
    w(6) =  0.273105609064246603352569187026D+00
    w(7) =  0.536405909712090149794921296776D+00
    w(8) =  0.536405909712090149794921296776D+00
    w(9) =  0.273105609064246603352569187026D+00
    w(10) = 0.685055342234652055387163312367D-01
    w(11) = 0.785005472645794431048644334608D-02
    w(12) = 0.355092613551923610483661076691D-03
    w(13) = 0.471648435501891674887688950105D-05
    w(14) = 0.862859116812515794532041783429D-08

  else if ( n == 15 ) then

    x(1) =  - 0.449999070730939155366438053053D+01
    x(2) =  - 0.366995037340445253472922383312D+01
    x(3) =  - 0.296716692790560324848896036355D+01
    x(4) =  - 0.232573248617385774545404479449D+01
    x(5) =  - 0.171999257518648893241583152515D+01
    x(6) =  - 0.113611558521092066631913490556D+01
    x(7) =  - 0.565069583255575748526020337198D+00
    x(8) =    0.0D+00
    x(9) =    0.565069583255575748526020337198D+00
    x(10) =   0.113611558521092066631913490556D+01
    x(11) =   0.171999257518648893241583152515D+01
    x(12) =   0.232573248617385774545404479449D+01
    x(13) =   0.296716692790560324848896036355D+01
    x(14) =   0.366995037340445253472922383312D+01
    x(15) =   0.449999070730939155366438053053D+01

    w(1) =  0.152247580425351702016062666965D-08
    w(2) =  0.105911554771106663577520791055D-05
    w(3) =  0.100004441232499868127296736177D-03
    w(4) =  0.277806884291277589607887049229D-02
    w(5) =  0.307800338725460822286814158758D-01
    w(6) =  0.158488915795935746883839384960D+00
    w(7) =  0.412028687498898627025891079568D+00
    w(8) =  0.564100308726417532852625797340D+00
    w(9) =  0.412028687498898627025891079568D+00
    w(10) = 0.158488915795935746883839384960D+00
    w(11) = 0.307800338725460822286814158758D-01
    w(12) = 0.277806884291277589607887049229D-02
    w(13) = 0.100004441232499868127296736177D-03
    w(14) = 0.105911554771106663577520791055D-05
    w(15) = 0.152247580425351702016062666965D-08

  else if ( n == 16 ) then

    x(1) =  - 0.468873893930581836468849864875D+01
    x(2) =  - 0.386944790486012269871942409801D+01
    x(3) =  - 0.317699916197995602681399455926D+01
    x(4) =  - 0.254620215784748136215932870545D+01
    x(5) =  - 0.195178799091625397743465541496D+01
    x(6) =  - 0.138025853919888079637208966969D+01
    x(7) =  - 0.822951449144655892582454496734D+00
    x(8) =  - 0.273481046138152452158280401965D+00
    x(9) =    0.273481046138152452158280401965D+00
    x(10) =   0.822951449144655892582454496734D+00
    x(11) =   0.138025853919888079637208966969D+01
    x(12) =   0.195178799091625397743465541496D+01
    x(13) =   0.254620215784748136215932870545D+01
    x(14) =   0.317699916197995602681399455926D+01
    x(15) =   0.386944790486012269871942409801D+01
    x(16) =   0.468873893930581836468849864875D+01

    w(1) =  0.265480747401118224470926366050D-09
    w(2) =  0.232098084486521065338749423185D-06
    w(3) =  0.271186009253788151201891432244D-04
    w(4) =  0.932284008624180529914277305537D-03
    w(5) =  0.128803115355099736834642999312D-01
    w(6) =  0.838100413989858294154207349001D-01
    w(7) =  0.280647458528533675369463335380D+00
    w(8) =  0.507929479016613741913517341791D+00
    w(9) =  0.507929479016613741913517341791D+00
    w(10) = 0.280647458528533675369463335380D+00
    w(11) = 0.838100413989858294154207349001D-01
    w(12) = 0.128803115355099736834642999312D-01
    w(13) = 0.932284008624180529914277305537D-03
    w(14) = 0.271186009253788151201891432244D-04
    w(15) = 0.232098084486521065338749423185D-06
    w(16) = 0.265480747401118224470926366050D-09

  else if ( n == 17 ) then

    x(1) =  - 0.487134519367440308834927655662D+01
    x(2) =  - 0.406194667587547430689245559698D+01
    x(3) =  - 0.337893209114149408338327069289D+01
    x(4) =  - 0.275776291570388873092640349574D+01
    x(5) =  - 0.217350282666662081927537907149D+01
    x(6) =  - 0.161292431422123133311288254454D+01
    x(7) =  - 0.106764872574345055363045773799D+01
    x(8) =  - 0.531633001342654731349086553718D+00
    x(9) =    0.0D+00
    x(10) =   0.531633001342654731349086553718D+00
    x(11) =   0.106764872574345055363045773799D+01
    x(12) =   0.161292431422123133311288254454D+01
    x(13) =   0.217350282666662081927537907149D+01
    x(14) =   0.275776291570388873092640349574D+01
    x(15) =   0.337893209114149408338327069289D+01
    x(16) =   0.406194667587547430689245559698D+01
    x(17) =   0.487134519367440308834927655662D+01

    w(1) =  0.458057893079863330580889281222D-10
    w(2) =  0.497707898163079405227863353715D-07
    w(3) =  0.711228914002130958353327376218D-05
    w(4) =  0.298643286697753041151336643059D-03
    w(5) =  0.506734995762753791170069495879D-02
    w(6) =  0.409200341495762798094994877854D-01
    w(7) =  0.172648297670097079217645196219D+00
    w(8) =  0.401826469470411956577635085257D+00
    w(9) =  0.530917937624863560331883103379D+00
    w(10) = 0.401826469470411956577635085257D+00
    w(11) = 0.172648297670097079217645196219D+00
    w(12) = 0.409200341495762798094994877854D-01
    w(13) = 0.506734995762753791170069495879D-02
    w(14) = 0.298643286697753041151336643059D-03
    w(15) = 0.711228914002130958353327376218D-05
    w(16) = 0.497707898163079405227863353715D-07
    w(17) = 0.458057893079863330580889281222D-10

  else if ( n == 18 ) then

    x(1) =  - 0.504836400887446676837203757885D+01
    x(2) =  - 0.424811787356812646302342016090D+01
    x(3) =  - 0.357376906848626607950067599377D+01
    x(4) =  - 0.296137750553160684477863254906D+01
    x(5) =  - 0.238629908916668600026459301424D+01
    x(6) =  - 0.183553160426162889225383944409D+01
    x(7) =  - 0.130092085838961736566626555439D+01
    x(8) =  - 0.776682919267411661316659462284D+00
    x(9) =  - 0.258267750519096759258116098711D+00
    x(10) =   0.258267750519096759258116098711D+00
    x(11) =   0.776682919267411661316659462284D+00
    x(12) =   0.130092085838961736566626555439D+01
    x(13) =   0.183553160426162889225383944409D+01
    x(14) =   0.238629908916668600026459301424D+01
    x(15) =   0.296137750553160684477863254906D+01
    x(16) =   0.357376906848626607950067599377D+01
    x(17) =   0.424811787356812646302342016090D+01
    x(18) =   0.504836400887446676837203757885D+01

    w(1) =  0.782819977211589102925147471012D-11
    w(2) =  0.104672057957920824443559608435D-07
    w(3) =  0.181065448109343040959702385911D-05
    w(4) =  0.918112686792940352914675407371D-04
    w(5) =  0.188852263026841789438175325426D-02
    w(6) =  0.186400423875446519219315221973D-01
    w(7) =  0.973017476413154293308537234155D-01
    w(8) =  0.284807285669979578595606820713D+00
    w(9) =  0.483495694725455552876410522141D+00
    w(10) = 0.483495694725455552876410522141D+00
    w(11) = 0.284807285669979578595606820713D+00
    w(12) = 0.973017476413154293308537234155D-01
    w(13) = 0.186400423875446519219315221973D-01
    w(14) = 0.188852263026841789438175325426D-02
    w(15) = 0.918112686792940352914675407371D-04
    w(16) = 0.181065448109343040959702385911D-05
    w(17) = 0.104672057957920824443559608435D-07
    w(18) = 0.782819977211589102925147471012D-11

  else if ( n == 19 ) then

    x(1) =  - 0.522027169053748216460967142500D+01
    x(2) =  - 0.442853280660377943723498532226D+01
    x(3) =  - 0.376218735196402009751489394104D+01
    x(4) =  - 0.315784881834760228184318034120D+01
    x(5) =  - 0.259113378979454256492128084112D+01
    x(6) =  - 0.204923170985061937575050838669D+01
    x(7) =  - 0.152417061939353303183354859367D+01
    x(8) =  - 0.101036838713431135136859873726D+01
    x(9) =  - 0.503520163423888209373811765050D+00
    x(10) =   0.0D+00
    x(11) =   0.503520163423888209373811765050D+00
    x(12) =   0.101036838713431135136859873726D+01
    x(13) =   0.152417061939353303183354859367D+01
    x(14) =   0.204923170985061937575050838669D+01
    x(15) =   0.259113378979454256492128084112D+01
    x(16) =   0.315784881834760228184318034120D+01
    x(17) =   0.376218735196402009751489394104D+01
    x(18) =   0.442853280660377943723498532226D+01
    x(19) =   0.522027169053748216460967142500D+01

    w(1) =  0.132629709449851575185289154385D-11
    w(2) =  0.216305100986355475019693077221D-08
    w(3) =  0.448824314722312295179447915594D-06
    w(4) =  0.272091977631616257711941025214D-04
    w(5) =  0.670877521407181106194696282100D-03
    w(6) =  0.798886677772299020922211491861D-02
    w(7) =  0.508103869090520673569908110358D-01
    w(8) =  0.183632701306997074156148485766D+00
    w(9) =  0.391608988613030244504042313621D+00
    w(10) = 0.502974888276186530840731361096D+00
    w(11) = 0.391608988613030244504042313621D+00
    w(12) = 0.183632701306997074156148485766D+00
    w(13) = 0.508103869090520673569908110358D-01
    w(14) = 0.798886677772299020922211491861D-02
    w(15) = 0.670877521407181106194696282100D-03
    w(16) = 0.272091977631616257711941025214D-04
    w(17) = 0.448824314722312295179447915594D-06
    w(18) = 0.216305100986355475019693077221D-08
    w(19) = 0.132629709449851575185289154385D-11

  else if ( n == 20 ) then

    x(1) =  - 0.538748089001123286201690041068D+01
    x(2) =  - 0.460368244955074427307767524898D+01
    x(3) =  - 0.394476404011562521037562880052D+01
    x(4) =  - 0.334785456738321632691492452300D+01
    x(5) =  - 0.278880605842813048052503375640D+01
    x(6) =  - 0.225497400208927552308233334473D+01
    x(7) =  - 0.173853771211658620678086566214D+01
    x(8) =  - 0.123407621539532300788581834696D+01
    x(9) =  - 0.737473728545394358705605144252D+00
    x(10) = - 0.245340708300901249903836530634D+00
    x(11) =   0.245340708300901249903836530634D+00
    x(12) =   0.737473728545394358705605144252D+00
    x(13) =   0.123407621539532300788581834696D+01
    x(14) =   0.173853771211658620678086566214D+01
    x(15) =   0.225497400208927552308233334473D+01
    x(16) =   0.278880605842813048052503375640D+01
    x(17) =   0.334785456738321632691492452300D+01
    x(18) =   0.394476404011562521037562880052D+01
    x(19) =   0.460368244955074427307767524898D+01
    x(20) =   0.538748089001123286201690041068D+01

    w(1) =  0.222939364553415129252250061603D-12
    w(2) =  0.439934099227318055362885145547D-09
    w(3) =  0.108606937076928169399952456345D-06
    w(4) =  0.780255647853206369414599199965D-05
    w(5) =  0.228338636016353967257145917963D-03
    w(6) =  0.324377334223786183218324713235D-02
    w(7) =  0.248105208874636108821649525589D-01
    w(8) =  0.109017206020023320013755033535D+00
    w(9) =  0.286675505362834129719659706228D+00
    w(10) = 0.462243669600610089650328639861D+00
    w(11) = 0.462243669600610089650328639861D+00
    w(12) = 0.286675505362834129719659706228D+00
    w(13) = 0.109017206020023320013755033535D+00
    w(14) = 0.248105208874636108821649525589D-01
    w(15) = 0.324377334223786183218324713235D-02
    w(16) = 0.228338636016353967257145917963D-03
    w(17) = 0.780255647853206369414599199965D-05
    w(18) = 0.108606937076928169399952456345D-06
    w(19) = 0.439934099227318055362885145547D-09
    w(20) = 0.222939364553415129252250061603D-12

  else if ( n == 31 ) then

    x(1) = -6.99568012371854027532485214732D+00
    x(2) = -6.27507870494286014270365678125D+00
    x(3) = -5.67396144461858832963325587893D+00
    x(4) = -5.13359557711238070458629689140D+00
    x(5) = -4.63155950631285994206679976543D+00
    x(6) = -4.15627175581814517248313523153D+00
    x(7) = -3.70074340323146942244971645897D+00
    x(8) = -3.26032073231354081046454015096D+00
    x(9) = -2.83168045339020545570156401514D+00
    x(10) = -2.41231770548042010517401845821D+00
    x(11) = -2.00025854893563896579755625986D+00
    x(12) = -1.59388586047213982613884194556D+00
    x(13) = -1.19182699835004642608213586492D+00
    x(14) = -0.792876976915308939685930329988D+00
    x(15) = -0.395942736471423110946700416634D+00
    x(16) = 0.0D+00
    x(17) = 0.395942736471423110946700416634D+00
    x(18) = 0.792876976915308939685930329988D+00
    x(19) = 1.19182699835004642608213586492D+00
    x(20) = 1.59388586047213982613884194556D+00
    x(21) = 2.00025854893563896579755625986D+00
    x(22) = 2.41231770548042010517401845821D+00
    x(23) = 2.83168045339020545570156401514D+00
    x(24) = 3.26032073231354081046454015096D+00
    x(25) = 3.70074340323146942244971645897D+00
    x(26) = 4.15627175581814517248313523153D+00
    x(27) = 4.63155950631285994206679976543D+00
    x(28) = 5.13359557711238070458629689140D+00
    x(29) = 5.67396144461858832963325587893D+00
    x(30) = 6.27507870494286014270365678125D+00
    x(31) = 6.99568012371854027532485214732D+00

    w(1) = 4.61896839446420502132944426974D-22
    w(2) = 5.11060900792715640739422641166D-18
    w(3) = 5.89955649875387299038431589378D-15
    w(4) = 1.86037352145214652437380892603D-12
    w(5) = 2.35249200320864163398597795323D-10
    w(6) = 1.46119883449105307352780323055D-08
    w(7) = 5.04371255893979974253745671633D-07
    w(8) = 0.0000104986027576756063228123279208D+00
    w(9) = 0.000139520903950470433823653754396D+00
    w(10) = 0.00123368330730688826551750402319D+00
    w(11) = 0.00748279991403519848345678003016D+00
    w(12) = 0.0318472307313003327772087235339D+00
    w(13) = 0.0967179481608704535580338478886D+00
    w(14) = 0.212132788668764779877735637343D+00
    w(15) = 0.338772657894107724675701919878D+00
    w(16) = 0.395778556098609545141783810611D+00
    w(17) = 0.338772657894107724675701919878D+00
    w(18) = 0.212132788668764779877735637343D+00
    w(19) = 0.0967179481608704535580338478886D+00
    w(20) = 0.0318472307313003327772087235339D+00
    w(21) = 0.00748279991403519848345678003016D+00
    w(22) = 0.00123368330730688826551750402319D+00
    w(23) = 0.000139520903950470433823653754396D+00
    w(24) = 0.0000104986027576756063228123279208D+00
    w(25) = 5.04371255893979974253745671633D-07
    w(26) = 1.46119883449105307352780323055D-08
    w(27) = 2.35249200320864163398597795323D-10
    w(28) = 1.86037352145214652437380892603D-12
    w(29) = 5.89955649875387299038431589378D-15
    w(30) = 5.11060900792715640739422641166D-18
    w(31) = 4.61896839446420502132944426974D-22

  else if ( n == 32 ) then

    x(1) = -7.12581390983072757279520760342D+00
    x(2) = -6.40949814926966041217376374153D+00
    x(3) = -5.81222594951591383276596615366D+00
    x(4) = -5.27555098651588012781906048140D+00
    x(5) = -4.77716450350259639303579405689D+00
    x(6) = -4.30554795335119844526348653193D+00
    x(7) = -3.85375548547144464388787292109D+00
    x(8) = -3.41716749281857073587392729564D+00
    x(9) = -2.99249082500237420628549407606D+00
    x(10) = -2.57724953773231745403092930114D+00
    x(11) = -2.16949918360611217330570559502D+00
    x(12) = -1.76765410946320160462767325853D+00
    x(13) = -1.37037641095287183816170564864D+00
    x(14) = -0.976500463589682838484704871982D+00
    x(15) = -0.584978765435932448466957544011D+00
    x(16) = -0.194840741569399326708741289532D+00
    x(17) = 0.194840741569399326708741289532D+00
    x(18) = 0.584978765435932448466957544011D+00
    x(19) = 0.976500463589682838484704871982D+00
    x(20) = 1.37037641095287183816170564864D+00
    x(21) = 1.76765410946320160462767325853D+00
    x(22) = 2.16949918360611217330570559502D+00
    x(23) = 2.57724953773231745403092930114D+00
    x(24) = 2.99249082500237420628549407606D+00
    x(25) = 3.41716749281857073587392729564D+00
    x(26) = 3.85375548547144464388787292109D+00
    x(27) = 4.30554795335119844526348653193D+00
    x(28) = 4.77716450350259639303579405689D+00
    x(29) = 5.27555098651588012781906048140D+00
    x(30) = 5.81222594951591383276596615366D+00
    x(31) = 6.40949814926966041217376374153D+00
    x(32) = 7.12581390983072757279520760342D+00

    w(1) = 7.31067642738416239327427845506D-23
    w(2) = 9.23173653651829223349442007207D-19
    w(3) = 1.19734401709284866582868189951D-15
    w(4) = 4.21501021132644757296944521183D-13
    w(5) = 5.93329146339663861451156821558D-11
    w(6) = 4.0988321647708966182350410138D-09
    w(7) = 1.57416779254559402926869257929D-07
    w(8) = 3.65058512956237605737032418746D-06
    w(9) = 0.0000541658406181998255800193939267D+00
    w(10) = 0.000536268365527972045970238101501D+00
    w(11) = 0.00365489032665442807912565712241D+00
    w(12) = 0.017553428831573430303437844611D+00
    w(13) = 0.0604581309559126141865857607833D+00
    w(14) = 0.151269734076642482575147114614D+00
    w(15) = 0.277458142302529898137698918542D+00
    w(16) = 0.375238352592802392866818388907D+00
    w(17) = 0.375238352592802392866818388907D+00
    w(18) = 0.277458142302529898137698918542D+00
    w(19) = 0.151269734076642482575147114614D+00
    w(20) = 0.0604581309559126141865857607833D+00
    w(21) = 0.017553428831573430303437844611D+00
    w(22) = 0.00365489032665442807912565712241D+00
    w(23) = 0.000536268365527972045970238101501D+00
    w(24) = 0.0000541658406181998255800193939267D+00
    w(25) = 3.65058512956237605737032418746D-06
    w(26) = 1.57416779254559402926869257929D-07
    w(27) = 4.0988321647708966182350410138D-09
    w(28) = 5.93329146339663861451156821558D-11
    w(29) = 4.21501021132644757296944521183D-13
    w(30) = 1.19734401709284866582868189951D-15
    w(31) = 9.23173653651829223349442007207D-19
    w(32) = 7.31067642738416239327427845506D-23

  else if ( n == 33 ) then

    x(1) = -7.25385182201520064607977244465D+00
    x(2) = -6.54165544573807726095826608811D+00
    x(3) = -5.94807118208714447981366477584D+00
    x(4) = -5.41492900261419253992709076454D+00
    x(5) = -4.92002852059500829241139910265D+00
    x(6) = -4.45191114883282719009473876206D+00
    x(7) = -4.00367160995693141451378357174D+00
    x(8) = -3.57072198023271828561890330658D+00
    x(9) = -3.14979668170382538461281438786D+00
    x(10) = -2.73844582435135490694887052899D+00
    x(11) = -2.33475115152951517708536069773D+00
    x(12) = -1.93715458182220661643452220908D+00
    x(13) = -1.54434826124312180914304288754D+00
    x(14) = -1.15520020412678961356412482063D+00
    x(15) = -0.768701379758868598107224561306D+00
    x(16) = -0.383926014508409083771145488401D+00
    x(17) = 0.0D+00
    x(18) = 0.383926014508409083771145488401D+00
    x(19) = 0.768701379758868598107224561306D+00
    x(20) = 1.15520020412678961356412482063D+00
    x(21) = 1.54434826124312180914304288754D+00
    x(22) = 1.93715458182220661643452220908D+00
    x(23) = 2.33475115152951517708536069773D+00
    x(24) = 2.73844582435135490694887052899D+00
    x(25) = 3.14979668170382538461281438786D+00
    x(26) = 3.57072198023271828561890330658D+00
    x(27) = 4.00367160995693141451378357174D+00
    x(28) = 4.45191114883282719009473876206D+00
    x(29) = 4.92002852059500829241139910265D+00
    x(30) = 5.41492900261419253992709076454D+00
    x(31) = 5.94807118208714447981366477584D+00
    x(32) = 6.54165544573807726095826608811D+00
    x(33) = 7.25385182201520064607977244465D+00

    w(1) = 1.15331621854588454082208890757D-23
    w(2) = 1.65709474153369051048226040291D-19
    w(3) = 2.40778567955799442824587707068D-16
    w(4) = 9.43481415901497503451931586527D-14
    w(5) = 1.47398093709248867676655543441D-11
    w(6) = 1.12892224710833129085848357165D-09
    w(7) = 4.8077456763231909801575826594D-08
    w(8) = 1.23769336720121013593677278301D-06
    w(9) = 0.000020423684051423773240532989618D+00
    w(10) = 0.000225442770596327415479556999963D+00
    w(11) = 0.00171845463776092445684897958755D+00
    w(12) = 0.00926568997068523330372581072023D+00
    w(13) = 0.0359879823185769744486866448437D+00
    w(14) = 0.102069079995541500792505520921D+00
    w(15) = 0.213493931150291836488258605194D+00
    w(16) = 0.331552000750741282288352789764D+00
    w(17) = 0.383785266519863801349608543622D+00
    w(18) = 0.331552000750741282288352789764D+00
    w(19) = 0.213493931150291836488258605194D+00
    w(20) = 0.102069079995541500792505520921D+00
    w(21) = 0.0359879823185769744486866448437D+00
    w(22) = 0.00926568997068523330372581072023D+00
    w(23) = 0.00171845463776092445684897958755D+00
    w(24) = 0.000225442770596327415479556999963D+00
    w(25) = 0.000020423684051423773240532989618D+00
    w(26) = 1.23769336720121013593677278301D-06
    w(27) = 4.8077456763231909801575826594D-08
    w(28) = 1.12892224710833129085848357165D-09
    w(29) = 1.47398093709248867676655543441D-11
    w(30) = 9.43481415901497503451931586527D-14
    w(31) = 2.40778567955799442824587707068D-16
    w(32) = 1.65709474153369051048226040291D-19
    w(33) = 1.15331621854588454082208890757D-23

  else if ( n == 63 ) then

    x(1) = -10.4354998778541680534681154273D+00
    x(2) = -9.80287599129749636352239352865D+00
    x(3) = -9.27920195430503913194047455065D+00
    x(4) = -8.81185814372845464425266282756D+00
    x(5) = -8.38076834518632193430106510438D+00
    x(6) = -7.97559508014203731815418062985D+00
    x(7) = -7.59013951986410667624797831945D+00
    x(8) = -7.22031670788896784611613242225D+00
    x(9) = -6.86325443317953685273532858761D+00
    x(10) = -6.51683481068211606052733958540D+00
    x(11) = -6.17943799227059698624184617873D+00
    x(12) = -5.84978840008106734625265829615D+00
    x(13) = -5.52685725264030314250475751228D+00
    x(14) = -5.20979798304083548615751364163D+00
    x(15) = -4.89790186449757423507450992149D+00
    x(16) = -4.59056657444351902292712945691D+00
    x(17) = -4.28727333528244040317276161995D+00
    x(18) = -3.98756991041971574852270520681D+00
    x(19) = -3.69105770009634651173228105598D+00
    x(20) = -3.39738177133039118527559418063D+00
    x(21) = -3.10622302792825663291386167460D+00
    x(22) = -2.81729196728379777507471356574D+00
    x(23) = -2.53032363047120109268552217185D+00
    x(24) = -2.24507346048120662989959181793D+00
    x(25) = -1.96131385830814852939220084113D+00
    x(26) = -1.67883127917201375208028006226D+00
    x(27) = -1.39742374860496251075707520637D+00
    x(28) = -1.11689870509964626905109702778D+00
    x(29) = -0.837071095589476159777377954613D+00
    x(30) = -0.557761664279082216687636652538D+00
    x(31) = -0.278795385671152239866876286272D+00
    x(32) = 0.0D+00
    x(33) = 0.278795385671152239866876286272D+00
    x(34) = 0.557761664279082216687636652538D+00
    x(35) = 0.837071095589476159777377954613D+00
    x(36) = 1.11689870509964626905109702778D+00
    x(37) = 1.39742374860496251075707520637D+00
    x(38) = 1.67883127917201375208028006226D+00
    x(39) = 1.96131385830814852939220084113D+00
    x(40) = 2.24507346048120662989959181793D+00
    x(41) = 2.53032363047120109268552217185D+00
    x(42) = 2.81729196728379777507471356574D+00
    x(43) = 3.10622302792825663291386167460D+00
    x(44) = 3.39738177133039118527559418063D+00
    x(45) = 3.69105770009634651173228105598D+00
    x(46) = 3.98756991041971574852270520681D+00
    x(47) = 4.28727333528244040317276161995D+00
    x(48) = 4.59056657444351902292712945691D+00
    x(49) = 4.89790186449757423507450992149D+00
    x(50) = 5.20979798304083548615751364163D+00
    x(51) = 5.52685725264030314250475751228D+00
    x(52) = 5.84978840008106734625265829615D+00
    x(53) = 6.17943799227059698624184617873D+00
    x(54) = 6.51683481068211606052733958540D+00
    x(55) = 6.86325443317953685273532858761D+00
    x(56) = 7.22031670788896784611613242225D+00
    x(57) = 7.59013951986410667624797831945D+00
    x(58) = 7.97559508014203731815418062985D+00
    x(59) = 8.38076834518632193430106510438D+00
    x(60) = 8.81185814372845464425266282756D+00
    x(61) = 9.27920195430503913194047455065D+00
    x(62) = 9.80287599129749636352239352865D+00
    x(63) = 10.4354998778541680534681154273D+00

    w(1) = 3.70992064349030055823376157823D-48
    w(2) = 1.04007786152246672212559599908D-42
    w(3) = 1.97968047083199197900260998813D-38
    w(4) = 8.46874781919035663281042885251D-35
    w(5) = 1.30713059308206243904769877879D-31
    w(6) = 9.34378371756582396450246862195D-29
    w(7) = 3.60274266352851638202340658522D-26
    w(8) = 8.29638631162099766157527065317D-24
    w(9) = 1.22666299091434557721622529775D-21
    w(10) = 1.22884356288353036990240371039D-19
    w(11) = 8.69255369584585252225619256428D-18
    w(12) = 4.48570586893158184069444097978D-16
    w(13) = 1.73358179557891044383064226749D-14
    w(14) = 5.1265062385197846998384009333D-13
    w(15) = 1.18089218445696923817995132237D-11
    w(16) = 2.15086982978749617679069862879D-10
    w(17) = 3.13719295353830786449435629291D-09
    w(18) = 3.70416259848969809883356560995D-08
    w(19) = 3.57347329499908777461505032558D-07
    w(20) = 2.83931144984692884712301165567D-06
    w(21) = 0.0000187091130037887216027832755405D+00
    w(22) = 0.000102848808006856425543062213642D+00
    w(23) = 0.000474117026103206754395975199216D+00
    w(24) = 0.0018409222622442103760124297917D+00
    w(25) = 0.00604360445513757113209247151533D+00
    w(26) = 0.0168292991996521044559098701555D+00
    w(27) = 0.0398582640278170328649908688578D+00
    w(28) = 0.0804670879942008323850873860195D+00
    w(29) = 0.138719508176584635072239096351D+00
    w(30) = 0.204486953468973988225911656103D+00
    w(31) = 0.25799889943138332612723393346D+00
    w(32) = 0.278766948849251654365527505911D+00
    w(33) = 0.25799889943138332612723393346D+00
    w(34) = 0.204486953468973988225911656103D+00
    w(35) = 0.138719508176584635072239096351D+00
    w(36) = 0.0804670879942008323850873860195D+00
    w(37) = 0.0398582640278170328649908688578D+00
    w(38) = 0.0168292991996521044559098701555D+00
    w(39) = 0.00604360445513757113209247151533D+00
    w(40) = 0.0018409222622442103760124297917D+00
    w(41) = 0.000474117026103206754395975199216D+00
    w(42) = 0.000102848808006856425543062213642D+00
    w(43) = 0.0000187091130037887216027832755405D+00
    w(44) = 2.83931144984692884712301165567D-06
    w(45) = 3.57347329499908777461505032558D-07
    w(46) = 3.70416259848969809883356560995D-08
    w(47) = 3.13719295353830786449435629291D-09
    w(48) = 2.15086982978749617679069862879D-10
    w(49) = 1.18089218445696923817995132237D-11
    w(50) = 5.1265062385197846998384009333D-13
    w(51) = 1.73358179557891044383064226749D-14
    w(52) = 4.48570586893158184069444097978D-16
    w(53) = 8.69255369584585252225619256428D-18
    w(54) = 1.22884356288353036990240371039D-19
    w(55) = 1.22666299091434557721622529775D-21
    w(56) = 8.29638631162099766157527065317D-24
    w(57) = 3.60274266352851638202340658522D-26
    w(58) = 9.34378371756582396450246862195D-29
    w(59) = 1.30713059308206243904769877879D-31
    w(60) = 8.46874781919035663281042885251D-35
    w(61) = 1.97968047083199197900260998813D-38
    w(62) = 1.04007786152246672212559599908D-42
    w(63) = 3.70992064349030055823376157823D-48

  else if ( n == 64 ) then

    x(1) = -10.5261231679605458833268262838D+00
    x(2) = -9.89528758682953902120446147716D+00
    x(3) = -9.37315954964672116254565243972D+00
    x(4) = -8.90724909996476975729597288564D+00
    x(5) = -8.47752908337986309056416634482D+00
    x(6) = -8.07368728501022522585879114076D+00
    x(7) = -7.68954016404049682844780422987D+00
    x(8) = -7.32101303278094920118956936372D+00
    x(9) = -6.96524112055110752924264219349D+00
    x(10) = -6.62011226263602737903666010894D+00
    x(11) = -6.28401122877482823541809319507D+00
    x(12) = -5.95566632679948604534456718098D+00
    x(13) = -5.63405216434997214724992048331D+00
    x(14) = -5.31832522463327085732364951520D+00
    x(15) = -5.00777960219876819644370262718D+00
    x(16) = -4.70181564740749981609753801581D+00
    x(17) = -4.39991716822813764776793253544D+00
    x(18) = -4.10163447456665671497098123846D+00
    x(19) = -3.80657151394536046116597200046D+00
    x(20) = -3.51437593574090621153995058647D+00
    x(21) = -3.22473129199203572584817111019D+00
    x(22) = -2.93735082300462180968533902619D+00
    x(23) = -2.65197243543063501100545778600D+00
    x(24) = -2.36835458863240140411151126534D+00
    x(25) = -2.08627287988176202083256330236D+00
    x(26) = -1.80551717146554491890377357419D+00
    x(27) = -1.52588914020986366294897013315D+00
    x(28) = -1.24720015694311794069356453069D+00
    x(29) = -0.969269423071178016743541489019D+00
    x(30) = -0.691922305810044577268219287596D+00
    x(31) = -0.414988824121078684576929129200D+00
    x(32) = -0.138302244987009724115049767967D+00
    x(33) = 0.138302244987009724115049767967D+00
    x(34) = 0.414988824121078684576929129200D+00
    x(35) = 0.691922305810044577268219287596D+00
    x(36) = 0.969269423071178016743541489019D+00
    x(37) = 1.24720015694311794069356453069D+00
    x(38) = 1.52588914020986366294897013315D+00
    x(39) = 1.80551717146554491890377357419D+00
    x(40) = 2.08627287988176202083256330236D+00
    x(41) = 2.36835458863240140411151126534D+00
    x(42) = 2.65197243543063501100545778600D+00
    x(43) = 2.93735082300462180968533902619D+00
    x(44) = 3.22473129199203572584817111019D+00
    x(45) = 3.51437593574090621153995058647D+00
    x(46) = 3.80657151394536046116597200046D+00
    x(47) = 4.10163447456665671497098123846D+00
    x(48) = 4.39991716822813764776793253544D+00
    x(49) = 4.70181564740749981609753801581D+00
    x(50) = 5.00777960219876819644370262718D+00
    x(51) = 5.31832522463327085732364951520D+00
    x(52) = 5.63405216434997214724992048331D+00
    x(53) = 5.95566632679948604534456718098D+00
    x(54) = 6.28401122877482823541809319507D+00
    x(55) = 6.62011226263602737903666010894D+00
    x(56) = 6.96524112055110752924264219349D+00
    x(57) = 7.32101303278094920118956936372D+00
    x(58) = 7.68954016404049682844780422987D+00
    x(59) = 8.07368728501022522585879114076D+00
    x(60) = 8.47752908337986309056416634482D+00
    x(61) = 8.90724909996476975729597288564D+00
    x(62) = 9.37315954964672116254565243972D+00
    x(63) = 9.89528758682953902120446147716D+00
    x(64) = 10.5261231679605458833268262838D+00

    w(1) = 5.53570653585694282057546330099D-49
    w(2) = 1.67974799010815921866628833063D-43
    w(3) = 3.42113801125574050432722182815D-39
    w(4) = 1.55739062462976380230933538026D-35
    w(5) = 2.54966089911299925660476658044D-32
    w(6) = 1.92910359546496685030196877907D-29
    w(7) = 7.86179778892591036909999149628D-27
    w(8) = 1.91170688330064282995845696553D-24
    w(9) = 2.98286278427985115447870070202D-22
    w(10) = 3.15225456650378141612134668341D-20
    w(11) = 2.35188471067581911695767591556D-18
    w(12) = 1.28009339132243804163956329526D-16
    w(13) = 5.21862372659084752295780851305D-15
    w(14) = 1.62834073070972036208430708124D-13
    w(15) = 3.95917776694772392723644586425D-12
    w(16) = 7.61521725014545135331529567532D-11
    w(17) = 1.17361674232154934354250646708D-09
    w(18) = 1.4651253164761093549266220038D-08
    w(19) = 1.49553293672724706110246169293D-07
    w(20) = 1.25834025103118457615784218002D-06
    w(21) = 8.7884992308503591814440474067D-06
    w(22) = 0.0000512592913578627466082191141274D+00
    w(23) = 0.000250983698513062486082362017982D+00
    w(24) = 0.00103632909950757766345674174628D+00
    w(25) = 0.00362258697853445876066812537162D+00
    w(26) = 0.0107560405098791370494651727867D+00
    w(27) = 0.0272031289536889184538348212615D+00
    w(28) = 0.0587399819640994345496889462518D+00
    w(29) = 0.108498349306186840633025845506D+00
    w(30) = 0.171685842349083702000727970124D+00
    w(31) = 0.232994786062678046650566029333D+00
    w(32) = 0.271377424941303977945606508418D+00
    w(33) = 0.271377424941303977945606508418D+00
    w(34) = 0.232994786062678046650566029333D+00
    w(35) = 0.171685842349083702000727970124D+00
    w(36) = 0.108498349306186840633025845506D+00
    w(37) = 0.0587399819640994345496889462518D+00
    w(38) = 0.0272031289536889184538348212615D+00
    w(39) = 0.0107560405098791370494651727867D+00
    w(40) = 0.00362258697853445876066812537162D+00
    w(41) = 0.00103632909950757766345674174628D+00
    w(42) = 0.000250983698513062486082362017982D+00
    w(43) = 0.0000512592913578627466082191141274D+00
    w(44) = 8.7884992308503591814440474067D-06
    w(45) = 1.25834025103118457615784218002D-06
    w(46) = 1.49553293672724706110246169293D-07
    w(47) = 1.4651253164761093549266220038D-08
    w(48) = 1.17361674232154934354250646708D-09
    w(49) = 7.61521725014545135331529567532D-11
    w(50) = 3.95917776694772392723644586425D-12
    w(51) = 1.62834073070972036208430708124D-13
    w(52) = 5.21862372659084752295780851305D-15
    w(53) = 1.28009339132243804163956329526D-16
    w(54) = 2.35188471067581911695767591556D-18
    w(55) = 3.15225456650378141612134668341D-20
    w(56) = 2.98286278427985115447870070202D-22
    w(57) = 1.91170688330064282995845696553D-24
    w(58) = 7.86179778892591036909999149628D-27
    w(59) = 1.92910359546496685030196877907D-29
    w(60) = 2.54966089911299925660476658044D-32
    w(61) = 1.55739062462976380230933538026D-35
    w(62) = 3.42113801125574050432722182815D-39
    w(63) = 1.67974799010815921866628833063D-43
    w(64) = 5.53570653585694282057546330099D-49

  else if ( n == 65 ) then

    x(1) = -10.6160229818782811890575602362D+00
    x(2) = -9.98694169167668475289516351866D+00
    x(3) = -9.46632932015538456209230219583D+00
    x(4) = -9.00182332295913301957373193078D+00
    x(5) = -8.57344474441790920512961590649D+00
    x(6) = -8.17090617805258532129873946832D+00
    x(7) = -7.78803908298957078257446963177D+00
    x(8) = -7.42077883436632423648540822493D+00
    x(9) = -7.06626794030689283695626648500D+00
    x(10) = -6.72239982016573443713269171234D+00
    x(11) = -6.38756373978709108766264932629D+00
    x(12) = -6.06049177883150521431844958114D+00
    x(13) = -5.74016182369022552144785845004D+00
    x(14) = -5.42573329769734933377879520007D+00
    x(15) = -5.11650300472141247403405938903D+00
    x(16) = -4.81187385202746476469375701958D+00
    x(17) = -4.51133211136821338520053781440D+00
    x(18) = -4.21443050997195460766227433230D+00
    x(19) = -3.92077540444472380734866143333D+00
    x(20) = -3.63001687763289533380743256995D+00
    x(21) = -3.34184096844683005348576537026D+00
    x(22) = -3.05596348432867100471992904351D+00
    x(23) = -2.77212500515709167382401243242D+00
    x(24) = -2.49008679530393550665251461372D+00
    x(25) = -2.20962741516918436363972079745D+00
    x(26) = -1.93053987597722550177796442890D+00
    x(27) = -1.65262921904032541019598694267D+00
    x(28) = -1.37571042772366833613088581439D+00
    x(29) = -1.09960660005694495033699221150D+00
    x(30) = -0.824147324402412861055989047706D+00
    x(31) = -0.549167211221599184571872835161D+00
    x(32) = -0.274504541753944755855051087074D+00
    x(33) = 0.0D+00
    x(34) = 0.274504541753944755855051087074D+00
    x(35) = 0.549167211221599184571872835161D+00
    x(36) = 0.824147324402412861055989047706D+00
    x(37) = 1.09960660005694495033699221150D+00
    x(38) = 1.37571042772366833613088581439D+00
    x(39) = 1.65262921904032541019598694267D+00
    x(40) = 1.93053987597722550177796442890D+00
    x(41) = 2.20962741516918436363972079745D+00
    x(42) = 2.49008679530393550665251461372D+00
    x(43) = 2.77212500515709167382401243242D+00
    x(44) = 3.05596348432867100471992904351D+00
    x(45) = 3.34184096844683005348576537026D+00
    x(46) = 3.63001687763289533380743256995D+00
    x(47) = 3.92077540444472380734866143333D+00
    x(48) = 4.21443050997195460766227433230D+00
    x(49) = 4.51133211136821338520053781440D+00
    x(50) = 4.81187385202746476469375701958D+00
    x(51) = 5.11650300472141247403405938903D+00
    x(52) = 5.42573329769734933377879520007D+00
    x(53) = 5.74016182369022552144785845004D+00
    x(54) = 6.06049177883150521431844958114D+00
    x(55) = 6.38756373978709108766264932629D+00
    x(56) = 6.72239982016573443713269171234D+00
    x(57) = 7.06626794030689283695626648500D+00
    x(58) = 7.42077883436632423648540822493D+00
    x(59) = 7.78803908298957078257446963177D+00
    x(60) = 8.17090617805258532129873946832D+00
    x(61) = 8.57344474441790920512961590649D+00
    x(62) = 9.00182332295913301957373193078D+00
    x(63) = 9.46632932015538456209230219583D+00
    x(64) = 9.98694169167668475289516351866D+00
    x(65) = 10.6160229818782811890575602362D+00

    w(1) = 8.25161081325244640518686536873D-50
    w(2) = 2.70767584528327632245086261566D-44
    w(3) = 5.89628446597893219238447711362D-40
    w(4) = 2.8541849032786262808377028501D-36
    w(5) = 4.95258625502059879210418105309D-33
    w(6) = 3.96328698707468682361835959189D-30
    w(7) = 1.70591158107580273148997822331D-27
    w(8) = 4.37697419487184691809226004173D-25
    w(9) = 7.20161078913500757836854034749D-23
    w(10) = 8.0222187354240312838311535001D-21
    w(11) = 6.30789104558609987896303941119D-19
    w(12) = 3.61819961904286485492939434525D-17
    w(13) = 1.55466357223809604941702812296D-15
    w(14) = 5.11391748171652449009988302839D-14
    w(15) = 1.31125161063902569430172028735D-12
    w(16) = 2.66086534779295548413319751434D-11
    w(17) = 4.32865615344850974821379264835D-10
    w(18) = 5.70758293277877491250362877931D-09
    w(19) = 6.15779622145053848599380659292D-08
    w(20) = 5.48045603501799498244047819842D-07
    w(21) = 4.05224939102373376093012342174D-06
    w(22) = 0.0000250453428904958321201946621231D+00
    w(23) = 0.000130082916298451204382435919638D+00
    w(24) = 0.000570398967523771524725931177791D+00
    w(25) = 0.00211998163203684165580510045255D+00
    w(26) = 0.00670140453800573713948633700424D+00
    w(27) = 0.018069433112703589006399924887D+00
    w(28) = 0.0416611087624784398909512089873D+00
    w(29) = 0.0823001633697352251543326980867D+00
    w(30) = 0.139526139482843953007755621004D+00
    w(31) = 0.203250574154441897747728738409D+00
    w(32) = 0.254628811852790103887643365928D+00
    w(33) = 0.274478226559263167375288621205D+00
    w(34) = 0.254628811852790103887643365928D+00
    w(35) = 0.203250574154441897747728738409D+00
    w(36) = 0.139526139482843953007755621004D+00
    w(37) = 0.0823001633697352251543326980867D+00
    w(38) = 0.0416611087624784398909512089873D+00
    w(39) = 0.018069433112703589006399924887D+00
    w(40) = 0.00670140453800573713948633700424D+00
    w(41) = 0.00211998163203684165580510045255D+00
    w(42) = 0.000570398967523771524725931177791D+00
    w(43) = 0.000130082916298451204382435919638D+00
    w(44) = 0.0000250453428904958321201946621231D+00
    w(45) = 4.05224939102373376093012342174D-06
    w(46) = 5.48045603501799498244047819842D-07
    w(47) = 6.15779622145053848599380659292D-08
    w(48) = 5.70758293277877491250362877931D-09
    w(49) = 4.32865615344850974821379264835D-10
    w(50) = 2.66086534779295548413319751434D-11
    w(51) = 1.31125161063902569430172028735D-12
    w(52) = 5.11391748171652449009988302839D-14
    w(53) = 1.55466357223809604941702812296D-15
    w(54) = 3.61819961904286485492939434525D-17
    w(55) = 6.30789104558609987896303941119D-19
    w(56) = 8.0222187354240312838311535001D-21
    w(57) = 7.20161078913500757836854034749D-23
    w(58) = 4.37697419487184691809226004173D-25
    w(59) = 1.70591158107580273148997822331D-27
    w(60) = 3.96328698707468682361835959189D-30
    w(61) = 4.95258625502059879210418105309D-33
    w(62) = 2.8541849032786262808377028501D-36
    w(63) = 5.89628446597893219238447711362D-40
    w(64) = 2.70767584528327632245086261566D-44
    w(65) = 8.25161081325244640518686536873D-50

  else if ( n == 127 ) then

    x(1) = -15.2283381481673509782469544335D+00
    x(2) = -14.6695951588339726327463541129D+00
    x(3) = -14.2090859952848707551682442509D+00
    x(4) = -13.7997222902116766346452467467D+00
    x(5) = -13.4235185900709500624382583219D+00
    x(6) = -13.0712086604746019015839954396D+00
    x(7) = -12.7372356524156863381380039241D+00
    x(8) = -12.4179393788697158054458796241D+00
    x(9) = -12.1107490209477476001321235081D+00
    x(10) = -11.8137721982677271951345841362D+00
    x(11) = -11.5255651125726965991678885886D+00
    x(12) = -11.2449945837855434451943841943D+00
    x(13) = -10.9711505698402474234230402639D+00
    x(14) = -10.7032882010274813476709407447D+00
    x(15) = -10.4407879577727728677425917980D+00
    x(16) = -10.1831274734503438886241264504D+00
    x(17) = -9.92986104951142507368470042737D+00
    x(18) = -9.68060444124747280381507127327D+00
    x(19) = -9.43502333898816501350195985063D+00
    x(20) = -9.19282449884603057157741950525D+00
    x(21) = -8.95374881085654043238078901700D+00
    x(22) = -8.71756580870763073638339995485D+00
    x(23) = -8.48406926898324733260971803400D+00
    x(24) = -8.25307364544571565796941242439D+00
    x(25) = -8.02441115147033755785947397968D+00
    x(26) = -7.79792935138701054208291204556D+00
    x(27) = -7.57348915560834540228349607633D+00
    x(28) = -7.35096313922690527019612580437D+00
    x(29) = -7.13023412203507106680640257134D+00
    x(30) = -6.91119396154657131974656331094D+00
    x(31) = -6.69374252087582941900744173817D+00
    x(32) = -6.47778678116453654481449038215D+00
    x(33) = -6.26324007427373543456097238571D+00
    x(34) = -6.05002141614198456944654744824D+00
    x(35) = -5.83805492487741873866016908078D+00
    x(36) = -5.62726931054648166594234557949D+00
    x(37) = -5.41759742592432407228484258729D+00
    x(38) = -5.20897586931539835875702583722D+00
    x(39) = -5.00134463203863600385208091074D+00
    x(40) = -4.79464678437649250097485099309D+00
    x(41) = -4.58882819476983729516064850312D+00
    x(42) = -4.38383727784647362942537444075D+00
    x(43) = -4.17962476753520313494211898924D+00
    x(44) = -3.97614351206733559160358141959D+00
    x(45) = -3.77334828812505267210046784001D+00
    x(46) = -3.57119563177821804471997564852D+00
    x(47) = -3.36964368417173978966436292400D+00
    x(48) = -3.16865205019536301918577982615D+00
    x(49) = -2.96818166859559102677616495215D+00
    x(50) = -2.76819469218240588012265459589D+00
    x(51) = -2.56865437694735017231440130224D+00
    x(52) = -2.36952497904904010800124746457D+00
    x(53) = -2.17077165874115068794984980837D+00
    x(54) = -1.97236039041950200793247432276D+00
    x(55) = -1.77425787805167915846764421037D+00
    x(56) = -1.57643147532678013155195976219D+00
    x(57) = -1.37884910992617780914415570537D+00
    x(58) = -1.18147921137006858486785835984D+00
    x(59) = -0.984290641940272777265689842138D+00
    x(60) = -0.787252630218250341515968318790D+00
    x(61) = -0.590334706809421021422304393461D+00
    x(62) = -0.393506641851301365680378262002D+00
    x(63) = -0.196738383924232519642722397371D+00
    x(64) = 0.0D+00
    x(65) = 0.196738383924232519642722397371D+00
    x(66) = 0.393506641851301365680378262002D+00
    x(67) = 0.590334706809421021422304393461D+00
    x(68) = 0.787252630218250341515968318790D+00
    x(69) = 0.984290641940272777265689842138D+00
    x(70) = 1.18147921137006858486785835984D+00
    x(71) = 1.37884910992617780914415570537D+00
    x(72) = 1.57643147532678013155195976219D+00
    x(73) = 1.77425787805167915846764421037D+00
    x(74) = 1.97236039041950200793247432276D+00
    x(75) = 2.17077165874115068794984980837D+00
    x(76) = 2.36952497904904010800124746457D+00
    x(77) = 2.56865437694735017231440130224D+00
    x(78) = 2.76819469218240588012265459589D+00
    x(79) = 2.96818166859559102677616495215D+00
    x(80) = 3.16865205019536301918577982615D+00
    x(81) = 3.36964368417173978966436292400D+00
    x(82) = 3.57119563177821804471997564852D+00
    x(83) = 3.77334828812505267210046784001D+00
    x(84) = 3.97614351206733559160358141959D+00
    x(85) = 4.17962476753520313494211898924D+00
    x(86) = 4.38383727784647362942537444075D+00
    x(87) = 4.58882819476983729516064850312D+00
    x(88) = 4.79464678437649250097485099309D+00
    x(89) = 5.00134463203863600385208091074D+00
    x(90) = 5.20897586931539835875702583722D+00
    x(91) = 5.41759742592432407228484258729D+00
    x(92) = 5.62726931054648166594234557949D+00
    x(93) = 5.83805492487741873866016908078D+00
    x(94) = 6.05002141614198456944654744824D+00
    x(95) = 6.26324007427373543456097238571D+00
    x(96) = 6.47778678116453654481449038215D+00
    x(97) = 6.69374252087582941900744173817D+00
    x(98) = 6.91119396154657131974656331094D+00
    x(99) = 7.13023412203507106680640257134D+00
    x(100) = 7.35096313922690527019612580437D+00
    x(101) = 7.57348915560834540228349607633D+00
    x(102) = 7.79792935138701054208291204556D+00
    x(103) = 8.02441115147033755785947397968D+00
    x(104) = 8.25307364544571565796941242439D+00
    x(105) = 8.48406926898324733260971803400D+00
    x(106) = 8.71756580870763073638339995485D+00
    x(107) = 8.95374881085654043238078901700D+00
    x(108) = 9.19282449884603057157741950525D+00
    x(109) = 9.43502333898816501350195985063D+00
    x(110) = 9.68060444124747280381507127327D+00
    x(111) = 9.92986104951142507368470042737D+00
    x(112) = 10.1831274734503438886241264504D+00
    x(113) = 10.4407879577727728677425917980D+00
    x(114) = 10.7032882010274813476709407447D+00
    x(115) = 10.9711505698402474234230402639D+00
    x(116) = 11.2449945837855434451943841943D+00
    x(117) = 11.5255651125726965991678885886D+00
    x(118) = 11.8137721982677271951345841362D+00
    x(119) = 12.1107490209477476001321235081D+00
    x(120) = 12.4179393788697158054458796241D+00
    x(121) = 12.7372356524156863381380039241D+00
    x(122) = 13.0712086604746019015839954396D+00
    x(123) = 13.4235185900709500624382583219D+00
    x(124) = 13.7997222902116766346452467467D+00
    x(125) = 14.2090859952848707551682442509D+00
    x(126) = 14.6695951588339726327463541129D+00
    x(127) = 15.2283381481673509782469544335D+00

    w(1) = 1.25044975770895101066558695394D-101
    w(2) = 1.72727980594728851329952877284D-94
    w(3) = 8.93216815722645216635320162557D-89
    w(4) = 7.7306185241134158744827181222D-84
    w(5) = 2.01439576527109443920782513994D-79
    w(6) = 2.15037147336771602203551878273D-75
    w(7) = 1.13419242086298913875376620343D-71
    w(8) = 3.34891390118992716444169809114D-68
    w(9) = 6.04865489642049179016214753843D-65
    w(10) = 7.13750929465743002965122123123D-62
    w(11) = 5.78845633750656959788340019085D-59
    w(12) = 3.3581166223962736386929935773D-56
    w(13) = 1.4394641949298720336141068619D-53
    w(14) = 4.68218083833618292793410025836D-51
    w(15) = 1.18170544407210392716367665268D-48
    w(16) = 2.35816591560823143778744566357D-46
    w(17) = 3.78144279409152203964384313149D-44
    w(18) = 4.9411031115925407477456893331D-42
    w(19) = 5.32553037755907921458489847863D-40
    w(20) = 4.78543906802804099967221020647D-38
    w(21) = 3.61918834460649868835433546523D-36
    w(22) = 2.3232083386415854084664074623D-34
    w(23) = 1.27533314110484056196532640642D-32
    w(24) = 6.02777538509463291699314327193D-31
    w(25) = 2.4679773241854004762148469348D-29
    w(26) = 8.8019567691972403392314252914D-28
    w(27) = 2.74824892121260880467531987939D-26
    w(28) = 7.54682189033203465872349657723D-25
    w(29) = 1.83031346363374264415878982576D-23
    w(30) = 3.93559908609832906838466602268D-22
    w(31) = 7.52931616388155067444192947319D-21
    w(32) = 1.28579977867628696999762170542D-19
    w(33) = 1.96593268885070384943390296306D-18
    w(34) = 2.69865119072980851232572568063D-17
    w(35) = 3.33444143033026256341061235315D-16
    w(36) = 3.71733031252663248624409938613D-15
    w(37) = 3.74739544729563577089986076081D-14
    w(38) = 3.42300944935037851188976963928D-13
    w(39) = 2.83853037250817094975750489262D-12
    w(40) = 2.14069202905212884993201956606D-11
    w(41) = 1.47063312734774830028408333227D-10
    w(42) = 9.21739409677215086782446989876D-10
    w(43) = 5.27816639371369729333040255118D-09
    w(44) = 2.76504970450371674155194812923D-08
    w(45) = 1.32678558425807549298485884004D-07
    w(46) = 5.83809442762947462901022315301D-07
    w(47) = 2.35815617248490159838145978859D-06
    w(48) = 8.75244680345528247507614056972D-06
    w(49) = 0.0000298767905360019901790649251988D+00
    w(50) = 0.0000938744357203646866361259710004D+00
    w(51) = 0.000271707626280157286781639661883D+00
    w(52) = 0.000724939297427239633212185817821D+00
    w(53) = 0.0017841208326818955520088211458D+00
    w(54) = 0.00405248551861722466559241860023D+00
    w(55) = 0.00850002630418086349941683729112D+00
    w(56) = 0.0164711422416609467530350356258D+00
    w(57) = 0.0294992962483054353948393364098D+00
    w(58) = 0.0488473871144520262535428484316D+00
    w(59) = 0.074807989768816537216026182806D+00
    w(60) = 0.10598520508123912472195529192D+00
    w(61) = 0.138939453090947794093360848265D+00
    w(62) = 0.168562360742603870987330592834D+00
    w(63) = 0.189278495801793364889704841035D+00
    w(64) = 0.196733406888845140995323677102D+00
    w(65) = 0.189278495801793364889704841035D+00
    w(66) = 0.168562360742603870987330592834D+00
    w(67) = 0.138939453090947794093360848265D+00
    w(68) = 0.10598520508123912472195529192D+00
    w(69) = 0.074807989768816537216026182806D+00
    w(70) = 0.0488473871144520262535428484316D+00
    w(71) = 0.0294992962483054353948393364098D+00
    w(72) = 0.0164711422416609467530350356258D+00
    w(73) = 0.00850002630418086349941683729112D+00
    w(74) = 0.00405248551861722466559241860023D+00
    w(75) = 0.0017841208326818955520088211458D+00
    w(76) = 0.000724939297427239633212185817821D+00
    w(77) = 0.000271707626280157286781639661883D+00
    w(78) = 0.0000938744357203646866361259710004D+00
    w(79) = 0.0000298767905360019901790649251988D+00
    w(80) = 8.75244680345528247507614056972D-06
    w(81) = 2.35815617248490159838145978859D-06
    w(82) = 5.83809442762947462901022315301D-07
    w(83) = 1.32678558425807549298485884004D-07
    w(84) = 2.76504970450371674155194812923D-08
    w(85) = 5.27816639371369729333040255118D-09
    w(86) = 9.21739409677215086782446989876D-10
    w(87) = 1.47063312734774830028408333227D-10
    w(88) = 2.14069202905212884993201956606D-11
    w(89) = 2.83853037250817094975750489262D-12
    w(90) = 3.42300944935037851188976963928D-13
    w(91) = 3.74739544729563577089986076081D-14
    w(92) = 3.71733031252663248624409938613D-15
    w(93) = 3.33444143033026256341061235315D-16
    w(94) = 2.69865119072980851232572568063D-17
    w(95) = 1.96593268885070384943390296306D-18
    w(96) = 1.28579977867628696999762170542D-19
    w(97) = 7.52931616388155067444192947319D-21
    w(98) = 3.93559908609832906838466602268D-22
    w(99) = 1.83031346363374264415878982576D-23
    w(100) = 7.54682189033203465872349657723D-25
    w(101) = 2.74824892121260880467531987939D-26
    w(102) = 8.8019567691972403392314252914D-28
    w(103) = 2.4679773241854004762148469348D-29
    w(104) = 6.02777538509463291699314327193D-31
    w(105) = 1.27533314110484056196532640642D-32
    w(106) = 2.3232083386415854084664074623D-34
    w(107) = 3.61918834460649868835433546523D-36
    w(108) = 4.78543906802804099967221020647D-38
    w(109) = 5.32553037755907921458489847863D-40
    w(110) = 4.9411031115925407477456893331D-42
    w(111) = 3.78144279409152203964384313149D-44
    w(112) = 2.35816591560823143778744566357D-46
    w(113) = 1.18170544407210392716367665268D-48
    w(114) = 4.68218083833618292793410025836D-51
    w(115) = 1.4394641949298720336141068619D-53
    w(116) = 3.3581166223962736386929935773D-56
    w(117) = 5.78845633750656959788340019085D-59
    w(118) = 7.13750929465743002965122123123D-62
    w(119) = 6.04865489642049179016214753843D-65
    w(120) = 3.34891390118992716444169809114D-68
    w(121) = 1.13419242086298913875376620343D-71
    w(122) = 2.15037147336771602203551878273D-75
    w(123) = 2.01439576527109443920782513994D-79
    w(124) = 7.7306185241134158744827181222D-84
    w(125) = 8.93216815722645216635320162557D-89
    w(126) = 1.72727980594728851329952877284D-94
    w(127) = 1.25044975770895101066558695394D-101

  else if ( n == 128 ) then

    x(1) = -15.2918197668827409717467886552D+00
    x(2) = -14.7338424735892990556131447177D+00
    x(3) = -14.2739813047878355625094356564D+00
    x(4) = -13.8652069844762415897768433203D+00
    x(5) = -13.4895564126231418263791177750D+00
    x(6) = -13.1377747880276511010650586719D+00
    x(7) = -12.8043120820671312950137141654D+00
    x(8) = -12.4855125853494481606990566084D+00
    x(9) = -12.1788086198312463132740693095D+00
    x(10) = -11.8823101188783115808359168093D+00
    x(11) = -11.5945750547414517467820845908D+00
    x(12) = -11.3144716442899779172120028451D+00
    x(13) = -11.0410909760196333842428986719D+00
    x(14) = -10.7736891151614406713116609896D+00
    x(15) = -10.5116473299148686173941369279D+00
    x(16) = -10.2544439284709307170245436605D+00
    x(17) = -10.0016337989301228460111363000D+00
    x(18) = -9.75283321343916867454942614151D+00
    x(19) = -9.50770832327905657695490182674D+00
    x(20) = -9.26596630029617592185364037517D+00
    x(21) = -9.02734841339478834482665573280D+00
    x(22) = -8.79162454488868640635040291427D+00
    x(23) = -8.55858879506450828896030380072D+00
    x(24) = -8.32805592079014664500802003672D+00
    x(25) = -8.09985842150789607545794348110D+00
    x(26) = -7.87384413353543446678710891884D+00
    x(27) = -7.64987422768100656113184995327D+00
    x(28) = -7.42782152995230111565739552073D+00
    x(29) = -7.20756910338733385441779947109D+00
    x(30) = -6.98900904264477401185223744438D+00
    x(31) = -6.77204144325592885820688621877D+00
    x(32) = -6.55657351526448288962578894289D+00
    x(33) = -6.34251881700177947172938955573D+00
    x(34) = -6.12979658942216202462059597877D+00
    x(35) = -5.91833117508581167511681743446D+00
    x(36) = -5.70805150876808626177490879113D+00
    x(37) = -5.49889066897390948452218926009D+00
    x(38) = -5.29078548147717957643674180866D+00
    x(39) = -5.08367616748933990673505368300D+00
    x(40) = -4.87750603026481441216755173491D+00
    x(41) = -4.67222117493263892214567470373D+00
    x(42) = -4.46777025714858268344631831723D+00
    x(43) = -4.26410425682551915674979043600D+00
    x(44) = -4.06117627374927282427754765790D+00
    x(45) = -3.85894134234428182659062673118D+00
    x(46) = -3.65735626323530809623058740618D+00
    x(47) = -3.45637944957173748220943445337D+00
    x(48) = -3.25597078635065934665290567700D+00
    x(49) = -3.05609150120268005595784091684D+00
    x(50) = -2.85670404529740528265184910544D+00
    x(51) = -2.65777198318948399631081621992D+00
    x(52) = -2.45925989056573940193677619199D+00
    x(53) = -2.26113325897306228028420817752D+00
    x(54) = -2.06335840670856597768175136750D+00
    x(55) = -1.86590239514059869664912407799D+00
    x(56) = -1.66873294980372363048660121191D+00
    x(57) = -1.47181838567448600067837560546D+00
    x(58) = -1.27512753608915832143251082623D+00
    x(59) = -1.07862968481090893047100757570D+00
    x(60) = -0.882294500792981406000508343227D+00
    x(61) = -0.686091975217334872045286432691D+00
    x(62) = -0.489992360415458918089044385637D+00
    x(63) = -0.293966110300295702813351867404D+00
    x(64) = -0.0979838219558189543137713246862D+00
    x(65) = 0.0979838219558189543137713246862D+00
    x(66) = 0.293966110300295702813351867404D+00
    x(67) = 0.489992360415458918089044385637D+00
    x(68) = 0.686091975217334872045286432691D+00
    x(69) = 0.882294500792981406000508343227D+00
    x(70) = 1.07862968481090893047100757570D+00
    x(71) = 1.27512753608915832143251082623D+00
    x(72) = 1.47181838567448600067837560546D+00
    x(73) = 1.66873294980372363048660121191D+00
    x(74) = 1.86590239514059869664912407799D+00
    x(75) = 2.06335840670856597768175136750D+00
    x(76) = 2.26113325897306228028420817752D+00
    x(77) = 2.45925989056573940193677619199D+00
    x(78) = 2.65777198318948399631081621992D+00
    x(79) = 2.85670404529740528265184910544D+00
    x(80) = 3.05609150120268005595784091684D+00
    x(81) = 3.25597078635065934665290567700D+00
    x(82) = 3.45637944957173748220943445337D+00
    x(83) = 3.65735626323530809623058740618D+00
    x(84) = 3.85894134234428182659062673118D+00
    x(85) = 4.06117627374927282427754765790D+00
    x(86) = 4.26410425682551915674979043600D+00
    x(87) = 4.46777025714858268344631831723D+00
    x(88) = 4.67222117493263892214567470373D+00
    x(89) = 4.87750603026481441216755173491D+00
    x(90) = 5.08367616748933990673505368300D+00
    x(91) = 5.29078548147717957643674180866D+00
    x(92) = 5.49889066897390948452218926009D+00
    x(93) = 5.70805150876808626177490879113D+00
    x(94) = 5.91833117508581167511681743446D+00
    x(95) = 6.12979658942216202462059597877D+00
    x(96) = 6.34251881700177947172938955573D+00
    x(97) = 6.55657351526448288962578894289D+00
    x(98) = 6.77204144325592885820688621877D+00
    x(99) = 6.98900904264477401185223744438D+00
    x(100) = 7.20756910338733385441779947109D+00
    x(101) = 7.42782152995230111565739552073D+00
    x(102) = 7.64987422768100656113184995327D+00
    x(103) = 7.87384413353543446678710891884D+00
    x(104) = 8.09985842150789607545794348110D+00
    x(105) = 8.32805592079014664500802003672D+00
    x(106) = 8.55858879506450828896030380072D+00
    x(107) = 8.79162454488868640635040291427D+00
    x(108) = 9.02734841339478834482665573280D+00
    x(109) = 9.26596630029617592185364037517D+00
    x(110) = 9.50770832327905657695490182674D+00
    x(111) = 9.75283321343916867454942614151D+00
    x(112) = 10.0016337989301228460111363000D+00
    x(113) = 10.2544439284709307170245436605D+00
    x(114) = 10.5116473299148686173941369279D+00
    x(115) = 10.7736891151614406713116609896D+00
    x(116) = 11.0410909760196333842428986719D+00
    x(117) = 11.3144716442899779172120028451D+00
    x(118) = 11.5945750547414517467820845908D+00
    x(119) = 11.8823101188783115808359168093D+00
    x(120) = 12.1788086198312463132740693095D+00
    x(121) = 12.4855125853494481606990566084D+00
    x(122) = 12.8043120820671312950137141654D+00
    x(123) = 13.1377747880276511010650586719D+00
    x(124) = 13.4895564126231418263791177750D+00
    x(125) = 13.8652069844762415897768433203D+00
    x(126) = 14.2739813047878355625094356564D+00
    x(127) = 14.7338424735892990556131447177D+00
    x(128) = 15.2918197668827409717467886552D+00

    w(1) = 1.79906598010928472082336338805D-102
    w(2) = 2.60817240240911107924885148459D-95
    w(3) = 1.40468977131508863479865725345D-89
    w(4) = 1.2612494833385383033093221663D-84
    w(5) = 3.4012300869366371268669286673D-80
    w(6) = 3.75121586880472499656274624235D-76
    w(7) = 2.04158579724398501580069247379D-72
    w(8) = 6.21424416183031366240930730224D-69
    w(9) = 1.15615516409637521334725409468D-65
    w(10) = 1.40446725774048726044186592003D-62
    w(11) = 1.17197850121298051738559888373D-59
    w(12) = 6.9930729240519559879874741506D-57
    w(13) = 3.08207738333929868710425541163D-54
    w(14) = 1.03048625205569473422672330856D-51
    w(15) = 2.67274375173606785452021989916D-49
    w(16) = 5.48021702897879649820616283051D-47
    w(17) = 9.02804013878664400917961564574D-45
    w(18) = 1.21177953413059190735434940091D-42
    w(19) = 1.34149748176436936696556841563D-40
    w(20) = 1.2380855579763680376188381669D-38
    w(21) = 9.6167080679675069775952182446D-37
    w(22) = 6.33991352636648906076753997388D-35
    w(23) = 3.57437889587942107216457034803D-33
    w(24) = 1.73510302028206120881601688138D-31
    w(25) = 7.29654500676840425381868704321D-30
    w(26) = 2.67292362005807324017266437183D-28
    w(27) = 8.5728304837693537445493254974D-27
    w(28) = 2.41840345964766496960390574396D-25
    w(29) = 6.02598403200645428864656987226D-24
    w(30) = 1.33136785903358960440599429474D-22
    w(31) = 2.61745758393481115586873166674D-21
    w(32) = 4.59400767732972159221172605389D-20
    w(33) = 7.22010731692829201964437734131D-19
    w(34) = 1.01893323042329252403658204469D-17
    w(35) = 1.2945481593393715343954569556D-16
    w(36) = 1.4842238375138564829118955689D-15
    w(37) = 1.53904973035354581424981070383D-14
    w(38) = 1.44634732119041656320590928428D-13
    w(39) = 1.23421448660055669081623604437D-12
    w(40) = 9.58031650873585770862066358548D-12
    w(41) = 6.77578048777455378630839649193D-11
    w(42) = 4.37318665984840344563217253619D-10
    w(43) = 2.57939722942639480114980569527D-9
    w(44) = 1.39219071529351788119578816175D-8
    w(45) = 6.88458112215435009064406266312D-8
    w(46) = 3.12287298617890308197944991751D-7
    w(47) = 1.30074700323819923351375586698D-6
    w(48) = 4.97992453259098701134099270598D-6
    w(49) = 0.0000175404858480939050383677619791D+00
    w(50) = 0.0000568874376004024109270187885882D+00
    w(51) = 0.000170014088262809409409897155763D+00
    w(52) = 0.000468551537808411365479802126842D+00
    w(53) = 0.00119156381445716723911680561041D+00
    w(54) = 0.00279783940160578927319080368252D+00
    w(55) = 0.00606886240692588762066801419927D+00
    w(56) = 0.0121669188644693394910166328856D+00
    w(57) = 0.0225543101678244224102498222492D+00
    w(58) = 0.0386739548106369026550248867136D+00
    w(59) = 0.061360721004490065664651069257D+00
    w(60) = 0.090108678376448919548057439804D+00
    w(61) = 0.122503273164135694618664605611D+00
    w(62) = 0.154210435298354383363527713284D+00
    w(63) = 0.179773083907799264988697956102D+00
    w(64) = 0.194097611864087756977697028723D+00
    w(65) = 0.194097611864087756977697028723D+00
    w(66) = 0.179773083907799264988697956102D+00
    w(67) = 0.154210435298354383363527713284D+00
    w(68) = 0.122503273164135694618664605611D+00
    w(69) = 0.090108678376448919548057439804D+00
    w(70) = 0.061360721004490065664651069257D+00
    w(71) = 0.0386739548106369026550248867136D+00
    w(72) = 0.0225543101678244224102498222492D+00
    w(73) = 0.0121669188644693394910166328856D+00
    w(74) = 0.00606886240692588762066801419927D+00
    w(75) = 0.00279783940160578927319080368252D+00
    w(76) = 0.00119156381445716723911680561041D+00
    w(77) = 0.000468551537808411365479802126842D+00
    w(78) = 0.000170014088262809409409897155763D+00
    w(79) = 0.0000568874376004024109270187885882D+00
    w(80) = 0.0000175404858480939050383677619791D+00
    w(81) = 4.97992453259098701134099270598D-06
    w(82) = 1.30074700323819923351375586698D-06
    w(83) = 3.12287298617890308197944991751D-07
    w(84) = 6.88458112215435009064406266312D-08
    w(85) = 1.39219071529351788119578816175D-08
    w(86) = 2.57939722942639480114980569527D-09
    w(87) = 4.37318665984840344563217253619D-10
    w(88) = 6.77578048777455378630839649193D-11
    w(89) = 9.58031650873585770862066358548D-12
    w(90) = 1.23421448660055669081623604437D-12
    w(91) = 1.44634732119041656320590928428D-13
    w(92) = 1.53904973035354581424981070383D-14
    w(93) = 1.4842238375138564829118955689D-15
    w(94) = 1.2945481593393715343954569556D-16
    w(95) = 1.01893323042329252403658204469D-17
    w(96) = 7.22010731692829201964437734131D-19
    w(97) = 4.59400767732972159221172605389D-20
    w(98) = 2.61745758393481115586873166674D-21
    w(99) = 1.33136785903358960440599429474D-22
    w(100) = 6.02598403200645428864656987226D-24
    w(101) = 2.41840345964766496960390574396D-25
    w(102) = 8.5728304837693537445493254974D-27
    w(103) = 2.67292362005807324017266437183D-28
    w(104) = 7.29654500676840425381868704321D-30
    w(105) = 1.73510302028206120881601688138D-31
    w(106) = 3.57437889587942107216457034803D-33
    w(107) = 6.33991352636648906076753997388D-35
    w(108) = 9.6167080679675069775952182446D-37
    w(109) = 1.2380855579763680376188381669D-38
    w(110) = 1.34149748176436936696556841563D-40
    w(111) = 1.21177953413059190735434940091D-42
    w(112) = 9.02804013878664400917961564574D-45
    w(113) = 5.48021702897879649820616283051D-47
    w(114) = 2.67274375173606785452021989916D-49
    w(115) = 1.03048625205569473422672330856D-51
    w(116) = 3.08207738333929868710425541163D-54
    w(117) = 6.9930729240519559879874741506D-57
    w(118) = 1.17197850121298051738559888373D-59
    w(119) = 1.40446725774048726044186592003D-62
    w(120) = 1.15615516409637521334725409468D-65
    w(121) = 6.21424416183031366240930730224D-69
    w(122) = 2.04158579724398501580069247379D-72
    w(123) = 3.75121586880472499656274624235D-76
    w(124) = 3.4012300869366371268669286673D-80
    w(125) = 1.2612494833385383033093221663D-84
    w(126) = 1.40468977131508863479865725345D-89
    w(127) = 2.60817240240911107924885148459D-95
    w(128) = 1.79906598010928472082336338805D-102

  else if ( n == 129 ) then

    x(1) = -15.3550496746831285549167746019D+00
    x(2) = -14.7978308964903080628845608050D+00
    x(3) = -14.3386115290089672811362217078D+00
    x(4) = -13.9304208664791805435265533989D+00
    x(5) = -13.5553179661308567022816946453D+00
    x(6) = -13.2040593596741921982903144147D+00
    x(7) = -12.8711017789036282758938040681D+00
    x(8) = -12.5527939524445397878411009214D+00
    x(9) = -12.2465713150240016840404064819D+00
    x(10) = -11.9505460927691823148587203418D+00
    x(11) = -11.6632780120689523111895976521D+00
    x(12) = -11.3836366735119364919041401601D+00
    x(13) = -11.1107142851416459382067369906D+00
    x(14) = -10.8437678377155441232588867872D+00
    x(15) = -10.5821793789735138638177686355D+00
    x(16) = -10.3254278845735933555309803756D+00
    x(17) = -10.0730688225840385168071109595D+00
    x(18) = -9.82471897583315292981163227664D+00
    x(19) = -9.58004495076523054718996925368D+00
    x(20) = -9.33875432946329683313144773753D+00
    x(21) = -9.10058875441877630705419698871D+00
    x(22) = -8.86531845144460445006059884245D+00
    x(23) = -8.63273783950843552405601767759D+00
    x(24) = -8.40266197362535012572889592790D+00
    x(25) = -8.17492363437398978801516910338D+00
    x(26) = -7.94937092512605027069441566240D+00
    x(27) = -7.72586527212128476858225880726D+00
    x(28) = -7.50427974726352240601473698475D+00
    x(29) = -7.28449765174015372725258835236D+00
    x(30) = -7.06641131216039069912389691607D+00
    x(31) = -6.84992105116014925339717178568D+00
    x(32) = -6.63493430223598850302862096202D+00
    x(33) = -6.42136484458510366813184121466D+00
    x(34) = -6.20913213839958724139275362341D+00
    x(35) = -5.99816074472179863235247556956D+00
    x(36) = -5.78837981685588271189500366573D+00
    x(37) = -5.57972265262736517721195076411D+00
    x(38) = -5.37212629862206810544406569124D+00
    x(39) = -5.16553119901808798749445925424D+00
    x(40) = -4.95988088282680494441019859192D+00
    x(41) = -4.75512168433945660698412105871D+00
    x(42) = -4.55120249237992751786552441506D+00
    x(43) = -4.34807452462720039287086617988D+00
    x(44) = -4.14569112381985885943227755134D+00
    x(45) = -3.94400757311174034591077592589D+00
    x(46) = -3.74298092822942697020355662873D+00
    x(47) = -3.54256986440235708810942555303D+00
    x(48) = -3.34273453630584653564135188865D+00
    x(49) = -3.14343644948499981855157130224D+00
    x(50) = -2.94463834192045889859029200431D+00
    x(51) = -2.74630407456093830886162508356D+00
    x(52) = -2.54839852978722422234239962660D+00
    x(53) = -2.35088751689161374336515868011D+00
    x(54) = -2.15373768375879465132813086025D+00
    x(55) = -1.95691643402151503173496953711D+00
    x(56) = -1.76039184903921457004681976795D+00
    x(57) = -1.56413261411186298656139617983D+00
    x(58) = -1.36810794839605047973791359079D+00
    x(59) = -1.17228753803712318264391216760D+00
    x(60) = -0.976641472070867557126534700881D+00
    x(61) = -0.781140180681760289238140546741D+00
    x(62) = -0.585754375432805697119652981369D+00
    x(63) = -0.390454991105046004780513867383D+00
    x(64) = -0.195213128803407573801607754230D+00
    x(65) = 0.0D+00
    x(66) = 0.195213128803407573801607754230D+00
    x(67) = 0.390454991105046004780513867383D+00
    x(68) = 0.585754375432805697119652981369D+00
    x(69) = 0.781140180681760289238140546741D+00
    x(70) = 0.976641472070867557126534700881D+00
    x(71) = 1.17228753803712318264391216760D+00
    x(72) = 1.36810794839605047973791359079D+00
    x(73) = 1.56413261411186298656139617983D+00
    x(74) = 1.76039184903921457004681976795D+00
    x(75) = 1.95691643402151503173496953711D+00
    x(76) = 2.15373768375879465132813086025D+00
    x(77) = 2.35088751689161374336515868011D+00
    x(78) = 2.54839852978722422234239962660D+00
    x(79) = 2.74630407456093830886162508356D+00
    x(80) = 2.94463834192045889859029200431D+00
    x(81) = 3.14343644948499981855157130224D+00
    x(82) = 3.34273453630584653564135188865D+00
    x(83) = 3.54256986440235708810942555303D+00
    x(84) = 3.74298092822942697020355662873D+00
    x(85) = 3.94400757311174034591077592589D+00
    x(86) = 4.14569112381985885943227755134D+00
    x(87) = 4.34807452462720039287086617988D+00
    x(88) = 4.55120249237992751786552441506D+00
    x(89) = 4.75512168433945660698412105871D+00
    x(90) = 4.95988088282680494441019859192D+00
    x(91) = 5.16553119901808798749445925424D+00
    x(92) = 5.37212629862206810544406569124D+00
    x(93) = 5.57972265262736517721195076411D+00
    x(94) = 5.78837981685588271189500366573D+00
    x(95) = 5.99816074472179863235247556956D+00
    x(96) = 6.20913213839958724139275362341D+00
    x(97) = 6.42136484458510366813184121466D+00
    x(98) = 6.63493430223598850302862096202D+00
    x(99) = 6.84992105116014925339717178568D+00
    x(100) = 7.06641131216039069912389691607D+00
    x(101) = 7.28449765174015372725258835236D+00
    x(102) = 7.50427974726352240601473698475D+00
    x(103) = 7.72586527212128476858225880726D+00
    x(104) = 7.94937092512605027069441566240D+00
    x(105) = 8.17492363437398978801516910338D+00
    x(106) = 8.40266197362535012572889592790D+00
    x(107) = 8.63273783950843552405601767759D+00
    x(108) = 8.86531845144460445006059884245D+00
    x(109) = 9.10058875441877630705419698871D+00
    x(110) = 9.33875432946329683313144773753D+00
    x(111) = 9.58004495076523054718996925368D+00
    x(112) = 9.82471897583315292981163227664D+00
    x(113) = 10.0730688225840385168071109595D+00
    x(114) = 10.3254278845735933555309803756D+00
    x(115) = 10.5821793789735138638177686355D+00
    x(116) = 10.8437678377155441232588867872D+00
    x(117) = 11.1107142851416459382067369906D+00
    x(118) = 11.3836366735119364919041401601D+00
    x(119) = 11.6632780120689523111895976521D+00
    x(120) = 11.9505460927691823148587203418D+00
    x(121) = 12.2465713150240016840404064819D+00
    x(122) = 12.5527939524445397878411009214D+00
    x(123) = 12.8711017789036282758938040681D+00
    x(124) = 13.2040593596741921982903144147D+00
    x(125) = 13.5553179661308567022816946453D+00
    x(126) = 13.9304208664791805435265533989D+00
    x(127) = 14.3386115290089672811362217078D+00
    x(128) = 14.7978308964903080628845608050D+00
    x(129) = 15.3550496746831285549167746019D+00

    w(1) = 2.58755395082114927399064139631D-103
    w(2) = 3.93601845908067608811461078697D-96
    w(3) = 2.20725529577484588586177997021D-90
    w(4) = 2.05563087297774646200941835216D-85
    w(5) = 5.73584763407311509769038083955D-81
    w(6) = 6.53456499014096713882711627986D-77
    w(7) = 3.66903606454555600244832281797D-73
    w(8) = 1.15105101975113879079427442365D-69
    w(9) = 2.20553774145133363585421051568D-66
    w(10) = 2.75763663311195533797446164671D-63
    w(11) = 2.36731747071610805241477009401D-60
    w(12) = 1.45257860403230704544333281907D-57
    w(13) = 6.58119121529392093666305170751D-55
    w(14) = 2.26137732951303228667152914802D-52
    w(15) = 6.02643011329776195986432204924D-50
    w(16) = 1.26938407638088455004457398255D-47
    w(17) = 2.14791778799787733305388107076D-45
    w(18) = 2.96092183537878053158423564486D-43
    w(19) = 3.36616090532826422441501486485D-41
    w(20) = 3.19014783528482307711547124192D-39
    w(21) = 2.54439796712780366695746038013D-37
    w(22) = 1.72239465322100711581154624691D-35
    w(23) = 9.97105538735197785176257533162D-34
    w(24) = 4.97009943352064894027841342072D-32
    w(25) = 2.14620630607238052957787041268D-30
    w(26) = 8.07377921555793000987040030256D-29
    w(27) = 2.65936924028161851577004868287D-27
    w(28) = 7.70515053183270746145645250031D-26
    w(29) = 1.97204966381589933881729892459D-24
    w(30) = 4.47579713475437089012921294273D-23
    w(31) = 9.0403370335874459959906960673D-22
    w(32) = 1.63036407035294103578729410788D-20
    w(33) = 2.63320491826449443345354482912D-19
    w(34) = 3.81944198027838553902522199764D-18
    w(35) = 4.98833273307808866457667338365D-17
    w(36) = 5.88022772755071728094452091844D-16
    w(37) = 6.27023947714728862011748531319D-15
    w(38) = 6.06072571359080078068964155295D-14
    w(39) = 5.32049070753884105044682362639D-13
    w(40) = 4.24955065877498808023415505556D-12
    w(41) = 3.09330203932473692244204789801D-11
    w(42) = 2.05524352987860630455845773203D-10
    w(43) = 1.2482251818681431772545606389D-09
    w(44) = 6.93896714453731562418029048785D-09
    w(45) = 3.53518234605234028369262582274D-08
    w(46) = 1.65252704577539544523562160076D-07
    w(47) = 7.09535030601389014268064639021D-07
    w(48) = 2.80106033677073567813925250808D-06
    w(49) = 0.0000101764715414468349837840217278D+00
    w(50) = 0.0000340541841724020078441933069804D+00
    w(51) = 0.000105047486997647220847004754607D+00
    w(52) = 0.000298922505941519029186629138321D+00
    w(53) = 0.000785197220610268688197653195861D+00
    w(54) = 0.00190506673927936544347937172051D+00
    w(55) = 0.00427162074179231114560048384305D+00
    w(56) = 0.00885609926394363549300290104701D+00
    w(57) = 0.0169845117091580731620255503875D+00
    w(58) = 0.0301436534848915408822025025241D+00
    w(59) = 0.0495245901368945546436264039895D+00
    w(60) = 0.0753454506416603410859342903275D+00
    w(61) = 0.106172669789632918045101372957D+00
    w(62) = 0.138604146980788427972651263139D+00
    w(63) = 0.167654732143619067997522798882D+00
    w(64) = 0.187923095463858179335367694936D+00
    w(65) = 0.195208341719164170910088609838D+00
    w(66) = 0.187923095463858179335367694936D+00
    w(67) = 0.167654732143619067997522798882D+00
    w(68) = 0.138604146980788427972651263139D+00
    w(69) = 0.106172669789632918045101372957D+00
    w(70) = 0.0753454506416603410859342903275D+00
    w(71) = 0.0495245901368945546436264039895D+00
    w(72) = 0.0301436534848915408822025025241D+00
    w(73) = 0.0169845117091580731620255503875D+00
    w(74) = 0.00885609926394363549300290104701D+00
    w(75) = 0.00427162074179231114560048384305D+00
    w(76) = 0.00190506673927936544347937172051D+00
    w(77) = 0.000785197220610268688197653195861D+00
    w(78) = 0.000298922505941519029186629138321D+00
    w(79) = 0.000105047486997647220847004754607D+00
    w(80) = 0.0000340541841724020078441933069804D+00
    w(81) = 0.0000101764715414468349837840217278D+00
    w(82) = 2.80106033677073567813925250808D-06
    w(83) = 7.09535030601389014268064639021D-07
    w(84) = 1.65252704577539544523562160076D-07
    w(85) = 3.53518234605234028369262582274D-08
    w(86) = 6.93896714453731562418029048785D-09
    w(87) = 1.2482251818681431772545606389D-09
    w(88) = 2.05524352987860630455845773203D-10
    w(89) = 3.09330203932473692244204789801D-11
    w(90) = 4.24955065877498808023415505556D-12
    w(91) = 5.32049070753884105044682362639D-13
    w(92) = 6.06072571359080078068964155295D-14
    w(93) = 6.27023947714728862011748531319D-15
    w(94) = 5.88022772755071728094452091844D-16
    w(95) = 4.98833273307808866457667338365D-17
    w(96) = 3.81944198027838553902522199764D-18
    w(97) = 2.63320491826449443345354482912D-19
    w(98) = 1.63036407035294103578729410788D-20
    w(99) = 9.0403370335874459959906960673D-22
    w(100) = 4.47579713475437089012921294273D-23
    w(101) = 1.97204966381589933881729892459D-24
    w(102) = 7.70515053183270746145645250031D-26
    w(103) = 2.65936924028161851577004868287D-27
    w(104) = 8.07377921555793000987040030256D-29
    w(105) = 2.14620630607238052957787041268D-30
    w(106) = 4.97009943352064894027841342072D-32
    w(107) = 9.97105538735197785176257533162D-34
    w(108) = 1.72239465322100711581154624691D-35
    w(109) = 2.54439796712780366695746038013D-37
    w(110) = 3.19014783528482307711547124192D-39
    w(111) = 3.36616090532826422441501486485D-41
    w(112) = 2.96092183537878053158423564486D-43
    w(113) = 2.14791778799787733305388107076D-45
    w(114) = 1.26938407638088455004457398255D-47
    w(115) = 6.02643011329776195986432204924D-50
    w(116) = 2.26137732951303228667152914802D-52
    w(117) = 6.58119121529392093666305170751D-55
    w(118) = 1.45257860403230704544333281907D-57
    w(119) = 2.36731747071610805241477009401D-60
    w(120) = 2.75763663311195533797446164671D-63
    w(121) = 2.20553774145133363585421051568D-66
    w(122) = 1.15105101975113879079427442365D-69
    w(123) = 3.66903606454555600244832281797D-73
    w(124) = 6.53456499014096713882711627986D-77
    w(125) = 5.73584763407311509769038083955D-81
    w(126) = 2.05563087297774646200941835216D-85
    w(127) = 2.20725529577484588586177997021D-90
    w(128) = 3.93601845908067608811461078697D-96
    w(129) = 2.58755395082114927399064139631D-103

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 to 20,'
    write ( *, '(a)' ) '  31/32/33, 63/64/65, 127/128/129.'
    stop

  end if

  return
end
subroutine hermite_ss_compute ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_SS_COMPUTE computes a Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The code uses an algorithm by Stroud and Secrest.
!
!    The abscissas are the zeros of the N-th order Hermite polynomial.
!
!    The integral:
!
!      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) s
  real ( kind = 8 ) temp
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  cc = sqrt ( pi ) * r8_gamma ( real ( n, kind = 8 ) ) / ( 2.0D+00**( n - 1 ) )

  s = ( 2.0D+00 * real ( n, kind = 8 ) + 1.0D+00 )**( 1.0D+00 / 6.0D+00 )

  do i = 1, ( n + 1 ) / 2

    if ( i == 1 ) then

      xval = s * s * s - 1.85575D+00 / s

    else if ( i == 2 ) then

      xval = xval - 1.14D+00 * ( ( real ( n, kind = 8 ) )**0.426D+00 ) / xval

    else if ( i == 3 ) then

      xval = 1.86D+00 * xval - 0.86D+00 * x(1)

    else if ( i == 4 ) then

      xval = 1.91D+00 * xval - 0.91D+00 * x(2)

    else

      xval = 2.0D+00 * xval - x(i-2)

    end if

    call hermite_ss_root ( xval, n, dp2, p1 )

    x(i) = xval
    w(i) = ( cc / dp2 ) / p1

    x(n-i+1) = - xval
    w(n-i+1) = w(i)

  end do
!
!  Reverse the abscissas.
!  Because of symmetry, the weights are unchanged,
!  and the abscissas simply change sign.
!
  x(1:n) = - x(1:n)

  return
end
subroutine hermite_ss_recur ( p2, dp2, p1, x, n )

!*****************************************************************************80
!
!! HERMITE_SS_RECUR finds the value and derivative of a Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of H(N)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of H'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of H(N-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) n
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, n

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2  = x * p1 - 0.5D+00 * ( real ( i, kind = 8 ) - 1.0D+00 ) * p0
    dp2 = x * dp1 + p1 - 0.5D+00 * ( real ( i, kind = 8 ) - 1.0D+00 ) * dp0

  end do

  return
end
subroutine hermite_ss_root ( x, n, dp2, p1 )

!*****************************************************************************80
!
!! HERMITE_SS_ROOT improves an approximate root of a Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) DP2, the value of H'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of H(N-1)(X).
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  integer ( kind = 4 ) n
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = epsilon ( eps )

  do step = 1, step_max

    call hermite_ss_recur ( p2, dp2, p1, x, n )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine imtql2 ( n, d, e, z, ierr )

!*****************************************************************************80
!
!! IMTQL2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This subroutine finds the eigenvalues and eigenvectors
!    of a symmetric tridiagonal matrix by the implicit QL method.
!    The eigenvectors of a full symmetric matrix can also
!    be found if TRED2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
!    the input matrix.  On output, the eigenvalues in ascending order.  If an
!    error exit is made, the eigenvalues are correct but
!    unordered for indices 1,2,...,IERR-1.
!
!    Input/output, real ( kind = 8 ) E(N).  On input, the subdiagonal elements
!    of the input matrix in E(2:N).  E(1) is arbitrary.  On output, E is
!    overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N,N).  On input, the transformation
!    matrix produced in the reduction by TRED2, if performed.  If the
!    eigenvectors of the tridiagonal matrix are desired, Z must contain the
!    identity matrix.  On output, Z contains orthonormal eigenvectors of the
!    symmetric tridiagonal (or full) matrix.  If an error exit is made, Z
!    contains the eigenvectors associated with the stored eigenvalues.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, for normal return,
!    J, if the J-th eigenvalue has not been determined after 30 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) tst1
  real ( kind = 8 ) tst2
  real ( kind = 8 ) z(n,n)

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do
  e(n) = 0.0D+00

  do l = 1, n

    j = 0
!
!  Look for a small sub-diagonal element.
!
105 continue

    do m = l, n

      if ( m == n ) then
        exit
      end if

      tst1 = abs ( d(m) ) + abs ( d(m+1) )
      tst2 = tst1 + abs ( e(m) )

      if ( tst2 == tst1 ) then
        exit
      end if

    end do

    p = d(l)

    if ( m == l ) then
      cycle
    end if

    if ( 30 <= j ) then
      ierr = l
      return
    end if

    j = j + 1
!
!  Form shift.
!
    g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
    r = pythag ( g, 1.0D+00 )
    g = d(m) - p + e(l) / ( g + sign ( r, g ) )
    s = 1.0D+00
    c = 1.0D+00
    p = 0.0D+00
    mml = m - l

    do ii = 1, mml

      i = m - ii
      f = s * e(i)
      b = c * e(i)
      r = pythag ( f, g )
      e(i+1) = r
!
!  Recover from underflow.
!
      if ( r == 0.0D+00 ) then
        d(i+1) = d(i+1) - p
        e(m) = 0.0D+00
        go to 105
      end if

      s = f / r
      c = g / r
      g = d(i+1) - p
      r = ( d(i) - g ) * s + 2.0D+00 * c * b
      p = s * r
      d(i+1) = g + p
      g = c * r - b
!
!  Form vector.
!
      do k = 1, n
        f = z(k,i+1)
        z(k,i+1) = s * z(k,i) + c * f
        z(k,i) = c * z(k,i) - s * f
      end do

    end do

    d(l) = d(l) - p
    e(l) = g
    e(m) = 0.0D+00
    go to 105

  end do
!
!  Order eigenvalues and eigenvectors.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then

      d(k) = d(i)
      d(i) = p

      t(1:n)   = z(1:n,i)
      z(1:n,i) = z(1:n,k)
      z(1:n,k) = t(1:n)

    end if

  end do

  return
end
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
subroutine jacobi_ek_compute ( n, alpha, beta, x, w )

!*****************************************************************************80
!
!! JACOBI_EK_COMPUTE: Elhay-Kautsky method for Gauss-Jacobi quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) (1-x)^alpha * (1+x)^beta * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
!    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) abi
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) i_r8
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = 2.0D+00**( alpha + beta + 1.0D+00 ) &
    * r8_gamma ( alpha + 1.0D+00 ) &
    * r8_gamma ( beta + 1.0D+00 ) &
    / r8_gamma ( 2.0D+00 + alpha + beta )
!
!  Define the Jacobi matrix.
!
  x(1) = ( beta - alpha ) / ( 2.0D+00 + alpha + beta )

  bj(1) = 4.0D+00 * ( 1.0D+00 + alpha ) * ( 1.0D+00 + beta ) &
    / ( ( 3.0D+00 + alpha + beta ) * ( 2.0D+00 + alpha + beta )**2 )

  do i = 2, n
    i_r8 = real ( i, kind = 8 )
    abi = 2.0D+00 * i_r8 + alpha + beta
    x(i) = ( beta + alpha ) * ( beta - alpha ) / ( ( abi - 2.0D+00 ) * abi )
    bj(i) = 4.0D+00 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) &
      * ( i_r8 + alpha + beta ) &
      / ( ( abi - 1.0D+00 ) * ( abi + 1.0D+00 ) * abi * abi )
  end do

  bj(1:n) = sqrt ( bj(1:n) )

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine jacobi_integral ( expon, alpha, beta, value )

!*****************************************************************************80
!
!! JACOBI_INTEGRAL evaluates the integral of a monomial with Jacobi weight.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= +1 ) x^expon (1-x)^alpha (1+x)^beta dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X) in the weight factor.
!
!    Input, real ( kind = 8 ) BETA, the exponent of (1+X) in the weight factor.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) s
  real ( kind = 8 ) value
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2

  c = real ( expon, kind = 8 )

  if ( mod ( expon, 2 ) == 0 ) then
    s = +1.0D+00
  else
    s = -1.0D+00
  end if

  arg1 = - alpha
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + beta + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

  arg1 = - beta
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + alpha + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value2 )

  value = r8_gamma ( 1.0D+00 + c ) * ( &
      s * r8_gamma ( 1.0D+00 + beta  ) * value1 &
    / r8_gamma ( 2.0D+00 + beta  + c ) &
    +     r8_gamma ( 1.0D+00 + alpha ) * value2 &
    / r8_gamma ( 2.0D+00 + alpha + c ) )

  return
end
subroutine jacobi_ss_compute ( n, alpha, beta, x, w )

!*****************************************************************************80
!
!! JACOBI_SS_COMPUTE computes a Gauss-Jacobi quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) (1-x)^alpha * (1+x)^beta * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    Thanks to Xu Xiang of Fudan University for pointing out that
!    an earlier implementation of this routine was incorrect!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2007
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
!    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) an
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  real ( kind = 8 ) bn
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) cc
  real ( kind = 8 ) delta
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval
!
!  Check ALPHA and BETA.
!
  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_SS_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  -1.0 < ALPHA is required.'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_SS_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  -1.0 < BETA is required.'
    stop
  end if
!
!  Set the recursion coefficients.
!
  do i = 1, n

    if ( alpha + beta == 0.0D+00 .or. beta - alpha == 0.0D+00 ) then

      b(i) = 0.0D+00

    else

      b(i) = ( alpha + beta ) * ( beta - alpha ) / &
            ( ( alpha + beta + real ( 2 * i, kind = 8 ) ) &
            * ( alpha + beta + real ( 2 * i - 2, kind = 8 ) ) )

    end if

    if ( i == 1 ) then

      c(i) = 0.0D+00

    else

      c(i) = 4.0D+00 * real ( i - 1, kind = 8 ) &
            * ( alpha + real ( i - 1, kind = 8 ) ) &
            * ( beta + real ( i - 1, kind = 8 ) ) &
            * ( alpha + beta + real ( i - 1, kind = 8 ) ) / &
            ( ( alpha + beta + real ( 2 * i - 1, kind = 8 ) ) &
            * ( alpha + beta + real ( 2 * i - 2, kind = 8 ) )**2 &
            * ( alpha + beta + real ( 2 * i - 3, kind = 8 ) ) )

    end if

  end do

  delta = r8_gamma ( alpha        + 1.0D+00 ) &
        * r8_gamma (         beta + 1.0D+00 ) &
        / r8_gamma ( alpha + beta + 2.0D+00 )

  cc = delta * 2.0D+00**( alpha + beta + 1.0D+00 ) * product ( c(2:n) )

  do i = 1, n

    if ( i == 1 ) then

      an = alpha / real ( n, kind = 8 )
      bn = beta / real ( n, kind = 8 )

      r1 = ( 1.0D+00 + alpha ) &
        * ( 2.78D+00 / ( 4.0D+00 + real ( n * n, kind = 8 ) ) &
        + 0.768D+00 * an / real ( n, kind = 8 ) )

      r2 = 1.0D+00 + 1.48D+00 * an + 0.96D+00 * bn &
        + 0.452D+00 * an * an + 0.83D+00 * an * bn

      xval = ( r2 - r1 ) / r2

    else if ( i == 2 ) then

      r1 = ( 4.1D+00 + alpha ) / &
        ( ( 1.0D+00 + alpha ) * ( 1.0D+00 + 0.156D+00 * alpha ) )

      r2 = 1.0D+00 + 0.06D+00 * ( real ( n, kind = 8 ) - 8.0D+00 ) * &
        ( 1.0D+00 + 0.12D+00 * alpha ) / real ( n, kind = 8 )

      r3 = 1.0D+00 + 0.012D+00 * beta * &
        ( 1.0D+00 + 0.25D+00 * abs ( alpha ) ) / real ( n, kind = 8 )

      xval = xval - r1 * r2 * r3 * ( 1.0D+00 - xval )

    else if ( i == 3 ) then

      r1 = ( 1.67D+00 + 0.28D+00 * alpha ) / ( 1.0D+00 + 0.37D+00 * alpha )

      r2 = 1.0D+00 + 0.22D+00 * ( real ( n, kind = 8 ) - 8.0D+00 ) &
        / real ( n, kind = 8 )

      r3 = 1.0D+00 + 8.0D+00 * beta / &
        ( ( 6.28D+00 + beta ) * real ( n * n, kind = 8 ) )

      xval = xval - r1 * r2 * r3 * ( x(1) - xval )

    else if ( i < n - 1 ) then

      xval = 3.0D+00 * x(i-1) - 3.0D+00 * x(i-2) + x(i-3)

    else if ( i == n - 1 ) then

      r1 = ( 1.0D+00 + 0.235D+00 * beta ) / ( 0.766D+00 + 0.119D+00 * beta )

      r2 = 1.0D+00 / ( 1.0D+00 + 0.639D+00 &
        * ( real ( n, kind = 8 ) - 4.0D+00 ) &
        / ( 1.0D+00 + 0.71D+00 * ( real ( n, kind = 8 ) - 4.0D+00 ) ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 20.0D+00 * alpha / ( ( 7.5D+00 + alpha ) * &
        real ( n * n, kind = 8 ) ) )

      xval = xval + r1 * r2 * r3 * ( xval - x(i-2) )

    else if ( i == n ) then

      r1 = ( 1.0D+00 + 0.37D+00 * beta ) / ( 1.67D+00 + 0.28D+00 * beta )

      r2 = 1.0D+00 / &
        ( 1.0D+00 + 0.22D+00 * ( real ( n, kind = 8 ) - 8.0D+00 ) &
        / real ( n, kind = 8 ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( n * n, kind = 8 ) ) )

      xval = xval + r1 * r2 * r3 * ( xval - x(i-2) )

    end if

    call jacobi_ss_root ( xval, n, alpha, beta, dp2, p1, b, c )

    x(i) = xval
    w(i) = cc / ( dp2 * p1 )

  end do
!
!  Reverse the data.
!
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

  return
end
subroutine jacobi_ss_recur ( p2, dp2, p1, x, n, alpha, beta, b, c )

!*****************************************************************************80
!
!! JACOBI_SS_RECUR finds the value and derivative of a Jacobi polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of J(N)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of J'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(N-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.
!
!    Input, real ( kind = 8 ) B(N), C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0D+00 )
  dp2 = 1.0D+00

  do i = 2, n

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine jacobi_ss_root ( x, n, alpha, beta, dp2, p1, b, c )

!*****************************************************************************80
!
!! JACOBI_SS_ROOT improves an approximate root of a Jacobi polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2000
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.
!
!    Output, real ( kind = 8 ) DP2, the value of J'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of J(N-1)(X).
!
!    Input, real ( kind = 8 ) B(N), C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = epsilon ( eps )

  do step = 1, step_max

    call jacobi_ss_recur ( p2, dp2, p1, x, n, alpha, beta, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine kronrod_set ( n, x, w )

!*****************************************************************************80
!
!! KRONROD_SET sets abscissas and weights for Gauss-Kronrod quadrature.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    A Kronrod rule is used in conjunction with a lower order
!    Gauss rule, and provides an efficient error estimation.
!
!    The error may be estimated as the difference in the two integral
!    approximations.
!
!    The efficiency comes about because the Kronrod rule uses the abscissas
!    of the Gauss rule, thus saving on the number of function evaluations
!    necessary.  If the Kronrod rule were replaced by a Gauss rule of
!    the same order, a higher precision integral estimate would be
!    made, but the function would have to be evaluated at many more
!    points.
!
!    The Gauss Kronrod pair of rules involves an ( N + 1 ) / 2
!    point Gauss-Legendre rule and an N point Kronrod rule.
!    Thus, the 15 point Kronrod rule should be paired with the
!    Gauss-Legendre 7 point rule.
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
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise deDoncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983,
!    ISBN: 3540125531.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N may be 15, 21, 31 or 41, corresponding to Gauss-Legendre rules of
!    order 7, 10, 15 or 20.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 15 ) then

    x(1) =  - 0.9914553711208126D+00
    x(2) =  - 0.9491079123427585D+00
    x(3) =  - 0.8648644233597691D+00
    x(4) =  - 0.7415311855993944D+00
    x(5) =  - 0.5860872354676911D+00
    x(6) =  - 0.4058451513773972D+00
    x(7) =  - 0.2077849550789850D+00
    x(8) =    0.0D+00
    x(9) =    0.2077849550789850D+00
    x(10) =   0.4058451513773972D+00
    x(11) =   0.5860872354676911D+00
    x(12) =   0.7415311855993944D+00
    x(13) =   0.8648644233597691D+00
    x(14) =   0.9491079123427585D+00
    x(15) =   0.9914553711208126D+00

    w(1) =  0.2293532201052922D-01
    w(2) =  0.6309209262997855D-01
    w(3) =  0.1047900103222502D+00
    w(4) =  0.1406532597155259D+00
    w(5) =  0.1690047266392679D+00
    w(6) =  0.1903505780647854D+00
    w(7) =  0.2044329400752989D+00
    w(8) =  0.2094821410847278D+00
    w(9) =  0.2044329400752989D+00
    w(10) = 0.1903505780647854D+00
    w(11) = 0.1690047266392679D+00
    w(12) = 0.1406532597155259D+00
    w(13) = 0.1047900103222502D+00
    w(14) = 0.6309209262997855D-01
    w(15) = 0.2293532201052922D-01

  else if ( n == 21 ) then

    x(1) =  - 0.9956571630258081D+00
    x(2) =  - 0.9739065285171717D+00
    x(3) =  - 0.9301574913557082D+00
    x(4) =  - 0.8650633666889845D+00
    x(5) =  - 0.7808177265864169D+00
    x(6) =  - 0.6794095682990244D+00
    x(7) =  - 0.5627571346686047D+00
    x(8) =  - 0.4333953941292472D+00
    x(9) =  - 0.2943928627014602D+00
    x(10) = - 0.1488743389816312D+00
    x(11) =   0.0D+00
    x(12) =   0.1488743389816312D+00
    x(13) =   0.2943928627014602D+00
    x(14) =   0.4333953941292472D+00
    x(15) =   0.5627571346686047D+00
    x(16) =   0.6794095682990244D+00
    x(17) =   0.7808177265864169D+00
    x(18) =   0.8650633666889845D+00
    x(19) =   0.9301574913557082D+00
    x(20) =   0.9739065285171717D+00
    x(21) =   0.9956571630258081D+00

    w(1) =  0.1169463886737187D-01
    w(2) =  0.3255816230796473D-01
    w(3) =  0.5475589657435200D-01
    w(4) =  0.7503967481091995D-01
    w(5) =  0.9312545458369761D-01
    w(6) =  0.1093871588022976D+00
    w(7) =  0.1234919762620659D+00
    w(8) =  0.1347092173114733D+00
    w(9) =  0.1427759385770601D+00
    w(10) = 0.1477391049013385D+00
    w(11) = 0.1494455540029169D+00
    w(12) = 0.1477391049013385D+00
    w(13) = 0.1427759385770601D+00
    w(14) = 0.1347092173114733D+00
    w(15) = 0.1234919762620659D+00
    w(16) = 0.1093871588022976D+00
    w(17) = 0.9312545458369761D-01
    w(18) = 0.7503967481091995D-01
    w(19) = 0.5475589657435200D-01
    w(20) = 0.3255816230796473D-01
    w(21) = 0.1169463886737187D-01

  else if ( n == 31 ) then

    x(1) =  - 0.9980022986933971D+00
    x(2) =  - 0.9879925180204854D+00
    x(3) =  - 0.9677390756791391D+00
    x(4) =  - 0.9372733924007059D+00
    x(5) =  - 0.8972645323440819D+00
    x(6) =  - 0.8482065834104272D+00
    x(7) =  - 0.7904185014424659D+00
    x(8) =  - 0.7244177313601700D+00
    x(9) =  - 0.6509967412974170D+00
    x(10) = - 0.5709721726085388D+00
    x(11) = - 0.4850818636402397D+00
    x(12) = - 0.3941513470775634D+00
    x(13) = - 0.2991800071531688D+00
    x(14) = - 0.2011940939974345D+00
    x(15) = - 0.1011420669187175D+00
    x(16) =   0.0D+00
    x(17) =   0.1011420669187175D+00
    x(18) =   0.2011940939974345D+00
    x(19) =   0.2991800071531688D+00
    x(20) =   0.3941513470775634D+00
    x(21) =   0.4850818636402397D+00
    x(22) =   0.5709721726085388D+00
    x(23) =   0.6509967412974170D+00
    x(24) =   0.7244177313601700D+00
    x(25) =   0.7904185014424659D+00
    x(26) =   0.8482065834104272D+00
    x(27) =   0.8972645323440819D+00
    x(28) =   0.9372733924007059D+00
    x(29) =   0.9677390756791391D+00
    x(30) =   0.9879925180204854D+00
    x(31) =   0.9980022986933971D+00

    w(1) =  0.5377479872923349D-02
    w(2) =  0.1500794732931612D-01
    w(3) =  0.2546084732671532D-01
    w(4) =  0.3534636079137585D-01
    w(5) =  0.4458975132476488D-01
    w(6) =  0.5348152469092809D-01
    w(7) =  0.6200956780067064D-01
    w(8) =  0.6985412131872826D-01
    w(9) =  0.7684968075772038D-01
    w(10) = 0.8308050282313302D-01
    w(11) = 0.8856444305621177D-01
    w(12) = 0.9312659817082532D-01
    w(13) = 0.9664272698362368D-01
    w(14) = 0.9917359872179196D-01
    w(15) = 0.1007698455238756D+00
    w(16) = 0.1013300070147915D+00
    w(17) = 0.1007698455238756D+00
    w(18) = 0.9917359872179196D-01
    w(19) = 0.9664272698362368D-01
    w(20) = 0.9312659817082532D-01
    w(21) = 0.8856444305621177D-01
    w(22) = 0.8308050282313302D-01
    w(23) = 0.7684968075772038D-01
    w(24) = 0.6985412131872826D-01
    w(25) = 0.6200956780067064D-01
    w(26) = 0.5348152469092809D-01
    w(27) = 0.4458975132476488D-01
    w(28) = 0.3534636079137585D-01
    w(29) = 0.2546084732671532D-01
    w(30) = 0.1500794732931612D-01
    w(31) = 0.5377479872923349D-02

  else if ( n == 41 ) then

    x(1) =  - 0.9988590315882777D+00
    x(2) =  - 0.9931285991850949D+00
    x(3) =  - 0.9815078774502503D+00
    x(4) =  - 0.9639719272779138D+00
    x(5) =  - 0.9408226338317548D+00
    x(6) =  - 0.9122344282513259D+00
    x(7) =  - 0.8782768112522820D+00
    x(8) =  - 0.8391169718222188D+00
    x(9) =  - 0.7950414288375512D+00
    x(10) = - 0.7463319064601508D+00
    x(11) = - 0.6932376563347514D+00
    x(12) = - 0.6360536807265150D+00
    x(13) = - 0.5751404468197103D+00
    x(14) = - 0.5108670019508271D+00
    x(15) = - 0.4435931752387251D+00
    x(16) = - 0.3737060887154196D+00
    x(17) = - 0.3016278681149130D+00
    x(18) = - 0.2277858511416451D+00
    x(19) = - 0.1526054652409227D+00
    x(20) = - 0.7652652113349733D-01
    x(21) =   0.0D+00
    x(22) =   0.7652652113349733D-01
    x(23) =   0.1526054652409227D+00
    x(24) =   0.2277858511416451D+00
    x(25) =   0.3016278681149130D+00
    x(26) =   0.3737060887154196D+00
    x(27) =   0.4435931752387251D+00
    x(28) =   0.5108670019508271D+00
    x(29) =   0.5751404468197103D+00
    x(30) =   0.6360536807265150D+00
    x(31) =   0.6932376563347514D+00
    x(32) =   0.7463319064601508D+00
    x(33) =   0.7950414288375512D+00
    x(34) =   0.8391169718222188D+00
    x(35) =   0.8782768112522820D+00
    x(36) =   0.9122344282513259D+00
    x(37) =   0.9408226338317548D+00
    x(38) =   0.9639719272779138D+00
    x(39) =   0.9815078774502503D+00
    x(40) =   0.9931285991850949D+00
    x(41) =   0.9988590315882777D+00

    w(1) =  0.3073583718520532D-02
    w(2) =  0.8600269855642942D-02
    w(3) =  0.1462616925697125D-01
    w(4) =  0.2038837346126652D-01
    w(5) =  0.2588213360495116D-01
    w(6) =  0.3128730677703280D-01
    w(7) =  0.3660016975820080D-01
    w(8) =  0.4166887332797369D-01
    w(9) =  0.4643482186749767D-01
    w(10) = 0.5094457392372869D-01
    w(11) = 0.5519510534828599D-01
    w(12) = 0.5911140088063957D-01
    w(13) = 0.6265323755478117D-01
    w(14) = 0.6583459713361842D-01
    w(15) = 0.6864867292852162D-01
    w(16) = 0.7105442355344407D-01
    w(17) = 0.7303069033278667D-01
    w(18) = 0.7458287540049919D-01
    w(19) = 0.7570449768455667D-01
    w(20) = 0.7637786767208074D-01
    w(21) = 0.7660071191799966D-01
    w(22) = 0.7637786767208074D-01
    w(23) = 0.7570449768455667D-01
    w(24) = 0.7458287540049919D-01
    w(25) = 0.7303069033278667D-01
    w(26) = 0.7105442355344407D-01
    w(27) = 0.6864867292852162D-01
    w(28) = 0.6583459713361842D-01
    w(29) = 0.6265323755478117D-01
    w(30) = 0.5911140088063957D-01
    w(31) = 0.5519510534828599D-01
    w(32) = 0.5094457392372869D-01
    w(33) = 0.4643482186749767D-01
    w(34) = 0.4166887332797369D-01
    w(35) = 0.3660016975820080D-01
    w(36) = 0.3128730677703280D-01
    w(37) = 0.2588213360495116D-01
    w(38) = 0.2038837346126652D-01
    w(39) = 0.1462616925697125D-01
    w(40) = 0.8600269855642942D-02
    w(41) = 0.3073583718520532D-02

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KRONROD_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 15, 21, 31 or 41.'
    stop

  end if

  return
end
subroutine laguerre_ek_compute ( n, x, w )

!*****************************************************************************80
!
!! LAGUERRE_EK_COMPUTE: Laguerre quadrature rule by the Elhay-Kautsky method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = 1.0D+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = 8 )
  end do

  do i = 1, n
    x(i) = real ( 2 * i - 1, kind = 8 )
  end do

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine laguerre_integral ( expon, exact )

!*****************************************************************************80
!
!! LAGUERRE_INTEGRAL evaluates a monomial Laguerre integral.
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_factorial

  exact = r8_factorial ( expon )

  return
end
subroutine laguerre_set ( n, x, w )

!*****************************************************************************80
!
!! LAGUERRE_SET sets abscissas and weights for Laguerre quadrature.
!
!  Discussion:
!
!    The abscissas are the zeroes of the Laguerre polynomial L(N)(X).
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) exp ( -x ) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * exp ( x(i) ) * f ( x(i) )
!
!    Mathematica can numerically estimate the abscissas for the
!    n-th order polynomial to p digits of precision by the command:
!
!      NSolve [ LaguerreL[n,x] == 0, x, p ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    Dover, 2006,
!    ISBN: 0486445798,
!    LC: QA311.K713.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 20, 31/32/33, 63/64/65, 127/128/129.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =  1.00000000000000000000000000000D+00

    w(1) =  1.00000000000000000000000000000D+00

  else if ( n == 2 ) then

    x(1) = 0.585786437626904951198311275790D+00
    x(2) = 3.41421356237309504880168872421D+00

    w(1) = 0.85355339059327376220042218105D+00
    w(2) = 0.146446609406726237799577818948D+00

  else if ( n == 3 ) then

    x(1) = 0.415774556783479083311533873128D+00
    x(2) = 2.29428036027904171982205036136D+00
    x(3) = 6.28994508293747919686641576551D+00

    w(1) = 0.71109300992917301544959019114D+00
    w(2) = 0.27851773356924084880144488846D+00
    w(3) = 0.010389256501586135748964920401D+00

  else if ( n == 4 ) then

    x(1) = 0.322547689619392311800361459104D+00
    x(2) = 1.74576110115834657568681671252D+00
    x(3) = 4.53662029692112798327928538496D+00
    x(4) = 9.39507091230113312923353644342D+00

    w(1) = 0.60315410434163360163596602382D+00
    w(2) = 0.35741869243779968664149201746D+00
    w(3) = 0.03888790851500538427243816816D+00
    w(4) = 0.0005392947055613274501037905676D+00

  else if ( n == 5 ) then

    x(1) = 0.263560319718140910203061943361D+00
    x(2) = 1.41340305910651679221840798019D+00
    x(3) = 3.59642577104072208122318658878D+00
    x(4) = 7.08581000585883755692212418111D+00
    x(5) = 12.6408008442757826594332193066D+00

    w(1) = 0.52175561058280865247586092879D+00
    w(2) = 0.3986668110831759274541333481D+00
    w(3) = 0.0759424496817075953876533114D+00
    w(4) = 0.00361175867992204845446126257D+00
    w(5) = 0.00002336997238577622789114908455D+00

  else if ( n == 6 ) then

    x(1) = 0.222846604179260689464354826787D+00
    x(2) = 1.18893210167262303074315092194D+00
    x(3) = 2.99273632605931407769132528451D+00
    x(4) = 5.77514356910451050183983036943D+00
    x(5) = 9.83746741838258991771554702994D+00
    x(6) = 15.9828739806017017825457915674D+00

    w(1) = 0.45896467394996359356828487771D+00
    w(2) = 0.4170008307721209941133775662D+00
    w(3) = 0.1133733820740449757387061851D+00
    w(4) = 0.01039919745314907489891330285D+00
    w(5) = 0.000261017202814932059479242860D+00
    w(6) = 8.98547906429621238825292053D-07

  else if ( n == 7 ) then

    x(1) = 0.193043676560362413838247885004D+00
    x(2) = 1.02666489533919195034519944317D+00
    x(3) = 2.56787674495074620690778622666D+00
    x(4) = 4.90035308452648456810171437810D+00
    x(5) = 8.18215344456286079108182755123D+00
    x(6) = 12.7341802917978137580126424582D+00
    x(7) = 19.3957278622625403117125820576D+00

    w(1) = 0.40931895170127390213043288002D+00
    w(2) = 0.4218312778617197799292810054D+00
    w(3) = 0.1471263486575052783953741846D+00
    w(4) = 0.0206335144687169398657056150D+00
    w(5) = 0.00107401014328074552213195963D+00
    w(6) = 0.0000158654643485642012687326223D+00
    w(7) = 3.17031547899558056227132215D-08

  else if ( n == 8 ) then

    x(1) = 0.170279632305100999788861856608D+00
    x(2) = 0.903701776799379912186020223555D+00
    x(3) = 2.25108662986613068930711836697D+00
    x(4) = 4.26670017028765879364942182690D+00
    x(5) = 7.04590540239346569727932548212D+00
    x(6) = 10.7585160101809952240599567880D+00
    x(7) = 15.7406786412780045780287611584D+00
    x(8) = 22.8631317368892641057005342974D+00

    w(1) = 0.36918858934163752992058283938D+00
    w(2) = 0.4187867808143429560769785813D+00
    w(3) = 0.175794986637171805699659867D+00
    w(4) = 0.033343492261215651522132535D+00
    w(5) = 0.0027945362352256725249389241D+00
    w(6) = 0.00009076508773358213104238501D+00
    w(7) = 8.4857467162725315448680183D-07
    w(8) = 1.04800117487151038161508854D-09

  else if ( n == 9 ) then

    x(1) = 0.152322227731808247428107073127D+00
    x(2) = 0.807220022742255847741419210952D+00
    x(3) = 2.00513515561934712298303324701D+00
    x(4) = 3.78347397333123299167540609364D+00
    x(5) = 6.20495677787661260697353521006D+00
    x(6) = 9.37298525168757620180971073215D+00
    x(7) = 13.4662369110920935710978818397D+00
    x(8) = 18.8335977889916966141498992996D+00
    x(9) = 26.3740718909273767961410072937D+00

    w(1) = 0.336126421797962519673467717606D+00
    w(2) = 0.411213980423984387309146942793D+00
    w(3) = 0.199287525370885580860575607212D+00
    w(4) = 0.0474605627656515992621163600479D+00
    w(5) = 0.00559962661079458317700419900556D+00
    w(6) = 0.000305249767093210566305412824291D+00
    w(7) = 6.59212302607535239225572284875D-06
    w(8) = 4.1107693303495484429024104033D-08
    w(9) = 3.29087403035070757646681380323D-11

  else if ( n == 10 ) then

    x(1) = 0.137793470540492430830772505653D+00
    x(2) = 0.729454549503170498160373121676D+00
    x(3) = 1.80834290174031604823292007575D+00
    x(4) = 3.40143369785489951448253222141D+00
    x(5) = 5.55249614006380363241755848687D+00
    x(6) = 8.33015274676449670023876719727D+00
    x(7) = 11.8437858379000655649185389191D+00
    x(8) = 16.2792578313781020995326539358D+00
    x(9) = 21.9965858119807619512770901956D+00
    x(10) = 29.9206970122738915599087933408D+00

    w(1) = 0.30844111576502014154747083468D+00
    w(2) = 0.4011199291552735515157803099D+00
    w(3) = 0.218068287611809421588648523D+00
    w(4) = 0.062087456098677747392902129D+00
    w(5) = 0.009501516975181100553839072D+00
    w(6) = 0.0007530083885875387754559644D+00
    w(7) = 0.00002825923349599565567422564D+00
    w(8) = 4.249313984962686372586577D-07
    w(9) = 1.839564823979630780921535D-09
    w(10) = 9.911827219609008558377547D-13

  else if ( n == 11 ) then

    x(1) = 0.125796442187967522675794577516D+00
    x(2) = 0.665418255839227841678127839420D+00
    x(3) = 1.64715054587216930958700321365D+00
    x(4) = 3.09113814303525495330195934259D+00
    x(5) = 5.02928440157983321236999508366D+00
    x(6) = 7.50988786380661681941099714450D+00
    x(7) = 10.6059509995469677805559216457D+00
    x(8) = 14.4316137580641855353200450349D+00
    x(9) = 19.1788574032146786478174853989D+00
    x(10) = 25.2177093396775611040909447797D+00
    x(11) = 33.4971928471755372731917259395D+00

    w(1) = 0.28493321289420060505605102472D+00
    w(2) = 0.3897208895278493779375535080D+00
    w(3) = 0.232781831848991333940223796D+00
    w(4) = 0.076564453546196686400854179D+00
    w(5) = 0.014393282767350695091863919D+00
    w(6) = 0.001518880846484873069847776D+00
    w(7) = 0.0000851312243547192259720424D+00
    w(8) = 2.29240387957450407857683D-06
    w(9) = 2.48635370276779587373391D-08
    w(10) = 7.71262693369132047028153D-11
    w(11) = 2.883775868323623861597778D-14

  else if ( n == 12 ) then

    x(1) = 0.115722117358020675267196428240D+00
    x(2) = 0.611757484515130665391630053042D+00
    x(3) = 1.51261026977641878678173792687D+00
    x(4) = 2.83375133774350722862747177657D+00
    x(5) = 4.59922763941834848460572922485D+00
    x(6) = 6.84452545311517734775433041849D+00
    x(7) = 9.62131684245686704391238234923D+00
    x(8) = 13.0060549933063477203460524294D+00
    x(9) = 17.1168551874622557281840528008D+00
    x(10) = 22.1510903793970056699218950837D+00
    x(11) = 28.4879672509840003125686072325D+00
    x(12) = 37.0991210444669203366389142764D+00

    w(1) = 0.26473137105544319034973889206D+00
    w(2) = 0.3777592758731379820244905567D+00
    w(3) = 0.244082011319877564254870818D+00
    w(4) = 0.09044922221168093072750549D+00
    w(5) = 0.02010238115463409652266129D+00
    w(6) = 0.002663973541865315881054158D+00
    w(7) = 0.000203231592662999392121433D+00
    w(8) = 8.3650558568197987453363D-06
    w(9) = 1.66849387654091026116990D-07
    w(10) = 1.34239103051500414552392D-09
    w(11) = 3.06160163503502078142408D-12
    w(12) = 8.148077467426241682473119D-16

  else if ( n == 13 ) then

    x(1) = 0.107142388472252310648493376977D+00
    x(2) = 0.566131899040401853406036347177D+00
    x(3) = 1.39856433645101971792750259921D+00
    x(4) = 2.61659710840641129808364008472D+00
    x(5) = 4.23884592901703327937303389926D+00
    x(6) = 6.29225627114007378039376523025D+00
    x(7) = 8.81500194118697804733348868036D+00
    x(8) = 11.8614035888112425762212021880D+00
    x(9) = 15.5107620377037527818478532958D+00
    x(10) = 19.8846356638802283332036594634D+00
    x(11) = 25.1852638646777580842970297823D+00
    x(12) = 31.8003863019472683713663283526D+00
    x(13) = 40.7230086692655795658979667001D+00

    w(1) = 0.24718870842996262134624918596D+00
    w(2) = 0.3656888229005219453067175309D+00
    w(3) = 0.252562420057658502356824289D+00
    w(4) = 0.10347075802418370511421863D+00
    w(5) = 0.02643275441556161577815877D+00
    w(6) = 0.00422039604025475276555209D+00
    w(7) = 0.000411881770472734774892473D+00
    w(8) = 0.0000235154739815532386882897D+00
    w(9) = 7.3173116202490991040105D-07
    w(10) = 1.10884162570398067979151D-08
    w(11) = 6.7708266922058988406462D-11
    w(12) = 1.15997995990507606094507D-13
    w(13) = 2.245093203892758415991872D-17

  else if ( n == 14 ) then

    x(1) = 0.0997475070325975745736829452514D+00
    x(2) = 0.526857648851902896404583451502D+00
    x(3) = 1.30062912125149648170842022116D+00
    x(4) = 2.43080107873084463616999751038D+00
    x(5) = 3.93210282229321888213134366778D+00
    x(6) = 5.82553621830170841933899983898D+00
    x(7) = 8.14024014156514503005978046052D+00
    x(8) = 10.9164995073660188408130510904D+00
    x(9) = 14.2108050111612886831059780825D+00
    x(10) = 18.1048922202180984125546272083D+00
    x(11) = 22.7233816282696248232280886985D+00
    x(12) = 28.2729817232482056954158923218D+00
    x(13) = 35.1494436605924265828643121364D+00
    x(14) = 44.3660817111174230416312423666D+00

    w(1) = 0.23181557714486497784077486110D+00
    w(2) = 0.3537846915975431518023313013D+00
    w(3) = 0.258734610245428085987320561D+00
    w(4) = 0.11548289355692321008730499D+00
    w(5) = 0.03319209215933736003874996D+00
    w(6) = 0.00619286943700661021678786D+00
    w(7) = 0.00073989037786738594242589D+00
    w(8) = 0.000054907194668416983785733D+00
    w(9) = 2.4095857640853774967578D-06
    w(10) = 5.801543981676495180886D-08
    w(11) = 6.819314692484974119616D-10
    w(12) = 3.2212077518948479398089D-12
    w(13) = 4.2213524405165873515980D-15
    w(14) = 6.05237502228918880839871D-19

  else if ( n == 15 ) then

    x(1) = 0.0933078120172818047629030383672D+00
    x(2) = 0.492691740301883908960101791412D+00
    x(3) = 1.21559541207094946372992716488D+00
    x(4) = 2.26994952620374320247421741375D+00
    x(5) = 3.66762272175143727724905959436D+00
    x(6) = 5.42533662741355316534358132596D+00
    x(7) = 7.56591622661306786049739555812D+00
    x(8) = 10.1202285680191127347927394568D+00
    x(9) = 13.1302824821757235640991204176D+00
    x(10) = 16.6544077083299578225202408430D+00
    x(11) = 20.7764788994487667729157175676D+00
    x(12) = 25.6238942267287801445868285977D+00
    x(13) = 31.4075191697539385152432196202D+00
    x(14) = 38.5306833064860094162515167595D+00
    x(15) = 48.0260855726857943465734308508D+00

    w(1) = 0.21823488594008688985641323645D+00
    w(2) = 0.3422101779228833296389489568D+00
    w(3) = 0.263027577941680097414812275D+00
    w(4) = 0.12642581810593053584303055D+00
    w(5) = 0.04020686492100091484158548D+00
    w(6) = 0.00856387780361183836391576D+00
    w(7) = 0.00121243614721425207621921D+00
    w(8) = 0.00011167439234425194199258D+00
    w(9) = 6.459926762022900924653D-06
    w(10) = 2.226316907096272630332D-07
    w(11) = 4.227430384979365007351D-09
    w(12) = 3.921897267041089290385D-11
    w(13) = 1.4565152640731264063327D-13
    w(14) = 1.4830270511133013354616D-16
    w(15) = 1.60059490621113323104998D-20

  else if ( n == 16 ) then

    x(1) = 0.0876494104789278403601980973401D+00
    x(2) = 0.462696328915080831880838260664D+00
    x(3) = 1.14105777483122685687794501811D+00
    x(4) = 2.12928364509838061632615907066D+00
    x(5) = 3.43708663389320664523510701675D+00
    x(6) = 5.07801861454976791292305830814D+00
    x(7) = 7.07033853504823413039598947080D+00
    x(8) = 9.43831433639193878394724672911D+00
    x(9) = 12.2142233688661587369391246088D+00
    x(10) = 15.4415273687816170767647741622D+00
    x(11) = 19.1801568567531348546631409497D+00
    x(12) = 23.5159056939919085318231872752D+00
    x(13) = 28.5787297428821403675206137099D+00
    x(14) = 34.5833987022866258145276871778D+00
    x(15) = 41.9404526476883326354722330252D+00
    x(16) = 51.7011603395433183643426971197D+00

    w(1) = 0.20615171495780099433427363674D+00
    w(2) = 0.3310578549508841659929830987D+00
    w(3) = 0.265795777644214152599502021D+00
    w(4) = 0.13629693429637753997554751D+00
    w(5) = 0.0473289286941252189780623D+00
    w(6) = 0.0112999000803394532312490D+00
    w(7) = 0.0018490709435263108642918D+00
    w(8) = 0.00020427191530827846012602D+00
    w(9) = 0.00001484458687398129877135D+00
    w(10) = 6.828319330871199564396D-07
    w(11) = 1.881024841079673213882D-08
    w(12) = 2.862350242973881619631D-10
    w(13) = 2.127079033224102967390D-12
    w(14) = 6.297967002517867787174D-15
    w(15) = 5.050473700035512820402D-18
    w(16) = 4.1614623703728551904265D-22

  else if ( n == 17 ) then

    x(1) = 0.0826382147089476690543986151980D+00
    x(2) = 0.436150323558710436375959029847D+00
    x(3) = 1.07517657751142857732980316755D+00
    x(4) = 2.00519353164923224070293371933D+00
    x(5) = 3.23425612404744376157380120696D+00
    x(6) = 4.77351351370019726480932076262D+00
    x(7) = 6.63782920536495266541643929703D+00
    x(8) = 8.84668551116980005369470571184D+00
    x(9) = 11.4255293193733525869726151469D+00
    x(10) = 14.4078230374813180021982874959D+00
    x(11) = 17.8382847307011409290658752412D+00
    x(12) = 21.7782682577222653261749080522D+00
    x(13) = 26.3153178112487997766149598369D+00
    x(14) = 31.5817716804567331343908517497D+00
    x(15) = 37.7960938374771007286092846663D+00
    x(16) = 45.3757165339889661829258363215D+00
    x(17) = 55.3897517898396106640900199790D+00

    w(1) = 0.19533220525177083214592729770D+00
    w(2) = 0.3203753572745402813366256320D+00
    w(3) = 0.267329726357171097238809604D+00
    w(4) = 0.14512985435875862540742645D+00
    w(5) = 0.0544369432453384577793806D+00
    w(6) = 0.0143572977660618672917767D+00
    w(7) = 0.0026628247355727725684324D+00
    w(8) = 0.0003436797271562999206118D+00
    w(9) = 0.00003027551783782870109437D+00
    w(10) = 1.768515053231676895381D-06
    w(11) = 6.57627288681043332199D-08
    w(12) = 1.469730932159546790344D-09
    w(13) = 1.81691036255544979555D-11
    w(14) = 1.095401388928687402976D-13
    w(15) = 2.617373882223370421551D-16
    w(16) = 1.6729356931461546908502D-19
    w(17) = 1.06562631627404278815253D-23

  else if ( n == 18 ) then

    x(1) = 0.0781691666697054712986747615334D+00
    x(2) = 0.412490085259129291039101536536D+00
    x(3) = 1.01652017962353968919093686187D+00
    x(4) = 1.89488850996976091426727831954D+00
    x(5) = 3.05435311320265975115241130719D+00
    x(6) = 4.50420553888989282633795571455D+00
    x(7) = 6.25672507394911145274209116326D+00
    x(8) = 8.32782515660563002170470261564D+00
    x(9) = 10.7379900477576093352179033397D+00
    x(10) = 13.5136562075550898190863812108D+00
    x(11) = 16.6893062819301059378183984163D+00
    x(12) = 20.3107676262677428561313764553D+00
    x(13) = 24.4406813592837027656442257980D+00
    x(14) = 29.1682086625796161312980677805D+00
    x(15) = 34.6279270656601721454012429438D+00
    x(16) = 41.0418167728087581392948614284D+00
    x(17) = 48.8339227160865227486586093290D+00
    x(18) = 59.0905464359012507037157810181D+00

    w(1) = 0.18558860314691880562333775228D+00
    w(2) = 0.3101817663702252936495975957D+00
    w(3) = 0.267866567148536354820854395D+00
    w(4) = 0.15297974746807490655384308D+00
    w(5) = 0.0614349178609616527076780D+00
    w(6) = 0.0176872130807729312772600D+00
    w(7) = 0.0036601797677599177980266D+00
    w(8) = 0.0005406227870077353231284D+00
    w(9) = 0.0000561696505121423113818D+00
    w(10) = 4.01530788370115755859D-06
    w(11) = 1.91466985667567497969D-07
    w(12) = 5.8360952686315941292D-09
    w(13) = 1.07171126695539012773D-10
    w(14) = 1.08909871388883385562D-12
    w(15) = 5.38666474837830887608D-15
    w(16) = 1.049865978035703408779D-17
    w(17) = 5.405398451631053643566D-21
    w(18) = 2.6916532692010286270838D-25

  else if ( n == 19 ) then

    x(1) = 0.0741587837572050877131369916024D+00
    x(2) = 0.391268613319994607337648350299D+00
    x(3) = 0.963957343997958058624878377130D+00
    x(4) = 1.79617558206832812557725825252D+00
    x(5) = 2.89365138187378399116494713237D+00
    x(6) = 4.26421553962776647436040018167D+00
    x(7) = 5.91814156164404855815360191408D+00
    x(8) = 7.86861891533473373105668358176D+00
    x(9) = 10.1324237168152659251627415800D+00
    x(10) = 12.7308814638423980045092979656D+00
    x(11) = 15.6912783398358885454136069861D+00
    x(12) = 19.0489932098235501532136429732D+00
    x(13) = 22.8508497608294829323930586693D+00
    x(14) = 27.1606693274114488789963947149D+00
    x(15) = 32.0691222518622423224362865906D+00
    x(16) = 37.7129058012196494770647508283D+00
    x(17) = 44.3173627958314961196067736013D+00
    x(18) = 52.3129024574043831658644222420D+00
    x(19) = 62.8024231535003758413504690673D+00

    w(1) = 0.17676847491591250225103547981D+00
    w(2) = 0.3004781436072543794821568077D+00
    w(3) = 0.267599547038175030772695441D+00
    w(4) = 0.15991337213558021678551215D+00
    w(5) = 0.0682493799761491134552355D+00
    w(6) = 0.0212393076065443249244062D+00
    w(7) = 0.0048416273511483959672501D+00
    w(8) = 0.0008049127473813667665946D+00
    w(9) = 0.0000965247209315350170843D+00
    w(10) = 8.20730525805103054409D-06
    w(11) = 4.8305667247307725394D-07
    w(12) = 1.90499136112328569994D-08
    w(13) = 4.8166846309280615577D-10
    w(14) = 7.3482588395511443768D-12
    w(15) = 6.2022753875726163989D-14
    w(16) = 2.54143084301542272372D-16
    w(17) = 4.07886129682571235007D-19
    w(18) = 1.707750187593837061004D-22
    w(19) = 6.715064649908189959990D-27

  else if ( n == 20 ) then

    x(1) = 0.0705398896919887533666890045842D+00
    x(2) = 0.372126818001611443794241388761D+00
    x(3) = 0.916582102483273564667716277074D+00
    x(4) = 1.70730653102834388068768966741D+00
    x(5) = 2.74919925530943212964503046049D+00
    x(6) = 4.04892531385088692237495336913D+00
    x(7) = 5.61517497086161651410453988565D+00
    x(8) = 7.45901745367106330976886021837D+00
    x(9) = 9.59439286958109677247367273428D+00
    x(10) = 12.0388025469643163096234092989D+00
    x(11) = 14.8142934426307399785126797100D+00
    x(12) = 17.9488955205193760173657909926D+00
    x(13) = 21.4787882402850109757351703696D+00
    x(14) = 25.4517027931869055035186774846D+00
    x(15) = 29.9325546317006120067136561352D+00
    x(16) = 35.0134342404790000062849359067D+00
    x(17) = 40.8330570567285710620295677078D+00
    x(18) = 47.6199940473465021399416271529D+00
    x(19) = 55.8107957500638988907507734445D+00
    x(20) = 66.5244165256157538186403187915D+00

    w(1) = 0.168746801851113862149223899689D+00
    w(2) = 0.291254362006068281716795323812D+00
    w(3) = 0.266686102867001288549520868998D+00
    w(4) = 0.166002453269506840031469127816D+00
    w(5) = 0.0748260646687923705400624639615D+00
    w(6) = 0.0249644173092832210728227383234D+00
    w(7) = 0.00620255084457223684744754785395D+00
    w(8) = 0.00114496238647690824203955356969D+00
    w(9) = 0.000155741773027811974779809513214D+00
    w(10) = 0.0000154014408652249156893806714048D+00
    w(11) = 1.08648636651798235147970004439D-06
    w(12) = 5.33012090955671475092780244305D-08
    w(13) = 1.7579811790505820035778763784D-09
    w(14) = 3.72550240251232087262924585338D-11
    w(15) = 4.76752925157819052449488071613D-13
    w(16) = 3.37284424336243841236506064991D-15
    w(17) = 1.15501433950039883096396247181D-17
    w(18) = 1.53952214058234355346383319667D-20
    w(19) = 5.28644272556915782880273587683D-24
    w(20) = 1.65645661249902329590781908529D-28

  else if ( n == 31 ) then

    x(1) = 0.0459019476211082907434960802752D+00
    x(2) = 0.241980163824772048904089741517D+00
    x(3) = 0.595253894222350737073301650054D+00
    x(4) = 1.10668949953299871621113087898D+00
    x(5) = 1.77759569287477272115937274827D+00
    x(6) = 2.60970341525668065038933759253D+00
    x(7) = 3.60519680234004426988058175542D+00
    x(8) = 4.76674708447176113136291272711D+00
    x(9) = 6.09755456718174092699254293285D+00
    x(10) = 7.60140094923313742293601069429D+00
    x(11) = 9.28271431347088941825366952977D+00
    x(12) = 11.1466497556192913589938156296D+00
    x(13) = 13.1991895762449985224649250286D+00
    x(14) = 15.4472683155493100758093258918D+00
    x(15) = 17.8989298266447576467257938178D+00
    x(16) = 20.5635263367158221707430489688D+00
    x(17) = 23.4519734820118585910502555759D+00
    x(18) = 26.5770813521182604599758769865D+00
    x(19) = 29.9539908723464455069519178400D+00
    x(20) = 33.6007595329022027354103138858D+00
    x(21) = 37.5391644073304408828879025580D+00
    x(22) = 41.7958308701822199813479458533D+00
    x(23) = 46.4038668064111231360292276044D+00
    x(24) = 51.4053144767977551618614610884D+00
    x(25) = 56.8549928687158436205119220557D+00
    x(26) = 62.8268559087863214536775233048D+00
    x(27) = 69.4252771910803456233222516564D+00
    x(28) = 76.8070477638627328376099722855D+00
    x(29) = 85.2303586075456691693870656070D+00
    x(30) = 95.1889398915256299813086068540D+00
    x(31) = 107.952243827578714750024401177D+00

    w(1) = 0.112527895503725838208477280828D+00
    w(2) = 0.21552760818089123795222505285D+00
    w(3) = 0.238308251645696547319057880892D+00
    w(4) = 0.195388309297902292499153033907D+00
    w(5) = 0.126982832893061901436352729046D+00
    w(6) = 0.0671861689238993006709294419935D+00
    w(7) = 0.029303224993879487404888669312D+00
    w(8) = 0.0105975699152957360895293803144D+00
    w(9) = 0.0031851272582386980320974842433D+00
    w(10) = 0.000795495483079403829220921490125D+00
    w(11) = 0.000164800521266366873178629671164D+00
    w(12) = 0.000028229237864310816393860971469D+00
    w(13) = 3.98029025510085803871161749001D-06
    w(14) = 4.59318398418010616737296945103D-07
    w(15) = 4.30755451877311009301314574659D-08
    w(16) = 3.25512499382715708551757492579D-09
    w(17) = 1.96202466754105949962471515931D-10
    w(18) = 9.31904990866175871295347164313D-12
    w(19) = 3.43775418194116205203125978983D-13
    w(20) = 9.67952471304467169974050357762D-15
    w(21) = 2.03680661101152473980106242193D-16
    w(22) = 3.12126872807135268317653586326D-18
    w(23) = 3.37295817041610524533956783083D-20
    w(24) = 2.46727963866166960110383632425D-22
    w(25) = 1.15822019045256436348345645766D-24
    w(26) = 3.2472922591425422434798022809D-27
    w(27) = 4.91430173080574327408200762597D-30
    w(28) = 3.45000711048083941322231359538D-33
    w(29) = 8.76637101171620414729327607329D-37
    w(30) = 5.03636439211614904112971723166D-41
    w(31) = 1.99099845825314564824395490803D-46

  else if ( n == 32 ) then

    x(1) = 0.0444893658332670184188501945244D+00
    x(2) = 0.234526109519618537452909561302D+00
    x(3) = 0.576884629301886426491552569378D+00
    x(4) = 1.07244875381781763304091397718D+00
    x(5) = 1.72240877644464544113093292797D+00
    x(6) = 2.52833670642579488112419990556D+00
    x(7) = 3.49221327302199448960880339977D+00
    x(8) = 4.61645676974976738776205229617D+00
    x(9) = 5.90395850417424394656152149158D+00
    x(10) = 7.35812673318624111322198973719D+00
    x(11) = 8.98294092421259610337824752677D+00
    x(12) = 10.7830186325399720675019491381D+00
    x(13) = 12.7636979867427251149690330822D+00
    x(14) = 14.9311397555225573197969646873D+00
    x(15) = 17.2924543367153147892357183836D+00
    x(16) = 19.8558609403360547397899445841D+00
    x(17) = 22.6308890131967744886757793394D+00
    x(18) = 25.6286360224592477674761768768D+00
    x(19) = 28.8621018163234747443426407115D+00
    x(20) = 32.3466291539647370032321654237D+00
    x(21) = 36.1004948057519738040171189479D+00
    x(22) = 40.1457197715394415362093439289D+00
    x(23) = 44.5092079957549379759066043775D+00
    x(24) = 49.2243949873086391767222218066D+00
    x(25) = 54.3337213333969073328671815512D+00
    x(26) = 59.8925091621340181961304753247D+00
    x(27) = 65.9753772879350527965630761193D+00
    x(28) = 72.6876280906627086386753490878D+00
    x(29) = 80.1874469779135230674916385687D+00
    x(30) = 88.7353404178923986893554495243D+00
    x(31) = 98.8295428682839725591844784095D+00
    x(32) = 111.751398097937695213664716539D+00

    w(1) = 0.109218341952384971136131337943D+00
    w(2) = 0.210443107938813232936062071392D+00
    w(3) = 0.235213229669848005394941106676D+00
    w(4) = 0.195903335972881043413247901182D+00
    w(5) = 0.129983786286071760607216822326D+00
    w(6) = 0.0705786238657174415601643320433D+00
    w(7) = 0.0317609125091750703058255211629D+00
    w(8) = 0.0119182148348385570565446505756D+00
    w(9) = 0.00373881629461152478966122847796D+00
    w(10) = 0.000980803306614955132230630611308D+00
    w(11) = 0.000214864918801364188023199483686D+00
    w(12) = 0.0000392034196798794720432695682782D+00
    w(13) = 5.93454161286863287835582893773D-06
    w(14) = 7.4164045786675522190708220213D-07
    w(15) = 7.60456787912078148111926545943D-08
    w(16) = 6.35060222662580674242777108552D-09
    w(17) = 4.28138297104092887881360582098D-10
    w(18) = 2.30589949189133607927336809618D-11
    w(19) = 9.79937928872709406333455225995D-13
    w(20) = 3.23780165772926646231042646142D-14
    w(21) = 8.17182344342071943320186059177D-16
    w(22) = 1.54213383339382337217855949129D-17
    w(23) = 2.11979229016361861204093474373D-19
    w(24) = 2.05442967378804542665570987602D-21
    w(25) = 1.34698258663739515580519340478D-23
    w(26) = 5.66129413039735937112634432382D-26
    w(27) = 1.41856054546303690595142933892D-28
    w(28) = 1.91337549445422430937127829683D-31
    w(29) = 1.19224876009822235654164532831D-34
    w(30) = 2.67151121924013698599867893958D-38
    w(31) = 1.33861694210625628271905701423D-42
    w(32) = 4.51053619389897423222342830132D-48

  else if ( n == 33 ) then

    x(1) = 0.0431611356173268921917334738206D+00
    x(2) = 0.227517802803371123850290226913D+00
    x(3) = 0.559616655851539887586282303916D+00
    x(4) = 1.04026850775100205382209621927D+00
    x(5) = 1.67055919607571519092562973257D+00
    x(6) = 2.45192079589763054651073898192D+00
    x(7) = 3.38615533758800483230187851832D+00
    x(8) = 4.47545949839977145702059137905D+00
    x(9) = 5.72245472027210352266790817933D+00
    x(10) = 7.13022434440010801631414039534D+00
    x(11) = 8.70235923062140624893696399459D+00
    x(12) = 10.4430136502059824268455293839D+00
    x(13) = 12.3569737593502859624441255236D+00
    x(14) = 14.4497416815855402377145121178D+00
    x(15) = 16.7276392186383223532615229942D+00
    x(16) = 19.1979365872124466372283088222D+00
    x(17) = 21.8690135249281898713512287043D+00
    x(18) = 24.7505629061577956433730931987D+00
    x(19) = 27.8538511114133567797747375537D+00
    x(20) = 31.1920555455751298677734295989D+00
    x(21) = 34.7807091535383377002292521853D+00
    x(22) = 38.6382967177740302250360622751D+00
    x(23) = 42.7870720782534794879639219927D+00
    x(24) = 47.2542066029932658172690829767D+00
    x(25) = 52.0734519015142202671640200482D+00
    x(26) = 57.2876345410929400754514841078D+00
    x(27) = 62.9525659469066302071906336861D+00
    x(28) = 69.1435133801098924457366348147D+00
    x(29) = 75.9666870142470623437939790250D+00
    x(30) = 83.5816372232708807614192336050D+00
    x(31) = 92.2511394441351012341481184391D+00
    x(32) = 102.477844336823322575825984750D+00
    x(33) = 115.554756448995807306876850793D+00

    w(1) = 0.106097745553686759448980241986D+00
    w(2) = 0.205582983661932603502389046086D+00
    w(3) = 0.232126523496060850848680143719D+00
    w(4) = 0.196207372769141916829837191484D+00
    w(5) = 0.132744856705171098099698375677D+00
    w(6) = 0.0738518038877138714048058524075D+00
    w(7) = 0.0342232334108640270351258175641D+00
    w(8) = 0.0132939751808086665861981468532D+00
    w(9) = 0.00434094309504623645941229723703D+00
    w(10) = 0.0011922509906686840510776728263D+00
    w(11) = 0.000275158225582396584420253012954D+00
    w(12) = 0.0000532433409922782412444424323192D+00
    w(13) = 8.60957132646331618369236329569D-06
    w(14) = 1.15837796102469195604266695945D-06
    w(15) = 1.28985114856525884538052927779D-07
    w(16) = 1.18096786980325241031580363325D-08
    w(17) = 8.82276640967020246770192412272D-10
    w(18) = 5.32961213410302701363055555216D-11
    w(19) = 2.57550403748317439431393144398D-12
    w(20) = 9.8314133225207825980561863437D-14
    w(21) = 2.92041495556546845392792035892D-15
    w(22) = 6.63077156752381637730149511056D-17
    w(23) = 1.12609863704995018019882580368D-18
    w(24) = 1.39311657122392009414616980902D-20
    w(25) = 1.21481009891544297673141523063D-22
    w(26) = 7.16158181142099840535743381278D-25
    w(27) = 2.70320712488116872172720734689D-27
    w(28) = 6.07192361286922243586897316027D-30
    w(29) = 7.32134211132579407517616095945D-33
    w(30) = 4.06131706145569795511645700604D-36
    w(31) = 8.04952284545203726871355981553D-40
    w(32) = 3.52902990360469937522008596417D-44
    w(33) = 1.01716656299412569799194166119D-49

  else if ( n == 63 ) then

    x(1) = 0.0227688937325761537859943302486D+00
    x(2) = 0.119983252427278247157714164264D+00
    x(3) = 0.294941854447701495774277385174D+00
    x(4) = 0.547790878962377253638650737759D+00
    x(5) = 0.878690611799319016738955670523D+00
    x(6) = 1.28784643359197063023092077886D+00
    x(7) = 1.77551238153885537639794632687D+00
    x(8) = 2.34199255670859892560556283377D+00
    x(9) = 2.98764232232464739399767310536D+00
    x(10) = 3.71286959920180003462996374134D+00
    x(11) = 4.51813633495035843911055685616D+00
    x(12) = 5.40396017818259462869025997827D+00
    x(13) = 6.37091637878653302203922508918D+00
    x(14) = 7.41963993393117111548884931990D+00
    x(15) = 8.55082800084033283125890487222D+00
    x(16) = 9.76524259992453668070045929780D+00
    x(17) = 11.0637136351406617362205504106D+00
    x(18) = 12.4471422623564927497986875693D+00
    x(19) = 13.9165046410578185629129670082D+00
    x(20) = 15.4728561100362964247771436078D+00
    x(21) = 17.1173358338635887531169003039D+00
    x(22) = 18.8511719741548568508734837875D+00
    x(23) = 20.6756874480565156603772656674D+00
    x(24) = 22.5923063463115283812922777600D+00
    x(25) = 24.6025610949726388837006427600D+00
    x(26) = 26.7081004587373439697790879988D+00
    x(27) = 28.9106985004513826401777181032D+00
    x(28) = 31.2122646311759128854777738208D+00
    x(29) = 33.6148549091011548365988428883D+00
    x(30) = 36.1206847744848230563063287408D+00
    x(31) = 38.7321434429335821456260416077D+00
    x(32) = 41.4518102223187411911147261814D+00
    x(33) = 44.2824730714792338393588571346D+00
    x(34) = 47.2271497842956868989350952315D+00
    x(35) = 50.2891122642406957617490218394D+00
    x(36) = 53.4719144567886528083482806195D+00
    x(37) = 56.7794246363420622130997810571D+00
    x(38) = 60.2158629090198628864175501144D+00
    x(39) = 63.7858450042359746317011396018D+00
    x(40) = 67.4944337022938858303743256950D+00
    x(41) = 71.3471996042952662866548033761D+00
    x(42) = 75.3502934256532342542905047443D+00
    x(43) = 79.5105326299863091495553913548D+00
    x(44) = 83.8355060808722578433398176585D+00
    x(45) = 88.3337015703543690861127663265D+00
    x(46) = 93.0146627285585474053033990371D+00
    x(47) = 97.8891841475781400433867276771D+00
    x(48) = 102.969556907413816507839527468D+00
    x(49) = 108.269881619615953922263509672D+00
    x(50) = 113.806473502874627389344859559D+00
    x(51) = 119.598395388304586669624529633D+00
    x(52) = 125.668172558561194312911963033D+00
    x(53) = 132.042772720911657465855905830D+00
    x(54) = 138.754984181037890781675905675D+00
    x(55) = 145.845413183135403582839942484D+00
    x(56) = 153.365484594978636237108159627D+00
    x(57) = 161.382151948137612435621726696D+00
    x(58) = 169.985706006658394387951753012D+00
    x(59) = 179.303662474015809102518278585D+00
    x(60) = 189.527895965324754736687213330D+00
    x(61) = 200.975211599246567416286718410D+00
    x(62) = 214.253685366387886426980562964D+00
    x(63) = 230.934657470897039712465629851D+00

    w(1) = 0.0571186332138689798115872833905D+00
    w(2) = 0.120674760906403952833199320363D+00
    w(3) = 0.159250010965818737238705610965D+00
    w(4) = 0.168751783275607992345961929636D+00
    w(5) = 0.153666419776689566961937113101D+00
    w(6) = 0.123687706147164816410866522619D+00
    w(7) = 0.0892750988548486715452791500574D+00
    w(8) = 0.0582584854461059449575718257252D+00
    w(9) = 0.0345466575459925808747170858125D+00
    w(10) = 0.0186756859857146567982865525912D+00
    w(11) = 0.00922334490440935365284900752416D+00
    w(12) = 0.00416712506848395927625826634702D+00
    w(13) = 0.0017238120299900582715386728542D+00
    w(14) = 0.00065320845029716311169340559359D+00
    w(15) = 0.000226776446709095869524051732075D+00
    w(16) = 0.0000721276741548106684107502702349D+00
    w(17) = 0.0000210112611804664845988115368512D+00
    w(18) = 5.60355008933572127491815360713D-06
    w(19) = 1.3673642785604888017836641282D-06
    w(20) = 3.050726393019581724073609719D-07
    w(21) = 6.21800618393097635599817754D-08
    w(22) = 1.1566529551931711260022449D-08
    w(23) = 1.9614588267565478081534782D-09
    w(24) = 3.028617119570941124433476D-10
    w(25) = 4.252134453940068676901296D-11
    w(26) = 5.42022205780738193346988D-12
    w(27) = 6.2627306838597672554167D-13
    w(28) = 6.5474443156573322992307D-14
    w(29) = 6.18155758087291818463D-15
    w(30) = 5.259272136350738140426D-16
    w(31) = 4.023092009264648401539D-17
    w(32) = 2.7600740511819536505D-18
    w(33) = 1.69369467569682960533D-19
    w(34) = 9.2689146872177087315D-21
    w(35) = 4.509373906036563294D-22
    w(36) = 1.9435162876132376574D-23
    w(37) = 7.392627089516920704D-25
    w(38) = 2.471436415443463262D-26
    w(39) = 7.228864944674159766D-28
    w(40) = 1.840761729261403936D-29
    w(41) = 4.058349856684196011D-31
    w(42) = 7.70004964164383681D-33
    w(43) = 1.248850576499933433D-34
    w(44) = 1.71850002267670107D-36
    w(45) = 1.989637263667239694D-38
    w(46) = 1.919967137880405827D-40
    w(47) = 1.527858828552216692D-42
    w(48) = 9.9054752688842143D-45
    w(49) = 5.159752367302921188D-47
    w(50) = 2.124984666408411125D-49
    w(51) = 6.790385276685291059D-52
    w(52) = 1.6466654148296177468D-54
    w(53) = 2.9509065402691055027D-57
    w(54) = 3.78384206475710519849D-60
    w(55) = 3.33581300685424318782D-63
    w(56) = 1.922346102227388098136D-66
    w(57) = 6.781269696108301687278D-70
    w(58) = 1.34047528024406046076205D-73
    w(59) = 1.3109745101805029757648D-77
    w(60) = 5.262486388140178738869458D-82
    w(61) = 6.3780013856587414257760666D-87
    w(62) = 1.299707894237292456634747392D-92
    w(63) = 1.00085114969687540634437401684D-99

  else if ( n == 64 ) then

    x(1) = 0.0224158741467052800228118848190D+00
    x(2) = 0.118122512096770479797466436710D+00
    x(3) = 0.290365744018036483999130066385D+00
    x(4) = 0.539286221227979039318144947812D+00
    x(5) = 0.865037004648113944619955074710D+00
    x(6) = 1.26781404077524139811570887769D+00
    x(7) = 1.74785962605943625282996395129D+00
    x(8) = 2.30546373930750871854807054389D+00
    x(9) = 2.94096515672525184067946815211D+00
    x(10) = 3.65475265020729052703539791209D+00
    x(11) = 4.44726634331309435674255098016D+00
    x(12) = 5.31899925449639034352210985919D+00
    x(13) = 6.27049904692365391291106464633D+00
    x(14) = 7.30237000258739574722349840952D+00
    x(15) = 8.41527523948302419449521859120D+00
    x(16) = 9.60993919279610803576288204955D+00
    x(17) = 10.8871503838863721425945504202D+00
    x(18) = 12.2477645042443016181623692907D+00
    x(19) = 13.6927078455475051527299325746D+00
    x(20) = 15.2229811115247288480082687834D+00
    x(21) = 16.8396636526487372105288380392D+00
    x(22) = 18.5439181708591905236196259711D+00
    x(23) = 20.3369959487302355011498158064D+00
    x(24) = 22.2202426659508765399221543471D+00
    x(25) = 24.1951048759332539898864438802D+00
    x(26) = 26.2631372271184857851260239548D+00
    x(27) = 28.4260105275010272994997715268D+00
    x(28) = 30.6855207675259717710485823984D+00
    x(29) = 33.0435992364378291255202106805D+00
    x(30) = 35.5023238911412095869787785351D+00
    x(31) = 38.0639321656464682603573179150D+00
    x(32) = 40.7308354444586263657318695132D+00
    x(33) = 43.5056354664215298527031849317D+00
    x(34) = 46.3911429786161920736053999424D+00
    x(35) = 49.3903990256246866792358008227D+00
    x(36) = 52.5066993413463016501792769805D+00
    x(37) = 55.7436224132783804633357112912D+00
    x(38) = 59.1050619190171066088487420918D+00
    x(39) = 62.5952644001513955960550179012D+00
    x(40) = 66.2188732512475643822137626710D+00
    x(41) = 69.9809803771468292285346579722D+00
    x(42) = 73.8871872324829632109574031135D+00
    x(43) = 77.9436774344631203136879758706D+00
    x(44) = 82.1573037783193042951958683422D+00
    x(45) = 86.5356933494565182102162783753D+00
    x(46) = 91.0873756131330901456367153493D+00
    x(47) = 95.8219400155207320947672154365D+00
    x(48) = 100.750231969513979629259261451D+00
    x(49) = 105.884599468799949356360427851D+00
    x(50) = 111.239207524439582063484736638D+00
    x(51) = 116.830445051306498463386669077D+00
    x(52) = 122.677460268538576577419690565D+00
    x(53) = 128.802878769237672512753623054D+00
    x(54) = 135.233787949525827833980498879D+00
    x(55) = 142.003121489931519025140038291D+00
    x(56) = 149.151665900049388587293462932D+00
    x(57) = 156.731075132671161233616960814D+00
    x(58) = 164.808602655150522993190109025D+00
    x(59) = 173.474946836424274522152844867D+00
    x(60) = 182.858204691431463646342794510D+00
    x(61) = 193.151136037072911479385527417D+00
    x(62) = 204.672028485059455949064433343D+00
    x(63) = 218.031851935328516332452384448D+00
    x(64) = 234.809579171326164713055529725D+00

    w(1) = 0.0562528423390298457410218545063D+00
    w(2) = 0.119023987312426027814903505889D+00
    w(3) = 0.157496403862144523820196434706D+00
    w(4) = 0.167547050415773947880904411659D+00
    w(5) = 0.153352855779236618085454792564D+00
    w(6) = 0.124221053609329744512613782193D+00
    w(7) = 0.0903423009864850577389741092016D+00
    w(8) = 0.0594777557683550242122469974397D+00
    w(9) = 0.0356275189040360718541657353369D+00
    w(10) = 0.0194804104311664060433373802715D+00
    w(11) = 0.00974359489938200224010796027138D+00
    w(12) = 0.00446431036416627529236482419054D+00
    w(13) = 0.00187535958132311482675012822252D+00
    w(14) = 0.000722646981575005122719108706619D+00
    w(15) = 0.00025548753283349670971444840218D+00
    w(16) = 0.0000828714353439694217906322403988D+00
    w(17) = 0.0000246568639678855874597337022172D+00
    w(18) = 6.7267138788296685276125455363D-06
    w(19) = 1.6817853699640888978212010865D-06
    w(20) = 3.850812981546684414827886111D-07
    w(21) = 8.06872804099049979041511092D-08
    w(22) = 1.54572370675768882800370393D-08
    w(23) = 2.7044801476174814099886074D-09
    w(24) = 4.316775475427200912314164D-10
    w(25) = 6.27775254176145220165296D-11
    w(26) = 8.30631737628895806387893D-12
    w(27) = 9.9840317872201640558973D-13
    w(28) = 1.0883538871166626853261D-13
    w(29) = 1.0740174034415901864828D-14
    w(30) = 9.57573723157444210559D-16
    w(31) = 7.69702802364858609886D-17
    w(32) = 5.56488113745402536653D-18
    w(33) = 3.6097564090104464983D-19
    w(34) = 2.0950953695489462348D-20
    w(35) = 1.0847933010975493612D-21
    w(36) = 4.994699486363804116D-23
    w(37) = 2.037836974598822311D-24
    w(38) = 7.339537564278837039D-26
    w(39) = 2.323783082198694261D-27
    w(40) = 6.43823470690876242D-29
    w(41) = 1.553121095788275271D-30
    w(42) = 3.24425009201953731D-32
    w(43) = 5.8323862678362015D-34
    w(44) = 8.96325483310285406D-36
    w(45) = 1.168703989550736241D-37
    w(46) = 1.282055984359980381D-39
    w(47) = 1.172094937405002292D-41
    w(48) = 8.83533967232860498D-44
    w(49) = 5.42495559030618659D-46
    w(50) = 2.675542666678893829D-48
    w(51) = 1.042917031411367078D-50
    w(52) = 3.152902351957772624D-53
    w(53) = 7.22954191064752234D-56
    w(54) = 1.2242353012300822645D-58
    w(55) = 1.48216850490191041178D-61
    w(56) = 1.23251934881451880806D-64
    w(57) = 6.69149900457126952681D-68
    w(58) = 2.220465941850448995507D-71
    w(59) = 4.1209460947388762499791D-75
    w(60) = 3.77439906189648917041762D-79
    w(61) = 1.414115052917619417463147D-83
    w(62) = 1.5918330640413679178610318D-88
    w(63) = 2.98948434886063430774131269D-94
    w(64) = 2.0890635084369527708281542544D-101

  else if ( n == 65 ) then

    x(1) = 0.0220736343882500875264595737093D+00
    x(2) = 0.116318612213376151729622328698D+00
    x(3) = 0.285929513070813951834551909860D+00
    x(4) = 0.531041775784488438911664842933D+00
    x(5) = 0.851801670809046586655299989571D+00
    x(6) = 1.24839628201831707727195703687D+00
    x(7) = 1.72105687981755734816623923359D+00
    x(8) = 2.27006018491690256109784837001D+00
    x(9) = 2.89572940756799471192249006458D+00
    x(10) = 3.59843535756478540694788870796D+00
    x(11) = 4.37859769744118464804519621154D+00
    x(12) = 5.23668636694320050319815958793D+00
    x(13) = 6.17322319658730631921440545744D+00
    x(14) = 7.18878372629467202131889616235D+00
    x(15) = 8.28399924588871831467573911308D+00
    x(16) = 9.45955907610065448076092041454D+00
    x(17) = 10.7162131111897048393914349801D+00
    x(18) = 12.0547746472322707150580323735D+00
    x(19) = 13.4761235235671154056092903386D+00
    x(20) = 14.9812096088541520753799758172D+00
    x(21) = 16.5710566677964491577962746244D+00
    x(22) = 18.2467666498987672118974219097D+00
    x(23) = 20.0095244478287060731090711532D+00
    x(24) = 21.8606031801774165914390699117D+00
    x(25) = 23.8013700618933285635853631176D+00
    x(26) = 25.8332929356397220962712840077D+00
    x(27) = 27.9579475491204073878989540109D+00
    x(28) = 30.1770256774182582795717815112D+00
    x(29) = 32.4923442060863420003362903722D+00
    x(30) = 34.9058553107319841822666865568D+00
    x(31) = 37.4196578929107315489849223258D+00
    x(32) = 40.0360104612770088389674875188D+00
    x(33) = 42.7573456823686148915804841094D+00
    x(34) = 45.5862868687358177811114530444D+00
    x(35) = 48.5256667254367437021219033875D+00
    x(36) = 51.5785487419130140104165726852D+00
    x(37) = 54.7482516984863756567192528045D+00
    x(38) = 58.0383778598867354584946224426D+00
    x(39) = 61.4528455586293273580430507845D+00
    x(40) = 64.9959270371996369214889214525D+00
    x(41) = 68.6722926314626183977917033408D+00
    x(42) = 72.4870626544527259157609330997D+00
    x(43) = 76.4458687019847672439750100021D+00
    x(44) = 80.5549265807842719390299008467D+00
    x(45) = 84.8211237010524840445780520727D+00
    x(46) = 89.2521246438779093867045988317D+00
    x(47) = 93.8564998060605244752886814683D+00
    x(48) = 98.6438836854008198627513004124D+00
    x(49) = 103.625171719648548339736534290D+00
    x(50) = 108.812767977808831403825388295D+00
    x(51) = 114.220900975902615729136419234D+00
    x(52) = 119.866032356536657109576887729D+00
    x(53) = 125.767394661236316522917056970D+00
    x(54) = 131.947712599442095787498603776D+00
    x(55) = 138.434191890473005752666966171D+00
    x(56) = 145.259909991527292147157412636D+00
    x(57) = 152.465831757831086161761248192D+00
    x(58) = 160.103837850755827943698932303D+00
    x(59) = 168.241478626055331550636477652D+00
    x(60) = 176.969855979856624954098184957D+00
    x(61) = 186.417642483351089954653964415D+00
    x(62) = 196.778474440876949259615139650D+00
    x(63) = 208.372107380940491095442543113D+00
    x(64) = 221.812376576320945562507410640D+00
    x(65) = 238.685811594674270102373709659D+00

    w(1) = 0.0554129011565536469555170462551D+00
    w(2) = 0.117417396564162014386769936525D+00
    w(3) = 0.15577793159527363526750345236D+00
    w(4) = 0.166347955884031811873697185854D+00
    w(5) = 0.153013207065446887512281359026D+00
    w(6) = 0.12471061536737329712442529837D+00
    w(7) = 0.0913671486268474804579363734334D+00
    w(8) = 0.0606693673532322224974724445698D+00
    w(9) = 0.0366985078337756899608575633553D+00
    w(10) = 0.0202884358923229233158787730215D+00
    w(11) = 0.0102732022687699783894639990081D+00
    w(12) = 0.00477128781081110626106879095602D+00
    w(13) = 0.00203437021965744474885076853755D+00
    w(14) = 0.00079674216708789273511811886929D+00
    w(15) = 0.000286683812097562728436084278246D+00
    w(16) = 0.0000947743836779584423250711498309D+00
    w(17) = 0.0000287808776386491122522878225936D+00
    w(18) = 8.025921887785674361327140519D-6
    w(19) = 2.0542534210210521080625753887D-6
    w(20) = 4.822974163227069271800228998D-7
    w(21) = 1.037906330975825689484858253D-7
    w(22) = 2.04556458252437904249066042D-8
    w(23) = 3.6885834829334325209249956D-9
    w(24) = 6.07893150020096839531226D-10
    w(25) = 9.14524981574057744811517D-11
    w(26) = 1.25427414380884616209197D-11
    w(27) = 1.56600212784699337922365D-12
    w(28) = 1.7771074191193507379193D-13
    w(29) = 1.829851898040443279934D-14
    w(30) = 1.70645929547049792703D-15
    w(31) = 1.438412537484673440005D-16
    w(32) = 1.09354404172789017674D-17
    w(33) = 7.4805627565827998507D-19
    w(34) = 4.5927490046001844919D-20
    w(35) = 2.5237976161210732972D-21
    w(36) = 1.2376030276622370811D-22
    w(37) = 5.398129612170324038D-24
    w(38) = 2.086925656272980914D-25
    w(39) = 7.12363601216912351D-27
    w(40) = 2.13797949495352507D-28
    w(41) = 5.61585996656691837D-30
    w(42) = 1.2845455347843963D-31
    w(43) = 2.5444410680246895D-33
    w(44) = 4.33793017530204526D-35
    w(45) = 6.3222072485073684D-37
    w(46) = 7.8174607001716233D-39
    w(47) = 8.1320382378781478D-41
    w(48) = 7.0491916339430245D-43
    w(49) = 5.03749289878460876D-45
    w(50) = 2.93162035252999008D-47
    w(51) = 1.370002490040754814D-49
    w(52) = 5.05827802154356236D-52
    w(53) = 1.447822639997111427D-54
    w(54) = 3.141466864703300069D-57
    w(55) = 5.03056294905086669D-60
    w(56) = 5.7547901892214244076D-63
    w(57) = 4.5172924325696051375D-66
    w(58) = 2.31224468889972999556D-69
    w(59) = 7.22311071834335770012D-73
    w(60) = 1.2595483022441380495635D-76
    w(61) = 1.08123622593537779961256D-80
    w(62) = 3.78397004144107162756714D-85
    w(63) = 3.9595563374024690284872814D-90
    w(64) = 6.85929105947869810086088696D-96
    w(65) = 4.3543678721167358517710196959D-103

  else if ( n == 127 ) then

    x(1) = 0.0113396352985186116918931696313D+00
    x(2) = 0.0597497534357266202813482370574D+00
    x(3) = 0.146850986907461676123882236874D+00
    x(4) = 0.272675907358595531313780082789D+00
    x(5) = 0.437246006441926655545770358699D+00
    x(6) = 0.640586882225669295335764164000D+00
    x(7) = 0.882729686390583644814876536500D+00
    x(8) = 1.16371141601665376615605847010D+00
    x(9) = 1.48357501528346138913135848610D+00
    x(10) = 1.84236943516135653806863208099D+00
    x(11) = 2.24014968395790242445133156565D+00
    x(12) = 2.67697687801413036921678699612D+00
    x(13) = 3.15291829570828255657715083088D+00
    x(14) = 3.66804743603047525402263399265D+00
    x(15) = 4.22244408233018884559778766674D+00
    x(16) = 4.81619437158705024756655350873D+00
    x(17) = 5.44939086945594167558621789084D+00
    x(18) = 6.12213265129972541939445847632D+00
    x(19) = 6.83452538941226681122379949733D+00
    x(20) = 7.58668144663674721742059868368D+00
    x(21) = 8.37871997659327252548421206595D+00
    x(22) = 9.21076703074265587779225061024D+00
    x(23) = 10.0829556725286438091664393536D+00
    x(24) = 10.9954260988581254298031473588D+00
    x(25) = 11.9483257691977259976106051279D+00
    x(26) = 12.9418095425855310537233810982D+00
    x(27) = 13.9760398228785065200144056687D+00
    x(28) = 15.0511867125795236315747963654D+00
    x(29) = 16.1674281756128529229773950518D+00
    x(30) = 17.3249502094436734465611637126D+00
    x(31) = 18.5239470269656885608117113093D+00
    x(32) = 19.7646212486115041040716693869D+00
    x(33) = 21.0471841051731836068770440201D+00
    x(34) = 22.3718556518555428176481239181D+00
    x(35) = 23.7388649941224971836523137887D+00
    x(36) = 25.1484505259373682340772783856D+00
    x(37) = 26.6008601810417496072533842798D+00
    x(38) = 28.0963516979646192017539612921D+00
    x(39) = 29.6351928995041789106102271386D+00
    x(40) = 31.2176619874797591442144671526D+00
    x(41) = 32.8440478536104304605229513413D+00
    x(42) = 34.5146504074411491491056359474D+00
    x(43) = 36.2297809223068040196153885089D+00
    x(44) = 37.9897624003999564359687801403D+00
    x(45) = 39.7949299580899617783964371417D+00
    x(46) = 41.6456312327301807051539908975D+00
    x(47) = 43.5422268122868595499508929938D+00
    x(48) = 45.4850906892287911379961513367D+00
    x(49) = 47.4746107402319647194687665991D+00
    x(50) = 49.5111892333790877167288845844D+00
    x(51) = 51.5952433646712444431827712669D+00
    x(52) = 53.7272058258193167582881400691D+00
    x(53) = 55.9075254054475533058306059917D+00
    x(54) = 58.1366676260224391970775260257D+00
    x(55) = 60.4151154190185902957071920538D+00
    x(56) = 62.7433698410518097002071267427D+00
    x(57) = 65.1219508339499963119560254171D+00
    x(58) = 67.5513980319978863144118724431D+00
    x(59) = 70.0322716198845845112298711920D+00
    x(60) = 72.5651532452068490908886694168D+00
    x(61) = 75.1506469897399352993543623251D+00
    x(62) = 77.7893804040858160006474054621D+00
    x(63) = 80.4820056107507292058039629268D+00
    x(64) = 83.2292004811959148867961200190D+00
    x(65) = 86.0316698929535829667982387326D+00
    x(66) = 88.8901470735120510996525185443D+00
    x(67) = 91.8053950383581779949712501705D+00
    x(68) = 94.7782081313315832053870310348D+00
    x(69) = 97.8094136763051164110541101154D+00
    x(70) = 100.899873750172859403719397622D+00
    x(71) = 104.050487088215989347040768450D+00
    x(72) = 107.262191134146004284231164014D+00
    x(73) = 110.535964248515005306027713513D+00
    x(74) = 113.872828090758394853483761877D+00
    x(75) = 117.273850191925177740954778864D+00
    x(76) = 120.740146737188801061739780027D+00
    x(77) = 124.272885579556983542595064469D+00
    x(78) = 127.873289508859426450938417454D+00
    x(79) = 131.542639803143669218093777421D+00
    x(80) = 135.282280093118369701327381064D+00
    x(81) = 139.093620574329700139644220870D+00
    x(82) = 142.978142606436017768082277536D+00
    x(83) = 146.937403744373665494410809691D+00
    x(84) = 150.973043252521871274925114375D+00
    x(85) = 155.086788160346125722296414206D+00
    x(86) = 159.280459926632882354019569899D+00
    x(87) = 163.555981789575711040159671821D+00
    x(88) = 167.915386891943601342455471847D+00
    x(89) = 172.360827284738125368381561917D+00
    x(90) = 176.894583929601921763116749935D+00
    x(91) = 181.519077840368130692275288340D+00
    x(92) = 186.236882528281123738612025304D+00
    x(93) = 191.050737944509291967908366108D+00
    x(94) = 195.963566148798798378390025430D+00
    x(95) = 200.978488976000251536964755261D+00
    x(96) = 206.098848024688711121272830428D+00
    x(97) = 211.328227356716552605723772570D+00
    x(98) = 216.670479376582303234770894658D+00
    x(99) = 222.129754459296872462673049638D+00
    x(100) = 227.710535020722324190891324313D+00
    x(101) = 233.417674882826024533677753226D+00
    x(102) = 239.256444988303086200187496671D+00
    x(103) = 245.232586778715671725312540190D+00
    x(104) = 251.352374887181280300055009918D+00
    x(105) = 257.622691237920614130761918823D+00
    x(106) = 264.051113229082405517543772418D+00
    x(107) = 270.646019457227967492991117186D+00
    x(108) = 277.416717501636510717983882181D+00
    x(109) = 284.373599742208703266744028731D+00
    x(110) = 291.528335213464957195812820216D+00
    x(111) = 298.894108370282486008788956154D+00
    x(112) = 306.485919782626113204181124239D+00
    x(113) = 314.320969864711774874000075076D+00
    x(114) = 322.419155891286796833494403613D+00
    x(115) = 330.803726638024056519338473349D+00
    x(116) = 339.502161278324337477353675960D+00
    x(117) = 348.547375594726973554807617874D+00
    x(118) = 357.979420280298454540490074431D+00
    x(119) = 367.847945200760045788583414229D+00
    x(120) = 378.215906231355328183329791889D+00
    x(121) = 389.165391412510041015794753252D+00
    x(122) = 400.807293314517025899963612864D+00
    x(123) = 413.298536817793844180082600819D+00
    x(124) = 426.875791536636755382885090171D+00
    x(125) = 441.930854853108414124603092718D+00
    x(126) = 459.218046398884299819712673132D+00
    x(127) = 480.693782633883738597042692293D+00

    w(1) = 0.0287732466920001243557700103008D+00
    w(2) = 0.0638174681751346493634809492651D+00
    w(3) = 0.0919196697215705713898641946531D+00
    w(4) = 0.11054167914413766381245463003D+00
    w(5) = 0.118797716333758501883283294227D+00
    w(6) = 0.117378185300526951488044516301D+00
    w(7) = 0.108193059841805514883351455812D+00
    w(8) = 0.0938270752904896280803772614011D+00
    w(9) = 0.0769664509605888439958224859284D+00
    w(10) = 0.0599349039129397143325707300635D+00
    w(11) = 0.0444177420738890013717083162729D+00
    w(12) = 0.0313850809662523209830093722151D+00
    w(13) = 0.021172316041924506411370709025D+00
    w(14) = 0.0136501453642305416521711855646D+00
    w(15) = 0.00841728527105991722793666573854D+00
    w(16) = 0.00496749900598827605159128586202D+00
    w(17) = 0.00280699038950018846319619574464D+00
    w(18) = 0.00151929510039419524604453410578D+00
    w(19) = 0.000787890287517960840862172871405D+00
    w(20) = 0.00039156751064868450584507324649D+00
    w(21) = 0.000186524342688258605500935662601D+00
    w(22) = 0.0000851731604155766219088098281602D+00
    w(23) = 0.0000372856391978530377121453215777D+00
    w(24) = 0.0000156484167917129939474478052968D+00
    w(25) = 6.2964340695224829035692735525D-06
    w(26) = 2.42889297113287245745413799382D-06
    w(27) = 8.9824607890051007201922871545D-07
    w(28) = 3.18441747407603537107429663281D-07
    w(29) = 1.08212729055668392118618075427D-07
    w(30) = 3.52450767506355360159027790853D-08
    w(31) = 1.10012243657193474070638397617D-08
    w(32) = 3.29040796167179321253293430033D-09
    w(33) = 9.4289145237889976419772700773D-10
    w(34) = 2.58825789046683181840501953093D-10
    w(35) = 6.80474371033707626309422590176D-11
    w(36) = 1.71313988051208378353995644756D-11
    w(37) = 4.12917445240528654694439223049D-12
    w(38) = 9.52641897188072732207076648735D-13
    w(39) = 2.10326044324424259329629420475D-13
    w(40) = 4.44271519387293528609404342858D-14
    w(41) = 8.97605003628337033233198464055D-15
    w(42) = 1.73415114077692870746279483468D-15
    w(43) = 3.20280995489883566314943798352D-16
    w(44) = 5.65313889507936820226607420952D-17
    w(45) = 9.53296727990265912345880440259D-18
    w(46) = 1.53534534773101425652885094376D-18
    w(47) = 2.36089621794673656860578421322D-19
    w(48) = 3.46487427944566113321938766532D-20
    w(49) = 4.85152418970864613201269576635D-21
    w(50) = 6.47862286335198134281373737907D-22
    w(51) = 8.2476020965403242936448553126D-23
    w(52) = 1.0005361880214719793491658283D-23
    w(53) = 1.1561395116207304954233181264D-24
    w(54) = 1.271934273116792265561213426D-25
    w(55) = 1.331658471416537296734000416D-26
    w(56) = 1.32612184546789440336461085D-27
    w(57) = 1.25549954476439498072860741D-28
    w(58) = 1.1294412178579462703240913D-29
    w(59) = 9.649102027956211922850061D-31
    w(60) = 7.82418467683020993967331D-32
    w(61) = 6.01815035422196266582499D-33
    w(62) = 4.38824827049617415515105D-34
    w(63) = 3.0314137647517256304036D-35
    w(64) = 1.9826016543944539545225D-36
    w(65) = 1.2267623373665926559014D-37
    w(66) = 7.176393169250888894381D-39
    w(67) = 3.965937883383696358411D-40
    w(68) = 2.068897055386804009958D-41
    w(69) = 1.017958701797951724527D-42
    w(70) = 4.72008277459863746257D-44
    w(71) = 2.06068289855533748257D-45
    w(72) = 8.4627575907305987246D-47
    w(73) = 3.2661123687088798658D-48
    w(74) = 1.1833939207883162381D-49
    w(75) = 4.021120912389501381D-51
    w(76) = 1.2799824394111125389D-52
    w(77) = 3.81238777475488465D-54
    w(78) = 1.061205754270115677D-55
    w(79) = 2.757144694720040359D-57
    w(80) = 6.67725442409284929D-59
    w(81) = 1.505243838386823495D-60
    w(82) = 3.15389868001137585D-62
    w(83) = 6.13266142994831808D-64
    w(84) = 1.10485100303248106D-65
    w(85) = 1.84105635380913481D-67
    w(86) = 2.83239265700528322D-69
    w(87) = 4.01544098437636555D-71
    w(88) = 5.23515302156837088D-73
    w(89) = 6.2634476665005101D-75
    w(90) = 6.861221053566653D-77
    w(91) = 6.8651298840956019D-79
    w(92) = 6.25813884337280849D-81
    w(93) = 5.1833271237514904D-83
    w(94) = 3.88936215719184435D-85
    w(95) = 2.63577113794769328D-87
    w(96) = 1.60788512939179797D-89
    w(97) = 8.79780420709689396D-92
    w(98) = 4.30134050774951099D-94
    w(99) = 1.871343588134283853D-96
    w(100) = 7.212574470806047168D-99
    w(101) = 2.450874606217787438D-101
    w(102) = 7.304209461947087578D-104
    w(103) = 1.8983290818383463538D-106
    w(104) = 4.2757400244246684123D-109
    w(105) = 8.2894681420515755691D-112
    w(106) = 1.37294322193244000131D-114
    w(107) = 1.926546412640497322204D-117
    w(108) = 2.269334450330135482614D-120
    w(109) = 2.2209290603717355061909D-123
    w(110) = 1.7851087685544512662857D-126
    w(111) = 1.16309319903871644674312D-129
    w(112) = 6.05244435846523922909528D-133
    w(113) = 2.472956911506352864762838D-136
    w(114) = 7.778906500648941036499721D-140
    w(115) = 1.8409738662712607039570678D-143
    w(116) = 3.1900921131079114970179072D-147
    w(117) = 3.917948713917419973761766608D-151
    w(118) = 3.2782158394188697053774429821D-155
    w(119) = 1.77935907131388880628196401287D-159
    w(120) = 5.88823534089326231574678353812D-164
    w(121) = 1.09572365090711698777472032739D-168
    w(122) = 1.02816211148670008982850769758D-173
    w(123) = 4.1704725557697758145816510854D-179
    w(124) = 5.8002877720316101774638319602D-185
    w(125) = 1.88735077458255171061716191011D-191
    w(126) = 6.91066018267309116827867059509D-199
    w(127) = 4.35068132011058556283833133344D-208

  else if ( n == 128 ) then

    x(1) = 0.0112513882636759629608518403162D+00
    x(2) = 0.0592847412690264542879220089614D+00
    x(3) = 0.145707966594312465141854059102D+00
    x(4) = 0.270553178758665066190760897100D+00
    x(5) = 0.433841407553836803056096580754D+00
    x(6) = 0.635597665781621938340867677969D+00
    x(7) = 0.875852384546520779155346013261D+00
    x(8) = 1.15464170197439795008153355708D+00
    x(9) = 1.47200756316673547554446633038D+00
    x(10) = 1.82799777831235028535984528718D+00
    x(11) = 2.22266607156190244817896452914D+00
    x(12) = 2.65607212988348119522885329309D+00
    x(13) = 3.12828165502791695498310369738D+00
    x(14) = 3.63936641985240321221074522169D+00
    x(15) = 4.18940432959404797478493079865D+00
    x(16) = 4.77847948843487609165724239213D+00
    x(17) = 5.40668227160049918527893820105D+00
    x(18) = 6.07410940319653684309155844506D+00
    x(19) = 6.78086403997562541104804929943D+00
    x(20) = 7.52705586122937588585512279842D+00
    x(21) = 8.31280116500777060337884381191D+00
    x(22) = 9.13822297088039239600262641969D+00
    x(23) = 10.0034511294682220892682761435D+00
    x(24) = 10.9086224389908825478488613010D+00
    x(25) = 11.8538807690918568038332644538D+00
    x(26) = 12.8393771922232496874551935673D+00
    x(27) = 13.8652701228920803029536799971D+00
    x(28) = 14.9317254650919274737473553133D+00
    x(29) = 16.0389167682670793509213783428D+00
    x(30) = 17.1870253921812651027585044591D+00
    x(31) = 18.3762406810896949333827523370D+00
    x(32) = 19.6067601476416467279054690989D+00
    x(33) = 20.8787896669713729158932014403D+00
    x(34) = 22.1925436814678369066763923182D+00
    x(35) = 23.5482454167489205609249730097D+00
    x(36) = 24.9461271094034886279510396640D+00
    x(37) = 26.3864302471052908976269305132D+00
    x(38) = 27.8694058217463902295696818564D+00
    x(39) = 29.3953145962849137215656288381D+00
    x(40) = 30.9644273860527540023220861317D+00
    x(41) = 32.5770253553237501781419486456D+00
    x(42) = 34.2334003300022426604794108753D+00
    x(43) = 35.9338551273561538107722924963D+00
    x(44) = 37.6787039037883744300655582016D+00
    x(45) = 39.4682725217157641271489033004D+00
    x(46) = 41.3028989367070896380080417637D+00
    x(47) = 43.1829336061203832438783635225D+00
    x(48) = 45.1087399205772317441506148507D+00
    x(49) = 47.0806946597172168560725351128D+00
    x(50) = 49.0991884737910268021535852860D+00
    x(51) = 51.1646263927766594446404335916D+00
    x(52) = 53.2774283648407739161085367944D+00
    x(53) = 55.4380298261178918683089638291D+00
    x(54) = 57.6468823039452288144249811220D+00
    x(55) = 59.9044540558720556965292635062D+00
    x(56) = 62.2112307469614582456791552962D+00
    x(57) = 64.5677161681212154290410515467D+00
    x(58) = 66.9744329984415610548027156195D+00
    x(59) = 69.4319236147834299557621097742D+00
    x(60) = 71.9407509521543751573018481062D+00
    x(61) = 74.5014994187340277930279831855D+00
    x(62) = 77.1147758697705705283198924354D+00
    x(63) = 79.7812106449685528544582124991D+00
    x(64) = 82.5014586744314529140391768845D+00
    x(65) = 85.2762006587153587377964042582D+00
    x(66) = 88.1061443290995036940317393258D+00
    x(67) = 90.9920257947926131560303030245D+00
    x(68) = 93.9346109844796944955244642925D+00
    x(69) = 96.9346971903819404516199551240D+00
    x(70) = 99.9931147238642715216213000267D+00
    x(71) = 103.110728692593987392319749158D+00
    x(72) = 106.288440910345442668426129659D+00
    x(73) = 109.527191951777550806618056918D+00
    x(74) = 112.827963365904193877333487264D+00
    x(75) = 116.191780063556780940235871708D+00
    x(76) = 119.619712895932010462348887420D+00
    x(77) = 123.112881443360190060911814509D+00
    x(78) = 126.672457035760183662338694957D+00
    x(79) = 130.299666028913462587217492864D+00
    x(80) = 133.995793363747964343582120836D+00
    x(81) = 137.762186439339380964302666578D+00
    x(82) = 141.600259334393040305789642722D+00
    x(83) = 145.511497416659592393640597008D+00
    x(84) = 149.497462385177707088173175451D+00
    x(85) = 153.559797796566440117982748261D+00
    x(86) = 157.700235133978105059095336546D+00
    x(87) = 161.920600485975634163753629031D+00
    x(88) = 166.222821912768875092875739160D+00
    x(89) = 170.608937589242234646550663310D+00
    x(90) = 175.081104828414604880617405502D+00
    x(91) = 179.641610105866994602634964639D+00
    x(92) = 184.292880225846805341834505020D+00
    x(93) = 189.037494793954109001292998345D+00
    x(94) = 193.878200190472967540802875940D+00
    x(95) = 198.817925273720602804745944994D+00
    x(96) = 203.859799085769844571664824897D+00
    x(97) = 209.007170885510867853387511181D+00
    x(98) = 214.263632898788021280527758492D+00
    x(99) = 219.633046255578174038387401024D+00
    x(100) = 225.119570684209027756659796566D+00
    x(101) = 230.727698658203619681658868680D+00
    x(102) = 236.462294850177665966018904158D+00
    x(103) = 242.328641949702698267864519866D+00
    x(104) = 248.332494162357178892016601780D+00
    x(105) = 254.480140044869131893543803358D+00
    x(106) = 260.778476773579736538560064538D+00
    x(107) = 267.235098528953836763992472029D+00
    x(108) = 273.858402462693609793414602648D+00
    x(109) = 280.657716776323492397504100977D+00
    x(110) = 287.643456899219330638473677900D+00
    x(111) = 294.827317787647739179806672104D+00
    x(112) = 302.222513246449465380535981711D+00
    x(113) = 309.844077326612663447772363643D+00
    x(114) = 317.709248954906289495678052340D+00
    x(115) = 325.837970121194949650401277931D+00
    x(116) = 334.253542067654135375184450174D+00
    x(117) = 342.983506273825316408508913329D+00
    x(118) = 352.060853546526185083043426984D+00
    x(119) = 361.525726392325047599066851839D+00
    x(120) = 371.427889214327523078517984867D+00
    x(121) = 381.830444119061080196207616882D+00
    x(122) = 392.815671240808098809377819898D+00
    x(123) = 404.494724750515074389666071660D+00
    x(124) = 417.024902977989015820197277594D+00
    x(125) = 430.643444166597381558323551668D+00
    x(126) = 445.743096973927989652171720726D+00
    x(127) = 463.080034109446258208013793406D+00
    x(128) = 484.615543986443976044063131110D+00

    w(1) = 0.0285518444532397286290731773612D+00
    w(2) = 0.0633502117845051187797978127259D+00
    w(3) = 0.0913083813661343144231616325903D+00
    w(4) = 0.109913900410911746101013915833D+00
    w(5) = 0.118274171034173698789809688874D+00
    w(6) = 0.117045739000406721566458439207D+00
    w(7) = 0.108089987545568415436473783125D+00
    w(8) = 0.0939428886389285017878088436356D+00
    w(9) = 0.0772536687978980532077800252359D+00
    w(10) = 0.0603270562656615705389303086003D+00
    w(11) = 0.0448473482471952140682424657998D+00
    w(12) = 0.0317969479368768461739632484821D+00
    w(13) = 0.0215301494537944439261107285438D+00
    w(14) = 0.0139369517338463483277576885975D+00
    w(15) = 0.00863158538020224714884473096489D+00
    w(16) = 0.00511777701366922852873936722845D+00
    w(17) = 0.00290634743648595585817980077219D+00
    w(18) = 0.00158143294331667939723416013489D+00
    w(19) = 0.000824738985098812150435438593253D+00
    w(20) = 0.000412326088539694730970290830804D+00
    w(21) = 0.000197649426442591498620529889783D+00
    w(22) = 0.0000908515788782451508022826306153D+00
    w(23) = 0.0000400484927835805298977887660442D+00
    w(24) = 0.0000169307623980817855755102888475D+00
    w(25) = 6.86452529111068208938636278412D-06
    w(26) = 2.66921659814210266015872228584D-06
    w(27) = 9.95364010286384477177483332196D-07
    w(28) = 3.55943575300306543988020166563D-07
    w(29) = 1.22053255194881194831205615734D-07
    w(30) = 4.01279192093563506890167766024D-08
    w(31) = 1.26481141474759786445650110908D-08
    w(32) = 3.82148972942657229023411372003D-09
    w(33) = 1.10664105922734169994480044024D-09
    w(34) = 3.07100923709742319582290034639D-10
    w(35) = 8.16549938415448956026437885004D-11
    w(36) = 2.07985363278137784234189612116D-11
    w(37) = 5.0739537708398704043296986402D-12
    w(38) = 1.1853143771796112305093733131D-12
    w(39) = 2.65092752372887358600565488195D-13
    w(40) = 5.67463221575765876681065606161D-14
    w(41) = 1.16237381470751589221529434901D-14
    w(42) = 2.27776629270238637919733104451D-15
    w(43) = 4.26883197029764927739172104126D-16
    w(44) = 7.64928879936327510525948457803D-17
    w(45) = 1.31013139198382464188082886821D-17
    w(46) = 2.1441452341246636343706788692D-18
    w(47) = 3.35194428720884780801470729044D-19
    w(48) = 5.00373308645947376823179365121D-20
    w(49) = 7.13003064195856212049702464626D-21
    w(50) = 9.6945407403972664035320905829D-22
    w(51) = 1.25728475563978459844059927432D-22
    w(52) = 1.5546610955630634482202731199D-23
    w(53) = 1.832109793253421778719084254D-24
    w(54) = 2.056797978136734920722781372D-25
    w(55) = 2.19866605262329119257657449D-26
    w(56) = 2.23691600732428936729406222D-27
    w(57) = 2.1649606446339054400256309D-28
    w(58) = 1.9922276806187937873877251D-29
    w(59) = 1.742153886325439585907653D-30
    w(60) = 1.446949786106284637699605D-31
    w(61) = 1.1407517061230822834189D-32
    w(62) = 8.5318050978102090722116D-34
    w(63) = 6.04970117793885843505D-35
    w(64) = 4.0643432432648003017795D-36
    w(65) = 2.585349374987909630703D-37
    w(66) = 1.556028762522623447585D-38
    w(67) = 8.85462584966333001103D-40
    w(68) = 4.76045751736458068032D-41
    w(69) = 2.416078510661232205D-42
    w(70) = 1.15664705033873749321D-43
    w(71) = 5.2185106194923759952D-45
    w(72) = 2.2169743353361803305D-46
    w(73) = 8.86010275661369606D-48
    w(74) = 3.327811159201095553D-49
    w(75) = 1.173490043078302544D-50
    w(76) = 3.880967726420921431D-52
    w(77) = 1.202426327933061418D-53
    w(78) = 3.48602304410554638D-55
    w(79) = 9.44554522159556681D-57
    w(80) = 2.38888427455968395D-58
    w(81) = 5.63188475075463052D-60
    w(82) = 1.23592861191216019D-61
    w(83) = 2.52100420237726743D-63
    w(84) = 4.7722246219998052D-65
    w(85) = 8.3700198919995783D-67
    w(86) = 1.35782434112020985D-68
    w(87) = 2.03368872715315416D-70
    w(88) = 2.8068384806953538D-72
    w(89) = 3.562567607062096D-74
    w(90) = 4.1494527492937706D-76
    w(91) = 4.4250079657663219D-78
    w(92) = 4.3100842612898497D-80
    w(93) = 3.8246610167617398D-82
    w(94) = 3.08354784259879275D-84
    w(95) = 2.25213982217062084D-86
    w(96) = 1.48551474064504312D-88
    w(97) = 8.8196354763726564D-91
    w(98) = 4.69641782212598507D-93
    w(99) = 2.23439382545477274D-95
    w(100) = 9.45878703822074032D-98
    w(101) = 3.546960831240672614D-100
    w(102) = 1.17253213003488723D-102
    w(103) = 3.399090555639915548D-105
    w(104) = 8.591907200623898045D-108
    w(105) = 1.8818913973535359647D-110
    w(106) = 3.5473586323062565237D-113
    w(107) = 5.7114822282836004745D-116
    w(108) = 7.78947378804446095611D-119
    w(109) = 8.91589869949126935148D-122
    w(110) = 8.476856358868403207418D-125
    w(111) = 6.617326935494900345408D-128
    w(112) = 4.1862163574157095190077D-131
    w(113) = 2.11438516898114207120093D-134
    w(114) = 8.38216350136786953641675D-138
    w(115) = 2.557202302197677884687798D-141
    w(116) = 5.8667686421912043720461236D-145
    w(117) = 9.8498610300648438019885689D-149
    w(118) = 1.171383943342068942456857274D-152
    w(119) = 9.483963265567383663821702301D-157
    w(120) = 4.9770963811238028116653976343D-161
    w(121) = 1.59089852775099765481980638695D-165
    w(122) = 2.85630382911900292320607568044D-170
    w(123) = 2.58225071969148999265031459122D-175
    w(124) = 1.00735025005079740952983187255D-180
    w(125) = 1.34425250044381631821772983363D-186
    w(126) = 4.18296221403683473389726627221D-193
    w(127) = 1.45716530772618631594481663188D-200
    w(128) = 8.64059169046870867692891422354D-210

  else if ( n == 129 ) then

    x(1) = 0.0111645041367687260935881187114D+00
    x(2) = 0.0588269115255121725669144777376D+00
    x(3) = 0.144582603939087375345544455104D+00
    x(4) = 0.268463250498790809142537571727D+00
    x(5) = 0.430489433028069665583513882755D+00
    x(6) = 0.630685596971157529700818698614D+00
    x(7) = 0.869081474989540465988995980646D+00
    x(8) = 1.14571237269034129786358037349D+00
    x(9) = 1.46061926689785560022252086358D+00
    x(10) = 1.81384886225287260620048305182D+00
    x(11) = 2.20545363849013952710373368048D+00
    x(12) = 2.63549189753739459262727316570D+00
    x(13) = 3.10402781353627480023526416641D+00
    x(14) = 3.61113148701933289479734007535D+00
    x(15) = 4.15687900382495881133416031205D+00
    x(16) = 4.74135249908325871733484826319D+00
    x(17) = 5.36464022650680264548807369539D+00
    x(18) = 6.02683663318167548105631862177D+00
    x(19) = 6.72804244004243132025609101021D+00
    x(20) = 7.46836472821534963467632383543D+00
    x(21) = 8.24791703142169723816558449856D+00
    x(22) = 9.06681943464370270026626900050D+00
    x(23) = 9.92519867926931734070041188408D+00
    x(24) = 10.8231882749469306495297612192D+00
    x(25) = 11.7609286183977310387181197615D+00
    x(26) = 12.7385671194512351722694605084D+00
    x(27) = 13.7562583345886271805335101149D+00
    x(28) = 14.8141641082989854857712290559D+00
    x(29) = 15.9124537225752979381236294324D+00
    x(30) = 17.0513040549004685335351914932D+00
    x(31) = 18.2308997450984136591429080617D+00
    x(32) = 19.4514333714519620150207362048D+00
    x(33) = 20.7131056365177555775299262985D+00
    x(34) = 22.0161255630988608706489924430D+00
    x(35) = 23.3607107008685190470514486749D+00
    x(36) = 24.7470873441735867407744490421D+00
    x(37) = 26.1754907615839641134296855243D+00
    x(38) = 27.6461654377949106206644801830D+00
    x(39) = 29.1593653285328756045576144321D+00
    x(40) = 30.7153541291626095441732915451D+00
    x(41) = 32.3144055577441922161871665693D+00
    x(42) = 33.9568036533435689847296719094D+00
    x(43) = 35.6428430904596160112634717165D+00
    x(44) = 37.3728295104950910213327545948D+00
    x(45) = 39.1470798712685466663878582421D+00
    x(46) = 40.9659228156399190364649448364D+00
    x(47) = 42.8296990604046437906422357564D+00
    x(48) = 44.7387618067004519884778950119D+00
    x(49) = 46.6934771732681867037686990052D+00
    x(50) = 48.6942246540138734219622380649D+00
    x(51) = 50.7413976014347803131042845818D+00
    x(52) = 52.8354037375983340937979164025D+00
    x(53) = 54.9766656945006500240481182310D+00
    x(54) = 57.1656215857823649284158070179D+00
    x(55) = 59.4027256119448606421881531943D+00
    x(56) = 61.6884487013914405461822377003D+00
    x(57) = 64.0232791898173852597210780437D+00
    x(58) = 66.4077235406921080587914699631D+00
    x(59) = 68.8423071098181639647332557636D+00
    x(60) = 71.3275749572182499453797757024D+00
    x(61) = 73.8640927098955268421782575110D+00
    x(62) = 76.4524474793379566181613942983D+00
    x(63) = 79.0932488379977030145472340597D+00
    x(64) = 81.7871298593763443093790704140D+00
    x(65) = 84.5347482267906647046323684124D+00
    x(66) = 87.3367874163878117910865422310D+00
    x(67) = 90.1939579605291478450570652459D+00
    x(68) = 93.1069987982766656767050611186D+00
    x(69) = 96.0766787204029972427806506124D+00
    x(70) = 99.1037979171157474398207782757D+00
    x(71) = 102.189189637550616978355114969D+00
    x(72) = 105.333721971058838012388514189D+00
    x(73) = 108.538299761408281119757506569D+00
    x(74) = 111.803866666252185387269569516D+00
    x(75) = 115.131407375615803792171876281D+00
    x(76) = 118.521950004733905726449958829D+00
    x(77) = 121.976568678369858697173472594D+00
    x(78) = 125.496386325793836628100280130D+00
    x(79) = 129.082577707933597477650969878D+00
    x(80) = 132.736372700883616552797038522D+00
    x(81) = 136.459059863023413416361147154D+00
    x(82) = 140.251990316520692590651584246D+00
    x(83) = 144.116581978059547282038191264D+00
    x(84) = 148.054324178334554730971189024D+00
    x(85) = 152.066782715303825347545842677D+00
    x(86) = 156.155605392537787829354492826D+00
    x(87) = 160.322528101405362530717313062D+00
    x(88) = 164.569381514511906899139962575D+00
    x(89) = 168.898098467996847713856358122D+00
    x(90) = 173.310722122324145053009369479D+00
    x(91) = 177.809415005439611927392788370D+00
    x(92) = 182.396469059102149766157602559D+00
    x(93) = 187.074316829415599837320402996D+00
    x(94) = 191.845543966839763114071401295D+00
    x(95) = 196.712903230183761859963922576D+00
    x(96) = 201.679330224475912387142876872D+00
    x(97) = 206.747961145685983577272619640D+00
    x(98) = 211.922152858007833039218477193D+00
    x(99) = 217.205505694330143211608701365D+00
    x(100) = 222.601889450939732907762241579D+00
    x(101) = 228.115473147766504188912042615D+00
    x(102) = 233.750759251359480867215068547D+00
    x(103) = 239.512623216994726852324857048D+00
    x(104) = 245.406359409276117119061170837D+00
    x(105) = 251.437734721509742305967783880D+00
    x(106) = 257.613051552607945102309836371D+00
    x(107) = 263.939222243647173894814944292D+00
    x(108) = 270.423857663083214051346532320D+00
    x(109) = 277.075373415313344577287378499D+00
    x(110) = 283.903118212107869887400929941D+00
    x(111) = 290.917530409009503510470042900D+00
    x(112) = 298.130330747241946479391511151D+00
    x(113) = 305.554762228622700877217556637D+00
    x(114) = 313.205892212538716296350101328D+00
    x(115) = 321.100997941634721519100026717D+00
    x(116) = 329.260065894410473958350680155D+00
    x(117) = 337.706449515634131989920236326D+00
    x(118) = 346.467752279350659621376501841D+00
    x(119) = 355.577039643063413893183224979D+00
    x(120) = 365.074545471124791778391196263D+00
    x(121) = 375.010148136708978052975802762D+00
    x(122) = 385.447095254054417308720464640D+00
    x(123) = 396.467858100744210106334636127D+00
    x(124) = 408.183851152492844798297769341D+00
    x(125) = 420.752744334742187498526476928D+00
    x(126) = 434.412341688764625555428748148D+00
    x(127) = 449.556338392256949417199002480D+00
    x(128) = 466.942750921706688536121321308D+00
    x(129) = 488.537715007400745716181291102D+00

    w(1) = 0.0283338232816188129433412493366D+00
    w(2) = 0.0628897352309939992519628028429D+00
    w(3) = 0.0907050560197830441591715791845D+00
    w(4) = 0.109292734964339745013347523543D+00
    w(5) = 0.117753891824430328742552706746D+00
    w(6) = 0.116712333575132760088854393741D+00
    w(7) = 0.1079821092277907522768638822D+00
    w(8) = 0.0940513886437790878162542877426D+00
    w(9) = 0.0775328171368385256641246588694D+00
    w(10) = 0.0607119801995722871258201910352D+00
    w(11) = 0.0452716214541695196710137988047D+00
    w(12) = 0.0322057586869443933590250840601D+00
    w(13) = 0.0218870093879284288723521418152D+00
    w(14) = 0.0142243242185532561642375502974D+00
    w(15) = 0.00884734285745239479408590424342D+00
    w(16) = 0.00526983370954167607842815218011D+00
    w(17) = 0.00300740619275414763773247784756D+00
    w(18) = 0.00164498171784021535901621253553D+00
    w(19) = 0.000862641473273809069700952476134D+00
    w(20) = 0.000433807488545501081264834235514D+00
    w(21) = 0.000209234988721404556453070968853D+00
    w(22) = 0.0000968044053231071525887634259114D+00
    w(23) = 0.0000429650601010182583779356860953D+00
    w(24) = 0.0000182943298240488545326843922155D+00
    w(25) = 7.47320473307839845584026474317D-06
    w(26) = 2.92876004890558731746712968433D-06
    w(27) = 1.10111937532188602299646730309D-06
    w(28) = 3.97133727854894494886436944708D-07
    w(29) = 1.37391737873739072964678053016D-07
    w(30) = 4.55898285044463980401770363171D-08
    w(31) = 1.45082031554226827387170770004D-08
    w(32) = 4.42736861865778798052557346184D-09
    w(33) = 1.29540549841465072618582643105D-09
    w(34) = 3.63353401896969889688016611161D-10
    w(35) = 9.76889957112077658988662065766D-11
    w(36) = 2.51697359198850123687093430919D-11
    w(37) = 6.2136427115425329244941688694D-12
    w(38) = 1.46947065273427272155102255087D-12
    w(39) = 3.3283536396381226771786168693D-13
    w(40) = 7.21860543546415622515782097245D-14
    w(41) = 1.49874700296546634758941598894D-14
    w(42) = 2.97813865190408297766537928957D-15
    w(43) = 5.66223500996744709699260363288D-16
    w(44) = 1.02976110977345161229212606736D-16
    w(45) = 1.79087076765055918801501255712D-17
    w(46) = 2.97741214327584722879794953728D-18
    w(47) = 4.73066849378813640521244315218D-19
    w(48) = 7.18076704552391091114386577815D-20
    w(49) = 1.0409591754013912471892470954D-20
    w(50) = 1.44063705945958837668569771815D-21
    w(51) = 1.9027009013059586477368991424D-22
    w(52) = 2.3972421860336028068385342016D-23
    w(53) = 2.880065029076382866335001882D-24
    w(54) = 3.2980570110683255202892323D-25
    w(55) = 3.59822818119059018987046195D-26
    w(56) = 3.7384843519427824153681456D-27
    w(57) = 3.6971969670644497346136084D-28
    w(58) = 3.478607942989822329014257D-29
    w(59) = 3.112229078360896126467536D-30
    w(60) = 2.64630166366922810478446D-31
    w(61) = 2.13731385180863223984415D-32
    w(62) = 1.63873356712820982018691D-33
    w(63) = 1.1920639048111247727415D-34
    w(64) = 8.221888191494076473793D-36
    w(65) = 5.373327742595686629791D-37
    w(66) = 3.325235584661609413228D-38
    w(67) = 1.94717317556033610096D-39
    w(68) = 1.07813497736466418105D-40
    w(69) = 5.6402504582069233692D-42
    w(70) = 2.785716667292756732D-43
    w(71) = 1.2978694111929463222D-44
    w(72) = 5.699117216622829387D-46
    w(73) = 2.356556045713220169D-47
    w(74) = 9.167179452095711245D-49
    w(75) = 3.351643630271094859D-50
    w(76) = 1.15053967361148792D-51
    w(77) = 3.70428664291287775D-53
    w(78) = 1.117334474142203311D-54
    w(79) = 3.15377989811063792D-56
    w(80) = 8.319920981942047D-58
    w(81) = 2.04876111892933112D-59
    w(82) = 4.7028955186049464D-61
    w(83) = 1.00491633674668433D-62
    w(84) = 1.9959187047623038D-64
    w(85) = 3.6789923736675531D-66
    w(86) = 6.2831482675040959D-68
    w(87) = 9.925201342288209D-70
    w(88) = 1.4475221077412768D-71
    w(89) = 1.945364935931307D-73
    w(90) = 2.4042822695448614D-75
    w(91) = 2.7267496829701407D-77
    w(92) = 2.8313374255297656D-79
    w(93) = 2.6851895059223692D-81
    w(94) = 2.3199549783717045D-83
    w(95) = 1.821032672647817D-85
    w(96) = 1.29486019972133753D-87
    w(97) = 8.3146100960594316D-90
    w(98) = 4.8053665090563748D-92
    w(99) = 2.49071240066108676D-94
    w(100) = 1.15335704284873844D-96
    w(101) = 4.75169815023478164D-99
    w(102) = 1.733951399870136754D-101
    w(103) = 5.57731896834145892D-104
    w(104) = 1.573010564351007982D-106
    w(105) = 3.867845242632879313D-109
    w(106) = 8.239883435606238718D-112
    w(107) = 1.5104570697877326124D-114
    w(108) = 2.3645657754433596259D-117
    w(109) = 3.1349053289923477642D-120
    w(110) = 3.48739145376585928069D-123
    w(111) = 3.22170074744057989255D-126
    w(112) = 2.443048415722317309221D-129
    w(113) = 1.5008657805760609578501D-132
    w(114) = 7.3592251345721592465131D-136
    w(115) = 2.83121162238276011127992D-139
    w(116) = 8.3785828758598937096069D-143
    w(117) = 1.8637689328976254234922931D-146
    w(118) = 3.0323700940390393081087066D-150
    w(119) = 3.49260330326226204565809172D-154
    w(120) = 2.736761201290944128360070077D-158
    w(121) = 1.3888959774881077581342370711D-162
    w(122) = 4.28912860126508716947322409477D-167
    w(123) = 7.43133882324715291928018394696D-172
    w(124) = 6.47421443374096511679045401121D-177
    w(125) = 2.42953692988216878005255824922D-182
    w(126) = 3.11143287762562176520260181694D-188
    w(127) = 9.26127289624597363219192415542D-195
    w(128) = 3.07023341560782650495387872798D-202
    w(129) = 1.71530871887294016615286222244D-211

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) &
      '  Legal values are 1 to 20, 31/32/33, 63/64/65 or 127/128/129.'
    stop

  end if

  return
end
subroutine laguerre_ss_compute ( n, x, w )

!*****************************************************************************80
!
!! LAGUERRE_SS_COMPUTE computes a Gauss-Laguerre quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) exp ( - x ) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * exp ( x(i) ) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be at least 1.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) ratio
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval
!
!  Set the recursion coefficients.
!
  do i = 1, n
    b(i) = ( real ( 2 * i - 1, kind = 8 ) )
  end do

  do i = 1, n
    c(i) = real ( i - 1, kind = 8 ) * ( real ( i - 1, kind = 8 ) )
  end do

  cc = product ( c(2:n) )

  do i = 1, n
!
!  Compute an estimate for the root.
!
    if ( i == 1 ) then

      xval = 3.0D+00 / ( 1.0D+00 + 2.4D+00 * real ( n, kind = 8 ) )

    else if ( i == 2 ) then

      xval = xval + 15.0D+00 / ( 1.0D+00 + 2.5D+00 * real ( n, kind = 8 ) )

    else

      r1 = ( 1.0D+00 + 2.55D+00 * real ( i - 2, kind = 8 ) ) &
        / ( 1.9D+00 * real ( i - 2, kind = 8 ) )

      xval = xval + r1 * ( xval - x(i-2) )

    end if
!
!  Use iteration to find the root.
!
    call laguerre_ss_root ( xval, n, dp2, p1, b, c )
!
!  Set the abscissa and weight.
!
    x(i) = xval
    w(i) = ( cc / dp2 ) / p1

  end do

  return
end
subroutine laguerre_ss_recur ( p2, dp2, p1, x, n, b, c )

!*****************************************************************************80
!
!! LAGUERRE_SS_RECUR finds the value and derivative of a Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of L(N)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of L'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(N-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) B(N), C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x - 1.0D+00
  dp2 = 1.0D+00

  do i = 2, n

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine laguerre_ss_root ( x, n, dp2, p1, b, c )

!*****************************************************************************80
!
!! LAGUERRE_SS_ROOT improves an approximate root of a Laguerre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) DP2, the value of L'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(N-1)(X).
!
!    Input, real ( kind = 8 ) B(N), C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(N)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = epsilon ( eps )

  do step = 1, step_max

    call laguerre_ss_recur ( p2, dp2, p1, x, n, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine laguerre_sum ( func, a, n, x, w, result )

!*****************************************************************************80
!
!! LAGUERRE_SUM carries out Laguerre quadrature over [ A, +oo ).
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x <= +oo ) exp ( -x ) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The integral:
!
!      integral ( a <= x <= +oo ) exp ( -x ) * f(x) dx
!
!    The quadrature rule:
!
!      exp ( - a ) * sum ( 1 <= i <= n ) w(i) * f ( x(i) + a )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, real ( kind = 8 ) A, the beginning of the integration interval.
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Input, real ( kind = 8 ) W(N), the weights.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_SUM - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop
  end if

  result = 0.0D+00
  do i = 1, n
    result = result + w(i) * func ( x(i) + a )
  end do
  result = exp ( - a ) * result

  return
end
subroutine legendre_cos2_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_COS2_SET sets a Gauss-Legendre rule for COS(X) * F(X) on [0,PI/2].
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x <= pi/2 ) cos(x) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The integral:
!
!      integral ( 0 <= x <= pi/2 ) sin(x) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( pi/2 - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gwynne Evans,
!    Practical Numerical Integration,
!    Wiley, 1993,
!    ISBN: 047193898X,
!    LC: QA299.3E93.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 2, 4, 8 or 16.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 2 ) then

    x(1) = 0.26587388056307823382D+00
    x(2) = 1.0351526093171315182D+00

    w(1) = 0.60362553280827113087D+00
    w(2) = 0.39637446719172886913D+00

  else if ( n == 4 ) then

    x(1) = 0.095669389196858636773D+00
    x(2) = 0.45240902327067096554D+00
    x(3) = 0.93185057672024082424D+00
    x(4) = 1.3564439599666466230D+00

    w( 1) = 0.23783071419515504517D+00
    w( 2) = 0.40265695523581253512D+00
    w( 3) = 0.28681737948564715225D+00
    w( 4) = 0.072694951083385267446D+00

  else if ( n == 8 ) then

    x(1) = 0.029023729768913933432D+00
    x(2) = 0.14828524404581819442D+00
    x(3) = 0.34531111151664787488D+00
    x(4) = 0.59447696797658360178D+00
    x(5) = 0.86538380686123504827D+00
    x(6) = 1.1263076093187456632D+00
    x(7) = 1.3470150460281258016D+00
    x(8) = 1.5015603622059195568D+00

    w( 1) = 0.073908998095117384985D+00
    w( 2) = 0.16002993702338006099D+00
    w( 3) = 0.21444434341803549108D+00
    w( 4) = 0.21979581268851903339D+00
    w( 5) = 0.17581164478209568886D+00
    w( 6) = 0.10560448025308322171D+00
    w( 7) = 0.042485497299217201089D+00
    w( 8) = 0.0079192864405519178899D+00

  else if ( n == 16 ) then

    x( 1) = 0.0080145034906295973494D+00
    x( 2) = 0.041893031354246254797D+00
    x( 3) = 0.10149954486757579459D+00
    x( 4) = 0.18463185923836617507D+00
    x( 5) = 0.28826388487760574589D+00
    x( 6) = 0.40870579076464794191D+00
    x( 7) = 0.54176054986913847463D+00
    x( 8) = 0.68287636658719416893D+00
    x( 9) = 0.82729287620416833520D+00
    x(10) = 0.97018212594829367065D+00
    x(11) = 1.1067865150286247873D+00
    x(12) = 1.2325555697227748824D+00
    x(13) = 1.3432821921580721861D+00
    x(14) = 1.4352370549295032923D+00
    x(15) = 1.5052970876794669248D+00
    x(16) = 1.5510586944086135769D+00

    w( 1) = 0.020528714977215248902D+00
    w( 2) = 0.046990919853597958123D+00
    w( 3) = 0.071441021312218541698D+00
    w( 4) = 0.092350338329243052271D+00
    w( 5) = 0.10804928026816236935D+00
    w( 6) = 0.11698241243306261791D+00
    w( 7) = 0.11812395361762037649D+00
    w( 8) = 0.11137584940420091049D+00
    w( 9) = 0.097778236145946543110D+00
    w(10) = 0.079418758985944482077D+00
    w(11) = 0.059039620053768691402D+00
    w(12) = 0.039458876783728165671D+00
    w(13) = 0.022987785677206847531D+00
    w(14) = 0.011010405600421536861D+00
    w(15) = 0.0038123928030499915653D+00
    w(16) = 0.00065143375461266656171D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COS2_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  end if

  return
end
subroutine legendre_cos_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_COS_SET: Gauss-Legendre rule for COS(X) * F(X) on [-PI/2,PI/2].
!
!  Discussion:
!
!    The integral:
!
!      integral ( -pi/2 <= x <= pi/2 ) cos(x) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The integral:
!
!      integral ( 0 <= x <= pi ) sin(x) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) + pi/2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gwynne Evans,
!    Practical Numerical Integration,
!    Wiley, 1993,
!    ISBN: 047193898X,
!    LC: QA299.3E93.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1, 2, 4, 8 or 16.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = 0.0D+00

    w(1) = 2.0D+00

  else if ( n == 2 ) then

    x(1) = - 0.68366739008990304094D+00
    x(2) =   0.68366739008990304094D+00

    w(1) = 1.0D+00
    w(2) = 1.0D+00

  else if ( n == 4 ) then

    x(1) = - 1.1906765638948557415D+00
    x(2) = - 0.43928746686001514756D+00
    x(3) =   0.43928746686001514756D+00
    x(4) =   1.1906765638948557415D+00

    w(1) = 0.22407061812762016065D+00
    w(2) = 0.77592938187237983935D+00
    w(3) = 0.77592938187237983935D+00
    w(4) = 0.22407061812762016065D+00

  else if ( n == 8 ) then

    x(1) = - 1.4414905401823575701D+00
    x(2) = - 1.1537256454567275850D+00
    x(3) = - 0.74346864787549244989D+00
    x(4) = - 0.25649650741623123020D+00
    x(5) =   0.25649650741623123020D+00
    x(6) =   0.74346864787549244989D+00
    x(7) =   1.1537256454567275850D+00
    x(8) =   1.4414905401823575701D+00

    w(1) = 0.027535633513767011149D+00
    w(2) = 0.14420409203022750950D+00
    w(3) = 0.33626447785280459621D+00
    w(4) = 0.49199579660320088314D+00
    w(5) = 0.49199579660320088314D+00
    w(6) = 0.33626447785280459621D+00
    w(7) = 0.14420409203022750950D+00
    w(8) = 0.027535633513767011149D+00

  else if ( n == 16 ) then

    x( 1) = - 1.5327507132362304779D+00
    x( 2) = - 1.4446014873666514608D+00
    x( 3) = - 1.3097818904452936698D+00
    x( 4) = - 1.1330068786005003695D+00
    x( 5) = - 0.92027786206637096497D+00
    x( 6) = - 0.67861108097560545347D+00
    x( 7) = - 0.41577197673418943962D+00
    x( 8) = - 0.14003444424696773778D+00
    x( 9) =   0.14003444424696773778D+00
    x(10) =   0.41577197673418943962D+00
    x(11) =   0.67861108097560545347D+00
    x(12) =   0.92027786206637096497D+00
    x(13) =   1.1330068786005003695D+00
    x(14) =   1.3097818904452936698D+00
    x(15) =   1.4446014873666514608D+00
    x(16) =   1.5327507132362304779D+00

    w( 1) = 0.0024194677567615628193D+00
    w( 2) = 0.014115268156854008264D+00
    w( 3) = 0.040437893946503669410D+00
    w( 4) = 0.083026647573217742131D+00
    w( 5) = 0.13834195526951273359D+00
    w( 6) = 0.19741148870253455567D+00
    w( 7) = 0.24763632094635522403D+00
    w( 8) = 0.27661095764826050408D+00
    w( 9) = 0.27661095764826050408D+00
    w(10) = 0.24763632094635522403D+00
    w(11) = 0.19741148870253455567D+00
    w(12) = 0.13834195526951273359D+00
    w(13) = 0.083026647573217742131D+00
    w(14) = 0.040437893946503669410D+00
    w(15) = 0.014115268156854008264D+00
    w(16) = 0.0024194677567615628193D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COS_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  end if

  return
end
subroutine legendre_dr_compute ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_DR_COMPUTE: Gauss-Legendre quadrature by Davis-Rabinowitz method.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
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
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer ( kind = 4 ) N, the order.
!    0 < N.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

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
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xtemp

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_DR_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop
  end if

  e1 = real ( n * ( n + 1 ), kind = 8 )

  m = ( n + 1 ) / 2

  do i = 1, m

    mp1mi = m + 1 - i

    t = real ( 4 * i - 1, kind = 8 ) * pi &
      / real ( 4 * n + 2, kind = 8 )

    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( n, kind = 8 ) ) &
      / real ( 8 * n * n, kind = 8 ) )

    pkm1 = 1.0D+00
    pk = x0

    do k = 2, n
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do

    d1 = real ( n, kind = 8 ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 ) &
      / ( 1.0D+00 + x0 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) &
      / ( 1.0D+00 - x0 ) / ( 1.0D+00 + x0 )

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

    x(mp1mi) = xtemp

    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

    w(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp ) * ( 1.0D+00 + xtemp ) / ( fx * fx )

  end do

  if ( mod ( n, 2 ) == 1 ) then
    x(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = ( n + 1 ) / 2
  ncopy = n - nmove

  do i = 1, nmove
    iback = n + 1 - i
    x(iback) = x(iback-ncopy)
    w(iback) = w(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, n - nmove
    x(i) = - x(n+1-i)
    w(i) = w(n+1-i)
  end do

  return
end
subroutine legendre_ek_compute ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_EK_COMPUTE: Legendre quadrature rule by the Elhay-Kautsky method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) zemu
!
!  Define the zero-th moment.
!
  zemu = 2.0D+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i * i, kind = 8 ) / real ( 4 * i * i - 1, kind = 8 )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0D+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine legendre_gw_compute ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_GW_COMPUTE: Legendre quadrature rule by the Golub-Welsch method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Gene Golub, John Welsch,
!    Calculation of Gaussian Quadrature Rules,
!    Mathematics of Computation,
!    Volume 23, Number 106, April 1969, pages 221-230.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) beta(n)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) z(n,n)

  do i = 1, n - 1
    beta(i) = 1.0D+00 / real ( 2 * i, kind = 8 )
  end do
  beta(1:n-1) = 0.5D+00 &
   / sqrt ( ( 1.0D+00 - beta(1:n-1) ) * ( 1.0D+00 + beta(1:n-1) ) )
!
!  Compute eigenvalues and eigenvectors.
!
  d(1:n) = 0.0D+00
  beta(2:n) = beta(1:n-1)
  beta(1) = 0.0D+00

  z(1:n,1:n) = 0.0D+00
  do i = 1, n
    z(i,i) = 1.0D+00
  end do

  call imtql2 ( n, d, beta, z, ierr )
!
!  X values are eigenvalues.
!
  x(1:n) = d(1:n)
!
!  W is related to first eigenvector.
!
  w(1:n) = 2.0D+00 * z(1,1:n)**2

  w(1:n) = 2.0D+00 * w(1:n) / sum ( w(1:n) )

  return
end
subroutine legendre_integral ( expon, exact )

!*****************************************************************************80
!
!! LEGENDRE_INTEGRAL evaluates a monomial Legendre integral.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= +1 ) x^n dx
!
!    This routine is given the value of the exponent, and returns the
!    exact value of the integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Output, real ( kind = 8 ) EXACT, the value of the exact integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
!
!  Get the exact value of the integral.
!
  if ( mod ( expon, 2 ) == 0 ) then

    exact = 2.0D+00 / real ( expon + 1, kind = 8 )
    
  else

    exact = 0.0D+00
    
  end if

  return
end
subroutine legendre_log_compute ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_LOG_COMPUTE: Gauss-Legendre rules for - LOG(X) * F(X) on [0,1].
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x <= 1 ) - log(x) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    Note that this rule is reasonably accurate for values of N up to 13.
!    However, for larger values of N, the weight calculation loses so much
!    accuracy that the formula is not reliable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2007
!
!  Author:
!
!    Original FORTRAN77 version by Federico Paris, Jose Canas.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Donald Anderson,
!    Gaussian Quadrature Formulae for the integral from 0 to 1 of
!    -ln(X) f(X) dx,
!    Mathematics of Computation,
!    Volume 19, Number 91, Number 3, July 1965, pages 477-481.
!
!    Federico Paris, Jose Canas,
!    Boundary Element Method: Fundamentals and Applications,
!    Oxford, 1997,
!    ISBN: 0-19-856543-7
!    LC: TA347.B69.P34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) coef(-1:n,3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ir0
  integer ( kind = 4 ) ir1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) rk
  real ( kind = 8 ) root(0:n+1,2)
  real ( kind = 8 ) sk
  real ( kind = 8 ) tk
  real ( kind = 8 ) tk1
  real ( kind = 8 ) tol
  real ( kind = 8 ) uk
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = 0.25D+00
    w(1) = 1.0D+00
    return

  end if

  coef(-1:n,1:3) = 0.0D+00
  root(0:n,1:2) = 0.0D+00

  coef(-1,1) = 0.0D+00
  coef( 0,1) = 1.0D+00

  coef(-1,2) =  0.0D+00
  coef( 0,2) = -0.5D+00
  coef( 1,2) =  2.0D+00

  tk1 = 1.0D+00

  tol = 1.0D-12
!  tol = 1.0D-10

  i0 = 3
  i1 = 2
  i2 = 1

  ir1 = 1
  ir0 = 2

  root(0,1) = 0.0D+00
  root(1,1) = 0.25D+00
  root(2,1) = 1.0D+00
!
!  Computation of polynomial coefficients and roots.
!
  do k = 2, n
    tk = 0.0D+00
    uk = 0.0D+00
    do i = 0, k - 1
      do j = 0, k - 1
        tk = tk + coef(i,i1) * coef(j,i1) / ( real ( i + j + 1, kind = 8 ) )**2
        uk = uk + coef(i,i1) * coef(j,i1) / ( real ( i + j + 2, kind = 8 ) )**2
      end do
    end do
    rk = uk / tk
    sk = tk / tk1
    do i = 0, k - 1
      coef(i,i0) = 2.0D+00 &
        * ( coef(i-1,i1) - rk * coef(i,i1) ) - sk * coef(i,i2)
    end do
    coef(k,i0) = 2.0D+00**k

    do i = 1, k
      call legendre_log_roots ( k, coef(0:k,i0), root(i-1,ir1), root(i,ir1), &
        tol, root(i,ir0) )
    end do
    root(  0,ir0) = 0.0D+00
    root(k+1,ir0) = 1.0D+00

    ii = i0
    i0 = i2
    i2 = i1
    i1 = ii

    ii = ir0
    ir0 = ir1
    ir1 = ii

    tk1 = tk

  end do
!
!  Compute the abscissas and weights.
!
  x(1:n) = root(1:n,ir1)

  do i = 1, n
    j = 0
    v1 = coef(j+1,i1)
    v2 = coef(j,  i2)
    do j = 1, n - 1
      v1 = v1 + real ( j + 1, kind = 8 ) * coef(j+1,i1) * x(i)**j
      v2 = v2 +                            coef(j,  i2) * x(i)**j
    end do
    w(i) = 2.0D+00 * tk / ( v1 * v2 )
  end do

  return
end
subroutine legendre_log_roots ( n, coef, x00, xff, tol, xi )

!*****************************************************************************80
!
!! LEGENDRE_LOG_ROOTS is a root finder for LEGENDRE_LOG_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2007
!
!  Author:
!
!    Original FORTRAN77 version by Federico Paris, Jose Canas.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Federico Paris, Jose Canas,
!    Boundary Element Method: Fundamentals and Applications,
!    Oxford, 1997,
!    ISBN: 0-19-856543-7
!    LC: TA347.B69.P34.
!
!  Parameters:
!
!    Input, integer N, the degree of the polynomial.
!
!    Input, real ( kind = 8 ) COEF(N+1), the coefficients of the polynomial.
!
!    Input, real ( kind = 8 ) X00, XFF, two points at which the polynomial
!    has opposite sign.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for the size of the
!    uncertainty interval in the value of XI.
!
!    Output, real ( kind = 8 ) XI, the approximate root.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) coef(n+1)
  real ( kind = 8 ) errac
  real ( kind = 8 ) tol
  real ( kind = 8 ) v0
  real ( kind = 8 ) vf
  real ( kind = 8 ) vi
  real ( kind = 8 ) x0
  real ( kind = 8 ) x00
  real ( kind = 8 ) xf
  real ( kind = 8 ) xff
  real ( kind = 8 ) xi

  x0 = x00
  call legendre_log_value ( n, coef, x0, v0 )

  xf = xff
  call legendre_log_value ( n, coef, xf, vf )

  errac = abs ( xf - x0 )

  do while ( tol .lt. errac )

    xi = 0.5D+00 * ( x0 + xf )

    call legendre_log_value ( n, coef, xi, vi )

    if ( abs ( vi ) .lt. tol ) then
      return
    end if

    if ( vi * v0 .lt. 0.0D+00 ) then
      xf = xi
      vf = vi
    else
      x0 = xi
      v0 = vi
    end if

    errac = abs ( xf - x0 )

  end do

  return
end
subroutine legendre_log_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_LOG_SET sets a Gauss-Legendre rule for - LOG(X) * F(X) on [0,1].
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x <= 1 ) - log(x) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Gwynne Evans,
!    Practical Numerical Integration,
!    Wiley, 1993,
!    ISBN: 047193898X,
!    LC: QA299.3E93.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 through 8, or 16.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = 0.25D+00

    w(1) = 1.0D+00

  else if ( n == 2 ) then

    x(1) = 0.112008806166976182957205488948D+00
    x(2) = 0.602276908118738102757080225338D+00

    w(1) = 0.718539319030384440665510200891D+00
    w(2) = 0.281460680969615559334489799109D+00

  else if ( n == 3 ) then

    x(1) = 0.0638907930873254049961166031363D+00
    x(2) = 0.368997063715618765546197645857D+00
    x(3) = 0.766880303938941455423682659817D+00

    w(1) = 0.513404552232363325129300497567D+00
    w(2) = 0.391980041201487554806287180966D+00
    w(3) = 0.0946154065661491200644123214672D+00

  else if ( n == 4 ) then

    x(1) = 0.0414484801993832208033213101564D+00
    x(2) = 0.245274914320602251939675759523D+00
    x(3) = 0.556165453560275837180184354376D+00
    x(4) = 0.848982394532985174647849188085D+00

    w(1) = 0.383464068145135124850046522343D+00
    w(2) = 0.386875317774762627336008234554D+00
    w(3) = 0.190435126950142415361360014547D+00
    w(4) = 0.0392254871299598324525852285552D+00

  else if ( n == 5 ) then

    x(1) = 0.0291344721519720533037267621154D+00
    x(2) = 0.173977213320897628701139710829D+00
    x(3) = 0.411702520284902043174931924646D+00
    x(4) = 0.677314174582820380701802667998D+00
    x(5) = 0.894771361031008283638886204455D+00

    w(1) = 0.297893471782894457272257877888D+00
    w(2) = 0.349776226513224180375071870307D+00
    w(3) = 0.234488290044052418886906857943D+00
    w(4) = 0.0989304595166331469761807114404D+00
    w(5) = 0.0189115521431957964895826824218D+00

  else if ( n == 6 ) then

    x(1) = 0.0216340058441169489956958558537D+00
    x(2) = 0.129583391154950796131158505009D+00
    x(3) = 0.314020449914765508798248188420D+00
    x(4) = 0.538657217351802144548941893993D+00
    x(5) = 0.756915337377402852164544156139D+00
    x(6) = 0.922668851372120237333873231507D+00

    w(1) = 0.238763662578547569722268303330D+00
    w(2) = 0.308286573273946792969383109211D+00
    w(3) = 0.245317426563210385984932540188D+00
    w(4) = 0.142008756566476685421345576030D+00
    w(5) = 0.0554546223248862900151353549662D+00
    w(6) = 0.0101689586929322758869351162755D+00

  else if ( n == 7 ) then

    x(1) = 0.0167193554082585159416673609320D+00
    x(2) = 0.100185677915675121586885031757D+00
    x(3) = 0.246294246207930599046668547239D+00
    x(4) = 0.433463493257033105832882482601D+00
    x(5) = 0.632350988047766088461805812245D+00
    x(6) = 0.811118626740105576526226796782D+00
    x(7) = 0.940848166743347721760134113379D+00

    w(1) = 0.196169389425248207525427377585D+00
    w(2) = 0.270302644247272982145271719533D+00
    w(3) = 0.239681873007690948308072785252D+00
    w(4) = 0.165775774810432906560869687736D+00
    w(5) = 0.0889432271376579644357238403458D+00
    w(6) = 0.0331943043565710670254494111034D+00
    w(7) = 0.00593278701512592399918517844468D+00

  else if ( n == 8 ) then

    x(1) = 0.0133202441608924650122526725243D+00
    x(2) = 0.0797504290138949384098277291424D+00
    x(3) = 0.197871029326188053794476159516D+00
    x(4) = 0.354153994351909419671463603538D+00
    x(5) = 0.529458575234917277706149699996D+00
    x(6) = 0.701814529939099963837152670310D+00
    x(7) = 0.849379320441106676048309202301D+00
    x(8) = 0.953326450056359788767379678514D+00

    w(1) = 0.164416604728002886831472568326D+00
    w(2) = 0.237525610023306020501348561960D+00
    w(3) = 0.226841984431919126368780402936D+00
    w(4) = 0.175754079006070244988056212006D+00
    w(5) = 0.112924030246759051855000442086D+00
    w(6) = 0.0578722107177820723985279672940D+00
    w(7) = 0.0209790737421329780434615241150D+00
    w(8) = 0.00368640710402761901335232127647D+00

  else if ( n == 16 ) then

    x( 1) = 0.00389783448711591592405360527037D+00
    x( 2) = 0.0230289456168732398203176309848D+00
    x( 3) = 0.0582803983062404123483532298394D+00
    x( 4) = 0.108678365091054036487713613051D+00
    x( 5) = 0.172609454909843937760843776232D+00
    x( 6) = 0.247937054470578495147671753047D+00
    x( 7) = 0.332094549129917155984755859320D+00
    x( 8) = 0.422183910581948600115088366745D+00
    x( 9) = 0.515082473381462603476277704052D+00
    x(10) = 0.607556120447728724086384921709D+00
    x(11) = 0.696375653228214061156318166581D+00
    x(12) = 0.778432565873265405203868167732D+00
    x(13) = 0.850850269715391083233822761319D+00
    x(14) = 0.911086857222271905418818994060D+00
    x(15) = 0.957025571703542157591520509383D+00
    x(16) = 0.987047800247984476758697436516D+00

    w( 1) = 0.0607917100435912328511733871235D+00
    w( 2) = 0.102915677517582144387691736210D+00
    w( 3) = 0.122355662046009193557547513197D+00
    w( 4) = 0.127569246937015988717042209239D+00
    w( 5) = 0.123013574600070915423123365137D+00
    w( 6) = 0.111847244855485722621848903429D+00
    w( 7) = 0.0965963851521243412529681650802D+00
    w( 8) = 0.0793566643514731387824416770520D+00
    w( 9) = 0.0618504945819652070951360853113D+00
    w(10) = 0.0454352465077266686288299526509D+00
    w(11) = 0.0310989747515818064092528627927D+00
    w(12) = 0.0194597659273608420780860268669D+00
    w(13) = 0.0107762549632055256455393162159D+00
    w(14) = 0.00497254289008764171250524951646D+00
    w(15) = 0.00167820111005119451503546419059D+00
    w(16) = 0.000282353764668436321778085987413D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_LOG_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 through 8 and 16.'
    stop

  end if

  return
end
subroutine legendre_log_value ( n, coef, x, value )

!*****************************************************************************80
!
!! LEGENDRE_LOG_VALUE evaluates a polynomial for LEGENDRE_LOG_ROOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 December 2007
!
!  Author:
!
!    Original FORTRAN77 version by Federico Paris, Jose Canas.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Federico Paris, Jose Canas,
!    Boundary Element Method: Fundamentals and Applications,
!    Oxford, 1997,
!    ISBN: 0-19-856543-7
!    LC: TA347.B69.P34.
!
!  Parameters:
!
!    Input, integer N, the degree of the polynomial.
!
!    Input, real ( kind = 8 ) COEF(N+1), the polynomial coefficients.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) V, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) coef(n+1)
  integer ( kind = 4 ) k
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = coef(n+1)
  do k = n, 1, -1
    value = value * x + coef(k)
  end do

  return
end
subroutine legendre_polynomial_value ( m, n, x, p )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_VALUE evaluates a Legendre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of values.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) P(M), the value of P(N)(X(1:M)).
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) dp0(m)
  real ( kind = 8 ) dp1(m)
  real ( kind = 8 ) dp2(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) p(m)
  real ( kind = 8 ) p0(m)
  real ( kind = 8 ) p1(m)
  real ( kind = 8 ) p2(m)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    p(1:m) = 0.0D+00
    return
  end if

  p1(1:m) = 1.0D+00
  dp1(1:m) = 0.0D+00

  if ( n == 0 ) then
    p(1:m) = p1(1:m)
    return
  end if

  p2(1:m) = x(1:m)
  dp2(1:m) = 1.0D+00

  if ( n == 1 ) then
    p(1:m) = p2(1:m)
    return
  end if

  do i = 2, n

    p0(1:m) = p1(1:m)
    dp0(1:m) = dp1(1:m)

    p1(1:m) = p2(1:m)
    dp1(1:m) = dp2(1:m)

    p2(1:m) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * p1(1:m) &
              + real (   - i + 1, kind = 8 )          * p0(1:m) ) &
              / real (     i,     kind = 8 )

    dp2(1:m) = &
      ( real ( 2 * i - 1, kind = 8 ) * ( p1(1:m) + x(1:m) * dp1(1:m) ) &
      - real (     i - 1, kind = 8 ) * dp0(1:m) ) &
      / real (     i,     kind = 8 )

  end do

  p(1:m) = p2(1:m)

  return
end
subroutine legendre_recur ( p2, dp2, p1, x, n )

!*****************************************************************************80
!
!! LEGENDRE_RECUR finds the value and derivative of a Legendre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of P(N)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of P'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of P(N-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, n

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( real ( 2 * i - 1, kind = 8 ) * x * p1 &
         + real (   - i + 1, kind = 8 )     * p0 ) &
         / real (     i,     kind = 8 )

    dp2 = ( real ( 2 * i - 1, kind = 8 ) * ( p1 + x * dp1 ) &
          - real (     i - 1, kind = 8 ) * dp0 ) &
          / real (     i,     kind = 8 )

  end do

  return
end
subroutine legendre_root ( x, n, dp2, p1 )

!*****************************************************************************80
!
!! LEGENDRE_ROOT improves an approximate root of a Legendre polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) DP2, the value of L'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(N-1)(X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = epsilon ( eps )

  do step = 1, step_max

    call legendre_recur ( p2, dp2, p1, x, n )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine legendre_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The quadrature rule is exact for polynomials through degree 2*N-1.
!
!    The abscissas are the zeroes of the Legendre polynomial P(N)(X).
!
!    Mathematica can compute the abscissas and weights of a Gauss-Legendre
!    rule of order N for the interval [A,B] with P digits of precision
!    by the commands:
!
!    Needs["NumericalDifferentialEquationAnalysis`"]
!    GaussianQuadratureWeights [ n, a, b, p ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    Dover, 2006,
!    ISBN: 0486445798,
!    LC: QA311.K713.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 33 or 63/64/65, 127/128/129,
!    255/256/257.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = 0.000000000000000000000000000000D+00

    w(1) = 2.000000000000000000000000000000D+00

  else if ( n == 2 ) then

    x(1) = -0.577350269189625764509148780502D+00
    x(2) = 0.577350269189625764509148780502D+00

    w(1) = 1.000000000000000000000000000000D+00
    w(2) = 1.000000000000000000000000000000D+00

  else if ( n == 3 ) then

    x(1) = -0.774596669241483377035853079956D+00
    x(2) = 0.000000000000000000000000000000D+00
    x(3) = 0.774596669241483377035853079956D+00

    w(1) = 0.555555555555555555555555555556D+00
    w(2) = 0.888888888888888888888888888889D+00
    w(3) = 0.555555555555555555555555555556D+00

  else if ( n == 4 ) then

    x(1) = -0.861136311594052575223946488893D+00
    x(2) = -0.339981043584856264802665759103D+00
    x(3) = 0.339981043584856264802665759103D+00
    x(4) = 0.861136311594052575223946488893D+00

    w(1) = 0.347854845137453857373063949222D+00
    w(2) = 0.652145154862546142626936050778D+00
    w(3) = 0.652145154862546142626936050778D+00
    w(4) = 0.347854845137453857373063949222D+00

  else if ( n == 5 ) then

    x(1) = -0.906179845938663992797626878299D+00
    x(2) = -0.538469310105683091036314420700D+00
    x(3) = 0.000000000000000000000000000000D+00
    x(4) = 0.538469310105683091036314420700D+00
    x(5) = 0.906179845938663992797626878299D+00

    w(1) = 0.236926885056189087514264040720D+00
    w(2) = 0.478628670499366468041291514836D+00
    w(3) = 0.568888888888888888888888888889D+00
    w(4) = 0.478628670499366468041291514836D+00
    w(5) = 0.236926885056189087514264040720D+00

  else if ( n == 6 ) then

    x(1) = -0.932469514203152027812301554494D+00
    x(2) = -0.661209386466264513661399595020D+00
    x(3) = -0.238619186083196908630501721681D+00
    x(4) = 0.238619186083196908630501721681D+00
    x(5) = 0.661209386466264513661399595020D+00
    x(6) = 0.932469514203152027812301554494D+00

    w(1) = 0.171324492379170345040296142173D+00
    w(2) = 0.360761573048138607569833513838D+00
    w(3) = 0.467913934572691047389870343990D+00
    w(4) = 0.467913934572691047389870343990D+00
    w(5) = 0.360761573048138607569833513838D+00
    w(6) = 0.171324492379170345040296142173D+00

  else if ( n == 7 ) then

    x(1) = -0.949107912342758524526189684048D+00
    x(2) = -0.741531185599394439863864773281D+00
    x(3) = -0.405845151377397166906606412077D+00
    x(4) = 0.000000000000000000000000000000D+00
    x(5) = 0.405845151377397166906606412077D+00
    x(6) = 0.741531185599394439863864773281D+00
    x(7) = 0.949107912342758524526189684048D+00

    w(1) = 0.129484966168869693270611432679D+00
    w(2) = 0.279705391489276667901467771424D+00
    w(3) = 0.381830050505118944950369775489D+00
    w(4) = 0.417959183673469387755102040816D+00
    w(5) = 0.381830050505118944950369775489D+00
    w(6) = 0.279705391489276667901467771424D+00
    w(7) = 0.129484966168869693270611432679D+00

  else if ( n == 8 ) then

    x(1) = -0.960289856497536231683560868569D+00
    x(2) = -0.796666477413626739591553936476D+00
    x(3) = -0.525532409916328985817739049189D+00
    x(4) = -0.183434642495649804939476142360D+00
    x(5) = 0.183434642495649804939476142360D+00
    x(6) = 0.525532409916328985817739049189D+00
    x(7) = 0.796666477413626739591553936476D+00
    x(8) = 0.960289856497536231683560868569D+00

    w(1) = 0.101228536290376259152531354310D+00
    w(2) = 0.222381034453374470544355994426D+00
    w(3) = 0.313706645877887287337962201987D+00
    w(4) = 0.362683783378361982965150449277D+00
    w(5) = 0.362683783378361982965150449277D+00
    w(6) = 0.313706645877887287337962201987D+00
    w(7) = 0.222381034453374470544355994426D+00
    w(8) = 0.101228536290376259152531354310D+00

  else if ( n == 9 ) then

    x(1) = -0.968160239507626089835576203D+00
    x(2) = -0.836031107326635794299429788D+00
    x(3) = -0.613371432700590397308702039D+00
    x(4) = -0.324253423403808929038538015D+00
    x(5) = 0.000000000000000000000000000D+00
    x(6) = 0.324253423403808929038538015D+00
    x(7) = 0.613371432700590397308702039D+00
    x(8) = 0.836031107326635794299429788D+00
    x(9) = 0.968160239507626089835576203D+00

    w(1) = 0.081274388361574411971892158111D+00
    w(2) = 0.18064816069485740405847203124D+00
    w(3) = 0.26061069640293546231874286942D+00
    w(4) = 0.31234707704000284006863040658D+00
    w(5) = 0.33023935500125976316452506929D+00
    w(6) = 0.31234707704000284006863040658D+00
    w(7) = 0.26061069640293546231874286942D+00
    w(8) = 0.18064816069485740405847203124D+00
    w(9) = 0.081274388361574411971892158111D+00

  else if ( n == 10 ) then

    x(1) = -0.973906528517171720077964012D+00
    x(2) = -0.865063366688984510732096688D+00
    x(3) = -0.679409568299024406234327365D+00
    x(4) = -0.433395394129247190799265943D+00
    x(5) = -0.148874338981631210884826001D+00
    x(6) = 0.148874338981631210884826001D+00
    x(7) = 0.433395394129247190799265943D+00
    x(8) = 0.679409568299024406234327365D+00
    x(9) = 0.865063366688984510732096688D+00
    x(10) = 0.973906528517171720077964012D+00

    w(1) = 0.066671344308688137593568809893D+00
    w(2) = 0.14945134915058059314577633966D+00
    w(3) = 0.21908636251598204399553493423D+00
    w(4) = 0.26926671930999635509122692157D+00
    w(5) = 0.29552422471475287017389299465D+00
    w(6) = 0.29552422471475287017389299465D+00
    w(7) = 0.26926671930999635509122692157D+00
    w(8) = 0.21908636251598204399553493423D+00
    w(9) = 0.14945134915058059314577633966D+00
    w(10) = 0.066671344308688137593568809893D+00

  else if ( n == 11 ) then

    x(1) = -0.978228658146056992803938001D+00
    x(2) = -0.887062599768095299075157769D+00
    x(3) = -0.730152005574049324093416252D+00
    x(4) = -0.519096129206811815925725669D+00
    x(5) = -0.269543155952344972331531985D+00
    x(6) = 0.000000000000000000000000000D+00
    x(7) = 0.269543155952344972331531985D+00
    x(8) = 0.519096129206811815925725669D+00
    x(9) = 0.730152005574049324093416252D+00
    x(10) = 0.887062599768095299075157769D+00
    x(11) = 0.978228658146056992803938001D+00

    w(1) = 0.055668567116173666482753720443D+00
    w(2) = 0.12558036946490462463469429922D+00
    w(3) = 0.18629021092773425142609764143D+00
    w(4) = 0.23319376459199047991852370484D+00
    w(5) = 0.26280454451024666218068886989D+00
    w(6) = 0.27292508677790063071448352834D+00
    w(7) = 0.26280454451024666218068886989D+00
    w(8) = 0.23319376459199047991852370484D+00
    w(9) = 0.18629021092773425142609764143D+00
    w(10) = 0.12558036946490462463469429922D+00
    w(11) = 0.055668567116173666482753720443D+00

  else if ( n == 12 ) then

    x(1) = -0.981560634246719250690549090D+00
    x(2) = -0.904117256370474856678465866D+00
    x(3) = -0.769902674194304687036893833D+00
    x(4) = -0.587317954286617447296702419D+00
    x(5) = -0.367831498998180193752691537D+00
    x(6) = -0.125233408511468915472441369D+00
    x(7) = 0.125233408511468915472441369D+00
    x(8) = 0.367831498998180193752691537D+00
    x(9) = 0.587317954286617447296702419D+00
    x(10) = 0.769902674194304687036893833D+00
    x(11) = 0.904117256370474856678465866D+00
    x(12) = 0.981560634246719250690549090D+00

    w(1) = 0.047175336386511827194615961485D+00
    w(2) = 0.10693932599531843096025471819D+00
    w(3) = 0.16007832854334622633465252954D+00
    w(4) = 0.20316742672306592174906445581D+00
    w(5) = 0.23349253653835480876084989892D+00
    w(6) = 0.24914704581340278500056243604D+00
    w(7) = 0.24914704581340278500056243604D+00
    w(8) = 0.23349253653835480876084989892D+00
    w(9) = 0.20316742672306592174906445581D+00
    w(10) = 0.16007832854334622633465252954D+00
    w(11) = 0.10693932599531843096025471819D+00
    w(12) = 0.047175336386511827194615961485D+00

  else if ( n == 13 ) then

    x(1) = -0.984183054718588149472829449D+00
    x(2) = -0.917598399222977965206547837D+00
    x(3) = -0.801578090733309912794206490D+00
    x(4) = -0.642349339440340220643984607D+00
    x(5) = -0.448492751036446852877912852D+00
    x(6) = -0.230458315955134794065528121D+00
    x(7) = 0.000000000000000000000000000D+00
    x(8) = 0.230458315955134794065528121D+00
    x(9) = 0.448492751036446852877912852D+00
    x(10) = 0.642349339440340220643984607D+00
    x(11) = 0.80157809073330991279420649D+00
    x(12) = 0.91759839922297796520654784D+00
    x(13) = 0.98418305471858814947282945D+00

    w(1) = 0.040484004765315879520021592201D+00
    w(2) = 0.092121499837728447914421775954D+00
    w(3) = 0.13887351021978723846360177687D+00
    w(4) = 0.17814598076194573828004669200D+00
    w(5) = 0.20781604753688850231252321931D+00
    w(6) = 0.22628318026289723841209018604D+00
    w(7) = 0.23255155323087391019458951527D+00
    w(8) = 0.22628318026289723841209018604D+00
    w(9) = 0.20781604753688850231252321931D+00
    w(10) = 0.17814598076194573828004669200D+00
    w(11) = 0.13887351021978723846360177687D+00
    w(12) = 0.092121499837728447914421775954D+00
    w(13) = 0.040484004765315879520021592201D+00

  else if ( n == 14 ) then

    x(1) = -0.986283808696812338841597267D+00
    x(2) = -0.928434883663573517336391139D+00
    x(3) = -0.827201315069764993189794743D+00
    x(4) = -0.687292904811685470148019803D+00
    x(5) = -0.515248636358154091965290719D+00
    x(6) = -0.319112368927889760435671824D+00
    x(7) = -0.108054948707343662066244650D+00
    x(8) = 0.108054948707343662066244650D+00
    x(9) = 0.31911236892788976043567182D+00
    x(10) = 0.51524863635815409196529072D+00
    x(11) = 0.68729290481168547014801980D+00
    x(12) = 0.82720131506976499318979474D+00
    x(13) = 0.92843488366357351733639114D+00
    x(14) = 0.98628380869681233884159727D+00

    w(1) = 0.035119460331751863031832876138D+00
    w(2) = 0.08015808715976020980563327706D+00
    w(3) = 0.12151857068790318468941480907D+00
    w(4) = 0.15720316715819353456960193862D+00
    w(5) = 0.18553839747793781374171659013D+00
    w(6) = 0.20519846372129560396592406566D+00
    w(7) = 0.21526385346315779019587644332D+00
    w(8) = 0.21526385346315779019587644332D+00
    w(9) = 0.20519846372129560396592406566D+00
    w(10) = 0.18553839747793781374171659013D+00
    w(11) = 0.15720316715819353456960193862D+00
    w(12) = 0.12151857068790318468941480907D+00
    w(13) = 0.08015808715976020980563327706D+00
    w(14) = 0.035119460331751863031832876138D+00

  else if ( n == 15 ) then

    x(1) = -0.987992518020485428489565719D+00
    x(2) = -0.937273392400705904307758948D+00
    x(3) = -0.848206583410427216200648321D+00
    x(4) = -0.724417731360170047416186055D+00
    x(5) = -0.570972172608538847537226737D+00
    x(6) = -0.394151347077563369897207371D+00
    x(7) = -0.201194093997434522300628303D+00
    x(8) = 0.00000000000000000000000000D+00
    x(9) = 0.20119409399743452230062830D+00
    x(10) = 0.39415134707756336989720737D+00
    x(11) = 0.57097217260853884753722674D+00
    x(12) = 0.72441773136017004741618605D+00
    x(13) = 0.84820658341042721620064832D+00
    x(14) = 0.93727339240070590430775895D+00
    x(15) = 0.98799251802048542848956572D+00

    w(1) = 0.030753241996117268354628393577D+00
    w(2) = 0.070366047488108124709267416451D+00
    w(3) = 0.107159220467171935011869546686D+00
    w(4) = 0.13957067792615431444780479451D+00
    w(5) = 0.16626920581699393355320086048D+00
    w(6) = 0.18616100001556221102680056187D+00
    w(7) = 0.19843148532711157645611832644D+00
    w(8) = 0.20257824192556127288062019997D+00
    w(9) = 0.19843148532711157645611832644D+00
    w(10) = 0.18616100001556221102680056187D+00
    w(11) = 0.16626920581699393355320086048D+00
    w(12) = 0.13957067792615431444780479451D+00
    w(13) = 0.107159220467171935011869546686D+00
    w(14) = 0.070366047488108124709267416451D+00
    w(15) = 0.030753241996117268354628393577D+00

  else if ( n == 16 ) then

    x(1) = -0.989400934991649932596154173D+00
    x(2) = -0.944575023073232576077988416D+00
    x(3) = -0.865631202387831743880467898D+00
    x(4) = -0.755404408355003033895101195D+00
    x(5) = -0.617876244402643748446671764D+00
    x(6) = -0.458016777657227386342419443D+00
    x(7) = -0.281603550779258913230460501D+00
    x(8) = -0.09501250983763744018531934D+00
    x(9) = 0.09501250983763744018531934D+00
    x(10) = 0.28160355077925891323046050D+00
    x(11) = 0.45801677765722738634241944D+00
    x(12) = 0.61787624440264374844667176D+00
    x(13) = 0.75540440835500303389510119D+00
    x(14) = 0.86563120238783174388046790D+00
    x(15) = 0.94457502307323257607798842D+00
    x(16) = 0.98940093499164993259615417D+00

    w(1) = 0.027152459411754094851780572456D+00
    w(2) = 0.062253523938647892862843836994D+00
    w(3) = 0.09515851168249278480992510760D+00
    w(4) = 0.12462897125553387205247628219D+00
    w(5) = 0.14959598881657673208150173055D+00
    w(6) = 0.16915651939500253818931207903D+00
    w(7) = 0.18260341504492358886676366797D+00
    w(8) = 0.18945061045506849628539672321D+00
    w(9) = 0.18945061045506849628539672321D+00
    w(10) = 0.18260341504492358886676366797D+00
    w(11) = 0.16915651939500253818931207903D+00
    w(12) = 0.14959598881657673208150173055D+00
    w(13) = 0.12462897125553387205247628219D+00
    w(14) = 0.09515851168249278480992510760D+00
    w(15) = 0.062253523938647892862843836994D+00
    w(16) = 0.027152459411754094851780572456D+00

  else if ( n == 17 ) then

    x(1) = -0.990575475314417335675434020D+00
    x(2) = -0.950675521768767761222716958D+00
    x(3) = -0.880239153726985902122955694D+00
    x(4) = -0.781514003896801406925230056D+00
    x(5) = -0.657671159216690765850302217D+00
    x(6) = -0.512690537086476967886246569D+00
    x(7) = -0.35123176345387631529718552D+00
    x(8) = -0.17848418149584785585067749D+00
    x(9) = 0.00000000000000000000000000D+00
    x(10) = 0.17848418149584785585067749D+00
    x(11) = 0.35123176345387631529718552D+00
    x(12) = 0.51269053708647696788624657D+00
    x(13) = 0.65767115921669076585030222D+00
    x(14) = 0.78151400389680140692523006D+00
    x(15) = 0.88023915372698590212295569D+00
    x(16) = 0.95067552176876776122271696D+00
    x(17) = 0.99057547531441733567543402D+00

    w(1) = 0.024148302868547931960110026288D+00
    w(2) = 0.055459529373987201129440165359D+00
    w(3) = 0.085036148317179180883535370191D+00
    w(4) = 0.111883847193403971094788385626D+00
    w(5) = 0.13513636846852547328631998170D+00
    w(6) = 0.15404576107681028808143159480D+00
    w(7) = 0.16800410215645004450997066379D+00
    w(8) = 0.17656270536699264632527099011D+00
    w(9) = 0.17944647035620652545826564426D+00
    w(10) = 0.17656270536699264632527099011D+00
    w(11) = 0.16800410215645004450997066379D+00
    w(12) = 0.15404576107681028808143159480D+00
    w(13) = 0.13513636846852547328631998170D+00
    w(14) = 0.111883847193403971094788385626D+00
    w(15) = 0.085036148317179180883535370191D+00
    w(16) = 0.055459529373987201129440165359D+00
    w(17) = 0.024148302868547931960110026288D+00

  else if ( n == 18 ) then

    x(1) = -0.991565168420930946730016005D+00
    x(2) = -0.955823949571397755181195893D+00
    x(3) = -0.892602466497555739206060591D+00
    x(4) = -0.803704958972523115682417455D+00
    x(5) = -0.691687043060353207874891081D+00
    x(6) = -0.55977083107394753460787155D+00
    x(7) = -0.41175116146284264603593179D+00
    x(8) = -0.25188622569150550958897285D+00
    x(9) = -0.08477501304173530124226185D+00
    x(10) = 0.08477501304173530124226185D+00
    x(11) = 0.25188622569150550958897285D+00
    x(12) = 0.41175116146284264603593179D+00
    x(13) = 0.55977083107394753460787155D+00
    x(14) = 0.69168704306035320787489108D+00
    x(15) = 0.80370495897252311568241746D+00
    x(16) = 0.89260246649755573920606059D+00
    x(17) = 0.95582394957139775518119589D+00
    x(18) = 0.99156516842093094673001600D+00

    w(1) = 0.021616013526483310313342710266D+00
    w(2) = 0.049714548894969796453334946203D+00
    w(3) = 0.07642573025488905652912967762D+00
    w(4) = 0.10094204410628716556281398492D+00
    w(5) = 0.12255520671147846018451912680D+00
    w(6) = 0.14064291467065065120473130375D+00
    w(7) = 0.15468467512626524492541800384D+00
    w(8) = 0.16427648374583272298605377647D+00
    w(9) = 0.16914238296314359184065647013D+00
    w(10) = 0.16914238296314359184065647013D+00
    w(11) = 0.16427648374583272298605377647D+00
    w(12) = 0.15468467512626524492541800384D+00
    w(13) = 0.14064291467065065120473130375D+00
    w(14) = 0.12255520671147846018451912680D+00
    w(15) = 0.10094204410628716556281398492D+00
    w(16) = 0.07642573025488905652912967762D+00
    w(17) = 0.049714548894969796453334946203D+00
    w(18) = 0.021616013526483310313342710266D+00

  else if ( n == 19 ) then

    x(1) = -0.992406843843584403189017670D+00
    x(2) = -0.960208152134830030852778841D+00
    x(3) = -0.903155903614817901642660929D+00
    x(4) = -0.822714656537142824978922487D+00
    x(5) = -0.72096617733522937861709586D+00
    x(6) = -0.60054530466168102346963816D+00
    x(7) = -0.46457074137596094571726715D+00
    x(8) = -0.31656409996362983199011733D+00
    x(9) = -0.16035864564022537586809612D+00
    x(10) = 0.00000000000000000000000000D+00
    x(11) = 0.16035864564022537586809612D+00
    x(12) = 0.31656409996362983199011733D+00
    x(13) = 0.46457074137596094571726715D+00
    x(14) = 0.60054530466168102346963816D+00
    x(15) = 0.72096617733522937861709586D+00
    x(16) = 0.82271465653714282497892249D+00
    x(17) = 0.90315590361481790164266093D+00
    x(18) = 0.96020815213483003085277884D+00
    x(19) = 0.99240684384358440318901767D+00

    w(1) = 0.019461788229726477036312041464D+00
    w(2) = 0.044814226765699600332838157402D+00
    w(3) = 0.069044542737641226580708258006D+00
    w(4) = 0.091490021622449999464462094124D+00
    w(5) = 0.111566645547333994716023901682D+00
    w(6) = 0.12875396253933622767551578486D+00
    w(7) = 0.14260670217360661177574610944D+00
    w(8) = 0.15276604206585966677885540090D+00
    w(9) = 0.15896884339395434764995643946D+00
    w(10) = 0.16105444984878369597916362532D+00
    w(11) = 0.15896884339395434764995643946D+00
    w(12) = 0.15276604206585966677885540090D+00
    w(13) = 0.14260670217360661177574610944D+00
    w(14) = 0.12875396253933622767551578486D+00
    w(15) = 0.111566645547333994716023901682D+00
    w(16) = 0.091490021622449999464462094124D+00
    w(17) = 0.069044542737641226580708258006D+00
    w(18) = 0.044814226765699600332838157402D+00
    w(19) = 0.019461788229726477036312041464D+00

  else if ( n == 20 ) then

    x(1) = -0.993128599185094924786122388D+00
    x(2) = -0.963971927277913791267666131D+00
    x(3) = -0.912234428251325905867752441D+00
    x(4) = -0.83911697182221882339452906D+00
    x(5) = -0.74633190646015079261430507D+00
    x(6) = -0.63605368072651502545283670D+00
    x(7) = -0.51086700195082709800436405D+00
    x(8) = -0.37370608871541956067254818D+00
    x(9) = -0.22778585114164507808049620D+00
    x(10) = -0.07652652113349733375464041D+00
    x(11) = 0.07652652113349733375464041D+00
    x(12) = 0.22778585114164507808049620D+00
    x(13) = 0.37370608871541956067254818D+00
    x(14) = 0.51086700195082709800436405D+00
    x(15) = 0.63605368072651502545283670D+00
    x(16) = 0.74633190646015079261430507D+00
    x(17) = 0.83911697182221882339452906D+00
    x(18) = 0.91223442825132590586775244D+00
    x(19) = 0.96397192727791379126766613D+00
    x(20) = 0.99312859918509492478612239D+00

    w(1) = 0.017614007139152118311861962352D+00
    w(2) = 0.040601429800386941331039952275D+00
    w(3) = 0.062672048334109063569506535187D+00
    w(4) = 0.08327674157670474872475814322D+00
    w(5) = 0.10193011981724043503675013548D+00
    w(6) = 0.11819453196151841731237737771D+00
    w(7) = 0.13168863844917662689849449975D+00
    w(8) = 0.14209610931838205132929832507D+00
    w(9) = 0.14917298647260374678782873700D+00
    w(10) = 0.15275338713072585069808433195D+00
    w(11) = 0.15275338713072585069808433195D+00
    w(12) = 0.14917298647260374678782873700D+00
    w(13) = 0.14209610931838205132929832507D+00
    w(14) = 0.13168863844917662689849449975D+00
    w(15) = 0.11819453196151841731237737771D+00
    w(16) = 0.10193011981724043503675013548D+00
    w(17) = 0.08327674157670474872475814322D+00
    w(18) = 0.062672048334109063569506535187D+00
    w(19) = 0.040601429800386941331039952275D+00
    w(20) = 0.017614007139152118311861962352D+00

  else if ( n == 21 ) then

    x( 1) =  -0.99375217062038950026024204D+00
    x( 2) =  -0.96722683856630629431662221D+00
    x( 3) =  -0.92009933415040082879018713D+00
    x( 4) =  -0.85336336458331728364725064D+00
    x( 5) =  -0.76843996347567790861587785D+00
    x( 6) =  -0.66713880419741231930596667D+00
    x( 7) =  -0.55161883588721980705901880D+00
    x( 8) =  -0.42434212020743878357366889D+00
    x( 9) =  -0.28802131680240109660079252D+00
    x(10) =  -0.14556185416089509093703098D+00
    x(11) =   0.00000000000000000000000000D+00
    x(12) =  +0.14556185416089509093703098D+00
    x(13) =  +0.28802131680240109660079252D+00
    x(14) =  +0.42434212020743878357366889D+00
    x(15) =  +0.55161883588721980705901880D+00
    x(16) =  +0.66713880419741231930596667D+00
    x(17) =  +0.76843996347567790861587785D+00
    x(18) =  +0.85336336458331728364725064D+00
    x(19) =  +0.92009933415040082879018713D+00
    x(20) =  +0.96722683856630629431662221D+00
    x(21) =  +0.99375217062038950026024204D+00

    w( 1) =   0.016017228257774333324224616858D+00
    w( 2) =   0.036953789770852493799950668299D+00
    w( 3) =   0.057134425426857208283635826472D+00
    w( 4) =   0.076100113628379302017051653300D+00
    w( 5) =   0.093444423456033861553289741114D+00
    w( 6) =   0.108797299167148377663474578070D+00
    w( 7) =   0.12183141605372853419536717713D+00
    w( 8) =   0.13226893863333746178105257450D+00
    w( 9) =   0.13988739479107315472213342387D+00
    w(10) =   0.14452440398997005906382716655D+00
    w(11) =   0.14608113364969042719198514768D+00
    w(12) =   0.14452440398997005906382716655D+00
    w(13) =   0.13988739479107315472213342387D+00
    w(14) =   0.13226893863333746178105257450D+00
    w(15) =   0.12183141605372853419536717713D+00
    w(16) =   0.108797299167148377663474578070D+00
    w(17) =   0.093444423456033861553289741114D+00
    w(18) =   0.076100113628379302017051653300D+00
    w(19) =   0.057134425426857208283635826472D+00
    w(20) =   0.036953789770852493799950668299D+00
    w(21) =   0.016017228257774333324224616858D+00

  else if ( n == 22 ) then

    x(1) = -0.99429458548239929207303142D+00
    x(2) = -0.97006049783542872712395099D+00
    x(3) = -0.92695677218717400052069294D+00
    x(4) = -0.86581257772030013653642564D+00
    x(5) = -0.78781680597920816200427796D+00
    x(6) = -0.69448726318668278005068984D+00
    x(7) = -0.58764040350691159295887693D+00
    x(8) = -0.46935583798675702640633071D+00
    x(9) = -0.34193582089208422515814742D+00
    x(10) = -0.20786042668822128547884653D+00
    x(11) = -0.06973927331972222121384180D+00
    x(12) = 0.06973927331972222121384180D+00
    x(13) = 0.20786042668822128547884653D+00
    x(14) = 0.34193582089208422515814742D+00
    x(15) = 0.46935583798675702640633071D+00
    x(16) = 0.58764040350691159295887693D+00
    x(17) = 0.69448726318668278005068984D+00
    x(18) = 0.78781680597920816200427796D+00
    x(19) = 0.86581257772030013653642564D+00
    x(20) = 0.92695677218717400052069294D+00
    x(21) = 0.97006049783542872712395099D+00
    x(22) = 0.99429458548239929207303142D+00

    w(1) = 0.014627995298272200684991098047D+00
    w(2) = 0.033774901584814154793302246866D+00
    w(3) = 0.052293335152683285940312051273D+00
    w(4) = 0.06979646842452048809496141893D+00
    w(5) = 0.08594160621706772741444368137D+00
    w(6) = 0.10041414444288096493207883783D+00
    w(7) = 0.11293229608053921839340060742D+00
    w(8) = 0.12325237681051242428556098615D+00
    w(9) = 0.13117350478706237073296499253D+00
    w(10) = 0.13654149834601517135257383123D+00
    w(11) = 0.13925187285563199337541024834D+00
    w(12) = 0.13925187285563199337541024834D+00
    w(13) = 0.13654149834601517135257383123D+00
    w(14) = 0.13117350478706237073296499253D+00
    w(15) = 0.12325237681051242428556098615D+00
    w(16) = 0.11293229608053921839340060742D+00
    w(17) = 0.10041414444288096493207883783D+00
    w(18) = 0.08594160621706772741444368137D+00
    w(19) = 0.06979646842452048809496141893D+00
    w(20) = 0.052293335152683285940312051273D+00
    w(21) = 0.033774901584814154793302246866D+00
    w(22) = 0.014627995298272200684991098047D+00

  else if ( n == 23 ) then

    x(1) = -0.99476933499755212352392572D+00
    x(2) = -0.97254247121811523195602408D+00
    x(3) = -0.93297108682601610234919699D+00
    x(4) = -0.87675235827044166737815689D+00
    x(5) = -0.80488840161883989215111841D+00
    x(6) = -0.71866136313195019446162448D+00
    x(7) = -0.61960987576364615638509731D+00
    x(8) = -0.50950147784600754968979305D+00
    x(9) = -0.39030103803029083142148887D+00
    x(10) = -0.26413568097034493053386954D+00
    x(11) = -0.13325682429846611093174268D+00
    x(12) = 0.00000000000000000000000000D+00
    x(13) = 0.13325682429846611093174268D+00
    x(14) = 0.26413568097034493053386954D+00
    x(15) = 0.39030103803029083142148887D+00
    x(16) = 0.50950147784600754968979305D+00
    x(17) = 0.61960987576364615638509731D+00
    x(18) = 0.71866136313195019446162448D+00
    x(19) = 0.80488840161883989215111841D+00
    x(20) = 0.87675235827044166737815689D+00
    x(21) = 0.93297108682601610234919699D+00
    x(22) = 0.97254247121811523195602408D+00
    x(23) = 0.99476933499755212352392572D+00

    w(1) = 0.013411859487141772081309493459D+00
    w(2) = 0.030988005856979444310694219642D+00
    w(3) = 0.048037671731084668571641071632D+00
    w(4) = 0.064232421408525852127169615159D+00
    w(5) = 0.079281411776718954922892524742D+00
    w(6) = 0.092915766060035147477018617370D+00
    w(7) = 0.104892091464541410074086185015D+00
    w(8) = 0.11499664022241136494164351293D+00
    w(9) = 0.12304908430672953046757840067D+00
    w(10) = 0.12890572218808214997859533940D+00
    w(11) = 0.13246203940469661737164246470D+00
    w(12) = 0.13365457218610617535145711055D+00
    w(13) = 0.13246203940469661737164246470D+00
    w(14) = 0.12890572218808214997859533940D+00
    w(15) = 0.12304908430672953046757840067D+00
    w(16) = 0.11499664022241136494164351293D+00
    w(17) = 0.104892091464541410074086185015D+00
    w(18) = 0.092915766060035147477018617370D+00
    w(19) = 0.079281411776718954922892524742D+00
    w(20) = 0.064232421408525852127169615159D+00
    w(21) = 0.048037671731084668571641071632D+00
    w(22) = 0.030988005856979444310694219642D+00
    w(23) = 0.013411859487141772081309493459D+00

  else if ( n == 24 ) then

    x(1) = -0.99518721999702136017999741D+00
    x(2) = -0.97472855597130949819839199D+00
    x(3) = -0.93827455200273275852364900D+00
    x(4) = -0.88641552700440103421315434D+00
    x(5) = -0.82000198597390292195394987D+00
    x(6) = -0.74012419157855436424382810D+00
    x(7) = -0.64809365193697556925249579D+00
    x(8) = -0.54542147138883953565837562D+00
    x(9) = -0.43379350762604513848708423D+00
    x(10) = -0.31504267969616337438679329D+00
    x(11) = -0.19111886747361630915863982D+00
    x(12) = -0.06405689286260562608504308D+00
    x(13) = 0.06405689286260562608504308D+00
    x(14) = 0.19111886747361630915863982D+00
    x(15) = 0.31504267969616337438679329D+00
    x(16) = 0.43379350762604513848708423D+00
    x(17) = 0.54542147138883953565837562D+00
    x(18) = 0.64809365193697556925249579D+00
    x(19) = 0.74012419157855436424382810D+00
    x(20) = 0.82000198597390292195394987D+00
    x(21) = 0.88641552700440103421315434D+00
    x(22) = 0.93827455200273275852364900D+00
    x(23) = 0.97472855597130949819839199D+00
    x(24) = 0.99518721999702136017999741D+00

    w(1) = 0.012341229799987199546805667070D+00
    w(2) = 0.028531388628933663181307815952D+00
    w(3) = 0.044277438817419806168602748211D+00
    w(4) = 0.059298584915436780746367758500D+00
    w(5) = 0.07334648141108030573403361525D+00
    w(6) = 0.08619016153195327591718520298D+00
    w(7) = 0.09761865210411388826988066446D+00
    w(8) = 0.10744427011596563478257734245D+00
    w(9) = 0.11550566805372560135334448391D+00
    w(10) = 0.12167047292780339120446315348D+00
    w(11) = 0.12583745634682829612137538251D+00
    w(12) = 0.12793819534675215697405616522D+00
    w(13) = 0.12793819534675215697405616522D+00
    w(14) = 0.12583745634682829612137538251D+00
    w(15) = 0.12167047292780339120446315348D+00
    w(16) = 0.11550566805372560135334448391D+00
    w(17) = 0.10744427011596563478257734245D+00
    w(18) = 0.09761865210411388826988066446D+00
    w(19) = 0.08619016153195327591718520298D+00
    w(20) = 0.07334648141108030573403361525D+00
    w(21) = 0.059298584915436780746367758500D+00
    w(22) = 0.044277438817419806168602748211D+00
    w(23) = 0.028531388628933663181307815952D+00
    w(24) = 0.012341229799987199546805667070D+00

  else if ( n == 25 ) then

    x(1) = -0.99555696979049809790878495D+00
    x(2) = -0.97666392145951751149831539D+00
    x(3) = -0.94297457122897433941401117D+00
    x(4) = -0.89499199787827536885104201D+00
    x(5) = -0.83344262876083400142102111D+00
    x(6) = -0.75925926303735763057728287D+00
    x(7) = -0.67356636847346836448512063D+00
    x(8) = -0.57766293024122296772368984D+00
    x(9) = -0.47300273144571496052218212D+00
    x(10) = -0.36117230580938783773582173D+00
    x(11) = -0.24386688372098843204519036D+00
    x(12) = -0.12286469261071039638735982D+00
    x(13) = 0.00000000000000000000000000D+00
    x(14) = 0.12286469261071039638735982D+00
    x(15) = 0.24386688372098843204519036D+00
    x(16) = 0.36117230580938783773582173D+00
    x(17) = 0.47300273144571496052218212D+00
    x(18) = 0.57766293024122296772368984D+00
    x(19) = 0.67356636847346836448512063D+00
    x(20) = 0.75925926303735763057728287D+00
    x(21) = 0.83344262876083400142102111D+00
    x(22) = 0.89499199787827536885104201D+00
    x(23) = 0.94297457122897433941401117D+00
    x(24) = 0.97666392145951751149831539D+00
    x(25) = 0.99555696979049809790878495D+00

    w(1) = 0.0113937985010262879479029641132D+00
    w(2) = 0.026354986615032137261901815295D+00
    w(3) = 0.040939156701306312655623487712D+00
    w(4) = 0.054904695975835191925936891541D+00
    w(5) = 0.068038333812356917207187185657D+00
    w(6) = 0.080140700335001018013234959669D+00
    w(7) = 0.091028261982963649811497220703D+00
    w(8) = 0.100535949067050644202206890393D+00
    w(9) = 0.108519624474263653116093957050D+00
    w(10) = 0.11485825914571164833932554587D+00
    w(11) = 0.11945576353578477222817812651D+00
    w(12) = 0.12224244299031004168895951895D+00
    w(13) = 0.12317605372671545120390287308D+00
    w(14) = 0.12224244299031004168895951895D+00
    w(15) = 0.11945576353578477222817812651D+00
    w(16) = 0.11485825914571164833932554587D+00
    w(17) = 0.108519624474263653116093957050D+00
    w(18) = 0.100535949067050644202206890393D+00
    w(19) = 0.091028261982963649811497220703D+00
    w(20) = 0.080140700335001018013234959669D+00
    w(21) = 0.068038333812356917207187185657D+00
    w(22) = 0.054904695975835191925936891541D+00
    w(23) = 0.040939156701306312655623487712D+00
    w(24) = 0.026354986615032137261901815295D+00
    w(25) = 0.0113937985010262879479029641132D+00

  else if ( n == 26 ) then

    x(1) = -0.99588570114561692900321696D+00
    x(2) = -0.97838544595647099110058035D+00
    x(3) = -0.94715906666171425013591528D+00
    x(4) = -0.90263786198430707421766560D+00
    x(5) = -0.84544594278849801879750706D+00
    x(6) = -0.77638594882067885619296725D+00
    x(7) = -0.69642726041995726486381391D+00
    x(8) = -0.60669229301761806323197875D+00
    x(9) = -0.50844071482450571769570306D+00
    x(10) = -0.40305175512348630648107738D+00
    x(11) = -0.29200483948595689514283538D+00
    x(12) = -0.17685882035689018396905775D+00
    x(13) = -0.05923009342931320709371858D+00
    x(14) = 0.05923009342931320709371858D+00
    x(15) = 0.17685882035689018396905775D+00
    x(16) = 0.29200483948595689514283538D+00
    x(17) = 0.40305175512348630648107738D+00
    x(18) = 0.50844071482450571769570306D+00
    x(19) = 0.60669229301761806323197875D+00
    x(20) = 0.69642726041995726486381391D+00
    x(21) = 0.77638594882067885619296725D+00
    x(22) = 0.84544594278849801879750706D+00
    x(23) = 0.90263786198430707421766560D+00
    x(24) = 0.94715906666171425013591528D+00
    x(25) = 0.97838544595647099110058035D+00
    x(26) = 0.99588570114561692900321696D+00

    w(1) = 0.010551372617343007155651187685D+00
    w(2) = 0.024417851092631908789615827520D+00
    w(3) = 0.037962383294362763950303141249D+00
    w(4) = 0.050975825297147811998319900724D+00
    w(5) = 0.063274046329574835539453689907D+00
    w(6) = 0.07468414976565974588707579610D+00
    w(7) = 0.08504589431348523921044776508D+00
    w(8) = 0.09421380035591414846366488307D+00
    w(9) = 0.10205916109442542323841407025D+00
    w(10) = 0.10847184052857659065657942673D+00
    w(11) = 0.11336181654631966654944071844D+00
    w(12) = 0.11666044348529658204466250754D+00
    w(13) = 0.11832141527926227651637108570D+00
    w(14) = 0.11832141527926227651637108570D+00
    w(15) = 0.11666044348529658204466250754D+00
    w(16) = 0.11336181654631966654944071844D+00
    w(17) = 0.10847184052857659065657942673D+00
    w(18) = 0.10205916109442542323841407025D+00
    w(19) = 0.09421380035591414846366488307D+00
    w(20) = 0.08504589431348523921044776508D+00
    w(21) = 0.07468414976565974588707579610D+00
    w(22) = 0.063274046329574835539453689907D+00
    w(23) = 0.050975825297147811998319900724D+00
    w(24) = 0.037962383294362763950303141249D+00
    w(25) = 0.024417851092631908789615827520D+00
    w(26) = 0.010551372617343007155651187685D+00

  else if ( n == 27 ) then

    x(1) = -0.99617926288898856693888721D+00
    x(2) = -0.97992347596150122285587336D+00
    x(3) = -0.95090055781470500685190803D+00
    x(4) = -0.90948232067749110430064502D+00
    x(5) = -0.85620790801829449030273722D+00
    x(6) = -0.79177163907050822714439734D+00
    x(7) = -0.71701347373942369929481621D+00
    x(8) = -0.63290797194649514092773464D+00
    x(9) = -0.54055156457945689490030094D+00
    x(10) = -0.44114825175002688058597416D+00
    x(11) = -0.33599390363850889973031903D+00
    x(12) = -0.22645936543953685885723911D+00
    x(13) = -0.11397258560952996693289498D+00
    x(14) = 0.00000000000000000000000000D+00
    x(15) = 0.11397258560952996693289498D+00
    x(16) = 0.22645936543953685885723911D+00
    x(17) = 0.33599390363850889973031903D+00
    x(18) = 0.44114825175002688058597416D+00
    x(19) = 0.54055156457945689490030094D+00
    x(20) = 0.63290797194649514092773464D+00
    x(21) = 0.71701347373942369929481621D+00
    x(22) = 0.79177163907050822714439734D+00
    x(23) = 0.85620790801829449030273722D+00
    x(24) = 0.90948232067749110430064502D+00
    x(25) = 0.95090055781470500685190803D+00
    x(26) = 0.97992347596150122285587336D+00
    x(27) = 0.99617926288898856693888721D+00

    w(1) = 0.0097989960512943602611500550912D+00
    w(2) = 0.022686231596180623196034206447D+00
    w(3) = 0.035297053757419711022578289305D+00
    w(4) = 0.047449412520615062704096710114D+00
    w(5) = 0.058983536859833599110300833720D+00
    w(6) = 0.069748823766245592984322888357D+00
    w(7) = 0.079604867773057771263074959010D+00
    w(8) = 0.088423158543756950194322802854D+00
    w(9) = 0.096088727370028507565652646558D+00
    w(10) = 0.102501637817745798671247711533D+00
    w(11) = 0.107578285788533187212162984427D+00
    w(12) = 0.111252488356845192672163096043D+00
    w(13) = 0.113476346108965148620369948092D+00
    w(14) = 0.11422086737895698904504573690D+00
    w(15) = 0.113476346108965148620369948092D+00
    w(16) = 0.111252488356845192672163096043D+00
    w(17) = 0.107578285788533187212162984427D+00
    w(18) = 0.102501637817745798671247711533D+00
    w(19) = 0.096088727370028507565652646558D+00
    w(20) = 0.088423158543756950194322802854D+00
    w(21) = 0.079604867773057771263074959010D+00
    w(22) = 0.069748823766245592984322888357D+00
    w(23) = 0.058983536859833599110300833720D+00
    w(24) = 0.047449412520615062704096710114D+00
    w(25) = 0.035297053757419711022578289305D+00
    w(26) = 0.022686231596180623196034206447D+00
    w(27) = 0.0097989960512943602611500550912D+00

  else if ( n == 28 ) then

    x(1) = -0.99644249757395444995043639D+00
    x(2) = -0.98130316537087275369455995D+00
    x(3) = -0.95425928062893819725410184D+00
    x(4) = -0.91563302639213207386968942D+00
    x(5) = -0.86589252257439504894225457D+00
    x(6) = -0.80564137091717917144788596D+00
    x(7) = -0.73561087801363177202814451D+00
    x(8) = -0.65665109403886496121989818D+00
    x(9) = -0.56972047181140171930800328D+00
    x(10) = -0.47587422495511826103441185D+00
    x(11) = -0.37625151608907871022135721D+00
    x(12) = -0.27206162763517807767682636D+00
    x(13) = -0.16456928213338077128147178D+00
    x(14) = -0.05507928988403427042651653D+00
    x(15) = 0.05507928988403427042651653D+00
    x(16) = 0.16456928213338077128147178D+00
    x(17) = 0.27206162763517807767682636D+00
    x(18) = 0.37625151608907871022135721D+00
    x(19) = 0.47587422495511826103441185D+00
    x(20) = 0.56972047181140171930800328D+00
    x(21) = 0.65665109403886496121989818D+00
    x(22) = 0.73561087801363177202814451D+00
    x(23) = 0.80564137091717917144788596D+00
    x(24) = 0.86589252257439504894225457D+00
    x(25) = 0.91563302639213207386968942D+00
    x(26) = 0.95425928062893819725410184D+00
    x(27) = 0.98130316537087275369455995D+00
    x(28) = 0.99644249757395444995043639D+00

    w(1) = 0.009124282593094517738816153923D+00
    w(2) = 0.021132112592771259751500380993D+00
    w(3) = 0.032901427782304379977630819171D+00
    w(4) = 0.044272934759004227839587877653D+00
    w(5) = 0.055107345675716745431482918227D+00
    w(6) = 0.06527292396699959579339756678D+00
    w(7) = 0.07464621423456877902393188717D+00
    w(8) = 0.08311341722890121839039649824D+00
    w(9) = 0.09057174439303284094218603134D+00
    w(10) = 0.09693065799792991585048900610D+00
    w(11) = 0.10211296757806076981421663851D+00
    w(12) = 0.10605576592284641791041643700D+00
    w(13) = 0.10871119225829413525357151930D+00
    w(14) = 0.11004701301647519628237626560D+00
    w(15) = 0.11004701301647519628237626560D+00
    w(16) = 0.10871119225829413525357151930D+00
    w(17) = 0.10605576592284641791041643700D+00
    w(18) = 0.10211296757806076981421663851D+00
    w(19) = 0.09693065799792991585048900610D+00
    w(20) = 0.09057174439303284094218603134D+00
    w(21) = 0.08311341722890121839039649824D+00
    w(22) = 0.07464621423456877902393188717D+00
    w(23) = 0.06527292396699959579339756678D+00
    w(24) = 0.055107345675716745431482918227D+00
    w(25) = 0.044272934759004227839587877653D+00
    w(26) = 0.032901427782304379977630819171D+00
    w(27) = 0.021132112592771259751500380993D+00
    w(28) = 0.009124282593094517738816153923D+00

  else if ( n == 29 ) then

    x(1) = -0.99667944226059658616319153D+00
    x(2) = -0.98254550526141317487092602D+00
    x(3) = -0.95728559577808772579820804D+00
    x(4) = -0.92118023295305878509375344D+00
    x(5) = -0.87463780492010279041779342D+00
    x(6) = -0.81818548761525244498957221D+00
    x(7) = -0.75246285173447713391261008D+00
    x(8) = -0.67821453760268651515618501D+00
    x(9) = -0.59628179713822782037958621D+00
    x(10) = -0.50759295512422764210262792D+00
    x(11) = -0.41315288817400866389070659D+00
    x(12) = -0.31403163786763993494819592D+00
    x(13) = -0.21135228616600107450637573D+00
    x(14) = -0.10627823013267923017098239D+00
    x(15) = 0.00000000000000000000000000D+00
    x(16) = 0.10627823013267923017098239D+00
    x(17) = 0.21135228616600107450637573D+00
    x(18) = 0.31403163786763993494819592D+00
    x(19) = 0.41315288817400866389070659D+00
    x(20) = 0.50759295512422764210262792D+00
    x(21) = 0.59628179713822782037958621D+00
    x(22) = 0.67821453760268651515618501D+00
    x(23) = 0.75246285173447713391261008D+00
    x(24) = 0.81818548761525244498957221D+00
    x(25) = 0.87463780492010279041779342D+00
    x(26) = 0.92118023295305878509375344D+00
    x(27) = 0.95728559577808772579820804D+00
    x(28) = 0.98254550526141317487092602D+00
    x(29) = 0.99667944226059658616319153D+00

    w(1) = 0.0085169038787464096542638133022D+00
    w(2) = 0.019732085056122705983859801640D+00
    w(3) = 0.030740492202093622644408525375D+00
    w(4) = 0.041402062518682836104830010114D+00
    w(5) = 0.051594826902497923912594381180D+00
    w(6) = 0.061203090657079138542109848024D+00
    w(7) = 0.070117933255051278569581486949D+00
    w(8) = 0.078238327135763783828144888660D+00
    w(9) = 0.085472257366172527545344849297D+00
    w(10) = 0.091737757139258763347966411077D+00
    w(11) = 0.096963834094408606301900074883D+00
    w(12) = 0.101091273759914966121820546907D+00
    w(13) = 0.104073310077729373913328471285D+00
    w(14) = 0.105876155097320941406591327852D+00
    w(15) = 0.10647938171831424424651112691D+00
    w(16) = 0.105876155097320941406591327852D+00
    w(17) = 0.104073310077729373913328471285D+00
    w(18) = 0.101091273759914966121820546907D+00
    w(19) = 0.096963834094408606301900074883D+00
    w(20) = 0.091737757139258763347966411077D+00
    w(21) = 0.085472257366172527545344849297D+00
    w(22) = 0.078238327135763783828144888660D+00
    w(23) = 0.070117933255051278569581486949D+00
    w(24) = 0.061203090657079138542109848024D+00
    w(25) = 0.051594826902497923912594381180D+00
    w(26) = 0.041402062518682836104830010114D+00
    w(27) = 0.030740492202093622644408525375D+00
    w(28) = 0.019732085056122705983859801640D+00
    w(29) = 0.0085169038787464096542638133022D+00

  else if ( n == 30 ) then

    x(1) = -0.99689348407464954027163005D+00
    x(2) = -0.98366812327974720997003258D+00
    x(3) = -0.96002186496830751221687103D+00
    x(4) = -0.92620004742927432587932428D+00
    x(5) = -0.88256053579205268154311646D+00
    x(6) = -0.82956576238276839744289812D+00
    x(7) = -0.76777743210482619491797734D+00
    x(8) = -0.69785049479331579693229239D+00
    x(9) = -0.62052618298924286114047756D+00
    x(10) = -0.53662414814201989926416979D+00
    x(11) = -0.44703376953808917678060990D+00
    x(12) = -0.35270472553087811347103721D+00
    x(13) = -0.25463692616788984643980513D+00
    x(14) = -0.15386991360858354696379467D+00
    x(15) = -0.05147184255531769583302521D+00
    x(16) = 0.05147184255531769583302521D+00
    x(17) = 0.15386991360858354696379467D+00
    x(18) = 0.25463692616788984643980513D+00
    x(19) = 0.35270472553087811347103721D+00
    x(20) = 0.44703376953808917678060990D+00
    x(21) = 0.53662414814201989926416979D+00
    x(22) = 0.62052618298924286114047756D+00
    x(23) = 0.69785049479331579693229239D+00
    x(24) = 0.76777743210482619491797734D+00
    x(25) = 0.82956576238276839744289812D+00
    x(26) = 0.88256053579205268154311646D+00
    x(27) = 0.92620004742927432587932428D+00
    x(28) = 0.96002186496830751221687103D+00
    x(29) = 0.98366812327974720997003258D+00
    x(30) = 0.99689348407464954027163005D+00

    w(1) = 0.007968192496166605615465883475D+00
    w(2) = 0.018466468311090959142302131912D+00
    w(3) = 0.028784707883323369349719179611D+00
    w(4) = 0.038799192569627049596801936446D+00
    w(5) = 0.048402672830594052902938140423D+00
    w(6) = 0.057493156217619066481721689402D+00
    w(7) = 0.06597422988218049512812851512D+00
    w(8) = 0.07375597473770520626824385002D+00
    w(9) = 0.08075589522942021535469493846D+00
    w(10) = 0.08689978720108297980238753072D+00
    w(11) = 0.09212252223778612871763270709D+00
    w(12) = 0.09636873717464425963946862635D+00
    w(13) = 0.09959342058679526706278028210D+00
    w(14) = 0.10176238974840550459642895217D+00
    w(15) = 0.10285265289355884034128563671D+00
    w(16) = 0.10285265289355884034128563671D+00
    w(17) = 0.10176238974840550459642895217D+00
    w(18) = 0.09959342058679526706278028210D+00
    w(19) = 0.09636873717464425963946862635D+00
    w(20) = 0.09212252223778612871763270709D+00
    w(21) = 0.08689978720108297980238753072D+00
    w(22) = 0.08075589522942021535469493846D+00
    w(23) = 0.07375597473770520626824385002D+00
    w(24) = 0.06597422988218049512812851512D+00
    w(25) = 0.057493156217619066481721689402D+00
    w(26) = 0.048402672830594052902938140423D+00
    w(27) = 0.038799192569627049596801936446D+00
    w(28) = 0.028784707883323369349719179611D+00
    w(29) = 0.018466468311090959142302131912D+00
    w(30) = 0.007968192496166605615465883475D+00

  else if ( n == 31 ) then

    x(1) = -0.99708748181947707405562655D+00
    x(2) = -0.98468590966515248400246517D+00
    x(3) = -0.96250392509294966178905240D+00
    x(4) = -0.93075699789664816495694576D+00
    x(5) = -0.88976002994827104337419201D+00
    x(6) = -0.83992032014626734008690454D+00
    x(7) = -0.78173314841662494040636002D+00
    x(8) = -0.71577678458685328390597087D+00
    x(9) = -0.64270672292426034618441820D+00
    x(10) = -0.56324916140714926272094492D+00
    x(11) = -0.47819378204490248044059404D+00
    x(12) = -0.38838590160823294306135146D+00
    x(13) = -0.29471806998170161661790390D+00
    x(14) = -0.19812119933557062877241300D+00
    x(15) = -0.09955531215234152032517479D+00
    x(16) = 0.00000000000000000000000000D+00
    x(17) = 0.09955531215234152032517479D+00
    x(18) = 0.19812119933557062877241300D+00
    x(19) = 0.29471806998170161661790390D+00
    x(20) = 0.38838590160823294306135146D+00
    x(21) = 0.47819378204490248044059404D+00
    x(22) = 0.56324916140714926272094492D+00
    x(23) = 0.64270672292426034618441820D+00
    x(24) = 0.71577678458685328390597087D+00
    x(25) = 0.78173314841662494040636002D+00
    x(26) = 0.83992032014626734008690454D+00
    x(27) = 0.88976002994827104337419201D+00
    x(28) = 0.93075699789664816495694576D+00
    x(29) = 0.96250392509294966178905240D+00
    x(30) = 0.98468590966515248400246517D+00
    x(31) = 0.99708748181947707405562655D+00

    w(1) = 0.0074708315792487758586968750322D+00
    w(2) = 0.017318620790310582463157996087D+00
    w(3) = 0.027009019184979421800608708092D+00
    w(4) = 0.036432273912385464024392010468D+00
    w(5) = 0.045493707527201102902315857895D+00
    w(6) = 0.054103082424916853711666259087D+00
    w(7) = 0.062174786561028426910343543687D+00
    w(8) = 0.069628583235410366167756126255D+00
    w(9) = 0.076390386598776616426357674901D+00
    w(10) = 0.082392991761589263903823367432D+00
    w(11) = 0.087576740608477876126198069695D+00
    w(12) = 0.091890113893641478215362871607D+00
    w(13) = 0.095290242912319512807204197488D+00
    w(14) = 0.097743335386328725093474010979D+00
    w(15) = 0.099225011226672307874875514429D+00
    w(16) = 0.09972054479342645142753383373D+00
    w(17) = 0.099225011226672307874875514429D+00
    w(18) = 0.097743335386328725093474010979D+00
    w(19) = 0.095290242912319512807204197488D+00
    w(20) = 0.091890113893641478215362871607D+00
    w(21) = 0.087576740608477876126198069695D+00
    w(22) = 0.082392991761589263903823367432D+00
    w(23) = 0.076390386598776616426357674901D+00
    w(24) = 0.069628583235410366167756126255D+00
    w(25) = 0.062174786561028426910343543687D+00
    w(26) = 0.054103082424916853711666259087D+00
    w(27) = 0.045493707527201102902315857895D+00
    w(28) = 0.036432273912385464024392010468D+00
    w(29) = 0.027009019184979421800608708092D+00
    w(30) = 0.017318620790310582463157996087D+00
    w(31) = 0.0074708315792487758586968750322D+00

  else if ( n == 32 ) then

    x(1) = -0.99726386184948156354498113D+00
    x(2) = -0.98561151154526833540017504D+00
    x(3) = -0.96476225558750643077381193D+00
    x(4) = -0.93490607593773968917091913D+00
    x(5) = -0.89632115576605212396530724D+00
    x(6) = -0.84936761373256997013369300D+00
    x(7) = -0.79448379596794240696309730D+00
    x(8) = -0.73218211874028968038742667D+00
    x(9) = -0.66304426693021520097511517D+00
    x(10) = -0.58771575724076232904074548D+00
    x(11) = -0.50689990893222939002374747D+00
    x(12) = -0.42135127613063534536411944D+00
    x(13) = -0.33186860228212764977991681D+00
    x(14) = -0.23928736225213707454460321D+00
    x(15) = -0.14447196158279649348518637D+00
    x(16) = -0.04830766568773831623481257D+00
    x(17) = 0.04830766568773831623481257D+00
    x(18) = 0.14447196158279649348518637D+00
    x(19) = 0.23928736225213707454460321D+00
    x(20) = 0.33186860228212764977991681D+00
    x(21) = 0.42135127613063534536411944D+00
    x(22) = 0.50689990893222939002374747D+00
    x(23) = 0.58771575724076232904074548D+00
    x(24) = 0.66304426693021520097511517D+00
    x(25) = 0.73218211874028968038742667D+00
    x(26) = 0.79448379596794240696309730D+00
    x(27) = 0.84936761373256997013369300D+00
    x(28) = 0.89632115576605212396530724D+00
    x(29) = 0.93490607593773968917091913D+00
    x(30) = 0.96476225558750643077381193D+00
    x(31) = 0.98561151154526833540017504D+00
    x(32) = 0.99726386184948156354498113D+00

    w(1) = 0.007018610009470096600407063739D+00
    w(2) = 0.016274394730905670605170562206D+00
    w(3) = 0.025392065309262059455752589789D+00
    w(4) = 0.034273862913021433102687732252D+00
    w(5) = 0.042835898022226680656878646606D+00
    w(6) = 0.050998059262376176196163244690D+00
    w(7) = 0.058684093478535547145283637300D+00
    w(8) = 0.06582222277636184683765006371D+00
    w(9) = 0.07234579410884850622539935648D+00
    w(10) = 0.07819389578707030647174091883D+00
    w(11) = 0.08331192422694675522219907460D+00
    w(12) = 0.08765209300440381114277146275D+00
    w(13) = 0.09117387869576388471286857711D+00
    w(14) = 0.09384439908080456563918023767D+00
    w(15) = 0.09563872007927485941908200220D+00
    w(16) = 0.09654008851472780056676483006D+00
    w(17) = 0.09654008851472780056676483006D+00
    w(18) = 0.09563872007927485941908200220D+00
    w(19) = 0.09384439908080456563918023767D+00
    w(20) = 0.09117387869576388471286857711D+00
    w(21) = 0.08765209300440381114277146275D+00
    w(22) = 0.08331192422694675522219907460D+00
    w(23) = 0.07819389578707030647174091883D+00
    w(24) = 0.07234579410884850622539935648D+00
    w(25) = 0.06582222277636184683765006371D+00
    w(26) = 0.058684093478535547145283637300D+00
    w(27) = 0.050998059262376176196163244690D+00
    w(28) = 0.042835898022226680656878646606D+00
    w(29) = 0.034273862913021433102687732252D+00
    w(30) = 0.025392065309262059455752589789D+00
    w(31) = 0.016274394730905670605170562206D+00
    w(32) = 0.007018610009470096600407063739D+00

  else if ( n == 33 ) then

    x(1) = -0.99742469424645521726616802D+00
    x(2) = -0.98645572623064248811037570D+00
    x(3) = -0.96682290968999276892837771D+00
    x(4) = -0.93869437261116835035583512D+00
    x(5) = -0.90231676774343358304053133D+00
    x(6) = -0.85800965267650406464306148D+00
    x(7) = -0.80616235627416658979620087D+00
    x(8) = -0.74723049644956215785905512D+00
    x(9) = -0.68173195996974278626821595D+00
    x(10) = -0.61024234583637902730728751D+00
    x(11) = -0.53338990478634764354889426D+00
    x(12) = -0.45185001727245069572599328D+00
    x(13) = -0.36633925774807334107022062D+00
    x(14) = -0.27760909715249702940324807D+00
    x(15) = -0.18643929882799157233579876D+00
    x(16) = -0.09363106585473338567074292D+00
    x(17) = 0.00000000000000000000000000D+00
    x(18) = 0.09363106585473338567074292D+00
    x(19) = 0.18643929882799157233579876D+00
    x(20) = 0.27760909715249702940324807D+00
    x(21) = 0.36633925774807334107022062D+00
    x(22) = 0.45185001727245069572599328D+00
    x(23) = 0.53338990478634764354889426D+00
    x(24) = 0.61024234583637902730728751D+00
    x(25) = 0.68173195996974278626821595D+00
    x(26) = 0.74723049644956215785905512D+00
    x(27) = 0.80616235627416658979620087D+00
    x(28) = 0.85800965267650406464306148D+00
    x(29) = 0.90231676774343358304053133D+00
    x(30) = 0.93869437261116835035583512D+00
    x(31) = 0.96682290968999276892837771D+00
    x(32) = 0.98645572623064248811037570D+00
    x(33) = 0.99742469424645521726616802D+00

    w(1) = 0.0066062278475873780586492352085D+00
    w(2) = 0.015321701512934676127945768534D+00
    w(3) = 0.023915548101749480350533257529D+00
    w(4) = 0.032300358632328953281561447250D+00
    w(5) = 0.040401541331669591563409790527D+00
    w(6) = 0.048147742818711695670146880138D+00
    w(7) = 0.055470846631663561284944495439D+00
    w(8) = 0.062306482530317480031627725771D+00
    w(9) = 0.068594572818656712805955073015D+00
    w(10) = 0.074279854843954149342472175919D+00
    w(11) = 0.079312364794886738363908384942D+00
    w(12) = 0.083647876067038707613928014518D+00
    w(13) = 0.087248287618844337607281670945D+00
    w(14) = 0.090081958660638577239743705500D+00
    w(15) = 0.092123986643316846213240977717D+00
    w(16) = 0.093356426065596116160999126274D+00
    w(17) = 0.09376844616020999656730454155D+00
    w(18) = 0.093356426065596116160999126274D+00
    w(19) = 0.092123986643316846213240977717D+00
    w(20) = 0.090081958660638577239743705500D+00
    w(21) = 0.087248287618844337607281670945D+00
    w(22) = 0.083647876067038707613928014518D+00
    w(23) = 0.079312364794886738363908384942D+00
    w(24) = 0.074279854843954149342472175919D+00
    w(25) = 0.068594572818656712805955073015D+00
    w(26) = 0.062306482530317480031627725771D+00
    w(27) = 0.055470846631663561284944495439D+00
    w(28) = 0.048147742818711695670146880138D+00
    w(29) = 0.040401541331669591563409790527D+00
    w(30) = 0.032300358632328953281561447250D+00
    w(31) = 0.023915548101749480350533257529D+00
    w(32) = 0.015321701512934676127945768534D+00
    w(33) = 0.0066062278475873780586492352085D+00

  else if ( n == 63 ) then

    x(1) = -0.99928298402912378037893614D+00
    x(2) = -0.99622401277797010860219336D+00
    x(3) = -0.99072854689218946681089467D+00
    x(4) = -0.98280881059372723486251141D+00
    x(5) = -0.97248403469757002280196068D+00
    x(6) = -0.95977944975894192707035417D+00
    x(7) = -0.94472613404100980296637532D+00
    x(8) = -0.92736092062184320544703138D+00
    x(9) = -0.90772630277853155803695313D+00
    x(10) = -0.88587032850785342629029846D+00
    x(11) = -0.86184648236412371953961184D+00
    x(12) = -0.83571355431950284347180777D+00
    x(13) = -0.80753549577345676005146599D+00
    x(14) = -0.7773812629903723355633302D+00
    x(15) = -0.7453246483178474178293217D+00
    x(16) = -0.7114440995848458078514315D+00
    x(17) = -0.6758225281149860901311033D+00
    x(18) = -0.6385471058213653850003070D+00
    x(19) = -0.5997090518776252357390089D+00
    x(20) = -0.5594034094862850132676978D+00
    x(21) = -0.5177288132900332481244776D+00
    x(22) = -0.4747872479948043999222123D+00
    x(23) = -0.4306837987951116006620889D+00
    x(24) = -0.3855263942122478924776150D+00
    x(25) = -0.3394255419745844024688344D+00
    x(26) = -0.2924940585862514400361572D+00
    x(27) = -0.2448467932459533627484046D+00
    x(28) = -0.1966003467915066845576275D+00
    x(29) = -0.1478727863578719685698391D+00
    x(30) = -0.0987833564469452795297037D+00
    x(31) = -0.0494521871161596272342338D+00
    x(32) = 0.0000000000000000000000000D+00
    x(33) = 0.0494521871161596272342338D+00
    x(34) = 0.0987833564469452795297037D+00
    x(35) = 0.1478727863578719685698391D+00
    x(36) = 0.1966003467915066845576275D+00
    x(37) = 0.2448467932459533627484046D+00
    x(38) = 0.2924940585862514400361572D+00
    x(39) = 0.3394255419745844024688344D+00
    x(40) = 0.3855263942122478924776150D+00
    x(41) = 0.4306837987951116006620889D+00
    x(42) = 0.4747872479948043999222123D+00
    x(43) = 0.5177288132900332481244776D+00
    x(44) = 0.5594034094862850132676978D+00
    x(45) = 0.5997090518776252357390089D+00
    x(46) = 0.6385471058213653850003070D+00
    x(47) = 0.6758225281149860901311033D+00
    x(48) = 0.7114440995848458078514315D+00
    x(49) = 0.7453246483178474178293217D+00
    x(50) = 0.7773812629903723355633302D+00
    x(51) = 0.8075354957734567600514660D+00
    x(52) = 0.8357135543195028434718078D+00
    x(53) = 0.8618464823641237195396118D+00
    x(54) = 0.8858703285078534262902985D+00
    x(55) = 0.9077263027785315580369531D+00
    x(56) = 0.9273609206218432054470314D+00
    x(57) = 0.9447261340410098029663753D+00
    x(58) = 0.9597794497589419270703542D+00
    x(59) = 0.9724840346975700228019607D+00
    x(60) = 0.9828088105937272348625114D+00
    x(61) = 0.9907285468921894668108947D+00
    x(62) = 0.9962240127779701086021934D+00
    x(63) = 0.9992829840291237803789361D+00

    w(1) = 0.0018398745955770841170924455540D+00
    w(2) = 0.0042785083468637618660784110826D+00
    w(3) = 0.0067102917659601362519069307298D+00
    w(4) = 0.0091259686763266563540586454218D+00
    w(5) = 0.011519376076880041750750606149D+00
    w(6) = 0.013884612616115610824866086368D+00
    w(7) = 0.016215878410338338882283672975D+00
    w(8) = 0.018507464460161270409260545805D+00
    w(9) = 0.020753761258039090775341953421D+00
    w(10) = 0.022949271004889933148942319562D+00
    w(11) = 0.025088620553344986618630138068D+00
    w(12) = 0.027166574359097933225189839439D+00
    w(13) = 0.029178047208280526945551502154D+00
    w(14) = 0.031118116622219817508215988557D+00
    w(15) = 0.032982034883779341765683179672D+00
    w(16) = 0.034765240645355877697180504643D+00
    w(17) = 0.036463370085457289630452409788D+00
    w(18) = 0.038072267584349556763638324928D+00
    w(19) = 0.039587995891544093984807928149D+00
    w(20) = 0.041006845759666398635110037009D+00
    w(21) = 0.042325345020815822982505485403D+00
    w(22) = 0.043540267083027590798964315704D+00
    w(23) = 0.044648638825941395370332669517D+00
    w(24) = 0.045647747876292608685885992609D+00
    w(25) = 0.046535149245383696510395418747D+00
    w(26) = 0.047308671312268919080604988339D+00
    w(27) = 0.047966421137995131411052756195D+00
    w(28) = 0.048506789097883847864090099146D+00
    w(29) = 0.048928452820511989944709361549D+00
    w(30) = 0.049230380423747560785043116988D+00
    w(31) = 0.049411833039918178967039646117D+00
    w(32) = 0.04947236662393102088866936042D+00
    w(33) = 0.049411833039918178967039646117D+00
    w(34) = 0.049230380423747560785043116988D+00
    w(35) = 0.048928452820511989944709361549D+00
    w(36) = 0.048506789097883847864090099146D+00
    w(37) = 0.047966421137995131411052756195D+00
    w(38) = 0.047308671312268919080604988339D+00
    w(39) = 0.046535149245383696510395418747D+00
    w(40) = 0.045647747876292608685885992609D+00
    w(41) = 0.044648638825941395370332669517D+00
    w(42) = 0.043540267083027590798964315704D+00
    w(43) = 0.042325345020815822982505485403D+00
    w(44) = 0.041006845759666398635110037009D+00
    w(45) = 0.039587995891544093984807928149D+00
    w(46) = 0.038072267584349556763638324928D+00
    w(47) = 0.036463370085457289630452409788D+00
    w(48) = 0.034765240645355877697180504643D+00
    w(49) = 0.032982034883779341765683179672D+00
    w(50) = 0.031118116622219817508215988557D+00
    w(51) = 0.029178047208280526945551502154D+00
    w(52) = 0.027166574359097933225189839439D+00
    w(53) = 0.025088620553344986618630138068D+00
    w(54) = 0.022949271004889933148942319562D+00
    w(55) = 0.020753761258039090775341953421D+00
    w(56) = 0.018507464460161270409260545805D+00
    w(57) = 0.016215878410338338882283672975D+00
    w(58) = 0.013884612616115610824866086368D+00
    w(59) = 0.011519376076880041750750606149D+00
    w(60) = 0.0091259686763266563540586454218D+00
    w(61) = 0.0067102917659601362519069307298D+00
    w(62) = 0.0042785083468637618660784110826D+00
    w(63) = 0.0018398745955770841170924455540D+00

  else if ( n == 64 ) then

    x(1) = -0.99930504173577213945690562D+00
    x(2) = -0.99634011677195527934692450D+00
    x(3) = -0.99101337147674432073938238D+00
    x(4) = -0.98333625388462595693129930D+00
    x(5) = -0.97332682778991096374185351D+00
    x(6) = -0.96100879965205371891861412D+00
    x(7) = -0.94641137485840281606248149D+00
    x(8) = -0.92956917213193957582149015D+00
    x(9) = -0.91052213707850280575638067D+00
    x(10) = -0.88931544599511410585340404D+00
    x(11) = -0.86599939815409281976078339D+00
    x(12) = -0.8406292962525803627516915D+00
    x(13) = -0.8132653151227975597419233D+00
    x(14) = -0.7839723589433414076102205D+00
    x(15) = -0.7528199072605318966118638D+00
    x(16) = -0.7198818501716108268489402D+00
    x(17) = -0.6852363130542332425635584D+00
    x(18) = -0.6489654712546573398577612D+00
    x(19) = -0.6111553551723932502488530D+00
    x(20) = -0.5718956462026340342838781D+00
    x(21) = -0.5312794640198945456580139D+00
    x(22) = -0.4894031457070529574785263D+00
    x(23) = -0.4463660172534640879849477D+00
    x(24) = -0.4022701579639916036957668D+00
    x(25) = -0.3572201583376681159504426D+00
    x(26) = -0.3113228719902109561575127D+00
    x(27) = -0.2646871622087674163739642D+00
    x(28) = -0.2174236437400070841496487D+00
    x(29) = -0.1696444204239928180373136D+00
    x(30) = -0.1214628192961205544703765D+00
    x(31) = -0.0729931217877990394495429D+00
    x(32) = -0.0243502926634244325089558D+00
    x(33) = 0.0243502926634244325089558D+00
    x(34) = 0.0729931217877990394495429D+00
    x(35) = 0.1214628192961205544703765D+00
    x(36) = 0.1696444204239928180373136D+00
    x(37) = 0.2174236437400070841496487D+00
    x(38) = 0.2646871622087674163739642D+00
    x(39) = 0.3113228719902109561575127D+00
    x(40) = 0.3572201583376681159504426D+00
    x(41) = 0.4022701579639916036957668D+00
    x(42) = 0.4463660172534640879849477D+00
    x(43) = 0.4894031457070529574785263D+00
    x(44) = 0.5312794640198945456580139D+00
    x(45) = 0.5718956462026340342838781D+00
    x(46) = 0.6111553551723932502488530D+00
    x(47) = 0.6489654712546573398577612D+00
    x(48) = 0.6852363130542332425635584D+00
    x(49) = 0.7198818501716108268489402D+00
    x(50) = 0.7528199072605318966118638D+00
    x(51) = 0.7839723589433414076102205D+00
    x(52) = 0.8132653151227975597419233D+00
    x(53) = 0.8406292962525803627516915D+00
    x(54) = 0.8659993981540928197607834D+00
    x(55) = 0.8893154459951141058534040D+00
    x(56) = 0.9105221370785028057563807D+00
    x(57) = 0.9295691721319395758214902D+00
    x(58) = 0.9464113748584028160624815D+00
    x(59) = 0.9610087996520537189186141D+00
    x(60) = 0.9733268277899109637418535D+00
    x(61) = 0.9833362538846259569312993D+00
    x(62) = 0.9910133714767443207393824D+00
    x(63) = 0.9963401167719552793469245D+00
    x(64) = 0.9993050417357721394569056D+00

    w(1) = 0.0017832807216964329472960791450D+00
    w(2) = 0.0041470332605624676352875357286D+00
    w(3) = 0.006504457968978362856117360400D+00
    w(4) = 0.008846759826363947723030914660D+00
    w(5) = 0.011168139460131128818590493019D+00
    w(6) = 0.013463047896718642598060766686D+00
    w(7) = 0.015726030476024719321965995298D+00
    w(8) = 0.017951715775697343085045302001D+00
    w(9) = 0.020134823153530209372340316729D+00
    w(10) = 0.022270173808383254159298330384D+00
    w(11) = 0.024352702568710873338177550409D+00
    w(12) = 0.026377469715054658671691792625D+00
    w(13) = 0.028339672614259483227511305200D+00
    w(14) = 0.030234657072402478867974059820D+00
    w(15) = 0.032057928354851553585467504348D+00
    w(16) = 0.033805161837141609391565482111D+00
    w(17) = 0.035472213256882383810693146715D+00
    w(18) = 0.037055128540240046040415101810D+00
    w(19) = 0.038550153178615629128962496947D+00
    w(20) = 0.039953741132720341386656926128D+00
    w(21) = 0.041262563242623528610156297474D+00
    w(22) = 0.042473515123653589007339767909D+00
    w(23) = 0.043583724529323453376827860974D+00
    w(24) = 0.044590558163756563060134710031D+00
    w(25) = 0.045491627927418144479770996971D+00
    w(26) = 0.046284796581314417295953249232D+00
    w(27) = 0.046968182816210017325326285755D+00
    w(28) = 0.047540165714830308662282206944D+00
    w(29) = 0.04799938859645830772812617987D+00
    w(30) = 0.04834476223480295716976952716D+00
    w(31) = 0.04857546744150342693479906678D+00
    w(32) = 0.04869095700913972038336539073D+00
    w(33) = 0.04869095700913972038336539073D+00
    w(34) = 0.04857546744150342693479906678D+00
    w(35) = 0.04834476223480295716976952716D+00
    w(36) = 0.04799938859645830772812617987D+00
    w(37) = 0.047540165714830308662282206944D+00
    w(38) = 0.046968182816210017325326285755D+00
    w(39) = 0.046284796581314417295953249232D+00
    w(40) = 0.045491627927418144479770996971D+00
    w(41) = 0.044590558163756563060134710031D+00
    w(42) = 0.043583724529323453376827860974D+00
    w(43) = 0.042473515123653589007339767909D+00
    w(44) = 0.041262563242623528610156297474D+00
    w(45) = 0.039953741132720341386656926128D+00
    w(46) = 0.038550153178615629128962496947D+00
    w(47) = 0.037055128540240046040415101810D+00
    w(48) = 0.035472213256882383810693146715D+00
    w(49) = 0.033805161837141609391565482111D+00
    w(50) = 0.032057928354851553585467504348D+00
    w(51) = 0.030234657072402478867974059820D+00
    w(52) = 0.028339672614259483227511305200D+00
    w(53) = 0.026377469715054658671691792625D+00
    w(54) = 0.024352702568710873338177550409D+00
    w(55) = 0.022270173808383254159298330384D+00
    w(56) = 0.020134823153530209372340316729D+00
    w(57) = 0.017951715775697343085045302001D+00
    w(58) = 0.015726030476024719321965995298D+00
    w(59) = 0.013463047896718642598060766686D+00
    w(60) = 0.011168139460131128818590493019D+00
    w(61) = 0.008846759826363947723030914660D+00
    w(62) = 0.006504457968978362856117360400D+00
    w(63) = 0.0041470332605624676352875357286D+00
    w(64) = 0.0017832807216964329472960791450D+00

  else if ( n == 65 ) then

    x(1) = -0.99932609707541287726569361D+00
    x(2) = -0.99645094806184916305579494D+00
    x(3) = -0.99128527617680166872182118D+00
    x(4) = -0.98383981218703494137763778D+00
    x(5) = -0.97413153983355116907496789D+00
    x(6) = -0.96218275471805523771198375D+00
    x(7) = -0.94802092816840750637376974D+00
    x(8) = -0.93167862822874933796567699D+00
    x(9) = -0.91319344054284626173654692D+00
    x(10) = -0.89260788050473893142328554D+00
    x(11) = -0.8699692949264070361941320D+00
    x(12) = -0.8453297528999302839424500D+00
    x(13) = -0.8187459259226514534339191D+00
    x(14) = -0.7902789574921218430473804D+00
    x(15) = -0.7599943224419997868739828D+00
    x(16) = -0.7279616763294246790119737D+00
    x(17) = -0.6942546952139916335526225D+00
    x(18) = -0.6589509061936251330409408D+00
    x(19) = -0.6221315090854002415825996D+00
    x(20) = -0.5838811896604873133271545D+00
    x(21) = -0.5442879248622271385455725D+00
    x(22) = -0.5034427804550068823410431D+00
    x(23) = -0.4614397015691450576978341D+00
    x(24) = -0.4183752966234090092641990D+00
    x(25) = -0.3743486151220660120087939D+00
    x(26) = -0.3294609198374864076452867D+00
    x(27) = -0.2838154539022487306176554D+00
    x(28) = -0.2375172033464168065707124D+00
    x(29) = -0.1906726556261427697749124D+00
    x(30) = -0.1433895546989751711312496D+00
    x(31) = -0.0957766532091975056522186D+00
    x(32) = -0.0479434623531718575225298D+00
    x(33) = 0.0000000000000000000000000D+00
    x(34) = 0.0479434623531718575225298D+00
    x(35) = 0.0957766532091975056522186D+00
    x(36) = 0.1433895546989751711312496D+00
    x(37) = 0.1906726556261427697749124D+00
    x(38) = 0.2375172033464168065707124D+00
    x(39) = 0.2838154539022487306176554D+00
    x(40) = 0.3294609198374864076452867D+00
    x(41) = 0.3743486151220660120087939D+00
    x(42) = 0.4183752966234090092641990D+00
    x(43) = 0.4614397015691450576978341D+00
    x(44) = 0.5034427804550068823410431D+00
    x(45) = 0.5442879248622271385455725D+00
    x(46) = 0.5838811896604873133271545D+00
    x(47) = 0.6221315090854002415825996D+00
    x(48) = 0.6589509061936251330409408D+00
    x(49) = 0.6942546952139916335526225D+00
    x(50) = 0.7279616763294246790119737D+00
    x(51) = 0.7599943224419997868739828D+00
    x(52) = 0.7902789574921218430473804D+00
    x(53) = 0.8187459259226514534339191D+00
    x(54) = 0.8453297528999302839424500D+00
    x(55) = 0.8699692949264070361941320D+00
    x(56) = 0.8926078805047389314232855D+00
    x(57) = 0.9131934405428462617365469D+00
    x(58) = 0.9316786282287493379656770D+00
    x(59) = 0.9480209281684075063737697D+00
    x(60) = 0.9621827547180552377119837D+00
    x(61) = 0.9741315398335511690749679D+00
    x(62) = 0.9838398121870349413776378D+00
    x(63) = 0.9912852761768016687218212D+00
    x(64) = 0.9964509480618491630557949D+00
    x(65) = 0.9993260970754128772656936D+00

    w(1) = 0.0017292582513002508983395851463D+00
    w(2) = 0.0040215241720037363470786599528D+00
    w(3) = 0.0063079425789717545501888719039D+00
    w(4) = 0.0085801482668814598936358121592D+00
    w(5) = 0.0108326787895979686215140551272D+00
    w(6) = 0.013060311639994846336168342922D+00
    w(7) = 0.015257912146448310349265388145D+00
    w(8) = 0.017420421997670248495365759969D+00
    w(9) = 0.019542865836750062826837429313D+00
    w(10) = 0.021620361284934062841654274667D+00
    w(11) = 0.023648129691287236698780978994D+00
    w(12) = 0.025621506938037758214084978694D+00
    w(13) = 0.027535954088450343942499722327D+00
    w(14) = 0.029387067789310668062644859210D+00
    w(15) = 0.031170590380189142464431845777D+00
    w(16) = 0.032882419676368574984049638008D+00
    w(17) = 0.034518618398549058625221276859D+00
    w(18) = 0.036075423225565273932166270524D+00
    w(19) = 0.037549253448257709809772223198D+00
    w(20) = 0.038936719204051197616673806364D+00
    w(21) = 0.040234629273005533815446337743D+00
    w(22) = 0.041439998417240293022686299233D+00
    w(23) = 0.042550054246755802719217150803D+00
    w(24) = 0.043562243595800486532284821661D+00
    w(25) = 0.044474238395082974427323504000D+00
    w(26) = 0.045283941026300230657128240574D+00
    w(27) = 0.045989489146651696963893390818D+00
    w(28) = 0.046589259972233498302255136790D+00
    w(29) = 0.047081874010454522246006808290D+00
    w(30) = 0.047466198232885503152644458740D+00
    w(31) = 0.047741348681240621559038972227D+00
    w(32) = 0.047906692500495862031347289176D+00
    w(33) = 0.04796184939446661812070762137D+00
    w(34) = 0.047906692500495862031347289176D+00
    w(35) = 0.047741348681240621559038972227D+00
    w(36) = 0.047466198232885503152644458740D+00
    w(37) = 0.047081874010454522246006808290D+00
    w(38) = 0.046589259972233498302255136790D+00
    w(39) = 0.045989489146651696963893390818D+00
    w(40) = 0.045283941026300230657128240574D+00
    w(41) = 0.044474238395082974427323504000D+00
    w(42) = 0.043562243595800486532284821661D+00
    w(43) = 0.042550054246755802719217150803D+00
    w(44) = 0.041439998417240293022686299233D+00
    w(45) = 0.040234629273005533815446337743D+00
    w(46) = 0.038936719204051197616673806364D+00
    w(47) = 0.037549253448257709809772223198D+00
    w(48) = 0.036075423225565273932166270524D+00
    w(49) = 0.034518618398549058625221276859D+00
    w(50) = 0.032882419676368574984049638008D+00
    w(51) = 0.031170590380189142464431845777D+00
    w(52) = 0.029387067789310668062644859210D+00
    w(53) = 0.027535954088450343942499722327D+00
    w(54) = 0.025621506938037758214084978694D+00
    w(55) = 0.023648129691287236698780978994D+00
    w(56) = 0.021620361284934062841654274667D+00
    w(57) = 0.019542865836750062826837429313D+00
    w(58) = 0.017420421997670248495365759969D+00
    w(59) = 0.015257912146448310349265388145D+00
    w(60) = 0.013060311639994846336168342922D+00
    w(61) = 0.0108326787895979686215140551272D+00
    w(62) = 0.0085801482668814598936358121592D+00
    w(63) = 0.0063079425789717545501888719039D+00
    w(64) = 0.0040215241720037363470786599528D+00
    w(65) = 0.0017292582513002508983395851463D+00

  else if ( n == 127 ) then

    x(1) = -0.9998221304153061462673512D+00
    x(2) = -0.9990629343553118951383159D+00
    x(3) = -0.9976975661898046210744170D+00
    x(4) = -0.9957265513520272266354334D+00
    x(5) = -0.9931510492545171473611308D+00
    x(6) = -0.9899726145914841576077867D+00
    x(7) = -0.9861931740169316667104383D+00
    x(8) = -0.9818150208038141100334631D+00
    x(9) = -0.9768408123430703268174439D+00
    x(10) = -0.9712735681615291922889469D+00
    x(11) = -0.9651166679452921210908251D+00
    x(12) = -0.9583738494252387711491029D+00
    x(13) = -0.9510492060778803105479076D+00
    x(14) = -0.9431471846248148273454496D+00
    x(15) = -0.9346725823247379685736349D+00
    x(16) = -0.9256305440562338491274647D+00
    x(17) = -0.9160265591914658093130886D+00
    x(18) = -0.9058664582618213828024613D+00
    x(19) = -0.8951564094170837089690438D+00
    x(20) = -0.8839029146800265699452579D+00
    x(21) = -0.8721128059985607114196375D+00
    x(22) = -0.8597932410977408098120313D+00
    x(23) = -0.8469516991340975984533393D+00
    x(24) = -0.8335959761548995143795572D+00
    x(25) = -0.8197341803650786741551191D+00
    x(26) = -0.8053747272046802146665608D+00
    x(27) = -0.7905263342398137999454500D+00
    x(28) = -0.7751980158702023824449628D+00
    x(29) = -0.7593990778565366715566637D+00
    x(30) = -0.7431391116709545129205669D+00
    x(31) = -0.7264279886740726855356929D+00
    x(32) = -0.7092758541221045609994446D+00
    x(33) = -0.6916931210077006701564414D+00
    x(34) = -0.6736904637382504853466825D+00
    x(35) = -0.6552788116554826302767651D+00
    x(36) = -0.6364693424002972413476082D+00
    x(37) = -0.6172734751268582838576392D+00
    x(38) = -0.5977028635700652293844120D+00
    x(39) = -0.5777693889706125800032517D+00
    x(40) = -0.5574851528619322329218619D+00
    x(41) = -0.5368624697233975674581664D+00
    x(42) = -0.5159138595042493572772773D+00
    x(43) = -0.4946520400227821173949402D+00
    x(44) = -0.4730899192454052416450999D+00
    x(45) = -0.4512405874502662273318986D+00
    x(46) = -0.4291173092801933762625441D+00
    x(47) = -0.4067335156897825634086729D+00
    x(48) = -0.3841027957915169357790778D+00
    x(49) = -0.3612388886058697060709248D+00
    x(50) = -0.3381556747203985013760003D+00
    x(51) = -0.3148671678628949814860148D+00
    x(52) = -0.2913875063937056207945188D+00
    x(53) = -0.2677309447223886208883435D+00
    x(54) = -0.2439118446539178579707132D+00
    x(55) = -0.2199446666696875424545234D+00
    x(56) = -0.1958439611486108515042816D+00
    x(57) = -0.1716243595336421650083449D+00
    x(58) = -0.1473005654490856693893293D+00
    x(59) = -0.1228873457740829717260337D+00
    x(60) = -0.0983995216776989707510918D+00
    x(61) = -0.0738519596210485452734404D+00
    x(62) = -0.0492595623319266303153793D+00
    x(63) = -0.0246372597574209446148971D+00
    x(64) = 0.0000000000000000000000000D+00
    x(65) = 0.0246372597574209446148971D+00
    x(66) = 0.0492595623319266303153793D+00
    x(67) = 0.0738519596210485452734404D+00
    x(68) = 0.0983995216776989707510918D+00
    x(69) = 0.1228873457740829717260337D+00
    x(70) = 0.1473005654490856693893293D+00
    x(71) = 0.1716243595336421650083449D+00
    x(72) = 0.1958439611486108515042816D+00
    x(73) = 0.2199446666696875424545234D+00
    x(74) = 0.2439118446539178579707132D+00
    x(75) = 0.2677309447223886208883435D+00
    x(76) = 0.2913875063937056207945188D+00
    x(77) = 0.3148671678628949814860148D+00
    x(78) = 0.3381556747203985013760003D+00
    x(79) = 0.3612388886058697060709248D+00
    x(80) = 0.3841027957915169357790778D+00
    x(81) = 0.4067335156897825634086729D+00
    x(82) = 0.4291173092801933762625441D+00
    x(83) = 0.4512405874502662273318986D+00
    x(84) = 0.4730899192454052416450999D+00
    x(85) = 0.4946520400227821173949402D+00
    x(86) = 0.5159138595042493572772773D+00
    x(87) = 0.5368624697233975674581664D+00
    x(88) = 0.5574851528619322329218619D+00
    x(89) = 0.5777693889706125800032517D+00
    x(90) = 0.5977028635700652293844120D+00
    x(91) = 0.6172734751268582838576392D+00
    x(92) = 0.6364693424002972413476082D+00
    x(93) = 0.6552788116554826302767651D+00
    x(94) = 0.6736904637382504853466825D+00
    x(95) = 0.6916931210077006701564414D+00
    x(96) = 0.7092758541221045609994446D+00
    x(97) = 0.7264279886740726855356929D+00
    x(98) = 0.7431391116709545129205669D+00
    x(99) = 0.7593990778565366715566637D+00
    x(100) = 0.7751980158702023824449628D+00
    x(101) = 0.7905263342398137999454500D+00
    x(102) = 0.8053747272046802146665608D+00
    x(103) = 0.8197341803650786741551191D+00
    x(104) = 0.8335959761548995143795572D+00
    x(105) = 0.8469516991340975984533393D+00
    x(106) = 0.8597932410977408098120313D+00
    x(107) = 0.8721128059985607114196375D+00
    x(108) = 0.8839029146800265699452579D+00
    x(109) = 0.8951564094170837089690438D+00
    x(110) = 0.9058664582618213828024613D+00
    x(111) = 0.9160265591914658093130886D+00
    x(112) = 0.9256305440562338491274647D+00
    x(113) = 0.9346725823247379685736349D+00
    x(114) = 0.9431471846248148273454496D+00
    x(115) = 0.9510492060778803105479076D+00
    x(116) = 0.9583738494252387711491029D+00
    x(117) = 0.965116667945292121090825D+00
    x(118) = 0.971273568161529192288947D+00
    x(119) = 0.976840812343070326817444D+00
    x(120) = 0.981815020803814110033463D+00
    x(121) = 0.986193174016931666710438D+00
    x(122) = 0.989972614591484157607787D+00
    x(123) = 0.993151049254517147361131D+00
    x(124) = 0.995726551352027226635433D+00
    x(125) = 0.997697566189804621074417D+00
    x(126) = 0.999062934355311895138316D+00
    x(127) = 0.999822130415306146267351D+00

    w(1) = 0.00045645726109586662791936519265D+00
    w(2) = 0.00106227668695384869596523598532D+00
    w(3) = 0.0016683488125171936761028862915D+00
    w(4) = 0.0022734860707492547802810840776D+00
    w(5) = 0.0028772587656289004082883197514D+00
    w(6) = 0.0034792893810051465908910894100D+00
    w(7) = 0.0040792095178254605327114733457D+00
    w(8) = 0.0046766539777779034772638165663D+00
    w(9) = 0.0052712596565634400891303815906D+00
    w(10) = 0.0058626653903523901033648343751D+00
    w(11) = 0.0064505120486899171845442463869D+00
    w(12) = 0.0070344427036681608755685893033D+00
    w(13) = 0.0076141028256526859356393930849D+00
    w(14) = 0.0081891404887415730817235884719D+00
    w(15) = 0.0087592065795403145773316804234D+00
    w(16) = 0.0093239550065309714787536985834D+00
    w(17) = 0.0098830429087554914716648010900D+00
    w(18) = 0.0104361308631410052256731719977D+00
    w(19) = 0.0109828830900689757887996573761D+00
    w(20) = 0.011522967656921087154811609735D+00
    w(21) = 0.012056056679400848183529562145D+00
    w(22) = 0.012581826520465013101514365424D+00
    w(23) = 0.013099957986718627426172681913D+00
    w(24) = 0.013610136522139249906034237534D+00
    w(25) = 0.014112052399003395774044161634D+00
    w(26) = 0.014605400905893418351737288079D+00
    w(27) = 0.015089882532666922992635733981D+00
    w(28) = 0.015565203152273955098532590263D+00
    w(29) = 0.016031074199309941802254151843D+00
    w(30) = 0.016487212845194879399346060358D+00
    w(31) = 0.016933342169871654545878815295D+00
    w(32) = 0.017369191329918731922164721250D+00
    w(33) = 0.017794495722974774231027912900D+00
    w(34) = 0.018208997148375106468721469154D+00
    w(35) = 0.018612443963902310429440419899D+00
    w(36) = 0.019004591238555646611148901045D+00
    w(37) = 0.019385200901246454628112623489D+00
    w(38) = 0.019754041885329183081815217323D+00
    w(39) = 0.020110890268880247225644623956D+00
    w(40) = 0.020455529410639508279497065713D+00
    w(41) = 0.020787750081531811812652137291D+00
    w(42) = 0.021107350591688713643523847922D+00
    w(43) = 0.021414136912893259295449693234D+00
    w(44) = 0.021707922796373466052301324695D+00
    w(45) = 0.021988529885872983756478409759D+00
    w(46) = 0.022255787825930280235631416460D+00
    w(47) = 0.022509534365300608085694429903D+00
    w(48) = 0.022749615455457959852242553241D+00
    w(49) = 0.022975885344117206754377437839D+00
    w(50) = 0.023188206663719640249922582982D+00
    w(51) = 0.023386450514828194170722043497D+00
    w(52) = 0.023570496544381716050033676844D+00
    w(53) = 0.023740233018760777777714726703D+00
    w(54) = 0.023895556891620665983864481754D+00
    w(55) = 0.024036373866450369675132086026D+00
    w(56) = 0.024162598453819584716522917711D+00
    w(57) = 0.024274154023278979833195063937D+00
    w(58) = 0.024370972849882214952813561907D+00
    w(59) = 0.024452996155301467956140198472D+00
    w(60) = 0.024520174143511508275183033290D+00
    w(61) = 0.024572466031020653286354137335D+00
    w(62) = 0.024609840071630254092545634003D+00
    w(63) = 0.024632273575707679066033370218D+00
    w(64) = 0.02463975292396109441957941748D+00
    w(65) = 0.024632273575707679066033370218D+00
    w(66) = 0.024609840071630254092545634003D+00
    w(67) = 0.024572466031020653286354137335D+00
    w(68) = 0.024520174143511508275183033290D+00
    w(69) = 0.024452996155301467956140198472D+00
    w(70) = 0.024370972849882214952813561907D+00
    w(71) = 0.024274154023278979833195063937D+00
    w(72) = 0.024162598453819584716522917711D+00
    w(73) = 0.024036373866450369675132086026D+00
    w(74) = 0.023895556891620665983864481754D+00
    w(75) = 0.023740233018760777777714726703D+00
    w(76) = 0.023570496544381716050033676844D+00
    w(77) = 0.023386450514828194170722043497D+00
    w(78) = 0.023188206663719640249922582982D+00
    w(79) = 0.022975885344117206754377437839D+00
    w(80) = 0.022749615455457959852242553241D+00
    w(81) = 0.022509534365300608085694429903D+00
    w(82) = 0.022255787825930280235631416460D+00
    w(83) = 0.021988529885872983756478409759D+00
    w(84) = 0.021707922796373466052301324695D+00
    w(85) = 0.021414136912893259295449693234D+00
    w(86) = 0.021107350591688713643523847922D+00
    w(87) = 0.020787750081531811812652137291D+00
    w(88) = 0.020455529410639508279497065713D+00
    w(89) = 0.020110890268880247225644623956D+00
    w(90) = 0.019754041885329183081815217323D+00
    w(91) = 0.019385200901246454628112623489D+00
    w(92) = 0.019004591238555646611148901045D+00
    w(93) = 0.018612443963902310429440419899D+00
    w(94) = 0.018208997148375106468721469154D+00
    w(95) = 0.017794495722974774231027912900D+00
    w(96) = 0.017369191329918731922164721250D+00
    w(97) = 0.016933342169871654545878815295D+00
    w(98) = 0.016487212845194879399346060358D+00
    w(99) = 0.016031074199309941802254151843D+00
    w(100) = 0.015565203152273955098532590263D+00
    w(101) = 0.015089882532666922992635733981D+00
    w(102) = 0.014605400905893418351737288079D+00
    w(103) = 0.014112052399003395774044161634D+00
    w(104) = 0.013610136522139249906034237534D+00
    w(105) = 0.013099957986718627426172681913D+00
    w(106) = 0.012581826520465013101514365424D+00
    w(107) = 0.012056056679400848183529562145D+00
    w(108) = 0.011522967656921087154811609735D+00
    w(109) = 0.0109828830900689757887996573761D+00
    w(110) = 0.0104361308631410052256731719977D+00
    w(111) = 0.0098830429087554914716648010900D+00
    w(112) = 0.0093239550065309714787536985834D+00
    w(113) = 0.0087592065795403145773316804234D+00
    w(114) = 0.0081891404887415730817235884719D+00
    w(115) = 0.0076141028256526859356393930849D+00
    w(116) = 0.0070344427036681608755685893033D+00
    w(117) = 0.0064505120486899171845442463869D+00
    w(118) = 0.0058626653903523901033648343751D+00
    w(119) = 0.0052712596565634400891303815906D+00
    w(120) = 0.0046766539777779034772638165663D+00
    w(121) = 0.0040792095178254605327114733457D+00
    w(122) = 0.0034792893810051465908910894100D+00
    w(123) = 0.0028772587656289004082883197514D+00
    w(124) = 0.0022734860707492547802810840776D+00
    w(125) = 0.0016683488125171936761028862915D+00
    w(126) = 0.00106227668695384869596523598532D+00
    w(127) = 0.00045645726109586662791936519265D+00

  else if ( n == 128 ) then

    x(1) = -0.9998248879471319144736081D+00
    x(2) = -0.9990774599773758950119878D+00
    x(3) = -0.9977332486255140198821574D+00
    x(4) = -0.9957927585349811868641612D+00
    x(5) = -0.9932571129002129353034372D+00
    x(6) = -0.9901278184917343833379303D+00
    x(7) = -0.9864067427245862088712355D+00
    x(8) = -0.9820961084357185360247656D+00
    x(9) = -0.9771984914639073871653744D+00
    x(10) = -0.9717168187471365809043384D+00
    x(11) = -0.9656543664319652686458290D+00
    x(12) = -0.9590147578536999280989185D+00
    x(13) = -0.9518019613412643862177963D+00
    x(14) = -0.9440202878302201821211114D+00
    x(15) = -0.9356743882779163757831268D+00
    x(16) = -0.9267692508789478433346245D+00
    x(17) = -0.9173101980809605370364836D+00
    x(18) = -0.9073028834017568139214859D+00
    x(19) = -0.8967532880491581843864474D+00
    x(20) = -0.8856677173453972174082924D+00
    x(21) = -0.8740527969580317986954180D+00
    x(22) = -0.8619154689395484605906323D+00
    x(23) = -0.8492629875779689691636001D+00
    x(24) = -0.8361029150609068471168753D+00
    x(25) = -0.8224431169556438424645942D+00
    x(26) = -0.8082917575079136601196422D+00
    x(27) = -0.7936572947621932902433329D+00
    x(28) = -0.7785484755064119668504941D+00
    x(29) = -0.7629743300440947227797691D+00
    x(30) = -0.7469441667970619811698824D+00
    x(31) = -0.7304675667419088064717369D+00
    x(32) = -0.7135543776835874133438599D+00
    x(33) = -0.6962147083695143323850866D+00
    x(34) = -0.6784589224477192593677557D+00
    x(35) = -0.6602976322726460521059468D+00
    x(36) = -0.6417416925623075571535249D+00
    x(37) = -0.6228021939105849107615396D+00
    x(38) = -0.6034904561585486242035732D+00
    x(39) = -0.5838180216287630895500389D+00
    x(40) = -0.5637966482266180839144308D+00
    x(41) = -0.5434383024128103634441936D+00
    x(42) = -0.5227551520511754784539479D+00
    x(43) = -0.5017595591361444642896063D+00
    x(44) = -0.4804640724041720258582757D+00
    x(45) = -0.4588814198335521954490891D+00
    x(46) = -0.4370245010371041629370429D+00
    x(47) = -0.4149063795522750154922739D+00
    x(48) = -0.3925402750332674427356482D+00
    x(49) = -0.3699395553498590266165917D+00
    x(50) = -0.3471177285976355084261628D+00
    x(51) = -0.3240884350244133751832523D+00
    x(52) = -0.3008654388776772026671541D+00
    x(53) = -0.2774626201779044028062316D+00
    x(54) = -0.2538939664226943208556180D+00
    x(55) = -0.2301735642266599864109866D+00
    x(56) = -0.2063155909020792171540580D+00
    x(57) = -0.1823343059853371824103826D+00
    x(58) = -0.1582440427142249339974755D+00
    x(59) = -0.1340591994611877851175753D+00
    x(60) = -0.1097942311276437466729747D+00
    x(61) = -0.0854636405045154986364980D+00
    x(62) = -0.0610819696041395681037870D+00
    x(63) = -0.0366637909687334933302153D+00
    x(64) = -0.0122236989606157641980521D+00
    x(65) = 0.0122236989606157641980521D+00
    x(66) = 0.0366637909687334933302153D+00
    x(67) = 0.0610819696041395681037870D+00
    x(68) = 0.0854636405045154986364980D+00
    x(69) = 0.1097942311276437466729747D+00
    x(70) = 0.1340591994611877851175753D+00
    x(71) = 0.1582440427142249339974755D+00
    x(72) = 0.1823343059853371824103826D+00
    x(73) = 0.2063155909020792171540580D+00
    x(74) = 0.2301735642266599864109866D+00
    x(75) = 0.2538939664226943208556180D+00
    x(76) = 0.2774626201779044028062316D+00
    x(77) = 0.3008654388776772026671541D+00
    x(78) = 0.3240884350244133751832523D+00
    x(79) = 0.3471177285976355084261628D+00
    x(80) = 0.3699395553498590266165917D+00
    x(81) = 0.3925402750332674427356482D+00
    x(82) = 0.4149063795522750154922739D+00
    x(83) = 0.4370245010371041629370429D+00
    x(84) = 0.4588814198335521954490891D+00
    x(85) = 0.4804640724041720258582757D+00
    x(86) = 0.5017595591361444642896063D+00
    x(87) = 0.5227551520511754784539479D+00
    x(88) = 0.5434383024128103634441936D+00
    x(89) = 0.5637966482266180839144308D+00
    x(90) = 0.5838180216287630895500389D+00
    x(91) = 0.6034904561585486242035732D+00
    x(92) = 0.6228021939105849107615396D+00
    x(93) = 0.6417416925623075571535249D+00
    x(94) = 0.6602976322726460521059468D+00
    x(95) = 0.6784589224477192593677557D+00
    x(96) = 0.6962147083695143323850866D+00
    x(97) = 0.7135543776835874133438599D+00
    x(98) = 0.7304675667419088064717369D+00
    x(99) = 0.7469441667970619811698824D+00
    x(100) = 0.7629743300440947227797691D+00
    x(101) = 0.7785484755064119668504941D+00
    x(102) = 0.7936572947621932902433329D+00
    x(103) = 0.8082917575079136601196422D+00
    x(104) = 0.8224431169556438424645942D+00
    x(105) = 0.8361029150609068471168753D+00
    x(106) = 0.8492629875779689691636001D+00
    x(107) = 0.8619154689395484605906323D+00
    x(108) = 0.8740527969580317986954180D+00
    x(109) = 0.8856677173453972174082924D+00
    x(110) = 0.8967532880491581843864474D+00
    x(111) = 0.9073028834017568139214859D+00
    x(112) = 0.9173101980809605370364836D+00
    x(113) = 0.926769250878947843334625D+00
    x(114) = 0.935674388277916375783127D+00
    x(115) = 0.944020287830220182121111D+00
    x(116) = 0.951801961341264386217796D+00
    x(117) = 0.959014757853699928098919D+00
    x(118) = 0.965654366431965268645829D+00
    x(119) = 0.971716818747136580904338D+00
    x(120) = 0.977198491463907387165374D+00
    x(121) = 0.982096108435718536024766D+00
    x(122) = 0.986406742724586208871236D+00
    x(123) = 0.990127818491734383337930D+00
    x(124) = 0.993257112900212935303437D+00
    x(125) = 0.995792758534981186864161D+00
    x(126) = 0.997733248625514019882157D+00
    x(127) = 0.999077459977375895011988D+00
    x(128) = 0.999824887947131914473608D+00

    w(1) = 0.00044938096029209037639429223999D+00
    w(2) = 0.0010458126793403487793128516001D+00
    w(3) = 0.0016425030186690295387908755948D+00
    w(4) = 0.0022382884309626187436220542727D+00
    w(5) = 0.0028327514714579910952857346468D+00
    w(6) = 0.0034255260409102157743377846601D+00
    w(7) = 0.0040162549837386423131943434863D+00
    w(8) = 0.0046045842567029551182905419803D+00
    w(9) = 0.0051901618326763302050707671348D+00
    w(10) = 0.0057726375428656985893346176261D+00
    w(11) = 0.006351663161707188787214327826D+00
    w(12) = 0.006926892566898813563426670360D+00
    w(13) = 0.007497981925634728687671962688D+00
    w(14) = 0.008064589890486057972928598698D+00
    w(15) = 0.008626377798616749704978843782D+00
    w(16) = 0.009183009871660874334478743688D+00
    w(17) = 0.009734153415006805863548266094D+00
    w(18) = 0.010279479015832157133215340326D+00
    w(19) = 0.010818660739503076247659646277D+00
    w(20) = 0.011351376324080416693281668453D+00
    w(21) = 0.011877307372740279575891106926D+00
    w(22) = 0.012396139543950922968821728197D+00
    w(23) = 0.012907562739267347220442834004D+00
    w(24) = 0.013411271288616332314488951616D+00
    w(25) = 0.013906964132951985244288007396D+00
    w(26) = 0.014394345004166846176823892009D+00
    w(27) = 0.014873122602147314252385498520D+00
    w(28) = 0.015343010768865144085990853741D+00
    w(29) = 0.015803728659399346858965631687D+00
    w(30) = 0.016255000909785187051657456477D+00
    w(31) = 0.016696557801589204589091507954D+00
    w(32) = 0.017128135423111376830680987619D+00
    w(33) = 0.017549475827117704648706925634D+00
    w(34) = 0.017960327185008685940196927525D+00
    w(35) = 0.018360443937331343221289290991D+00
    w(36) = 0.018749586940544708650919548474D+00
    w(37) = 0.019127523609950945486518531668D+00
    w(38) = 0.019494028058706602823021918681D+00
    w(39) = 0.019848881232830862219944413265D+00
    w(40) = 0.020191871042130041180673158406D+00
    w(41) = 0.020522792486960069432284967788D+00
    w(42) = 0.020841447780751149113583948423D+00
    w(43) = 0.021147646468221348537019535180D+00
    w(44) = 0.021441205539208460137111853878D+00
    w(45) = 0.021721949538052075375260957768D+00
    w(46) = 0.021989710668460491434122106599D+00
    w(47) = 0.022244328893799765104629133607D+00
    w(48) = 0.022485652032744966871824603941D+00
    w(49) = 0.022713535850236461309712635923D+00
    w(50) = 0.022927844143686846920410987209D+00
    w(51) = 0.023128448824387027879297902403D+00
    w(52) = 0.023315229994062760122415671273D+00
    w(53) = 0.023488076016535913153025273282D+00
    w(54) = 0.023646883584447615143651392303D+00
    w(55) = 0.023791557781003400638780709885D+00
    w(56) = 0.023922012136703455672450408817D+00
    w(57) = 0.024038168681024052637587316820D+00
    w(58) = 0.024139957989019284997716653890D+00
    w(59) = 0.024227319222815248120093308442D+00
    w(60) = 0.024300200167971865323442606364D+00
    w(61) = 0.024358557264690625853268520246D+00
    w(62) = 0.024402355633849582093297989694D+00
    w(63) = 0.02443156909785004505484856143D+00
    w(64) = 0.02444618019626251821132585261D+00
    w(65) = 0.02444618019626251821132585261D+00
    w(66) = 0.02443156909785004505484856143D+00
    w(67) = 0.024402355633849582093297989694D+00
    w(68) = 0.024358557264690625853268520246D+00
    w(69) = 0.024300200167971865323442606364D+00
    w(70) = 0.024227319222815248120093308442D+00
    w(71) = 0.024139957989019284997716653890D+00
    w(72) = 0.024038168681024052637587316820D+00
    w(73) = 0.023922012136703455672450408817D+00
    w(74) = 0.023791557781003400638780709885D+00
    w(75) = 0.023646883584447615143651392303D+00
    w(76) = 0.023488076016535913153025273282D+00
    w(77) = 0.023315229994062760122415671273D+00
    w(78) = 0.023128448824387027879297902403D+00
    w(79) = 0.022927844143686846920410987209D+00
    w(80) = 0.022713535850236461309712635923D+00
    w(81) = 0.022485652032744966871824603941D+00
    w(82) = 0.022244328893799765104629133607D+00
    w(83) = 0.021989710668460491434122106599D+00
    w(84) = 0.021721949538052075375260957768D+00
    w(85) = 0.021441205539208460137111853878D+00
    w(86) = 0.021147646468221348537019535180D+00
    w(87) = 0.020841447780751149113583948423D+00
    w(88) = 0.020522792486960069432284967788D+00
    w(89) = 0.020191871042130041180673158406D+00
    w(90) = 0.019848881232830862219944413265D+00
    w(91) = 0.019494028058706602823021918681D+00
    w(92) = 0.019127523609950945486518531668D+00
    w(93) = 0.018749586940544708650919548474D+00
    w(94) = 0.018360443937331343221289290991D+00
    w(95) = 0.017960327185008685940196927525D+00
    w(96) = 0.017549475827117704648706925634D+00
    w(97) = 0.017128135423111376830680987619D+00
    w(98) = 0.016696557801589204589091507954D+00
    w(99) = 0.016255000909785187051657456477D+00
    w(100) = 0.015803728659399346858965631687D+00
    w(101) = 0.015343010768865144085990853741D+00
    w(102) = 0.014873122602147314252385498520D+00
    w(103) = 0.014394345004166846176823892009D+00
    w(104) = 0.013906964132951985244288007396D+00
    w(105) = 0.013411271288616332314488951616D+00
    w(106) = 0.012907562739267347220442834004D+00
    w(107) = 0.012396139543950922968821728197D+00
    w(108) = 0.011877307372740279575891106926D+00
    w(109) = 0.011351376324080416693281668453D+00
    w(110) = 0.010818660739503076247659646277D+00
    w(111) = 0.010279479015832157133215340326D+00
    w(112) = 0.009734153415006805863548266094D+00
    w(113) = 0.009183009871660874334478743688D+00
    w(114) = 0.008626377798616749704978843782D+00
    w(115) = 0.008064589890486057972928598698D+00
    w(116) = 0.007497981925634728687671962688D+00
    w(117) = 0.006926892566898813563426670360D+00
    w(118) = 0.006351663161707188787214327826D+00
    w(119) = 0.0057726375428656985893346176261D+00
    w(120) = 0.0051901618326763302050707671348D+00
    w(121) = 0.0046045842567029551182905419803D+00
    w(122) = 0.0040162549837386423131943434863D+00
    w(123) = 0.0034255260409102157743377846601D+00
    w(124) = 0.0028327514714579910952857346468D+00
    w(125) = 0.0022382884309626187436220542727D+00
    w(126) = 0.0016425030186690295387908755948D+00
    w(127) = 0.0010458126793403487793128516001D+00
    w(128) = 0.00044938096029209037639429223999D+00

  else if ( n == 129 ) then

    x(1) = -0.9998275818477487191077441D+00
    x(2) = -0.9990916504696409986514389D+00
    x(3) = -0.9977681080525852721429460D+00
    x(4) = -0.9958574393142831982149111D+00
    x(5) = -0.9933607326210712814854011D+00
    x(6) = -0.9902794486488178389207689D+00
    x(7) = -0.9866153978313475022005761D+00
    x(8) = -0.9823707352517413115507418D+00
    x(9) = -0.9775479582993672474447814D+00
    x(10) = -0.9721499048427034297274163D+00
    x(11) = -0.9661797514202097197778763D+00
    x(12) = -0.9596410113101918904168119D+00
    x(13) = -0.9525375324342090471027732D+00
    x(14) = -0.9448734950776734726784764D+00
    x(15) = -0.9366534094216514605284616D+00
    x(16) = -0.9278821128840036204317296D+00
    x(17) = -0.9185647672698286252225115D+00
    x(18) = -0.9087068557320696331245539D+00
    x(19) = -0.8983141795436338850435985D+00
    x(20) = -0.8873928546826803665034968D+00
    x(21) = -0.8759493082329433892035217D+00
    x(22) = -0.8639902746011257878940216D+00
    x(23) = -0.8515227915535356930243826D+00
    x(24) = -0.8385541960742664442975407D+00
    x(25) = -0.8250921200473358809210133D+00
    x(26) = -0.8111444857653120742087717D+00
    x(27) = -0.7967195012670592680339606D+00
    x(28) = -0.7818256555073413245387500D+00
    x(29) = -0.7664717133611208816717785D+00
    x(30) = -0.7506667104654910227632368D+00
    x(31) = -0.7344199479022727047791516D+00
    x(32) = -0.7177409867244055767721220D+00
    x(33) = -0.7006396423293521790044710D+00
    x(34) = -0.6831259786828258512462248D+00
    x(35) = -0.6652103023962409818802202D+00
    x(36) = -0.6469031566613704719753373D+00
    x(37) = -0.6282153150457794374886895D+00
    x(38) = -0.6091577751526861909563306D+00
    x(39) = -0.5897417521489813916767844D+00
    x(40) = -0.5699786721652138894754096D+00
    x(41) = -0.5498801655714271702189358D+00
    x(42) = -0.5294580601328034000099406D+00
    x(43) = -0.5087243740491428186199463D+00
    x(44) = -0.4876913088822746111853066D+00
    x(45) = -0.4663712423755613514331869D+00
    x(46) = -0.4447767211697226217818454D+00
    x(47) = -0.4229204534192644388475065D+00
    x(48) = -0.4008153013138596117693121D+00
    x(49) = -0.3784742735090801012801265D+00
    x(50) = -0.3559105174709357969672656D+00
    x(51) = -0.3331373117387248575049982D+00
    x(52) = -0.3101680581107488341147318D+00
    x(53) = -0.2870162737574911929568755D+00
    x(54) = -0.2636955832669005409666949D+00
    x(55) = -0.2402197106264598167721148D+00
    x(56) = -0.2166024711467599103221439D+00
    x(57) = -0.1928577633313305998663880D+00
    x(58) = -0.1689995606975133227390302D+00
    x(59) = -0.1450419035531891084328306D+00
    x(60) = -0.1209988907342009817690539D+00
    x(61) = -0.0968846713073332753086909D+00
    x(62) = -0.0727134362437305599118207D+00
    x(63) = -0.0484994100676562986191764D+00
    x(64) = -0.0242568424855058415749954D+00
    x(65) = 0.0000000000000000000000000D+00
    x(66) = 0.0242568424855058415749954D+00
    x(67) = 0.0484994100676562986191764D+00
    x(68) = 0.0727134362437305599118207D+00
    x(69) = 0.0968846713073332753086909D+00
    x(70) = 0.1209988907342009817690539D+00
    x(71) = 0.1450419035531891084328306D+00
    x(72) = 0.1689995606975133227390302D+00
    x(73) = 0.1928577633313305998663880D+00
    x(74) = 0.2166024711467599103221439D+00
    x(75) = 0.2402197106264598167721148D+00
    x(76) = 0.2636955832669005409666949D+00
    x(77) = 0.2870162737574911929568755D+00
    x(78) = 0.3101680581107488341147318D+00
    x(79) = 0.3331373117387248575049982D+00
    x(80) = 0.3559105174709357969672656D+00
    x(81) = 0.3784742735090801012801265D+00
    x(82) = 0.4008153013138596117693121D+00
    x(83) = 0.4229204534192644388475065D+00
    x(84) = 0.4447767211697226217818454D+00
    x(85) = 0.4663712423755613514331869D+00
    x(86) = 0.4876913088822746111853066D+00
    x(87) = 0.5087243740491428186199463D+00
    x(88) = 0.5294580601328034000099406D+00
    x(89) = 0.5498801655714271702189358D+00
    x(90) = 0.5699786721652138894754096D+00
    x(91) = 0.5897417521489813916767844D+00
    x(92) = 0.6091577751526861909563306D+00
    x(93) = 0.6282153150457794374886895D+00
    x(94) = 0.6469031566613704719753373D+00
    x(95) = 0.6652103023962409818802202D+00
    x(96) = 0.6831259786828258512462248D+00
    x(97) = 0.7006396423293521790044710D+00
    x(98) = 0.7177409867244055767721220D+00
    x(99) = 0.7344199479022727047791516D+00
    x(100) = 0.7506667104654910227632368D+00
    x(101) = 0.7664717133611208816717785D+00
    x(102) = 0.7818256555073413245387500D+00
    x(103) = 0.7967195012670592680339606D+00
    x(104) = 0.8111444857653120742087717D+00
    x(105) = 0.8250921200473358809210133D+00
    x(106) = 0.8385541960742664442975407D+00
    x(107) = 0.8515227915535356930243826D+00
    x(108) = 0.8639902746011257878940216D+00
    x(109) = 0.875949308232943389203522D+00
    x(110) = 0.887392854682680366503497D+00
    x(111) = 0.898314179543633885043599D+00
    x(112) = 0.908706855732069633124554D+00
    x(113) = 0.918564767269828625222511D+00
    x(114) = 0.927882112884003620431730D+00
    x(115) = 0.936653409421651460528462D+00
    x(116) = 0.944873495077673472678476D+00
    x(117) = 0.952537532434209047102773D+00
    x(118) = 0.959641011310191890416812D+00
    x(119) = 0.966179751420209719777876D+00
    x(120) = 0.972149904842703429727416D+00
    x(121) = 0.977547958299367247444781D+00
    x(122) = 0.982370735251741311550742D+00
    x(123) = 0.986615397831347502200576D+00
    x(124) = 0.990279448648817838920769D+00
    x(125) = 0.993360732621071281485401D+00
    x(126) = 0.995857439314283198214911D+00
    x(127) = 0.997768108052585272142946D+00
    x(128) = 0.999091650469640998651439D+00
    x(129) = 0.999827581847748719107744D+00

    w(1) = 0.00044246794182939296923668005717D+00
    w(2) = 0.00102972844619622394463273519315D+00
    w(3) = 0.0016172530556785534682413679271D+00
    w(4) = 0.0022039015180966937075786419741D+00
    w(5) = 0.0027892681877797554940944677057D+00
    w(6) = 0.0033729979506246246117755709288D+00
    w(7) = 0.0039547444682113562172392974765D+00
    w(8) = 0.0045341644298525434513226874954D+00
    w(9) = 0.0051109164669246267289761565766D+00
    w(10) = 0.0056846609912469045788016012203D+00
    w(11) = 0.0062550602724461408889348709586D+00
    w(12) = 0.0068217785893519121070498527769D+00
    w(13) = 0.0073844824072454014447165055698D+00
    w(14) = 0.0079428405646668029041114107832D+00
    w(15) = 0.0084965244635723279730542832506D+00
    w(16) = 0.0090452082602137316404219313819D+00
    w(17) = 0.0095885690555104190787301294510D+00
    w(18) = 0.0101262870842733548093160774580D+00
    w(19) = 0.0106580459029055185304204093001D+00
    w(20) = 0.0111835325753305049735380697538D+00
    w(21) = 0.011702437856964778185746436834D+00
    w(22) = 0.012214456376582979416221105914D+00
    w(23) = 0.012719286815944623465099036330D+00
    w(24) = 0.013216632087061724231482387345D+00
    w(25) = 0.013706199506993971244060563234D+00
    w(26) = 0.014187700970062900419317230938D+00
    w(27) = 0.014660853117380060971041027493D+00
    w(28) = 0.015125377503587024690403432771D+00
    w(29) = 0.015581000760707523415881287558D+00
    w(30) = 0.016027454759014214436403950465D+00
    w(31) = 0.016464476764814667467169189640D+00
    w(32) = 0.016891809595063204177526208819D+00
    w(33) = 0.017309201768707240731293596444D+00
    w(34) = 0.017716407654678809269702031810D+00
    w(35) = 0.018113187616443980503999783812D+00
    w(36) = 0.018499308153024985727791918518D+00
    w(37) = 0.018874542036411948181617592169D+00
    w(38) = 0.019238668445283284085199492202D+00
    w(39) = 0.019591473094956024580283987216D+00
    w(40) = 0.019932748363489542089706675388D+00
    w(41) = 0.020262293413868438317104423081D+00
    w(42) = 0.020579914312192665948185517085D+00
    w(43) = 0.020885424141805311409990024684D+00
    w(44) = 0.021178643113290860912881038703D+00
    w(45) = 0.021459398670279205389981598196D+00
    w(46) = 0.021727525590993110687305178710D+00
    w(47) = 0.021982866085479386179554968899D+00
    w(48) = 0.022225269888466526554736910919D+00
    w(49) = 0.022454594347794176432066564511D+00
    w(50) = 0.022670704508362374313093970958D+00
    w(51) = 0.022873473191551169638592083492D+00
    w(52) = 0.023062781070063872924670495006D+00
    w(53) = 0.023238516738149892544490435771D+00
    w(54) = 0.023400576777165831146714346635D+00
    w(55) = 0.023548865816436258377269094263D+00
    w(56) = 0.023683296589378342897341543485D+00
    w(57) = 0.023803789984857314051325299744D+00
    w(58) = 0.023910275093742530302367230296D+00
    w(59) = 0.024002689250636756075547029720D+00
    w(60) = 0.024080978070754089272959634041D+00
    w(61) = 0.024145095481924836783843156014D+00
    w(62) = 0.024195003751708503129818111597D+00
    w(63) = 0.024230673509598936275508460625D+00
    w(64) = 0.024252083764308562906498864071D+00
    w(65) = 0.02425922191612154143202867472D+00
    w(66) = 0.024252083764308562906498864071D+00
    w(67) = 0.024230673509598936275508460625D+00
    w(68) = 0.024195003751708503129818111597D+00
    w(69) = 0.024145095481924836783843156014D+00
    w(70) = 0.024080978070754089272959634041D+00
    w(71) = 0.024002689250636756075547029720D+00
    w(72) = 0.023910275093742530302367230296D+00
    w(73) = 0.023803789984857314051325299744D+00
    w(74) = 0.023683296589378342897341543485D+00
    w(75) = 0.023548865816436258377269094263D+00
    w(76) = 0.023400576777165831146714346635D+00
    w(77) = 0.023238516738149892544490435771D+00
    w(78) = 0.023062781070063872924670495006D+00
    w(79) = 0.022873473191551169638592083492D+00
    w(80) = 0.022670704508362374313093970958D+00
    w(81) = 0.022454594347794176432066564511D+00
    w(82) = 0.022225269888466526554736910919D+00
    w(83) = 0.021982866085479386179554968899D+00
    w(84) = 0.021727525590993110687305178710D+00
    w(85) = 0.021459398670279205389981598196D+00
    w(86) = 0.021178643113290860912881038703D+00
    w(87) = 0.020885424141805311409990024684D+00
    w(88) = 0.020579914312192665948185517085D+00
    w(89) = 0.020262293413868438317104423081D+00
    w(90) = 0.019932748363489542089706675388D+00
    w(91) = 0.019591473094956024580283987216D+00
    w(92) = 0.019238668445283284085199492202D+00
    w(93) = 0.018874542036411948181617592169D+00
    w(94) = 0.018499308153024985727791918518D+00
    w(95) = 0.018113187616443980503999783812D+00
    w(96) = 0.017716407654678809269702031810D+00
    w(97) = 0.017309201768707240731293596444D+00
    w(98) = 0.016891809595063204177526208819D+00
    w(99) = 0.016464476764814667467169189640D+00
    w(100) = 0.016027454759014214436403950465D+00
    w(101) = 0.015581000760707523415881287558D+00
    w(102) = 0.015125377503587024690403432771D+00
    w(103) = 0.014660853117380060971041027493D+00
    w(104) = 0.014187700970062900419317230938D+00
    w(105) = 0.013706199506993971244060563234D+00
    w(106) = 0.013216632087061724231482387345D+00
    w(107) = 0.012719286815944623465099036330D+00
    w(108) = 0.012214456376582979416221105914D+00
    w(109) = 0.011702437856964778185746436834D+00
    w(110) = 0.0111835325753305049735380697538D+00
    w(111) = 0.0106580459029055185304204093001D+00
    w(112) = 0.0101262870842733548093160774580D+00
    w(113) = 0.0095885690555104190787301294510D+00
    w(114) = 0.0090452082602137316404219313819D+00
    w(115) = 0.0084965244635723279730542832506D+00
    w(116) = 0.0079428405646668029041114107832D+00
    w(117) = 0.0073844824072454014447165055698D+00
    w(118) = 0.0068217785893519121070498527769D+00
    w(119) = 0.0062550602724461408889348709586D+00
    w(120) = 0.0056846609912469045788016012203D+00
    w(121) = 0.0051109164669246267289761565766D+00
    w(122) = 0.0045341644298525434513226874954D+00
    w(123) = 0.0039547444682113562172392974765D+00
    w(124) = 0.0033729979506246246117755709288D+00
    w(125) = 0.0027892681877797554940944677057D+00
    w(126) = 0.0022039015180966937075786419741D+00
    w(127) = 0.0016172530556785534682413679271D+00
    w(128) = 0.00102972844619622394463273519315D+00
    w(129) = 0.00044246794182939296923668005717D+00

  else if ( n == 255 ) then

    x(1) = -0.999955705317563751730191D+00
    x(2) = -0.999766621312000569367063D+00
    x(3) = -0.999426474680169959344386D+00
    x(4) = -0.998935241284654635142155D+00
    x(5) = -0.998292986136967889228248D+00
    x(6) = -0.997499804126615814044844D+00
    x(7) = -0.996555814435198617028738D+00
    x(8) = -0.995461159480026294089975D+00
    x(9) = -0.994216004616630164799381D+00
    x(10) = -0.992820538021989138984811D+00
    x(11) = -0.991274970630385567164523D+00
    x(12) = -0.989579536085920123498574D+00
    x(13) = -0.987734490699732356281248D+00
    x(14) = -0.985740113407419277752900D+00
    x(15) = -0.983596705724776358640192D+00
    x(16) = -0.981304591701017185126565D+00
    x(17) = -0.978864117869068155239121D+00
    x(18) = -0.976275653192735980815246D+00
    x(19) = -0.973539589010643617645393D+00
    x(20) = -0.970656338976880365477697D+00
    x(21) = -0.967626338998338798105523D+00
    x(22) = -0.964450047168726298761719D+00
    x(23) = -0.961127943699247839572910D+00
    x(24) = -0.957660530845962076295490D+00
    x(25) = -0.954048332833816317950921D+00
    x(26) = -0.950291895777368285733522D+00
    x(27) = -0.946391787598204251752103D+00
    x(28) = -0.942348597939064408301480D+00
    x(29) = -0.938162938074687317626793D+00
    x(30) = -0.933835440819386124349338D+00
    x(31) = -0.929366760431369935739045D+00
    x(32) = -0.924757572513824425220425D+00
    x(33) = -0.920008573912766315142721D+00
    x(34) = -0.915120482611686961035103D+00
    x(35) = -0.910094037623000801254172D+00
    x(36) = -0.904929998876314959753358D+00
    x(37) = -0.899629147103536800144342D+00
    x(38) = -0.894192283720836729335637D+00
    x(39) = -0.888620230707484040924981D+00
    x(40) = -0.882913830481574073645470D+00
    x(41) = -0.877073945772665439532627D+00
    x(42) = -0.871101459491346550796200D+00
    x(43) = -0.864997274595751144137121D+00
    x(44) = -0.858762313955042966785823D+00
    x(45) = -0.852397520209890250084237D+00
    x(46) = -0.845903855629951054143931D+00
    x(47) = -0.839282301968391021084600D+00
    x(48) = -0.832533860313455524647230D+00
    x(49) = -0.825659550937118650611534D+00
    x(50) = -0.818660413140831885432406D+00
    x(51) = -0.811537505098395829833580D+00
    x(52) = -0.804291903695978689734633D+00
    x(53) = -0.796924704369305728807154D+00
    x(54) = -0.789437020938044295117764D+00
    x(55) = -0.781829985437409458675147D+00
    x(56) = -0.774104747947015717207115D+00
    x(57) = -0.766262476417000644100858D+00
    x(58) = -0.758304356491446765092016D+00
    x(59) = -0.750231591329128358931528D+00
    x(60) = -0.742045401421610281838045D+00
    x(61) = -0.733747024408726316001889D+00
    x(62) = -0.725337714891464938687812D+00
    x(63) = -0.716818744242290800531501D+00
    x(64) = -0.708191400412930589382399D+00
    x(65) = -0.699456987739652339456557D+00
    x(66) = -0.690616826746067624571761D+00
    x(67) = -0.681672253943486448787259D+00
    x(68) = -0.672624621628855017806731D+00
    x(69) = -0.663475297680306939970658D+00
    x(70) = -0.654225665350358766508700D+00
    x(71) = -0.644877123056781136890077D+00
    x(72) = -0.635431084171177146547142D+00
    x(73) = -0.625888976805299900901619D+00
    x(74) = -0.616252243595141561442344D+00
    x(75) = -0.606522341482826526536576D+00
    x(76) = -0.596700741496341721653202D+00
    x(77) = -0.586788928527137300685706D+00
    x(78) = -0.576788401105631382036211D+00
    x(79) = -0.566700671174652760010815D+00
    x(80) = -0.556527263860855843833077D+00
    x(81) = -0.546269717244142383159817D+00
    x(82) = -0.535929582125124840335150D+00
    x(83) = -0.525508421790666565699453D+00
    x(84) = -0.515007811777534223035005D+00
    x(85) = -0.504429339634198197635551D+00
    x(86) = -0.493774604680816999489812D+00
    x(87) = -0.483045217767441948626854D+00
    x(88) = -0.472242801030478698742627D+00
    x(89) = -0.461368987647442418771401D+00
    x(90) = -0.450425421590043710043279D+00
    x(91) = -0.439413757375642589040685D+00
    x(92) = -0.428335659817108112494341D+00
    x(93) = -0.417192803771121462605751D+00
    x(94) = -0.405986873884960545511889D+00
    x(95) = -0.394719564341804385683361D+00
    x(96) = -0.383392578604595822734854D+00
    x(97) = -0.372007629158501235092510D+00
    x(98) = -0.360566437252006227074021D+00
    x(99) = -0.349070732636686422161576D+00
    x(100) = -0.337522253305692705554261D+00
    x(101) = -0.325922745230990453444769D+00
    x(102) = -0.314273962099392474845918D+00
    x(103) = -0.302577665047425574167140D+00
    x(104) = -0.290835622395070819082047D+00
    x(105) = -0.279049609378417768508970D+00
    x(106) = -0.267221407881273079721012D+00
    x(107) = -0.255352806165764071686080D+00
    x(108) = -0.243445598601977973686482D+00
    x(109) = -0.231501585396677734059116D+00
    x(110) = -0.219522572321135403508985D+00
    x(111) = -0.207510370438124240859625D+00
    x(112) = -0.195466795828110816293869D+00
    x(113) = -0.183393669314688508087976D+00
    x(114) = -0.171292816189293903533225D+00
    x(115) = -0.159166065935247723154292D+00
    x(116) = -0.147015251951161989456661D+00
    x(117) = -0.134842211273755257250625D+00
    x(118) = -0.122648784300117812092492D+00
    x(119) = -0.110436814509468826540991D+00
    x(120) = -0.098208148184447540736015D+00
    x(121) = -0.085964634131980604256000D+00
    x(122) = -0.073708123403767780288977D+00
    x(123) = -0.061440469016428270850728D+00
    x(124) = -0.049163525671349973093019D+00
    x(125) = -0.036879149474284021657652D+00
    x(126) = -0.024589197654727010541405D+00
    x(127) = -0.012295528285133320036860D+00
    x(128) = 0.000000000000000000000000D+00
    x(129) = 0.012295528285133320036860D+00
    x(130) = 0.024589197654727010541405D+00
    x(131) = 0.036879149474284021657652D+00
    x(132) = 0.049163525671349973093019D+00
    x(133) = 0.061440469016428270850728D+00
    x(134) = 0.073708123403767780288977D+00
    x(135) = 0.085964634131980604256000D+00
    x(136) = 0.098208148184447540736015D+00
    x(137) = 0.110436814509468826540991D+00
    x(138) = 0.122648784300117812092492D+00
    x(139) = 0.134842211273755257250625D+00
    x(140) = 0.147015251951161989456661D+00
    x(141) = 0.159166065935247723154292D+00
    x(142) = 0.171292816189293903533225D+00
    x(143) = 0.183393669314688508087976D+00
    x(144) = 0.195466795828110816293869D+00
    x(145) = 0.207510370438124240859625D+00
    x(146) = 0.219522572321135403508985D+00
    x(147) = 0.231501585396677734059116D+00
    x(148) = 0.243445598601977973686482D+00
    x(149) = 0.255352806165764071686080D+00
    x(150) = 0.267221407881273079721012D+00
    x(151) = 0.279049609378417768508970D+00
    x(152) = 0.290835622395070819082047D+00
    x(153) = 0.302577665047425574167140D+00
    x(154) = 0.314273962099392474845918D+00
    x(155) = 0.325922745230990453444769D+00
    x(156) = 0.337522253305692705554261D+00
    x(157) = 0.349070732636686422161576D+00
    x(158) = 0.360566437252006227074021D+00
    x(159) = 0.372007629158501235092510D+00
    x(160) = 0.383392578604595822734854D+00
    x(161) = 0.394719564341804385683361D+00
    x(162) = 0.405986873884960545511889D+00
    x(163) = 0.417192803771121462605751D+00
    x(164) = 0.428335659817108112494341D+00
    x(165) = 0.439413757375642589040685D+00
    x(166) = 0.450425421590043710043279D+00
    x(167) = 0.461368987647442418771401D+00
    x(168) = 0.472242801030478698742627D+00
    x(169) = 0.483045217767441948626854D+00
    x(170) = 0.493774604680816999489812D+00
    x(171) = 0.504429339634198197635551D+00
    x(172) = 0.515007811777534223035005D+00
    x(173) = 0.525508421790666565699453D+00
    x(174) = 0.535929582125124840335150D+00
    x(175) = 0.546269717244142383159817D+00
    x(176) = 0.556527263860855843833077D+00
    x(177) = 0.566700671174652760010815D+00
    x(178) = 0.576788401105631382036211D+00
    x(179) = 0.586788928527137300685706D+00
    x(180) = 0.596700741496341721653202D+00
    x(181) = 0.606522341482826526536576D+00
    x(182) = 0.616252243595141561442344D+00
    x(183) = 0.625888976805299900901619D+00
    x(184) = 0.635431084171177146547142D+00
    x(185) = 0.644877123056781136890077D+00
    x(186) = 0.654225665350358766508700D+00
    x(187) = 0.663475297680306939970658D+00
    x(188) = 0.672624621628855017806731D+00
    x(189) = 0.681672253943486448787259D+00
    x(190) = 0.690616826746067624571761D+00
    x(191) = 0.699456987739652339456557D+00
    x(192) = 0.708191400412930589382399D+00
    x(193) = 0.716818744242290800531501D+00
    x(194) = 0.725337714891464938687812D+00
    x(195) = 0.733747024408726316001889D+00
    x(196) = 0.742045401421610281838045D+00
    x(197) = 0.750231591329128358931528D+00
    x(198) = 0.758304356491446765092016D+00
    x(199) = 0.766262476417000644100858D+00
    x(200) = 0.774104747947015717207115D+00
    x(201) = 0.781829985437409458675147D+00
    x(202) = 0.789437020938044295117764D+00
    x(203) = 0.796924704369305728807154D+00
    x(204) = 0.804291903695978689734633D+00
    x(205) = 0.811537505098395829833580D+00
    x(206) = 0.818660413140831885432406D+00
    x(207) = 0.825659550937118650611534D+00
    x(208) = 0.832533860313455524647230D+00
    x(209) = 0.839282301968391021084600D+00
    x(210) = 0.845903855629951054143931D+00
    x(211) = 0.852397520209890250084237D+00
    x(212) = 0.858762313955042966785823D+00
    x(213) = 0.864997274595751144137121D+00
    x(214) = 0.871101459491346550796200D+00
    x(215) = 0.877073945772665439532627D+00
    x(216) = 0.882913830481574073645470D+00
    x(217) = 0.888620230707484040924981D+00
    x(218) = 0.894192283720836729335637D+00
    x(219) = 0.899629147103536800144342D+00
    x(220) = 0.904929998876314959753358D+00
    x(221) = 0.910094037623000801254172D+00
    x(222) = 0.915120482611686961035103D+00
    x(223) = 0.920008573912766315142721D+00
    x(224) = 0.924757572513824425220425D+00
    x(225) = 0.929366760431369935739045D+00
    x(226) = 0.933835440819386124349338D+00
    x(227) = 0.938162938074687317626793D+00
    x(228) = 0.942348597939064408301480D+00
    x(229) = 0.946391787598204251752103D+00
    x(230) = 0.950291895777368285733522D+00
    x(231) = 0.954048332833816317950921D+00
    x(232) = 0.957660530845962076295490D+00
    x(233) = 0.961127943699247839572910D+00
    x(234) = 0.964450047168726298761719D+00
    x(235) = 0.967626338998338798105523D+00
    x(236) = 0.970656338976880365477697D+00
    x(237) = 0.973539589010643617645393D+00
    x(238) = 0.976275653192735980815246D+00
    x(239) = 0.978864117869068155239121D+00
    x(240) = 0.981304591701017185126565D+00
    x(241) = 0.983596705724776358640192D+00
    x(242) = 0.985740113407419277752900D+00
    x(243) = 0.987734490699732356281248D+00
    x(244) = 0.989579536085920123498574D+00
    x(245) = 0.991274970630385567164523D+00
    x(246) = 0.992820538021989138984811D+00
    x(247) = 0.994216004616630164799381D+00
    x(248) = 0.995461159480026294089975D+00
    x(249) = 0.996555814435198617028738D+00
    x(250) = 0.997499804126615814044844D+00
    x(251) = 0.998292986136967889228248D+00
    x(252) = 0.998935241284654635142155D+00
    x(253) = 0.999426474680169959344386D+00
    x(254) = 0.999766621312000569367063D+00
    x(255) = 0.999955705317563751730191D+00

    w(1) = 0.00011367361999142272115645954414D+00
    w(2) = 0.00026459387119083065532790838855D+00
    w(3) = 0.00041569762526823913616284210066D+00
    w(4) = 0.00056675794564824918946626058353D+00
    w(5) = 0.00071773647800611087798371518325D+00
    w(6) = 0.00086860766611945667949717690640D+00
    w(7) = 0.00101934797642732530281229369360D+00
    w(8) = 0.0011699343729388079886897709773D+00
    w(9) = 0.0013203439900221692090523602144D+00
    w(10) = 0.0014705540427783843160097204304D+00
    w(11) = 0.0016205417990415653896921100325D+00
    w(12) = 0.0017702845706603213070421243905D+00
    w(13) = 0.0019197597117132050055085980675D+00
    w(14) = 0.0020689446195015801533643667413D+00
    w(15) = 0.0022178167367540171700373764020D+00
    w(16) = 0.0023663535543962867157201855305D+00
    w(17) = 0.0025145326145997073931298921370D+00
    w(18) = 0.0026623315139717112732749157331D+00
    w(19) = 0.0028097279068204407457332299361D+00
    w(20) = 0.0029566995084575002760043344138D+00
    w(21) = 0.0031032240985191112621977893133D+00
    w(22) = 0.0032492795242943133198690930777D+00
    w(23) = 0.0033948437040533928255056951665D+00
    w(24) = 0.0035398946303722552150296713510D+00
    w(25) = 0.0036844103734499176530742235517D+00
    w(26) = 0.0038283690844171626400743524999D+00
    w(27) = 0.0039717489986349171988699773906D+00
    w(28) = 0.0041145284389812475901826468094D+00
    w(29) = 0.0042566858191260658425395494472D+00
    w(30) = 0.0043981996467927779838546384780D+00
    w(31) = 0.0045390485270061921259394035112D+00
    w(32) = 0.0046792111653260640506279893190D+00
    w(33) = 0.0048186663710656988918572043815D+00
    w(34) = 0.0049573930604950563104281084148D+00
    w(35) = 0.0050953702600278273039420404117D+00
    w(36) = 0.0052325771093919661294970523234D+00
    w(37) = 0.0053689928647831724787741258653D+00
    w(38) = 0.0055045969020008281904902120813D+00
    w(39) = 0.0056393687195659001929970994675D+00
    w(40) = 0.0057732879418203275712033691864D+00
    w(41) = 0.0059063343220074160130475409466D+00
    w(42) = 0.0060384877453327676663371666884D+00
    w(43) = 0.0061697282320052788060812561217D+00
    w(44) = 0.0063000359402577418025981070425D+00
    w(45) = 0.0064293911693465917826140832500D+00
    w(46) = 0.0065577743625303421548456356354D+00
    w(47) = 0.0066851661100262568757892743568D+00
    w(48) = 0.0068115471519448109954345674817D+00
    w(49) = 0.0069368983812014946719507501243D+00
    w(50) = 0.0070612008464055194979848418291D+00
    w(51) = 0.0071844357547249896530757997058D+00
    w(52) = 0.0073065844747281040972736443146D+00
    w(53) = 0.0074276285391999597581348419714D+00
    w(54) = 0.0075475496479345294426435656724D+00
    w(55) = 0.0076663296705013920315933272426D+00
    w(56) = 0.0077839506489867963897419914623D+00
    w(57) = 0.0079003948007086443529587296692D+00
    w(58) = 0.0080156445209049821352946484008D+00
    w(59) = 0.0081296823853955935356080649925D+00
    w(60) = 0.0082424911532162924158504385939D+00
    w(61) = 0.0083540537692255160718568405530D+00
    w(62) = 0.0084643533666828253227353760036D+00
    w(63) = 0.0085733732697989214067758505840D+00
    w(64) = 0.0086810969962567940901133439612D+00
    w(65) = 0.0087875082597036197689825483144D+00
    w(66) = 0.0088925909722130327769834298578D+00
    w(67) = 0.0089963292467173975949700110383D+00
    w(68) = 0.0090987073994097142025303711406D+00
    w(69) = 0.0091997099521147934060534414075D+00
    w(70) = 0.0092993216346293436285393234867D+00
    w(71) = 0.0093975273870306153500305317074D+00
    w(72) = 0.0094943123619532541442165010292D+00
    w(73) = 0.0095896619268340180657610209655D+00
    w(74) = 0.0096835616661240200035669970076D+00
    w(75) = 0.0097759973834681605268499842249D+00
    w(76) = 0.0098669551038514217128483481814D+00
    w(77) = 0.0099564210757116974565448593910D+00
    w(78) = 0.0100443817730188408231888789497D+00
    w(79) = 0.0101308238973196141129538950955D+00
    w(80) = 0.0102157343797482324629939488415D+00
    w(81) = 0.0102991003830021970147153502911D+00
    w(82) = 0.0103809093032831189224876935085D+00
    w(83) = 0.0104611487722022407735015844669D+00
    w(84) = 0.0105398066586503673262517188088D+00
    w(85) = 0.0106168710706319228563864391054D+00
    w(86) = 0.0106923303570628578226139809571D+00
    w(87) = 0.0107661731095321330311788312990D+00
    w(88) = 0.0108383881640265149842990798832D+00
    w(89) = 0.0109089646026184216450603134401D+00
    w(90) = 0.0109778917551165634377595759712D+00
    w(91) = 0.0110451592006791299277436662993D+00
    w(92) = 0.0111107567693892782875426356195D+00
    w(93) = 0.0111746745437926853557086684962D+00
    w(94) = 0.0112369028603969308303734810332D+00
    w(95) = 0.0112974323111324849102690558722D+00
    w(96) = 0.0113562537447750795009464486204D+00
    w(97) = 0.011413358268329247942299599697D+00
    w(98) = 0.011468737248372824084374355981D+00
    w(99) = 0.011522382312362197440930930031D+00
    w(100) = 0.011574285349898127083439539046D+00
    w(101) = 0.011624438513951922901227922331D+00
    w(102) = 0.011672834222051808845465154244D+00
    w(103) = 0.011719465157429288794653489478D+00
    w(104) = 0.011764324270125341726399410909D+00
    w(105) = 0.011807404778056278953532930501D+00
    w(106) = 0.011848700168039102281222824051D+00
    w(107) = 0.011888204196776208064673282076D+00
    w(108) = 0.011925910891799288293359117699D+00
    w(109) = 0.011961814552372285996633285380D+00
    w(110) = 0.011995909750353268455989686823D+00
    w(111) = 0.012028191331015087920350431142D+00
    w(112) = 0.012058654413824705751531083631D+00
    w(113) = 0.012087294393181062176578184854D+00
    w(114) = 0.012114106939111380091025793650D+00
    w(115) = 0.012139087997925797641334635250D+00
    w(116) = 0.012162233792830230614908682534D+00
    w(117) = 0.012183540824497371981177306326D+00
    w(118) = 0.012203005871595742256331865516D+00
    w(119) = 0.012220625991276710706457005806D+00
    w(120) = 0.012236398519619413758040249691D+00
    w(121) = 0.012250321072033503350218104906D+00
    w(122) = 0.012262391543619664338660618398D+00
    w(123) = 0.012272608109487846445745237751D+00
    w(124) = 0.012280969225033162644659793962D+00
    w(125) = 0.012287473626169412265336919908D+00
    w(126) = 0.012292120329520193516690694701D+00
    w(127) = 0.012294908632567576531532225710D+00
    w(128) = 0.01229583811375831445681490730D+00
    w(129) = 0.012294908632567576531532225710D+00
    w(130) = 0.012292120329520193516690694701D+00
    w(131) = 0.012287473626169412265336919908D+00
    w(132) = 0.012280969225033162644659793962D+00
    w(133) = 0.012272608109487846445745237751D+00
    w(134) = 0.012262391543619664338660618398D+00
    w(135) = 0.012250321072033503350218104906D+00
    w(136) = 0.012236398519619413758040249691D+00
    w(137) = 0.012220625991276710706457005806D+00
    w(138) = 0.012203005871595742256331865516D+00
    w(139) = 0.012183540824497371981177306326D+00
    w(140) = 0.012162233792830230614908682534D+00
    w(141) = 0.012139087997925797641334635250D+00
    w(142) = 0.012114106939111380091025793650D+00
    w(143) = 0.012087294393181062176578184854D+00
    w(144) = 0.012058654413824705751531083631D+00
    w(145) = 0.012028191331015087920350431142D+00
    w(146) = 0.011995909750353268455989686823D+00
    w(147) = 0.011961814552372285996633285380D+00
    w(148) = 0.011925910891799288293359117699D+00
    w(149) = 0.011888204196776208064673282076D+00
    w(150) = 0.011848700168039102281222824051D+00
    w(151) = 0.011807404778056278953532930501D+00
    w(152) = 0.011764324270125341726399410909D+00
    w(153) = 0.011719465157429288794653489478D+00
    w(154) = 0.011672834222051808845465154244D+00
    w(155) = 0.011624438513951922901227922331D+00
    w(156) = 0.011574285349898127083439539046D+00
    w(157) = 0.011522382312362197440930930031D+00
    w(158) = 0.011468737248372824084374355981D+00
    w(159) = 0.011413358268329247942299599697D+00
    w(160) = 0.0113562537447750795009464486204D+00
    w(161) = 0.0112974323111324849102690558722D+00
    w(162) = 0.0112369028603969308303734810332D+00
    w(163) = 0.0111746745437926853557086684962D+00
    w(164) = 0.0111107567693892782875426356195D+00
    w(165) = 0.0110451592006791299277436662993D+00
    w(166) = 0.0109778917551165634377595759712D+00
    w(167) = 0.0109089646026184216450603134401D+00
    w(168) = 0.0108383881640265149842990798832D+00
    w(169) = 0.0107661731095321330311788312990D+00
    w(170) = 0.0106923303570628578226139809571D+00
    w(171) = 0.0106168710706319228563864391054D+00
    w(172) = 0.0105398066586503673262517188088D+00
    w(173) = 0.0104611487722022407735015844669D+00
    w(174) = 0.0103809093032831189224876935085D+00
    w(175) = 0.0102991003830021970147153502911D+00
    w(176) = 0.0102157343797482324629939488415D+00
    w(177) = 0.0101308238973196141129538950955D+00
    w(178) = 0.0100443817730188408231888789497D+00
    w(179) = 0.0099564210757116974565448593910D+00
    w(180) = 0.0098669551038514217128483481814D+00
    w(181) = 0.0097759973834681605268499842249D+00
    w(182) = 0.0096835616661240200035669970076D+00
    w(183) = 0.0095896619268340180657610209655D+00
    w(184) = 0.0094943123619532541442165010292D+00
    w(185) = 0.0093975273870306153500305317074D+00
    w(186) = 0.0092993216346293436285393234867D+00
    w(187) = 0.0091997099521147934060534414075D+00
    w(188) = 0.0090987073994097142025303711406D+00
    w(189) = 0.0089963292467173975949700110383D+00
    w(190) = 0.0088925909722130327769834298578D+00
    w(191) = 0.0087875082597036197689825483144D+00
    w(192) = 0.0086810969962567940901133439612D+00
    w(193) = 0.0085733732697989214067758505840D+00
    w(194) = 0.0084643533666828253227353760036D+00
    w(195) = 0.0083540537692255160718568405530D+00
    w(196) = 0.0082424911532162924158504385939D+00
    w(197) = 0.0081296823853955935356080649925D+00
    w(198) = 0.0080156445209049821352946484008D+00
    w(199) = 0.0079003948007086443529587296692D+00
    w(200) = 0.0077839506489867963897419914623D+00
    w(201) = 0.0076663296705013920315933272426D+00
    w(202) = 0.0075475496479345294426435656724D+00
    w(203) = 0.0074276285391999597581348419714D+00
    w(204) = 0.0073065844747281040972736443146D+00
    w(205) = 0.0071844357547249896530757997058D+00
    w(206) = 0.0070612008464055194979848418291D+00
    w(207) = 0.0069368983812014946719507501243D+00
    w(208) = 0.0068115471519448109954345674817D+00
    w(209) = 0.0066851661100262568757892743568D+00
    w(210) = 0.0065577743625303421548456356354D+00
    w(211) = 0.0064293911693465917826140832500D+00
    w(212) = 0.0063000359402577418025981070425D+00
    w(213) = 0.0061697282320052788060812561217D+00
    w(214) = 0.0060384877453327676663371666884D+00
    w(215) = 0.0059063343220074160130475409466D+00
    w(216) = 0.0057732879418203275712033691864D+00
    w(217) = 0.0056393687195659001929970994675D+00
    w(218) = 0.0055045969020008281904902120813D+00
    w(219) = 0.0053689928647831724787741258653D+00
    w(220) = 0.0052325771093919661294970523234D+00
    w(221) = 0.0050953702600278273039420404117D+00
    w(222) = 0.0049573930604950563104281084148D+00
    w(223) = 0.0048186663710656988918572043815D+00
    w(224) = 0.0046792111653260640506279893190D+00
    w(225) = 0.0045390485270061921259394035112D+00
    w(226) = 0.0043981996467927779838546384780D+00
    w(227) = 0.0042566858191260658425395494472D+00
    w(228) = 0.0041145284389812475901826468094D+00
    w(229) = 0.0039717489986349171988699773906D+00
    w(230) = 0.0038283690844171626400743524999D+00
    w(231) = 0.0036844103734499176530742235517D+00
    w(232) = 0.0035398946303722552150296713510D+00
    w(233) = 0.0033948437040533928255056951665D+00
    w(234) = 0.0032492795242943133198690930777D+00
    w(235) = 0.0031032240985191112621977893133D+00
    w(236) = 0.0029566995084575002760043344138D+00
    w(237) = 0.0028097279068204407457332299361D+00
    w(238) = 0.0026623315139717112732749157331D+00
    w(239) = 0.0025145326145997073931298921370D+00
    w(240) = 0.0023663535543962867157201855305D+00
    w(241) = 0.0022178167367540171700373764020D+00
    w(242) = 0.0020689446195015801533643667413D+00
    w(243) = 0.0019197597117132050055085980675D+00
    w(244) = 0.0017702845706603213070421243905D+00
    w(245) = 0.0016205417990415653896921100325D+00
    w(246) = 0.0014705540427783843160097204304D+00
    w(247) = 0.0013203439900221692090523602144D+00
    w(248) = 0.0011699343729388079886897709773D+00
    w(249) = 0.00101934797642732530281229369360D+00
    w(250) = 0.00086860766611945667949717690640D+00
    w(251) = 0.00071773647800611087798371518325D+00
    w(252) = 0.00056675794564824918946626058353D+00
    w(253) = 0.00041569762526823913616284210066D+00
    w(254) = 0.00026459387119083065532790838855D+00
    w(255) = 0.00011367361999142272115645954414D+00

  else if ( n == 256 ) then

    x(1) = -0.999956050018992230734801D+00
    x(2) = -0.999768437409263186104879D+00
    x(3) = -0.999430937466261408240854D+00
    x(4) = -0.998943525843408856555026D+00
    x(5) = -0.998306266473006444055500D+00
    x(6) = -0.997519252756720827563409D+00
    x(7) = -0.996582602023381540430504D+00
    x(8) = -0.995496454481096356592647D+00
    x(9) = -0.994260972922409664962878D+00
    x(10) = -0.992876342608822117143534D+00
    x(11) = -0.991342771207583086922189D+00
    x(12) = -0.989660488745065218319244D+00
    x(13) = -0.987829747564860608916488D+00
    x(14) = -0.985850822286125956479245D+00
    x(15) = -0.983724009760315496166686D+00
    x(16) = -0.981449629025464405769303D+00
    x(17) = -0.979028021257622038824238D+00
    x(18) = -0.976459549719234155621011D+00
    x(19) = -0.973744599704370405266079D+00
    x(20) = -0.970883578480743029320923D+00
    x(21) = -0.967876915228489454909004D+00
    x(22) = -0.964725060975706430932612D+00
    x(23) = -0.961428488530732144006407D+00
    x(24) = -0.957987692411178129365790D+00
    x(25) = -0.954403188769716241764448D+00
    x(26) = -0.950675515316628276363852D+00
    x(27) = -0.946805231239127481372052D+00
    x(28) = -0.942792917117462443183076D+00
    x(29) = -0.938639174837814804981926D+00
    x(30) = -0.934344627502003094292477D+00
    x(31) = -0.929909919334005641180246D+00
    x(32) = -0.925335715583316202872730D+00
    x(33) = -0.920622702425146495505047D+00
    x(34) = -0.915771586857490384526670D+00
    x(35) = -0.910783096595065011890907D+00
    x(36) = -0.905657979960144647082682D+00
    x(37) = -0.900397005770303544771620D+00
    x(38) = -0.895000963223084577441223D+00
    x(39) = -0.889470661777610888828677D+00
    x(40) = -0.883806931033158284859826D+00
    x(41) = -0.878010620604706543986435D+00
    x(42) = -0.872082599995488289130046D+00
    x(43) = -0.866023758466554519297515D+00
    x(44) = -0.859835004903376350696173D+00
    x(45) = -0.853517267679502965073036D+00
    x(46) = -0.847071494517296207187072D+00
    x(47) = -0.840498652345762713895068D+00
    x(48) = -0.833799727155504894348444D+00
    x(49) = -0.826975723850812514289093D+00
    x(50) = -0.820027666098917067403478D+00
    x(51) = -0.812956596176431543136410D+00
    x(52) = -0.805763574812998623257389D+00
    x(53) = -0.798449681032170758782543D+00
    x(54) = -0.791016011989545994546707D+00
    x(55) = -0.783463682808183820750670D+00
    x(56) = -0.775793826411325739132053D+00
    x(57) = -0.768007593352445635975891D+00
    x(58) = -0.760106151642655454941907D+00
    x(59) = -0.752090686575492059587530D+00
    x(60) = -0.743962400549111568455683D+00
    x(61) = -0.735722512885917834620373D+00
    x(62) = -0.727372259649652126586894D+00
    x(63) = -0.718912893459971448372640D+00
    x(64) = -0.710345683304543313394566D+00
    x(65) = -0.701671914348685159406084D+00
    x(66) = -0.692892887742576960105342D+00
    x(67) = -0.684009920426075953124877D+00
    x(68) = -0.675024344931162763855919D+00
    x(69) = -0.665937509182048559906408D+00
    x(70) = -0.656750776292973221887500D+00
    x(71) = -0.647465524363724862617016D+00
    x(72) = -0.638083146272911368668689D+00
    x(73) = -0.628605049469014975432210D+00
    x(74) = -0.619032655759261219430968D+00
    x(75) = -0.609367401096333939522311D+00
    x(76) = -0.599610735362968321730388D+00
    x(77) = -0.589764122154454300785786D+00
    x(78) = -0.579829038559082944921832D+00
    x(79) = -0.569806974936568759057668D+00
    x(80) = -0.559699434694481145136907D+00
    x(81) = -0.549507934062718557042427D+00
    x(82) = -0.539234001866059181127936D+00
    x(83) = -0.528879179294822261951476D+00
    x(84) = -0.518445019673674476221662D+00
    x(85) = -0.507933088228616036231925D+00
    x(86) = -0.497344961852181477119512D+00
    x(87) = -0.486682228866890350103621D+00
    x(88) = -0.475946488786983306390738D+00
    x(89) = -0.465139352078479313645570D+00
    x(90) = -0.454262439917589998774455D+00
    x(91) = -0.443317383947527357216926D+00
    x(92) = -0.432305826033741309953441D+00
    x(93) = -0.421229418017623824976812D+00
    x(94) = -0.410089821468716550006434D+00
    x(95) = -0.398888707435459127713463D+00
    x(96) = -0.387627756194515583637985D+00
    x(97) = -0.376308656998716390283056D+00
    x(98) = -0.364933107823654018533465D+00
    x(99) = -0.353502815112969989537790D+00
    x(100) = -0.342019493522371636480730D+00
    x(101) = -0.330484865662416976229187D+00
    x(102) = -0.318900661840106275631683D+00
    x(103) = -0.307268619799319076258610D+00
    x(104) = -0.295590484460135614563787D+00
    x(105) = -0.283868007657081741799766D+00
    x(106) = -0.272102947876336609505245D+00
    x(107) = -0.260297069991942541978561D+00
    x(108) = -0.248452145001056666833243D+00
    x(109) = -0.236569949758284018477508D+00
    x(110) = -0.224652266709131967147878D+00
    x(111) = -0.212700883622625957937040D+00
    x(112) = -0.200717593323126670068001D+00
    x(113) = -0.188704193421388826461504D+00
    x(114) = -0.176662486044901997403722D+00
    x(115) = -0.164594277567553849829285D+00
    x(116) = -0.152501378338656395374607D+00
    x(117) = -0.140385602411375885913025D+00
    x(118) = -0.128248767270607094742050D+00
    x(119) = -0.116092693560332804940735D+00
    x(120) = -0.103919204810509403639197D+00
    x(121) = -0.091730127163519552031146D+00
    x(122) = -0.079527289100232965903227D+00
    x(123) = -0.067312521165716400242290D+00
    x(124) = -0.055087655694633984104561D+00
    x(125) = -0.042854526536379098381242D+00
    x(126) = -0.030614968779979029366279D+00
    x(127) = -0.018370818478813665117926D+00
    x(128) = -0.006123912375189529501170D+00
    x(129) = 0.006123912375189529501170D+00
    x(130) = 0.018370818478813665117926D+00
    x(131) = 0.030614968779979029366279D+00
    x(132) = 0.042854526536379098381242D+00
    x(133) = 0.055087655694633984104561D+00
    x(134) = 0.067312521165716400242290D+00
    x(135) = 0.079527289100232965903227D+00
    x(136) = 0.091730127163519552031146D+00
    x(137) = 0.103919204810509403639197D+00
    x(138) = 0.116092693560332804940735D+00
    x(139) = 0.128248767270607094742050D+00
    x(140) = 0.140385602411375885913025D+00
    x(141) = 0.152501378338656395374607D+00
    x(142) = 0.164594277567553849829285D+00
    x(143) = 0.176662486044901997403722D+00
    x(144) = 0.188704193421388826461504D+00
    x(145) = 0.200717593323126670068001D+00
    x(146) = 0.212700883622625957937040D+00
    x(147) = 0.224652266709131967147878D+00
    x(148) = 0.236569949758284018477508D+00
    x(149) = 0.248452145001056666833243D+00
    x(150) = 0.260297069991942541978561D+00
    x(151) = 0.272102947876336609505245D+00
    x(152) = 0.283868007657081741799766D+00
    x(153) = 0.295590484460135614563787D+00
    x(154) = 0.307268619799319076258610D+00
    x(155) = 0.318900661840106275631683D+00
    x(156) = 0.330484865662416976229187D+00
    x(157) = 0.342019493522371636480730D+00
    x(158) = 0.353502815112969989537790D+00
    x(159) = 0.364933107823654018533465D+00
    x(160) = 0.376308656998716390283056D+00
    x(161) = 0.387627756194515583637985D+00
    x(162) = 0.398888707435459127713463D+00
    x(163) = 0.410089821468716550006434D+00
    x(164) = 0.421229418017623824976812D+00
    x(165) = 0.432305826033741309953441D+00
    x(166) = 0.443317383947527357216926D+00
    x(167) = 0.454262439917589998774455D+00
    x(168) = 0.465139352078479313645570D+00
    x(169) = 0.475946488786983306390738D+00
    x(170) = 0.486682228866890350103621D+00
    x(171) = 0.497344961852181477119512D+00
    x(172) = 0.507933088228616036231925D+00
    x(173) = 0.518445019673674476221662D+00
    x(174) = 0.528879179294822261951476D+00
    x(175) = 0.539234001866059181127936D+00
    x(176) = 0.549507934062718557042427D+00
    x(177) = 0.559699434694481145136907D+00
    x(178) = 0.569806974936568759057668D+00
    x(179) = 0.579829038559082944921832D+00
    x(180) = 0.589764122154454300785786D+00
    x(181) = 0.599610735362968321730388D+00
    x(182) = 0.609367401096333939522311D+00
    x(183) = 0.619032655759261219430968D+00
    x(184) = 0.628605049469014975432210D+00
    x(185) = 0.638083146272911368668689D+00
    x(186) = 0.647465524363724862617016D+00
    x(187) = 0.656750776292973221887500D+00
    x(188) = 0.665937509182048559906408D+00
    x(189) = 0.675024344931162763855919D+00
    x(190) = 0.684009920426075953124877D+00
    x(191) = 0.692892887742576960105342D+00
    x(192) = 0.701671914348685159406084D+00
    x(193) = 0.710345683304543313394566D+00
    x(194) = 0.718912893459971448372640D+00
    x(195) = 0.727372259649652126586894D+00
    x(196) = 0.735722512885917834620373D+00
    x(197) = 0.743962400549111568455683D+00
    x(198) = 0.752090686575492059587530D+00
    x(199) = 0.760106151642655454941907D+00
    x(200) = 0.768007593352445635975891D+00
    x(201) = 0.775793826411325739132053D+00
    x(202) = 0.783463682808183820750670D+00
    x(203) = 0.791016011989545994546707D+00
    x(204) = 0.798449681032170758782543D+00
    x(205) = 0.805763574812998623257389D+00
    x(206) = 0.812956596176431543136410D+00
    x(207) = 0.820027666098917067403478D+00
    x(208) = 0.826975723850812514289093D+00
    x(209) = 0.833799727155504894348444D+00
    x(210) = 0.840498652345762713895068D+00
    x(211) = 0.847071494517296207187072D+00
    x(212) = 0.853517267679502965073036D+00
    x(213) = 0.859835004903376350696173D+00
    x(214) = 0.866023758466554519297515D+00
    x(215) = 0.872082599995488289130046D+00
    x(216) = 0.878010620604706543986435D+00
    x(217) = 0.883806931033158284859826D+00
    x(218) = 0.889470661777610888828677D+00
    x(219) = 0.895000963223084577441223D+00
    x(220) = 0.900397005770303544771620D+00
    x(221) = 0.905657979960144647082682D+00
    x(222) = 0.910783096595065011890907D+00
    x(223) = 0.915771586857490384526670D+00
    x(224) = 0.920622702425146495505047D+00
    x(225) = 0.925335715583316202872730D+00
    x(226) = 0.929909919334005641180246D+00
    x(227) = 0.934344627502003094292477D+00
    x(228) = 0.938639174837814804981926D+00
    x(229) = 0.942792917117462443183076D+00
    x(230) = 0.946805231239127481372052D+00
    x(231) = 0.950675515316628276363852D+00
    x(232) = 0.954403188769716241764448D+00
    x(233) = 0.957987692411178129365790D+00
    x(234) = 0.961428488530732144006407D+00
    x(235) = 0.964725060975706430932612D+00
    x(236) = 0.967876915228489454909004D+00
    x(237) = 0.970883578480743029320923D+00
    x(238) = 0.973744599704370405266079D+00
    x(239) = 0.976459549719234155621011D+00
    x(240) = 0.979028021257622038824238D+00
    x(241) = 0.981449629025464405769303D+00
    x(242) = 0.983724009760315496166686D+00
    x(243) = 0.985850822286125956479245D+00
    x(244) = 0.987829747564860608916488D+00
    x(245) = 0.989660488745065218319244D+00
    x(246) = 0.991342771207583086922189D+00
    x(247) = 0.992876342608822117143534D+00
    x(248) = 0.994260972922409664962878D+00
    x(249) = 0.995496454481096356592647D+00
    x(250) = 0.996582602023381540430504D+00
    x(251) = 0.997519252756720827563409D+00
    x(252) = 0.998306266473006444055500D+00
    x(253) = 0.998943525843408856555026D+00
    x(254) = 0.999430937466261408240854D+00
    x(255) = 0.999768437409263186104879D+00
    x(256) = 0.999956050018992230734801D+00

    w(1) = 0.00011278901782227217551253887725D+00
    w(2) = 0.00026253494429644590628745756250D+00
    w(3) = 0.00041246325442617632843218583774D+00
    w(4) = 0.00056234895403140980281523674759D+00
    w(5) = 0.0007121541634733206669089891511D+00
    w(6) = 0.0008618537014200890378140934163D+00
    w(7) = 0.0010114243932084404526058128414D+00
    w(8) = 0.0011608435575677247239705981135D+00
    w(9) = 0.0013100886819025044578316804271D+00
    w(10) = 0.0014591373333107332010883864996D+00
    w(11) = 0.0016079671307493272424499395690D+00
    w(12) = 0.0017565557363307299936069145295D+00
    w(13) = 0.0019048808534997184044191411746D+00
    w(14) = 0.0020529202279661431745487818492D+00
    w(15) = 0.0022006516498399104996848834189D+00
    w(16) = 0.0023480529563273120170064609087D+00
    w(17) = 0.0024951020347037068508395354372D+00
    w(18) = 0.0026417768254274905641208292516D+00
    w(19) = 0.0027880553253277068805747610763D+00
    w(20) = 0.0029339155908297166460123254142D+00
    w(21) = 0.0030793357411993375832053528316D+00
    w(22) = 0.0032242939617941981570107134269D+00
    w(23) = 0.0033687685073155510120191062489D+00
    w(24) = 0.0035127377050563073309710549844D+00
    w(25) = 0.0036561799581425021693892413052D+00
    w(26) = 0.0037990737487662579981170192082D+00
    w(27) = 0.0039413976414088336277290349840D+00
    w(28) = 0.0040831302860526684085997759212D+00
    w(29) = 0.0042242504213815362723565049060D+00
    w(30) = 0.0043647368779680566815684200621D+00
    w(31) = 0.0045045685814478970686417923159D+00
    w(32) = 0.0046437245556800603139790923525D+00
    w(33) = 0.0047821839258926913729317340448D+00
    w(34) = 0.0049199259218138656695587765655D+00
    w(35) = 0.0050569298807868423875578160762D+00
    w(36) = 0.0051931752508692809303287536296D+00
    w(37) = 0.0053286415939159303170811114788D+00
    w(38) = 0.0054633085886443102775705318566D+00
    w(39) = 0.0055971560336829100775514452572D+00
    w(40) = 0.005730163850601437177384417555D+00
    w(41) = 0.005862312086922653060661598801D+00
    w(42) = 0.005993580919115338221127696870D+00
    w(43) = 0.006123950655567932542389081187D+00
    w(44) = 0.006253401739542401272063645975D+00
    w(45) = 0.006381914752107880570375164275D+00
    w(46) = 0.006509470415053660267809899951D+00
    w(47) = 0.006636049593781065044590038355D+00
    w(48) = 0.006761633300173798780927861108D+00
    w(49) = 0.006886202695446320346713323775D+00
    w(50) = 0.007009739092969822621234436194D+00
    w(51) = 0.007132223961075390071672422986D+00
    w(52) = 0.007253638925833913783829137214D+00
    w(53) = 0.007373965773812346437572440695D+00
    w(54) = 0.007493186454805883358599761133D+00
    w(55) = 0.007611283084545659461618719618D+00
    w(56) = 0.007728237947381555631110194958D+00
    w(57) = 0.007844033498939711866810316151D+00
    w(58) = 0.007958652368754348353613161227D+00
    w(59) = 0.008072077362873499500946974804D+00
    w(60) = 0.008184291466438269935619761004D+00
    w(61) = 0.008295277846235225425171412553D+00
    w(62) = 0.008405019853221535756180301698D+00
    w(63) = 0.008513501025022490693838354790D+00
    w(64) = 0.008620705088401014305368838410D+00
    w(65) = 0.008726615961698807140336632217D+00
    w(66) = 0.008831217757248750025318272685D+00
    w(67) = 0.008934494783758207548408417085D+00
    w(68) = 0.009036431548662873680227775572D+00
    w(69) = 0.009137012760450806402000472219D+00
    w(70) = 0.009236223330956302687378716714D+00
    w(71) = 0.009334048377623269712466014486D+00
    w(72) = 0.009430473225737752747352764482D+00
    w(73) = 0.009525483410629284811829685754D+00
    w(74) = 0.009619064679840727857162164401D+00
    w(75) = 0.009711202995266279964249670496D+00
    w(76) = 0.009801884535257327825498800250D+00
    w(77) = 0.009891095696695828602630683809D+00
    w(78) = 0.009978823097034910124733949495D+00
    w(79) = 0.010065053576306383309460978930D+00
    w(80) = 0.010149774199094865654634066042D+00
    w(81) = 0.010232972256478219656954857160D+00
    w(82) = 0.010314635267934015068260713997D+00
    w(83) = 0.010394750983211728997101725205D+00
    w(84) = 0.010473307384170403003569566927D+00
    w(85) = 0.010550292686581481517533575536D+00
    w(86) = 0.010625695341896561133961681801D+00
    w(87) = 0.010699504038979785603048200583D+00
    w(88) = 0.010771707705804626636653631927D+00
    w(89) = 0.010842295511114795995293477058D+00
    w(90) = 0.010911256866049039700796847788D+00
    w(91) = 0.010978581425729570637988203448D+00
    w(92) = 0.011044259090813901263517571044D+00
    w(93) = 0.011108280009009843630460815451D+00
    w(94) = 0.011170634576553449462710881938D+00
    w(95) = 0.011231313439649668572656802083D+00
    w(96) = 0.011290307495875509508367594121D+00
    w(97) = 0.011347607895545491941625714297D+00
    w(98) = 0.011403206043039185964847059552D+00
    w(99) = 0.011457093598090639152334392298D+00
    w(100) = 0.011509262477039497958586392439D+00
    w(101) = 0.011559704854043635772668656950D+00
    w(102) = 0.011608413162253105722084706677D+00
    w(103) = 0.011655380094945242121298939730D+00
    w(104) = 0.011700598606620740288189823359D+00
    w(105) = 0.011744061914060550305376732759D+00
    w(106) = 0.011785763497343426181690117627D+00
    w(107) = 0.011825697100823977771160737958D+00
    w(108) = 0.011863856734071078731904572908D+00
    w(109) = 0.011900236672766489754287204237D+00
    w(110) = 0.011934831459563562255873201696D+00
    w(111) = 0.011967635904905893729007282670D+00
    w(112) = 0.011998645087805811934536710071D+00
    w(113) = 0.012027854356582571161267533498D+00
    w(114) = 0.012055259329560149814347085327D+00
    w(115) = 0.012080855895724544655975183976D+00
    w(116) = 0.012104640215340463097757829736D+00
    w(117) = 0.012126608720527321034718492205D+00
    w(118) = 0.012146758115794459815559837664D+00
    w(119) = 0.012165085378535502061307291839D+00
    w(120) = 0.012181587759481772174047585032D+00
    w(121) = 0.012196262783114713518180974196D+00
    w(122) = 0.012209108248037240407514094371D+00
    w(123) = 0.012220122227303969191708737227D+00
    w(124) = 0.012229303068710278904146266083D+00
    w(125) = 0.012236649395040158109242574767D+00
    w(126) = 0.012242160104272800769728083260D+00
    w(127) = 0.012245834369747920142463857550D+00
    w(128) = 0.01224767164028975590407032649D+00
    w(129) = 0.01224767164028975590407032649D+00
    w(130) = 0.012245834369747920142463857550D+00
    w(131) = 0.012242160104272800769728083260D+00
    w(132) = 0.012236649395040158109242574767D+00
    w(133) = 0.012229303068710278904146266083D+00
    w(134) = 0.012220122227303969191708737227D+00
    w(135) = 0.012209108248037240407514094371D+00
    w(136) = 0.012196262783114713518180974196D+00
    w(137) = 0.012181587759481772174047585032D+00
    w(138) = 0.012165085378535502061307291839D+00
    w(139) = 0.012146758115794459815559837664D+00
    w(140) = 0.012126608720527321034718492205D+00
    w(141) = 0.012104640215340463097757829736D+00
    w(142) = 0.012080855895724544655975183976D+00
    w(143) = 0.012055259329560149814347085327D+00
    w(144) = 0.012027854356582571161267533498D+00
    w(145) = 0.011998645087805811934536710071D+00
    w(146) = 0.011967635904905893729007282670D+00
    w(147) = 0.011934831459563562255873201696D+00
    w(148) = 0.011900236672766489754287204237D+00
    w(149) = 0.011863856734071078731904572908D+00
    w(150) = 0.011825697100823977771160737958D+00
    w(151) = 0.011785763497343426181690117627D+00
    w(152) = 0.011744061914060550305376732759D+00
    w(153) = 0.011700598606620740288189823359D+00
    w(154) = 0.011655380094945242121298939730D+00
    w(155) = 0.011608413162253105722084706677D+00
    w(156) = 0.011559704854043635772668656950D+00
    w(157) = 0.011509262477039497958586392439D+00
    w(158) = 0.011457093598090639152334392298D+00
    w(159) = 0.011403206043039185964847059552D+00
    w(160) = 0.011347607895545491941625714297D+00
    w(161) = 0.011290307495875509508367594121D+00
    w(162) = 0.011231313439649668572656802083D+00
    w(163) = 0.011170634576553449462710881938D+00
    w(164) = 0.011108280009009843630460815451D+00
    w(165) = 0.011044259090813901263517571044D+00
    w(166) = 0.010978581425729570637988203448D+00
    w(167) = 0.010911256866049039700796847788D+00
    w(168) = 0.010842295511114795995293477058D+00
    w(169) = 0.010771707705804626636653631927D+00
    w(170) = 0.010699504038979785603048200583D+00
    w(171) = 0.010625695341896561133961681801D+00
    w(172) = 0.010550292686581481517533575536D+00
    w(173) = 0.010473307384170403003569566927D+00
    w(174) = 0.010394750983211728997101725205D+00
    w(175) = 0.010314635267934015068260713997D+00
    w(176) = 0.010232972256478219656954857160D+00
    w(177) = 0.010149774199094865654634066042D+00
    w(178) = 0.010065053576306383309460978930D+00
    w(179) = 0.009978823097034910124733949495D+00
    w(180) = 0.009891095696695828602630683809D+00
    w(181) = 0.009801884535257327825498800250D+00
    w(182) = 0.009711202995266279964249670496D+00
    w(183) = 0.009619064679840727857162164401D+00
    w(184) = 0.009525483410629284811829685754D+00
    w(185) = 0.009430473225737752747352764482D+00
    w(186) = 0.009334048377623269712466014486D+00
    w(187) = 0.009236223330956302687378716714D+00
    w(188) = 0.009137012760450806402000472219D+00
    w(189) = 0.009036431548662873680227775572D+00
    w(190) = 0.008934494783758207548408417085D+00
    w(191) = 0.008831217757248750025318272685D+00
    w(192) = 0.008726615961698807140336632217D+00
    w(193) = 0.008620705088401014305368838410D+00
    w(194) = 0.008513501025022490693838354790D+00
    w(195) = 0.008405019853221535756180301698D+00
    w(196) = 0.008295277846235225425171412553D+00
    w(197) = 0.008184291466438269935619761004D+00
    w(198) = 0.008072077362873499500946974804D+00
    w(199) = 0.007958652368754348353613161227D+00
    w(200) = 0.007844033498939711866810316151D+00
    w(201) = 0.007728237947381555631110194958D+00
    w(202) = 0.007611283084545659461618719618D+00
    w(203) = 0.007493186454805883358599761133D+00
    w(204) = 0.007373965773812346437572440695D+00
    w(205) = 0.007253638925833913783829137214D+00
    w(206) = 0.007132223961075390071672422986D+00
    w(207) = 0.007009739092969822621234436194D+00
    w(208) = 0.006886202695446320346713323775D+00
    w(209) = 0.006761633300173798780927861108D+00
    w(210) = 0.006636049593781065044590038355D+00
    w(211) = 0.006509470415053660267809899951D+00
    w(212) = 0.006381914752107880570375164275D+00
    w(213) = 0.006253401739542401272063645975D+00
    w(214) = 0.006123950655567932542389081187D+00
    w(215) = 0.005993580919115338221127696870D+00
    w(216) = 0.005862312086922653060661598801D+00
    w(217) = 0.005730163850601437177384417555D+00
    w(218) = 0.0055971560336829100775514452572D+00
    w(219) = 0.0054633085886443102775705318566D+00
    w(220) = 0.0053286415939159303170811114788D+00
    w(221) = 0.0051931752508692809303287536296D+00
    w(222) = 0.0050569298807868423875578160762D+00
    w(223) = 0.0049199259218138656695587765655D+00
    w(224) = 0.0047821839258926913729317340448D+00
    w(225) = 0.0046437245556800603139790923525D+00
    w(226) = 0.0045045685814478970686417923159D+00
    w(227) = 0.0043647368779680566815684200621D+00
    w(228) = 0.0042242504213815362723565049060D+00
    w(229) = 0.0040831302860526684085997759212D+00
    w(230) = 0.0039413976414088336277290349840D+00
    w(231) = 0.0037990737487662579981170192082D+00
    w(232) = 0.0036561799581425021693892413052D+00
    w(233) = 0.0035127377050563073309710549844D+00
    w(234) = 0.0033687685073155510120191062489D+00
    w(235) = 0.0032242939617941981570107134269D+00
    w(236) = 0.0030793357411993375832053528316D+00
    w(237) = 0.0029339155908297166460123254142D+00
    w(238) = 0.0027880553253277068805747610763D+00
    w(239) = 0.0026417768254274905641208292516D+00
    w(240) = 0.0024951020347037068508395354372D+00
    w(241) = 0.0023480529563273120170064609087D+00
    w(242) = 0.0022006516498399104996848834189D+00
    w(243) = 0.0020529202279661431745487818492D+00
    w(244) = 0.0019048808534997184044191411746D+00
    w(245) = 0.0017565557363307299936069145295D+00
    w(246) = 0.0016079671307493272424499395690D+00
    w(247) = 0.0014591373333107332010883864996D+00
    w(248) = 0.0013100886819025044578316804271D+00
    w(249) = 0.0011608435575677247239705981135D+00
    w(250) = 0.0010114243932084404526058128414D+00
    w(251) = 0.0008618537014200890378140934163D+00
    w(252) = 0.0007121541634733206669089891511D+00
    w(253) = 0.00056234895403140980281523674759D+00
    w(254) = 0.00041246325442617632843218583774D+00
    w(255) = 0.00026253494429644590628745756250D+00
    w(256) = 0.00011278901782227217551253887725D+00

  else if ( n == 257 ) then

    x(1) = -0.999956390712330402472857D+00
    x(2) = -0.999770232390338019056053D+00
    x(3) = -0.999435348366365078441838D+00
    x(4) = -0.998951714093223210129834D+00
    x(5) = -0.998319392445383847808766D+00
    x(6) = -0.997538475365520218731818D+00
    x(7) = -0.996609078365487004512326D+00
    x(8) = -0.995531339486830143483750D+00
    x(9) = -0.994305419008553630362377D+00
    x(10) = -0.992931499332908653172844D+00
    x(11) = -0.991409784923101705201254D+00
    x(12) = -0.989740502257507526030375D+00
    x(13) = -0.987923899788618253106809D+00
    x(14) = -0.985960247902290665366669D+00
    x(15) = -0.983849838875444644048531D+00
    x(16) = -0.981592986831381877693095D+00
    x(17) = -0.979190027692327124191591D+00
    x(18) = -0.976641319128992592610888D+00
    x(19) = -0.973947240507062326750976D+00
    x(20) = -0.971108192830542793021113D+00
    x(21) = -0.968124598681952354372943D+00
    x(22) = -0.964996902159337170373447D+00
    x(23) = -0.961725568810109767190665D+00
    x(24) = -0.958311085561711847074814D+00
    x(25) = -0.954753960649106318830855D+00
    x(26) = -0.951054723539105826691801D+00
    x(27) = -0.947213924851546682950881D+00
    x(28) = -0.943232136277318328151464D+00
    x(29) = -0.939109950493259404355123D+00
    x(30) = -0.934847981073932324370129D+00
    x(31) = -0.930446862400288909805510D+00
    x(32) = -0.925907249565240289235888D+00
    x(33) = -0.921229818276144817520964D+00
    x(34) = -0.916415264754228313295468D+00
    x(35) = -0.911464305630951423630955D+00
    x(36) = -0.906377677841339419411308D+00
    x(37) = -0.901156138514290206476301D+00
    x(38) = -0.895800464859876809085345D+00
    x(39) = -0.890311454053661045810287D+00
    x(40) = -0.884689923118035575018750D+00
    x(41) = -0.878936708800611938658765D+00
    x(42) = -0.873052667449672679799858D+00
    x(43) = -0.867038674886706051812473D+00
    x(44) = -0.860895626276042275514686D+00
    x(45) = -0.854624435991610735314055D+00
    x(46) = -0.848226037480837936478636D+00
    x(47) = -0.841701383125706473284556D+00
    x(48) = -0.835051444100995681967937D+00
    x(49) = -0.828277210229725073186687D+00
    x(50) = -0.821379689835822056081139D+00
    x(51) = -0.814359909594035880004229D+00
    x(52) = -0.807218914377120130552073D+00
    x(53) = -0.799957767100306523636066D+00
    x(54) = -0.792577548563093144962574D+00
    x(55) = -0.785079357288370682385816D+00
    x(56) = -0.777464309358910595129671D+00
    x(57) = -0.769733538251239556788216D+00
    x(58) = -0.761888194666924898264210D+00
    x(59) = -0.753929446361296162339238D+00
    x(60) = -0.745858477969628263337895D+00
    x(61) = -0.737676490830812123299244D+00
    x(62) = -0.729384702808539030149808D+00
    x(63) = -0.720984348110025333531072D+00
    x(64) = -0.712476677102304460118510D+00
    x(65) = -0.703862956126113592426171D+00
    x(66) = -0.695144467307402713168813D+00
    x(67) = -0.686322508366494071200553D+00
    x(68) = -0.677398392424920474813593D+00
    x(69) = -0.668373447809971163711735D+00
    x(70) = -0.659249017856974352220492D+00
    x(71) = -0.650026460709345873208532D+00
    x(72) = -0.640707149116433684724434D+00
    x(73) = -0.631292470229188329449219D+00
    x(74) = -0.621783825393689760680446D+00
    x(75) = -0.612182629942561267650033D+00
    x(76) = -0.602490312984301547488097D+00
    x(77) = -0.592708317190566281032495D+00
    x(78) = -0.582838098581430874902446D+00
    x(79) = -0.572881126308666332759406D+00
    x(80) = -0.562838882437060514424546D+00
    x(81) = -0.552712861723817332466074D+00
    x(82) = -0.542504571396066721967792D+00
    x(83) = -0.532215530926518500400434D+00
    x(84) = -0.521847271807293510797499D+00
    x(85) = -0.511401337321965712746629D+00
    x(86) = -0.500879282315849152005553D+00
    x(87) = -0.490282672964564000798817D+00
    x(88) = -0.479613086540916117008992D+00
    x(89) = -0.468872111180124821505728D+00
    x(90) = -0.458061345643433838720630D+00
    x(91) = -0.447182399080140586238810D+00
    x(92) = -0.436236890788079234603398D+00
    x(93) = -0.425226449972593188682213D+00
    x(94) = -0.414152715504032866791986D+00
    x(95) = -0.403017335673814873281489D+00
    x(96) = -0.391821967949078874408131D+00
    x(97) = -0.380568278725978696070941D+00
    x(98) = -0.369257943081644365255611D+00
    x(99) = -0.357892644524852014873858D+00
    x(100) = -0.346474074745438764010632D+00
    x(101) = -0.335003933362499872399782D+00
    x(102) = -0.323483927671405649204085D+00
    x(103) = -0.311915772389675771851948D+00
    x(104) = -0.300301189401748840754520D+00
    x(105) = -0.288641907502685160168097D+00
    x(106) = -0.276939662140840894253032D+00
    x(107) = -0.265196195159551900488370D+00
    x(108) = -0.253413254537865690008131D+00
    x(109) = -0.241592594130360106108882D+00
    x(110) = -0.229735973406087448117604D+00
    x(111) = -0.217845157186682897983880D+00
    x(112) = -0.205921915383676231351599D+00
    x(113) = -0.193968022735045913454182D+00
    x(114) = -0.181985258541054792946197D+00
    x(115) = -0.169975406399406713716337D+00
    x(116) = -0.157940253939763465806087D+00
    x(117) = -0.145881592557661591770148D+00
    x(118) = -0.133801217147868654144405D+00
    x(119) = -0.121700925837218653121859D+00
    x(120) = -0.109582519716966361063898D+00
    x(121) = -0.097447802574700412082119D+00
    x(122) = -0.085298580625855050603929D+00
    x(123) = -0.073136662244860502573600D+00
    x(124) = -0.060963857695971986730406D+00
    x(125) = -0.048781978863817431238958D+00
    x(126) = -0.036592838983704002816750D+00
    x(127) = -0.024398252371723591403953D+00
    x(128) = -0.012200034154697423345412D+00
    x(129) = 0.000000000000000000000000D+00
    x(130) = 0.012200034154697423345412D+00
    x(131) = 0.024398252371723591403953D+00
    x(132) = 0.036592838983704002816750D+00
    x(133) = 0.048781978863817431238958D+00
    x(134) = 0.060963857695971986730406D+00
    x(135) = 0.073136662244860502573600D+00
    x(136) = 0.085298580625855050603929D+00
    x(137) = 0.097447802574700412082119D+00
    x(138) = 0.109582519716966361063898D+00
    x(139) = 0.121700925837218653121859D+00
    x(140) = 0.133801217147868654144405D+00
    x(141) = 0.145881592557661591770148D+00
    x(142) = 0.157940253939763465806087D+00
    x(143) = 0.169975406399406713716337D+00
    x(144) = 0.181985258541054792946197D+00
    x(145) = 0.193968022735045913454182D+00
    x(146) = 0.205921915383676231351599D+00
    x(147) = 0.217845157186682897983880D+00
    x(148) = 0.229735973406087448117604D+00
    x(149) = 0.241592594130360106108882D+00
    x(150) = 0.253413254537865690008131D+00
    x(151) = 0.265196195159551900488370D+00
    x(152) = 0.276939662140840894253032D+00
    x(153) = 0.288641907502685160168097D+00
    x(154) = 0.300301189401748840754520D+00
    x(155) = 0.311915772389675771851948D+00
    x(156) = 0.323483927671405649204085D+00
    x(157) = 0.335003933362499872399782D+00
    x(158) = 0.346474074745438764010632D+00
    x(159) = 0.357892644524852014873858D+00
    x(160) = 0.369257943081644365255611D+00
    x(161) = 0.380568278725978696070941D+00
    x(162) = 0.391821967949078874408131D+00
    x(163) = 0.403017335673814873281489D+00
    x(164) = 0.414152715504032866791986D+00
    x(165) = 0.425226449972593188682213D+00
    x(166) = 0.436236890788079234603398D+00
    x(167) = 0.447182399080140586238810D+00
    x(168) = 0.458061345643433838720630D+00
    x(169) = 0.468872111180124821505728D+00
    x(170) = 0.479613086540916117008992D+00
    x(171) = 0.490282672964564000798817D+00
    x(172) = 0.500879282315849152005553D+00
    x(173) = 0.511401337321965712746629D+00
    x(174) = 0.521847271807293510797499D+00
    x(175) = 0.532215530926518500400434D+00
    x(176) = 0.542504571396066721967792D+00
    x(177) = 0.552712861723817332466074D+00
    x(178) = 0.562838882437060514424546D+00
    x(179) = 0.572881126308666332759406D+00
    x(180) = 0.582838098581430874902446D+00
    x(181) = 0.592708317190566281032495D+00
    x(182) = 0.602490312984301547488097D+00
    x(183) = 0.612182629942561267650033D+00
    x(184) = 0.621783825393689760680446D+00
    x(185) = 0.631292470229188329449219D+00
    x(186) = 0.640707149116433684724434D+00
    x(187) = 0.650026460709345873208532D+00
    x(188) = 0.659249017856974352220492D+00
    x(189) = 0.668373447809971163711735D+00
    x(190) = 0.677398392424920474813593D+00
    x(191) = 0.686322508366494071200553D+00
    x(192) = 0.695144467307402713168813D+00
    x(193) = 0.703862956126113592426171D+00
    x(194) = 0.712476677102304460118510D+00
    x(195) = 0.720984348110025333531072D+00
    x(196) = 0.729384702808539030149808D+00
    x(197) = 0.737676490830812123299244D+00
    x(198) = 0.745858477969628263337895D+00
    x(199) = 0.753929446361296162339238D+00
    x(200) = 0.761888194666924898264210D+00
    x(201) = 0.769733538251239556788216D+00
    x(202) = 0.777464309358910595129671D+00
    x(203) = 0.785079357288370682385816D+00
    x(204) = 0.792577548563093144962574D+00
    x(205) = 0.799957767100306523636066D+00
    x(206) = 0.807218914377120130552073D+00
    x(207) = 0.814359909594035880004229D+00
    x(208) = 0.821379689835822056081139D+00
    x(209) = 0.828277210229725073186687D+00
    x(210) = 0.835051444100995681967937D+00
    x(211) = 0.841701383125706473284556D+00
    x(212) = 0.848226037480837936478636D+00
    x(213) = 0.854624435991610735314055D+00
    x(214) = 0.860895626276042275514686D+00
    x(215) = 0.867038674886706051812473D+00
    x(216) = 0.873052667449672679799858D+00
    x(217) = 0.878936708800611938658765D+00
    x(218) = 0.884689923118035575018750D+00
    x(219) = 0.890311454053661045810287D+00
    x(220) = 0.895800464859876809085345D+00
    x(221) = 0.901156138514290206476301D+00
    x(222) = 0.906377677841339419411308D+00
    x(223) = 0.911464305630951423630955D+00
    x(224) = 0.916415264754228313295468D+00
    x(225) = 0.921229818276144817520964D+00
    x(226) = 0.925907249565240289235888D+00
    x(227) = 0.930446862400288909805510D+00
    x(228) = 0.934847981073932324370129D+00
    x(229) = 0.939109950493259404355123D+00
    x(230) = 0.943232136277318328151464D+00
    x(231) = 0.947213924851546682950881D+00
    x(232) = 0.951054723539105826691801D+00
    x(233) = 0.954753960649106318830855D+00
    x(234) = 0.958311085561711847074814D+00
    x(235) = 0.961725568810109767190665D+00
    x(236) = 0.964996902159337170373447D+00
    x(237) = 0.968124598681952354372943D+00
    x(238) = 0.971108192830542793021113D+00
    x(239) = 0.973947240507062326750976D+00
    x(240) = 0.976641319128992592610888D+00
    x(241) = 0.979190027692327124191591D+00
    x(242) = 0.981592986831381877693095D+00
    x(243) = 0.983849838875444644048531D+00
    x(244) = 0.985960247902290665366669D+00
    x(245) = 0.987923899788618253106809D+00
    x(246) = 0.989740502257507526030375D+00
    x(247) = 0.991409784923101705201254D+00
    x(248) = 0.992931499332908653172844D+00
    x(249) = 0.994305419008553630362377D+00
    x(250) = 0.995531339486830143483750D+00
    x(251) = 0.996609078365487004512326D+00
    x(252) = 0.997538475365520218731818D+00
    x(253) = 0.998319392445383847808766D+00
    x(254) = 0.998951714093223210129834D+00
    x(255) = 0.999435348366365078441838D+00
    x(256) = 0.999770232390338019056053D+00
    x(257) = 0.999956390712330402472857D+00

    w(1) = 0.00011191470145601756450862287886D+00
    w(2) = 0.00026049995580176964436806680831D+00
    w(3) = 0.00040926648283531339591138751432D+00
    w(4) = 0.00055799120546880640169677292533D+00
    w(5) = 0.00070663671051592291949335494247D+00
    w(6) = 0.00085517818446696565626595950963D+00
    w(7) = 0.00100359280467969441299468763292D+00
    w(8) = 0.0011518582377826677880963146741D+00
    w(9) = 0.0012999523174235227389668643832D+00
    w(10) = 0.0014478529559255120065233994722D+00
    w(11) = 0.0015955381166175133369701690235D+00
    w(12) = 0.0017429858051468299509941139300D+00
    w(13) = 0.0018901740676190104269878470891D+00
    w(14) = 0.0020370809914723626741694800322D+00
    w(15) = 0.0021836847075455253317921866057D+00
    w(16) = 0.0023299633927021828561308282641D+00
    w(17) = 0.0024758952727301488651840215879D+00
    w(18) = 0.0026214586253808109266552781372D+00
    w(19) = 0.0027666317834818283552560256501D+00
    w(20) = 0.0029113931380877846359302447381D+00
    w(21) = 0.0030557211416493711130936102459D+00
    w(22) = 0.0031995943111899437356540290142D+00
    w(23) = 0.0033429912314827618499065991316D+00
    w(24) = 0.0034858905582247143702551557840D+00
    w(25) = 0.0036282710212037760873102463983D+00
    w(26) = 0.0037701114274582873548537007645D+00
    w(27) = 0.0039113906644266662571543468015D+00
    w(28) = 0.0040520877030864825223229951262D+00
    w(29) = 0.0041921816010820254766367595011D+00
    w(30) = 0.0043316515058396297504806208252D+00
    w(31) = 0.0044704766576701092218388764046D+00
    w(32) = 0.0046086363928577081326523656522D+00
    w(33) = 0.0047461101467350184936945641585D+00
    w(34) = 0.0048828774567433411142588306018D+00
    w(35) = 0.0050189179654779878773297516544D+00
    w(36) = 0.0051542114237180378340642003713D+00
    w(37) = 0.0052887376934400710240953933529D+00
    w(38) = 0.0054224767508154127788846727083D+00
    w(39) = 0.0055554086891904284012033890901D+00
    w(40) = 0.0056875137220494140577838938236D+00
    w(41) = 0.0058187721859596348346566361185D+00
    w(42) = 0.0059491645434980654366600347567D+00
    w(43) = 0.0060786713861593931405204596709D+00
    w(44) = 0.0062072734372448464599330978665D+00
    w(45) = 0.0063349515547314166407936938524D+00
    w(46) = 0.0064616867341210426397202932350D+00
    w(47) = 0.0065874601112693336961737372300D+00
    w(48) = 0.0067122529651934070221351960200D+00
    w(49) = 0.0068360467208584215286561508406D+00
    w(50) = 0.0069588229519423919043121805236D+00
    w(51) = 0.0070805633835788707705149901066D+00
    w(52) = 0.0072012498950770900730828552207D+00
    w(53) = 0.0073208645226191563361371026044D+00
    w(54) = 0.0074393894619338979090297315972D+00
    w(55) = 0.0075568070709469658838993300454D+00
    w(56) = 0.0076730998724067939537782250476D+00
    w(57) = 0.0077882505564860261212726654404D+00
    w(58) = 0.0079022419833580248574070864277D+00
    w(59) = 0.0080150571857480760504667455353D+00
    w(60) = 0.0081266793714589108764118189068D+00
    w(61) = 0.0082370919258701685661946145361D+00
    w(62) = 0.0083462784144114279413811886655D+00
    w(63) = 0.0084542225850084395379670551258D+00
    w(64) = 0.0085609083705021941391459209280D+00
    w(65) = 0.0086663198910404675908861979240D+00
    w(66) = 0.0087704414564414858792445834744D+00
    w(67) = 0.0088732575685293586050755892934D+00
    w(68) = 0.0089747529234409331997949023068D+00
    w(69) = 0.0090749124139037264846862498962D+00
    w(70) = 0.0091737211314845944854270065178D+00
    w(71) = 0.0092711643688088057725325917169D+00
    w(72) = 0.0093672276217491880067391857021D+00
    w(73) = 0.0094618965915850218253881576301D+00
    w(74) = 0.0095551571871303607110514249099D+00
    w(75) = 0.0096469955268314600363329731559D+00
    w(76) = 0.0097373979408330030783691793250D+00
    w(77) = 0.0098263509730128164423854701706D+00
    w(78) = 0.0099138413829847720250916955489D+00
    w(79) = 0.0099998561480695773850435626986D+00
    w(80) = 0.0100843824652331611676814627839D+00
    w(81) = 0.0101674077529923650568895461852D+00
    w(82) = 0.0102489196532876585918958554047D+00
    w(83) = 0.0103289060333225980974485876288D+00
    w(84) = 0.0104073549873697559257355517893D+00
    w(85) = 0.0104842548385428511997370260353D+00
    w(86) = 0.0105595941405348182788823332058D+00
    w(87) = 0.0106333616793215542382761147904D+00
    w(88) = 0.0107055464748310917616231511294D+00
    w(89) = 0.0107761377825779489945556541150D+00
    w(90) = 0.0108451250952624130885928632830D+00
    w(91) = 0.0109124981443345193856719616965D+00
    w(92) = 0.0109782469015224934483083029166D+00
    w(93) = 0.0110423615803254284301924654946D+00
    w(94) = 0.0111048326374699756056269264803D+00
    w(95) = 0.0111656507743308312328559850485D+00
    w(96) = 0.0112248069383148083152535688671D+00
    w(97) = 0.0112822923242082872447042603128D+00
    w(98) = 0.0113380983754878447625379269120D+00
    w(99) = 0.011392216785593866154247619654D+00
    w(100) = 0.011444639499166951104119199270D+00
    w(101) = 0.011495358713246929174010288914D+00
    w(102) = 0.011544366878434306436012137033D+00
    w(103) = 0.011591656700013970380783131035D+00
    w(104) = 0.011637221139040985841125311445D+00
    w(105) = 0.011681053413388320313049670635D+00
    w(106) = 0.011723146998756342723302879656D+00
    w(107) = 0.011763495629643945382264331878D+00
    w(108) = 0.011802093300281144573421477037D+00
    w(109) = 0.011838934265523020964443424791D+00
    w(110) = 0.011874013041704866779344562066D+00
    w(111) = 0.011907324407458412445505183140D+00
    w(112) = 0.011938863404489011222535627643D+00
    w(113) = 0.011968625338313666131272065445D+00
    w(114) = 0.011996605778959789329711050159D+00
    w(115) = 0.012022800561624589927558893338D+00
    w(116) = 0.012047205787294992091420946532D+00
    w(117) = 0.012069817823327991167612855626D+00
    w(118) = 0.012090633303991361438266420912D+00
    w(119) = 0.012109649130964635027950450318D+00
    w(120) = 0.012126862473800277391553601370D+00
    w(121) = 0.012142270770344990738801546574D+00
    w(122) = 0.012155871727121082685623083829D+00
    w(123) = 0.012167663319667843366755737416D+00
    w(124) = 0.012177643792842880196606249581D+00
    w(125) = 0.012185811661083365425569178819D+00
    w(126) = 0.012192165708627157605870499188D+00
    w(127) = 0.012196704989693764053654538465D+00
    w(128) = 0.012199428828625117371582840212D+00
    w(129) = 0.01220033681998614507777289232D+00
    w(130) = 0.012199428828625117371582840212D+00
    w(131) = 0.012196704989693764053654538465D+00
    w(132) = 0.012192165708627157605870499188D+00
    w(133) = 0.012185811661083365425569178819D+00
    w(134) = 0.012177643792842880196606249581D+00
    w(135) = 0.012167663319667843366755737416D+00
    w(136) = 0.012155871727121082685623083829D+00
    w(137) = 0.012142270770344990738801546574D+00
    w(138) = 0.012126862473800277391553601370D+00
    w(139) = 0.012109649130964635027950450318D+00
    w(140) = 0.012090633303991361438266420912D+00
    w(141) = 0.012069817823327991167612855626D+00
    w(142) = 0.012047205787294992091420946532D+00
    w(143) = 0.012022800561624589927558893338D+00
    w(144) = 0.011996605778959789329711050159D+00
    w(145) = 0.011968625338313666131272065445D+00
    w(146) = 0.011938863404489011222535627643D+00
    w(147) = 0.011907324407458412445505183140D+00
    w(148) = 0.011874013041704866779344562066D+00
    w(149) = 0.011838934265523020964443424791D+00
    w(150) = 0.011802093300281144573421477037D+00
    w(151) = 0.011763495629643945382264331878D+00
    w(152) = 0.011723146998756342723302879656D+00
    w(153) = 0.011681053413388320313049670635D+00
    w(154) = 0.011637221139040985841125311445D+00
    w(155) = 0.011591656700013970380783131035D+00
    w(156) = 0.011544366878434306436012137033D+00
    w(157) = 0.011495358713246929174010288914D+00
    w(158) = 0.011444639499166951104119199270D+00
    w(159) = 0.011392216785593866154247619654D+00
    w(160) = 0.0113380983754878447625379269120D+00
    w(161) = 0.0112822923242082872447042603128D+00
    w(162) = 0.0112248069383148083152535688671D+00
    w(163) = 0.0111656507743308312328559850485D+00
    w(164) = 0.0111048326374699756056269264803D+00
    w(165) = 0.0110423615803254284301924654946D+00
    w(166) = 0.0109782469015224934483083029166D+00
    w(167) = 0.0109124981443345193856719616965D+00
    w(168) = 0.0108451250952624130885928632830D+00
    w(169) = 0.0107761377825779489945556541150D+00
    w(170) = 0.0107055464748310917616231511294D+00
    w(171) = 0.0106333616793215542382761147904D+00
    w(172) = 0.0105595941405348182788823332058D+00
    w(173) = 0.0104842548385428511997370260353D+00
    w(174) = 0.0104073549873697559257355517893D+00
    w(175) = 0.0103289060333225980974485876288D+00
    w(176) = 0.0102489196532876585918958554047D+00
    w(177) = 0.0101674077529923650568895461852D+00
    w(178) = 0.0100843824652331611676814627839D+00
    w(179) = 0.0099998561480695773850435626986D+00
    w(180) = 0.0099138413829847720250916955489D+00
    w(181) = 0.0098263509730128164423854701706D+00
    w(182) = 0.0097373979408330030783691793250D+00
    w(183) = 0.0096469955268314600363329731559D+00
    w(184) = 0.0095551571871303607110514249099D+00
    w(185) = 0.0094618965915850218253881576301D+00
    w(186) = 0.0093672276217491880067391857021D+00
    w(187) = 0.0092711643688088057725325917169D+00
    w(188) = 0.0091737211314845944854270065178D+00
    w(189) = 0.0090749124139037264846862498962D+00
    w(190) = 0.0089747529234409331997949023068D+00
    w(191) = 0.0088732575685293586050755892934D+00
    w(192) = 0.0087704414564414858792445834744D+00
    w(193) = 0.0086663198910404675908861979240D+00
    w(194) = 0.0085609083705021941391459209280D+00
    w(195) = 0.0084542225850084395379670551258D+00
    w(196) = 0.0083462784144114279413811886655D+00
    w(197) = 0.0082370919258701685661946145361D+00
    w(198) = 0.0081266793714589108764118189068D+00
    w(199) = 0.0080150571857480760504667455353D+00
    w(200) = 0.0079022419833580248574070864277D+00
    w(201) = 0.0077882505564860261212726654404D+00
    w(202) = 0.0076730998724067939537782250476D+00
    w(203) = 0.0075568070709469658838993300454D+00
    w(204) = 0.0074393894619338979090297315972D+00
    w(205) = 0.0073208645226191563361371026044D+00
    w(206) = 0.0072012498950770900730828552207D+00
    w(207) = 0.0070805633835788707705149901066D+00
    w(208) = 0.0069588229519423919043121805236D+00
    w(209) = 0.0068360467208584215286561508406D+00
    w(210) = 0.0067122529651934070221351960200D+00
    w(211) = 0.0065874601112693336961737372300D+00
    w(212) = 0.0064616867341210426397202932350D+00
    w(213) = 0.0063349515547314166407936938524D+00
    w(214) = 0.0062072734372448464599330978665D+00
    w(215) = 0.0060786713861593931405204596709D+00
    w(216) = 0.0059491645434980654366600347567D+00
    w(217) = 0.0058187721859596348346566361185D+00
    w(218) = 0.0056875137220494140577838938236D+00
    w(219) = 0.0055554086891904284012033890901D+00
    w(220) = 0.0054224767508154127788846727083D+00
    w(221) = 0.0052887376934400710240953933529D+00
    w(222) = 0.0051542114237180378340642003713D+00
    w(223) = 0.0050189179654779878773297516544D+00
    w(224) = 0.0048828774567433411142588306018D+00
    w(225) = 0.0047461101467350184936945641585D+00
    w(226) = 0.0046086363928577081326523656522D+00
    w(227) = 0.0044704766576701092218388764046D+00
    w(228) = 0.0043316515058396297504806208252D+00
    w(229) = 0.0041921816010820254766367595011D+00
    w(230) = 0.0040520877030864825223229951262D+00
    w(231) = 0.0039113906644266662571543468015D+00
    w(232) = 0.0037701114274582873548537007645D+00
    w(233) = 0.0036282710212037760873102463983D+00
    w(234) = 0.0034858905582247143702551557840D+00
    w(235) = 0.0033429912314827618499065991316D+00
    w(236) = 0.0031995943111899437356540290142D+00
    w(237) = 0.0030557211416493711130936102459D+00
    w(238) = 0.0029113931380877846359302447381D+00
    w(239) = 0.0027666317834818283552560256501D+00
    w(240) = 0.0026214586253808109266552781372D+00
    w(241) = 0.0024758952727301488651840215879D+00
    w(242) = 0.0023299633927021828561308282641D+00
    w(243) = 0.0021836847075455253317921866057D+00
    w(244) = 0.0020370809914723626741694800322D+00
    w(245) = 0.0018901740676190104269878470891D+00
    w(246) = 0.0017429858051468299509941139300D+00
    w(247) = 0.0015955381166175133369701690235D+00
    w(248) = 0.0014478529559255120065233994722D+00
    w(249) = 0.0012999523174235227389668643832D+00
    w(250) = 0.0011518582377826677880963146741D+00
    w(251) = 0.00100359280467969441299468763292D+00
    w(252) = 0.00085517818446696565626595950963D+00
    w(253) = 0.00070663671051592291949335494247D+00
    w(254) = 0.00055799120546880640169677292533D+00
    w(255) = 0.00040926648283531339591138751432D+00
    w(256) = 0.00026049995580176964436806680831D+00
    w(257) = 0.00011191470145601756450862287886D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) &
      '  Legal values are 1 through 33, 63/64/65, 127/128/129 and 255/256/257.'
    stop

  end if

  return
end
subroutine legendre_sqrtx2_01_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_SQRTX2_01_SET: Gauss-Legendre rule for F(X) / SQRT(X) on [0,1].
!
!  Discussion:
!
!    The integral:
!
!      Integral ( 0 <= X <= 1 ) F(X) / SQRT ( X ) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) n2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) w2(2*n+1)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x2(2*n+1)

  n2 = 2 * n + 1

  call legendre_set ( n2, x2, w2 )

  x(1:n) = x2(n+2:2*n+1)**2
  w(1:n) = 2.0D+00 * w2(n+2:2*n+1)

  return
end
subroutine legendre_sqrtx_01_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_SQRTX_01_SET sets a Gauss-Legendre rule for SQRT(X) * F(X) on [0,1].
!
!  Discussion:
!
!    The integral:
!
!      Integral ( 0 <= X <= 1 ) SQRT ( X ) * F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) n2

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) w2(2*n+1)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x2(2*n+1)


  n2 = 2 * n + 1

  call legendre_set ( n2, x2, w2 )

  x(1:n) = x2(n+2:2*n+1)**2
  w(1:n) = 2.0D+00 * w2(n+2:2*n+1) * x(1:n)

  return
end
subroutine legendre_ss_compute ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_SS_COMPUTE: Gauss-Legendre quadrature by Stroud-Secrest method.
!
!  Discussion:
!
!    The Stroud and Secrest reference did not print a specific computer program
!    for the Gauss-Legendre case.  Therefore, this code is based on the
!    printed code for the Gauss-Jacobi case, with ALPHA and BETA set to 0.
!    This means that the LEGENDRE_SS_ROOT and LEGENDRE_SS_RECUR routines,
!    while appropriate for this computation, do not use the standard
!    normalization for the Legendre polynomials in which Pn(1) = 1.
!    The unusual scaling does not, however, affect the location of the
!    roots, which is the primary thing of interest.
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2009
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) cc
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p1
  real ( kind = 8 ) r
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtemp
!
!  Set the recursion coefficients.
!
  do i = 1, n

    c(i) = real ( ( i - 1 ) * ( i - 1 ), kind = 8 ) / &
           real ( ( 2 * i - 1 ) * ( 2 * i - 3 ), kind = 8 )

  end do

  cc = 2.0D+00 * product ( c(2:n) )

  do i = 1, n

    if ( i == 1 ) then

      r = 2.78D+00 / ( 4.0D+00 + real ( n * n, kind = 8 ) )

      xtemp = 1.0D+00 - r

    else if ( i == 2 ) then

      r = 1.0D+00 + 0.06D+00 * real ( n - 8, kind = 8 ) &
        / real ( n, kind = 8 )

      xtemp = xtemp - 4.1D+00 * r * ( 1.0D+00 - xtemp )

    else if ( i == 3 ) then

      r = 1.0D+00 + 0.22D+00 * real ( n - 8, kind = 8 )  &
        / real ( n, kind = 8 )

      xtemp = xtemp - 1.67D+00 * r * ( x(1) - xtemp )

    else if ( i < n - 1 ) then

      xtemp = 3.0D+00 * x(i-1) - 3.0D+00 * x(i-2) + x(i-3)

    else if ( i == n - 1 ) then

      r = 1.0D+00 / ( 1.0D+00 + 0.639D+00 * real ( n - 4, kind = 8 ) &
        / ( 1.0D+00 + 0.71D+00 * real ( n - 4, kind = 8 ) ) )

      xtemp = xtemp + r * ( xtemp - x(i-2) ) / 0.766D+00

    else if ( i == n ) then

      r = 1.0D+00 / ( 1.0D+00 + 0.22D+00 * real ( n - 8, kind = 8 ) &
        / real ( n, kind = 8 ) )

      xtemp = xtemp + r * ( xtemp - x(i-2) ) / 1.67D+00

    end if

    call legendre_ss_root ( xtemp, n, dp2, p1, c )

    x(i) = xtemp
    w(i) = cc / dp2 / p1

  end do
!
!  Reverse the data.
!
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

  return
end
subroutine legendre_ss_recur ( p2, dp2, p1, x, n, c )

!*****************************************************************************80
!
!! LEGENDRE_SS_RECUR: value and derivative of a scaled Legendre polynomial.
!
!  Discussion:
!
!    The Stroud and Secrest reference did not print a specific computer program
!    for the Gauss-Legendre case.  Therefore, this code is based on the
!    printed code for the Gauss-Jacobi case, with ALPHA and BETA set to 0.
!    This means that the LEGENDRE_SS_ROOT and LEGENDRE_SS_RECUR routines,
!    while appropriate for this computation, do not use the standard
!    normalization for the Legendre polynomials in which Pn(1) = 1.
!    The unusual scaling does not, however, affect the location of the
!    roots, which is the primary thing of interest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P2, the value of L(N)(X).
!
!    Output, real ( kind = 8 ) DP2, the value of L'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(N-1)(X).
!
!    Input, real ( kind = 8 ) X, the point at which polynomials are evaluated.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) dp0
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dp2
  integer ( kind = 4 ) i
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, n

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = x * p1 - c(i) * p0
    dp2 = x * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine legendre_ss_root ( x, n, dp2, p1, c )

!*****************************************************************************80
!
!! LEGENDRE_SS_ROOT: improve approximate root of scaled Legendre polynomial.
!
!  Discussion:
!
!    The Stroud and Secrest reference did not print a specific computer program
!    for the Gauss-Legendre case.  Therefore, this code is based on the
!    printed code for the Gauss-Jacobi case, with ALPHA and BETA set to 0.
!    This means that the LEGENDRE_SS_ROOT and LEGENDRE_SS_RECUR routines,
!    while appropriate for this computation, do not use the standard
!    normalization for the Legendre polynomials in which Pn(1) = 1.
!    The unusual scaling does not, however, affect the location of the
!    roots, which is the primary thing of interest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) DP2, the value of L'(N)(X).
!
!    Output, real ( kind = 8 ) P1, the value of L(N-1)(X).
!
!    Input, real ( kind = 8 ) C(N), the recursion coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d
  real ( kind = 8 ) dp2
  real ( kind = 8 ) eps
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  integer ( kind = 4 ) step
  integer ( kind = 4 ), parameter :: step_max = 10
  real ( kind = 8 ) x

  eps = epsilon ( eps )

  do step = 1, step_max

    call legendre_ss_recur ( p2, dp2, p1, x, n, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine legendre_x0_01_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_X0_01_SET sets a Gauss-Legendre rule for F(X) on [0,1].
!
!  Discussion:
!
!    The integral:
!
!      Integral ( 0 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 8.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = 0.5D+00

    w(1) = 1.0D+00

  else if ( n == 2 ) then

    x(1) = 0.2113248654D+00
    x(2) = 0.7886751346D+00

    w(1) = 0.5D+00
    w(2) = 0.5D+00

  else if ( n == 3 ) then

    x(1) = 0.1127016654D+00
    x(2) = 0.5000000000D+00
    x(3) = 0.8872983346D+00

    w(1) = 5.0D+00 / 18.0D+00
    w(2) = 8.0D+00 / 18.0D+00
    w(3) = 5.0D+00 / 18.0D+00

  else if ( n == 4 ) then

    x(1) = 0.0694318442D+00
    x(2) = 0.3300094782D+00
    x(3) = 0.6699905218D+00
    x(4) = 0.9305681558D+00

    w(1) = 0.1739274226D+00
    w(2) = 0.3260725774D+00
    w(3) = 0.3260725774D+00
    w(4) = 0.1739274226D+00

  else if ( n == 5 ) then

    x(1) = 0.0469100770D+00
    x(2) = 0.2307653449D+00
    x(3) = 0.5000000000D+00
    x(4) = 0.7692346551D+00
    x(5) = 0.9530899230D+00

    w(1) = 0.1184634425D+00
    w(2) = 0.2393143352D+00
    w(3) = 0.2844444444D+00
    w(4) = 0.2393143352D+00
    w(5) = 0.1184634425D+00

  else if ( n == 6 ) then

    x(1) = 0.0337652429D+00
    x(2) = 0.1693953068D+00
    x(3) = 0.3806904070D+00
    x(4) = 0.6193095930D+00
    x(5) = 0.8306046932D+00
    x(6) = 0.9662347571D+00

    w(1) = 0.0856622462D+00
    w(2) = 0.1803807865D+00
    w(3) = 0.2339569673D+00
    w(4) = 0.2339569673D+00
    w(5) = 0.1803807865D+00
    w(6) = 0.0856622462D+00

  else if ( n == 7 ) then

    x(1) = 0.0254460438D+00
    x(2) = 0.1292344072D+00
    x(3) = 0.2970774243D+00
    x(4) = 0.5000000000D+00
    x(5) = 0.7029225757D+00
    x(6) = 0.8707655928D+00
    x(7) = 0.9745539562D+00

    w(1) = 0.0647424831D+00
    w(2) = 0.1398526957D+00
    w(3) = 0.1909150253D+00
    w(4) = 0.2089795918D+00
    w(5) = 0.1909150253D+00
    w(6) = 0.1398526957D+00
    w(7) = 0.0647424831D+00

  else if ( n == 8 ) then

    x(1) = 0.0198550718D+00
    x(2) = 0.1016667613D+00
    x(3) = 0.2372337950D+00
    x(4) = 0.4082826788D+00
    x(5) = 0.5917173212D+00
    x(6) = 0.7627662050D+00
    x(7) = 0.8983332387D+00
    x(8) = 0.9801449282D+00

    w(1) = 0.0506142681D+00
    w(2) = 0.1111905172D+00
    w(3) = 0.1568533229D+00
    w(4) = 0.1813418917D+00
    w(5) = 0.1813418917D+00
    w(6) = 0.1568533229D+00
    w(7) = 0.1111905172D+00
    w(8) = 0.0506142681D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_X0_01_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  end if

  return
end
subroutine legendre_x1_01_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_X1_01_SET sets a Gauss-Legendre rule for X * F(X) on [0,1].
!
!  Discussion:
!
!    The integral:
!
!      Integral ( 0 <= X <= 1 ) X * F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 8.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =   0.6666666667D+00

    w(1) = 0.5000000000D+00

  else if ( n == 2 ) then

    x(1) = 0.3550510257D+00
    x(2) = 0.8449489743D+00

    w(1) = 0.1819586183D+00
    w(2) = 0.3180413817D+00

  else if ( n == 3 ) then

    x(1) = 0.2123405382D+00
    x(2) = 0.5905331356D+00
    x(3) = 0.9114120405D+00

    w(1) = 0.0698269799D+00
    w(2) = 0.2292411064D+00
    w(3) = 0.2009319137D+00

  else if ( n == 4 ) then

    x(1) = 0.1397598643D+00
    x(2) = 0.4164095676D+00
    x(3) = 0.7231569864D+00
    x(4) = 0.9428958039D+00

    w(1) = 0.0311809710D+00
    w(2) = 0.1298475476D+00
    w(3) = 0.2034645680D+00
    w(4) = 0.1355069134D+00

  else if ( n == 5 ) then

    x(1) = 0.0985350858D+00
    x(2) = 0.3045357266D+00
    x(3) = 0.5620251898D+00
    x(4) = 0.8019865821D+00
    x(5) = 0.9601901429D+00

    w(1) = 0.0157479145D+00
    w(2) = 0.0739088701D+00
    w(3) = 0.1463869871D+00
    w(4) = 0.1671746381D+00
    w(5) = 0.0967815902D+00

  else if ( n == 6 ) then

    x(1) = 0.0730543287D+00
    x(2) = 0.2307661380D+00
    x(3) = 0.4413284812D+00
    x(4) = 0.6630153097D+00
    x(5) = 0.8519214003D+00
    x(6) = 0.9706835728D+00

    w(1) = 0.0087383108D+00
    w(2) = 0.0439551656D+00
    w(3) = 0.0986611509D+00
    w(4) = 0.1407925538D+00
    w(5) = 0.1355424972D+00
    w(6) = 0.0723103307D+00

  else if ( n == 7 ) then

    x(1) = 0.0562625605D+00
    x(2) = 0.1802406917D+00
    x(3) = 0.3526247171D+00
    x(4) = 0.5471536263D+00
    x(5) = 0.7342101772D+00
    x(6) = 0.8853209468D+00
    x(7) = 0.9775206136D+00

    w(1) = 0.0052143622D+00
    w(2) = 0.0274083567D+00
    w(3) = 0.0663846965D+00
    w(4) = 0.1071250657D+00
    w(5) = 0.1273908973D+00
    w(6) = 0.1105092582D+00
    w(7) = 0.0559673634D+00

  else if ( n == 8 ) then

    x(1) = 0.0446339553D+00
    x(2) = 0.1443662570D+00
    x(3) = 0.2868247571D+00
    x(4) = 0.4548133152D+00
    x(5) = 0.6280678354D+00
    x(6) = 0.7856915206D+00
    x(7) = 0.9086763921D+00
    x(8) = 0.9822200849D+00

    w(1) = 0.0032951914D+00
    w(2) = 0.0178429027D+00
    w(3) = 0.0454393195D+00
    w(4) = 0.0791995995D+00
    w(5) = 0.1060473594D+00
    w(6) = 0.1125057995D+00
    w(7) = 0.0911190236D+00
    w(8) = 0.0445508044D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_X1_01_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  end if

  return
end
subroutine legendre_x1_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_X1_SET sets a Gauss-Legendre rule for ( 1 + X ) * F(X) on [-1,1].
!
!  Discussion:
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) ( 1 + X ) * F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 9.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =  0.333333333333333333333333333333D+00

    w(1) = 2.0D+00

  else if ( n == 2 ) then

    x(1) = -0.289897948556635619639456814941D+00
    x(2) =  0.689897948556635619639456814941D+00

    w(1) =  0.727834473024091322422523991699D+00
    w(2) =  1.27216552697590867757747600830D+00

  else if ( n == 3 ) then

    x(1) = -0.575318923521694112050483779752D+00
    x(2) =  0.181066271118530578270147495862D+00
    x(3) =  0.822824080974592105208907712461D+00

    w(1) =  0.279307919605816490135525088716D+00
    w(2) =  0.916964425438344986775682378225D+00
    w(3) =  0.803727654955838523088792533058D+00

  else if ( n == 4 ) then

    x(1) = -0.720480271312438895695825837750D+00
    x(2) = -0.167180864737833640113395337326D+00
    x(3) =  0.446313972723752344639908004629D+00
    x(4) =  0.885791607770964635613757614892D+00

    w(1) =  0.124723883800032328695500588386D+00
    w(2) =  0.519390190432929763305824811559D+00
    w(3) =  0.813858272041085443165617903743D+00
    w(4) =  0.542027653725952464833056696312D+00

  else if ( n == 5 ) then

    x(1) = -0.802929828402347147753002204224D+00
    x(2) = -0.390928546707272189029229647442D+00
    x(3) =  0.124050379505227711989974959990D+00
    x(4) =  0.603973164252783654928415726409D+00
    x(5) =  0.920380285897062515318386619813D+00

    w(1) =  0.0629916580867691047411692662740D+00
    w(2) =  0.295635480290466681402532877367D+00
    w(3) =  0.585547948338679234792151477424D+00
    w(4) =  0.668698552377478261966702492391D+00
    w(5) =  0.387126360906606717097443886545D+00

  else if ( n == 6 ) then

    x(1) = -0.853891342639482229703747931639D+00
    x(2) = -0.538467724060109001833766720231D+00
    x(3) = -0.117343037543100264162786683611D+00
    x(4) =  0.326030619437691401805894055838D+00
    x(5) =  0.703842800663031416300046295008D+00
    x(6) =  0.941367145680430216055899446174D+00

    w(1) =  0.0349532072544381270240692132496D+00
    w(2) =  0.175820662202035902032706497222D+00
    w(3) =  0.394644603562621056482338042193D+00
    w(4) =  0.563170215152795712476307356284D+00
    w(5) =  0.542169988926074467362761586552D+00
    w(6) =  0.289241322902034734621817304499D+00

  else if ( n == 7 ) then

    x(1) = -0.887474878926155707068695617935D+00
    x(2) = -0.639518616526215270024840114382D+00
    x(3) = -0.294750565773660725252184459658D+00
    x(4) =  0.0943072526611107660028971153047D+00
    x(5) =  0.468420354430821063046421216613D+00
    x(6) =  0.770641893678191536180719525865D+00
    x(7) =  0.955041227122575003782349000858D+00

    w(1) =  0.0208574488112296163587654972151D+00
    w(2) =  0.109633426887493901777324193433D+00
    w(3) =  0.265538785861965879934591955055D+00
    w(4) =  0.428500262783494679963649011999D+00
    w(5) =  0.509563589198353307674937943100D+00
    w(6) =  0.442037032763498409684482945478D+00
    w(7) =  0.223869453693964204606248453720D+00

  else if ( n == 8 ) then

    x(1) = -0.910732089420060298533757956283D+00
    x(2) = -0.711267485915708857029562959544D+00
    x(3) = -0.426350485711138962102627520502D+00
    x(4) = -0.0903733696068532980645444599064D+00
    x(5) =  0.256135670833455395138292079035D+00
    x(6) =  0.571383041208738483284917464837D+00
    x(7) =  0.817352784200412087992517083851D+00
    x(8) =  0.964440169705273096373589797925D+00

    w(1) =  0.0131807657689951954189692640444D+00
    w(2) =  0.0713716106239448335742111888042D+00
    w(3) =  0.181757278018795592332221684383D+00
    w(4) =  0.316798397969276640481632757440D+00
    w(5) =  0.424189437743720042818124385645D+00
    w(6) =  0.450023197883549464687088394417D+00
    w(7) =  0.364476094545494505382889847132D+00
    w(8) =  0.178203217446223725304862478136D+00

  else if ( n == 9 ) then

    x(1) = -0.927484374233581078117671398464D+00
    x(2) = -0.763842042420002599615429776011D+00
    x(3) = -0.525646030370079229365386614293D+00
    x(4) = -0.236234469390588049278459503207D+00
    x(5) =  0.0760591978379781302337137826389D+00
    x(6) =  0.380664840144724365880759065541D+00
    x(7) =  0.647766687674009436273648507855D+00
    x(8) =  0.851225220581607910728163628088D+00
    x(9) =  0.971175180702246902734346518378D+00

    w(1) =  0.00872338834309252349019620448007D+00
    w(2) =  0.0482400171391415162069086091476D+00
    w(3) =  0.127219285964216005046760427743D+00
    w(4) =  0.233604781180660442262926091607D+00
    w(5) =  0.337433287379681397577000079834D+00
    w(6) =  0.401235236773473158616600898930D+00
    w(7) =  0.394134968689382820640692081477D+00
    w(8) =  0.304297020437232650320317215016D+00
    w(9) =  0.145112014093119485838598391765D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_X1_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of N = ', n
    stop

  end if

  return
end
subroutine legendre_x2_01_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_X2_01_SET sets a Gauss-Legendre rule for X * X * F(X) on [0,1].
!
!  Discussion:
!
!    The integral:
!
!      Integral ( 0 <= X <= 1 ) X * X * F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 8.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =   0.75D+00

    w(1) = 1.0D+00 / 3.0D+00

  else if ( n == 2 ) then

    x(1) = 0.4558481560D+00
    x(2) = 0.8774851773D+00

    w(1) = 0.1007858821D+00
    w(2) = 0.2325474513D+00

  else if ( n == 3 ) then

    x(1) = 0.2949977901D+00
    x(2) = 0.6529962340D+00
    x(3) = 0.9270059759D+00

    w(1) = 0.0299507030D+00
    w(2) = 0.1462462693D+00
    w(3) = 0.1571363611D+00

  else if ( n == 4 ) then

    x(1) = 0.2041485821D+00
    x(2) = 0.4829527049D+00
    x(3) = 0.7613992624D+00
    x(4) = 0.9514994506D+00

    w(1) = 0.0103522408D+00
    w(2) = 0.0686338872D+00
    w(3) = 0.1434587898D+00
    w(4) = 0.1108884156D+00

  else if ( n == 5 ) then

    x(1) = 0.1489457871D+00
    x(2) = 0.3656665274D+00
    x(3) = 0.6101136129D+00
    x(4) = 0.8265196792D+00
    x(5) = 0.9654210601D+00

    w(1) = 0.0041138252D+00
    w(2) = 0.0320556007D+00
    w(3) = 0.0892001612D+00
    w(4) = 0.1261989619D+00
    w(5) = 0.0817647843D+00

  else if ( n == 6 ) then

    x(1) = 0.1131943838D+00
    x(2) = 0.2843188727D+00
    x(3) = 0.4909635868D+00
    x(4) = 0.6975630820D+00
    x(5) = 0.8684360583D+00
    x(6) = 0.9740954449D+00

    w(1) = 0.0018310758D+00
    w(2) = 0.0157202972D+00
    w(3) = 0.0512895711D+00
    w(4) = 0.0945771867D+00
    w(5) = 0.1073764997D+00
    w(6) = 0.0625387027D+00

  else if ( n == 7 ) then

    x(1) = 0.0888168334D+00
    x(2) = 0.2264827534D+00
    x(3) = 0.3999784867D+00
    x(4) = 0.5859978554D+00
    x(5) = 0.7594458740D+00
    x(6) = 0.8969109709D+00
    x(7) = 0.9798672262D+00

    w(1) = 0.0008926880D+00
    w(2) = 0.0081629256D+00
    w(3) = 0.0294222113D+00
    w(4) = 0.0631463787D+00
    w(5) = 0.0917338033D+00
    w(6) = 0.0906988246D+00
    w(7) = 0.0492765018D+00

  else if ( n == 8 ) then

    x(1) = 0.0714910350D+00
    x(2) = 0.1842282964D+00
    x(3) = 0.3304477282D+00
    x(4) = 0.4944029218D+00
    x(5) = 0.6583480085D+00
    x(6) = 0.8045248315D+00
    x(7) = 0.9170993825D+00
    x(8) = 0.9839022404D+00

    w(1) = 0.0004685178D+00
    w(2) = 0.0044745217D+00
    w(3) = 0.0172468638D+00
    w(4) = 0.0408144264D+00
    w(5) = 0.0684471834D+00
    w(6) = 0.0852847692D+00
    w(7) = 0.0768180933D+00
    w(8) = 0.0397789578D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_X2_01_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  end if

  return
end
subroutine legendre_x2_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_X2_SET sets a Gauss-Legendre rule for ( 1 + X )^2 * F(X) on [-1,1].
!
!  Discussion:
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) ( 1 + X ) * ( 1 + X ) * F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 9.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =  0.5D+00

    w(1) =  2.66666666666666666666666666666D+00

  else if ( n == 2 ) then

    x(1) = -0.0883036880224505775998524725910D+00
    x(2) =  0.754970354689117244266519139258D+00

    w(1) =  0.806287056638603444666851075928D+00
    w(2) =  1.86037961002806322199981559074D+00

  else if ( n == 3 ) then

    x(1) = -0.410004419776996766244796955168D+00
    x(2) =  0.305992467923296230556472913192D+00
    x(3) =  0.854011951853700535688324041976D+00

    w(1) =  0.239605624068645584091811926047D+00
    w(2) =  1.16997015407892817602809616291D+00
    w(3) =  1.25709088851909290654675857771D+00

  else if ( n == 4 ) then

    x(1) = -0.591702835793545726606755921586D+00
    x(2) = -0.0340945902087350046811467387661D+00
    x(3) =  0.522798524896275389882037174551D+00
    x(4) =  0.902998901106005341405865485802D+00

    w(1) =  0.0828179259993445222751812523731D+00
    w(2) =  0.549071097383384602539010760334D+00
    w(3) =  1.14767031839371367238662411421D+00
    w(4) =  0.887107324890223869465850539752D+00

  else if ( n == 5 ) then

    x(1) = -0.702108425894032836232448374820D+00
    x(2) = -0.268666945261773544694327777841D+00
    x(3) =  0.220227225868961343518209179230D+00
    x(4) =  0.653039358456608553790815164028D+00
    x(5) =  0.930842120163569816951085142737D+00

    w(1) =  0.0329106016247920636689299329544D+00
    w(2) =  0.256444805783695354037991444453D+00
    w(3) =  0.713601289772720001490035944563D+00
    w(4) =  1.00959169519929190423066348132D+00
    w(5) =  0.654118274286167343239045863379D+00

  else if ( n == 6 ) then

    x(1) = -0.773611232355123732602532012021D+00
    x(2) = -0.431362254623427837535325249187D+00
    x(3) = -0.0180728263295041680220798103354D+00
    x(4) =  0.395126163954217534500188844163D+00
    x(5) =  0.736872116684029732026178298518D+00
    x(6) =  0.948190889812665614490712786006D+00

    w(1) =  0.0146486064549543818622276447204D+00
    w(2) =  0.125762377479560410622810097040D+00
    w(3) =  0.410316569036929681761034600615D+00
    w(4) =  0.756617493988329628546336413760D+00
    w(5) =  0.859011997894245060846045458784D+00
    w(6) =  0.500309621812647503028212451747D+00

  else if ( n == 7 ) then

    x(1) = -0.822366333126005527278634734418D+00
    x(2) = -0.547034493182875002223997992852D+00
    x(3) = -0.200043026557985860387937545780D+00
    x(4) =  0.171995710805880507163425502299D+00
    x(5) =  0.518891747903884926692601716998D+00
    x(6) =  0.793821941703901970495546427988D+00
    x(7) =  0.959734452453198985538996625765D+00

    w(1) =  0.00714150426951365443207221475404D+00
    w(2) =  0.0653034050584375560578544725498D+00
    w(3) =  0.235377690316228918725962815880D+00
    w(4) =  0.505171029671130381676271523850D+00
    w(5) =  0.733870426238362032891332767175D+00
    w(6) =  0.725590596901489156295739839779D+00
    w(7) =  0.394212014211504966587433032679D+00

  else if ( n == 8 ) then

    x(1) = -0.857017929919813794402037235698D+00
    x(2) = -0.631543407166567521509503573952D+00
    x(3) = -0.339104543648722903660229021109D+00
    x(4) = -0.0111941563689783438801237300122D+00
    x(5) =  0.316696017045595559454075475675D+00
    x(6) =  0.609049663022520165351466780939D+00
    x(7) =  0.834198765028697794599267293239D+00
    x(8) =  0.967804480896157932935972899807D+00

    w(1) =  0.00374814227227757804631954025851D+00
    w(2) =  0.0357961737041152639660521680263D+00
    w(3) =  0.137974910241879862433949246199D+00
    w(4) =  0.326515411108352185491692769217D+00
    w(5) =  0.547577467373226177976217604887D+00
    w(6) =  0.682278153375510121675529810121D+00
    w(7) =  0.614544746137780998436053880546D+00
    w(8) =  0.318231662453524478640851647411D+00

  else if ( n == 9 ) then

    x(1) = -0.882491728426548422828684254270D+00
    x(2) = -0.694873684026474640346360850039D+00
    x(3) = -0.446537143480670863635920316400D+00
    x(4) = -0.159388112702326252531544826624D+00
    x(5) =  0.141092709224374414981503995427D+00
    x(6) =  0.428217823321559204544020866175D+00
    x(7) =  0.676480966471850715860378175342D+00
    x(8) =  0.863830940812464825046988286026D+00
    x(9) =  0.973668228805771018909618924364D+00

    w(1) =  0.00209009877215570354392734918986D+00
    w(2) =  0.0205951891648697848186537272448D+00
    w(3) =  0.0832489326348178964194106978875D+00
    w(4) =  0.210746247220398685903797568021D+00
    w(5) =  0.388325022916052063676224499399D+00
    w(6) =  0.554275165518437673725822282791D+00
    w(7) =  0.621388553284444032628761363828D+00
    w(8) =  0.523916296267173054255512857631D+00
    w(9) =  0.262081160888317771694556320674D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_X2_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of N = ', n
    stop

  end if

  return
end
subroutine lobatto_compute ( n, x, w )

!*****************************************************************************80
!
!! LOBATTO_COMPUTE computes a Lobatto quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*N-3).
!
!    The Lobatto rule is distinguished by the fact that both endpoints
!    (-1 and 1) are always abscissas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2007
!
!  Author:
!
!    Original MATLAB version by Greg von Winckel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
!    Spectral Methods in Fluid Dynamics,
!    Springer, 1993,
!    ISNB13: 978-3540522058,
!    LC: QA377.S676.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be at least 2.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(n,n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) tolerance
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xold(n)

  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOBATTO_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) ' N must be at least 2.'
    stop
  end if

  tolerance = 100.0D+00 * epsilon ( tolerance )
!
!  Initial estimate for the abscissas is the Chebyshev-Gauss-Lobatto nodes.
!
  do i = 1, n
    x(i) = cos ( pi * real ( i - 1, kind = 8 ) / real ( n - 1, kind = 8 ) )
  end do

  xold(1:n) = 2.0D+00

  do while ( tolerance < maxval ( abs ( x(1:n) - xold(1:n) ) ) )

    xold(1:n) = x(1:n)

    p(1:n,1) = 1.0D+00
    p(1:n,2) = x(1:n)
    do j = 2, n-1
      p(1:n,j+1) = ( real ( 2 * j - 1, kind = 8 ) * x(1:n) * p(1:n,j)     &
                   + real (   - j + 1, kind = 8 ) *          p(1:n,j-1) ) &
                   / real (     j,     kind = 8 )
    end do

    x(1:n) = xold(1:n) - ( x(1:n) * p(1:n,n) - p(1:n,n-1) ) &
             / ( real ( n, kind = 8 ) * p(1:n,n) )
  end do

  x(1:n) = x(n:1:-1)
  w(1:n) = 2.0D+00 / ( real ( ( n - 1 ) * n, kind = 8 ) * p(1:n,n)**2 )

  return
end
subroutine lobatto_set ( n, x, w )

!*****************************************************************************80
!
!! LOBATTO_SET sets abscissas and weights for Lobatto quadrature.
!
!  Discussion:
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*N-3).
!
!    The Lobatto rule is distinguished by the fact that both endpoints
!    (-1 and 1) are always abscissas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 2 and 20.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOBATTO_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are between 2 and 20.'
    stop

  else if ( n == 2 ) then

    x(1) =  - 1.0D+00
    x(2) =    1.0D+00

    w(1) =  1.0D+00
    w(2) =  1.0D+00

  else if ( n == 3 ) then

    x(1) =  - 1.0D+00
    x(2) =    0.0D+00
    x(3) =    1.0D+00

    w(1) =  1.0D+00 / 3.0D+00
    w(2) =  4.0D+00 / 3.0D+00
    w(3) =  1.0D+00 / 3.0D+00

  else if ( n == 4 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.447213595499957939281834733746D+00
    x(3) =    0.447213595499957939281834733746D+00
    x(4) =    1.0D+00

    w(1) =  1.0D+00 / 6.0D+00
    w(2) =  5.0D+00 / 6.0D+00
    w(3) =  5.0D+00 / 6.0D+00
    w(4) =  1.0D+00 / 6.0D+00

  else if ( n == 5 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.654653670707977143798292456247D+00
    x(3) =    0.0D+00
    x(4) =    0.654653670707977143798292456247D+00
    x(5) =    1.0D+00

    w(1) =  9.0D+00 / 90.0D+00
    w(2) = 49.0D+00 / 90.0D+00
    w(3) = 64.0D+00 / 90.0D+00
    w(4) = 49.0D+00 / 90.0D+00
    w(5) =  9.0D+00 / 90.0D+00

  else if ( n == 6 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.765055323929464692851002973959D+00
    x(3) =  - 0.285231516480645096314150994041D+00
    x(4) =    0.285231516480645096314150994041D+00
    x(5) =    0.765055323929464692851002973959D+00
    x(6) =    1.0D+00

    w(1) =  0.066666666666666666666666666667D+00
    w(2) =  0.378474956297846980316612808212D+00
    w(3) =  0.554858377035486353016720525121D+00
    w(4) =  0.554858377035486353016720525121D+00
    w(5) =  0.378474956297846980316612808212D+00
    w(6) =  0.066666666666666666666666666667D+00

  else if ( n == 7 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.830223896278566929872032213967D+00
    x(3) =  - 0.468848793470714213803771881909D+00
    x(4) =    0.0D+00
    x(5) =    0.468848793470714213803771881909D+00
    x(6) =    0.830223896278566929872032213967D+00
    x(7) =    1.0D+00

    w(1) =  0.476190476190476190476190476190D-01
    w(2) =  0.276826047361565948010700406290D+00
    w(3) =  0.431745381209862623417871022281D+00
    w(4) =  0.487619047619047619047619047619D+00
    w(5) =  0.431745381209862623417871022281D+00
    w(6) =  0.276826047361565948010700406290D+00
    w(7) =  0.476190476190476190476190476190D-01

  else if ( n == 8 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.871740148509606615337445761221D+00
    x(3) =  - 0.591700181433142302144510731398D+00
    x(4) =  - 0.209299217902478868768657260345D+00
    x(5) =    0.209299217902478868768657260345D+00
    x(6) =    0.591700181433142302144510731398D+00
    x(7) =    0.871740148509606615337445761221D+00
    x(8) =    1.0D+00

    w(1) =  0.357142857142857142857142857143D-01
    w(2) =  0.210704227143506039382991065776D+00
    w(3) =  0.341122692483504364764240677108D+00
    w(4) =  0.412458794658703881567052971402D+00
    w(5) =  0.412458794658703881567052971402D+00
    w(6) =  0.341122692483504364764240677108D+00
    w(7) =  0.210704227143506039382991065776D+00
    w(8) =  0.357142857142857142857142857143D-01

  else if ( n == 9 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.899757995411460157312345244418D+00
    x(3) =  - 0.677186279510737753445885427091D+00
    x(4) =  - 0.363117463826178158710752068709D+00
    x(5) =    0.0D+00
    x(6) =    0.363117463826178158710752068709D+00
    x(7) =    0.677186279510737753445885427091D+00
    x(8) =    0.899757995411460157312345244418D+00
    x(9) =    1.0D+00

    w(1) =  0.277777777777777777777777777778D-01
    w(2) =  0.165495361560805525046339720029D+00
    w(3) =  0.274538712500161735280705618579D+00
    w(4) =  0.346428510973046345115131532140D+00
    w(5) =  0.371519274376417233560090702948D+00
    w(6) =  0.346428510973046345115131532140D+00
    w(7) =  0.274538712500161735280705618579D+00
    w(8) =  0.165495361560805525046339720029D+00
    w(9) =  0.277777777777777777777777777778D-01

  else if ( n == 10 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.919533908166458813828932660822D+00
    x(3) =  - 0.738773865105505075003106174860D+00
    x(4) =  - 0.477924949810444495661175092731D+00
    x(5) =  - 0.165278957666387024626219765958D+00
    x(6) =    0.165278957666387024626219765958D+00
    x(7) =    0.477924949810444495661175092731D+00
    x(8) =    0.738773865105505075003106174860D+00
    x(9) =    0.919533908166458813828932660822D+00
    x(10) =   1.0D+00

    w(1) =  0.222222222222222222222222222222D-01
    w(2) =  0.133305990851070111126227170755D+00
    w(3) =  0.224889342063126452119457821731D+00
    w(4) =  0.292042683679683757875582257374D+00
    w(5) =  0.327539761183897456656510527917D+00
    w(6) =  0.327539761183897456656510527917D+00
    w(7) =  0.292042683679683757875582257374D+00
    w(8) =  0.224889342063126452119457821731D+00
    w(9) =  0.133305990851070111126227170755D+00
    w(10) = 0.222222222222222222222222222222D-01

  else if ( n == 11 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.934001430408059134332274136099D+00
    x(3) =  - 0.784483473663144418622417816108D+00
    x(4) =  - 0.565235326996205006470963969478D+00
    x(5) =  - 0.295758135586939391431911515559D+00
    x(6) =    0.0D+00
    x(7) =    0.295758135586939391431911515559D+00
    x(8) =    0.565235326996205006470963969478D+00
    x(9) =    0.784483473663144418622417816108D+00
    x(10) =   0.934001430408059134332274136099D+00
    x(11) =   1.0D+00

    w(1) =  0.181818181818181818181818181818D-01
    w(2) =  0.109612273266994864461403449580D+00
    w(3) =  0.187169881780305204108141521899D+00
    w(4) =  0.248048104264028314040084866422D+00
    w(5) =  0.286879124779008088679222403332D+00
    w(6) =  0.300217595455690693785931881170D+00
    w(7) =  0.286879124779008088679222403332D+00
    w(8) =  0.248048104264028314040084866422D+00
    w(9) =  0.187169881780305204108141521899D+00
    w(10) = 0.109612273266994864461403449580D+00
    w(11) = 0.181818181818181818181818181818D-01

  else if ( n == 12 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.944899272222882223407580138303D+00
    x(3) =  - 0.819279321644006678348641581717D+00
    x(4) =  - 0.632876153031869677662404854444D+00
    x(5) =  - 0.399530940965348932264349791567D+00
    x(6) =  - 0.136552932854927554864061855740D+00
    x(7) =    0.136552932854927554864061855740D+00
    x(8) =    0.399530940965348932264349791567D+00
    x(9) =    0.632876153031869677662404854444D+00
    x(10) =   0.819279321644006678348641581717D+00
    x(11) =   0.944899272222882223407580138303D+00
    x(12) =   1.0D+00

    w(1) =  0.151515151515151515151515151515D-01
    w(2) =  0.916845174131961306683425941341D-01
    w(3) =  0.157974705564370115164671062700D+00
    w(4) =  0.212508417761021145358302077367D+00
    w(5) =  0.251275603199201280293244412148D+00
    w(6) =  0.271405240910696177000288338500D+00
    w(7) =  0.271405240910696177000288338500D+00
    w(8) =  0.251275603199201280293244412148D+00
    w(9) =  0.212508417761021145358302077367D+00
    w(10) = 0.157974705564370115164671062700D+00
    w(11) = 0.916845174131961306683425941341D-01
    w(12) = 0.151515151515151515151515151515D-01

  else if ( n == 13 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.953309846642163911896905464755D+00
    x(3) =  - 0.846347564651872316865925607099D+00
    x(4) =  - 0.686188469081757426072759039566D+00
    x(5) =  - 0.482909821091336201746937233637D+00
    x(6) =  - 0.249286930106239992568673700374D+00
    x(7) =    0.0D+00
    x(8) =    0.249286930106239992568673700374D+00
    x(9) =    0.482909821091336201746937233637D+00
    x(10) =   0.686188469081757426072759039566D+00
    x(11) =   0.846347564651872316865925607099D+00
    x(12) =   0.953309846642163911896905464755D+00
    x(13) =   1.0D+00

    w(1) =  0.128205128205128205128205128205D-01
    w(2) =  0.778016867468189277935889883331D-01
    w(3) =  0.134981926689608349119914762589D+00
    w(4) =  0.183646865203550092007494258747D+00
    w(5) =  0.220767793566110086085534008379D+00
    w(6) =  0.244015790306676356458578148360D+00
    w(7) =  0.251930849333446736044138641541D+00
    w(8) =  0.244015790306676356458578148360D+00
    w(9) =  0.220767793566110086085534008379D+00
    w(10) = 0.183646865203550092007494258747D+00
    w(11) = 0.134981926689608349119914762589D+00
    w(12) = 0.778016867468189277935889883331D-01
    w(13) = 0.128205128205128205128205128205D-01

  else if ( n == 14 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.959935045267260901355100162015D+00
    x(3) =  - 0.867801053830347251000220202908D+00
    x(4) =  - 0.728868599091326140584672400521D+00
    x(5) =  - 0.550639402928647055316622705859D+00
    x(6) =  - 0.342724013342712845043903403642D+00
    x(7) =  - 0.116331868883703867658776709736D+00
    x(8) =    0.116331868883703867658776709736D+00
    x(9) =    0.342724013342712845043903403642D+00
    x(10) =   0.550639402928647055316622705859D+00
    x(11) =   0.728868599091326140584672400521D+00
    x(12) =   0.867801053830347251000220202908D+00
    x(13) =   0.959935045267260901355100162015D+00
    x(14) =   1.0D+00

    w(1) =  0.109890109890109890109890109890D-01
    w(2) =  0.668372844976812846340706607461D-01
    w(3) =  0.116586655898711651540996670655D+00
    w(4) =  0.160021851762952142412820997988D+00
    w(5) =  0.194826149373416118640331778376D+00
    w(6) =  0.219126253009770754871162523954D+00
    w(7) =  0.231612794468457058889628357293D+00
    w(8) =  0.231612794468457058889628357293D+00
    w(9) =  0.219126253009770754871162523954D+00
    w(10) = 0.194826149373416118640331778376D+00
    w(11) = 0.160021851762952142412820997988D+00
    w(12) = 0.116586655898711651540996670655D+00
    w(13) = 0.668372844976812846340706607461D-01
    w(14) = 0.109890109890109890109890109890D-01

  else if ( n == 15 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.965245926503838572795851392070D+00
    x(3) =  - 0.885082044222976298825401631482D+00
    x(4) =  - 0.763519689951815200704118475976D+00
    x(5) =  - 0.606253205469845711123529938637D+00
    x(6) =  - 0.420638054713672480921896938739D+00
    x(7) =  - 0.215353955363794238225679446273D+00
    x(8) =    0.0D+00
    x(9) =    0.215353955363794238225679446273D+00
    x(10) =   0.420638054713672480921896938739D+00
    x(11) =   0.606253205469845711123529938637D+00
    x(12) =   0.763519689951815200704118475976D+00
    x(13) =   0.885082044222976298825401631482D+00
    x(14) =   0.965245926503838572795851392070D+00
    x(15) =   1.0D+00

    w(1) =  0.952380952380952380952380952381D-02
    w(2) =  0.580298930286012490968805840253D-01
    w(3) =  0.101660070325718067603666170789D+00
    w(4) =  0.140511699802428109460446805644D+00
    w(5) =  0.172789647253600949052077099408D+00
    w(6) =  0.196987235964613356092500346507D+00
    w(7) =  0.211973585926820920127430076977D+00
    w(8) =  0.217048116348815649514950214251D+00
    w(9) =  0.211973585926820920127430076977D+00
    w(10) = 0.196987235964613356092500346507D+00
    w(11) = 0.172789647253600949052077099408D+00
    w(12) = 0.140511699802428109460446805644D+00
    w(13) = 0.101660070325718067603666170789D+00
    w(14) = 0.580298930286012490968805840253D-01
    w(15) = 0.952380952380952380952380952381D-02

  else if ( n == 16 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.969568046270217932952242738367D+00
    x(3) =  - 0.899200533093472092994628261520D+00
    x(4) =  - 0.792008291861815063931088270963D+00
    x(5) =  - 0.652388702882493089467883219641D+00
    x(6) =  - 0.486059421887137611781890785847D+00
    x(7) =  - 0.299830468900763208098353454722D+00
    x(8) =  - 0.101326273521949447843033005046D+00
    x(9) =    0.101326273521949447843033005046D+00
    x(10) =   0.299830468900763208098353454722D+00
    x(11) =   0.486059421887137611781890785847D+00
    x(12) =   0.652388702882493089467883219641D+00
    x(13) =   0.792008291861815063931088270963D+00
    x(14) =   0.899200533093472092994628261520D+00
    x(15) =   0.969568046270217932952242738367D+00
    x(16) =   1.0D+00

    w(1) =  0.833333333333333333333333333333D-02
    w(2) =  0.508503610059199054032449195655D-01
    w(3) =  0.893936973259308009910520801661D-01
    w(4) =  0.124255382132514098349536332657D+00
    w(5) =  0.154026980807164280815644940485D+00
    w(6) =  0.177491913391704125301075669528D+00
    w(7) =  0.193690023825203584316913598854D+00
    w(8) =  0.201958308178229871489199125411D+00
    w(9) =  0.201958308178229871489199125411D+00
    w(10) = 0.193690023825203584316913598854D+00
    w(11) = 0.177491913391704125301075669528D+00
    w(12) = 0.154026980807164280815644940485D+00
    w(13) = 0.124255382132514098349536332657D+00
    w(14) = 0.893936973259308009910520801661D-01
    w(15) = 0.508503610059199054032449195655D-01
    w(16) = 0.833333333333333333333333333333D-02

  else if ( n == 17 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.973132176631418314156979501874D+00
    x(3) =  - 0.910879995915573595623802506398D+00
    x(4) =  - 0.815696251221770307106750553238D+00
    x(5) =  - 0.691028980627684705394919357372D+00
    x(6) =  - 0.541385399330101539123733407504D+00
    x(7) =  - 0.372174433565477041907234680735D+00
    x(8) =  - 0.189511973518317388304263014753D+00
    x(9) =    0.0D+00
    x(10) =   0.189511973518317388304263014753D+00
    x(11) =   0.372174433565477041907234680735D+00
    x(12) =   0.541385399330101539123733407504D+00
    x(13) =   0.691028980627684705394919357372D+00
    x(14) =   0.815696251221770307106750553238D+00
    x(15) =   0.910879995915573595623802506398D+00
    x(16) =   0.973132176631418314156979501874D+00
    x(17) =   1.0D+00

    w(1) =  0.735294117647058823529411764706D-02
    w(2) =  0.449219405432542096474009546232D-01
    w(3) =  0.791982705036871191902644299528D-01
    w(4) =  0.110592909007028161375772705220D+00
    w(5) =  0.137987746201926559056201574954D+00
    w(6) =  0.160394661997621539516328365865D+00
    w(7) =  0.177004253515657870436945745363D+00
    w(8) =  0.187216339677619235892088482861D+00
    w(9) =  0.190661874753469433299407247028D+00
    w(10) = 0.187216339677619235892088482861D+00
    w(11) = 0.177004253515657870436945745363D+00
    w(12) = 0.160394661997621539516328365865D+00
    w(13) = 0.137987746201926559056201574954D+00
    w(14) = 0.110592909007028161375772705220D+00
    w(15) = 0.791982705036871191902644299528D-01
    w(16) = 0.449219405432542096474009546232D-01
    w(17) = 0.735294117647058823529411764706D-02

  else if ( n == 18 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.976105557412198542864518924342D+00
    x(3) =  - 0.920649185347533873837854625431D+00
    x(4) =  - 0.835593535218090213713646362328D+00
    x(5) =  - 0.723679329283242681306210365302D+00
    x(6) =  - 0.588504834318661761173535893194D+00
    x(7) =  - 0.434415036912123975342287136741D+00
    x(8) =  - 0.266362652878280984167665332026D+00
    x(9) =  - 0.897490934846521110226450100886D-01
    x(10) =   0.897490934846521110226450100886D-01
    x(11) =   0.266362652878280984167665332026D+00
    x(12) =   0.434415036912123975342287136741D+00
    x(13) =   0.588504834318661761173535893194D+00
    x(14) =   0.723679329283242681306210365302D+00
    x(15) =   0.835593535218090213713646362328D+00
    x(16) =   0.920649185347533873837854625431D+00
    x(17) =   0.976105557412198542864518924342D+00
    x(18) =   1.0D+00

    w(1) =  0.653594771241830065359477124183D-02
    w(2) =  0.399706288109140661375991764101D-01
    w(3) =  0.706371668856336649992229601678D-01
    w(4) =  0.990162717175028023944236053187D-01
    w(5) =  0.124210533132967100263396358897D+00
    w(6) =  0.145411961573802267983003210494D+00
    w(7) =  0.161939517237602489264326706700D+00
    w(8) =  0.173262109489456226010614403827D+00
    w(9) =  0.179015863439703082293818806944D+00
    w(10) = 0.179015863439703082293818806944D+00
    w(11) = 0.173262109489456226010614403827D+00
    w(12) = 0.161939517237602489264326706700D+00
    w(13) = 0.145411961573802267983003210494D+00
    w(14) = 0.124210533132967100263396358897D+00
    w(15) = 0.990162717175028023944236053187D-01
    w(16) = 0.706371668856336649992229601678D-01
    w(17) = 0.399706288109140661375991764101D-01
    w(18) = 0.653594771241830065359477124183D-02

  else if ( n == 19 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.978611766222080095152634063110D+00
    x(3) =  - 0.928901528152586243717940258797D+00
    x(4) =  - 0.852460577796646093085955970041D+00
    x(5) =  - 0.751494202552613014163637489634D+00
    x(6) =  - 0.628908137265220497766832306229D+00
    x(7) =  - 0.488229285680713502777909637625D+00
    x(8) =  - 0.333504847824498610298500103845D+00
    x(9) =  - 0.169186023409281571375154153445D+00
    x(10) =   0.0D+00
    x(11) =   0.169186023409281571375154153445D+00
    x(12) =   0.333504847824498610298500103845D+00
    x(13) =   0.488229285680713502777909637625D+00
    x(14) =   0.628908137265220497766832306229D+00
    x(15) =   0.751494202552613014163637489634D+00
    x(16) =   0.852460577796646093085955970041D+00
    x(17) =   0.928901528152586243717940258797D+00
    x(18) =   0.978611766222080095152634063110D+00
    x(19) =   1.0D+00

    w(1) =  0.584795321637426900584795321637D-02
    w(2) =  0.357933651861764771154255690351D-01
    w(3) =  0.633818917626297368516956904183D-01
    w(4) =  0.891317570992070844480087905562D-01
    w(5) =  0.112315341477305044070910015464D+00
    w(6) =  0.132267280448750776926046733910D+00
    w(7) =  0.148413942595938885009680643668D+00
    w(8) =  0.160290924044061241979910968184D+00
    w(9) =  0.167556584527142867270137277740D+00
    w(10) = 0.170001919284827234644672715617D+00
    w(11) = 0.167556584527142867270137277740D+00
    w(12) = 0.160290924044061241979910968184D+00
    w(13) = 0.148413942595938885009680643668D+00
    w(14) = 0.132267280448750776926046733910D+00
    w(15) = 0.112315341477305044070910015464D+00
    w(16) = 0.891317570992070844480087905562D-01
    w(17) = 0.633818917626297368516956904183D-01
    w(18) = 0.357933651861764771154255690351D-01
    w(19) = 0.584795321637426900584795321637D-02

  else if ( n == 20 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.980743704893914171925446438584D+00
    x(3) =  - 0.935934498812665435716181584931D+00
    x(4) =  - 0.866877978089950141309847214616D+00
    x(5) =  - 0.775368260952055870414317527595D+00
    x(6) =  - 0.663776402290311289846403322971D+00
    x(7) =  - 0.534992864031886261648135961829D+00
    x(8) =  - 0.392353183713909299386474703816D+00
    x(9) =  - 0.239551705922986495182401356927D+00
    x(10) = - 0.805459372388218379759445181596D-01
    x(11) =   0.805459372388218379759445181596D-01
    x(12) =   0.239551705922986495182401356927D+00
    x(13) =   0.392353183713909299386474703816D+00
    x(14) =   0.534992864031886261648135961829D+00
    x(15) =   0.663776402290311289846403322971D+00
    x(16) =   0.775368260952055870414317527595D+00
    x(17) =   0.866877978089950141309847214616D+00
    x(18) =   0.935934498812665435716181584931D+00
    x(19) =   0.980743704893914171925446438584D+00
    x(20) =   1.0D+00

    w(1) =  0.526315789473684210526315789474D-02
    w(2) =  0.322371231884889414916050281173D-01
    w(3) =  0.571818021275668260047536271732D-01
    w(4) =  0.806317639961196031447768461137D-01
    w(5) =  0.101991499699450815683781205733D+00
    w(6) =  0.120709227628674725099429705002D+00
    w(7) =  0.136300482358724184489780792989D+00
    w(8) =  0.148361554070916825814713013734D+00
    w(9) =  0.156580102647475487158169896794D+00
    w(10) = 0.160743286387845749007726726449D+00
    w(11) = 0.160743286387845749007726726449D+00
    w(12) = 0.156580102647475487158169896794D+00
    w(13) = 0.148361554070916825814713013734D+00
    w(14) = 0.136300482358724184489780792989D+00
    w(15) = 0.120709227628674725099429705002D+00
    w(16) = 0.101991499699450815683781205733D+00
    w(17) = 0.806317639961196031447768461137D-01
    w(18) = 0.571818021275668260047536271732D-01
    w(19) = 0.322371231884889414916050281173D-01
    w(20) = 0.526315789473684210526315789474D-02

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOBATTO_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are between 2 and 20.'
    stop

  end if

  return
end
subroutine moulton_set ( n, x, w )

!*****************************************************************************80
!
!! MOULTON_SET sets weights for Adams-Moulton quadrature.
!
!  Discussion:
!
!    Adams-Moulton quadrature formulas are normally used in solving
!    ordinary differential equations, and are not suitable for general
!    quadrature computations.  However, an Adams-Moulton formula is
!    equivalent to approximating the integral of F(Y(X)) between X(M)
!    and X(M+1), using an implicit formula that relies on known values
!    of F(Y(X)) at X(M-N+1) through X(M), plus the unknown value at X(M+1).
!
!    Suppose the unknown function is denoted by Y(X), with derivative F(Y(X)),
!    and that approximate values of the function are known at a series of
!    X values, which we write as X(1), X(2), ..., X(M).  We write the value
!    Y(X(1)) as Y(1) and so on.
!
!    Then the solution of the ODE Y' = F(X,Y) at the next point X(M+1) is
!    computed by:
!
!      Y(M+1) = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dx
!             = Y(M) + H * Sum ( 1 <= I <= N ) W(I) * F(Y(M+2-I)) approximately.
!
!    Note that this formula is implicit, since the unknown value Y(M+1)
!    appears on the right hand side.  Hence, in ODE applications, this
!    equation must be solved via a nonlinear equation solver.  For
!    quadrature problems, where the function to be integrated is known
!    beforehand, this is not a problem, and the calculation is explicit.
!
!    In the documentation that follows, we replace F(Y(X)) by F(X).
!
!
!    The Adams-Moulton formulas require equally spaced data.
!
!    Here is how the formula is applied in the case with non-unit spacing:
!
!      Integral ( A <= X <= A+H ) F(X) dx =
!      H * Sum ( 1 <= I <= N ) W(I) * F ( A - (I-2)*H ),
!      approximately.
!
!    The integral:
!
!      Integral ( 0 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( 2 - I )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Leon Lapidus, John Seinfeld,
!    Numerical Solution of Ordinary Differential Equations,
!    Academic Press, 1971,
!    ISBN: 0124366503.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 10, 12, 14, 16, 18 or 20.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!    W(1) is the weight at X = 1, W(2) the weight at X = 0, and so on.
!    The weights are rational.  The weights are not symmetric, and
!    some weights may be negative.  They should sum to 1.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    w(1) =  1.0D+00

  else if ( n == 2 ) then

    d = 2.0D+00

    w(1) =  1.0D+00 / d
    w(2) =  1.0D+00 / d

  else if ( n == 3 ) then

    d = 12.0D+00

    w(1) =    5.0D+00 / d
    w(2) =    8.0D+00 / d
    w(3) =  - 1.0D+00 / d

  else if ( n == 4 ) then

    d = 24.0D+00

    w(1) =    9.0D+00 / d
    w(2) =   19.0D+00 / d
    w(3) =  - 5.0D+00 / d
    w(4) =    1.0D+00 / d

  else if ( n == 5 ) then

    d = 720.0D+00

    w(1) =    251.0D+00 / d
    w(2) =    646.0D+00 / d
    w(3) =  - 264.0D+00 / d
    w(4) =    106.0D+00 / d
    w(5) =   - 19.0D+00 / d

  else if ( n == 6 ) then

    d = 1440.0D+00

    w(1) =    475.0D+00 / d
    w(2) =   1427.0D+00 / d
    w(3) =  - 798.0D+00 / d
    w(4) =    482.0D+00 / d
    w(5) =  - 173.0D+00 / d
    w(6) =     27.0D+00 / d

  else if ( n == 7 ) then

    d = 60480.0D+00

    w(1) =    19087.0D+00 / d
    w(2) =    65112.0D+00 / d
    w(3) =  - 46461.0D+00 / d
    w(4) =    37504.0D+00 / d
    w(5) =  - 20211.0D+00 / d
    w(6) =     6312.0D+00 / d
    w(7) =    - 863.0D+00 / d

  else if ( n == 8 ) then

    d = 120960.0D+00

    w(1) =    36799.0D+00 / d
    w(2) =   139849.0D+00 / d
    w(3) = - 121797.0D+00 / d
    w(4) =   123133.0D+00 / d
    w(5) =  - 88547.0D+00 / d
    w(6) =    41499.0D+00 / d
    w(7) =  - 11351.0D+00 / d
    w(8) =     1375.0D+00 / d

  else if ( n == 9 ) then

    d = 3628800.0D+00

    w(1) =   1070017.0D+00 / d
    w(2) =   4467094.0D+00 / d
    w(3) = - 4604594.0D+00 / d
    w(4) =   5595358.0D+00 / d
    w(5) = - 5033120.0D+00 / d
    w(6) =   3146338.0D+00 / d
    w(7) = - 1291214.0D+00 / d
    w(8) =    312874.0D+00 / d
    w(9) =   - 33953.0D+00 / d

  else if ( n == 10 ) then

    d = 7257600.0D+00

    w(1) =    2082753.0D+00 / d
    w(2) =    9449717.0D+00 / d
    w(3) = - 11271304.0D+00 / d
    w(4) =   16002320.0D+00 / d
    w(5) = - 17283646.0D+00 / d
    w(6) =   13510082.0D+00 / d
    w(7) =  - 7394032.0D+00 / d
    w(8) =    2687864.0D+00 / d
    w(9) =   - 583435.0D+00 / d
    w(10) =     57281.0D+00 / d

  else if ( n == 12 ) then

    d = 958003200.0D+00

    w(1) =    262747265.0D+00 / d
    w(2) =   1374799219.0D+00 / d
    w(3) =  -2092490673.0D+00 / d
    w(4) =   3828828885.0D+00 / d
    w(5) =  -5519460582.0D+00 / d
    w(6) =   6043521486.0D+00 / d
    w(7) =  -4963166514.0D+00 / d
    w(8) =   3007739418.0D+00 / d
    w(9) =  -1305971115.0D+00 / d
    w(10) =   384709327.0D+00 / d
    w(11) =   -68928781.0D+00 / d
    w(12) =     5675265.0D+00 / d

  else if ( n == 14 ) then

    d = 5230697472000.0D+00

    w(1) =    1382741929621.0D+00 / d
    w(2) =    8153167962181.0D+00 / d
    w(3) =  -15141235084110.0D+00 / d
    w(4) =   33928990133618.0D+00 / d
    w(5) =  -61188680131285.0D+00 / d
    w(6) =   86180228689563.0D+00 / d
    w(7) =  -94393338653892.0D+00 / d
    w(8) =   80101021029180.0D+00 / d
    w(9) =  -52177910882661.0D+00 / d
    w(10) =  25620259777835.0D+00 / d
    w(11) =  -9181635605134.0D+00 / d
    w(12) =   2268078814386.0D+00 / d
    w(13) =   -345457086395.0D+00 / d
    w(14) =     24466579093.0D+00 / d

  else if ( n == 16 ) then

    d = 62768369664000.0D+00

    w(1) =     16088129229375.0D+00 / d
    w(2) =    105145058757073.0D+00 / d
    w(3) =   -230992163723849.0D+00 / d
    w(4) =    612744541065337.0D+00 / d
    w(5) =  -1326978663058069.0D+00 / d
    w(6) =   2285168598349733.0D+00 / d
    w(7) =  -3129453071993581.0D+00 / d
    w(8) =   3414941728852893.0D+00 / d
    w(9) =  -2966365730265699.0D+00 / d
    w(10) =  2039345879546643.0D+00 / d
    w(11) = -1096355235402331.0D+00 / d
    w(12) =   451403108933483.0D+00 / d
    w(13) =  -137515713789319.0D+00 / d
    w(14) =    29219384284087.0D+00 / d
    w(15) =    -3867689367599.0D+00 / d
    w(16) =      240208245823.0D+00 / d

  else if ( n == 18 ) then

    d = 64023737057280000.0D+00

    w(1) =      15980174332775873.0D+00 / d
    w(2) =     114329243705491117.0D+00 / d
    w(3) =    -290470969929371220.0D+00 / d
    w(4) =     890337710266029860.0D+00 / d
    w(5) =   -2250854333681641520.0D+00 / d
    w(6) =    4582441343348851896.0D+00 / d
    w(7) =   -7532171919277411636.0D+00 / d
    w(8) =   10047287575124288740.0D+00 / d
    w(9) =  -10910555637627652470.0D+00 / d
    w(10) =   9644799218032932490.0D+00 / d
    w(11) =  -6913858539337636636.0D+00 / d
    w(12) =   3985516155854664396.0D+00 / d
    w(13) =  -1821304040326216520.0D+00 / d
    w(14) =    645008976643217360.0D+00 / d
    w(15) =   -170761422500096220.0D+00 / d
    w(16) =     31816981024600492.0D+00 / d
    w(17) =     -3722582669836627.0D+00 / d
    w(18) =       205804074290625.0D+00 / d

  else if ( n == 20 ) then

    d = 102181884343418880000.0D+00

    w(1) =       24919383499187492303.0D+00 / d
    w(2) =      193280569173472261637.0D+00 / d
    w(3) =     -558160720115629395555.0D+00 / d
    w(4) =     1941395668950986461335.0D+00 / d
    w(5) =    -5612131802364455926260.0D+00 / d
    w(6) =    13187185898439270330756.0D+00 / d
    w(7) =   -25293146116627869170796.0D+00 / d
    w(8) =    39878419226784442421820.0D+00 / d
    w(9) =   -51970649453670274135470.0D+00 / d
    w(10) =   56154678684618739939910.0D+00 / d
    w(11) =  -50320851025594566473146.0D+00 / d
    w(12) =   37297227252822858381906.0D+00 / d
    w(13) =  -22726350407538133839300.0D+00 / d
    w(14) =   11268210124987992327060.0D+00 / d
    w(15) =   -4474886658024166985340.0D+00 / d
    w(16) =    1389665263296211699212.0D+00 / d
    w(17) =    -325187970422032795497.0D+00 / d
    w(18) =      53935307402575440285.0D+00 / d
    w(19) =      -5652892248087175675.0D+00 / d
    w(20) =        281550972898020815.0D+00 / d

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOULTON_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 through 10, 12, 14, 16, 18 or 20.'
    stop

  end if

  do i = 1, n
    x(i) = real ( 2 - i, kind = 8 )
  end do

  return
end
subroutine nc_compute ( n, x_min, x_max, x, w )

!*****************************************************************************80
!
!! NC_COMPUTE computes a Newton-Cotes quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      Integral ( X_MIN <= X <= X_MAX ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
!
!    For the CLOSED rule, the abscissas include the end points.
!    For an OPEN rule, the abscissas do not include the end points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X_MIN, X_MAX, the endpoints of the interval.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  real ( kind = 8 ) yvala
  real ( kind = 8 ) yvalb

  do i = 1, n
!
!  Compute the Lagrange basis polynomial which is 1 at X(I),
!  and zero at the other nodes.
!
    d(1:n) = 0.0D+00
    d(i) = 1.0D+00

    do j = 2, n
      do k = j, n
        d(n+j-k) = ( d(n+j-k-1) - d(n+j-k) ) / ( x(n+1-k) - x(n+j-k) )
      end do
    end do

    do j = 1, n - 1
      do k = 1, n - j
        d(n-k) = d(n-k) - x(n-k-j+1) * d(n-k+1)
      end do
    end do
!
!  Evaluate the antiderivative of the polynomial at the endpoints.
!
    yvala = d(n) / real ( n, kind = 8 )
    do j = n - 1, 1, -1
      yvala = yvala * x_min + d(j) / real ( j, kind = 8 )
    end do
    yvala = yvala * x_min

    yvalb = d(n) / real ( n, kind = 8 )
    do j = n - 1, 1, -1
      yvalb = yvalb * x_max + d(j) / real ( j, kind = 8 )
    end do
    yvalb = yvalb * x_max

    w(i) = yvalb - yvala

  end do

  return
end
subroutine ncc_compute ( n, x, w )

!*****************************************************************************80
!
!! NCC_COMPUTE: Newton-Cotes Closed quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      Integral ( -1 <= X <= +1 ) F(X) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
!
!    For the CLOSED rule, the abscissas are equally spaced and include
!    the end points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call ncc_compute_points ( n, x )

  call ncc_compute_weights ( n, w )

  return
end
subroutine ncc_compute_points ( n, x )

!*****************************************************************************80
!
!! NCC_COMPUTE_POINTS: Newton-Cotes Closed points
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  if ( n == 1 ) then

    x(1) = ( x_max + x_min ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = 8 ) * x_min   &
             + real (     i - 1, kind = 8 ) * x_max ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine ncc_compute_weights ( n, w )

!*****************************************************************************80
!
!! NCC_COMPUTE_WEIGHTS: Newton-Cotes Closed weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  if ( n == 1 ) then

    w(1) = x_max - x_min

  else

    call ncc_compute_points ( n, x )

    call nc_compute ( n, x_min, x_max, x, w )

  end if

  return
end
subroutine ncc_set ( n, x, w )

!*****************************************************************************80
!
!! NCC_SET sets abscissas and weights for Newton-Cotes closed quadrature.
!
!  Discussion:
!
!    The Newton-Cotes closed rules use equally spaced abscissas, and
!    hence may be used with tabulated function data.
!
!    The rules are called "closed" because they include the endpoints.
!    As a favor, we include an order 1 rule, the midpoint rule, even
!    though this does not satisfy the requirement that the endpoints
!    be included!
!
!    The higher order rules involve negative weights.  These can produce
!    loss of accuracy due to the subtraction of large, nearly equal quantities.
!
!    N = 1 is the midpoint rule (and is not really an NCC rule!)
!    N = 2 is the trapezoidal rule.
!    N = 3 is Simpson's rule.
!    N = 4 is Simpson's 3/8 rule.
!    N = 5 is Bode's rule.
!
!    The Kopal reference for N = 12 lists
!      W(6) = 15494566.0D+00 / 43545600.0D+00
!    but this results in a set of coeffients that don't add up to 2.
!    The correct value is
!      W(6) = 15493566.0D+00 / 43545600.0.
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    In Mathematica, the Newton-Cotes closed weights of order N in
!    the interval [-1,+1] can be computed by:
!
!      Needs["NumericalDifferentialEquationAnalysis`"]
!      NewtonCotesWeights [ n, -1, 1, QuadratureType -> Closed ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Johnson,
!    Quarterly Journal of Mathematics,
!    Volume 46, Number 52, 1915.
!
!    Zdenek Kopal,
!    Numerical Analysis,
!    John Wiley, 1955,
!    LC: QA297.K6.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N <= 21.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then
!
!  2
!
    x(1) = 0.00000000000000000000D+00

    w(1) = 2.00000000000000000000D+00

  else if ( n == 2 ) then
!
!  1
!  1
!
    x(1) = -1.00000000000000000000D+00
    x(2) =  1.00000000000000000000D+00

    w(1) = 1.00000000000000000000D+00
    w(2) = 1.00000000000000000000D+00

  else if ( n == 3 ) then
!
!  1 / 3
!  4 / 3
!  1 / 3
!
    x(1) = -1.00000000000000000000D+00
    x(2) =  0.00000000000000000000D+00
    x(3) =  1.00000000000000000000D+00

    w(1) = 0.33333333333333333333D+00
    w(2) = 1.33333333333333333333D+00
    w(3) = 0.33333333333333333333D+00

  else if ( n == 4 ) then
!
!  1 / 4
!  3 / 4
!  3 / 4
!  1 / 4
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.33333333333333333333D+00
    x(3) =  0.33333333333333333333D+00
    x(4) =  1.00000000000000000000D+00

    w(1) = 0.25000000000000000000D+00
    w(2) = 0.75000000000000000000D+00
    w(3) = 0.75000000000000000000D+00
    w(4) = 0.25000000000000000000D+00

  else if ( n == 5 ) then
!
!   7 / 45
!  32 / 45
!  12 / 45
!  32 / 45
!   7 / 45
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.50000000000000000000D+00
    x(3) =  0.00000000000000000000D+00
    x(4) =  0.50000000000000000000D+00
    x(5) =  1.00000000000000000000D+00

    w(1) = 0.15555555555555555556D+00
    w(2) = 0.71111111111111111111D+00
    w(3) = 0.26666666666666666667D+00
    w(4) = 0.71111111111111111111D+00
    w(5) = 0.15555555555555555556D+00

  else if ( n == 6 ) then
!
!  19 / 144
!  75 / 144
!  50 / 144
!  50 / 144
!  75 / 144
!  19 / 144
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.60000000000000000000D+00
    x(3) = -0.20000000000000000000D+00
    x(4) =  0.20000000000000000000D+00
    x(5) =  0.60000000000000000000D+00
    x(6) =  1.00000000000000000000D+00

    w(1) = 0.13194444444444444444D+00
    w(2) = 0.52083333333333333333D+00
    w(3) = 0.34722222222222222222D+00
    w(4) = 0.34722222222222222222D+00
    w(5) = 0.52083333333333333333D+00
    w(6) = 0.13194444444444444444D+00

  else if ( n == 7 ) then
!
!   41 / 420
!  216 / 420
!   27 / 420
!  272 / 420
!   27 / 420
!  216 / 420
!   41 / 420
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.66666666666666666667D+00
    x(3) = -0.33333333333333333333D+00
    x(4) =  0.00000000000000000000D+00
    x(5) =  0.33333333333333333333D+00
    x(6) =  0.66666666666666666667D+00
    x(7) =  1.00000000000000000000D+00

    w(1) = 0.097619047619047619048D+00
    w(2) = 0.51428571428571428571D+00
    w(3) = 0.064285714285714285714D+00
    w(4) = 0.64761904761904761905D+00
    w(5) = 0.064285714285714285714D+00
    w(6) = 0.51428571428571428571D+00
    w(7) = 0.097619047619047619048D+00

  else if ( n == 8 ) then
!
!   751 / 8640
!  3577 / 8640
!  1323 / 8640
!  2989 / 8640
!  2989 / 8640
!  1323 / 8640
!  3577 / 8640
!   751 / 8640
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.71428571428571428571D+00
    x(3) = -0.42857142857142857143D+00
    x(4) = -0.14285714285714285714D+00
    x(5) =  0.14285714285714285714D+00
    x(6) =  0.42857142857142857143D+00
    x(7) =  0.71428571428571428571D+00
    x(8) =  1.00000000000000000000D+00

    w(1) = 0.086921296296296296296D+00
    w(2) = 0.41400462962962962963D+00
    w(3) = 0.15312500000000000000D+00
    w(4) = 0.34594907407407407407D+00
    w(5) = 0.34594907407407407407D+00
    w(6) = 0.15312500000000000000D+00
    w(7) = 0.41400462962962962963D+00
    w(8) = 0.086921296296296296296D+00

  else if ( n == 9 ) then
!
!    989 / 14175
!   5888 / 14175
!   -928 / 14175
!  10496 / 14175
!  -4540 / 14175
!  10496 / 14175
!   -928 / 14175
!   5888 / 14175
!    989 / 14175
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.75000000000000000000D+00
    x(3) = -0.50000000000000000000D+00
    x(4) = -0.25000000000000000000D+00
    x(5) =  0.00000000000000000000D+00
    x(6) =  0.25000000000000000000D+00
    x(7) =  0.50000000000000000000D+00
    x(8) =  0.75000000000000000000D+00
    x(9) =  1.00000000000000000000D+00

    w(1) =  0.069770723104056437390D+00
    w(2) =  0.41537918871252204586D+00
    w(3) = -0.065467372134038800705D+00
    w(4) =  0.74045855379188712522D+00
    w(5) = -0.32028218694885361552D+00
    w(6) =  0.74045855379188712522D+00
    w(7) = -0.065467372134038800705D+00
    w(8) =  0.41537918871252204586D+00
    w(9) =  0.069770723104056437390D+00

  else if ( n == 10 ) then
!
!   2857 / 44800
!  15741 / 44800
!   1080 / 44800
!  19344 / 44800
!   5778 / 44800
!   5778 / 44800
!  19344 / 44800
!   1080 / 44800
!  15741 / 44800
!   2857 / 44800
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.77777777777777777778D+00
    x(3) = -0.55555555555555555556D+00
    x(4) = -0.33333333333333333333D+00
    x(5) = -0.11111111111111111111D+00
    x(6) =  0.11111111111111111111D+00
    x(7) =  0.33333333333333333333D+00
    x(8) =  0.55555555555555555556D+00
    x(9) =  0.77777777777777777778D+00
    x(10) = 1.00000000000000000000D+00

    w(1) =  0.063772321428571428571D+00
    w(2) =  0.35136160714285714286D+00
    w(3) =  0.024107142857142857143D+00
    w(4) =  0.43178571428571428571D+00
    w(5) =  0.12897321428571428571D+00
    w(6) =  0.12897321428571428571D+00
    w(7) =  0.43178571428571428571D+00
    w(8) =  0.024107142857142857143D+00
    w(9) =  0.35136160714285714286D+00
    w(10) = 0.063772321428571428571D+00

  else if ( n == 11 ) then
!
!     16067 / 299376
!    106300 / 299376
!   - 48525 / 299376
!    272400 / 299376
!  - 260550 / 299376
!    427368 / 299376
!  - 260550 / 299376
!    272400 / 299376
!   - 48525 / 299376
!    106300 / 299376
!     16067 / 299376
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.80000000000000000000D+00
    x(3) = -0.60000000000000000000D+00
    x(4) = -0.40000000000000000000D+00
    x(5) = -0.20000000000000000000D+00
    x(6) =  0.00000000000000000000D+00
    x(7) =  0.20000000000000000000D+00
    x(8) =  0.40000000000000000000D+00
    x(9) =  0.60000000000000000000D+00
    x(10) = 0.80000000000000000000D+00
    x(11) = 1.00000000000000000000D+00

    w(1) =  0.053668296723852279408D+00
    w(2) =  0.35507188284966062744D+00
    w(3) = -0.16208714125380792047D+00
    w(4) =  0.90989257655924322591D+00
    w(5) = -0.87031024531024531025D+00
    w(6) =  1.4275292608625941959D+00
    w(7) = -0.87031024531024531025D+00
    w(8) =  0.90989257655924322591D+00
    w(9) = -0.16208714125380792047D+00
    w(10) = 0.35507188284966062744D+00
    w(11) = 0.053668296723852279408D+00

  else if ( n == 12 ) then
!
!     2171465 / 43545600
!    13486539 / 43545600
!   - 3237113 / 43545600
!    25226685 / 43545600
!   - 9595542 / 43545600
!    15493566 / 43545600
!    15493566 / 43545600
!   - 9595542 / 43545600
!    25226685 / 43545600
!   - 3237113 / 43545600
!    13486539 / 43545600
!     2171465 / 43545600
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.81818181818181818182D+00
    x(3) = -0.63636363636363636364D+00
    x(4) = -0.45454545454545454545D+00
    x(5) = -0.27272727272727272727D+00
    x(6) = -0.090909090909090909091D+00
    x(7) =  0.090909090909090909091D+00
    x(8) =  0.27272727272727272727D+00
    x(9) =  0.45454545454545454545D+00
    x(10) = 0.63636363636363636364D+00
    x(11) = 0.81818181818181818182D+00
    x(12) = 1.00000000000000000000D+00

    w(1) =   0.049866461823927101705D+00
    w(2) =   0.30971071704144620811D+00
    w(3) =  -0.074338463587595532040D+00
    w(4) =   0.57931650958994708995D+00
    w(5) =  -0.22035617835097001764D+00
    w(6) =   0.35580095348324514991D+00
    w(7) =   0.35580095348324514991D+00
    w(8) =  -0.22035617835097001764D+00
    w(9) =   0.57931650958994708995D+00
    w(10) = -0.074338463587595532040D+00
    w(11) =  0.30971071704144620811D+00
    w(12) =  0.049866461823927101705D+00

  else if ( n == 13 ) then
!
!      1364651 / 31531500
!      9903168 / 31531500
!    - 7587864 / 31531500
!     35725120 / 31531500
!   - 51491295 / 31531500
!     87516288 / 31531500
!   - 87797136 / 31531500
!     87516288 / 31531500
!   - 51491295 / 31531500
!     35725120 / 31531500
!    - 7587864 / 31531500
!      9903168 / 31531500
!      1364651 / 31531500
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.83333333333333333333D+00
    x(3) = -0.66666666666666666667D+00
    x(4) = -0.50000000000000000000D+00
    x(5) = -0.33333333333333333333D+00
    x(6) = -0.16666666666666666667D+00
    x(7) =  0.00000000000000000000D+00
    x(8) =  0.16666666666666666667D+00
    x(9) =  0.33333333333333333333D+00
    x(10) = 0.50000000000000000000D+00
    x(11) = 0.66666666666666666667D+00
    x(12) = 0.83333333333333333333D+00
    x(13) = 1.00000000000000000000D+00

    w(1) =   0.043278974993260707546D+00
    w(2) =   0.31407221350078492936D+00
    w(3) =  -0.24064392750107035821D+00
    w(4) =   1.1329977958549387121D+00
    w(5) =  -1.6330112744398458684D+00
    w(6) =   2.7755193378050520908D+00
    w(7) =  -2.7844262404262404262D+00
    w(8) =   2.7755193378050520908D+00
    w(9) =  -1.6330112744398458684D+00
    w(10) =  1.1329977958549387121D+00
    w(11) = -0.24064392750107035821D+00
    w(12) =  0.31407221350078492936D+00
    w(13) =  0.043278974993260707546D+00

  else if ( n == 14 ) then
!
!      6137698213 / 150885504000
!     42194238652 / 150885504000
!   - 23361540993 / 150885504000
!    116778274403 / 150885504000
!  - 113219777650 / 150885504000
!    154424590209 / 150885504000
!   - 32067978834 / 150885504000
!   - 32067978834 / 150885504000
!    154424590209 / 150885504000
!  - 113219777650 / 150885504000
!    116778274403 / 150885504000
!   - 23361540993 / 150885504000
!     42194238652 / 150885504000
!      6137698213 / 150885504000
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.84615384615384615385D+00
    x(3) = -0.69230769230769230769D+00
    x(4) = -0.53846153846153846154D+00
    x(5) = -0.38461538461538461538D+00
    x(6) = -0.23076923076923076923D+00
    x(7) = -0.076923076923076923077D+00
    x(8) =  0.076923076923076923077D+00
    x(9) =  0.23076923076923076923D+00
    x(10) = 0.38461538461538461538D+00
    x(11) = 0.53846153846153846154D+00
    x(12) = 0.69230769230769230769D+00
    x(13) = 0.84615384615384615385D+00
    x(14) = 1.00000000000000000000D+00

    w(1) =   0.040669438210247155353D+00
    w(2) =   0.27975217053157074652D+00
    w(3) =  -0.15542374057682837445D+00
    w(4) =   0.77579230848776566369D+00
    w(5) =  -0.75384763266423526013D+00
    w(6) =   1.0273523591123107492D+00
    w(7) =  -0.21429490310083068020D+00
    w(8) =  -0.21429490310083068020D+00
    w(9) =   1.0273523591123107492D+00
    w(10) = -0.75384763266423526013D+00
    w(11) =  0.77579230848776566369D+00
    w(12) = -0.15542374057682837445D+00
    w(13) =  0.27975217053157074652D+00
    w(14) =  0.040669438210247155353D+00

  else if ( n == 15 ) then
!
!       90241897 / 2501928000
!      710986864 / 2501928000
!    - 770720657 / 2501928000
!     3501442784 / 2501928000
!   - 6625093363 / 2501928000
!    12630121616 / 2501928000
!  - 16802270373 / 2501928000
!    19534438464 / 2501928000
!  - 16802270373 / 2501928000
!    12630121616 / 2501928000
!   - 6625093363 / 2501928000
!     3501442784 / 2501928000
!    - 770720657 / 2501928000
!      710986864 / 2501928000
!       90241897 / 2501928000
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.85714285714285714286D+00
    x(3) = -0.71428571428571428571D+00
    x(4) = -0.57142857142857142857D+00
    x(5) = -0.42857142857142857143D+00
    x(6) = -0.28571428571428571429D+00
    x(7) = -0.14285714285714285714D+00
    x(8) =  0.00000000000000000000D+00
    x(9) =  0.14285714285714285714D+00
    x(10) = 0.28571428571428571429D+00
    x(11) = 0.42857142857142857143D+00
    x(12) = 0.57142857142857142857D+00
    x(13) = 0.71428571428571428571D+00
    x(14) = 0.85714285714285714286D+00
    x(15) = 1.00000000000000000000D+00

    w(1) =   0.036068942431596752584D+00
    w(2) =   0.28417558938546592868D+00
    w(3) =  -0.30805069410470645039D+00
    w(4) =   1.3994978208805369299D+00
    w(5) =  -2.6479952112930507992D+00
    w(6) =   5.0481555088715582543D+00
    w(7) =  -6.7157289790113864188D+00
    w(8) =   7.8077540456799716059D+00
    w(9) =  -6.7157289790113864188D+00
    w(10) =  5.0481555088715582543D+00
    w(11) = -2.6479952112930507992D+00
    w(12) =  1.3994978208805369299D+00
    w(13) = -0.30805069410470645039D+00
    w(14) =  0.28417558938546592868D+00
    w(15) =  0.036068942431596752584D+00

  else if ( n == 16 ) then
!
!     105930069 / 3099672576
!     796661595 / 3099672576
!   - 698808195 / 3099672576
!    3143332755 / 3099672576
!  - 4688522055 / 3099672576
!    7385654007 / 3099672576
!  - 6000998415 / 3099672576
!    3056422815 / 3099672576
!    3056422815 / 3099672576
!  - 6000998415 / 3099672576
!    7385654007 / 3099672576
!  - 4688522055 / 3099672576
!    3143332755 / 3099672576
!   - 698808195 / 3099672576
!     796661595 / 3099672576
!     105930069 / 3099672576
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.86666666666666666667D+00
    x(3) = -0.73333333333333333333D+00
    x(4) = -0.60000000000000000000D+00
    x(5) = -0.46666666666666666667D+00
    x(6) = -0.33333333333333333333D+00
    x(7) = -0.20000000000000000000D+00
    x(8) = -0.066666666666666666667D+00
    x(9) =  0.066666666666666666667D+00
    x(10) = 0.20000000000000000000D+00
    x(11) = 0.33333333333333333333D+00
    x(12) = 0.46666666666666666667D+00
    x(13) = 0.60000000000000000000D+00
    x(14) = 0.73333333333333333333D+00
    x(15) = 0.86666666666666666667D+00
    x(16) = 1.00000000000000000000D+00

    w(1) =   0.034174599543251887002D+00
    w(2) =   0.25701475735481036820D+00
    w(3) =  -0.22544581011901045383D+00
    w(4) =   1.0140854164204471124D+00
    w(5) =  -1.5125862296882804695D+00
    w(6) =   2.3827206990135980091D+00
    w(7) =  -1.9360104229924960952D+00
    w(8) =   0.98604699046767964179D+00
    w(9) =   0.98604699046767964179D+00
    w(10) = -1.9360104229924960952D+00
    w(11) =  2.3827206990135980091D+00
    w(12) = -1.5125862296882804695D+00
    w(13) =  1.0140854164204471124D+00
    w(14) = -0.22544581011901045383D+00
    w(15) =  0.25701475735481036820D+00
    w(16) =  0.034174599543251887002D+00

  else if ( n == 17 ) then
!
!       15043611773 / 488462349375
!      127626606592 / 488462349375
!    - 179731134720 / 488462349375
!      832211855360 / 488462349375
!   - 1929498607520 / 488462349375
!     4177588893696 / 488462349375
!   - 6806534407936 / 488462349375
!     9368875018240 / 488462349375
!  - 10234238972220 / 488462349375
!     9368875018240 / 488462349375
!   - 6806534407936 / 488462349375
!     4177588893696 / 488462349375
!   - 1929498607520 / 488462349375
!      832211855360 / 488462349375
!    - 179731134720 / 488462349375
!      127626606592 / 488462349375
!       15043611773 / 488462349375
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.87500000000000000000D+00
    x(3) = -0.75000000000000000000D+00
    x(4) = -0.62500000000000000000D+00
    x(5) = -0.50000000000000000000D+00
    x(6) = -0.37500000000000000000D+00
    x(7) = -0.25000000000000000000D+00
    x(8) = -0.12500000000000000000D+00
    x(9) =  0.00000000000000000000D+00
    x(10) = 0.12500000000000000000D+00
    x(11) = 0.25000000000000000000D+00
    x(12) = 0.37500000000000000000D+00
    x(13) = 0.50000000000000000000D+00
    x(14) = 0.62500000000000000000D+00
    x(15) = 0.75000000000000000000D+00
    x(16) = 0.87500000000000000000D+00
    x(17) = 1.00000000000000000000D+00

    w(1)  =   0.030797894233299012495D+00
    w(2)  =   0.26128238288028031086D+00
    w(3)  =  -0.36795289329867605622D+00
    w(4)  =   1.7037379778090086905D+00
    w(5)  =  -3.9501480717783930427D+00
    w(6)  =   8.5525299934402953388D+00
    w(7)  = -13.934614237197880038D+00
    w(8)  =  19.180342211078732848D+00
    w(9)  = -20.951950514333334128D+00
    w(10) =  19.180342211078732848D+00
    w(11) = -13.934614237197880038D+00
    w(12) =   8.5525299934402953388D+00
    w(13) =  -3.9501480717783930427D+00
    w(14) =   1.7037379778090086905D+00
    w(15) =  -0.36795289329867605622D+00
    w(16) =   0.26128238288028031086D+00
    w(17) =   0.030797894233299012495D+00

  else if ( n == 18 ) then
!
!       55294720874657 / 1883051089920000
!      450185515446285 / 1883051089920000
!    - 542023437008852 / 1883051089920000
!     2428636525764260 / 1883051089920000
!   - 4768916800123440 / 1883051089920000
!     8855416648684984 / 1883051089920000
!  - 10905371859796660 / 1883051089920000
!    10069615750132836 / 1883051089920000
!   - 3759785974054070 / 1883051089920000
!   - 3759785974054070 / 1883051089920000
!    10069615750132836 / 1883051089920000
!  - 10905371859796660 / 1883051089920000
!     8855416648684984 / 1883051089920000
!   - 4768916800123440 / 1883051089920000
!     2428636525764260 / 1883051089920000
!    - 542023437008852 / 1883051089920000
!      450185515446285 / 1883051089920000
!       55294720874657 / 1883051089920000
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.88235294117647058824D+00
    x(3) = -0.76470588235294117647D+00
    x(4) = -0.64705882352941176471D+00
    x(5) = -0.52941176470588235294D+00
    x(6) = -0.41176470588235294118D+00
    x(7) = -0.29411764705882352941D+00
    x(8) = -0.17647058823529411765D+00
    x(9) = -0.058823529411764705882D+00
    x(10) = 0.058823529411764705882D+00
    x(11) = 0.17647058823529411765D+00
    x(12) = 0.29411764705882352941D+00
    x(13) = 0.41176470588235294118D+00
    x(14) = 0.52941176470588235294D+00
    x(15) = 0.64705882352941176471D+00
    x(16) = 0.76470588235294117647D+00
    x(17) = 0.88235294117647058824D+00
    x(18) = 1.00000000000000000000D+00

    w(1) =   0.029364429446790078519D+00
    w(2) =   0.23907238516051669677D+00
    w(3) =  -0.28784319231183443641D+00
    w(4) =   1.2897348026109258587D+00
    w(5) =  -2.5325477495812627261D+00
    w(6) =   4.7026959045817496499D+00
    w(7) =  -5.7913308450170443690D+00
    w(8) =   5.3475000248456540826D+00
    w(9) =  -1.9966457597354948350D+00
    w(10) = -1.9966457597354948350D+00
    w(11) =  5.3475000248456540826D+00
    w(12) = -5.7913308450170443690D+00
    w(13) =  4.7026959045817496499D+00
    w(14) = -2.5325477495812627261D+00
    w(15) =  1.2897348026109258587D+00
    w(16) = -0.28784319231183443641D+00
    w(17) =  0.23907238516051669677D+00
    w(18) =  0.029364429446790078519D+00

  else if ( n == 19 ) then
!
!       203732352169 / 7604556960000
!      1848730221900 / 7604556960000
!    - 3212744374395 / 7604556960000
!     15529830312096 / 7604556960000
!   - 42368630685840 / 7604556960000
!    103680563465808 / 7604556960000
!  - 198648429867720 / 7604556960000
!    319035784479840 / 7604556960000
!  - 419127951114198 / 7604556960000
!    461327344340680 / 7604556960000
!  - 419127951114198 / 7604556960000
!    319035784479840 / 7604556960000
!  - 198648429867720 / 7604556960000
!    103680563465808 / 7604556960000
!   - 42368630685840 / 7604556960000
!     15529830312096 / 7604556960000
!    - 3212744374395 / 7604556960000
!      1848730221900 / 7604556960000
!       203732352169 / 7604556960000
!
    x(1) = -1.00000000000000000000D+00
    x(2) = -0.88888888888888888889D+00
    x(3) = -0.77777777777777777778D+00
    x(4) = -0.66666666666666666667D+00
    x(5) = -0.55555555555555555556D+00
    x(6) = -0.44444444444444444444D+00
    x(7) = -0.33333333333333333333D+00
    x(8) = -0.22222222222222222222D+00
    x(9) = -0.11111111111111111111D+00
    x(10) = 0.00000000000000000000D+00
    x(11) = 0.11111111111111111111D+00
    x(12) = 0.22222222222222222222D+00
    x(13) = 0.33333333333333333333D+00
    x(14) = 0.44444444444444444444D+00
    x(15) = 0.55555555555555555556D+00
    x(16) = 0.66666666666666666667D+00
    x(17) = 0.77777777777777777778D+00
    x(18) = 0.88888888888888888889D+00
    x(19) = 1.00000000000000000000D+00

    w(1) =    0.026790824664820447344D+00
    w(2) =    0.24310820888374278151D+00
    w(3) =   -0.42247620621346493274D+00
    w(4) =    2.0421742376029227612D+00
    w(5) =   -5.5714791681749728126D+00
    w(6) =   13.634004454324976218D+00
    w(7) =  -26.122288374274995239D+00
    w(8) =   41.953237533490708445D+00
    w(9) =  -55.115367445968607749D+00
    w(10) =  60.664591871329740161D+00
    w(11) = -55.115367445968607749D+00
    w(12) =  41.953237533490708445D+00
    w(13) = -26.122288374274995239D+00
    w(14) =  13.634004454324976218D+00
    w(15) =  -5.5714791681749728126D+00
    w(16) =   2.0421742376029227612D+00
    w(17) =  -0.42247620621346493274D+00
    w(18) =   0.24310820888374278151D+00
    w(19) =   0.026790824664820447344D+00

  else if ( n == 20 ) then
!
!       69028763155644023 / 2688996956405760000
!      603652082270808125 / 2688996956405760000
!    - 926840515700222955 / 2688996956405760000
!     4301581538450500095 / 2688996956405760000
!  - 10343692234243192788 / 2688996956405760000
!    22336420328479961316 / 2688996956405760000
!  - 35331888421114781580 / 2688996956405760000
!    43920768370565135580 / 2688996956405760000
!  - 37088370261379851390 / 2688996956405760000
!    15148337305921759574 / 2688996956405760000
!    15148337305921759574 / 2688996956405760000
!  - 37088370261379851390 / 2688996956405760000
!    43920768370565135580 / 2688996956405760000
!  - 35331888421114781580 / 2688996956405760000
!    22336420328479961316 / 2688996956405760000
!  - 10343692234243192788 / 2688996956405760000
!     4301581538450500095 / 2688996956405760000
!    - 926840515700222955 / 2688996956405760000
!      603652082270808125 / 2688996956405760000
!       69028763155644023 / 2688996956405760000
!
    x(1) =  -1.00000000000000000000D+00
    x(2) =  -0.89473684210526315789D+00
    x(3) =  -0.78947368421052631579D+00
    x(4) =  -0.68421052631578947368D+00
    x(5) =  -0.57894736842105263158D+00
    x(6) =  -0.47368421052631578947D+00
    x(7) =  -0.36842105263157894737D+00
    x(8) =  -0.26315789473684210526D+00
    x(9) =  -0.15789473684210526316D+00
    x(10) = -0.052631578947368421053D+00
    x(11) =  0.052631578947368421053D+00
    x(12) =  0.15789473684210526316D+00
    x(13) =  0.26315789473684210526D+00
    x(14) =  0.36842105263157894737D+00
    x(15) =  0.47368421052631578947D+00
    x(16) =  0.57894736842105263158D+00
    x(17) =  0.68421052631578947368D+00
    x(18) =  0.78947368421052631579D+00
    x(19) =  0.89473684210526315789D+00
    x(20) =  1.00000000000000000000D+00

    w(1) =    0.025670822345560078100D+00
    w(2) =    0.22448968595251886556D+00
    w(3) =   -0.34467890099030890987D+00
    w(4) =    1.5996974366978074270D+00
    w(5) =   -3.8466730910952978835D+00
    w(6) =    8.3065993344729824120D+00
    w(7) =  -13.139430424771119113D+00
    w(8) =   16.333513604742678295D+00
    w(9) =  -13.792641220001198577D+00
    w(10) =   5.6334527526463774045D+00
    w(11) =   5.6334527526463774045D+00
    w(12) = -13.792641220001198577D+00
    w(13) =  16.333513604742678295D+00
    w(14) = -13.139430424771119113D+00
    w(15) =   8.3065993344729824120D+00
    w(16) =  -3.8466730910952978835D+00
    w(17) =   1.5996974366978074270D+00
    w(18) =  -0.34467890099030890987D+00
    w(19) =   0.22448968595251886556D+00
    w(20) =   0.025670822345560078100D+00

  else if ( n == 21 ) then

    x(1) =  -1.00000000000000000000D+00
    x(2) =  -0.90000000000000000000D+00
    x(3) =  -0.80000000000000000000D+00
    x(4) =  -0.70000000000000000000D+00
    x(5) =  -0.60000000000000000000D+00
    x(6) =  -0.50000000000000000000D+00
    x(7) =  -0.40000000000000000000D+00
    x(8) =  -0.30000000000000000000D+00
    x(9) =  -0.20000000000000000000D+00
    x(10) = -0.10000000000000000000D+00
    x(11) =  0.00000000000000000000D+00
    x(12) =  0.10000000000000000000D+00
    x(13) =  0.20000000000000000000D+00
    x(14) =  0.30000000000000000000D+00
    x(15) =  0.40000000000000000000D+00
    x(16) =  0.50000000000000000000D+00
    x(17) =  0.60000000000000000000D+00
    x(18) =  0.70000000000000000000D+00
    x(19) =  0.80000000000000000000D+00
    x(20) =  0.90000000000000000000D+00
    x(21) =  1.00000000000000000000D+00

    w(1) =     0.023650546498063206389D+00
    w(2) =     0.22827543528921394997D+00
    w(3) =    -0.47295674102285392846D+00
    w(4) =     2.4123737869637513288D+00
    w(5) =    -7.5420634534306609355D+00
    w(6) =    20.673596439879602287D+00
    w(7) =   -45.417631687959024596D+00
    w(8) =    83.656114844387109207D+00
    w(9) =  -128.15055898030800930D+00
    w(10) =  165.59456694494570344D+00
    w(11) = -180.01073427048578932D+00
    w(12) =  165.59456694494570344D+00
    w(13) = -128.15055898030800930D+00
    w(14) =   83.656114844387109207D+00
    w(15) =  -45.417631687959024596D+00
    w(16) =   20.673596439879602287D+00
    w(17) =   -7.5420634534306609355D+00
    w(18) =    2.4123737869637513288D+00
    w(19) =   -0.47295674102285392846D+00
    w(20) =    0.22827543528921394997D+00
    w(21) =    0.023650546498063206389D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCC_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 through 21.'
    stop

  end if

  return
end
subroutine nco_compute ( n, x, w )

!*****************************************************************************80
!
!! NCO_COMPUTE computes a Newton-Cotes Open quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      Integral ( -1 <= X <= +1 ) F(X) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
!
!    For the OPEN rule, the abscissas do not include the end points.
!
!    For this common type of Newton-Cotes Open rule, the abscissas for
!    rule N are found by computing the abscissas for the closed rule of
!    order N+2 and dropping the two endpoints.
!
!    (An alternative "open half" rule of order N can be defined, which divides
!    the interval into N equal subintervals, and uses the midpoint
!    of each subinterval as the abscissa.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  call nco_compute_points ( n, x )

  call nco_compute_weights ( n, w )

  return
end
subroutine nco_compute_points ( n, x )

!*****************************************************************************80
!
!! NCO_COMPUTE_POINTS: points for a Newton-Cotes Open quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00
  real ( kind = 8 ) x(n)

  do i = 1, n
    x(i) = ( real ( n - i + 1, kind = 8 ) * x_min   &
           + real (     i,     kind = 8 ) * x_max ) &
           / real ( n     + 1, kind = 8 )
  end do

  return
end
subroutine nco_compute_weights ( n, w )

!*****************************************************************************80
!
!! NCO_COMPUTE_WEIGHTS: weights for a Newton-Cotes Open quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  call nco_compute_points ( n, x )

  call nc_compute ( n, x_min, x_max, x, w )

  return
end
subroutine nco_set ( n, x, w )

!*****************************************************************************80
!
!! NCO_SET sets abscissas and weights for open Newton-Cotes quadrature.
!
!  Discussion:
!
!    The open Newton-Cotes rules use equally spaced abscissas, and
!    hence may be used with equally spaced data.
!
!    The rules are called "open" because the abscissas do not include
!    the interval endpoints.
!
!    For this common type of open Newton-Cotes rule, the abscissas for
!    rule N are found by computing the abscissas for the closed rule of
!    order N+2 and dropping the two endpoints.
!
!    (An alternative "open half" rule of order N can be defined, which divides
!    the interval into N equal subintervals, and uses the midpoint
!    of each subinterval as the abscissa.)
!
!    Most of the rules involve negative weights.  These can produce loss
!    of accuracy due to the subtraction of large, nearly equal quantities.
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 7, and 9.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    w(1) = 2.0D+00

  else if ( n == 2 ) then

    w(1) = 1.0D+00
    w(2) = 1.0D+00

  else if ( n == 3 ) then

    d = 3.0D+00

    w(1) =   4.0D+00 / d
    w(2) = - 2.0D+00 / d
    w(3) =   4.0D+00 / d

  else if ( n == 4 ) then

    d = 12.0D+00

    w(1) = 11.0D+00 / d
    w(2) =  1.0D+00 / d
    w(3) =  1.0D+00 / d
    w(4) = 11.0D+00 / d

  else if ( n == 5 ) then

    d = 10.0D+00

    w(1) =   11.0D+00 / d
    w(2) = - 14.0D+00 / d
    w(3) =   26.0D+00 / d
    w(4) = - 14.0D+00 / d
    w(5) =   11.0D+00 / d

  else if ( n == 6 ) then

    d = 1440.0D+00

    w(1) =  1222.0D+00 / d
    w(2) = - 906.0D+00 / d
    w(3) =  1124.0D+00 / d
    w(4) =  1124.0D+00 / d
    w(5) = - 906.0D+00 / d
    w(6) =  1222.0D+00 / d

  else if ( n == 7 ) then

    d = 945.0D+00

    w(1) =    920.0D+00 / d
    w(2) = - 1908.0D+00 / d
    w(3) =   4392.0D+00 / d
    w(4) = - 4918.0D+00 / d
    w(5) =   4392.0D+00 / d
    w(6) = - 1908.0D+00 / d
    w(7) =    920.0D+00 / d

  else if ( n == 9 ) then

    d = 4536.0D+00

    w(1) =    4045.0D+00 / d
    w(2) = - 11690.0D+00 / d
    w(3) =   33340.0D+00 / d
    w(4) = - 55070.0D+00 / d
    w(5) =   67822.0D+00 / d
    w(6) = - 55070.0D+00 / d
    w(7) =   33340.0D+00 / d
    w(8) = - 11690.0D+00 / d
    w(9) =    4045.0D+00 / d

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCO_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 to 7, and 9.'
    stop

  end if
!
!  Set the abscissas.
!
  a = -1.0D+00
  b =  1.0D+00

  do i = 1, n
    x(i) = ( real ( n + 1 - i, kind = 8 ) * a   &
           + real (         i, kind = 8 ) * b ) &
           / real ( n + 1,     kind = 8 )
  end do

  return
end
subroutine ncoh_compute ( n, x, w )

!*****************************************************************************80
!
!! NCOH_COMPUTE computes a Newton-Cotes Open Half quadrature rule.
!
!  Discussion:
!
!    The input value N is used to define N equal subintervals of [-1,+1].
!    The I-th abscissa is the center of the I-th subinterval.
!
!    The integral:
!
!      Integral ( X_MIN <= X <= X_MAX ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  do i = 1, n
    x(i) = ( real ( 2 * n - 2 * i + 1, kind = 8 ) * x_min   &
           + real (         2 * i - 1, kind = 8 ) * x_max ) &
           / real ( 2 * n,             kind = 8 )
  end do

  call nc_compute ( n, x_min, x_max, x, w )

  return
end
subroutine ncoh_compute_points ( n, x )

!*****************************************************************************80
!
!! NCOH_COMPUTE_POINTS: points for a Newton-Cotes Open Half quadrature rule.
!
!  Discussion:
!
!    The input value N is used to define N equal subintervals of [-1,+1].
!    The I-th abscissa is the center of the I-th subinterval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  do i = 1, n
    x(i) = ( real ( 2 * n - 2 * i + 1, kind = 8 ) * x_min   &
           + real (         2 * i - 1, kind = 8 ) * x_max ) &
           / real ( 2 * n,             kind = 8 )
  end do

  return
end
subroutine ncoh_compute_weights ( n, w )

!*****************************************************************************80
!
!! NCOH_COMPUTE_WEIGHTS: weights for a Newton-Cotes Open Half quadrature rule.
!
!  Discussion:
!
!    The input value N is used to define N equal subintervals of [-1,+1].
!    The I-th abscissa is the center of the I-th subinterval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  call ncoh_compute_points ( n, x )

  call nc_compute ( n, x_min, x_max, x, w )

  return
end
subroutine ncoh_set ( n, x, w )

!*****************************************************************************80
!
!! NCOH_SET sets abscissas and weights for Newton-Cotes "open half" quadrature.
!
!  Discussion:
!
!    The Newton-Cotes rules use equally spaced abscissas, and
!    hence may be used with equally spaced data.
!
!    The rules are called "open" because the abscissas do not include
!    the interval endpoints.
!
!    For this uncommon type of Newton-Cotes open rule, the abscissas for
!    rule N are found by dividing the interval into N equal subintervals,
!    and using the midpoint of each subinterval as the abscissa.
!
!    Most of the rules involve negative weights.  These can produce loss
!    of accuracy due to the subtraction of large, nearly equal quantities.
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    In Mathematica, the Newton-Cotes "open half" weights and abscissas
!    can be computed by the commands:
!
!      Needs["NumericalDifferentialEquationAnalysis`"]
!      NewtonCotesWeights [ n, -1, 1, QuadratureType -> Open ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N <= 17.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    w(1) = 2.0D+00

  else if ( n == 2 ) then

    w(1) = 1.0D+00
    w(2) = 1.0D+00

  else if ( n == 3 ) then

    d = 4.0D+00

    w(1) =   3.0D+00 / d
    w(2) =   2.0D+00 / d
    w(3) =   3.0D+00 / d

  else if ( n == 4 ) then

    d = 24.0D+00

    w(1) = 13.0D+00 / d
    w(2) = 11.0D+00 / d
    w(3) = 11.0D+00 / d
    w(4) = 13.0D+00 / d

  else if ( n == 5 ) then

    d = 576.0D+00

    w(1) =  275.0D+00 / d
    w(2) =  100.0D+00 / d
    w(3) =  402.0D+00 / d
    w(4) =  100.0D+00 / d
    w(5) =  275.0D+00 / d

  else if ( n == 6 ) then

    d = 640.0D+00

    w(1) =   247.0D+00 / d
    w(2) =   139.0D+00 / d
    w(3) =   254.0D+00 / d
    w(4) =   254.0D+00 / d
    w(5) =   139.0D+00 / d
    w(6) =   247.0D+00 / d

  else if ( n == 7 ) then

    d = 138240.0D+00

    w(1) =   49490.0D+00 / d
    w(2) =    1764.0D+00 / d
    w(3) =  112014.0D+00 / d
    w(4) =  -50056.0D+00 / d
    w(5) =  112014.0D+00 / d
    w(6) =    1764.0D+00 / d
    w(7) =   49490.0D+00 / d

  else if ( n == 8 ) then

    d = 967680.0D+00

    w(1) =  295627.0D+00 / d
    w(2) =   71329.0D+00 / d
    w(3) =  471771.0D+00 / d
    w(4) =  128953.0D+00 / d
    w(5) =  128953.0D+00 / d
    w(6) =  471771.0D+00 / d
    w(7) =   71329.0D+00 / d
    w(8) =  295627.0D+00 / d

  else if ( n == 9 ) then

    d = 2867200.0D+00

    w(1) =    832221.0D+00 / d
    w(2) =   -260808.0D+00 / d
    w(3) =   2903148.0D+00 / d
    w(4) =  -3227256.0D+00 / d
    w(5) =   5239790.0D+00 / d
    w(6) =  -3227256.0D+00 / d
    w(7) =   2903148.0D+00 / d
    w(8) =   -260808.0D+00 / d
    w(9) =    832221.0D+00 / d

  else if ( n == 10 ) then

    w( 1) = + 0.255727885681905864197530864198D+00
    w( 2) = - 0.0265214977230764991181657848325D+00
    w( 3) = + 0.660404481164572310405643738977D+00
    w( 4) = - 0.337696647307649911816578483245D+00
    w( 5) = + 0.448085778184248236331569664903D+00
    w( 6) = + 0.448085778184248236331569664903D+00
    w( 7) = - 0.337696647307649911816578483245D+00
    w( 8) = + 0.660404481164572310405643738977D+00
    w( 9) = - 0.0265214977230764991181657848325D+00
    w(10) = + 0.255727885681905864197530864198D+00

  else if ( n == 11 ) then

    w( 1) = + 0.246271364278193434009406231628D+00
    w( 2) = - 0.167027133984260177836566725456D+00
    w( 3) = + 1.27129728179339588844797178131D+00
    w( 4) = - 2.19004533609595458553791887125D+00
    w( 5) = + 3.91748917836991567460317460317D+00
    w( 6) = - 4.15597070872258046737213403880D+00
    w( 7) = + 3.91748917836991567460317460317D+00
    w( 8) = - 2.19004533609595458553791887125D+00
    w( 9) = + 1.27129728179339588844797178131D+00
    w(10) = - 0.167027133984260177836566725456D+00
    w(11) = + 0.246271364278193434009406231628D+00

  else if ( n == 12 ) then

    w( 1) = + 0.221603210581329477813852813853D+00
    w( 2) = - 0.103156166902352205086580086580D+00
    w( 3) = + 0.889254983348763866341991341991D+00
    w( 4) = - 1.08160728355506797889610389610D+00
    w( 5) = + 1.49180546087620062229437229437D+00
    w( 6) = - 0.417900204348873782467532467532D+00
    w( 7) = - 0.417900204348873782467532467532D+00
    w( 8) = + 1.49180546087620062229437229437D+00
    w( 9) = - 1.08160728355506797889610389610D+00
    w(10) = + 0.889254983348763866341991341991D+00
    w(11) = - 0.103156166902352205086580086580D+00
    w(12) = + 0.221603210581329477813852813853D+00

  else if ( n == 13 ) then

    w(1) = 0.215232356419153566228270676022D+00
    w(2) = -0.227154289276070155983970468097D+00
    w(3) = 1.57154640756958579127322926926D+00
    w(4) = -3.60188931740556785445074962271D+00
    w(5) = 7.51615534838963020202557032914D+00
    w(6) = -10.7785343238762239297023523214D+00
    w(7) = 12.6092876363589847612200042756D+00
    w(8) = -10.7785343238762239297023523214D+00
    w(9) = 7.51615534838963020202557032914D+00
    w(10) = -3.60188931740556785445074962271D+00
    w(11) = 1.57154640756958579127322926926D+00
    w(12) = -0.227154289276070155983970468097D+00
    w(13) = 0.215232356419153566228270676022D+00

  else if ( n == 14 ) then

    w(1) = 0.196600731862944474955289480752D+00
    w(2) = -0.165179242362168469504443173425D+00
    w(3) = 1.16085790162743923526801130968D+00
    w(4) = -2.14582271238684154514413484321D+00
    w(5) = 3.66584923423684682693019643251D+00
    w(6) = -3.34045051168652382743365816282D+00
    w(7) = 1.62814459870830330492873895652D+00
    w(8) = 1.62814459870830330492873895652D+00
    w(9) = -3.34045051168652382743365816282D+00
    w(10) = 3.66584923423684682693019643251D+00
    w(11) = -2.14582271238684154514413484321D+00
    w(12) = 1.16085790162743923526801130968D+00
    w(13) = -0.165179242362168469504443173425D+00
    w(14) = 0.196600731862944474955289480752D+00

  else if ( n == 15 ) then

    w(1) = 0.192053656112251156523443074782D+00
    w(2) = -0.277042941258250537556131864168D+00
    w(3) = 1.90509434600895135399617123947D+00
    w(4) = -5.39701622989083452471078029114D+00
    w(5) = 13.1085281753546466152623727959D+00
    w(6) = -23.3466945206436771323681898459D+00
    w(7) = 33.4478422682091702199516443378D+00
    w(8) = -37.2655295077845143021970588935D+00
    w(9) = 33.4478422682091702199516443378D+00
    w(10) = -23.3466945206436771323681898459D+00
    w(11) = 13.1085281753546466152623727959D+00
    w(12) = -5.39701622989083452471078029114D+00
    w(13) = 1.90509434600895135399617123947D+00
    w(14) = -0.277042941258250537556131864168D+00
    w(15) = 0.192053656112251156523443074782D+00

  else if ( n == 16 ) then

    w(1) = 0.177408479879589716830780293564D+00
    w(2) = -0.217359399771056183974705378131D+00
    w(3) = 1.46740967914595726066296033468D+00
    w(4) = -3.56820982596198712548876407280D+00
    w(5) = 7.42429624597447227175662173974D+00
    w(6) = -10.1614344802943189309930887295D+00
    w(7) = 9.74825566269696996625640284529D+00
    w(8) = -3.87036636166962697505020703289D+00
    w(9) = -3.87036636166962697505020703289D+00
    w(10) = 9.74825566269696996625640284529D+00
    w(11) = -10.1614344802943189309930887295D+00
    w(12) = 7.42429624597447227175662173974D+00
    w(13) = -3.56820982596198712548876407280D+00
    w(14) = 1.46740967914595726066296033468D+00
    w(15) = -0.217359399771056183974705378131D+00
    w(16) = 0.177408479879589716830780293564D+00

  else if ( n == 17 ) then

    w(1) = 0.174021728363203743659784786159D+00
    w(2) = -0.319844636797863622878597303396D+00
    w(3) = 2.26685253417620917889510925819D+00
    w(4) = -7.60565246092744264614795072379D+00
    w(5) = 21.2205863313196208783745036601D+00
    w(6) = -44.9326914054546061828308816595D+00
    w(7) = 76.6598740687724224896458863733D+00
    w(8) = -104.621704713086021393464459433D+00
    w(9) = 116.317117107268955109493210084D+00
    w(10) = -104.621704713086021393464459433D+00
    w(11) = 76.6598740687724224896458863733D+00
    w(12) = -44.9326914054546061828308816595D+00
    w(13) = 21.2205863313196208783745036601D+00
    w(14) = -7.60565246092744264614795072379D+00
    w(15) = 2.26685253417620917889510925819D+00
    w(16) = -0.319844636797863622878597303396D+00
    w(17) = 0.174021728363203743659784786159D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCOH_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 to 17.'
    stop

  end if
!
!  Set the abscissas.
!
  a = - 1.0D+00
  b = + 1.0D+00

  do i = 1, n
    x(i) = ( real ( 2 * n - 2 * i + 1, kind = 8 ) * a   &
           + real (         2 * i - 1, kind = 8 ) * b ) &
           / real ( 2 * n,             kind = 8 )
  end do

  return
end
subroutine patterson_set ( n, x, w )

!*****************************************************************************80
!
!! PATTERSON_SET sets abscissas and weights for Gauss-Patterson quadrature.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The zeroth rule, of order 1, is the standard Gauss-Legendre rule.
!
!    The first rule, of order 3, is the standard Gauss-Legendre rule.
!
!    The second rule, of order 7, includes the abscissas of the previous
!    rule.
!
!    Rules are available of orders 1, 3, 7, 15, 31, 63, 127, 255, and 511.
!
!    These rules constitute a nested family.  The rules can integrate exactly
!    any polynomial of degree 1, 5, 11, 23, 47, 95, 191, 383 or 767, 
!    respectively.
!
!    The data for N = 511 was supplied by Dirk Laurie, and is derived
!    from a NAG Library function d01arf.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!    NAG Library Documentation,
!    D01ARF,
!    The Numerical Algorithms Group.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be 1, 3, 7, 15, 31, 63, 127, 255 or 511.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =   0.0D+00

    w(1) = 2.0D+00

  else if ( n == 3 ) then

    x(1) = -0.77459666924148337704D+00
    x(2) =  0.0D+00
    x(3) =  0.77459666924148337704D+00

    w(1) = 0.555555555555555555556D+00
    w(2) = 0.888888888888888888889D+00
    w(3) = 0.555555555555555555556D+00

  else if ( n == 7 ) then

    x(1) = -0.96049126870802028342D+00
    x(2) = -0.77459666924148337704D+00
    x(3) = -0.43424374934680255800D+00
    x(4) =  0.0D+00
    x(5) =  0.43424374934680255800D+00
    x(6) =  0.77459666924148337704D+00
    x(7) =  0.96049126870802028342D+00

    w(1) = 0.104656226026467265194D+00
    w(2) = 0.268488089868333440729D+00
    w(3) = 0.401397414775962222905D+00
    w(4) = 0.450916538658474142345D+00
    w(5) = 0.401397414775962222905D+00
    w(6) = 0.268488089868333440729D+00
    w(7) = 0.104656226026467265194D+00

  else if ( n == 15 ) then

    x( 1) = -0.99383196321275502221D+00
    x( 2) = -0.96049126870802028342D+00
    x( 3) = -0.88845923287225699889D+00
    x( 4) = -0.77459666924148337704D+00
    x( 5) = -0.62110294673722640294D+00
    x( 6) = -0.43424374934680255800D+00
    x( 7) = -0.22338668642896688163D+00
    x( 8) =  0.0D+00
    x( 9) =  0.22338668642896688163D+00
    x(10) =  0.43424374934680255800D+00
    x(11) =  0.62110294673722640294D+00
    x(12) =  0.77459666924148337704D+00
    x(13) =  0.88845923287225699889D+00
    x(14) =  0.96049126870802028342D+00
    x(15) =  0.99383196321275502221D+00

    w( 1) = 0.0170017196299402603390D+00
    w( 2) = 0.0516032829970797396969D+00
    w( 3) = 0.0929271953151245376859D+00
    w( 4) = 0.134415255243784220360D+00
    w( 5) = 0.171511909136391380787D+00
    w( 6) = 0.200628529376989021034D+00
    w( 7) = 0.219156858401587496404D+00
    w( 8) = 0.225510499798206687386D+00
    w( 9) = 0.219156858401587496404D+00
    w(10) = 0.200628529376989021034D+00
    w(11) = 0.171511909136391380787D+00
    w(12) = 0.134415255243784220360D+00
    w(13) = 0.0929271953151245376859D+00
    w(14) = 0.0516032829970797396969D+00
    w(15) = 0.0170017196299402603390D+00

  else if ( n == 31 ) then

    x( 1) = -0.99909812496766759766D+00
    x( 2) = -0.99383196321275502221D+00
    x( 3) = -0.98153114955374010687D+00
    x( 4) = -0.96049126870802028342D+00
    x( 5) = -0.92965485742974005667D+00
    x( 6) = -0.88845923287225699889D+00
    x( 7) = -0.83672593816886873550D+00
    x( 8) = -0.77459666924148337704D+00
    x( 9) = -0.70249620649152707861D+00
    x(10) = -0.62110294673722640294D+00
    x(11) = -0.53131974364437562397D+00
    x(12) = -0.43424374934680255800D+00
    x(13) = -0.33113539325797683309D+00
    x(14) = -0.22338668642896688163D+00
    x(15) = -0.11248894313318662575D+00
    x(16) =  0.0D+00
    x(17) =  0.11248894313318662575D+00
    x(18) =  0.22338668642896688163D+00
    x(19) =  0.33113539325797683309D+00
    x(20) =  0.43424374934680255800D+00
    x(21) =  0.53131974364437562397D+00
    x(22) =  0.62110294673722640294D+00
    x(23) =  0.70249620649152707861D+00
    x(24) =  0.77459666924148337704D+00
    x(25) =  0.83672593816886873550D+00
    x(26) =  0.88845923287225699889D+00
    x(27) =  0.92965485742974005667D+00
    x(28) =  0.96049126870802028342D+00
    x(29) =  0.98153114955374010687D+00
    x(30) =  0.99383196321275502221D+00
    x(31) =  0.99909812496766759766D+00

    w( 1) = 0.00254478079156187441540D+00
    w( 2) = 0.00843456573932110624631D+00
    w( 3) = 0.0164460498543878109338D+00
    w( 4) = 0.0258075980961766535646D+00
    w( 5) = 0.0359571033071293220968D+00
    w( 6) = 0.0464628932617579865414D+00
    w( 7) = 0.0569795094941233574122D+00
    w( 8) = 0.0672077542959907035404D+00
    w( 9) = 0.0768796204990035310427D+00
    w(10) = 0.0857559200499903511542D+00
    w(11) = 0.0936271099812644736167D+00
    w(12) = 0.100314278611795578771D+00
    w(13) = 0.105669893580234809744D+00
    w(14) = 0.109578421055924638237D+00
    w(15) = 0.111956873020953456880D+00
    w(16) = 0.112755256720768691607D+00
    w(17) = 0.111956873020953456880D+00
    w(18) = 0.109578421055924638237D+00
    w(19) = 0.105669893580234809744D+00
    w(20) = 0.100314278611795578771D+00
    w(21) = 0.0936271099812644736167D+00
    w(22) = 0.0857559200499903511542D+00
    w(23) = 0.0768796204990035310427D+00
    w(24) = 0.0672077542959907035404D+00
    w(25) = 0.0569795094941233574122D+00
    w(26) = 0.0464628932617579865414D+00
    w(27) = 0.0359571033071293220968D+00
    w(28) = 0.0258075980961766535646D+00
    w(29) = 0.0164460498543878109338D+00
    w(30) = 0.00843456573932110624631D+00
    w(31) = 0.00254478079156187441540D+00

  else if ( n == 63 ) then

    x( 1) = -0.99987288812035761194D+00
    x( 2) = -0.99909812496766759766D+00
    x( 3) = -0.99720625937222195908D+00
    x( 4) = -0.99383196321275502221D+00
    x( 5) = -0.98868475754742947994D+00
    x( 6) = -0.98153114955374010687D+00
    x( 7) = -0.97218287474858179658D+00
    x( 8) = -0.96049126870802028342D+00
    x( 9) = -0.94634285837340290515D+00
    x(10) = -0.92965485742974005667D+00
    x(11) = -0.91037115695700429250D+00
    x(12) = -0.88845923287225699889D+00
    x(13) = -0.86390793819369047715D+00
    x(14) = -0.83672593816886873550D+00
    x(15) = -0.80694053195021761186D+00
    x(16) = -0.77459666924148337704D+00
    x(17) = -0.73975604435269475868D+00
    x(18) = -0.70249620649152707861D+00
    x(19) = -0.66290966002478059546D+00
    x(20) = -0.62110294673722640294D+00
    x(21) = -0.57719571005204581484D+00
    x(22) = -0.53131974364437562397D+00
    x(23) = -0.48361802694584102756D+00
    x(24) = -0.43424374934680255800D+00
    x(25) = -0.38335932419873034692D+00
    x(26) = -0.33113539325797683309D+00
    x(27) = -0.27774982202182431507D+00
    x(28) = -0.22338668642896688163D+00
    x(29) = -0.16823525155220746498D+00
    x(30) = -0.11248894313318662575D+00
    x(31) = -0.056344313046592789972D+00
    x(32) =  0.0D+00
    x(33) =  0.056344313046592789972D+00
    x(34) =  0.11248894313318662575D+00
    x(35) =  0.16823525155220746498D+00
    x(36) =  0.22338668642896688163D+00
    x(37) =  0.27774982202182431507D+00
    x(38) =  0.33113539325797683309D+00
    x(39) =  0.38335932419873034692D+00
    x(40) =  0.43424374934680255800D+00
    x(41) =  0.48361802694584102756D+00
    x(42) =  0.53131974364437562397D+00
    x(43) =  0.57719571005204581484D+00
    x(44) =  0.62110294673722640294D+00
    x(45) =  0.66290966002478059546D+00
    x(46) =  0.70249620649152707861D+00
    x(47) =  0.73975604435269475868D+00
    x(48) =  0.77459666924148337704D+00
    x(49) =  0.80694053195021761186D+00
    x(50) =  0.83672593816886873550D+00
    x(51) =  0.86390793819369047715D+00
    x(52) =  0.88845923287225699889D+00
    x(53) =  0.91037115695700429250D+00
    x(54) =  0.92965485742974005667D+00
    x(55) =  0.94634285837340290515D+00
    x(56) =  0.96049126870802028342D+00
    x(57) =  0.97218287474858179658D+00
    x(58) =  0.98153114955374010687D+00
    x(59) =  0.98868475754742947994D+00
    x(60) =  0.99383196321275502221D+00
    x(61) =  0.99720625937222195908D+00
    x(62) =  0.99909812496766759766D+00
    x(63) =  0.99987288812035761194D+00

    w( 1) = 0.000363221481845530659694D+00
    w( 2) = 0.00126515655623006801137D+00
    w( 3) = 0.00257904979468568827243D+00
    w( 4) = 0.00421763044155885483908D+00
    w( 5) = 0.00611550682211724633968D+00
    w( 6) = 0.00822300795723592966926D+00
    w( 7) = 0.0104982469096213218983D+00
    w( 8) = 0.0129038001003512656260D+00
    w( 9) = 0.0154067504665594978021D+00
    w(10) = 0.0179785515681282703329D+00
    w(11) = 0.0205942339159127111492D+00
    w(12) = 0.0232314466399102694433D+00
    w(13) = 0.0258696793272147469108D+00
    w(14) = 0.0284897547458335486125D+00
    w(15) = 0.0310735511116879648799D+00
    w(16) = 0.0336038771482077305417D+00
    w(17) = 0.0360644327807825726401D+00
    w(18) = 0.0384398102494555320386D+00
    w(19) = 0.0407155101169443189339D+00
    w(20) = 0.0428779600250077344929D+00
    w(21) = 0.0449145316536321974143D+00
    w(22) = 0.0468135549906280124026D+00
    w(23) = 0.0485643304066731987159D+00
    w(24) = 0.0501571393058995374137D+00
    w(25) = 0.0515832539520484587768D+00
    w(26) = 0.0528349467901165198621D+00
    w(27) = 0.0539054993352660639269D+00
    w(28) = 0.0547892105279628650322D+00
    w(29) = 0.0554814043565593639878D+00
    w(30) = 0.0559784365104763194076D+00
    w(31) = 0.0562776998312543012726D+00
    w(32) = 0.0563776283603847173877D+00
    w(33) = 0.0562776998312543012726D+00
    w(34) = 0.0559784365104763194076D+00
    w(35) = 0.0554814043565593639878D+00
    w(36) = 0.0547892105279628650322D+00
    w(37) = 0.0539054993352660639269D+00
    w(38) = 0.0528349467901165198621D+00
    w(39) = 0.0515832539520484587768D+00
    w(40) = 0.0501571393058995374137D+00
    w(41) = 0.0485643304066731987159D+00
    w(42) = 0.0468135549906280124026D+00
    w(43) = 0.0449145316536321974143D+00
    w(44) = 0.0428779600250077344929D+00
    w(45) = 0.0407155101169443189339D+00
    w(46) = 0.0384398102494555320386D+00
    w(47) = 0.0360644327807825726401D+00
    w(48) = 0.0336038771482077305417D+00
    w(49) = 0.0310735511116879648799D+00
    w(50) = 0.0284897547458335486125D+00
    w(51) = 0.0258696793272147469108D+00
    w(52) = 0.0232314466399102694433D+00
    w(53) = 0.0205942339159127111492D+00
    w(54) = 0.0179785515681282703329D+00
    w(55) = 0.0154067504665594978021D+00
    w(56) = 0.0129038001003512656260D+00
    w(57) = 0.0104982469096213218983D+00
    w(58) = 0.00822300795723592966926D+00
    w(59) = 0.00611550682211724633968D+00
    w(60) = 0.00421763044155885483908D+00
    w(61) = 0.00257904979468568827243D+00
    w(62) = 0.00126515655623006801137D+00
    w(63) = 0.000363221481845530659694D+00

  else if ( n == 127 ) then

    x(  1) = -0.99998243035489159858D+00
    x(  2) = -0.99987288812035761194D+00
    x(  3) = -0.99959879967191068325D+00
    x(  4) = -0.99909812496766759766D+00
    x(  5) = -0.99831663531840739253D+00
    x(  6) = -0.99720625937222195908D+00
    x(  7) = -0.99572410469840718851D+00
    x(  8) = -0.99383196321275502221D+00
    x(  9) = -0.99149572117810613240D+00
    x( 10) = -0.98868475754742947994D+00
    x( 11) = -0.98537149959852037111D+00
    x( 12) = -0.98153114955374010687D+00
    x( 13) = -0.97714151463970571416D+00
    x( 14) = -0.97218287474858179658D+00
    x( 15) = -0.96663785155841656709D+00
    x( 16) = -0.96049126870802028342D+00
    x( 17) = -0.95373000642576113641D+00
    x( 18) = -0.94634285837340290515D+00
    x( 19) = -0.93832039777959288365D+00
    x( 20) = -0.92965485742974005667D+00
    x( 21) = -0.92034002547001242073D+00
    x( 22) = -0.91037115695700429250D+00
    x( 23) = -0.89974489977694003664D+00
    x( 24) = -0.88845923287225699889D+00
    x( 25) = -0.87651341448470526974D+00
    x( 26) = -0.86390793819369047715D+00
    x( 27) = -0.85064449476835027976D+00
    x( 28) = -0.83672593816886873550D+00
    x( 29) = -0.82215625436498040737D+00
    x( 30) = -0.80694053195021761186D+00
    x( 31) = -0.79108493379984836143D+00
    x( 32) = -0.77459666924148337704D+00
    x( 33) = -0.75748396638051363793D+00
    x( 34) = -0.73975604435269475868D+00
    x( 35) = -0.72142308537009891548D+00
    x( 36) = -0.70249620649152707861D+00
    x( 37) = -0.68298743109107922809D+00
    x( 38) = -0.66290966002478059546D+00
    x( 39) = -0.64227664250975951377D+00
    x( 40) = -0.62110294673722640294D+00
    x( 41) = -0.59940393024224289297D+00
    x( 42) = -0.57719571005204581484D+00
    x( 43) = -0.55449513263193254887D+00
    x( 44) = -0.53131974364437562397D+00
    x( 45) = -0.50768775753371660215D+00
    x( 46) = -0.48361802694584102756D+00
    x( 47) = -0.45913001198983233287D+00
    x( 48) = -0.43424374934680255800D+00
    x( 49) = -0.40897982122988867241D+00
    x( 50) = -0.38335932419873034692D+00
    x( 51) = -0.35740383783153215238D+00
    x( 52) = -0.33113539325797683309D+00
    x( 53) = -0.30457644155671404334D+00
    x( 54) = -0.27774982202182431507D+00
    x( 55) = -0.25067873030348317661D+00
    x( 56) = -0.22338668642896688163D+00
    x( 57) = -0.19589750271110015392D+00
    x( 58) = -0.16823525155220746498D+00
    x( 59) = -0.14042423315256017459D+00
    x( 60) = -0.11248894313318662575D+00
    x( 61) = -0.084454040083710883710D+00
    x( 62) = -0.056344313046592789972D+00
    x( 63) = -0.028184648949745694339D+00
    x( 64) =  0.0D+00
    x( 65) =  0.028184648949745694339D+00
    x( 66) =  0.056344313046592789972D+00
    x( 67) =  0.084454040083710883710D+00
    x( 68) =  0.11248894313318662575D+00
    x( 69) =  0.14042423315256017459D+00
    x( 70) =  0.16823525155220746498D+00
    x( 71) =  0.19589750271110015392D+00
    x( 72) =  0.22338668642896688163D+00
    x( 73) =  0.25067873030348317661D+00
    x( 74) =  0.27774982202182431507D+00
    x( 75) =  0.30457644155671404334D+00
    x( 76) =  0.33113539325797683309D+00
    x( 77) =  0.35740383783153215238D+00
    x( 78) =  0.38335932419873034692D+00
    x( 79) =  0.40897982122988867241D+00
    x( 80) =  0.43424374934680255800D+00
    x( 81) =  0.45913001198983233287D+00
    x( 82) =  0.48361802694584102756D+00
    x( 83) =  0.50768775753371660215D+00
    x( 84) =  0.53131974364437562397D+00
    x( 85) =  0.55449513263193254887D+00
    x( 86) =  0.57719571005204581484D+00
    x( 87) =  0.59940393024224289297D+00
    x( 88) =  0.62110294673722640294D+00
    x( 89) =  0.64227664250975951377D+00
    x( 90) =  0.66290966002478059546D+00
    x( 91) =  0.68298743109107922809D+00
    x( 92) =  0.70249620649152707861D+00
    x( 93) =  0.72142308537009891548D+00
    x( 94) =  0.73975604435269475868D+00
    x( 95) =  0.75748396638051363793D+00
    x( 96) =  0.77459666924148337704D+00
    x( 97) =  0.79108493379984836143D+00
    x( 98) =  0.80694053195021761186D+00
    x( 99) =  0.82215625436498040737D+00
    x(100) =  0.83672593816886873550D+00
    x(101) =  0.85064449476835027976D+00
    x(102) =  0.86390793819369047715D+00
    x(103) =  0.87651341448470526974D+00
    x(104) =  0.88845923287225699889D+00
    x(105) =  0.89974489977694003664D+00
    x(106) =  0.91037115695700429250D+00
    x(107) =  0.92034002547001242073D+00
    x(108) =  0.92965485742974005667D+00
    x(109) =  0.93832039777959288365D+00
    x(110) =  0.94634285837340290515D+00
    x(111) =  0.95373000642576113641D+00
    x(112) =  0.96049126870802028342D+00
    x(113) =  0.96663785155841656709D+00
    x(114) =  0.97218287474858179658D+00
    x(115) =  0.97714151463970571416D+00
    x(116) =  0.98153114955374010687D+00
    x(117) =  0.98537149959852037111D+00
    x(118) =  0.98868475754742947994D+00
    x(119) =  0.99149572117810613240D+00
    x(120) =  0.99383196321275502221D+00
    x(121) =  0.99572410469840718851D+00
    x(122) =  0.99720625937222195908D+00
    x(123) =  0.99831663531840739253D+00
    x(124) =  0.99909812496766759766D+00
    x(125) =  0.99959879967191068325D+00
    x(126) =  0.99987288812035761194D+00
    x(127) =  0.99998243035489159858D+00

    w(  1) = 0.0000505360952078625176247D+00
    w(  2) = 0.000180739564445388357820D+00
    w(  3) = 0.000377746646326984660274D+00
    w(  4) = 0.000632607319362633544219D+00
    w(  5) = 0.000938369848542381500794D+00
    w(  6) = 0.00128952408261041739210D+00
    w(  7) = 0.00168114286542146990631D+00
    w(  8) = 0.00210881524572663287933D+00
    w(  9) = 0.00256876494379402037313D+00
    w( 10) = 0.00305775341017553113613D+00
    w( 11) = 0.00357289278351729964938D+00
    w( 12) = 0.00411150397865469304717D+00
    w( 13) = 0.00467105037211432174741D+00
    w( 14) = 0.00524912345480885912513D+00
    w( 15) = 0.00584344987583563950756D+00
    w( 16) = 0.00645190005017573692280D+00
    w( 17) = 0.00707248999543355546805D+00
    w( 18) = 0.00770337523327974184817D+00
    w( 19) = 0.00834283875396815770558D+00
    w( 20) = 0.00898927578406413572328D+00
    w( 21) = 0.00964117772970253669530D+00
    w( 22) = 0.0102971169579563555237D+00
    w( 23) = 0.0109557333878379016480D+00
    w( 24) = 0.0116157233199551347270D+00
    w( 25) = 0.0122758305600827700870D+00
    w( 26) = 0.0129348396636073734547D+00
    w( 27) = 0.0135915710097655467896D+00
    w( 28) = 0.0142448773729167743063D+00
    w( 29) = 0.0148936416648151820348D+00
    w( 30) = 0.0155367755558439824399D+00
    w( 31) = 0.0161732187295777199419D+00
    w( 32) = 0.0168019385741038652709D+00
    w( 33) = 0.0174219301594641737472D+00
    w( 34) = 0.0180322163903912863201D+00
    w( 35) = 0.0186318482561387901863D+00
    w( 36) = 0.0192199051247277660193D+00
    w( 37) = 0.0197954950480974994880D+00
    w( 38) = 0.0203577550584721594669D+00
    w( 39) = 0.0209058514458120238522D+00
    w( 40) = 0.0214389800125038672465D+00
    w( 41) = 0.0219563663053178249393D+00
    w( 42) = 0.0224572658268160987071D+00
    w( 43) = 0.0229409642293877487608D+00
    w( 44) = 0.0234067774953140062013D+00
    w( 45) = 0.0238540521060385400804D+00
    w( 46) = 0.0242821652033365993580D+00
    w( 47) = 0.0246905247444876769091D+00
    w( 48) = 0.0250785696529497687068D+00
    w( 49) = 0.0254457699654647658126D+00
    w( 50) = 0.0257916269760242293884D+00
    w( 51) = 0.0261156733767060976805D+00
    w( 52) = 0.0264174733950582599310D+00
    w( 53) = 0.0266966229274503599062D+00
    w( 54) = 0.0269527496676330319634D+00
    w( 55) = 0.0271855132296247918192D+00
    w( 56) = 0.0273946052639814325161D+00
    w( 57) = 0.0275797495664818730349D+00
    w( 58) = 0.0277407021782796819939D+00
    w( 59) = 0.0278772514766137016085D+00
    w( 60) = 0.0279892182552381597038D+00
    w( 61) = 0.0280764557938172466068D+00
    w( 62) = 0.0281388499156271506363D+00
    w( 63) = 0.0281763190330166021307D+00
    w( 64) = 0.0281888141801923586938D+00
    w( 65) = 0.0281763190330166021307D+00
    w( 66) = 0.0281388499156271506363D+00
    w( 67) = 0.0280764557938172466068D+00
    w( 68) = 0.0279892182552381597038D+00
    w( 69) = 0.0278772514766137016085D+00
    w( 70) = 0.0277407021782796819939D+00
    w( 71) = 0.0275797495664818730349D+00
    w( 72) = 0.0273946052639814325161D+00
    w( 73) = 0.0271855132296247918192D+00
    w( 74) = 0.0269527496676330319634D+00
    w( 75) = 0.0266966229274503599062D+00
    w( 76) = 0.0264174733950582599310D+00
    w( 77) = 0.0261156733767060976805D+00
    w( 78) = 0.0257916269760242293884D+00
    w( 79) = 0.0254457699654647658126D+00
    w( 80) = 0.0250785696529497687068D+00
    w( 81) = 0.0246905247444876769091D+00
    w( 82) = 0.0242821652033365993580D+00
    w( 83) = 0.0238540521060385400804D+00
    w( 84) = 0.0234067774953140062013D+00
    w( 85) = 0.0229409642293877487608D+00
    w( 86) = 0.0224572658268160987071D+00
    w( 87) = 0.0219563663053178249393D+00
    w( 88) = 0.0214389800125038672465D+00
    w( 89) = 0.0209058514458120238522D+00
    w( 90) = 0.0203577550584721594669D+00
    w( 91) = 0.0197954950480974994880D+00
    w( 92) = 0.0192199051247277660193D+00
    w( 93) = 0.0186318482561387901863D+00
    w( 94) = 0.0180322163903912863201D+00
    w( 95) = 0.0174219301594641737472D+00
    w( 96) = 0.0168019385741038652709D+00
    w( 97) = 0.0161732187295777199419D+00
    w( 98) = 0.0155367755558439824399D+00
    w( 99) = 0.0148936416648151820348D+00
    w(100) = 0.0142448773729167743063D+00
    w(101) = 0.0135915710097655467896D+00
    w(102) = 0.0129348396636073734547D+00
    w(103) = 0.0122758305600827700870D+00
    w(104) = 0.0116157233199551347270D+00
    w(105) = 0.0109557333878379016480D+00
    w(106) = 0.0102971169579563555237D+00
    w(107) = 0.00964117772970253669530D+00
    w(108) = 0.00898927578406413572328D+00
    w(109) = 0.00834283875396815770558D+00
    w(110) = 0.00770337523327974184817D+00
    w(111) = 0.00707248999543355546805D+00
    w(112) = 0.00645190005017573692280D+00
    w(113) = 0.00584344987583563950756D+00
    w(114) = 0.00524912345480885912513D+00
    w(115) = 0.00467105037211432174741D+00
    w(116) = 0.00411150397865469304717D+00
    w(117) = 0.00357289278351729964938D+00
    w(118) = 0.00305775341017553113613D+00
    w(119) = 0.00256876494379402037313D+00
    w(120) = 0.00210881524572663287933D+00
    w(121) = 0.00168114286542146990631D+00
    w(122) = 0.00128952408261041739210D+00
    w(123) = 0.000938369848542381500794D+00
    w(124) = 0.000632607319362633544219D+00
    w(125) = 0.000377746646326984660274D+00
    w(126) = 0.000180739564445388357820D+00
    w(127) = 0.0000505360952078625176247D+00

  else if ( n == 255 ) then

    x(  1) = -0.99999759637974846462D+00
    x(  2) = -0.99998243035489159858D+00
    x(  3) = -0.99994399620705437576D+00
    x(  4) = -0.99987288812035761194D+00
    x(  5) = -0.99976049092443204733D+00
    x(  6) = -0.99959879967191068325D+00
    x(  7) = -0.99938033802502358193D+00
    x(  8) = -0.99909812496766759766D+00
    x(  9) = -0.99874561446809511470D+00
    x( 10) = -0.99831663531840739253D+00
    x( 11) = -0.99780535449595727456D+00
    x( 12) = -0.99720625937222195908D+00
    x( 13) = -0.99651414591489027385D+00
    x( 14) = -0.99572410469840718851D+00
    x( 15) = -0.99483150280062100052D+00
    x( 16) = -0.99383196321275502221D+00
    x( 17) = -0.99272134428278861533D+00
    x( 18) = -0.99149572117810613240D+00
    x( 19) = -0.99015137040077015918D+00
    x( 20) = -0.98868475754742947994D+00
    x( 21) = -0.98709252795403406719D+00
    x( 22) = -0.98537149959852037111D+00
    x( 23) = -0.98351865757863272876D+00
    x( 24) = -0.98153114955374010687D+00
    x( 25) = -0.97940628167086268381D+00
    x( 26) = -0.97714151463970571416D+00
    x( 27) = -0.97473445975240266776D+00
    x( 28) = -0.97218287474858179658D+00
    x( 29) = -0.96948465950245923177D+00
    x( 30) = -0.96663785155841656709D+00
    x( 31) = -0.96364062156981213252D+00
    x( 32) = -0.96049126870802028342D+00
    x( 33) = -0.95718821610986096274D+00
    x( 34) = -0.95373000642576113641D+00
    x( 35) = -0.95011529752129487656D+00
    x( 36) = -0.94634285837340290515D+00
    x( 37) = -0.94241156519108305981D+00
    x( 38) = -0.93832039777959288365D+00
    x( 39) = -0.93406843615772578800D+00
    x( 40) = -0.92965485742974005667D+00
    x( 41) = -0.92507893290707565236D+00
    x( 42) = -0.92034002547001242073D+00
    x( 43) = -0.91543758715576504064D+00
    x( 44) = -0.91037115695700429250D+00
    x( 45) = -0.90514035881326159519D+00
    x( 46) = -0.89974489977694003664D+00
    x( 47) = -0.89418456833555902286D+00
    x( 48) = -0.88845923287225699889D+00
    x( 49) = -0.88256884024734190684D+00
    x( 50) = -0.87651341448470526974D+00
    x( 51) = -0.87029305554811390585D+00
    x( 52) = -0.86390793819369047715D+00
    x( 53) = -0.85735831088623215653D+00
    x( 54) = -0.85064449476835027976D+00
    x( 55) = -0.84376688267270860104D+00
    x( 56) = -0.83672593816886873550D+00
    x( 57) = -0.82952219463740140018D+00
    x( 58) = -0.82215625436498040737D+00
    x( 59) = -0.81462878765513741344D+00
    x( 60) = -0.80694053195021761186D+00
    x( 61) = -0.79909229096084140180D+00
    x( 62) = -0.79108493379984836143D+00
    x( 63) = -0.78291939411828301639D+00
    x( 64) = -0.77459666924148337704D+00
    x( 65) = -0.76611781930376009072D+00
    x( 66) = -0.75748396638051363793D+00
    x( 67) = -0.74869629361693660282D+00
    x( 68) = -0.73975604435269475868D+00
    x( 69) = -0.73066452124218126133D+00
    x( 70) = -0.72142308537009891548D+00
    x( 71) = -0.71203315536225203459D+00
    x( 72) = -0.70249620649152707861D+00
    x( 73) = -0.69281376977911470289D+00
    x( 74) = -0.68298743109107922809D+00
    x( 75) = -0.67301883023041847920D+00
    x( 76) = -0.66290966002478059546D+00
    x( 77) = -0.65266166541001749610D+00
    x( 78) = -0.64227664250975951377D+00
    x( 79) = -0.63175643771119423041D+00
    x( 80) = -0.62110294673722640294D+00
    x( 81) = -0.61031811371518640016D+00
    x( 82) = -0.59940393024224289297D+00
    x( 83) = -0.58836243444766254143D+00
    x( 84) = -0.57719571005204581484D+00
    x( 85) = -0.56590588542365442262D+00
    x( 86) = -0.55449513263193254887D+00
    x( 87) = -0.54296566649831149049D+00
    x( 88) = -0.53131974364437562397D+00
    x( 89) = -0.51955966153745702199D+00
    x( 90) = -0.50768775753371660215D+00
    x( 91) = -0.49570640791876146017D+00
    x( 92) = -0.48361802694584102756D+00
    x( 93) = -0.47142506587165887693D+00
    x( 94) = -0.45913001198983233287D+00
    x( 95) = -0.44673538766202847374D+00
    x( 96) = -0.43424374934680255800D+00
    x( 97) = -0.42165768662616330006D+00
    x( 98) = -0.40897982122988867241D+00
    x( 99) = -0.39621280605761593918D+00
    x(100) = -0.38335932419873034692D+00
    x(101) = -0.37042208795007823014D+00
    x(102) = -0.35740383783153215238D+00
    x(103) = -0.34430734159943802278D+00
    x(104) = -0.33113539325797683309D+00
    x(105) = -0.31789081206847668318D+00
    x(106) = -0.30457644155671404334D+00
    x(107) = -0.29119514851824668196D+00
    x(108) = -0.27774982202182431507D+00
    x(109) = -0.26424337241092676194D+00
    x(110) = -0.25067873030348317661D+00
    x(111) = -0.23705884558982972721D+00
    x(112) = -0.22338668642896688163D+00
    x(113) = -0.20966523824318119477D+00
    x(114) = -0.19589750271110015392D+00
    x(115) = -0.18208649675925219825D+00
    x(116) = -0.16823525155220746498D+00
    x(117) = -0.15434681148137810869D+00
    x(118) = -0.14042423315256017459D+00
    x(119) = -0.12647058437230196685D+00
    x(120) = -0.11248894313318662575D+00
    x(121) = -0.098482396598119202090D+00
    x(122) = -0.084454040083710883710D+00
    x(123) = -0.070406976042855179063D+00
    x(124) = -0.056344313046592789972D+00
    x(125) = -0.042269164765363603212D+00
    x(126) = -0.028184648949745694339D+00
    x(127) = -0.014093886410782462614D+00
    x(128) =  0.0D+00
    x(129) =  0.014093886410782462614D+00
    x(130) =  0.028184648949745694339D+00
    x(131) =  0.042269164765363603212D+00
    x(132) =  0.056344313046592789972D+00
    x(133) =  0.070406976042855179063D+00
    x(134) =  0.084454040083710883710D+00
    x(135) =  0.098482396598119202090D+00
    x(136) =  0.11248894313318662575D+00
    x(137) =  0.12647058437230196685D+00
    x(138) =  0.14042423315256017459D+00
    x(139) =  0.15434681148137810869D+00
    x(140) =  0.16823525155220746498D+00
    x(141) =  0.18208649675925219825D+00
    x(142) =  0.19589750271110015392D+00
    x(143) =  0.20966523824318119477D+00
    x(144) =  0.22338668642896688163D+00
    x(145) =  0.23705884558982972721D+00
    x(146) =  0.25067873030348317661D+00
    x(147) =  0.26424337241092676194D+00
    x(148) =  0.27774982202182431507D+00
    x(149) =  0.29119514851824668196D+00
    x(150) =  0.30457644155671404334D+00
    x(151) =  0.31789081206847668318D+00
    x(152) =  0.33113539325797683309D+00
    x(153) =  0.34430734159943802278D+00
    x(154) =  0.35740383783153215238D+00
    x(155) =  0.37042208795007823014D+00
    x(156) =  0.38335932419873034692D+00
    x(157) =  0.39621280605761593918D+00
    x(158) =  0.40897982122988867241D+00
    x(159) =  0.42165768662616330006D+00
    x(160) =  0.43424374934680255800D+00
    x(161) =  0.44673538766202847374D+00
    x(162) =  0.45913001198983233287D+00
    x(163) =  0.47142506587165887693D+00
    x(164) =  0.48361802694584102756D+00
    x(165) =  0.49570640791876146017D+00
    x(166) =  0.50768775753371660215D+00
    x(167) =  0.51955966153745702199D+00
    x(168) =  0.53131974364437562397D+00
    x(169) =  0.54296566649831149049D+00
    x(170) =  0.55449513263193254887D+00
    x(171) =  0.56590588542365442262D+00
    x(172) =  0.57719571005204581484D+00
    x(173) =  0.58836243444766254143D+00
    x(174) =  0.59940393024224289297D+00
    x(175) =  0.61031811371518640016D+00
    x(176) =  0.62110294673722640294D+00
    x(177) =  0.63175643771119423041D+00
    x(178) =  0.64227664250975951377D+00
    x(179) =  0.65266166541001749610D+00
    x(180) =  0.66290966002478059546D+00
    x(181) =  0.67301883023041847920D+00
    x(182) =  0.68298743109107922809D+00
    x(183) =  0.69281376977911470289D+00
    x(184) =  0.70249620649152707861D+00
    x(185) =  0.71203315536225203459D+00
    x(186) =  0.72142308537009891548D+00
    x(187) =  0.73066452124218126133D+00
    x(188) =  0.73975604435269475868D+00
    x(189) =  0.74869629361693660282D+00
    x(190) =  0.75748396638051363793D+00
    x(191) =  0.76611781930376009072D+00
    x(192) =  0.77459666924148337704D+00
    x(193) =  0.78291939411828301639D+00
    x(194) =  0.79108493379984836143D+00
    x(195) =  0.79909229096084140180D+00
    x(196) =  0.80694053195021761186D+00
    x(197) =  0.81462878765513741344D+00
    x(198) =  0.82215625436498040737D+00
    x(199) =  0.82952219463740140018D+00
    x(200) =  0.83672593816886873550D+00
    x(201) =  0.84376688267270860104D+00
    x(202) =  0.85064449476835027976D+00
    x(203) =  0.85735831088623215653D+00
    x(204) =  0.86390793819369047715D+00
    x(205) =  0.87029305554811390585D+00
    x(206) =  0.87651341448470526974D+00
    x(207) =  0.88256884024734190684D+00
    x(208) =  0.88845923287225699889D+00
    x(209) =  0.89418456833555902286D+00
    x(210) =  0.89974489977694003664D+00
    x(211) =  0.90514035881326159519D+00
    x(212) =  0.91037115695700429250D+00
    x(213) =  0.91543758715576504064D+00
    x(214) =  0.92034002547001242073D+00
    x(215) =  0.92507893290707565236D+00
    x(216) =  0.92965485742974005667D+00
    x(217) =  0.93406843615772578800D+00
    x(218) =  0.93832039777959288365D+00
    x(219) =  0.94241156519108305981D+00
    x(220) =  0.94634285837340290515D+00
    x(221) =  0.95011529752129487656D+00
    x(222) =  0.95373000642576113641D+00
    x(223) =  0.95718821610986096274D+00
    x(224) =  0.96049126870802028342D+00
    x(225) =  0.96364062156981213252D+00
    x(226) =  0.96663785155841656709D+00
    x(227) =  0.96948465950245923177D+00
    x(228) =  0.97218287474858179658D+00
    x(229) =  0.97473445975240266776D+00
    x(230) =  0.97714151463970571416D+00
    x(231) =  0.97940628167086268381D+00
    x(232) =  0.98153114955374010687D+00
    x(233) =  0.98351865757863272876D+00
    x(234) =  0.98537149959852037111D+00
    x(235) =  0.98709252795403406719D+00
    x(236) =  0.98868475754742947994D+00
    x(237) =  0.99015137040077015918D+00
    x(238) =  0.99149572117810613240D+00
    x(239) =  0.99272134428278861533D+00
    x(240) =  0.99383196321275502221D+00
    x(241) =  0.99483150280062100052D+00
    x(242) =  0.99572410469840718851D+00
    x(243) =  0.99651414591489027385D+00
    x(244) =  0.99720625937222195908D+00
    x(245) =  0.99780535449595727456D+00
    x(246) =  0.99831663531840739253D+00
    x(247) =  0.99874561446809511470D+00
    x(248) =  0.99909812496766759766D+00
    x(249) =  0.99938033802502358193D+00
    x(250) =  0.99959879967191068325D+00
    x(251) =  0.99976049092443204733D+00
    x(252) =  0.99987288812035761194D+00
    x(253) =  0.99994399620705437576D+00
    x(254) =  0.99998243035489159858D+00
    x(255) =  0.99999759637974846462D+00

    w(  1) = 0.69379364324108267170D-05
    w(  2) = 0.25157870384280661489D-04
    w(  3) = 0.53275293669780613125D-04
    w(  4) = 0.90372734658751149261D-04
    w(  5) = 0.13575491094922871973D-03
    w(  6) = 0.18887326450650491366D-03
    w(  7) = 0.24921240048299729402D-03
    w(  8) = 0.31630366082226447689D-03
    w(  9) = 0.38974528447328229322D-03
    w( 10) = 0.46918492424785040975D-03
    w( 11) = 0.55429531493037471492D-03
    w( 12) = 0.64476204130572477933D-03
    w( 13) = 0.74028280424450333046D-03
    w( 14) = 0.84057143271072246365D-03
    w( 15) = 0.94536151685852538246D-03
    w( 16) = 0.10544076228633167722D-02
    w( 17) = 0.11674841174299594077D-02
    w( 18) = 0.12843824718970101768D-02
    w( 19) = 0.14049079956551446427D-02
    w( 20) = 0.15288767050877655684D-02
    w( 21) = 0.16561127281544526052D-02
    w( 22) = 0.17864463917586498247D-02
    w( 23) = 0.19197129710138724125D-02
    w( 24) = 0.20557519893273465236D-02
    w( 25) = 0.21944069253638388388D-02
    w( 26) = 0.23355251860571608737D-02
    w( 27) = 0.24789582266575679307D-02
    w( 28) = 0.26245617274044295626D-02
    w( 29) = 0.27721957645934509940D-02
    w( 30) = 0.29217249379178197538D-02
    w( 31) = 0.30730184347025783234D-02
    w( 32) = 0.32259500250878684614D-02
    w( 33) = 0.33803979910869203823D-02
    w( 34) = 0.35362449977167777340D-02
    w( 35) = 0.36933779170256508183D-02
    w( 36) = 0.38516876166398709241D-02
    w( 37) = 0.40110687240750233989D-02
    w( 38) = 0.41714193769840788528D-02
    w( 39) = 0.43326409680929828545D-02
    w( 40) = 0.44946378920320678616D-02
    w( 41) = 0.46573172997568547773D-02
    w( 42) = 0.48205888648512683476D-02
    w( 43) = 0.49843645647655386012D-02
    w( 44) = 0.51485584789781777618D-02
    w( 45) = 0.53130866051870565663D-02
    w( 46) = 0.54778666939189508240D-02
    w( 47) = 0.56428181013844441585D-02
    w( 48) = 0.58078616599775673635D-02
    w( 49) = 0.59729195655081658049D-02
    w( 50) = 0.61379152800413850435D-02
    w( 51) = 0.63027734490857587172D-02
    w( 52) = 0.64674198318036867274D-02
    w( 53) = 0.66317812429018878941D-02
    w( 54) = 0.67957855048827733948D-02
    w( 55) = 0.69593614093904229394D-02
    w( 56) = 0.71224386864583871532D-02
    w( 57) = 0.72849479805538070639D-02
    w( 58) = 0.74468208324075910174D-02
    w( 59) = 0.76079896657190565832D-02
    w( 60) = 0.77683877779219912200D-02
    w( 61) = 0.79279493342948491103D-02
    w( 62) = 0.80866093647888599710D-02
    w( 63) = 0.82443037630328680306D-02
    w( 64) = 0.84009692870519326354D-02
    w( 65) = 0.85565435613076896192D-02
    w( 66) = 0.87109650797320868736D-02
    w( 67) = 0.88641732094824942641D-02
    w( 68) = 0.90161081951956431600D-02
    w( 69) = 0.91667111635607884067D-02
    w( 70) = 0.93159241280693950932D-02
    w( 71) = 0.94636899938300652943D-02
    w( 72) = 0.96099525623638830097D-02
    w( 73) = 0.97546565363174114611D-02
    w( 74) = 0.98977475240487497440D-02
    w( 75) = 0.10039172044056840798D-01
    w( 76) = 0.10178877529236079733D-01
    w( 77) = 0.10316812330947621682D-01
    w( 78) = 0.10452925722906011926D-01
    w( 79) = 0.10587167904885197931D-01
    w( 80) = 0.10719490006251933623D-01
    w( 81) = 0.10849844089337314099D-01
    w( 82) = 0.10978183152658912470D-01
    w( 83) = 0.11104461134006926537D-01
    w( 84) = 0.11228632913408049354D-01
    w( 85) = 0.11350654315980596602D-01
    w( 86) = 0.11470482114693874380D-01
    w( 87) = 0.11588074033043952568D-01
    w( 88) = 0.11703388747657003101D-01
    w( 89) = 0.11816385890830235763D-01
    w( 90) = 0.11927026053019270040D-01
    w( 91) = 0.12035270785279562630D-01
    w( 92) = 0.12141082601668299679D-01
    w( 93) = 0.12244424981611985899D-01
    w( 94) = 0.12345262372243838455D-01
    w( 95) = 0.12443560190714035263D-01
    w( 96) = 0.12539284826474884353D-01
    w( 97) = 0.12632403643542078765D-01
    w( 98) = 0.12722884982732382906D-01
    w( 99) = 0.12810698163877361967D-01
    w(100) = 0.12895813488012114694D-01
    w(101) = 0.12978202239537399286D-01
    w(102) = 0.13057836688353048840D-01
    w(103) = 0.13134690091960152836D-01
    w(104) = 0.13208736697529129966D-01
    w(105) = 0.13279951743930530650D-01
    w(106) = 0.13348311463725179953D-01
    w(107) = 0.13413793085110098513D-01
    w(108) = 0.13476374833816515982D-01
    w(109) = 0.13536035934956213614D-01
    w(110) = 0.13592756614812395910D-01
    w(111) = 0.13646518102571291428D-01
    w(112) = 0.13697302631990716258D-01
    w(113) = 0.13745093443001896632D-01
    w(114) = 0.13789874783240936517D-01
    w(115) = 0.13831631909506428676D-01
    w(116) = 0.13870351089139840997D-01
    w(117) = 0.13906019601325461264D-01
    w(118) = 0.13938625738306850804D-01
    w(119) = 0.13968158806516938516D-01
    w(120) = 0.13994609127619079852D-01
    w(121) = 0.14017968039456608810D-01
    w(122) = 0.14038227896908623303D-01
    w(123) = 0.14055382072649964277D-01
    w(124) = 0.14069424957813575318D-01
    w(125) = 0.14080351962553661325D-01
    w(126) = 0.14088159516508301065D-01
    w(127) = 0.14092845069160408355D-01
    w(128) = 0.14094407090096179347D-01
    w(129) = 0.14092845069160408355D-01
    w(130) = 0.14088159516508301065D-01
    w(131) = 0.14080351962553661325D-01
    w(132) = 0.14069424957813575318D-01
    w(133) = 0.14055382072649964277D-01
    w(134) = 0.14038227896908623303D-01
    w(135) = 0.14017968039456608810D-01
    w(136) = 0.13994609127619079852D-01
    w(137) = 0.13968158806516938516D-01
    w(138) = 0.13938625738306850804D-01
    w(139) = 0.13906019601325461264D-01
    w(140) = 0.13870351089139840997D-01
    w(141) = 0.13831631909506428676D-01
    w(142) = 0.13789874783240936517D-01
    w(143) = 0.13745093443001896632D-01
    w(144) = 0.13697302631990716258D-01
    w(145) = 0.13646518102571291428D-01
    w(146) = 0.13592756614812395910D-01
    w(147) = 0.13536035934956213614D-01
    w(148) = 0.13476374833816515982D-01
    w(149) = 0.13413793085110098513D-01
    w(150) = 0.13348311463725179953D-01
    w(151) = 0.13279951743930530650D-01
    w(152) = 0.13208736697529129966D-01
    w(153) = 0.13134690091960152836D-01
    w(154) = 0.13057836688353048840D-01
    w(155) = 0.12978202239537399286D-01
    w(156) = 0.12895813488012114694D-01
    w(157) = 0.12810698163877361967D-01
    w(158) = 0.12722884982732382906D-01
    w(159) = 0.12632403643542078765D-01
    w(160) = 0.12539284826474884353D-01
    w(161) = 0.12443560190714035263D-01
    w(162) = 0.12345262372243838455D-01
    w(163) = 0.12244424981611985899D-01
    w(164) = 0.12141082601668299679D-01
    w(165) = 0.12035270785279562630D-01
    w(166) = 0.11927026053019270040D-01
    w(167) = 0.11816385890830235763D-01
    w(168) = 0.11703388747657003101D-01
    w(169) = 0.11588074033043952568D-01
    w(170) = 0.11470482114693874380D-01
    w(171) = 0.11350654315980596602D-01
    w(172) = 0.11228632913408049354D-01
    w(173) = 0.11104461134006926537D-01
    w(174) = 0.10978183152658912470D-01
    w(175) = 0.10849844089337314099D-01
    w(176) = 0.10719490006251933623D-01
    w(177) = 0.10587167904885197931D-01
    w(178) = 0.10452925722906011926D-01
    w(179) = 0.10316812330947621682D-01
    w(180) = 0.10178877529236079733D-01
    w(181) = 0.10039172044056840798D-01
    w(182) = 0.98977475240487497440D-02
    w(183) = 0.97546565363174114611D-02
    w(184) = 0.96099525623638830097D-02
    w(185) = 0.94636899938300652943D-02
    w(186) = 0.93159241280693950932D-02
    w(187) = 0.91667111635607884067D-02
    w(188) = 0.90161081951956431600D-02
    w(189) = 0.88641732094824942641D-02
    w(190) = 0.87109650797320868736D-02
    w(191) = 0.85565435613076896192D-02
    w(192) = 0.84009692870519326354D-02
    w(193) = 0.82443037630328680306D-02
    w(194) = 0.80866093647888599710D-02
    w(195) = 0.79279493342948491103D-02
    w(196) = 0.77683877779219912200D-02
    w(197) = 0.76079896657190565832D-02
    w(198) = 0.74468208324075910174D-02
    w(199) = 0.72849479805538070639D-02
    w(200) = 0.71224386864583871532D-02
    w(201) = 0.69593614093904229394D-02
    w(202) = 0.67957855048827733948D-02
    w(203) = 0.66317812429018878941D-02
    w(204) = 0.64674198318036867274D-02
    w(205) = 0.63027734490857587172D-02
    w(206) = 0.61379152800413850435D-02
    w(207) = 0.59729195655081658049D-02
    w(208) = 0.58078616599775673635D-02
    w(209) = 0.56428181013844441585D-02
    w(210) = 0.54778666939189508240D-02
    w(211) = 0.53130866051870565663D-02
    w(212) = 0.51485584789781777618D-02
    w(213) = 0.49843645647655386012D-02
    w(214) = 0.48205888648512683476D-02
    w(215) = 0.46573172997568547773D-02
    w(216) = 0.44946378920320678616D-02
    w(217) = 0.43326409680929828545D-02
    w(218) = 0.41714193769840788528D-02
    w(219) = 0.40110687240750233989D-02
    w(220) = 0.38516876166398709241D-02
    w(221) = 0.36933779170256508183D-02
    w(222) = 0.35362449977167777340D-02
    w(223) = 0.33803979910869203823D-02
    w(224) = 0.32259500250878684614D-02
    w(225) = 0.30730184347025783234D-02
    w(226) = 0.29217249379178197538D-02
    w(227) = 0.27721957645934509940D-02
    w(228) = 0.26245617274044295626D-02
    w(229) = 0.24789582266575679307D-02
    w(230) = 0.23355251860571608737D-02
    w(231) = 0.21944069253638388388D-02
    w(232) = 0.20557519893273465236D-02
    w(233) = 0.19197129710138724125D-02
    w(234) = 0.17864463917586498247D-02
    w(235) = 0.16561127281544526052D-02
    w(236) = 0.15288767050877655684D-02
    w(237) = 0.14049079956551446427D-02
    w(238) = 0.12843824718970101768D-02
    w(239) = 0.11674841174299594077D-02
    w(240) = 0.10544076228633167722D-02
    w(241) = 0.94536151685852538246D-03
    w(242) = 0.84057143271072246365D-03
    w(243) = 0.74028280424450333046D-03
    w(244) = 0.64476204130572477933D-03
    w(245) = 0.55429531493037471492D-03
    w(246) = 0.46918492424785040975D-03
    w(247) = 0.38974528447328229322D-03
    w(248) = 0.31630366082226447689D-03
    w(249) = 0.24921240048299729402D-03
    w(250) = 0.18887326450650491366D-03
    w(251) = 0.13575491094922871973D-03
    w(252) = 0.90372734658751149261D-04
    w(253) = 0.53275293669780613125D-04
    w(254) = 0.25157870384280661489D-04
    w(255) = 0.69379364324108267170D-05

  else if ( n == 511 ) then

    x(  1) = -0.999999672956734384381D+00
    x(  2) = -0.999997596379748464620D+00
    x(  3) = -0.999992298136257588028D+00
    x(  4) = -0.999982430354891598580D+00
    x(  5) = -0.999966730098486276883D+00
    x(  6) = -0.999943996207054375764D+00
    x(  7) = -0.999913081144678282800D+00
    x(  8) = -0.999872888120357611938D+00
    x(  9) = -0.999822363679787739196D+00
    x( 10) = -0.999760490924432047330D+00
    x( 11) = -0.999686286448317731776D+00
    x( 12) = -0.999598799671910683252D+00
    x( 13) = -0.999497112467187190535D+00
    x( 14) = -0.999380338025023581928D+00
    x( 15) = -0.999247618943342473599D+00
    x( 16) = -0.999098124967667597662D+00
    x( 17) = -0.998931050830810562236D+00
    x( 18) = -0.998745614468095114704D+00
    x( 19) = -0.998541055697167906027D+00
    x( 20) = -0.998316635318407392531D+00
    x( 21) = -0.998071634524930323302D+00
    x( 22) = -0.997805354495957274562D+00
    x( 23) = -0.997517116063472399965D+00
    x( 24) = -0.997206259372221959076D+00
    x( 25) = -0.996872143485260161299D+00
    x( 26) = -0.996514145914890273849D+00
    x( 27) = -0.996131662079315037786D+00
    x( 28) = -0.995724104698407188509D+00
    x( 29) = -0.995290903148810302261D+00
    x( 30) = -0.994831502800621000519D+00
    x( 31) = -0.994345364356723405931D+00
    x( 32) = -0.993831963212755022209D+00
    x( 33) = -0.993290788851684966211D+00
    x( 34) = -0.992721344282788615328D+00
    x( 35) = -0.992123145530863117683D+00
    x( 36) = -0.991495721178106132399D+00
    x( 37) = -0.990838611958294243677D+00
    x( 38) = -0.990151370400770159181D+00
    x( 39) = -0.989433560520240838716D+00
    x( 40) = -0.988684757547429479939D+00
    x( 41) = -0.987904547695124280467D+00
    x( 42) = -0.987092527954034067190D+00
    x( 43) = -0.986248305913007552681D+00
    x( 44) = -0.985371499598520371114D+00
    x( 45) = -0.984461737328814534596D+00
    x( 46) = -0.983518657578632728762D+00
    x( 47) = -0.982541908851080604251D+00
    x( 48) = -0.981531149553740106867D+00
    x( 49) = -0.980486047876721339416D+00
    x( 50) = -0.979406281670862683806D+00
    x( 51) = -0.978291538324758539526D+00
    x( 52) = -0.977141514639705714156D+00
    x( 53) = -0.975955916702011753129D+00
    x( 54) = -0.974734459752402667761D+00
    x( 55) = -0.973476868052506926773D+00
    x( 56) = -0.972182874748581796578D+00
    x( 57) = -0.970852221732792443256D+00
    x( 58) = -0.969484659502459231771D+00
    x( 59) = -0.968079947017759947964D+00
    x( 60) = -0.966637851558416567092D+00
    x( 61) = -0.965158148579915665979D+00
    x( 62) = -0.963640621569812132521D+00
    x( 63) = -0.962085061904651475741D+00
    x( 64) = -0.960491268708020283423D+00
    x( 65) = -0.958859048710200221356D+00
    x( 66) = -0.957188216109860962736D+00
    x( 67) = -0.955478592438183697574D+00
    x( 68) = -0.953730006425761136415D+00
    x( 69) = -0.951942293872573589498D+00
    x( 70) = -0.950115297521294876558D+00
    x( 71) = -0.948248866934137357063D+00
    x( 72) = -0.946342858373402905148D+00
    x( 73) = -0.944397134685866648591D+00
    x( 74) = -0.942411565191083059813D+00
    x( 75) = -0.940386025573669721370D+00
    x( 76) = -0.938320397779592883655D+00
    x( 77) = -0.936214569916450806625D+00
    x( 78) = -0.934068436157725787999D+00
    x( 79) = -0.931881896650953639345D+00
    x( 80) = -0.929654857429740056670D+00
    x( 81) = -0.927387230329536696843D+00
    x( 82) = -0.925078932907075652364D+00
    x( 83) = -0.922729888363349241523D+00
    x( 84) = -0.920340025470012420730D+00
    x( 85) = -0.917909278499077501636D+00
    x( 86) = -0.915437587155765040644D+00
    x( 87) = -0.912924896514370590080D+00
    x( 88) = -0.910371156957004292498D+00
    x( 89) = -0.907776324115058903624D+00
    x( 90) = -0.905140358813261595189D+00
    x( 91) = -0.902463227016165675048D+00
    x( 92) = -0.899744899776940036639D+00
    x( 93) = -0.896985353188316590376D+00
    x( 94) = -0.894184568335559022859D+00
    x( 95) = -0.891342531251319871666D+00
    x( 96) = -0.888459232872256998890D+00
    x( 97) = -0.885534668997285008926D+00
    x( 98) = -0.882568840247341906842D+00
    x( 99) = -0.879561752026556262568D+00
    x(100) = -0.876513414484705269742D+00
    x(101) = -0.873423842480859310192D+00
    x(102) = -0.870293055548113905851D+00
    x(103) = -0.867121077859315215614D+00
    x(104) = -0.863907938193690477146D+00
    x(105) = -0.860653669904299969802D+00
    x(106) = -0.857358310886232156525D+00
    x(107) = -0.854021903545468625813D+00
    x(108) = -0.850644494768350279758D+00
    x(109) = -0.847226135891580884381D+00
    x(110) = -0.843766882672708601038D+00
    x(111) = -0.840266795261030442350D+00
    x(112) = -0.836725938168868735503D+00
    x(113) = -0.833144380243172624728D+00
    x(114) = -0.829522194637401400178D+00
    x(115) = -0.825859458783650001088D+00
    x(116) = -0.822156254364980407373D+00
    x(117) = -0.818412667287925807395D+00
    x(118) = -0.814628787655137413436D+00
    x(119) = -0.810804709738146594361D+00
    x(120) = -0.806940531950217611856D+00
    x(121) = -0.803036356819268687782D+00
    x(122) = -0.799092290960841401800D+00
    x(123) = -0.795108445051100526780D+00
    x(124) = -0.791084933799848361435D+00
    x(125) = -0.787021875923539422170D+00
    x(126) = -0.782919394118283016385D+00
    x(127) = -0.778777615032822744702D+00
    x(128) = -0.774596669241483377036D+00
    x(129) = -0.770376691217076824278D+00
    x(130) = -0.766117819303760090717D+00
    x(131) = -0.761820195689839149173D+00
    x(132) = -0.757483966380513637926D+00
    x(133) = -0.753109281170558142523D+00
    x(134) = -0.748696293616936602823D+00
    x(135) = -0.744245161011347082309D+00
    x(136) = -0.739756044352694758677D+00
    x(137) = -0.735229108319491547663D+00
    x(138) = -0.730664521242181261329D+00
    x(139) = -0.726062455075389632685D+00
    x(140) = -0.721423085370098915485D+00
    x(141) = -0.716746591245747095767D+00
    x(142) = -0.712033155362252034587D+00
    x(143) = -0.707282963891961103412D+00
    x(144) = -0.702496206491527078610D+00
    x(145) = -0.697673076273711232906D+00
    x(146) = -0.692813769779114702895D+00
    x(147) = -0.687918486947839325756D+00
    x(148) = -0.682987431091079228087D+00
    x(149) = -0.678020808862644517838D+00
    x(150) = -0.673018830230418479199D+00
    x(151) = -0.667981708447749702165D+00
    x(152) = -0.662909660024780595461D+00
    x(153) = -0.657802904699713735422D+00
    x(154) = -0.652661665410017496101D+00
    x(155) = -0.647486168263572388782D+00
    x(156) = -0.642276642509759513774D+00
    x(157) = -0.637033320510492495071D+00
    x(158) = -0.631756437711194230414D+00
    x(159) = -0.626446232611719746542D+00
    x(160) = -0.621102946737226402941D+00
    x(161) = -0.615726824608992638014D+00
    x(162) = -0.610318113715186400156D+00
    x(163) = -0.604877064481584353319D+00
    x(164) = -0.599403930242242892974D+00
    x(165) = -0.593898967210121954393D+00
    x(166) = -0.588362434447662541434D+00
    x(167) = -0.582794593837318850840D+00
    x(168) = -0.577195710052045814844D+00
    x(169) = -0.571566050525742833992D+00
    x(170) = -0.565905885423654422623D+00
    x(171) = -0.560215487612728441818D+00
    x(172) = -0.554495132631932548866D+00
    x(173) = -0.548745098662529448608D+00
    x(174) = -0.542965666498311490492D+00
    x(175) = -0.537157119515795115982D+00
    x(176) = -0.531319743644375623972D+00
    x(177) = -0.525453827336442687395D+00
    x(178) = -0.519559661537457021993D+00
    x(179) = -0.513637539655988578507D+00
    x(180) = -0.507687757533716602155D+00
    x(181) = -0.501710613415391878251D+00
    x(182) = -0.495706407918761460170D+00
    x(183) = -0.489675444004456155436D+00
    x(184) = -0.483618026945841027562D+00
    x(185) = -0.477534464298829155284D+00
    x(186) = -0.471425065871658876934D+00
    x(187) = -0.465290143694634735858D+00
    x(188) = -0.459130011989832332874D+00
    x(189) = -0.452944987140767283784D+00
    x(190) = -0.446735387662028473742D+00
    x(191) = -0.440501534168875795783D+00
    x(192) = -0.434243749346802558002D+00
    x(193) = -0.427962357921062742583D+00
    x(194) = -0.421657686626163300056D+00
    x(195) = -0.415330064175321663764D+00
    x(196) = -0.408979821229888672409D+00
    x(197) = -0.402607290368737092671D+00
    x(198) = -0.396212806057615939183D+00
    x(199) = -0.389796704618470795479D+00
    x(200) = -0.383359324198730346916D+00
    x(201) = -0.376901004740559344802D+00
    x(202) = -0.370422087950078230138D+00
    x(203) = -0.363922917266549655269D+00
    x(204) = -0.357403837831532152376D+00
    x(205) = -0.350865196458001209011D+00
    x(206) = -0.344307341599438022777D+00
    x(207) = -0.337730623318886219621D+00
    x(208) = -0.331135393257976833093D+00
    x(209) = -0.324522004605921855207D+00
    x(210) = -0.317890812068476683182D+00
    x(211) = -0.311242171836871800300D+00
    x(212) = -0.304576441556714043335D+00
    x(213) = -0.297893980296857823437D+00
    x(214) = -0.291195148518246681964D+00
    x(215) = -0.284480308042725577496D+00
    x(216) = -0.277749822021824315065D+00
    x(217) = -0.271004054905512543536D+00
    x(218) = -0.264243372410926761945D+00
    x(219) = -0.257468141491069790481D+00
    x(220) = -0.250678730303483176613D+00
    x(221) = -0.243875508178893021593D+00
    x(222) = -0.237058845589829727213D+00
    x(223) = -0.230229114119222177156D+00
    x(224) = -0.223386686428966881628D+00
    x(225) = -0.216531936228472628081D+00
    x(226) = -0.209665238243181194766D+00
    x(227) = -0.202786968183064697557D+00
    x(228) = -0.195897502711100153915D+00
    x(229) = -0.188997219411721861059D+00
    x(230) = -0.182086496759252198246D+00
    x(231) = -0.175165714086311475707D+00
    x(232) = -0.168235251552207464982D+00
    x(233) = -0.161295490111305257361D+00
    x(234) = -0.154346811481378108692D+00
    x(235) = -0.147389598111939940054D+00
    x(236) = -0.140424233152560174594D+00
    x(237) = -0.133451100421161601344D+00
    x(238) = -0.126470584372301966851D+00
    x(239) = -0.119483070065440005133D+00
    x(240) = -0.112488943133186625746D+00
    x(241) = -0.105488589749541988533D+00
    x(242) = -0.984823965981192020903D-01
    x(243) = -0.914707508403553909095D-01
    x(244) = -0.844540400837108837102D-01
    x(245) = -0.774326523498572825675D-01
    x(246) = -0.704069760428551790633D-01
    x(247) = -0.633773999173222898797D-01
    x(248) = -0.563443130465927899720D-01
    x(249) = -0.493081047908686267156D-01
    x(250) = -0.422691647653636032124D-01
    x(251) = -0.352278828084410232603D-01
    x(252) = -0.281846489497456943394D-01
    x(253) = -0.211398533783310883350D-01
    x(254) = -0.140938864107824626142D-01
    x(255) = -0.704713845933674648514D-02
    x(256) = +0.000000000000000000000D+00
    x(257) = +0.704713845933674648514D-02
    x(258) = +0.140938864107824626142D-01
    x(259) = +0.211398533783310883350D-01
    x(260) = +0.281846489497456943394D-01
    x(261) = +0.352278828084410232603D-01
    x(262) = +0.422691647653636032124D-01
    x(263) = +0.493081047908686267156D-01
    x(264) = +0.563443130465927899720D-01
    x(265) = +0.633773999173222898797D-01
    x(266) = +0.704069760428551790633D-01
    x(267) = +0.774326523498572825675D-01
    x(268) = +0.844540400837108837102D-01
    x(269) = +0.914707508403553909095D-01
    x(270) = +0.984823965981192020903D-01
    x(271) = +0.105488589749541988533D+00
    x(272) = +0.112488943133186625746D+00
    x(273) = +0.119483070065440005133D+00
    x(274) = +0.126470584372301966851D+00
    x(275) = +0.133451100421161601344D+00
    x(276) = +0.140424233152560174594D+00
    x(277) = +0.147389598111939940054D+00
    x(278) = +0.154346811481378108692D+00
    x(279) = +0.161295490111305257361D+00
    x(280) = +0.168235251552207464982D+00
    x(281) = +0.175165714086311475707D+00
    x(282) = +0.182086496759252198246D+00
    x(283) = +0.188997219411721861059D+00
    x(284) = +0.195897502711100153915D+00
    x(285) = +0.202786968183064697557D+00
    x(286) = +0.209665238243181194766D+00
    x(287) = +0.216531936228472628081D+00
    x(288) = +0.223386686428966881628D+00
    x(289) = +0.230229114119222177156D+00
    x(290) = +0.237058845589829727213D+00
    x(291) = +0.243875508178893021593D+00
    x(292) = +0.250678730303483176613D+00
    x(293) = +0.257468141491069790481D+00
    x(294) = +0.264243372410926761945D+00
    x(295) = +0.271004054905512543536D+00
    x(296) = +0.277749822021824315065D+00
    x(297) = +0.284480308042725577496D+00
    x(298) = +0.291195148518246681964D+00
    x(299) = +0.297893980296857823437D+00
    x(300) = +0.304576441556714043335D+00
    x(301) = +0.311242171836871800300D+00
    x(302) = +0.317890812068476683182D+00
    x(303) = +0.324522004605921855207D+00
    x(304) = +0.331135393257976833093D+00
    x(305) = +0.337730623318886219621D+00
    x(306) = +0.344307341599438022777D+00
    x(307) = +0.350865196458001209011D+00
    x(308) = +0.357403837831532152376D+00
    x(309) = +0.363922917266549655269D+00
    x(310) = +0.370422087950078230138D+00
    x(311) = +0.376901004740559344802D+00
    x(312) = +0.383359324198730346916D+00
    x(313) = +0.389796704618470795479D+00
    x(314) = +0.396212806057615939183D+00
    x(315) = +0.402607290368737092671D+00
    x(316) = +0.408979821229888672409D+00
    x(317) = +0.415330064175321663764D+00
    x(318) = +0.421657686626163300056D+00
    x(319) = +0.427962357921062742583D+00
    x(320) = +0.434243749346802558002D+00
    x(321) = +0.440501534168875795783D+00
    x(322) = +0.446735387662028473742D+00
    x(323) = +0.452944987140767283784D+00
    x(324) = +0.459130011989832332874D+00
    x(325) = +0.465290143694634735858D+00
    x(326) = +0.471425065871658876934D+00
    x(327) = +0.477534464298829155284D+00
    x(328) = +0.483618026945841027562D+00
    x(329) = +0.489675444004456155436D+00
    x(330) = +0.495706407918761460170D+00
    x(331) = +0.501710613415391878251D+00
    x(332) = +0.507687757533716602155D+00
    x(333) = +0.513637539655988578507D+00
    x(334) = +0.519559661537457021993D+00
    x(335) = +0.525453827336442687395D+00
    x(336) = +0.531319743644375623972D+00
    x(337) = +0.537157119515795115982D+00
    x(338) = +0.542965666498311490492D+00
    x(339) = +0.548745098662529448608D+00
    x(340) = +0.554495132631932548866D+00
    x(341) = +0.560215487612728441818D+00
    x(342) = +0.565905885423654422623D+00
    x(343) = +0.571566050525742833992D+00
    x(344) = +0.577195710052045814844D+00
    x(345) = +0.582794593837318850840D+00
    x(346) = +0.588362434447662541434D+00
    x(347) = +0.593898967210121954393D+00
    x(348) = +0.599403930242242892974D+00
    x(349) = +0.604877064481584353319D+00
    x(350) = +0.610318113715186400156D+00
    x(351) = +0.615726824608992638014D+00
    x(352) = +0.621102946737226402941D+00
    x(353) = +0.626446232611719746542D+00
    x(354) = +0.631756437711194230414D+00
    x(355) = +0.637033320510492495071D+00
    x(356) = +0.642276642509759513774D+00
    x(357) = +0.647486168263572388782D+00
    x(358) = +0.652661665410017496101D+00
    x(359) = +0.657802904699713735422D+00
    x(360) = +0.662909660024780595461D+00
    x(361) = +0.667981708447749702165D+00
    x(362) = +0.673018830230418479199D+00
    x(363) = +0.678020808862644517838D+00
    x(364) = +0.682987431091079228087D+00
    x(365) = +0.687918486947839325756D+00
    x(366) = +0.692813769779114702895D+00
    x(367) = +0.697673076273711232906D+00
    x(368) = +0.702496206491527078610D+00
    x(369) = +0.707282963891961103412D+00
    x(370) = +0.712033155362252034587D+00
    x(371) = +0.716746591245747095767D+00
    x(372) = +0.721423085370098915485D+00
    x(373) = +0.726062455075389632685D+00
    x(374) = +0.730664521242181261329D+00
    x(375) = +0.735229108319491547663D+00
    x(376) = +0.739756044352694758677D+00
    x(377) = +0.744245161011347082309D+00
    x(378) = +0.748696293616936602823D+00
    x(379) = +0.753109281170558142523D+00
    x(380) = +0.757483966380513637926D+00
    x(381) = +0.761820195689839149173D+00
    x(382) = +0.766117819303760090717D+00
    x(383) = +0.770376691217076824278D+00
    x(384) = +0.774596669241483377036D+00
    x(385) = +0.778777615032822744702D+00
    x(386) = +0.782919394118283016385D+00
    x(387) = +0.787021875923539422170D+00
    x(388) = +0.791084933799848361435D+00
    x(389) = +0.795108445051100526780D+00
    x(390) = +0.799092290960841401800D+00
    x(391) = +0.803036356819268687782D+00
    x(392) = +0.806940531950217611856D+00
    x(393) = +0.810804709738146594361D+00
    x(394) = +0.814628787655137413436D+00
    x(395) = +0.818412667287925807395D+00
    x(396) = +0.822156254364980407373D+00
    x(397) = +0.825859458783650001088D+00
    x(398) = +0.829522194637401400178D+00
    x(399) = +0.833144380243172624728D+00
    x(400) = +0.836725938168868735503D+00
    x(401) = +0.840266795261030442350D+00
    x(402) = +0.843766882672708601038D+00
    x(403) = +0.847226135891580884381D+00
    x(404) = +0.850644494768350279758D+00
    x(405) = +0.854021903545468625813D+00
    x(406) = +0.857358310886232156525D+00
    x(407) = +0.860653669904299969802D+00
    x(408) = +0.863907938193690477146D+00
    x(409) = +0.867121077859315215614D+00
    x(410) = +0.870293055548113905851D+00
    x(411) = +0.873423842480859310192D+00
    x(412) = +0.876513414484705269742D+00
    x(413) = +0.879561752026556262568D+00
    x(414) = +0.882568840247341906842D+00
    x(415) = +0.885534668997285008926D+00
    x(416) = +0.888459232872256998890D+00
    x(417) = +0.891342531251319871666D+00
    x(418) = +0.894184568335559022859D+00
    x(419) = +0.896985353188316590376D+00
    x(420) = +0.899744899776940036639D+00
    x(421) = +0.902463227016165675048D+00
    x(422) = +0.905140358813261595189D+00
    x(423) = +0.907776324115058903624D+00
    x(424) = +0.910371156957004292498D+00
    x(425) = +0.912924896514370590080D+00
    x(426) = +0.915437587155765040644D+00
    x(427) = +0.917909278499077501636D+00
    x(428) = +0.920340025470012420730D+00
    x(429) = +0.922729888363349241523D+00
    x(430) = +0.925078932907075652364D+00
    x(431) = +0.927387230329536696843D+00
    x(432) = +0.929654857429740056670D+00
    x(433) = +0.931881896650953639345D+00
    x(434) = +0.934068436157725787999D+00
    x(435) = +0.936214569916450806625D+00
    x(436) = +0.938320397779592883655D+00
    x(437) = +0.940386025573669721370D+00
    x(438) = +0.942411565191083059813D+00
    x(439) = +0.944397134685866648591D+00
    x(440) = +0.946342858373402905148D+00
    x(441) = +0.948248866934137357063D+00
    x(442) = +0.950115297521294876558D+00
    x(443) = +0.951942293872573589498D+00
    x(444) = +0.953730006425761136415D+00
    x(445) = +0.955478592438183697574D+00
    x(446) = +0.957188216109860962736D+00
    x(447) = +0.958859048710200221356D+00
    x(448) = +0.960491268708020283423D+00
    x(449) = +0.962085061904651475741D+00
    x(450) = +0.963640621569812132521D+00
    x(451) = +0.965158148579915665979D+00
    x(452) = +0.966637851558416567092D+00
    x(453) = +0.968079947017759947964D+00
    x(454) = +0.969484659502459231771D+00
    x(455) = +0.970852221732792443256D+00
    x(456) = +0.972182874748581796578D+00
    x(457) = +0.973476868052506926773D+00
    x(458) = +0.974734459752402667761D+00
    x(459) = +0.975955916702011753129D+00
    x(460) = +0.977141514639705714156D+00
    x(461) = +0.978291538324758539526D+00
    x(462) = +0.979406281670862683806D+00
    x(463) = +0.980486047876721339416D+00
    x(464) = +0.981531149553740106867D+00
    x(465) = +0.982541908851080604251D+00
    x(466) = +0.983518657578632728762D+00
    x(467) = +0.984461737328814534596D+00
    x(468) = +0.985371499598520371114D+00
    x(469) = +0.986248305913007552681D+00
    x(470) = +0.987092527954034067190D+00
    x(471) = +0.987904547695124280467D+00
    x(472) = +0.988684757547429479939D+00
    x(473) = +0.989433560520240838716D+00
    x(474) = +0.990151370400770159181D+00
    x(475) = +0.990838611958294243677D+00
    x(476) = +0.991495721178106132399D+00
    x(477) = +0.992123145530863117683D+00
    x(478) = +0.992721344282788615328D+00
    x(479) = +0.993290788851684966211D+00
    x(480) = +0.993831963212755022209D+00
    x(481) = +0.994345364356723405931D+00
    x(482) = +0.994831502800621000519D+00
    x(483) = +0.995290903148810302261D+00
    x(484) = +0.995724104698407188509D+00
    x(485) = +0.996131662079315037786D+00
    x(486) = +0.996514145914890273849D+00
    x(487) = +0.996872143485260161299D+00
    x(488) = +0.997206259372221959076D+00
    x(489) = +0.997517116063472399965D+00
    x(490) = +0.997805354495957274562D+00
    x(491) = +0.998071634524930323302D+00
    x(492) = +0.998316635318407392531D+00
    x(493) = +0.998541055697167906027D+00
    x(494) = +0.998745614468095114704D+00
    x(495) = +0.998931050830810562236D+00
    x(496) = +0.999098124967667597662D+00
    x(497) = +0.999247618943342473599D+00
    x(498) = +0.999380338025023581928D+00
    x(499) = +0.999497112467187190535D+00
    x(500) = +0.999598799671910683252D+00
    x(501) = +0.999686286448317731776D+00
    x(502) = +0.999760490924432047330D+00
    x(503) = +0.999822363679787739196D+00
    x(504) = +0.999872888120357611938D+00
    x(505) = +0.999913081144678282800D+00
    x(506) = +0.999943996207054375764D+00
    x(507) = +0.999966730098486276883D+00
    x(508) = +0.999982430354891598580D+00
    x(509) = +0.999992298136257588028D+00
    x(510) = +0.999997596379748464620D+00
    x(511) = +0.999999672956734384381D+00

    w(  1) = 0.945715933950007048827D-06
    w(  2) = 0.345456507169149134898D-05
    w(  3) = 0.736624069102321668857D-05
    w(  4) = 0.125792781889592743525D-04
    w(  5) = 0.190213681905875816679D-04
    w(  6) = 0.266376412339000901358D-04
    w(  7) = 0.353751372055189588628D-04
    w(  8) = 0.451863674126296143105D-04
    w(  9) = 0.560319507856164252140D-04
    w( 10) = 0.678774554733972416227D-04
    w( 11) = 0.806899228014035293851D-04
    w( 12) = 0.944366322532705527066D-04
    w( 13) = 0.109085545645741522051D-03
    w( 14) = 0.124606200241498368482D-03
    w( 15) = 0.140970302204104791413D-03
    w( 16) = 0.158151830411132242924D-03
    w( 17) = 0.176126765545083195474D-03
    w( 18) = 0.194872642236641146532D-03
    w( 19) = 0.214368090034216937149D-03
    w( 20) = 0.234592462123925204879D-03
    w( 21) = 0.255525589595236862014D-03
    w( 22) = 0.277147657465187357459D-03
    w( 23) = 0.299439176850911730874D-03
    w( 24) = 0.322381020652862389664D-03
    w( 25) = 0.345954492129903871350D-03
    w( 26) = 0.370141402122251665232D-03
    w( 27) = 0.394924138246873704434D-03
    w( 28) = 0.420285716355361231823D-03
    w( 29) = 0.446209810101403247488D-03
    w( 30) = 0.472680758429262691232D-03
    w( 31) = 0.499683553312800484519D-03
    w( 32) = 0.527203811431658386125D-03
    w( 33) = 0.555227733977307579715D-03
    w( 34) = 0.583742058714979703847D-03
    w( 35) = 0.612734008012225209294D-03
    w( 36) = 0.642191235948505088403D-03
    w( 37) = 0.672101776960108194646D-03
    w( 38) = 0.702453997827572321358D-03
    w( 39) = 0.733236554224767912055D-03
    w( 40) = 0.764438352543882784191D-03
    w( 41) = 0.796048517297550871506D-03
    w( 42) = 0.828056364077226302608D-03
    w( 43) = 0.860451377808527848128D-03
    w( 44) = 0.893223195879324912340D-03
    w( 45) = 0.926361595613111283368D-03
    w( 46) = 0.959856485506936206261D-03
    w( 47) = 0.993697899638760857945D-03
    w( 48) = 0.102787599466367326179D-02
    w( 49) = 0.106238104885340071375D-02
    w( 50) = 0.109720346268191941940D-02
    w( 51) = 0.113233376051597664917D-02
    w( 52) = 0.116776259302858043685D-02
    w( 53) = 0.120348074001265964881D-02
    w( 54) = 0.123947911332878396534D-02
    w( 55) = 0.127574875977346947345D-02
    w( 56) = 0.131228086370221478128D-02
    w( 57) = 0.134906674928353113127D-02
    w( 58) = 0.138609788229672549700D-02
    w( 59) = 0.142336587141720519900D-02
    w( 60) = 0.146086246895890987689D-02
    w( 61) = 0.149857957106456636214D-02
    w( 62) = 0.153650921735128916170D-02
    w( 63) = 0.157464359003212166189D-02
    w( 64) = 0.161297501254393423070D-02
    w( 65) = 0.165149594771914570655D-02
    w( 66) = 0.169019899554346019117D-02
    w( 67) = 0.172907689054461607168D-02
    w( 68) = 0.176812249885838886701D-02
    w( 69) = 0.180732881501808930079D-02
    w( 70) = 0.184668895851282540913D-02
    w( 71) = 0.188619617015808475394D-02
    w( 72) = 0.192584380831993546204D-02
    w( 73) = 0.196562534503150547732D-02
    w( 74) = 0.200553436203751169944D-02
    w( 75) = 0.204556454679958293446D-02
    w( 76) = 0.208570968849203942640D-02
    w( 77) = 0.212596367401472533045D-02
    w( 78) = 0.216632048404649142727D-02
    w( 79) = 0.220677418916003329194D-02
    w( 80) = 0.224731894601603393082D-02
    w( 81) = 0.228794899365195972378D-02
    w( 82) = 0.232865864987842738864D-02
    w( 83) = 0.236944230779380495146D-02
    w( 84) = 0.241029443242563417382D-02
    w( 85) = 0.245120955750556483923D-02
    w( 86) = 0.249218228238276930060D-02
    w( 87) = 0.253320726907925325750D-02
    w( 88) = 0.257427923948908888092D-02
    w( 89) = 0.261539297272236109225D-02
    w( 90) = 0.265654330259352828314D-02
    w( 91) = 0.269772511525294586667D-02
    w( 92) = 0.273893334695947541201D-02
    w( 93) = 0.278016298199139435045D-02
    w( 94) = 0.282140905069222207923D-02
    w( 95) = 0.286266662764757868253D-02
    w( 96) = 0.290393082998878368175D-02
    w( 97) = 0.294519681581857582284D-02
    w( 98) = 0.298645978275408290247D-02
    w( 99) = 0.302771496658198544480D-02
    w(100) = 0.306895764002069252174D-02
    w(101) = 0.311018311158427546158D-02
    w(102) = 0.315138672454287935858D-02
    w(103) = 0.319256385597434736790D-02
    w(104) = 0.323370991590184336368D-02
    w(105) = 0.327482034651233969564D-02
    w(106) = 0.331589062145094394706D-02
    w(107) = 0.335691624518616761342D-02
    w(108) = 0.339789275244138669739D-02
    w(109) = 0.343881570768790591876D-02
    w(110) = 0.347968070469521146972D-02
    w(111) = 0.352048336613417922682D-02
    w(112) = 0.356121934322919357659D-02
    w(113) = 0.360188431545532431869D-02
    w(114) = 0.364247399027690353194D-02
    w(115) = 0.368298410292403911967D-02
    w(116) = 0.372341041620379550870D-02
    w(117) = 0.376374872034296338241D-02
    w(118) = 0.380399483285952829161D-02
    w(119) = 0.384414459846013158917D-02
    w(120) = 0.388419388896099560998D-02
    w(121) = 0.392413860322995774660D-02
    w(122) = 0.396397466714742455513D-02
    w(123) = 0.400369803358421688562D-02
    w(124) = 0.404330468239442998549D-02
    w(125) = 0.408279062042157838350D-02
    w(126) = 0.412215188151643401528D-02
    w(127) = 0.416138452656509745764D-02
    w(128) = 0.420048464352596631772D-02
    w(129) = 0.423944834747438184434D-02
    w(130) = 0.427827178065384480959D-02
    w(131) = 0.431695111253279479928D-02
    w(132) = 0.435548253986604343679D-02
    w(133) = 0.439386228676004195260D-02
    w(134) = 0.443208660474124713206D-02
    w(135) = 0.447015177282692726900D-02
    w(136) = 0.450805409759782158001D-02
    w(137) = 0.454578991327213285488D-02
    w(138) = 0.458335558178039420335D-02
    w(139) = 0.462074749284080687482D-02
    w(140) = 0.465796206403469754658D-02
    w(141) = 0.469499574088179046532D-02
    w(142) = 0.473184499691503264714D-02
    w(143) = 0.476850633375474925263D-02
    w(144) = 0.480497628118194150483D-02
    w(145) = 0.484125139721057135214D-02
    w(146) = 0.487732826815870573054D-02
    w(147) = 0.491320350871841897367D-02
    w(148) = 0.494887376202437487201D-02
    w(149) = 0.498433569972103029914D-02
    w(150) = 0.501958602202842039909D-02
    w(151) = 0.505462145780650125058D-02
    w(152) = 0.508943876461803986674D-02
    w(153) = 0.512403472879005351831D-02
    w(154) = 0.515840616547381084096D-02
    w(155) = 0.519254991870341614863D-02
    w(156) = 0.522646286145300596306D-02
    w(157) = 0.526014189569259311205D-02
    w(158) = 0.529358395244259896547D-02
    w(159) = 0.532678599182711857974D-02
    w(160) = 0.535974500312596681161D-02
    w(161) = 0.539245800482555593606D-02
    w(162) = 0.542492204466865704951D-02
    w(163) = 0.545713419970309863995D-02
    w(164) = 0.548909157632945623482D-02
    w(165) = 0.552079131034778706457D-02
    w(166) = 0.555223056700346326850D-02
    w(167) = 0.558340654103215637610D-02
    w(168) = 0.561431645670402467678D-02
    w(169) = 0.564495756786715368885D-02
    w(170) = 0.567532715799029830087D-02
    w(171) = 0.570542254020497332312D-02
    w(172) = 0.573524105734693719020D-02
    w(173) = 0.576478008199711142954D-02
    w(174) = 0.579403701652197628421D-02
    w(175) = 0.582300929311348057702D-02
    w(176) = 0.585169437382850155033D-02
    w(177) = 0.588008975062788803205D-02
    w(178) = 0.590819294541511788161D-02
    w(179) = 0.593600151007459827614D-02
    w(180) = 0.596351302650963502011D-02
    w(181) = 0.599072510668009471472D-02
    w(182) = 0.601763539263978131522D-02
    w(183) = 0.604424155657354634589D-02
    w(184) = 0.607054130083414983949D-02
    w(185) = 0.609653235797888692923D-02
    w(186) = 0.612221249080599294931D-02
    w(187) = 0.614757949239083790214D-02
    w(188) = 0.617263118612191922727D-02
    w(189) = 0.619736542573665996342D-02
    w(190) = 0.622178009535701763157D-02
    w(191) = 0.624587310952490748541D-02
    w(192) = 0.626964241323744217671D-02
    w(193) = 0.629308598198198836688D-02
    w(194) = 0.631620182177103938227D-02
    w(195) = 0.633898796917690165912D-02
    w(196) = 0.636144249136619145314D-02
    w(197) = 0.638356348613413709795D-02
    w(198) = 0.640534908193868098342D-02
    w(199) = 0.642679743793437438922D-02
    w(200) = 0.644790674400605734710D-02
    w(201) = 0.646867522080231481688D-02
    w(202) = 0.648910111976869964292D-02
    w(203) = 0.650918272318071200827D-02
    w(204) = 0.652891834417652442012D-02
    w(205) = 0.654830632678944064054D-02
    w(206) = 0.656734504598007641819D-02
    w(207) = 0.658603290766824937794D-02
    w(208) = 0.660436834876456498276D-02
    w(209) = 0.662234983720168509457D-02
    w(210) = 0.663997587196526532519D-02
    w(211) = 0.665724498312454708217D-02
    w(212) = 0.667415573186258997654D-02
    w(213) = 0.669070671050613006584D-02
    w(214) = 0.670689654255504925648D-02
    w(215) = 0.672272388271144108036D-02
    w(216) = 0.673818741690825799086D-02
    w(217) = 0.675328586233752529078D-02
    w(218) = 0.676801796747810680683D-02
    w(219) = 0.678238251212300746082D-02
    w(220) = 0.679637830740619795480D-02
    w(221) = 0.681000419582894688374D-02
    w(222) = 0.682325905128564571420D-02
    w(223) = 0.683614177908911221841D-02
    w(224) = 0.684865131599535812903D-02
    w(225) = 0.686078663022780697951D-02
    w(226) = 0.687254672150094831613D-02
    w(227) = 0.688393062104341470995D-02
    w(228) = 0.689493739162046825872D-02
    w(229) = 0.690556612755588354803D-02
    w(230) = 0.691581595475321433825D-02
    w(231) = 0.692568603071643155621D-02
    w(232) = 0.693517554456992049848D-02
    w(233) = 0.694428371707782549438D-02
    w(234) = 0.695300980066273063177D-02
    w(235) = 0.696135307942366551493D-02
    w(236) = 0.696931286915342540213D-02
    w(237) = 0.697688851735519545845D-02
    w(238) = 0.698407940325846925786D-02
    w(239) = 0.699088493783425207545D-02
    w(240) = 0.699730456380953992594D-02
    w(241) = 0.700333775568106572820D-02
    w(242) = 0.700898401972830440494D-02
    w(243) = 0.701424289402572916425D-02
    w(244) = 0.701911394845431165171D-02
    w(245) = 0.702359678471225911031D-02
    w(246) = 0.702769103632498213858D-02
    w(247) = 0.703139636865428709508D-02
    w(248) = 0.703471247890678765907D-02
    w(249) = 0.703763909614153052319D-02
    w(250) = 0.704017598127683066242D-02
    w(251) = 0.704232292709631209597D-02
    w(252) = 0.704407975825415053266D-02
    w(253) = 0.704544633127951476780D-02
    w(254) = 0.704642253458020417748D-02
    w(255) = 0.704700828844548013730D-02
    w(256) = 0.704720354504808967346D-02
    w(257) = 0.704700828844548013730D-02
    w(258) = 0.704642253458020417748D-02
    w(259) = 0.704544633127951476780D-02
    w(260) = 0.704407975825415053266D-02
    w(261) = 0.704232292709631209597D-02
    w(262) = 0.704017598127683066242D-02
    w(263) = 0.703763909614153052319D-02
    w(264) = 0.703471247890678765907D-02
    w(265) = 0.703139636865428709508D-02
    w(266) = 0.702769103632498213858D-02
    w(267) = 0.702359678471225911031D-02
    w(268) = 0.701911394845431165171D-02
    w(269) = 0.701424289402572916425D-02
    w(270) = 0.700898401972830440494D-02
    w(271) = 0.700333775568106572820D-02
    w(272) = 0.699730456380953992594D-02
    w(273) = 0.699088493783425207545D-02
    w(274) = 0.698407940325846925786D-02
    w(275) = 0.697688851735519545845D-02
    w(276) = 0.696931286915342540213D-02
    w(277) = 0.696135307942366551493D-02
    w(278) = 0.695300980066273063177D-02
    w(279) = 0.694428371707782549438D-02
    w(280) = 0.693517554456992049848D-02
    w(281) = 0.692568603071643155621D-02
    w(282) = 0.691581595475321433825D-02
    w(283) = 0.690556612755588354803D-02
    w(284) = 0.689493739162046825872D-02
    w(285) = 0.688393062104341470995D-02
    w(286) = 0.687254672150094831613D-02
    w(287) = 0.686078663022780697951D-02
    w(288) = 0.684865131599535812903D-02
    w(289) = 0.683614177908911221841D-02
    w(290) = 0.682325905128564571420D-02
    w(291) = 0.681000419582894688374D-02
    w(292) = 0.679637830740619795480D-02
    w(293) = 0.678238251212300746082D-02
    w(294) = 0.676801796747810680683D-02
    w(295) = 0.675328586233752529078D-02
    w(296) = 0.673818741690825799086D-02
    w(297) = 0.672272388271144108036D-02
    w(298) = 0.670689654255504925648D-02
    w(299) = 0.669070671050613006584D-02
    w(300) = 0.667415573186258997654D-02
    w(301) = 0.665724498312454708217D-02
    w(302) = 0.663997587196526532519D-02
    w(303) = 0.662234983720168509457D-02
    w(304) = 0.660436834876456498276D-02
    w(305) = 0.658603290766824937794D-02
    w(306) = 0.656734504598007641819D-02
    w(307) = 0.654830632678944064054D-02
    w(308) = 0.652891834417652442012D-02
    w(309) = 0.650918272318071200827D-02
    w(310) = 0.648910111976869964292D-02
    w(311) = 0.646867522080231481688D-02
    w(312) = 0.644790674400605734710D-02
    w(313) = 0.642679743793437438922D-02
    w(314) = 0.640534908193868098342D-02
    w(315) = 0.638356348613413709795D-02
    w(316) = 0.636144249136619145314D-02
    w(317) = 0.633898796917690165912D-02
    w(318) = 0.631620182177103938227D-02
    w(319) = 0.629308598198198836688D-02
    w(320) = 0.626964241323744217671D-02
    w(321) = 0.624587310952490748541D-02
    w(322) = 0.622178009535701763157D-02
    w(323) = 0.619736542573665996342D-02
    w(324) = 0.617263118612191922727D-02
    w(325) = 0.614757949239083790214D-02
    w(326) = 0.612221249080599294931D-02
    w(327) = 0.609653235797888692923D-02
    w(328) = 0.607054130083414983949D-02
    w(329) = 0.604424155657354634589D-02
    w(330) = 0.601763539263978131522D-02
    w(331) = 0.599072510668009471472D-02
    w(332) = 0.596351302650963502011D-02
    w(333) = 0.593600151007459827614D-02
    w(334) = 0.590819294541511788161D-02
    w(335) = 0.588008975062788803205D-02
    w(336) = 0.585169437382850155033D-02
    w(337) = 0.582300929311348057702D-02
    w(338) = 0.579403701652197628421D-02
    w(339) = 0.576478008199711142954D-02
    w(340) = 0.573524105734693719020D-02
    w(341) = 0.570542254020497332312D-02
    w(342) = 0.567532715799029830087D-02
    w(343) = 0.564495756786715368885D-02
    w(344) = 0.561431645670402467678D-02
    w(345) = 0.558340654103215637610D-02
    w(346) = 0.555223056700346326850D-02
    w(347) = 0.552079131034778706457D-02
    w(348) = 0.548909157632945623482D-02
    w(349) = 0.545713419970309863995D-02
    w(350) = 0.542492204466865704951D-02
    w(351) = 0.539245800482555593606D-02
    w(352) = 0.535974500312596681161D-02
    w(353) = 0.532678599182711857974D-02
    w(354) = 0.529358395244259896547D-02
    w(355) = 0.526014189569259311205D-02
    w(356) = 0.522646286145300596306D-02
    w(357) = 0.519254991870341614863D-02
    w(358) = 0.515840616547381084096D-02
    w(359) = 0.512403472879005351831D-02
    w(360) = 0.508943876461803986674D-02
    w(361) = 0.505462145780650125058D-02
    w(362) = 0.501958602202842039909D-02
    w(363) = 0.498433569972103029914D-02
    w(364) = 0.494887376202437487201D-02
    w(365) = 0.491320350871841897367D-02
    w(366) = 0.487732826815870573054D-02
    w(367) = 0.484125139721057135214D-02
    w(368) = 0.480497628118194150483D-02
    w(369) = 0.476850633375474925263D-02
    w(370) = 0.473184499691503264714D-02
    w(371) = 0.469499574088179046532D-02
    w(372) = 0.465796206403469754658D-02
    w(373) = 0.462074749284080687482D-02
    w(374) = 0.458335558178039420335D-02
    w(375) = 0.454578991327213285488D-02
    w(376) = 0.450805409759782158001D-02
    w(377) = 0.447015177282692726900D-02
    w(378) = 0.443208660474124713206D-02
    w(379) = 0.439386228676004195260D-02
    w(380) = 0.435548253986604343679D-02
    w(381) = 0.431695111253279479928D-02
    w(382) = 0.427827178065384480959D-02
    w(383) = 0.423944834747438184434D-02
    w(384) = 0.420048464352596631772D-02
    w(385) = 0.416138452656509745764D-02
    w(386) = 0.412215188151643401528D-02
    w(387) = 0.408279062042157838350D-02
    w(388) = 0.404330468239442998549D-02
    w(389) = 0.400369803358421688562D-02
    w(390) = 0.396397466714742455513D-02
    w(391) = 0.392413860322995774660D-02
    w(392) = 0.388419388896099560998D-02
    w(393) = 0.384414459846013158917D-02
    w(394) = 0.380399483285952829161D-02
    w(395) = 0.376374872034296338241D-02
    w(396) = 0.372341041620379550870D-02
    w(397) = 0.368298410292403911967D-02
    w(398) = 0.364247399027690353194D-02
    w(399) = 0.360188431545532431869D-02
    w(400) = 0.356121934322919357659D-02
    w(401) = 0.352048336613417922682D-02
    w(402) = 0.347968070469521146972D-02
    w(403) = 0.343881570768790591876D-02
    w(404) = 0.339789275244138669739D-02
    w(405) = 0.335691624518616761342D-02
    w(406) = 0.331589062145094394706D-02
    w(407) = 0.327482034651233969564D-02
    w(408) = 0.323370991590184336368D-02
    w(409) = 0.319256385597434736790D-02
    w(410) = 0.315138672454287935858D-02
    w(411) = 0.311018311158427546158D-02
    w(412) = 0.306895764002069252174D-02
    w(413) = 0.302771496658198544480D-02
    w(414) = 0.298645978275408290247D-02
    w(415) = 0.294519681581857582284D-02
    w(416) = 0.290393082998878368175D-02
    w(417) = 0.286266662764757868253D-02
    w(418) = 0.282140905069222207923D-02
    w(419) = 0.278016298199139435045D-02
    w(420) = 0.273893334695947541201D-02
    w(421) = 0.269772511525294586667D-02
    w(422) = 0.265654330259352828314D-02
    w(423) = 0.261539297272236109225D-02
    w(424) = 0.257427923948908888092D-02
    w(425) = 0.253320726907925325750D-02
    w(426) = 0.249218228238276930060D-02
    w(427) = 0.245120955750556483923D-02
    w(428) = 0.241029443242563417382D-02
    w(429) = 0.236944230779380495146D-02
    w(430) = 0.232865864987842738864D-02
    w(431) = 0.228794899365195972378D-02
    w(432) = 0.224731894601603393082D-02
    w(433) = 0.220677418916003329194D-02
    w(434) = 0.216632048404649142727D-02
    w(435) = 0.212596367401472533045D-02
    w(436) = 0.208570968849203942640D-02
    w(437) = 0.204556454679958293446D-02
    w(438) = 0.200553436203751169944D-02
    w(439) = 0.196562534503150547732D-02
    w(440) = 0.192584380831993546204D-02
    w(441) = 0.188619617015808475394D-02
    w(442) = 0.184668895851282540913D-02
    w(443) = 0.180732881501808930079D-02
    w(444) = 0.176812249885838886701D-02
    w(445) = 0.172907689054461607168D-02
    w(446) = 0.169019899554346019117D-02
    w(447) = 0.165149594771914570655D-02
    w(448) = 0.161297501254393423070D-02
    w(449) = 0.157464359003212166189D-02
    w(450) = 0.153650921735128916170D-02
    w(451) = 0.149857957106456636214D-02
    w(452) = 0.146086246895890987689D-02
    w(453) = 0.142336587141720519900D-02
    w(454) = 0.138609788229672549700D-02
    w(455) = 0.134906674928353113127D-02
    w(456) = 0.131228086370221478128D-02
    w(457) = 0.127574875977346947345D-02
    w(458) = 0.123947911332878396534D-02
    w(459) = 0.120348074001265964881D-02
    w(460) = 0.116776259302858043685D-02
    w(461) = 0.113233376051597664917D-02
    w(462) = 0.109720346268191941940D-02
    w(463) = 0.106238104885340071375D-02
    w(464) = 0.102787599466367326179D-02
    w(465) = 0.993697899638760857945D-03
    w(466) = 0.959856485506936206261D-03
    w(467) = 0.926361595613111283368D-03
    w(468) = 0.893223195879324912340D-03
    w(469) = 0.860451377808527848128D-03
    w(470) = 0.828056364077226302608D-03
    w(471) = 0.796048517297550871506D-03
    w(472) = 0.764438352543882784191D-03
    w(473) = 0.733236554224767912055D-03
    w(474) = 0.702453997827572321358D-03
    w(475) = 0.672101776960108194646D-03
    w(476) = 0.642191235948505088403D-03
    w(477) = 0.612734008012225209294D-03
    w(478) = 0.583742058714979703847D-03
    w(479) = 0.555227733977307579715D-03
    w(480) = 0.527203811431658386125D-03
    w(481) = 0.499683553312800484519D-03
    w(482) = 0.472680758429262691232D-03
    w(483) = 0.446209810101403247488D-03
    w(484) = 0.420285716355361231823D-03
    w(485) = 0.394924138246873704434D-03
    w(486) = 0.370141402122251665232D-03
    w(487) = 0.345954492129903871350D-03
    w(488) = 0.322381020652862389664D-03
    w(489) = 0.299439176850911730874D-03
    w(490) = 0.277147657465187357459D-03
    w(491) = 0.255525589595236862014D-03
    w(492) = 0.234592462123925204879D-03
    w(493) = 0.214368090034216937149D-03
    w(494) = 0.194872642236641146532D-03
    w(495) = 0.176126765545083195474D-03
    w(496) = 0.158151830411132242924D-03
    w(497) = 0.140970302204104791413D-03
    w(498) = 0.124606200241498368482D-03
    w(499) = 0.109085545645741522051D-03
    w(500) = 0.944366322532705527066D-04
    w(501) = 0.806899228014035293851D-04
    w(502) = 0.678774554733972416227D-04
    w(503) = 0.560319507856164252140D-04
    w(504) = 0.451863674126296143105D-04
    w(505) = 0.353751372055189588628D-04
    w(506) = 0.266376412339000901358D-04
    w(507) = 0.190213681905875816679D-04
    w(508) = 0.125792781889592743525D-04
    w(509) = 0.736624069102321668857D-05
    w(510) = 0.345456507169149134898D-05
    w(511) = 0.945715933950007048827D-06

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PATTERSON_SET - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of N.'
    write ( *, '(a)' ) '  N must be 1, 3, 7, 15, 31, 63, 127, 255 or 511.'
    stop

  end if

  return
end
function pythag ( a, b )

!*****************************************************************************80
!
!! PYTHAG computes SQRT ( A * A + B * B ) carefully.
!
!  Discussion:
!
!    The formula
!
!      pythag = sqrt ( a * a + b * b )
!
!    is reasonably accurate, but can fail if, for example, A*A is larger
!    than the machine overflow.  The formula can lose most of its accuracy
!    if the sum of the squares is very large or very small.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the two legs of a right triangle.
!
!    Output, real ( kind = 8 ) PYTHAG, the length of the hypotenuse.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) p
  real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u

  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0D+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

    do

      t = 4.0D+00 + r

      if ( t == 4.0D+00 ) then
        exit
      end if

      s = r / t
      u = 1.0D+00 + 2.0D+00 * s
      p = u * p
      r = ( s / u ) * ( s / u ) * r

    end do

  end if

  pythag = p

  return
end
function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EPSILON, the R8 round-off unit.
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) d_test
  real ( kind = 8 ) r8_epsilon

  d = 1.0D+00
  d_test = 1.0D+00 + d / 2.0D+00

  do while ( 1.0D+00 < d_test )
    d = d / 2.0D+00
    d_test = 1.0D+00 + d / 2.0D+00
  end do

  r8_epsilon = d

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( n ) = product ( 1 <= i <= n ) i
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function of N.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    N!!
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial
!    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL2, the value of the double factorial.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial2
  real ( kind = 8 ) r8_n

  if ( n < 1 ) then
    r8_factorial2 = 1.0D+00
    return
  end if

  r8_n = real ( n, kind = 8 )
  r8_factorial2 = 1.0D+00

  do while ( 1.0D+00 < r8_n )
    r8_factorial2 = r8_factorial2 * r8_n
    r8_n = r8_n - 2.0D+00
  end do

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
!
!  Coefficients for minimax approximation over (12, INF).
!
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
subroutine r8_hyper_2f1 ( a_input, b_input, c_input, x_input, hf )

!*****************************************************************************80
!
!! R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
!
!  Discussion:
!
!    A minor bug was corrected.  The HW variable, used in several places as
!    the "old" value of a quantity being iteratively improved, was not
!    being initialized.  JVB, 11 February 2008.
!
!    The original version of this program allowed the input arguments to
!    be modified, although they were restored to their input values before exit.
!    This is unacceptable if the input arguments are allowed to be constants.
!    The code has been modified so that the input arguments are never modified.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2008
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program provided that the copyright
!    is acknowledged.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_INPUT, B_INPUT, C_INPUT, X_INPUT,
!    the arguments of the function.  The user is allowed to pass these
!    values as constants or variables.
!    C_INPUT must not be equal to a nonpositive integer.
!    X_INPUT < 1.
!
!    Output, real ( kind = 8 ) HF, the value of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a_input
  real ( kind = 8 ) a0
  real ( kind = 8 ) aa
  real ( kind = 8 ) b
  real ( kind = 8 ) b_input
  real ( kind = 8 ) bb
  real ( kind = 8 ) c
  real ( kind = 8 ) c_input
  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ), parameter :: el = 0.5772156649015329D+00
  real ( kind = 8 ) eps
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) g0
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) g3
  real ( kind = 8 ) ga
  real ( kind = 8 ) gabc
  real ( kind = 8 ) gam
  real ( kind = 8 ) gb
  real ( kind = 8 ) gbm
  real ( kind = 8 ) gc
  real ( kind = 8 ) gca
  real ( kind = 8 ) gcab
  real ( kind = 8 ) gcb
  real ( kind = 8 ) gm
  real ( kind = 8 ) hf
  real ( kind = 8 ) hw
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical l0
  logical l1
  logical l2
  logical l3
  logical l4
  logical l5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nm
  real ( kind = 8 ) pa
  real ( kind = 8 ) pb
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r0
  real ( kind = 8 ) r1
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) rm
  real ( kind = 8 ) rp
  real ( kind = 8 ) sm
  real ( kind = 8 ) sp
  real ( kind = 8 ) sp0
  real ( kind = 8 ) x
  real ( kind = 8 ) x_input
  real ( kind = 8 ) x1
!
!  Immediately copy the input arguments!
!
  a = a_input
  b = b_input
  c = c_input
  x = x_input

  l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
  l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
  l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
  l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
  l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
  l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )

  if ( l0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  Integral C < 0.'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    stop
  end if

  if ( l1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    write ( *, '(a)' ) '  1 - X < 0, C - A - B < 0.'
    stop
  end if

  if ( 0.95D+00 < x ) then
    eps = 1.0D-08
  else
    eps = 1.0D-15
  end if

  if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then

    hf = 1.0D+00
    return

  else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then

    gc = r8_gamma ( c )
    gcab = r8_gamma ( c - a - b )
    gca = r8_gamma ( c - a )
    gcb = r8_gamma ( c - b )
    hf = gc * gcab / ( gca * gcb )
    return

  else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then

    g0 = sqrt ( pi ) * 2.0D+00**( - a )
    g1 = r8_gamma ( c )
    g2 = r8_gamma ( 1.0D+00 + a / 2.0D+00 - b )
    g3 = r8_gamma ( 0.5D+00 + 0.5D+00 * a )
    hf = g0 * g1 / ( g2 * g3 )
    return

  else if ( l2 .or. l3 ) then

    if ( l2 ) then
      nm = int ( abs ( a ) )
    end if

    if ( l3 ) then
      nm = int ( abs ( b ) )
    end if

    hf = 1.0D+00
    r = 1.0D+00

    do k = 1, nm
      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do

    return

  else if ( l4 .or. l5 ) then

    if ( l4 ) then
      nm = int ( abs ( c - a ) )
    end if

    if ( l5 ) then
      nm = int ( abs ( c - b ) )
    end if

    hf = 1.0D+00
    r  = 1.0D+00
    do k = 1, nm
      r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do
    hf = ( 1.0D+00 - x )**( c - a - b ) * hf
    return

  end if

  aa = a
  bb = b
  x1 = x

  if ( x < 0.0D+00 ) then
    x = x / ( x - 1.0D+00 )
    if ( a < c .and. b < a .and. 0.0D+00 < b ) then
      a = bb
      b = aa
    end if
    b = c - b
  end if

  if ( 0.75D+00 <= x ) then

    gm = 0.0D+00

    if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then

      m = int ( c - a - b )
      ga = r8_gamma ( a )
      gb = r8_gamma ( b )
      gc = r8_gamma ( c )
      gam = r8_gamma ( a + m )
      gbm = r8_gamma ( b + m )

      pa = r8_psi ( a )
      pb = r8_psi ( b )

      if ( m /= 0 ) then
        gm = 1.0D+00
      end if

      do j = 1, abs ( m ) - 1
        gm = gm * j
      end do

      rm = 1.0D+00
      do j = 1, abs ( m )
        rm = rm * j
      end do

      f0 = 1.0D+00
      r0 = 1.0D+00
      r1 = 1.0D+00
      sp0 = 0.0D+00
      sp = 0.0D+00

      if ( 0 <= m ) then

        c0 = gm * gc / ( gam * gbm )
        c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

        do k = 1, m - 1
          r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
            + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / real ( k, kind = 8 )
        end do

        f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + ( 1.0D+00 - a ) &
              / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
              + 1.0D+00 / ( b + j + k - 1.0D+00 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      else if ( m < 0 ) then

        m = - m
        c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
        c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

        do k = 1, m - 1
          r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / real ( k, kind = 8 )
        end do

        f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) &
            / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + 1.0D+00 / real ( j + k, kind = 8 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      end if

    else

      ga = r8_gamma ( a )
      gb = r8_gamma ( b )
      gc = r8_gamma ( c )
      gca = r8_gamma ( c - a )
      gcb = r8_gamma ( c - b )
      gcab = r8_gamma ( c - a - b )
      gabc = r8_gamma ( a + b - c )
      c0 = gc * gcab / ( gca * gcb )
      c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
      hf = 0.0D+00
      hw = hf
      r0 = c0
      r1 = c1

      do k = 1, 250

        r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
          / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

        r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
          / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

        hf = hf + r0 + r1

        if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
          exit
        end if

        hw = hf

      end do

      hf = hf + c0 + c1

    end if

  else

    a0 = 1.0D+00

    if ( a < c .and. c < 2.0D+00 * a .and. b < c .and. c < 2.0D+00 * b ) then

      a0 = ( 1.0D+00 - x )**( c - a - b )
      a = c - a
      b = c - b

    end if

    hf = 1.0D+00
    hw = hf
    r = 1.0D+00

    do k = 1, 250

      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x

      hf = hf + r

      if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
        exit
      end if

      hw = hf

    end do

    hf = a0 * hf

  end if

  if ( x1 < 0.0D+00 ) then
    x = x1
    c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
    hf = c0 * hf
  end if

  a = aa
  b = bb

  if ( 120 < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Warning!'
    write ( *, '(a)' ) '  A large number of iterations were needed.'
    write ( *, '(a)' ) '  The accuracy of the results should be checked.'
  end if

  return
end
function r8_psi ( xx )

!*****************************************************************************80
!
!! R8_PSI evaluates the function Psi(X).
!
!  Discussion:
!
!    This routine evaluates the logarithmic derivative of the
!    Gamma function,
!
!      PSI(X) = d/dx ( GAMMA(X) ) / GAMMA(X)
!             = d/dx LN ( GAMMA(X) )
!
!    for real X, where either
!
!      - XMAX1 < X < - XMIN, and X is not a negative integer,
!
!    or
!
!      XMIN < X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Anthony Strecok, Henry Thacher,
!    Chebyshev Approximations for the Psi Function,
!    Mathematics of Computation,
!    Volume 27, Number 121, January 1973, pages 123-127.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XX, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_PSI, the value of the function.
!
  implicit none

  real ( kind = 8 ) aug
  real ( kind = 8 ) den
  real ( kind = 8 ), parameter :: four = 4.0D+00
  real ( kind = 8 ), parameter :: fourth = 0.25D+00
  real ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nq
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ), dimension ( 9 ) :: p1 = (/ &
   4.5104681245762934160D-03, &
   5.4932855833000385356D+00, &
   3.7646693175929276856D+02, &
   7.9525490849151998065D+03, &
   7.1451595818951933210D+04, &
   3.0655976301987365674D+05, &
   6.3606997788964458797D+05, &
   5.8041312783537569993D+05, &
   1.6585695029761022321D+05 /)
  real ( kind = 8 ), dimension ( 7 ) :: p2 = (/ &
  -2.7103228277757834192D+00, &
  -1.5166271776896121383D+01, &
  -1.9784554148719218667D+01, &
  -8.8100958828312219821D+00, &
  -1.4479614616899842986D+00, &
  -7.3689600332394549911D-02, &
  -6.5135387732718171306D-21 /)
  real ( kind = 8 ), parameter :: piov4 = 0.78539816339744830962D+00
  real ( kind = 8 ), dimension ( 8 ) :: q1 = (/ &
   9.6141654774222358525D+01, &
   2.6287715790581193330D+03, &
   2.9862497022250277920D+04, &
   1.6206566091533671639D+05, &
   4.3487880712768329037D+05, &
   5.4256384537269993733D+05, &
   2.4242185002017985252D+05, &
   6.4155223783576225996D-08 /)
  real ( kind = 8 ), dimension ( 6 ) :: q2 = (/ &
   4.4992760373789365846D+01, &
   2.0240955312679931159D+02, &
   2.4736979003315290057D+02, &
   1.0742543875702278326D+02, &
   1.7463965060678569906D+01, &
   8.8427520398873480342D-01 /)
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) sgn
  real ( kind = 8 ), parameter :: three = 3.0D+00
  real ( kind = 8 ) upper
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: x01 = 187.0D+00
  real ( kind = 8 ), parameter :: x01d = 128.0D+00
  real ( kind = 8 ), parameter :: x02 = 6.9464496836234126266D-04
  real ( kind = 8 ), parameter :: xinf = 1.70D+38
  real ( kind = 8 ), parameter :: xlarge = 2.04D+15
  real ( kind = 8 ), parameter :: xmax1 = 3.60D+16
  real ( kind = 8 ), parameter :: xmin1 = 5.89D-39
  real ( kind = 8 ), parameter :: xsmall = 2.05D-09
  real ( kind = 8 ) xx
  real ( kind = 8 ) z
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  x = xx
  w = abs ( x )
  aug = zero
!
!  Check for valid arguments, then branch to appropriate algorithm.
!
  if ( xmax1 <= - x .or. w < xmin1 ) then

    if ( zero < x ) then
      r8_psi = - xinf
    else
      r8_psi = xinf
    end if

    return
  end if

  if ( x < half ) then
!
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!
    if ( w <= xsmall ) then

      aug = - one / x
!
!  Argument reduction for cotangent.
!
    else

      if ( x < zero ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - real ( int ( w ), kind = 8 )
      nq = int ( w * four )
      w = four * ( w - real ( nq, kind = 8 ) * fourth )
!
!  W is now related to the fractional part of 4.0 * X.
!  Adjust argument to correspond to values in the first
!  quadrant and determine the sign.
!
      n = nq / 2

      if ( n + n /= nq ) then
        w = one - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) /= 0 ) then
        sgn = - sgn
      end if
!
!  Determine the final value for  -pi * cotan(pi*x).
!
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) == 0 ) then
!
!  Check for singularity.
!
        if ( z == zero ) then

          if ( zero < x ) then
            r8_psi = -xinf
          else
            r8_psi = xinf
          end if

          return
        end if

        aug = sgn * ( four / tan ( z ) )

      else

        aug = sgn * ( four * tan ( z ) )

      end if

    end if

    x = one - x

  end if
!
!  0.5 <= X <= 3.0.
!
  if ( x <= three ) then

    den = x
    upper = p1(1) * x
    do i = 1, 7
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do
    den = ( upper + p1(9) ) / ( den + q1(8) )
    x = ( x - x01 / x01d ) - x02
    r8_psi = den * x + aug
    return

  end if
!
!  3.0 < X.
!
  if ( x < xlarge ) then
    w = one / ( x * x )
    den = w
    upper = p2(1) * w
    do i = 1, 5
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do
    aug = ( upper + p2(7) ) / ( den + q2(6) ) - half / x + aug
  end if

  r8_psi = aug + log ( x )

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_reverse ( n, a )

!*****************************************************************************80
!
!! R8VEC_REVERSE reverses the elements of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of double precision values.
!
!    Input:
!
!      N = 5,
!      A = ( 11.0, 12.0, 13.0, 14.0, 15.0 )
!
!    Output:
!
!      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t

  do i = 1, n / 2
    t        = a(i)
    a(i)     = a(n+1-i)
    a(n+1-i) = t
  end do

  return
end
subroutine radau_compute ( n, x, w )

!*****************************************************************************80
!
!! RADAU_COMPUTE computes a Radau quadrature rule.
!
!  Discussion:
!
!    The Radau rule is distinguished by the fact that the left endpoint
!    (-1) is always an abscissa.
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*N-2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2007
!
!  Author:
!
!    Original MATLAB version by Greg von Winckel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
!    Spectral Methods in Fluid Dynamics,
!    Springer, 1993,
!    ISNB13: 978-3540522058,
!    LC: QA377.S676.
!
!    Francis Hildebrand,
!    Section 8.11,
!    Introduction to Numerical Analysis,
!    Dover, 1987,
!    ISBN13: 978-0486653631,
!    LC: QA300.H5.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be at least 1.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(n,n+1)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) tolerance
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xold(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RADAU_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) ' N must be at least 1.'
    stop
  end if

  tolerance = 100.0D+00 * epsilon ( tolerance )
!
!  Initial estimate for the abscissas is the Chebyshev-Gauss-Radau nodes.
!
  do i = 1, n
    x(i) = - cos ( 2.0D+00 * pi * real (     i - 1, kind = 8 ) &
                                / real ( 2 * n - 1, kind = 8 ) )
  end do

  xold(1:n) = 2.0D+00

  do while ( tolerance < maxval ( abs ( x(1:n) - xold(1:n) ) ) )

    xold(1:n) = x(1:n)

    do j = 1, n+1
      p(1,j) = ( -1.0D+00 ) **( j - 1 )
    end do

    p(2:n,1) = 1.0D+00
    p(2:n,2) = x(2:n)

    do j = 2, n
      p(2:n,j+1) = ( real ( 2 * j - 1, kind = 8 ) * x(2:n) * p(2:n,j)     &
                   + real (   - j + 1, kind = 8 ) *          p(2:n,j-1) ) &
                   / real (     j,     kind = 8 )
    end do

    x(2:n) = xold(2:n) - ( ( 1.0D+00 - xold(2:n) ) / real ( n, kind = 8 ) ) &
      * ( p(2:n,n) + p(2:n,n+1) ) / ( p(2:n,n) - p(2:n,n+1) )

  end do

  w(1) = 2.0D+00 / real ( n * n, kind = 8 )
  w(2:n) = ( 1.0D+00 - x(2:n) ) / ( real ( n, kind = 8 ) * p(2:n,n) )**2

  return
end
subroutine radau_set ( n, x, w )

!*****************************************************************************80
!
!! RADAU_SET sets abscissas and weights for Radau quadrature.
!
!  Discussion:
!
!    The Radau rule is distinguished by the fact that the left endpoint
!    (-1) is always an abscissa.
!
!    The integral:
!
!      Integral ( -1 <= X <= 1 ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*N-2).
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
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be between 1 and 15.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =   - 1.0D+00
    w(1) =   2.0D+00

  else if ( n == 2 ) then

    x(1) =  - 1.0D+00
    x(2) =    1.0D+00 / 3.0D+00

    w(1) =  0.5D+00
    w(2) =  1.5D+00

  else if ( n == 3 ) then

    x(1) =   - 1.0D+00
    x(2) =   - 0.289897948556635619639456814941D+00
    x(3) =     0.689897948556635619639456814941D+00

    w(1) =  0.222222222222222222222222222222D+00
    w(2) =  0.102497165237684322767762689304D+01
    w(3) =  0.752806125400934550100150884739D+00

  else if ( n == 4 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.575318923521694112050483779752D+00
    x(3) =    0.181066271118530578270147495862D+00
    x(4) =    0.822824080974592105208907712461D+00

    w(1) =  0.125D+00
    w(2) =  0.657688639960119487888578442146D+00
    w(3) =  0.776386937686343761560464613780D+00
    w(4) =  0.440924422353536750550956944074D+00

  else if ( n == 5 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.720480271312438895695825837750D+00
    x(3) =  - 0.167180864737833640113395337326D+00
    x(4) =    0.446313972723752344639908004629D+00
    x(5) =    0.885791607770964635613757614892D+00

    w(1) =  0.08D+00
    w(2) =  0.446207802167141488805120436457D+00
    w(3) =  0.623653045951482508163709823153D+00
    w(4) =  0.562712030298924120384345300681D+00
    w(5) =  0.287427121582451882646824439708D+00

  else if ( n == 6 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.802929828402347147753002204224D+00
    x(3) =  - 0.390928546707272189029229647442D+00
    x(4) =    0.124050379505227711989974959990D+00
    x(5) =    0.603973164252783654928415726409D+00
    x(6) =    0.920380285897062515318386619813D+00

    w(1) =  0.555555555555555555555555555556D-01
    w(2) =  0.319640753220510966545779983796D+00
    w(3) =  0.485387188468969916159827915587D+00
    w(4) =  0.520926783189574982570229406570D+00
    w(5) =  0.416901334311907738959406382743D+00
    w(6) =  0.201588385253480840209200755749D+00

  else if ( n == 7 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.853891342639482229703747931639D+00
    x(3) =  - 0.538467724060109001833766720231D+00
    x(4) =  - 0.117343037543100264162786683611D+00
    x(5) =    0.326030619437691401805894055838D+00
    x(6) =    0.703842800663031416300046295008D+00
    x(7) =    0.941367145680430216055899446174D+00

    w(1) =  0.408163265306122448979591836735D-01
    w(2) =  0.239227489225312405787077480770D+00
    w(3) =  0.380949873644231153805938347876D+00
    w(4) =  0.447109829014566469499348953642D+00
    w(5) =  0.424703779005955608398308039150D+00
    w(6) =  0.318204231467301481744870434470D+00
    w(7) =  0.148988471112020635866497560418D+00

  else if ( n == 8 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.887474878926155707068695617935D+00
    x(3) =  - 0.639518616526215270024840114382D+00
    x(4) =  - 0.294750565773660725252184459658D+00
    x(5) =    0.943072526611107660028971153047D-01
    x(6) =    0.468420354430821063046421216613D+00
    x(7) =    0.770641893678191536180719525865D+00
    x(8) =    0.955041227122575003782349000858D+00

    w(1) =  0.03125D+00
    w(2) =  0.185358154802979278540728972699D+00
    w(3) =  0.304130620646785128975743291400D+00
    w(4) =  0.376517545389118556572129261442D+00
    w(5) =  0.391572167452493593082499534004D+00
    w(6) =  0.347014795634501280228675918422D+00
    w(7) =  0.249647901329864963257869293513D+00
    w(8) =  0.114508814744257199342353728520D+00

  else if ( n == 9 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.910732089420060298533757956283D+00
    x(3) =  - 0.711267485915708857029562959544D+00
    x(4) =  - 0.426350485711138962102627520502D+00
    x(5) =  - 0.903733696068532980645444599064D-01
    x(6) =    0.256135670833455395138292079035D+00
    x(7) =    0.571383041208738483284917464837D+00
    x(8) =    0.817352784200412087992517083851D+00
    x(9) =    0.964440169705273096373589797925D+00

    w(1) =  0.246913580246913580246913580247D-01
    w(2) =  0.147654019046315385819588499802D+00
    w(3) =  0.247189378204593052361239794969D+00
    w(4) =  0.316843775670437978338000849642D+00
    w(5) =  0.348273002772966594071991031186D+00
    w(6) =  0.337693966975929585803724239792D+00
    w(7) =  0.286386696357231171146705637752D+00
    w(8) =  0.200553298024551957421165090417D+00
    w(9) =  0.907145049232829170128934984159D-01

  else if ( n == 10 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.927484374233581078117671398464D+00
    x(3) =  - 0.763842042420002599615429776011D+00
    x(4) =  - 0.525646030370079229365386614293D+00
    x(5) =  - 0.236234469390588049278459503207D+00
    x(6) =    0.760591978379781302337137826389D-01
    x(7) =    0.380664840144724365880759065541D+00
    x(8) =    0.647766687674009436273648507855D+00
    x(9) =    0.851225220581607910728163628088D+00
    x(10) =   0.971175180702246902734346518378D+00

    w(1) =  0.02D+00
    w(2) =  0.120296670557481631517310522702D+00
    w(3) =  0.204270131879000675555788672223D+00
    w(4) =  0.268194837841178696058554475262D+00
    w(5) =  0.305859287724422621016275475401D+00
    w(6) =  0.313582457226938376695902847302D+00
    w(7) =  0.290610164832918311146863077963D+00
    w(8) =  0.239193431714379713376571966160D+00
    w(9) =  0.164376012736921475701681668908D+00
    w(10) = 0.736170054867584989310512940790D-01

  else if ( n == 11 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.939941935677027005913871284731D+00
    x(3) =  - 0.803421975580293540697597956820D+00
    x(4) =  - 0.601957842073797690275892603234D+00
    x(5) =  - 0.351888923353330214714301017870D+00
    x(6) =  - 0.734775314313212657461903554238D-01
    x(7) =    0.210720306228426314076095789845D+00
    x(8) =    0.477680647983087519467896683890D+00
    x(9) =    0.705777100713859519144801128840D+00
    x(10) =   0.876535856245703748954741265611D+00
    x(11) =   0.976164773135168806180508826082D+00

    w(1) =  0.165289256198347107438016528926D-01
    w(2) =  0.998460819079680638957534695802D-01
    w(3) =  0.171317619206659836486712649042D+00
    w(4) =  0.228866123848976624401683231126D+00
    w(5) =  0.267867086189684177806638163355D+00
    w(6) =  0.285165563941007337460004408915D+00
    w(7) =  0.279361333103383045188962195720D+00
    w(8) =  0.250925377697128394649140267633D+00
    w(9) =  0.202163108540024418349931754266D+00
    w(10) = 0.137033682133202256310153880580D+00
    w(11) = 0.609250978121311347072183268883D-01

  else if ( n == 12 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.949452759204959300493337627077D+00
    x(3) =  - 0.833916773105189706586269254036D+00
    x(4) =  - 0.661649799245637148061133087811D+00
    x(5) =  - 0.444406569781935851126642615609D+00
    x(6) =  - 0.196994559534278366455441427346D+00
    x(7) =    0.637247738208319158337792384845D-01
    x(8) =    0.319983684170669623532789532206D+00
    x(9) =    0.554318785912324288984337093085D+00
    x(10) =   0.750761549711113852529400825472D+00
    x(11) =   0.895929097745638894832914608454D+00
    x(12) =   0.979963439076639188313950540264D+00

    w(1) =  0.138888888888888888888888888888D-01
    w(2) =  0.841721349386809762415796536813D-01
    w(3) =  0.145563668853995128522547654706D+00
    w(4) =  0.196998534826089634656049637969D+00
    w(5) =  0.235003115144985839348633985940D+00
    w(6) =  0.256991338152707776127974253598D+00
    w(7) =  0.261465660552133103438074715743D+00
    w(8) =  0.248121560804009959403073107079D+00
    w(9) =  0.217868879026192438848747482023D+00
    w(10) = 0.172770639313308564306065766966D+00
    w(11) = 0.115907480291738392750341908272D+00
    w(12) = 0.512480992072692974680229451351D-01

  else if ( n == 13 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.956875873668299278183813833834D+00
    x(3) =  - 0.857884202528822035697620310269D+00
    x(4) =  - 0.709105087529871761580423832811D+00
    x(5) =  - 0.519197779050454107485205148087D+00
    x(6) =  - 0.299201300554509985532583446686D+00
    x(7) =  - 0.619016986256353412578604857936D-01
    x(8) =    0.178909837597084635021931298881D+00
    x(9) =    0.409238231474839556754166331248D+00
    x(10) =   0.615697890940291918017885487543D+00
    x(11) =   0.786291018233046684731786459135D+00
    x(12) =   0.911107073689184553949066402429D+00
    x(13) =   0.982921890023145161262671078244D+00

    w(1) =  0.118343195266272189349112426036D-01
    w(2) =  0.719024162924955289397537405641D-01
    w(3) =  0.125103834331152358133769287976D+00
    w(4) =  0.171003460470616642463758674512D+00
    w(5) =  0.206960611455877074631132560829D+00
    w(6) =  0.230888862886995434012203758668D+00
    w(7) =  0.241398342287691148630866924129D+00
    w(8) =  0.237878547660712031342685189180D+00
    w(9) =  0.220534229288451464691077164199D+00
    w(10) = 0.190373715559631732254759820746D+00
    w(11) = 0.149150950090000205151491864242D+00
    w(12) = 0.992678068818470859847363877478D-01
    w(13) = 0.437029032679020748288533846051D-01

  else if ( n == 14 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.962779269978024297120561244319D+00
    x(3) =  - 0.877048918201462024795266773531D+00
    x(4) =  - 0.747389642613378838735429134263D+00
    x(5) =  - 0.580314056546874971105726664999D+00
    x(6) =  - 0.384202003439203313794083903375D+00
    x(7) =  - 0.168887928042680911008441695622D+00
    x(8) =    0.548312279917645496498107146428D-01
    x(9) =    0.275737205435522399182637403545D+00
    x(10) =   0.482752918588474966820418534355D+00
    x(11) =   0.665497977216884537008955042481D+00
    x(12) =   0.814809550601994729434217249123D+00
    x(13) =   0.923203722520643299246334950272D+00
    x(14) =   0.985270697947821356698617003172D+00

    w(1) =  0.102040816326530612244897959184D-01
    w(2) =  0.621220169077714601661329164668D-01
    w(3) =  0.108607722744362826826720935229D+00
    w(4) =  0.149620539353121355950520836946D+00
    w(5) =  0.183127002125729654123867302103D+00
    w(6) =  0.207449763335175672668082886489D+00
    w(7) =  0.221369811499570948931671683021D+00
    w(8) =  0.224189348002707794238414632220D+00
    w(9) =  0.215767100604618851381187446115D+00
    w(10) = 0.196525518452982430324613091930D+00
    w(11) = 0.167429727891086278990102277038D+00
    w(12) = 0.129939668737342347807425737146D+00
    w(13) = 0.859405354429804030893077310866D-01
    w(14) = 0.377071632698969142774627282919D-01

  else if ( n == 15 ) then

    x(1) =  - 1.0D+00
    x(2) =  - 0.967550468197200476562456018282D+00
    x(3) =  - 0.892605400120550767066811886849D+00
    x(4) =  - 0.778685617639031079381743321893D+00
    x(5) =  - 0.630779478886949283946148437224D+00
    x(6) =  - 0.455352905778529370872053455981D+00
    x(7) =  - 0.260073376740807915768961188263D+00
    x(8) =  - 0.534757226797460641074538896258D-01
    x(9) =    0.155410685384859484319182024964D+00
    x(10) =   0.357456512022127651195319205174D+00
    x(11) =   0.543831458701484016930711802760D+00
    x(12) =   0.706390264637572540152679669478D+00
    x(13) =   0.838029000636089631215097384520D+00
    x(14) =   0.932997190935973719928072142859D+00
    x(15) =   0.987166478414363086378359071811D+00

    w(1) =  0.888888888888888888888888888889D-02
    w(2) =  0.542027800486444943382142368018D-01
    w(3) =  0.951295994604808992038477266346D-01
    w(4) =  0.131875462504951632186262157944D+00
    w(5) =  0.162854477303832629448732245828D+00
    w(6) =  0.186715145839450908083795103799D+00
    w(7) =  0.202415187030618429872703310435D+00
    w(8) =  0.209268608147694581430889790306D+00
    w(9) =  0.206975960249553755479027321787D+00
    w(10) = 0.195637503045116116473556617575D+00
    w(11) = 0.175748872642447685670310440476D+00
    w(12) = 0.148179527003467253924682058743D+00
    w(13) = 0.114135203489752753013075582569D+00
    w(14) = 0.751083927605064397329716653914D-01
    w(15) = 0.328643915845935322530428528231D-01

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RADAU_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) '  Legal values are 1 to 15.'
    stop

  end if

  return
end
subroutine rule_adjust ( a, b, c, d, n, x, w )

!*****************************************************************************80
!
!! RULE_ADJUST maps a quadrature rule from [A,B] to [C,D].
!
!  Discussion:
!
!    Most quadrature rules are defined on a special interval, like
!    [-1,1] or [0,1].  To integrate over an interval, the abscissas
!    and weights must be adjusted.  This can be done on the fly,
!    or by calling this routine.
!
!    If the weight function W(X) is not 1, then the weight vector W will
!    require further adjustment by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the definition interval.
!
!    Input, real ( kind = 8 ) C, D, the endpoints of the integration interval.
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input/output, real ( kind = 8 ) X(N), the abscissas.
!
!    Input/output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  x(1:n) = ( ( b - x(1:n)     ) * c   &
           + (     x(1:n) - a ) * d ) &
           / ( b          - a )

  w(1:n) = ( ( d - c ) / ( b - a ) ) * w(1:n)

  return
end
subroutine sum_sub ( func, a, b, nsub, n, xlo, xhi, x, w, result )

!*****************************************************************************80
!
!! SUM_SUB carries out a composite quadrature rule.
!
!  Discussion:
!
!    SUM_SUB assumes the original rule was written for [XLO,XHI].
!
!    The integral:
!
!      Integral ( A <= X <= B ) F(X) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of integration.
!
!    Input, integer ( kind = 4 ) NSUB, the number of equal subintervals into
!    which the finite interval (A,B) is to be subdivided for
!    higher accuracy.  NSUB must be at least 1.
!
!    Input, integer ( kind = 4 ) N, the order.
!    N must be at least 1.
!
!    Input, real ( kind = 8 ) XLO, XHI, the left and right endpoints of
!    the interval over which the quadrature rule was defined.
!
!    Input, real ( kind = 8 ) X(N), the abscissas of a quadrature
!    rule for the interval [XLO,XHI].
!
!    Input, real ( kind = 8 ) W(N), the weights.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) a_sub
  real ( kind = 8 ) b
  real ( kind = 8 ) b_sub
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) quad_sub
  real ( kind = 8 ) result
  real ( kind = 8 ) result_sub
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_sub
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xval

  if ( a == b ) then
    result = 0.0D+00
    return
  end if

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUM_SUB - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive value of N = ', n
    stop
  end if

  if ( nsub < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUM_SUB - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive value of NSUB = ', nsub
    stop
  end if

  if ( xlo == xhi ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUM_SUB - Fatal error!'
    write ( *, '(a)' ) '  XLO = XHI.'
    stop
  end if

  volume = 0.0D+00
  result = 0.0D+00

  do j = 1, nsub

    a_sub = ( real ( nsub - j + 1, kind = 8 ) * a   &
            + real (        j - 1, kind = 8 ) * b ) &
            / real ( nsub,         kind = 8 )

    b_sub = ( real ( nsub - j, kind = 8 )     * a   &
            + real (        j, kind = 8 )     * b ) &
            / real ( nsub,     kind = 8 )

    quad_sub = 0.0D+00
    do i = 1, n
      xval = ( ( xhi - x(i)       ) * a_sub   &
             + (       x(i) - xlo ) * b_sub ) &
             / ( xhi        - xlo )
      quad_sub = quad_sub + w(i) * func ( xval )
    end do

    volume_sub = ( b - a ) / ( ( xhi - xlo ) * real ( nsub, kind = 8 ) )
    result_sub = quad_sub * volume_sub

    volume = volume + volume_sub
    result = result + result_sub

  end do

  return
end
subroutine sum_sub_gk ( func, a, b, nsub, ng, wg, resultg, nk, &
  xk, wk, resultk, error )

!*****************************************************************************80
!
!! SUM_SUB_GK carries out a composite Gauss-Kronrod rule.
!
!  Discussion:
!
!    The integral:
!
!      Integral ( A <= X <= B ) F(X) dx
!
!    The quadrature rule:
!
!      H = ( B - A ) / NSUB
!      XMID(J) = A + 0.5 * H * ( 2 * J - 1 )
!
!      Sum ( 1 <= J <= NSUB )
!        Sum ( 1 <= I <= NK )
!          WK(I) * F ( XMID(J) + 0.5 * H * XK(I) )
!
!    The orders of the Gauss-Legendre and Kronrod rules must satisfy
!    NK = 2 * NG + 1.
!
!    The Kronrod rule uses the abscissas of the Gauss-Legendre rule,
!    plus more points, resulting in an efficient and higher order estimate.
!
!    The difference between the Gauss-Legendre and Kronrod estimates
!    is taken as an estimate of the error in the approximation to the
!    integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of integration.
!
!    Input, integer ( kind = 4 ) NSUB, the number of equal subintervals into
!    which the finite interval (A,B) is to be subdivided for
!    higher accuracy.  NSUB must be at least 1.
!
!    Input, integer ( kind = 4 ) NG, the order of the Gauss-Legendre rule.
!    NG must be at least 1.
!
!    Input, real ( kind = 8 ) WG(NG), the weights of the
!    Gauss-Legendre rule.
!
!    Output, real ( kind = 8 ) RESULTG, the approximate value of the
!    integral based on the Gauss-Legendre rule.
!
!    Input, integer ( kind = 4 ) NK, the order of the Kronrod rule.
!    NK must be at least 1.
!
!    Input, real ( kind = 8 ) XK(NK), the abscissas of the
!    Kronrod rule.
!
!    Input, real ( kind = 8 ) WK(NK), the weights of the
!    Kronrod rule.
!
!    Output, real ( kind = 8 ) RESULTK, the approximate value of the
!    integral based on the Kronrod rule.
!
!    Output, real ( kind = 8 ) ERROR, an estimate of the approximation
!    error.  This is computed by taking the sum of the absolute values of
!    the differences between the Gauss-Legendre and Kronrod rules
!    over each subinterval.  This is usually a good estimate of
!    the error in the value RESULTG.  The error in the Kronrod
!    estimate RESULTK is usually much smaller.
!
  implicit none

  integer ( kind = 4 ) ng
  integer ( kind = 4 ) nk

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) error
  real ( kind = 8 ) fk
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) partg
  real ( kind = 8 ) partk
  real ( kind = 8 ) resultg
  real ( kind = 8 ) resultk
  real ( kind = 8 ) wg(ng)
  real ( kind = 8 ) wk(nk)
  real ( kind = 8 ) xk(nk)
  real ( kind = 8 ) xmid
  real ( kind = 8 ) xval

  resultg = 0.0D+00
  resultk = 0.0D+00
  error = 0.0D+00

  if ( a == b ) then
    return
  end if

  if ( nk /= 2 * ng + 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUM_SUB_GK - Fatal error!'
    write ( *, '(a)' ) '  NK must equal 2 * NG + 1.'
    write ( *, '(a,i8)' ) '  The input value was NG = ', ng
    write ( *, '(a,i8)' ) '  The input value was NK = ', nk
    stop
  end if

  h = ( b - a ) / real ( nsub, kind = 8 )

  do j = 1, nsub

    xmid = a + 0.5D+00 * h * real ( 2 * j - 1, kind = 8 )

    partg = 0.0D+00
    partk = 0.0D+00

    do i = 1, nk

      xval = xmid + 0.5D+00 * h * xk(i)
      fk = func ( xval )
      partk = partk + 0.5D+00 * h * wk(i) * fk

      if ( mod ( i, 2 ) == 0 ) then
        partg = partg + 0.5D+00 * h * wg(i/2) * fk
      end if

    end do

    resultg = resultg + partg
    resultk = resultk + partk
    error = error + abs ( partk - partg )

  end do

  return
end
subroutine summer ( func, n, x, w, result )

!*****************************************************************************80
!
!! SUMMER carries out a quadrature rule over a single interval.
!
!  Formula:
!
!    RESULT = sum ( 1 <= I <= N ) W(I) * FUNC ( X(I) )
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
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Input, real ( kind = 8 ) W(N), the weights.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUMMER - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 1.'
    write ( *, '(a,i8)' ) '  The input value was N = ', n
    stop
  end if

  result = 0.0D+00
  do i = 1, n
    result = result + w(i) * func ( x(i) )
  end do

  return
end
subroutine summer_gk ( func, ng, wg, resultg, nk, xk, wk, resultk )

!*****************************************************************************80
!
!! SUMMER_GK carries out Gauss-Kronrod quadrature over a single interval.
!
!  Discussion:
!
!    The abscissas for the Gauss-Legendre rule of order NG are
!    not required, since they are assumed to be the even-indexed
!    entries of the corresponding Kronrod rule.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      function func ( x ).
!
!    Input, integer ( kind = 4 ) NG, the order of the Gauss-Legendre rule.
!
!    Input, real ( kind = 8 ) WG(NG), the weights of the
!    Gauss-Legendre rule.
!
!    Output, real ( kind = 8 ) RESULTG, the approximate value of the
!    integral, based on the Gauss-Legendre rule.
!
!    Input, integer ( kind = 4 ) NK, the order of the Kronrod rule.  NK
!    must equal 2 * NG + 1.
!
!    Input, real ( kind = 8 ) XK(NK), the abscissas of the Kronrod rule.
!
!    Input, real ( kind = 8 ) WK(NK), the weights of the Kronrod rule.
!
!    Output, real ( kind = 8 ) RESULTK, the approximate value of the integral,
!    based on the Kronrod rule.
!
  implicit none

  integer ( kind = 4 ) ng
  integer ( kind = 4 ) nk

  real ( kind = 8 ) fk
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) resultg
  real ( kind = 8 ) resultk
  real ( kind = 8 ) wg(ng)
  real ( kind = 8 ) wk(nk)
  real ( kind = 8 ) xk(nk)

  if ( nk /= 2 * ng + 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUMMER_GK - Fatal error!'
    write ( *, '(a)' ) '  NK must equal 2 * NG + 1.'
    write ( *, '(a,i8)' ) '  The input value was NG = ', ng
    write ( *, '(a,i8)' ) '  The input value was NK = ', nk
    stop
  end if

  resultg = 0.0D+00
  resultk = 0.0D+00

  do i = 1, nk

    fk = func ( xk(i) )

    resultk = resultk + wk(i) * fk

    if ( mod ( i, 2 ) == 0 ) then
      resultg = resultg + wg(i/2) * fk
    end if

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2005
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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
