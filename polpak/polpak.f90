function agm ( a, b )

!*****************************************************************************80
!
!! AGM computes the arithmetic-geometric mean of A and B.
!
!  Discussion:
!
!    The AGM is defined for nonnegative A and B.
!
!    The AGM of numbers A and B is defined by setting
!
!      A(0) = A,
!      B(0) = B
!
!      A(N+1) = ( A(N) + B(N) ) / 2
!      B(N+1) = sqrt ( A(N) * B(N) )
!
!    The two sequences both converge to AGM(A,B).
!
!    In Mathematica, the AGM can be evaluated by
!
!      ArithmeticGeometricMean [ a, b ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the arguments whose AGM is to be computed.
!
!    Output, real ( kind = 8 ) AGM, the arithmetic-geometric mean of A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) agm
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 1000
  real ( kind = 8 ) tol

  if ( a < 0.0D+00 ) then
    stop
  end if

  if ( b < 0.0D+00 ) then
    stop
  end if

  if ( a == 0.0D+00 .or. b == 0.0D+00 ) then
    agm = 0.0D+00
    return
  end if

  it = 0
  tol = 100.0D+00 * epsilon ( tol )

  a2 = a
  b2 = b

  do

    it = it + 1

    c = ( a2 + b2 ) / 2.0D+00
    d = sqrt ( a2 * b2 )

    if ( abs ( c - d ) <= tol * ( c + d ) ) then
      exit
    end if

    if ( it_max < it ) then
      exit
    end if

    a2 = c
    b2 = d

  end do

  agm = c

  return
end
subroutine agm_values ( n_data, a, b, fx )

!*****************************************************************************80
!
!! AGM_VALUES returns some values of the AGM.
!
!  Discussion:
!
!    The AGM is defined for nonnegative A and B.
!
!    The AGM of numbers A and B is defined by setting
!
!      A(0) = A,
!      B(0) = B
!
!      A(N+1) = ( A(N) + B(N) ) / 2
!      B(N+1) = sqrt ( A(N) * B(N) )
!
!    The two sequences both converge to AGM(A,B).
!
!    In Mathematica, the AGM can be evaluated by
!
!      ArithmeticGeometricMean [ a, b ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the argument ofs the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     22.0D+00, &
     83.0D+00, &
     42.0D+00, &
     26.0D+00, &
      4.0D+00, &
      6.0D+00, &
     40.0D+00, &
     80.0D+00, &
     90.0D+00, &
      9.0D+00, &
     53.0D+00, &
      1.0D+00, &
      1.0D+00, &
      1.0D+00, &
      1.5D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     96.0D+00, &
     56.0D+00, &
      7.0D+00, &
     11.0D+00, &
     63.0D+00, &
     45.0D+00, &
     75.0D+00, &
      0.0D+00, &
     35.0D+00, &
      1.0D+00, &
     53.0D+00, &
      2.0D+00, &
      4.0D+00, &
      8.0D+00, &
      8.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     52.274641198704240049D+00, &
     68.836530059858524345D+00, &
     20.659301196734009322D+00, &
     17.696854873743648823D+00, &
     23.867049721753300163D+00, &
     20.717015982805991662D+00, &
     56.127842255616681863D+00, &
      0.000000000000000000D+00, &
     59.269565081229636528D+00, &
     3.9362355036495554780D+00, &
     53.000000000000000000D+00, &
     1.4567910310469068692D+00, &
     2.2430285802876025701D+00, &
     3.6157561775973627487D+00, &
     4.0816924080221632670D+00 /)
  integer ( kind = 4 ) n_data

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function agud ( gamma )

!*****************************************************************************80
!
!! AGUD evaluates the inverse Gudermannian function.
!
!  Discussion:
!
!    The Gudermannian function relates the hyperbolic and trigonometric
!    functions.  For any argument X, there is a corresponding value
!    GAMMA so that
!
!      SINH(X) = TAN(GAMMA).
!
!    This value GAMMA(X) is called the Gudermannian of X.  The inverse
!    Gudermannian function is given as input a value GAMMA and computes
!    the corresponding value X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) GAMMA, the value of the Gudermannian.
!
!    Output, real ( kind = 8 ) AGUD, the argument of the Gudermannian.
!
  implicit none

  real ( kind = 8 ) agud
  real ( kind = 8 ) gamma
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  agud = log ( tan ( 0.25D+00 * pi + 0.5D+00 * gamma ) )

  return
end
function align_enum ( m, n )

!*****************************************************************************80
!
!! ALIGN_ENUM counts the alignments of two sequences of M and N elements.
!
!  Discussion:
!
!    We assume that we have sequences A and B of M and N characters each.
!    An alignment of the two sequences is a rule matching corresponding
!    elements of one sequence to another.  Some elements of either sequence
!    can be matched to a null element.  If A(I1) and A(I2) are matched
!    to B(J1) and B(J2), and I1 < I2, then it must be the case that J1 < J2.
!
!    The 5 alignments of a sequence of 1 to a sequence of 2 are:
!
!          _1_   _2_   __3__   __4__   __5__
!
!      A:  1 -   - 1   - 1 -   - - 1   1 - -
!      B:  1 2   1 2   1 - 2   1 2 -   - 1 2
!
!    The formula is:
!
!      F(0,0) = 1
!      F(1,0) = 1
!      F(0,1) = 1
!      F(M,N) = F(M-1,N) + F(M-1,N-1) + F(M,N-1)
!
!    To compute F(M,N), it is not necessary to keep an M+1 by N+1
!    array in memory.  A vector of length N will do.
!
!    F(N,N) is approximately ( 1 + sqrt(2) )^(2*N+1) / sqrt ( N )
!
!  Example:
!
!    The initial portion of the table is:
!
!  
!  M/N   0    1    2    3    4       5       6       7       8       9      10
!  
!   0    1    1    1    1    1       1       1       1       1       1       1
!   1    1    3    5    7    9      11      13      15      17      19      21
!   2    1    5   13   25   41      61      85     113     145     181     221
!   3    1    7   25   63  129     231     377     575     833    1159    1561
!   4    1    9   41  129  321     681    1289    2241    3649    5641    8361
!   5    1   11   61  231  681    1683    3653    7183   13073   22363   36365
!   6    1   13   85  377 1289    3653    8989   19825   40081   75517  134245
!   7    1   15  113  575 2241    7183   19825   48639  108545  224143  433905
!   8    1   17  145  833 3649   13073   40081  108545  265729  598417 1256465
!   9    1   19  181 1159 5641   22363   75517  224143  598417 1462563 3317445
!  10    1   21  221 1561 8361   36365  134245  433905 1256465 3317445 8097453
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995,
!    ISBN: 0412993910,
!    LC: QH438.4.M33.W38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of elements of the 
!    two sequences.
!
!    Output, integer ( kind = 4 ) ALIGN_ENUM, the number of possible
!    alignments of the sequences.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) align_enum
  integer ( kind = 4 ) fi(0:n)
  integer ( kind = 4 ) fim1j
  integer ( kind = 4 ) fim1jm1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m

  if ( m < 0 ) then
    align_enum = 0
    return
  else if ( n < 0 ) then
    align_enum = 0
    return
  else if ( m == 0 ) then
    align_enum = 1
    return
  else if ( n == 0 ) then
    align_enum = 1
    return
  end if

  fi(0:n) = 1

  do i = 1, m

    fim1jm1 = 1

    do j = 1, n

      fim1j = fi(j)

      fi(j) = fi(j) + fi(j-1) + fim1jm1

      fim1jm1 = fim1j

    end do
  end do

  align_enum = fi(n)

  return
end
function arc_cosine ( c )

!*****************************************************************************80
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
!
!    In particular, you may get the value NaN returned.
!
!    This routine truncates arguments outside the range, avoiding the problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, the argument.
!
!    Output, real ( kind = 8 ) ARC_COSINE, an angle whose cosine is C.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) c
  real ( kind = 8 ) c2

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  arc_cosine = acos ( c2 )

  return
end
function arc_sine ( s )

!*****************************************************************************80
!
!! ARC_SINE computes the arc sine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ASIN routine with an input argument that is
!    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
!
!    In particular, you may get the value NaN returned.
!
!    This routine truncates arguments outside the range, avoiding the problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) S, the argument.
!
!    Output, real ( kind = 8 ) ARC_SINE, an angle whose sine is S.
!
  implicit none

  real ( kind = 8 ) arc_sine
  real ( kind = 8 ) s
  real ( kind = 8 ) s2

  s2 = s
  s2 = max ( s2, -1.0D+00 )
  s2 = min ( s2, +1.0D+00 )

  arc_sine = asin ( s2 )

  return
end
function atan4 ( y, x )

!*****************************************************************************80
!
!! ATAN4 computes the inverse tangent of the ratio Y / X.
!
!  Discussion:
!
!    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * ATAN4 always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Y, X, two quantities which represent 
!    the tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ( kind = 8 ) ATAN4, an angle between 0 and 2 * PI, 
!    whose tangent is (Y/X), and which lies in the appropriate quadrant 
!    so that the signs of its cosine and sine match those of X and Y.
!
  implicit none

  real ( kind = 8 ) abs_x
  real ( kind = 8 ) abs_y
  real ( kind = 8 ) atan4
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_0
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Special cases:
!
  if ( x == 0.0D+00 ) then

    if ( 0.0D+00 < y ) then
      theta = pi / 2.0D+00
    else if ( y < 0.0D+00 ) then
      theta = 3.0D+00 * pi / 2.0D+00
    else if ( y == 0.0D+00 ) then
      theta = 0.0D+00
    end if

  else if ( y == 0.0D+00 ) then

    if ( 0.0D+00 < x ) then
      theta = 0.0D+00
    else if ( x < 0.0D+00 ) then
      theta = pi
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0D+00 < x .and. 0.0D+00 < y ) then
      theta = theta_0
    else if ( x < 0.0D+00 .and. 0.0D+00 < y ) then
      theta = pi - theta_0
    else if ( x < 0.0D+00 .and. y < 0.0D+00 ) then
      theta = pi + theta_0
    else if ( 0.0D+00 < x .and. y < 0.0D+00 ) then
      theta = 2.0D+00 * pi - theta_0
    end if

  end if

  atan4 = theta

  return
end
subroutine bell ( n, b )

!*****************************************************************************80
!
!! BELL returns the Bell numbers from 0 to N.
!
!  Discussion:
!
!    The Bell number B(N) is the number of restricted growth functions on N.
!
!    Note that the Stirling numbers of the second kind, S^m_n, count the
!    number of partitions of N objects into M classes, and so it is
!    true that
!
!      B(N) = S^1_N + S^2_N + ... + S^N_N.
!
!    The Bell numbers were named for Eric Temple Bell.
!
!    The Bell number B(N) is defined as the number of partitions (of
!    any size) of a set of N distinguishable objects.  
!
!    A partition of a set is a division of the objects of the set into 
!    subsets.
!
!  Example:
!
!    There are 15 partitions of a set of 4 objects:
!
!      (1234), 
!      (123) (4), 
!      (124) (3), 
!      (12) (34), 
!      (12) (3) (4), 
!      (134) (2), 
!      (13) (24), 
!      (13) (2) (4), 
!      (14) (23), 
!      (1) (234),
!      (1) (23) (4), 
!      (14) (2) (3), 
!      (1) (24) (3), 
!      (1) (2) (34), 
!      (1) (2) (3) (4).
!
!    and so B(4) = 15.
!
!  First values:
!
!     N         B(N)
!     0           1
!     1           1
!     2           2
!     3           5
!     4          15
!     5          52
!     6         203
!     7         877
!     8        4140
!     9       21147
!    10      115975
!
!  Recursion:
!
!    B(I) = sum ( 1 <= J <= I ) Binomial ( I-1, J-1 ) * B(I-J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of Bell numbers desired.
!
!    Output, integer ( kind = 4 ) B(0:N), the Bell numbers from 0 to N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) b(0:n)
  integer ( kind = 4 ) combo
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    return
  end if

  b(0) = 1

  do i = 1, n
    b(i) = 0
    do j = 1, i
      combo = i4_choose ( i - 1, j - 1 )
      b(i) = b(i) + combo * b(i-j)
    end do
  end do

  return
end
subroutine bell_values ( n_data, n, c )

!*****************************************************************************80
!
!! BELL_VALUES returns some values of the Bell numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the order of the Bell number.
!
!    Output, integer ( kind = 4 ) C, the value of the Bell number.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 11

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
    1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
function benford ( ival )

!*****************************************************************************80
!
!! BENFORD returns the Benford probability of one or more significant digits.
!
!  Discussion:
!
!    Benford's law is an empirical formula explaining the observed
!    distribution of initial digits in lists culled from newspapers,
!    tax forms, stock market prices, and so on.  It predicts the observed
!    high frequency of the initial digit 1, for instance.
!
!    The probabilities of digits 1 through 9 are guaranteed
!    to add up to 1, since
!
!      log10 ( 2/1 ) + log10 ( 3/2) + log10 ( 4/3 ) + ... + log10 ( 10/9 )
!      = log10 ( 2/1 * 3/2 * 4/3 * ... * 10/9 ) = log10 ( 10 ) = 1.
!
!   The formula is:
!
!      Prob ( First significant digits are IVAL ) =
!        log10 ( ( IVAL + 1 ) / IVAL ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Frank Benford,
!    The Law of Anomalous Numbers,
!    Proceedings of the American Philosophical Society,
!    Volume 78, pages 551-572, 1938.
!
!    Ted Hill,
!    The First Digit Phenomenon,
!    American Scientist, 
!    Volume 86, July/August 1998, pages 358 - 363.
!
!    Ralph Raimi,
!    The Peculiar Distribution of First Digits,
!    Scientific American, 
!    December 1969, pages 109-119.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, the string of significant digits to 
!    be checked.  If IVAL is 1, then we are asking for the Benford probability
!    that a value will have first digit 1.  If IVAL is 123, we are asking for
!    the probability that the first three digits will be 123, and so on.
!    Note that IVAL must not be 0 or negative.
!
!    Output, real ( kind = 8 ) BENFORD, the Benford probability that an 
!    item taken from a real world distribution will have the initial 
!    digits IVAL.
!
  implicit none

  real ( kind = 8 ) benford
  integer ( kind = 4 ) ival

  if ( ival <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BENFORD - Fatal error!'
    write ( *, '(a)' ) '  The input argument must be positive.'
    write ( *, '(a,i8)' ) '  Your value was ', ival
    stop
  end if

  benford = log10 ( real ( ival + 1, kind = 8 ) / real ( ival, kind = 8 ) )

  return
end
subroutine bernoulli_number ( n, b )

!*****************************************************************************80
!
!! BERNOULLI_NUMBER computes the Bernoulli numbers B(0) through B(N).
!
!  Discussion:
!
!    The Bernoulli numbers are rational.
!
!    If we define the sum of the M-th powers of the first N integers as:
!
!      SIGMA(M,N) = sum ( 0 <= I <= N ) I^M
!
!    and let C(I,J) be the combinatorial coefficient:
!
!      C(I,J) = I! / ( ( I - J )! * J! )
!
!    then the Bernoulli numbers B(J) satisfy:
!
!      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)^(M+1-J)
!
!  First values:
!
!   B0  1                   =         1.00000000000
!   B1 -1/2                 =        -0.50000000000
!   B2  1/6                 =         1.66666666666
!   B3  0                   =         0
!   B4 -1/30                =        -0.03333333333
!   B5  0                   =         0
!   B6  1/42                =         0.02380952380
!   B7  0                   =         0
!   B8 -1/30                =        -0.03333333333
!   B9  0                   =         0
!  B10  5/66                =         0.07575757575
!  B11  0                   =         0
!  B12 -691/2730            =        -0.25311355311
!  B13  0                   =         0
!  B14  7/6                 =         1.16666666666
!  B15  0                   =         0
!  B16 -3617/510            =        -7.09215686274
!  B17  0                   =         0
!  B18  43867/798           =        54.97117794486
!  B19  0                   =         0
!  B20 -174611/330          =      -529.12424242424
!  B21  0                   =         0
!  B22  854,513/138         =      6192.123
!  B23  0                   =         0
!  B24 -236364091/2730      =    -86580.257
!  B25  0                   =         0
!  B26  8553103/6           =   1425517.16666
!  B27  0                   =         0
!  B28 -23749461029/870     = -27298231.0678
!  B29  0                   =         0
!  B30  8615841276005/14322 = 601580873.901
!
!  Recursion:
!
!    With C(N+1,K) denoting the standard binomial coefficient,
!
!    B(0) = 1.0
!    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
!
!  Warning:
!
!    This recursion, which is used in this routine, rapidly results
!    in significant errors.
!
!  Special Values:
!
!    Except for B(1), all Bernoulli numbers of odd index are 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the highest Bernoulli 
!    number to compute.
!
!    Output, real ( kind = 8 ) B(0:N), B(I) contains the I-th Bernoulli number.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(0:n)
  real ( kind = 8 ) b_sum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) iwork(0:n+1)
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    return
  end if

  b(0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  b(1) = -0.5D+00

  ido = 0
 
  do i = 2, n

    call comb_row ( ido, i+1, iwork )
    ido = 1
 
    if ( mod ( i, 2 ) == 1 ) then
 
      b(i) = 0.0D+00
 
    else
 
      b_sum = 0.0D+00
      do j = 0, i - 1
        b_sum = b_sum + b(j) * real ( iwork(j), kind = 8 )
      end do
 
      b(i) = -b_sum / real ( iwork(i), kind = 8 )
 
    end if

  end do
 
  return
end
subroutine bernoulli_number2 ( n, b )

!*****************************************************************************80
!
!! BERNOULLI_NUMBER2 evaluates the Bernoulli numbers.
!
!  Discussion:
!
!    The Bernoulli numbers are rational.
!
!    If we define the sum of the M-th powers of the first N integers as:
!
!      SIGMA(M,N) = sum ( 0 <= I <= N ) I^M
!
!    and let C(I,J) be the combinatorial coefficient:
!
!      C(I,J) = I! / ( ( I - J )! * J! )
!
!    then the Bernoulli numbers B(J) satisfy:
!
!      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)^(M+1-J)
!
!    Note that the Bernoulli numbers grow rapidly.  Bernoulli number
!    62 is probably the last that can be computed on the VAX without
!    overflow.
!
!    A different method than that used in BERN is employed.
!
!  First values:
!
!   B0  1                   =         1.00000000000
!   B1 -1/2                 =        -0.50000000000
!   B2  1/6                 =         1.66666666666
!   B3  0                   =         0
!   B4 -1/30                =        -0.03333333333
!   B5  0                   =         0
!   B6  1/42                =         0.02380952380
!   B7  0                   =         0
!   B8 -1/30                =        -0.03333333333
!   B9  0                   =         0
!  B10  5/66                =         0.07575757575
!  B11  0                   =         0
!  B12 -691/2730            =        -0.25311355311
!  B13  0                   =         0
!  B14  7/6                 =         1.16666666666
!  B15  0                   =         0
!  B16 -3617/510            =        -7.09215686274
!  B17  0                   =         0
!  B18  43867/798           =        54.97117794486
!  B19  0                   =         0
!  B20 -174611/330          =      -529.12424242424
!  B21  0                   =         0
!  B22  854,513/138         =      6192.123
!  B23  0                   =         0
!  B24 -236364091/2730      =    -86580.257
!  B25  0                   =         0
!  B26  8553103/6           =   1425517.16666
!  B27  0                   =         0
!  B28 -23749461029/870     = -27298231.0678
!  B29  0                   =         0
!  B30  8615841276005/14322 = 601580873.901
!
!  Recursion:
!
!    With C(N+1,K) denoting the standard binomial coefficient,
!
!    B(0) = 1.0
!    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
!
!  Special Values:
!
!    Except for B(1), all Bernoulli numbers of odd index are 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order Bernoulli number 
!    to compute.
!
!    Output, real ( kind = 8 ) B(0:N), the requested Bernoulli numbers.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) altpi
  real ( kind = 8 ) b(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: kmax = 400
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sgn
  real ( kind = 8 ) sum2
  real ( kind = 8 ) t
  real ( kind = 8 ) term
  real ( kind = 8 ), parameter :: tol = 1.0D-06

  if ( n < 0 ) then
    return
  end if

  b(0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  b(1) = -0.5D+00

  if ( n < 2 ) then
    return
  end if

  altpi = log ( 2.0D+00 * pi )
!
!  Initial estimates for B(I), I = 2 to N
!
  b(2) = log ( 2.0D+00 )
  do i = 3, n
    if ( mod ( i, 2 ) == 1 ) then
      b(i) = 0.0D+00
    else
      b(i) = log ( real ( i * ( i - 1 ), kind = 8 ) ) + b(i-2)
    end if
  end do

  b(2) = 1.0D+00 / 6.0D+00

  if ( n <= 3 ) then
    return
  end if

  b(4) = -1.0D+00 / 30.0D+00

  sgn = -1.0D+00
 
  do i = 6, n, 2
 
    sgn = -sgn
    t = 2.0D+00 * sgn * exp ( b(i) - real ( i, kind = 8 ) * altpi )
 
    sum2 = 1.0D+00

    do k = 2, kmax

      term = real ( k, kind = 8 )**(-i)
      sum2 = sum2 + term

      if ( term <= tol * sum2 ) then
        exit
      end if

    end do
 
    b(i) = t * sum2
 
  end do
 
  return
end
subroutine bernoulli_number3 ( n, b )

!*****************************************************************************80
!
!! BERNOULLI_NUMBER3 computes the value of the Bernoulli number B(N).
!
!  Discussion:
!
!    The Bernoulli numbers are rational.
!
!    If we define the sum of the M-th powers of the first N integers as:
!
!      SIGMA(M,N) = sum ( 0 <= I <= N ) I^M
!
!    and let C(I,J) be the combinatorial coefficient:
!
!      C(I,J) = I! / ( ( I - J )! * J! )
!
!    then the Bernoulli numbers B(J) satisfy:
!
!      SIGMA(M,N) = 
!        1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)^(M+1-J)
!
!  First values:
!
!     B0  1                   =         1.00000000000
!     B1 -1/2                 =        -0.50000000000
!     B2  1/6                 =         1.66666666666
!     B3  0                   =         0
!     B4 -1/30                =        -0.03333333333
!     B5  0                   =         0
!     B6  1/42                =         0.02380952380
!     B7  0                   =         0
!     B8 -1/30                =        -0.03333333333
!     B9  0                   =         0
!    B10  5/66                =         0.07575757575
!    B11  0                   =         0
!    B12 -691/2730            =        -0.25311355311
!    B13  0                   =         0
!    B14  7/6                 =         1.16666666666
!    B15  0                   =         0
!    B16 -3617/510            =        -7.09215686274
!    B17  0                   =         0
!    B18  43867/798           =        54.97117794486
!    B19  0                   =         0
!    B20 -174611/330          =      -529.12424242424
!    B21  0                   =         0
!    B22  854513/138          =      6192.123
!    B23  0                   =         0
!    B24 -236364091/2730      =    -86580.257
!    B25  0                   =         0
!    B26  8553103/6           =   1425517.16666
!    B27  0                   =         0
!    B28 -23749461029/870     = -27298231.0678
!    B29  0                   =         0
!    B30  8615841276005/14322 = 601580873.901
!
!  Recursion:
!
!    With C(N+1,K) denoting the standard binomial coefficient,
!
!    B(0) = 1.0
!    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
!
!  Special Values:
!
!    Except for B(1), all Bernoulli numbers of odd index are 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Bernoulli number 
!    to compute.
!
!    Output, real ( kind = 8 ) B, the desired Bernoulli number.
!
  implicit none

  real ( kind = 8 ) b
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 1000
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) sum2
  real ( kind = 8 ) term
  real ( kind = 8 ), parameter :: tol = 5.0D-07

  if ( n < 0 ) then

    b = 0.0D+00

  else if ( n == 0 ) then

    b = 1.0D+00

  else if ( n == 1 ) then

    b = -0.5D+00

  else if ( n == 2 ) then

    b = 1.0D+00 / 6.0D+00

  else if ( mod ( n, 2 ) == 1 ) then

    b = 0.0D+00

  else

    sum2 = 0.0D+00

    do it = 1, it_max

      term = 1.0D+00 / real ( it**n, kind = 8 )
      sum2 = sum2 + term

      if ( abs ( term ) < tol .or. abs ( term ) < tol * abs ( sum2 ) ) then
        exit
      end if

    end do

    b = 2.0D+00 * sum2 * r8_factorial ( n ) / ( 2.0D+00 * pi )**n

    if ( mod ( n, 4 ) == 0 ) then
      b = - b
    end if

  end if

  return
end
subroutine bernoulli_number_values ( n_data, n, c )

!*****************************************************************************80
!
!! BERNOULLI_NUMBER_VALUES returns some values of the Bernoulli numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the order of the Bernoulli number.
!
!    Output, real ( kind = 8 ) C, the value of the Bernoulli number.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 10

  real ( kind = 8 ) c
  real ( kind = 8 ), save, dimension ( nmax ) :: c_vec = (/ &
    1.0000000000D+00, -0.5000000000D+00,  0.1666666667D+00, &
    0.0000000000D+00, -0.0333333333D+00, -0.02380952381D+00, &
   -0.0333333333D+00,  0.0757575757D+00, -529.1242424D+00, &
    601580873.9D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 6,  8, 10, 20, 30 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine bernoulli_poly ( n, x, bx )

!*****************************************************************************80
!
!! BERNOULLI_POLY evaluates the Bernoulli polynomial of order N at X.
!
!  Discussion:
!
!    Thanks to Bart Vandewoestyne for pointing out an error in the previous
!    documentation, 31 January 2008.
!
!    Special values of the Bernoulli polynomial include:
!
!      B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
!
!      B'(N,X) = N * B(N-1,X)
!
!      B(N,X+1) - B(N,X) = N * X^(N-1)
!      B(N,X) = (-1)^N * B(N,1-X)
!
!    A formula for the Bernoulli polynomial in terms of the Bernoulli
!    numbers is:
!
!      B(N,X) = sum ( 0 <= K <= N ) B(K) * C(N,K) * X^(N-K)
!
!    The first few polynomials include:
!
!      B(0,X) = 1
!      B(1,X) = X    - 1/2
!      B(2,X) = X^2 -   X      +  1/6
!      B(3,X) = X^3 - 3/2*X^2 +  1/2*X
!      B(4,X) = X^4 - 2*X^3   +      X^2 - 1/30
!      B(5,X) = X^5 - 5/2*X^4 +  5/3*X^3 - 1/6*X
!      B(6,X) = X^6 - 3*X^5   +  5/2*X^4 - 1/2*X^2 + 1/42
!      B(7,X) = X^7 - 7/2*X^6 +  7/2*X^5 - 7/6*X^3 + 1/6*X
!      B(8,X) = X^8 - 4*X^7   + 14/3*X^6 - 7/3*X^4 + 2/3*X^2 - 1/30
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Bernoulli polynomial to
!    be evaluated.  N must be 0 or greater.
!
!    Input, real ( kind = 8 ) X, the value of X at which the polynomial is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) BX, the value of B(N,X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) iwork(0:n)
  real ( kind = 8 ) work(0:n)
  real ( kind = 8 ) x

  call bernoulli_number ( n, work )
 
  ido = 0
  call comb_row ( ido, n, iwork )
 
  bx = 1.0D+00
  do i = 1, n
    bx = bx * x + work(i) * real ( iwork(i), kind = 8 )
  end do
 
  return
end
subroutine bernoulli_poly2 ( n, x, bx )

!*****************************************************************************80
!
!! BERNOULLI_POLY2 evaluates the N-th Bernoulli polynomial at X.
!
!  Discussion:
!
!    Thanks to Bart Vandewoestyne for pointing out an error in the previous
!    documentation, 31 January 2008.
!
!    Special values of the Bernoulli polynomial include:
!
!      B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
!
!      B'(N,X) = N * B(N-1,X)
!
!      B(N,X+1) - B(N,X) = N * X^(N-1)
!      B(N,X) = (-1)^N * B(N,1-X)
!
!    A formula for the Bernoulli polynomial in terms of the Bernoulli
!    numbers is:
!
!      B(N,X) = sum ( 0 <= K <= N ) B(K) * C(N,K) * X^(N-K)
!
!    The first few polynomials include:
!
!      B(0,X) = 1
!      B(1,X) = X    - 1/2
!      B(2,X) = X^2 -   X      +  1/6
!      B(3,X) = X^3 - 3/2*X^2 +  1/2*X
!      B(4,X) = X^4 - 2*X^3   +      X^2 - 1/30
!      B(5,X) = X^5 - 5/2*X^4 +  5/3*X^3 - 1/6*X
!      B(6,X) = X^6 - 3*X^5   +  5/2*X^4 - 1/2*X^2 + 1/42
!      B(7,X) = X^7 - 7/2*X^6 +  7/2*X^5 - 7/6*X^3 + 1/6*X
!      B(8,X) = X^8 - 4*X^7   + 14/3*X^6 - 7/3*X^4 + 2/3*X^2 - 1/30
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Bernoulli polynomial to
!    be evaluated.  N must be 0 or greater.
!
!    Input, real ( kind = 8 ) X, the value at which the polynomial is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) BX, the value of B(N,X).
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bx
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  fact = 1.0D+00

  call bernoulli_number3 ( 0, b )

  bx = b

  do i = 1, n
    fact = fact * real ( n + 1 - i, kind = 8 ) / real ( i, kind = 8 )
    call bernoulli_number3 ( i, b )
    bx = bx * x + fact * b
  end do

  return
end
subroutine bernstein_poly ( n, x, bern )

!*****************************************************************************80
!
!! BERNSTEIN_POLY evaluates the Bernstein polynomials at a point X.
!
!  Discussion:
!
!    The Bernstein polynomials are assumed to be based on [0,1].
!
!    The formula is:
!
!      B(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
!
!  First values:
!
!    B(0,0,X) = 1
!
!    B(1,0,X) =      1-X
!    B(1,1,X) =               X
!
!    B(2,0,X) =     (1-X)^2
!    B(2,1,X) = 2 * (1-X)   * X
!    B(2,2,X) =               X^2
!
!    B(3,0,X) =     (1-X)^3
!    B(3,1,X) = 3 * (1-X)^2 * X
!    B(3,2,X) = 3 * (1-X)   * X^2
!    B(3,3,X) =               X^3
!
!    B(4,0,X) =     (1-X)^4
!    B(4,1,X) = 4 * (1-X)^3 * X
!    B(4,2,X) = 6 * (1-X)^2 * X^2
!    B(4,3,X) = 4 * (1-X)   * X^3
!    B(4,4,X) =               X^4
!
!  Special values:
!
!    B(N,I,X) has a unique maximum value at X = I/N.
!
!    B(N,I,X) has an I-fold zero at 0 and and N-I fold zero at 1.
!
!    B(N,I,1/2) = C(N,K) / 2^N
!
!    For a fixed X and N, the polynomials add up to 1:
!
!      Sum ( 0 <= I <= N ) B(N,I,X) = 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the Bernstein polynomials 
!    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
!    each of degree N, which form a basis for polynomials on [0,1].
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) BERN(0:N), the values of the N+1 
!    Bernstein polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bern(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  if ( n == 0 ) then
 
    bern(0) = 1.0D+00
 
  else if ( 0 < n ) then
 
    bern(0) = 1.0D+00 - x
    bern(1) = x
 
    do i = 2, n
      bern(i) = x * bern(i-1)
      do j = i - 1, 1, -1
        bern(j) =             x   * bern(j-1) &
                + ( 1.0D+00 - x ) * bern(j)
      end do
      bern(0) = ( 1.0D+00 - x ) * bern(0)
    end do
 
  end if
 
  return
end
subroutine bernstein_poly_values ( n_data, n, k, x, b )

!*****************************************************************************80
!
!! BERNSTEIN_POLY_VALUES returns some values of the Bernstein polynomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and 
!    N_DATA is set to the index of the test data.  On each subsequent 
!    call, N_DATA is incremented and that test data is returned.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Output, integer ( kind = 4 ) K, the index of the polynomial.
! 
!    Output, real ( kind = 8 ) X, the argument of the polynomial.
!
!    Output, real ( kind = 8 ) B, the value of the polynomial B(N,K,X).
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 15

  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( nmax ) :: b_vec = (/ &
    1.0D+00, &
    0.75D+00,       0.25D+00, &
    0.5625D+00,     0.3750D+00,   0.0625D+00, &
    0.421875D+00,   0.421875D+00, 0.140625D+00,  0.015625D+00, &
    0.31640625D+00, 0.421875D+00, 0.2109375D+00, 0.046875D+00, 0.00390625D+00 /)
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save, dimension ( nmax ) :: k_vec = (/ &
    0, &
    0, 1, &
    0, 1, 2, &
    0, 1, 2, 3, &
    0, 1, 2, 3, 4 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
    0, &
    1, 1, &
    2, 2, 2, &
    3, 3, 3, 3, &
    4, 4, 4, 4, 4 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( nmax ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, 0.25D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    k = 0
    x = 0.0D+00
    b = 0.0D+00
  else
    n = n_vec(n_data)
    k = k_vec(n_data)
    x = x_vec(n_data)
    b = b_vec(n_data)
  end if

  return
end
function beta ( x, y )

!*****************************************************************************80
!
!! BETA returns the value of the Beta function.
!
!  Discussion:
!
!    The Beta function can be defined in terms of the Gamma function:
!
!      BETA(X,Y) = ( GAMMA(X) * GAMMA(Y) ) / GAMMA(X+Y)
!
!      Both X and Y must be greater than 0.
!
!    The function has the following properties:
!
!      BETA(X,Y) = BETA(Y,X).
!      BETA(X,Y) = Integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
!      BETA(X,Y) = GAMMA(X) * GAMMA(Y) / GAMMA(X+Y)
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
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the two parameters that define 
!    the Beta function.  X and Y must be greater than 0.
!
!    Output, real ( kind = 8 ) BETA, the value of the Beta function.
!
  implicit none

  real ( kind = 8 ) beta
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( x <= 0.0D+00 .or. y <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA - Fatal error!'
    write ( *, '(a)' ) '  Both X and Y must be greater than 0.'
    stop
  end if

  beta = exp ( r8_gamma_log ( x ) &
             + r8_gamma_log ( y ) &
             - r8_gamma_log ( x + y ) )

  return
end
subroutine beta_values ( n_data, x, y, fxy )

!*****************************************************************************80
!
!! BETA_VALUES returns some values of the Beta function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2001
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and 
!    N_DATA is set to the index of the test data.  On each subsequent 
!    call, N_DATA is incremented and that test data is returned.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, real ( kind = 8 ) X, Y, the arguments of the function.
!
!    Output, real ( kind = 8 ) FXY, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 17

  real ( kind = 8 ), save, dimension ( nmax ) :: b_vec = (/ &
    5.000000D+00, 2.500000D+00, 1.666667D+00, 1.250000D+00, &
    5.000000D+00, 2.500000D+00, 1.000000D+00, 1.666667D-01, &
    0.333333D-01, 7.142857D-03, 1.587302D-03, 0.238095D-01, &
    5.952381D-03, 1.984127D-03, 7.936508D-04, 3.607504D-04, &
    8.325008D-05 /)
  real ( kind = 8 ) fxy
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( nmax ) :: x_vec = (/ &
    0.2D+00, 0.4D+00, 0.6D+00, 0.8D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, &
    3.0D+00, 4.0D+00, 5.0D+00, 6.0D+00, &
    6.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
    7.0D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( nmax ) :: y_vec = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    0.2D+00, 0.4D+00, 1.0D+00, 2.0D+00, &
    3.0D+00, 4.0D+00, 5.0D+00, 2.0D+00, &
    3.0D+00, 4.0D+00, 5.0D+00, 6.0D+00, &
    7.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    x = 0.0D+00
    y = 0.0D+00
    fxy = 0.0D+00
  else
    x = x_vec(n_data)
    y = y_vec(n_data)
    fxy = b_vec(n_data)
  end if

  return
end
subroutine bpab ( n, a, b, x, bern )

!*****************************************************************************80
!
!! BPAB evaluates at X the Bernstein polynomials based in [A,B].
!
!  Discussion:
!
!    The formula is:
!
!      BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
!
!  First values:
!
!    B(0,0,X) =   1
!
!    B(1,0,X) = (      B-X                ) / (B-A)
!    B(1,1,X) = (                 X-A     ) / (B-A)
!
!    B(2,0,X) = (     (B-X)^2             ) / (B-A)^2
!    B(2,1,X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
!    B(2,2,X) = (                (X-A)^2  ) / (B-A)^2
!
!    B(3,0,X) = (     (B-X)^3             ) / (B-A)^3
!    B(3,1,X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
!    B(3,2,X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
!    B(3,3,X) = (                (X-A)^3  ) / (B-A)^3
!
!    B(4,0,X) = (     (B-X)^4             ) / (B-A)^4
!    B(4,1,X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
!    B(4,2,X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4
!    B(4,3,X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4
!    B(4,4,X) = (                (X-A)^4  ) / (B-A)^4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the Bernstein polynomials 
!    to be used.  For any N, there is a set of N+1 Bernstein polynomials, 
!    each of degree N, which form a basis for polynomials on [A,B].
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval on which the
!    polynomials are to be based.  A and B should not be equal.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials 
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) BERN(0:N), the values of the N+1
!    Bernstein polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bern(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  if ( b == a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BPAB - Fatal error!'
    write ( *, '(a,g14.6)' ) '  A = B = ', a
    stop
  end if

  if ( n == 0 ) then
 
    bern(0) = 1.0D+00
 
  else if ( 0 < n ) then
 
    bern(0) = ( b - x ) / ( b - a )
    bern(1) = ( x - a ) / ( b - a )
 
    do i = 2, n
      bern(i) = ( x - a ) * bern(i-1) / ( b - a )
      do j = i - 1, 1, -1
        bern(j) = ( ( b - x     ) * bern(j)     &
                  + (     x - a ) * bern(j-1) ) &
                  / ( b     - a )
      end do
      bern(0) = ( b - x ) * bern(0) / ( b - a )
    end do
 
  end if
 
  return
end
subroutine cardan ( n, x, s, cx )

!*****************************************************************************80
!
!! CARDAN evaluates the Cardan polynomials.
!
!  Discussion:
!
!    If we write the N-th polynomial in terms of its coefficients:
!
!      C(N,S,X) = sum ( 0 <= I <= N ) D(N,I) * S^(N-I)/2 * X^I
!
!    then we have the rcursion:
!
!      D(0,0) = 1
!
!      D(1,1) = 1
!      D(1,0) = 0
!
!      D(N,N) = 1
!      D(N,K) = D(N-1,K-1) - D(N-2,K)
!
!  Example:
!
!     N  C(N,S,X)
!
!     0  2
!     1  X
!     2  X^2  -  2 S
!     3  X^3  -  3 S X
!     4  X^4  -  4 S X^2 +  2 S^2
!     5  X^5  -  5 S X^3 +  5 S^2 X
!     6  X^6  -  6 S X^4 +  9 S^2 X^2 -  2 S^3
!     7  X^7  -  7 S X^5 + 14 S^2 X^3 -  7 S^3 X
!     8  X^8  -  8 S X^6 + 20 S^2 X^4 - 16 S^3 X^3 +  2 S^4
!     9  X^9  -  9 S X^7 + 27 S^2 X^5 - 30 S^3 X^3 +  9 S^4 X
!    10  X^10 - 10 S X^8 + 35 S^2 X^6 - 50 S^3 X^4 + 25 S^4 X^2 -  2 S^5
!    11  X^11 - 11 S X^9 + 44 S^2 X^7 - 77 S^3 X^5 + 55 S^4 X^3 - 11 S^5 X
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Osler,
!    Cardan Polynomials and the Reduction of Radicals,
!    Mathematics Magazine, 
!    Volume 74, Number 1, February 2001, pages 26-32.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest polynomial to compute.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials 
!    are to be computed.
!
!    Input, real ( kind = 8 ) S, the value of the parameter, which 
!    must be positive.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the Cardan 
!    polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  real ( kind = 8 ) s
  real ( kind = 8 ) s2
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  s2 = sqrt ( s )
  x2 = 0.5D+00 * x / s2

  call cheby_t_poly ( 1, n, x2, cx )

  fact = 1.0D+00

  do i = 0, n
    cx(i) = 2.0D+00 * fact * cx(i)
    fact = fact * s2
  end do
 
  return
end
subroutine cardan_poly_coef ( n, s, c )

!*****************************************************************************80
!
!! CARDAN_POLY_COEF computes the coefficients of the N-th Cardan polynomial.
!
!  First terms:
!
!    2
!    0      1
!   -2 S    0      1
!    0     -3 S    0      1
!    2 S^2  0     -4 S    0      1
!    0      5 S^2  0     -5 S    0      1
!   -2 S^3  0      9 S^2  0     -6 S    0      1
!    0      7 S^3  0     14 S^2  0     -7 S    0      1
!    2 S^4  0    -16 S^3  0     20 S^2  0     -8 S    0       1
!    0      9 S^4  0    -30 S^3  0     27 S^2  0     -9 S     0     1
!   -2 S^5  0     25 S^4  0    -50 S^3  0     35 S^2  0     -10 S   0   1
!    0    -11 S^5  0     55 S^4  0    -77 S^3  0    +44 S^2   0   -11 S 0 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Osler,
!    Cardan Polynomials and the Reduction of Radicals,
!    Mathematics Magazine, 
!    Volume 74, Number 1, February 2001, pages 26-32.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial
!
!    Input, real ( kind = 8 ) S, the value of the parameter, which 
!    must be positive.
!
!    Output, real ( kind = 8 ) C(0:N), the coefficients.  C(0) is the 
!    constant term, and C(N) is the coefficient of X^N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) cm1(0:n)
  real ( kind = 8 ) cm2(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) s

  if ( n < 0 ) then
    return
  end if

  c(0) = 2.0D+00
  c(1:n) = 0.0D+00

  if ( n == 0 ) then
    return
  end if

  cm1(0:n) = c(0:n)

  c(0) = 0.0D+00
  c(1) = 1.0D+00
  c(2:n) = 0.0D+00

  do i = 2, n

    cm2(0:i-2) = cm1(0:i-2)
    cm1(0:i-1) = c(0:i-1)

    c(0) = 0.0D+00
    c(1:i) = cm1(0:i-1)
    c(0:i-2) = c(0:i-2) - s * cm2(0:i-2)

  end do

  return
end
subroutine catalan ( n, c )

!*****************************************************************************80
!
!! CATALAN computes the Catalan numbers, from C(0) to C(N).
!
!  Discussion:
!
!    The Catalan number C(N) counts:
!
!    1) the number of binary trees on N vertices;
!    2) the number of ordered trees on N+1 vertices;
!    3) the number of full binary trees on 2N+1 vertices;
!    4) the number of well formed sequences of 2N parentheses;
!    5) the number of ways 2N ballots can be counted, in order,
!       with N positive and N negative, so that the running sum
!       is never negative;
!    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
!    7) the number of monotone functions from [1..N} to [1..N} which 
!       satisfy f(i) <= i for all i;
!    8) the number of ways to triangulate a polygon with N+2 vertices.
!
!    The formula is:
!
!      C(N) = (2*N)! / ( (N+1) * (N!) * (N!) ) 
!           = 1 / (N+1) * COMB ( 2N, N )
!           = 1 / (2N+1) * COMB ( 2N+1, N+1).
!
!  First values:
!
!     C(0)     1
!     C(1)     1
!     C(2)     2
!     C(3)     5
!     C(4)    14
!     C(5)    42
!     C(6)   132
!     C(7)   429
!     C(8)  1430
!     C(9)  4862
!    C(10) 16796
!
!  Recursion:
!
!    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
!    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
!
!  Example:
!
!    N = 3
!
!    ()()()
!    ()(())
!    (()())
!    (())()
!    ((()))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of Catalan numbers desired.
!
!    Output, integer ( kind = 4 ) C(0:N), the Catalan numbers from C(0) to C(N).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0) = 1
!
!  The extra parentheses ensure that the integer division is
!  done AFTER the integer multiplication.
!
  do i = 1, n
    c(i) = ( c(i-1) * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
  end do
 
  return
end
function catalan_constant ( )

!*****************************************************************************80
!
!! CATALAN_CONSTANT returns the value of Catalan's constant.
!
!  Discussion:
!
!    Catalan's constant, which may be denoted by G, is defined as
!
!      G = sum ( 0 <= K ) ( -1 )^K / ( 2 * K + 1 )^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Output, real ( kind = 8 ) CATALAN_CONSTANT, the value of Catalan's
!    constant.
!
  implicit none

  real ( kind = 8 ) catalan_constant

  catalan_constant = 0.915965594177D+00

  return
end
subroutine catalan_row_next ( ido, n, irow )

!*****************************************************************************80
!
!! CATALAN_ROW_NEXT computes row N of Catalan's triangle.
!
!  Example:
!
!    I\J 0   1   2   3   4   5   6
!
!    0   1
!    1   1   1
!    2   1   2   2
!    3   1   3   5   5
!    4   1   4   9  14  14
!    5   1   5  14  28  42  42
!    6   1   6  20  48  90 132 132
!
!  Recursion:
!
!    C(0,0) = 1
!    C(I,0) = 1
!    C(I,J) = 0 for I < J
!    C(I,J) = C(I,J-1) + C(I-1,J)
!    C(I,I) is the I-th Catalan number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IDO, indicates whether this is a call for
!    the 'next' row of the triangle.
!    IDO = 0, this is a startup call.  Row N is desired, but
!    presumably this is a first call, or row N-1 was not computed
!    on the previous call.
!    IDO = 1, this is not the first call, and row N-1 was computed
!    on the previous call.  In this case, much work can be saved
!    by using the information from the previous values of IROW
!    to build the next values.
!
!    Input, integer ( kind = 4 ) N, the index of the row of the triangle 
!    desired.  
!
!    Input/output, integer ( kind = 4 ) IROW(0:N), the row of coefficients.
!    If IDO = 0, then IROW is not required to be set on input.
!    If IDO = 1, then IROW must be set on input to the value of
!    row N-1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) irow(0:n)
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    return
  end if

  if ( ido == 0 ) then
 
    irow(0) = 1
    irow(1:n) = 0
 
    do i = 1, n

      irow(0) = 1

      do j = 1, i - 1
        irow(j) = irow(j) + irow(j-1)
      end do

      irow(i) = irow(i-1)

    end do
 
  else
 
    irow(0) = 1

    do j = 1, n - 1
      irow(j) = irow(j) + irow(j-1)
    end do

    if ( 1 <= n ) then
      irow(n) = irow(n-1)
    end if
 
  end if
 
  return
end
subroutine catalan_values ( n_data, n, c )

!*****************************************************************************80
!
!! CATALAN_VALUES returns some values of the Catalan numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the order of the Catalan number.
!
!    Output, integer ( kind = 4 ) C, the value of the Catalan number.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 11

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
    1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine charlier ( n, a, x, value )

!*****************************************************************************80
!
!! CHARLIER evaluates Charlier polynomials at a point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J Simoes Pereira,
!    Algorithm 234: Poisson-Charliers Polynomials,
!    Communications of the ACM,
!    Volume 7, Number 7, page 420, July 1964.
!
!    Walter Gautschi,
!    Orthogonal Polynomials: Computation and Approximation,
!    Oxford, 2004,
!    ISBN: 0-19-850672-4,
!    LC: QA404.5 G3555.
!
!    Gabor Szego,
!    Orthogonal Polynomials,
!    American Mathematical Society, 1975,
!    ISBN: 0821810235,
!    LC: QA3.A5.v23.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum order of the polynomial.  
!    N must be at least 0.
!
!    Input, real ( kind = 8 ) A, the parameter.  A must not be 0.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) VALUE(0:N), the value of the polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  real ( kind = 8 ) value(0:n)
  real ( kind = 8 ) x

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHARLIER - Fatal error!'
    write ( *, '(a)' ) '  Parameter A cannot be zero.'
    stop
  end if

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHARLIER - Fatal error!'
    write ( *, '(a)' ) '  Parameter N must be nonnegative.'
    stop
  end if

  value(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  value(1) = - x / a

  if ( n == 1 ) then 
    return
  end if
  
  do i = 1, n - 1
    value(i+1) = ( ( real ( i, kind = 8 ) + a - x ) * value(i) &
                   - real ( i, kind = 8 ) * value(i-1) ) / a
  end do

  return
end
subroutine cheby_t_poly ( m, n, x, v )

!*****************************************************************************80
!
!! CHEBY_T_POLY evaluates Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    Chebyshev polynomials are useful as a basis for representing the
!    approximation of functions since they are well conditioned, in the sense
!    that in the interval [-1,1] they each have maximum absolute value 1.
!    Hence an error in the value of a coefficient of the approximation, of
!    size epsilon, is exactly reflected in an error of size epsilon between
!    the computed approximation and the theoretical approximation.
!
!    Typical usage is as follows, where we assume for the moment
!    that the interval of approximation is [-1,1].  The value
!    of N is chosen, the highest polynomial to be used in the
!    approximation.  Then the function to be approximated is
!    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
!    Chebyshev polynomial.  Let these values be denoted by F(XJ).
!
!    The coefficients of the approximation are now defined by
!
!      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I,XJ)
!
!    except that C(0) is given a value which is half that assigned
!    to it by the above formula,
!
!    and the representation is
!
!    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
!
!    Now note that, again because of the fact that the Chebyshev polynomials
!    have maximum absolute value 1, if the higher order terms of the
!    coefficients C are small, then we have the option of truncating
!    the approximation by dropping these terms, and we will have an
!    exact value for maximum perturbation to the approximation that
!    this will cause.
!
!    It should be noted that typically the error in approximation
!    is dominated by the first neglected basis function (some multiple of
!    T(N+1,X) in the example above).  If this term were the exact error,
!    then we would have found the minimax polynomial, the approximating
!    polynomial of smallest maximum deviation from the original function.
!    The minimax polynomial is hard to compute, and another important
!    feature of the Chebyshev approximation is that it tends to behave
!    like the minimax polynomial while being easy to compute.
!
!    To evaluate a sum like 
!
!      sum ( 0 <= J <= N ) C(J) T(J,X), 
!
!    Clenshaw's recurrence formula is recommended instead of computing the
!    polynomial values, forming the products and summing.
!
!    Assuming that the coefficients C(J) have been computed
!    for J = 0 to N, then the coefficients of the representation of the
!    indefinite integral of the function may be computed by
!
!      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1, 
!
!    with
! 
!      C(N+1)=0
!      B(0) arbitrary.  
!
!    Also, the coefficients of the representation of the derivative of the 
!    function may be computed by:
!
!      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0, 
!
!    with
!
!      D(N+1) = D(N)=0.
!
!    Some of the above may have to adjusted because of the irregularity of C(0).
!
!    The formula is:
!
!      T(N,X) = COS(N*ARCCOS(X))
!
!  Differential equation:
!
!    (1-X*X) Y'' - X Y' + N N Y = 0
!
!  First terms:
!
!    T(0,X) =  1
!    T(1,X) =  1 X
!    T(2,X) =  2 X^2 -   1
!    T(3,X) =  4 X^3 -   3 X
!    T(4,X) =  8 X^4 -   8 X^2 +  1
!    T(5,X) = 16 X^5 -  20 X^3 +  5 X
!    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
!    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
!
!  Inequality:
!
!    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
!
!  Orthogonality:
!
!    For integration over [-1,1] with weight
!
!      W(X) = 1 / sqrt(1-X*X), 
!
!    if we write the inner product of T(I,X) and T(J,X) as
!
!      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
!
!    then the result is:
!
!      < T(I,X), T(J,X) > = 0    if I /= J
!      < T(I,X), T(J,X) > = PI/2 if I == J /= 0
!      < T(I,X), T(J,X) > = PI   if I == J == 0
!
!    A discrete orthogonality relation is also satisfied at each of
!    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
!                              = 0 if I /= J
!                              = N/2 if I == J /= 0
!                              = N if I == J == 0
!
!  Recursion:
!
!    T(0,X) = 1,
!    T(1,X) = X,
!    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
!
!    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
!
!  Special values:
!
!    T(N,1) = 1
!    T(N,-1) = (-1)^N
!    T(2N,0) = (-1)^N
!    T(2N+1,0) = 0
!    T(N,X) = (-1)^N * T(N,-X)
!
!  Zeroes:
!
!    M-th zero of T(N,X) is X = cos((2*M-1)*PI/(2*N)), M = 1 to N.
!
!  Extrema:
!
!    M-th extremum of T(N,X) is X = cos(PI*M/N), M = 0 to N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2001
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
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest polynomial to compute.
!
!    Input, real ( kind = 8 ) X(1:M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(1:M,0:N), the values of the polynomials.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) x(1:m)
  real ( kind = 8 ) v(1:m,0:n)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
 
  do j = 2, n
    v(1:m,j) = 2.0D+00 * x(1:m) * v(1:m,j-1) - v(1:m,j-2)
  end do
 
  return
end
subroutine cheby_t_poly_coef ( n, c )

!*****************************************************************************80
!
!! CHEBY_T_POLY_COEF evaluates coefficients of Chebyshev polynomials T(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     1
!     2     -1     0      2
!     3      0    -3      0      4
!     4      1     0     -8      0       8
!     5      0     5      0    -20       0    16
!     6     -1     0     18      0     -48     0     32
!     7      0    -7      0     56       0  -112      0    64
!
!  Recursion:
!
!    T(0,X) = 1,
!    T(1,X) = X,
!    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 May 2003
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the Chebyshev T
!    polynomials.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 1.0D+00
 
  do i = 2, n
    c(i,0)     =                        - c(i-2,0)
    c(i,1:i-2) = 2.0D+00 * c(i-1,0:i-3) - c(i-2,1:i-2)
    c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
    c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
  end do
 
  return
end
subroutine cheby_t_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! CHEBY_T_POLY_VALUES returns values of Chebyshev polynomials T(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( nmax ) :: fx_vec = (/ &
     1.0000000000D+00,  0.8000000000D+00,  0.2800000000D+00, &
    -0.3520000000D+00, -0.8432000000D+00, -0.9971200000D+00, &
    -0.7521920000D+00, -0.2063872000D+00,  0.4219724800D+00, &
     0.8815431680D+00,  0.9884965888D+00,  0.7000513741D+00, &
     0.1315856097D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( nmax ) :: x_vec = (/ &
    0.8D+00,  0.8D+00,  0.8D+00, &
    0.8D+00,  0.8D+00,  0.8D+00, &
    0.8D+00,  0.8D+00,  0.8D+00, &
    0.8D+00,  0.8D+00,  0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cheby_t_poly_zero ( n, z )

!*****************************************************************************80
!
!! CHEBY_T_POLY_ZERO returns zeroes of Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    The I-th zero of T(N,X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) Z(N), the zeroes of T(N,X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) z(n)

  do i = 1, n
    angle = real ( 2 * i - 1, kind = 8 ) * pi / real ( 2 * n, kind = 8 )
    z(i) = cos ( angle )
  end do

  return
end
subroutine cheby_u_poly ( n, x, cx )

!*****************************************************************************80
!
!! CHEBY_U_POLY evaluates Chebyshev polynomials U(n,x).
!
!  Discussion:
!
!    The formula is:
!
!      If |X| <= 1, then
!
!        U(N,X) = sin ( (N+1) * arccos(X) ) / sqrt ( 1 - X^2 )
!               = sin ( (N+1) * arccos(X) ) / sin ( arccos(X) )
!
!      else
!
!        U(N,X) = sinh ( (N+1) * arccosh(X) ) / sinh ( arccosh(X) )
!
!  Differential equation:
!
!    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
!
!  First terms:
!
!    U(0,X) =   1
!    U(1,X) =   2 X
!    U(2,X) =   4 X^2 -   1
!    U(3,X) =   8 X^3 -   4 X
!    U(4,X) =  16 X^4 -  12 X^2 +  1
!    U(5,X) =  32 X^5 -  32 X^3 +  6 X
!    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
!    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
!
!  Orthogonality:
!
!    For integration over [-1,1] with weight
!
!      W(X) = sqrt(1-X*X), 
!
!    we have
!
!      < U(I,X), U(J,X) > = integral ( -1 <= X <= 1 ) W(X) U(I,X) U(J,X) dX 
!
!    then the result is:
!
!      < U(I,X), U(J,X) >  =  0    if I /= J
!      < U(I,X), U(J,X) >  =  PI/2 if I == J
!
!  Recursion:
!
!    U(0,X) = 1,
!    U(1,X) = 2 * X,
!    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
!
!  Special values:
!
!    U(N,1) = N + 1
!    U(2N,0) = (-1)^N
!    U(2N+1,0) = 0
!    U(N,X) = (-1)^N * U(N,-X)
!
!  Zeroes:
!
!    M-th zero of U(N,X) is X = cos( M*PI/(N+1)), M = 1 to N
!
!  Extrema:
!
!    M-th extremum of U(N,X) is X = cos( M*PI/N), M = 0 to N
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = PI/2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2002
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
!    Input, integer ( kind = 4 ) N, the highest polynomial to compute.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials 
!    are to be computed.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the N+1 Chebyshev
!    polynomials.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  cx(1) = 2.0D+00 * x

  do i = 2, n
    cx(i) = 2.0D+00 * x * cx(i-1) - cx(i-2)
  end do
 
  return
end
subroutine cheby_u_poly_coef ( n, c )

!*****************************************************************************80
!
!! CHEBY_U_POLY_COEF evaluates coefficients of Chebyshev polynomials U(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     2
!     2     -1     0      4
!     3      0    -4      0      8
!     4      1     0    -12      0      16
!     5      0     6      0    -32       0    32
!     6     -1     0     24      0     -80     0     64
!     7      0    -8      0     80       0  -192      0   128
!
!  Recursion:
!
!    U(0,X) = 1,
!    U(1,X) = 2*X,
!    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2003
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the Chebyshev U
!    polynomials.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 2.0D+00
 
  do i = 2, n
    c(i,0)     =                        - c(i-2,0)
    c(i,1:i-2) = 2.0D+00 * c(i-1,0:i-3) - c(i-2,1:i-2)
    c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
    c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
  end do
 
  return
end
subroutine cheby_u_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! CHEBY_U_POLY_VALUES returns values of the Chebyshev polynomial U(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ChebyshevU[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2012
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.1600000000000000D+01, &
     0.1560000000000000D+01, &
     0.8960000000000000D+00, &
    -0.1264000000000000D+00, &
    -0.1098240000000000D+01, &
    -0.1630784000000000D+01, &
    -0.1511014400000000D+01, &
    -0.7868390400000000D+00, &
     0.2520719360000000D+00, &
     0.1190154137600000D+01, &
     0.1652174684160000D+01, &
     0.1453325357056000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine cheby_u_poly_zero ( n, z )

!*****************************************************************************80
!
!! CHEBY_U_POLY_ZERO returns zeroes of Chebyshev polynomials U(n,x).
!
!  Discussion:
!
!    The I-th zero of U(N,X) is cos((I-1)*PI/(N-1)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) Z(N), the zeroes of U(N,X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) z(n)

  do i = 1, n
    angle = real ( i, kind = 8 ) * pi / real ( n + 1, kind = 8 )
    z(i) = cos ( angle )
  end do

  return
end
subroutine chebyshev_discrete ( n, m, x, v )

!*****************************************************************************80
!
!! CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials at a point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Walter Gautschi,
!    Orthogonal Polynomials: Computation and Approximation,
!    Oxford, 2004,
!    ISBN: 0-19-850672-4,
!    LC: QA404.5 G3555.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order of the polynomials to 
!    be evaluated.  0 <= N <= M.
!
!    Input, integer ( kind = 4 ) M, the maximum order of the polynomials.
!    0 <= M.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) V(0:N), the value of the polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) x
  real ( kind = 8 ) v(0:n)

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEBYSHEV_DISCRETE - Fatal error!'
    write ( *, '(a)' ) '  Parameter M must be nonnegative.'
    stop
  end if

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEBYSHEV_DISCRETE - Fatal error!'
    write ( *, '(a)' ) '  Parameter N must be nonnegative.'
    stop
  end if

  if ( m < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHEBYSHEV_DISCRETE - Fatal error!'
    write ( *, '(a)' ) '  Parameter N must be no greater than M.'
    stop
  end if

  v(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  v(1) = 2.0D+00 * x + real ( 1 - m, kind = 8 )

  if ( n == 1 ) then
    return
  end if

  do i = 1, n - 1
    v(i+1) = ( &
      real ( 2 * i + 1, kind = 8 ) &
      * ( 2.0D+00 * x + real ( 1 - m, kind = 8 ) ) * v(i) &
      - real ( i * ( m + i ) * ( m - i ), kind = 8 ) * v(i-1) &
    ) / real ( i + 1, kind = 8 )
  end do

  return
end
function collatz_count ( n )

!*****************************************************************************80
!
!! COLLATZ_COUNT counts the number of terms in a Collatz sequence.
!
!  Discussion:
!
!    The rules for generation of the Collatz sequence are recursive.
!    If T is the current entry of the sequence, (T is
!    assumed to be a positive integer), then the next
!    entry, U is determined as follows:
!    
!      if T is 1 (or less)
!        terminate the sequence;
!      else if T is even
!        U = T/2.
!      else (if T is odd and not 1)
!        U = 3*T+1;
!
!     N  Sequence                                                Length
!
!     1                                                               1
!     2   1                                                           2
!     3  10,  5, 16,  8,  4,  2,  1                                   8
!     4   2   1                                                       3
!     5  16,  8,  4,  2,  1                                           6
!     6   3, 10,  5, 16,  8,  4,  2,  1                               9
!     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
!     8   4,  2,  1                                                   4
!     9  28, 14,  7, ...                                             20
!    10   5, 16,  8,  4,  2,  1                                       7
!    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
!    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the first element of the sequence.
!
!    Output, integer ( kind = 4 ) COLLATZ_COUNT, the number of elements in
!    the Collatz sequence that begins with N.
!
  implicit none

  integer ( kind = 4 ) collatz_count
  integer ( kind = 4 ) count
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_local

  count = 1
  n_local = n

  do

    if ( n_local <= 1 ) then
      exit
    else if ( mod ( n_local, 2 ) == 0 ) then
      n_local = n_local / 2
    else
      n_local = 3 * n_local + 1
    end if

    count = count + 1

  end do

  collatz_count = count
 
  return
end
subroutine collatz_count_max ( n, i_max, j_max )

!*****************************************************************************80
!
!! COLLATZ_COUNT_MAX seeks the maximum Collatz count for 1 through N.
!
!  Discussion:
!
!    For each integer I, we compute a sequence of values that 
!    terminate when we reach 1.  The number of steps required to
!    reach 1 is the "rank" of I, and we are searching the numbers
!    from 1 to N for the number with maximum rank.
!
!    For a given I, the sequence is produced by:
!
!    1) J = 1, X(J) = I;
!    2) If X(J) = 1, stop.
!    3) J = J + 1; 
!       if X(J-1) was even, X(J) = X(J-1)/2;
!       else                X(J) = 3 * X(J-1) + 1;
!    4) Go to 3
!
!  Example:
!
!            N     I_MAX J_MAX
!
!           10         9    20
!          100        97   119
!        1,000       871   179
!       10,000     6,171   262
!      100,000    77,031   351
!    1,000,000   837,799   525
!   10,000,000 8,400,511   686
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum integer to check.
!
!    Output, integer ( kind = 4 ) I_MAX, J_MAX, an integer I with the maximum 
!    rank, and the value of the maximum rank.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x

  i_max = -1
  j_max = -1

  do i = 1, n

    j = 1
    x = i

    do while ( x /= 1 )
      j = j + 1
      if ( mod ( x, 2 ) == 0 ) then
        x = x / 2
      else
        x = 3 * x + 1
      end if
    end do

    if ( j_max < j ) then
      i_max = i
      j_max = j
    end if

  end do

  return
end
subroutine collatz_count_values ( n_data, n, count )

!*****************************************************************************80
!
!! COLLATZ_COUNT_VALUES returns some values of the Collatz count function.
!
!  Discussion:
!
!    The rules for generation of the Collatz sequence are recursive.
!    If T is the current entry of the sequence, (T is
!    assumed to be a positive integer), then the next
!    entry, U is determined as follows:
!
!      if T is 1 (or less)
!        terminate the sequence;
!      else if T is even
!        U = T/2.
!      else (if T is odd and not 1)
!        U = 3*T+1;
!
!    The Collatz count is the length of the Collatz sequence for a given
!    starting value.  By convention, we include the initial value in the
!    count, so the minimum value of the count is 1.
!
!     N  Sequence                                                 Count
!
!     1                                                               1
!     2   1                                                           2
!     3  10,  5, 16,  8,  4,  2,  1                                   8
!     4   2   1                                                       3
!     5  16,  8,  4,  2,  1                                           6
!     6   3, 10,  5, 16,  8,  4,  2,  1                               9
!     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
!     8   4,  2,  1                                                   4
!     9  28, 14,  7, ...                                             20
!    10   5, 16,  8,  4,  2,  1                                       7
!    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
!    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the initial value of a Collatz sequence.
!
!    Output, integer ( kind = 4 ) COUNT, the length of the Collatz sequence
!    starting with N.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) count
  integer ( kind = 4 ), save, dimension ( n_max ) :: count_vec = (/ &
      1,   2,   8,   3,   6,   9,   17,   4,  20,   7, &
    112,  25,  26,  27,  17,  28,  111,  18,  83,  29 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     27,  50, 100, 200, 300, 400, 500, 600, 700, 800 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    count = 0
  else
    n = n_vec(n_data)
    count = count_vec(n_data)
  end if

  return
end
subroutine comb_row ( ido, n, irow )

!*****************************************************************************80
!
!! COMB_ROW computes row N of Pascal's triangle.
!
!  Discussion:
!
!    Row N contains the combinatorial coefficients
!
!      C(N,0), C(N,1), C(N,2), ... C(N,N)
!
!    The sum of the elements of row N is equal to 2**N.
!
!    The formula is
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  First terms:
!
!     N K:0  1   2   3   4   5   6   7  8  9 10
!
!     0   1
!     1   1  1
!     2   1  2   1
!     3   1  3   3   1
!     4   1  4   6   4   1
!     5   1  5  10  10   5   1
!     6   1  6  15  20  15   6   1
!     7   1  7  21  35  35  21   7   1
!     8   1  8  28  56  70  56  28   8  1
!     9   1  9  36  84 126 126  84  36  9  1
!    10   1 10  45 120 210 252 210 120 45 10  1
!
!  Recursion:
!
!    C(N,K) = C(N-1,K-1)+C(N-1,K)
!
!  Special values:
!
!    C(N,0) = C(N,N) = 1
!    C(N,1) = C(N,N-1) = N
!    C(N,N-2) = sum ( 1 <= I <= N ) N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IDO, indicates whether this is a call for
!    the 'next' row of the triangle.
!    * 0 means this is a startup call.  Row N is desired, but
!      presumably this is a first call, or row N-1 was not computed
!      on the previous call.
!    * 1 means this is not the first call, and row N-1 was computed
!      on the previous call.  In this case, much work can be saved
!      by using the information from the previous values of IROW
!      to build the next values.
!
!    Input, integer ( kind = 4 ) N, the row of the triangle desired.  The 
!    triangle begins with row 0.
!
!    Output, integer ( kind = 4 ) IROW(N+1), the row of coefficients.
!    IROW(I) = C(N,I-1).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) irow(n+1)
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    return
  end if
 
  if ( ido == 1 ) then
 
    do i = n, 2, -1
      irow(i) = irow(i) + irow(i-1)
    end do
 
    irow(n+1) = 1
 
  else
 
    irow(1) = 1
    irow(2:n+1) = 0
 
    do j = 1, n
      do i = j + 1, 2, -1
        irow(i) = irow(i) + irow(i-1)
      end do
    end do
 
  end if
 
  return
end
subroutine commul ( n, nfactor, factor, ncomb )

!*****************************************************************************80
!
!! COMMUL computes a multinomial combinatorial coefficient.
!
!  Discussion:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where FACTOR(1) objects are indistinguishable of type 1,
!    ... and FACTOR(K) are indistinguishable of type NFACTOR.
!
!    The formula is:
!
!      NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, determines the numerator.
!
!    Input, integer ( kind = 4 ) NFACTOR, the number of factors in the 
!    numerator.
!
!    Input, integer ( kind = 4 ) FACTOR(NFACTOR).
!    FACTOR contains the NFACTOR values used in the denominator.
!    Note that the sum of these entries should be N,
!    and that all entries should be nonnegative.
!
!    Output, integer ( kind = 4 ) NCOMB, the value of the multinomial 
!    coefficient.
!
  implicit none

  integer ( kind = 4 ) nfactor

  real ( kind = 8 ) arg
  real ( kind = 8 ) fack
  real ( kind = 8 ) facn
  integer ( kind = 4 ) factor(nfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isum
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncomb
  real ( kind = 8 ) r8_gamma_log

  if ( nfactor < 1 ) then
    ncomb = 1
    return
  end if

  do i = 1, nfactor

    if ( factor(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMMUL - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Entry ', I, ' of FACTOR = ', factor(i)
      write ( *, '(a)' ) '  But this value must be nonnegative.'
      stop
    end if

  end do
 
  isum = sum ( factor(1:nfactor) )

  if ( isum /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMMUL - Fatal error!'
    write ( *, '(a,i8)' ) '  The sum of the FACTOR entries is ', isum
    write ( *, '(a,i8)' ) '  But it must equal N = ', n
    stop
  end if
 
  arg = real ( n + 1, kind = 8 )
  facn = r8_gamma_log ( arg )
 
  do i = 1, nfactor
 
    arg = real ( factor(i) + 1, kind = 8 )
    fack = r8_gamma_log ( arg )
    facn = facn - fack
 
  end do
 
  ncomb = nint ( exp ( facn ) )
 
  return
end
function cos_deg ( angle )

!*****************************************************************************80
!
!! COS_DEG returns the cosine of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees.
!
!    Output, real ( kind = 8 ) COS_DEG, the cosine of the angle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) cos_deg
  real ( kind = 8 ), parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00

  cos_deg = cos ( degrees_to_radians * angle )

  return
end
function cos_power_int ( a, b, n )

!*****************************************************************************80
!
!! COS_POWER_INT evaluates the cosine power integral.
!
!  Discussion:
!
!    The function is defined by
!
!      COS_POWER_INT(A,B,N) = Integral ( A <= t <= B ) ( cos ( t ))^n dt
!
!    The algorithm uses the following fact:
!
!      Integral cos^n ( t ) = - (1/n) * (
!        cos^(n-1)(t) * sin(t) + ( n-1 ) * Integral cos^(n-2) ( t ) dt )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, integer ( kind = 4 ) N, the power of the sine function.
!
!    Output, real ( kind = 8 ) COS_POWER_INT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  real ( kind = 8 ) cos_power_int
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mlo
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb

  real ( kind = 8 ) value

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COS_POWER_INT - Fatal error!'
    write ( *, '(a)' ) '  Power N < 0.'
    value = 0.0D+00
    stop
  end if

  sa = sin ( a )
  sb = sin ( b )
  ca = cos ( a )
  cb = cos ( b )

  if ( mod ( n, 2 ) == 0 ) then
    value = b - a
    mlo = 2
  else
    value = sb - sa
    mlo = 3
  end if

  do m = mlo, n, 2
    value = ( real ( m - 1, kind = 8 ) * value &
              - ca**(m-1) * sa + cb**(m-1) * sb ) &
      / real ( m, kind = 8 )
  end do

  cos_power_int = value

  return
end
subroutine cos_power_int_values ( n_data, a, b, n, fx )

!*****************************************************************************80
!
!! COS_POWER_INT_VALUES returns some values of the cosine power integral.
!
!  Discussion:
!
!    The function has the form
!
!      COS_POWER_INT(A,B,N) = Integral ( A <= t <= B ) ( cos(T) )^N dt
!
!    In Mathematica, the function can be evaluated by:
!
!      Integrate [ ( Cos[x] )^n, { x, a, b } ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, integer ( kind = 4 ) N, the power.
!
!    Output, real ( kind = 8 ) FX, the function value.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00, &
     3.141592653589793D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     3.141592653589793D+00, &
     0.0D+00, &
     1.570796326794897D+00, &
     0.0D+00, &
     1.178097245096172D+00, &
     0.0D+00, &
     0.9817477042468104D+00, &
     0.0D+00, &
     0.8590292412159591D+00, &
     0.0D+00, &
     0.7731263170943632D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, &
     1, &
     2, &
     3, &
     4, &
     5, &
     6, &
     7, &
     8, &
     9, &
    10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    n = 0
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine delannoy ( m, n, a )

!*****************************************************************************80
!
!! DELANNOY returns the Delannoy numbers up to orders (M,N).
!
!  Discussion:
!
!    The Delannoy number A(M,N) counts the number of distinct paths
!    from (0,0) to (M,N) in which the only steps used are
!    (1,1), (1,0) and (0,1).
!
!  First values:
!
!      \N 0  1   2   3    4     5     6      7      8
!     M-+--------------------------------------------
!     0 | 1  1   1   1    1     1     1      1      1
!     1 | 1  3   5   7    9    11    13     15     17
!     2 | 1  5  13  25   41    61    85    113    145
!     3 | 1  7  25  63  129   231   377    575    833
!     4 | 1  9  41 129  321   681  1289   2241   3649
!     5 | 1 11  61 231  681  1683  3653   7183  13073
!     6 | 1 13  85 377 1289  3653  8989  19825  40081
!     7 | 1 15 113 575 2241  7183 19825  48639 108545
!     8 | 1 17 145 833 3649 13073 40081 108545 265729
!
!  Recursion:
!
!    A(0,0) = 1
!    A(M,0) = 1
!    A(0,N) = 1
!    A(M,N) = A(M-1,N) + A(M,N-1) + A(M-1,N-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, define the highest order number to 
!    compute.
!
!    Output, integer ( kind = 4 ) A(0:M,0:N), the Delannoy numbers.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(0:m,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( m < 0 ) then
    return
  end if

  if ( n < 0 ) then
    return
  end if

  a(0,0) = 1

  a(1:m,0) = 1
  a(0,1:n) = 1

  do i = 1, m
    do j = 1, n
      a(i,j) = a(i-1,j) + a(i,j-1) + a(i-1,j-1)
    end do
  end do

  return
end
function e_constant ( )

!*****************************************************************************80
!
!! E_CONSTANT returns the value of the base of the natural logarithm system.
!
!  Discussion:
!
!    E = Limit ( N -> oo ) ( 1 + 1 / N )^N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) E_CONSTANT, the base of the natural 
!    logarithm system.
!
  implicit none

  real ( kind = 8 ) e_constant

  e_constant = 2.718281828459045235360287D+00
 
  return
end
subroutine erf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ERF_VALUES returns some values of the ERF or "error" function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and 
!    N_DATA is set to the index of the test data.  On each subsequent 
!    call, N_DATA is incremented and that test data is returned.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 21

  real ( kind = 8 ), save, dimension ( nmax ) :: bvec = (/ &
    0.0000000000D+00, 0.1124629160D+00, 0.2227025892D+00, 0.3286267595D+00, &
    0.4283923550D+00, 0.5204998778D+00, 0.6038560908D+00, 0.6778011938D+00, &
    0.7421009647D+00, 0.7969082124D+00, 0.8427007929D+00, 0.8802050696D+00, &
    0.9103139782D+00, 0.9340079449D+00, 0.9522851198D+00, 0.9661051465D+00, &
    0.9763483833D+00, 0.9837904586D+00, 0.9890905016D+00, 0.9927904292D+00, &
    0.9953222650D+00 /)
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( nmax ) :: xvec = (/ &
    0.0D+00, 0.1D+00, 0.2D+00, 0.3D+00, &
    0.4D+00, 0.5D+00, 0.6D+00, 0.7D+00, &
    0.8D+00, 0.9D+00, 1.0D+00, 1.1D+00, &
    1.2D+00, 1.3D+00, 1.4D+00, 1.5D+00, &
    1.6D+00, 1.7D+00, 1.8D+00, 1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = xvec(n_data)
    fx = bvec(n_data)
  end if

  return
end
function error_f ( x )

!*****************************************************************************80
!
!! ERROR_F evaluates the error function ERF.
!
!  Discussion:
!
!    Since some compilers already supply a routine named ERF which evaluates
!    the error function, this routine has been given a distinct, if
!    somewhat unnatural, name.
!    
!    The function is defined by:
!
!      ERF(X) = ( 2 / sqrt ( PI ) ) * Integral ( 0 <= t <= X ) exp ( - t^2 ) dt
!
!    Properties of the function include:
!
!      Limit ( X -> -oo ) ERF(X) =          -1.0;
!                               ERF(0) =           0.0;
!                               ERF(0.476936...) = 0.5;
!      Limit ( X -> +oo ) ERF(X) =          +1.0.
!
!      0.5D+00 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev Approximations for the Error Function,
!    Mathematics of Computation, 
!    1969, pages 631-638.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) ERROR_F, the value of the error function.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
    3.16112374387056560D+00, &
    1.13864154151050156D+02, &
    3.77485237685302021D+02, &
    3.20937758913846947D+03, &
    1.85777706184603153D-01 /)
  real ( kind = 8 ), parameter, dimension ( 4 ) :: b = (/ &
    2.36012909523441209D+01, &
    2.44024637934444173D+02, &
    1.28261652607737228D+03, &
    2.84423683343917062D+03 /)
  real ( kind = 8 ), parameter, dimension ( 9 ) :: c = (/ &
    5.64188496988670089D-01, &
    8.88314979438837594D+00, &
    6.61191906371416295D+01, &
    2.98635138197400131D+02, &
    8.81952221241769090D+02, &
    1.71204761263407058D+03, &
    2.05107837782607147D+03, &
    1.23033935479799725D+03, &
    2.15311535474403846D-08 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.57449261107098347D+01, &
    1.17693950891312499D+02, &
    5.37181101862009858D+02, &
    1.62138957456669019D+03, &
    3.29079923573345963D+03, &
    4.36261909014324716D+03, &
    3.43936767414372164D+03, &
    1.23033935480374942D+03 /)
  real ( kind = 8 ) del
  real ( kind = 8 ) error_f
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter, dimension ( 6 ) :: p = (/ &
    3.05326634961232344D-01, &
    3.60344899949804439D-01, &
    1.25781726111229246D-01, &
    1.60837851487422766D-02, &
    6.58749161529837803D-04, &
    1.63153871373020978D-02 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    2.56852019228982242D+00, &
    1.87295284992346047D+00, &
    5.27905102951428412D-01, &
    6.05183413124413191D-02, &
    2.33520497626869185D-03 /)
  real ( kind = 8 ), parameter :: sqrpi = 0.56418958354775628695D+00
  real ( kind = 8 ), parameter :: thresh = 0.46875D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xabs
  real ( kind = 8 ), parameter :: xbig = 26.543D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsq

  xabs = abs ( ( x ) )
!
!  Evaluate ERF(X) for |X| <= 0.46875.
!
  if ( xabs <= thresh ) then

    if ( epsilon ( xabs ) < xabs ) then
      xsq = xabs * xabs
    else
      xsq = 0.0D+00
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do

    error_f = x * ( xnum + a(4) ) / ( xden + b(4) )
!
!  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
!
  else if ( xabs <= 4.0D+00 ) then

    xnum = c(9) * xabs
    xden = xabs
    do i = 1, 7
      xnum = ( xnum + c(i) ) * xabs
      xden = ( xden + d(i) ) * xabs
    end do

    error_f = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = real ( int ( xabs * 16.0D+00 ), kind = 8 ) / 16.0D+00
    del = ( xabs - xsq ) * ( xabs + xsq )
    error_f = exp ( - xsq * xsq ) * exp ( - del ) * error_f

    error_f = ( 0.5D+00 - error_f ) + 0.5D+00

    if ( x < 0.0D+00 ) then
      error_f = - error_f
    end if
!
!  Evaluate ERFC(X) for 4.0D+00 < |X|.
!
  else

    if ( xbig <= xabs ) then

      if ( 0.0D+00 < x ) then
        error_f = 1.0D+00
      else
        error_f = - 1.0D+00
      end if

    else

      xsq = 1.0D+00 / ( xabs * xabs )

      xnum = p(6) * xsq
      xden = xsq
      do i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
      end do

      error_f = xsq * ( xnum + p(5) ) / ( xden + q(5) )
      error_f = ( sqrpi - error_f ) / xabs
      xsq = real ( int ( xabs * 16.0D+00 ), kind = 8 ) / 16.0D+00
      del = ( xabs - xsq ) * ( xabs + xsq )
      error_f = exp ( - xsq * xsq ) * exp ( - del ) * error_f

      error_f = ( 0.5D+00 - error_f ) + 0.5D+00

      if ( x < 0.0D+00 ) then
        error_f = - error_f
      end if

    end if

  end if

  return
end
function error_f_inverse ( y )

!*****************************************************************************80
!
!! ERROR_F_INVERSE inverts the error function ERF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Y, the value of the error function.
!
!    Output, real ( kind = 8 ) ERROR_F_INVERSE, the value X such that
!    ERROR_F(X) = Y.
!
  implicit none

  real ( kind = 8 ) error_f_inverse
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = ( y + 1.0D+00 ) / 2.0D+00

  call normal_01_cdf_inv ( z, x )

  error_f_inverse = x / sqrt ( 2.0D+00 )

  return
end
function euler_constant ( )

!*****************************************************************************80
!
!! EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
!
!  Discussion:
!
!    The Euler-Mascheroni constant is often denoted by a lower-case
!    Gamma.  Gamma is defined as
!
!      Gamma = limit ( M -> +oo )
!        ( sum ( 1 <= N <= M ) 1 / N ) - log ( M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EULER_CONSTANT, the value of the 
!    Euler-Mascheroni constant.
!
  implicit none

  real ( kind = 8 ) euler_constant

  euler_constant = 0.577215664901532860606512090082402431042D+00

  return
end
subroutine euler_number ( n, e )

!*****************************************************************************80
!
!! EULER_NUMBER computes the Euler numbers.
!
!  Discussion:
!
!    The Euler numbers can be evaluated in Mathematica by:
!
!      EulerE[n]
!
!    These numbers rapidly get too big to store in an ordinary integer!
!
!    The terms of odd index are 0.
!
!    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
!
!  First terms:
!
!    E0  = 1
!    E1  = 0
!    E2  = -1
!    E3  = 0
!    E4  = 5
!    E5  = 0
!    E6  = -61
!    E7  = 0
!    E8  = 1385
!    E9  = 0
!    E10 = -50521
!    E11 = 0
!    E12 = 2702765
!    E13 = 0
!    E14 = -199360981
!    E15 = 0
!    E16 = 19391512145
!    E17 = 0
!    E18 = -2404879675441
!    E19 = 0
!    E20 = 370371188237525
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the last Euler number
!    to compute.
!
!    Output, integer ( kind = 4 ) E(0:N), the Euler numbers.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) e(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    return
  end if

  e(0) = 1

  if ( n == 0 ) then
    return
  end if

  e(1) = 0
 
  if ( n == 1 ) then
    return
  end if

  e(2) = -1

  do i = 3, n

    e(i) = 0

    if ( mod ( i, 2 ) == 0 ) then

      do j = 2, i, 2
        e(i) = e(i) - i4_choose ( i, j ) * e(i-j)
      end do

    end if

  end do
 
  return
end
function euler_number2 ( n )

!*****************************************************************************80
!
!! EULER_NUMBER2 computes the Euler numbers.
!
!  Discussion:
!
!    The Euler numbers can be evaluated in Mathematica by:
!
!      EulerE[n]
!
!  First terms:
!
!    E0  = 1
!    E1  = 0
!    E2  = -1
!    E3  = 0
!    E4  = 5
!    E5  = 0
!    E6  = -61
!    E7  = 0
!    E8  = 1385
!    E9  = 0
!    E10 = -50521
!    E11 = 0
!    E12 = 2702765
!    E13 = 0
!    E14 = -199360981
!    E15 = 0
!    E16 = 19391512145
!    E17 = 0
!    E18 = -2404879675441
!    E19 = 0
!    E20 = 370371188237525
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the Euler number to compute.
!
!    Output, real ( kind = 8 ) EULER_NUMBER2, the value of E(N).
!
  implicit none

  real ( kind = 8 ) euler_number2
  real ( kind = 8 ), save, dimension ( 0:6 ) :: e = &
    (/ 1.0D+00, -1.0D+00, 5.0D+00, -61.0D+00, 1385.0D+00, &
       -50521.0D+00, 2702765.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: itmax = 1000
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) sum1
  real ( kind = 8 ) term

  if ( n < 0 ) then
    euler_number2 = 0.0D+00
    return
  end if

  if ( n == 0 ) then
    euler_number2 = e(0)
    return
  end if

  if ( mod ( n, 2 ) == 1 ) then
    euler_number2 = 0.0D+00
    return
  end if

  if ( n <= 12 ) then
    euler_number2 = e(n/2)
    return
  end if

  sum1 = 0.0D+00
  do i = 1, itmax

    term = 1.0D+00 / real ( ( 2 * i - 1 )**( n + 1 ), kind = 8 )

    if ( mod ( i, 2 ) == 1 ) then
      sum1 = sum1 + term
    else
      sum1 = sum1 - term
    end if

    if ( abs ( term ) < 1.0D-10 ) then
      exit
    else if ( abs ( term ) < 1.0D-08 * abs ( sum1 ) ) then
      exit
    end if

  end do

  euler_number2 = 2.0D+00**( n + 2 ) * sum1 * r8_factorial ( n ) / pi**( n + 1 )

  if ( mod ( n, 4 ) /= 0 ) then
    euler_number2 = - euler_number2
  end if

  return
end
subroutine euler_number_values ( n_data, n, c )

!*****************************************************************************80
!
!! EULER_NUMBER_VALUES returns some values of the Euler numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the order of the Euler number.
!
!    Output, integer ( kind = 4 ) C, the value of the Euler number.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 8

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
    1, 0, -1, 5, 61, 1385, -50521, 2702765 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0, 1, 2, 4, 6, 8, 10, 12 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
function euler_poly ( n, x )

!*****************************************************************************80
!
!! EULER_POLY evaluates the N-th Euler polynomial at X.
!
!  First values:
!
!    E(0,X) = 1
!    E(1,X) = X - 1/2
!    E(2,X) = X^2 - X 
!    E(3,X) = X^3 - 3/2 X^2 + 1/4
!    E(4,X) = X^4 - 2*X^3 + X
!    E(5,X) = X^5 - 5/2 X^4 + 5/2 X^2 - 1/2
!    E(6,X) = X^6 - 3 X^5 + 5 X^3 - 3 X
!    E(7,X) = X^7 - 7/2 X^6 + 35/4 X^4 - 21/2 X^2 + 17/8
!    E(8,X) = X^8 - 4 X^7 + 14 X^5 - 28 X^3 + 17 X
!
!  Special values:
!
!    E'(N,X) = N * E(N-1,X)
!
!    E(N,1/2) = E(N) / 2^N, where E(N) is the N-th Euler number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Euler polynomial to
!    be evaluated.  N must be 0 or greater.
!
!    Input, real ( kind = 8 ) X, the value at which the polynomial is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) EULER_POLY, the value of E(N,X).
!
  implicit none

  real ( kind = 8 ) bx1
  real ( kind = 8 ) bx2
  real ( kind = 8 ) euler_poly
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  call bernoulli_poly2 ( n+1, x, bx1 )
  call bernoulli_poly2 ( n+1, 0.5D+00 * x, bx2 )

  euler_poly = 2.0D+00 * ( bx1 - bx2 * 2.0D+00**( n + 1 ) ) &
    / real ( n + 1, kind = 8 )

  return
end
subroutine eulerian ( n, e )

!*****************************************************************************80
!
!! EULERIAN computes the Eulerian number E(N,K).
!
!  Discussion:
!
!    A run in a permutation is a sequence of consecutive ascending values.
!
!    E(N,K) is the number of permutations of N objects which contain
!    exactly K runs.
!
!  Examples:
!
!     N = 7
!
!     1     0     0     0     0     0     0
!     1     1     0     0     0     0     0
!     1     4     1     0     0     0     0
!     1    11    11     1     0     0     0
!     1    26    66    26     1     0     0
!     1    57   302   302    57     1     0
!     1   120  1191  2416  1191   120     1
!
!  Recursion:
!
!    E(N,K) = K * E(N-1,K) + (N-K+1) * E(N-1,K-1).
!
!  Properties:
!
!    E(N,1) = E(N,N) = 1.
!    E(N,K) = 0 if K <= 0 or N < K.
!    sum ( 1 <= K <= N ) E(N,K) = N!.
!    X^N = sum ( 0 <= K <= N ) COMB(X+K-1, N ) E(N,K)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, 1986
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows desired.
!
!    Output, integer ( kind = 4 ) E(N,N), the first N rows of Eulerian numbers.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) e(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( n < 1 ) then
    return
  end if
!
!  Construct rows 1, 2, ..., N of the Eulerian triangle.
!
  e(1,1) = 1
  e(1,2:n) = 0

  do i = 2, n
    e(i,1) = 1
    do j = 2, n
      e(i,j) = j * e(i-1,j) + ( i - j + 1 ) * e(i-1,j-1)
    end do
  end do

  return
end
recursive function f_hofstadter ( n ) result ( value )

!*****************************************************************************80
!
!! F_HOFSTADTER computes the Hofstadter F sequence.
!
!  Discussion:
!
!    F(N) = 0                if N = 0
!         = N - F ( N - 1 ), otherwise.
!
!    F(N) is defined for all nonnegative integers, and turns out
!    to be equal to int ( ( N + 1 ) / 2 ).
!
!  Table:
!
!     N  F(N)
!    --  ----
!
!     0     0
!     1     1
!     2     1
!     3     2
!     4     2
!     5     3
!     6     3
!     7     4
!     8     4
!     9     5
!    10     5
!    11     6
!    12     6
!    13     7
!    14     7
!    15     8
!    16     8
!    17     9
!    18     9
!    19    10
!    20    10
!    21    11
!    22    11
!    23    12
!    24    12
!    25    13
!    26    13
!    27    14
!    28    14
!    29    15
!    30    15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Douglas Hofstadter,
!    Goedel, Escher, Bach,
!    Basic Books, 1979,
!    ISBN: 0465026567,
!    LC: QA9.8H63.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, integer ( kind = 4 ) F_HOFSTADTER, the value of the function.
!
  implicit none

  integer ( kind = 4 ) arg
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  if ( n <= 0 ) then
    value = 0
  else
    arg = n - 1
    value = n - f_hofstadter ( arg )
  end if

  return
end
subroutine fibonacci_direct ( n, f )

!*****************************************************************************80
!
!! FIBONACCI_DIRECT computes the N-th Fibonacci number directly.
!
!  Discussion:
!
!    A direct formula for the N-th Fibonacci number is:
!
!      F(N) = ( PHIP^N - PHIM^N ) / sqrt(5)
!
!    where 
!
!      PHIP = ( 1 + sqrt(5) ) / 2, 
!      PHIM = ( 1 - sqrt(5) ) / 2.
!
!  Example:
!
!     N   F
!    --  --
!     0   0
!     1   1
!     2   1
!     3   2
!     4   3
!     5   5
!     6   8
!     7  13
!     8  21
!     9  34
!    10  55
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the Fibonacci number 
!    to compute.  N should be nonnegative.
!
!    Output, integer ( kind = 4 ) F, the value of the N-th Fibonacci number.
!
  implicit none

  integer ( kind = 4 ) f
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: sqrt5 = 2.236068D+00
  real ( kind = 8 ), parameter :: phim = ( 1.0D+00 - sqrt5 ) / 2.0D+00
  real ( kind = 8 ), parameter :: phip = ( 1.0D+00 + sqrt5 ) / 2.0D+00

  if ( n < 0 ) then
    f = 0
  else
    f = nint ( ( phip**n - phim**n ) / sqrt ( 5.0D+00 ) )
  end if
 
  return
end
subroutine fibonacci_floor ( n, f, i )

!*****************************************************************************80
!
!! FIBONACCI_FLOOR returns the largest Fibonacci number less than or equal to N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the positive integer whose Fibonacci 
!    "floor" is desired.
!
!    Output, integer ( kind = 4 ) F, the largest Fibonacci number less 
!    than or equal to N.
!
!    Output, integer ( kind = 4 ) I, the index of the F.
!
  implicit none

  integer ( kind = 4 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  if ( n <= 0 ) then

    i = 0
    f = 0

  else

    i = int ( &
        log ( 0.5D+00 * real ( 2 * n + 1, kind = 8 ) * sqrt ( 5.0D+00 ) ) &
      / log ( 0.5D+00 * ( 1.0D+00 + sqrt ( 5.0D+00 ) ) ) )

    call fibonacci_direct ( i, f )

    if ( n < f ) then
      i = i - 1
      call fibonacci_direct ( i, f )
    end if

  end if

  return
end
subroutine fibonacci_recursive ( n, f )

!*****************************************************************************80
!
!! FIBONACCI_RECURSIVE computes the first N Fibonacci numbers.
!
!  Discussion:
!
!    The 'golden ratio' 
!
!      PHI = (1+sqrt(5))/2 
!
!    satisfies the algebraic equation:
!
!      X*X-X-1=0
!
!    which is often written as:
!
!       X        1
!      --- =  ------
!       1      X - 1
!
!    expressing the fact that a rectangle, whose sides are in proportion X:1,
!    is similar to the rotated rectangle after a square of side 1 is removed.
!
!      <----X---->
! 
!      +-----*---*
!      |     |   |  1
!      |     |   | 
!      +-----*---+
!      <--1-><X-1>
!
!    A direct formula for the N-th Fibonacci number can be found.
!
!    Let
!
!      PHIP = ( 1 + sqrt(5) ) / 2
!      PHIM = ( 1 - sqrt(5) ) / 2
!
!    Then
!
!      F(N) = ( PHIP^N + PHIM^N ) / sqrt(5)
!
!    Moreover, F(N) can be computed by computing PHIP^N / sqrt(5) and rounding
!    to the nearest whole number.
!
!    The function 
!
!      F(X) = X / ( 1 - X - X^2 )
!
!    has a power series whose coefficients are the Fibonacci numbers:
!
!      F(X) = 0 + 1*X + 1*X^2 + 2*X^3 + 3*X^4 + 5*X^5+...
!
!  First terms:
!
!      0
!      1
!      1
!      2
!      3
!      5
!      8
!     13
!     21
!     34
!     55
!     89
!    144
!
!    The 40th number is                  102,334,155.
!    The 50th number is               12,586,269,025.
!    The 100th number is 354,224,848,179,261,915,075.
!
!  Recursion:
!
!    F(0) = 0
!    F(1) = 1
!
!    F(N) = F(N-1) + F(N-2)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest Fibonacci number to compute.
!
!    Output, integer ( kind = 4 ) F(N), the first N Fibonacci numbers.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) f(n)
  integer ( kind = 4 ) i

  if ( n <= 0 ) then
    return
  end if

  f(1) = 1

  if ( n <= 1 ) then
    return
  end if

  f(2) = 1

  do i = 3, n
    f(i) = f(i-1) + f(i-2)
  end do
 
  return
end
recursive function g_hofstadter ( n ) result ( value )

!*****************************************************************************80
!
!! G_HOFSTADTER computes the Hofstadter G sequence.
!
!  Discussion:
!
!    G(N) = 0                      if N = 0
!         = N - G ( G ( N - 1 ) ), otherwise.
!
!    G(N) is defined for all nonnegative integers.
!
!    The value of G(N) turns out to be related to the Zeckendorf 
!    representation of N as a sum of non-consecutive Fibonacci numbers.  
!    To compute G(N), determine the Zeckendorf representation:
!
!      N = sum ( 1 <= I <= M ) F(I)
!
!    and reduce the index of each Fibonacci number by 1:
!
!      G(N) = sum ( 1 <= I <= M ) F(I-1)
! 
!    However, this is NOT how the computation is done in this routine.
!    Instead, a straightforward recursive function call is defined
!    to correspond to the definition of the mathematical function.
!
!  Table:
!
!     N  G(N)  Zeckendorf   Decremented
!    --  ----  ----------   -----------
!
!     1   1    1            1
!     2   1    2            1             
!     3   2    3            2
!     4   3    3 + 1        2 + 1
!     5   3    5            3
!     6   4    5 + 1        3 + 1
!     7   4    5 + 2        3 + 1
!     8   5    8            5
!     9   6    8 + 1        5 + 1
!    10   6    8 + 2        5 + 1
!    11   7    8 + 3        5 + 2
!    12   8    8 + 3 + 1    5 + 2 + 1
!    13   8    13           8
!    14   9    13 + 1       8 + 1
!    15   9    13 + 2       8 + 1
!    16  10    13 + 3       8 + 2
!    17  11    13 + 3 + 1   8 + 2 + 1
!    18  11    13 + 5       8 + 3
!    19  12    13 + 5 + 1   8 + 3 + 1
!    20  12    13 + 5 + 2   8 + 3 + 1
!    21  13    21           13
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Douglas Hofstadter,
!    Goedel, Escher, Bach,
!    Basic Books, 1979,
!    ISBN: 0465026567,
!    LC: QA9.8H63.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, integer ( kind = 4 ) G_HOFSTADTER, the value of the function.
!
  implicit none

  integer ( kind = 4 ) arg
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  if ( n <= 0 ) then
    value = 0
  else
    arg = n - 1
    value = n - g_hofstadter ( g_hofstadter ( arg ) )
  end if

  return
end
subroutine gamma_log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_LOG_VALUES returns some values of the Log Gamma function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 April 2001
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and
!    N_DATA is set to the index of the test data.  On each subsequent
!    call, N_DATA is incremented and that test data is returned.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 18

  real ( kind = 8 ), save, dimension ( nmax ) :: bvec = (/ &
     1.524064183D+00,    0.7966780066D+00,   0.3982337117D+00,  &
     0.1520599127D+00,   0.000000000D+00,   -0.04987246543D+00, &
    -0.08537410945D+00, -0.1081747934D+00,  -0.1196128950D+00,  & 
    -0.1207822040D+00,  -0.1125917658D+00,  -0.09580771625D+00, &
    -0.07108385116D+00, -0.03898428380D+00,  0.000000000D+00,   &
    12.80182743D+00,    39.33988571D+00,    71.25704193D+00 /)
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( nmax ) :: xvec = (/ &
    0.2D+00,  0.4D+00,  0.6D+00,  0.8D+00, &
    1.0D+00,  1.1D+00,  1.2D+00,  1.3D+00, &
    1.4D+00,  1.5D+00,  1.6D+00,  1.7D+00, &
    1.8D+00,  1.9D+00,  2.0D+00, 10.0D+00, &
   20.0D+00, 30.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = xvec(n_data)
    fx = bvec(n_data)
  end if

  return
end
subroutine gamma_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_VALUES returns some values of the Gamma function.
!
!  Discussion:
!
!    The Gamma function is defined as:
!
!      Gamma(Z) = Integral ( 0 <= t < +oo) t^(Z-1) exp(-t) dt
!
!    It satisfies the recursion:
!
!      Gamma(X+1) = X * Gamma(X)
!
!    Gamma is undefined for nonpositive integral X.
!    Gamma(0.5) = sqrt(PI)
!    For N a positive integer, Gamma(N+1) = N!, the standard factorial.
!
!    In Mathematica, the function can be evaluated by:
!
!      Gamma[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.3544907701811032D+01, &
    -0.1005871979644108D+03, &
     0.9943258511915060D+02, &
     0.9513507698668732D+01, &
     0.4590843711998803D+01, &
     0.2218159543757688D+01, &
     0.1772453850905516D+01, &
     0.1489192248812817D+01, &
     0.1164229713725303D+01, &
     0.1000000000000000D+01, &
     0.9513507698668732D+00, &
     0.9181687423997606D+00, &
     0.8974706963062772D+00, &
     0.8872638175030753D+00, &
     0.8862269254527580D+00, &
     0.8935153492876903D+00, &
     0.9086387328532904D+00, &
     0.9313837709802427D+00, &
     0.9617658319073874D+00, &
     0.1000000000000000D+01, &
     0.2000000000000000D+01, &
     0.6000000000000000D+01, &
     0.3628800000000000D+06, &
     0.1216451004088320D+18, &
     0.8841761993739702D+31 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.50D+00, &
    -0.01D+00, &
     0.01D+00, &
     0.10D+00, &
     0.20D+00, &
     0.40D+00, &
     0.50D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /) 

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gegenbauer_poly ( n, alpha, x, cx )

!*****************************************************************************80
!
!! GEGENBAUER_POLY computes the Gegenbauer polynomials C(I,ALPHA,X).
!
!  Discussion:
!
!    The Gegenbauer polynomial can be evaluated in Mathematica with
!    the command 
!
!      GegenbauerC[n,m,x]
!
!    ALPHA must be greater than -0.5.
!
!    If ALPHA = 1, the Gegenbauer polynomials reduce to the Chebyshev
!    polynomials of the second kind.
!
!  Differential equation:
!
!    (1-X*X) Y'' - (2 ALPHA + 1) X Y' + N (N + 2 ALPHA) Y = 0
!
!  Recursion:
!
!    C(0,ALPHA,X) = 1,
!    C(1,ALPHA,X) = 2*ALPHA*X
!    C(N,ALPHA,X) = ( (2*N-2+2*ALPHA) * X * C(N-1,ALPHA,X) 
!                   + ( -N+2-2*ALPHA)     * C(N-2,ALPHA,X) ) / N
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 ) 
!      ( 1 - X^2 )^( ALPHA - 0.5 ) * C(N,ALPHA,X)^2 dX
!
!    = PI * 2^( 1 - 2 * ALPHA ) * Gamma ( N + 2 * ALPHA ) 
!      / ( N! * ( N + ALPHA ) * ( Gamma ( ALPHA ) )^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) ALPHA, a parameter which is part of the 
!    definition of the Gegenbauer polynomials.  It must be greater than -0.5.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials 
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 Gegenbauer
!    polynomials at the point X.  
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) cx(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( alpha <= -0.5D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEGENBAUER_POLY - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal value of ALPHA = ', alpha
    write ( *, '(a)' ) '  but ALPHA must be greater than -0.5.'
    stop
  end if

  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = 2.0D+00 * alpha * x

  do i = 2, n
    cx(i) = &
      ( ( real ( 2 * i - 2, kind = 8 ) + 2.0D+00 * alpha ) * x * cx(i-1)   &
      + ( real (   - i + 2, kind = 8 ) - 2.0D+00 * alpha )     * cx(i-2) ) &
      /   real (     i,     kind = 8 )
  end do
 
  return
end
subroutine gegenbauer_poly_values ( n_data, n, a, x, fx )

!*****************************************************************************80
!
!! GEGENBAUER_POLY_VALUES returns some values of the Gegenbauer polynomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2004
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order parameter of the function.
!
!    Output, real ( kind = 8 ) A, the real parameter of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 38

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.0D+00, &
     1.0D+00,  2.0D+00,  3.0D+00, &
     4.0D+00,  5.0D+00,  6.0D+00, &
     7.0D+00,  8.0D+00,  9.0D+00, &
    10.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00,  3.0D+00, &
     3.0D+00,  3.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.0000000000D+00,   0.2000000000D+00,  -0.4400000000D+00, &
   -0.2800000000D+00,   0.2320000000D+00,   0.3075200000D+00, &
   -0.0805760000D+00,  -0.2935168000D+00,  -0.0395648000D+00, &
    0.2459712000D+00,   0.1290720256D+00,   0.0000000000D+00, &
   -0.3600000000D+00,  -0.0800000000D+00,   0.8400000000D+00, &
    2.4000000000D+00,   4.6000000000D+00,   7.4400000000D+00, &
   10.9200000000D+00,  15.0400000000D+00,  19.8000000000D+00, &
   25.2000000000D+00,  -9.0000000000D+00,  -0.1612800000D+00, &
   -6.6729600000D+00,  -8.3750400000D+00,  -5.5267200000D+00, &
    0.0000000000D+00,   5.5267200000D+00,   8.3750400000D+00, &
    6.6729600000D+00,   0.1612800000D+00,  -9.0000000000D+00, &
  -15.4252800000D+00,  -9.6969600000D+00,  22.4409600000D+00, &
  100.8892800000D+00, 252.0000000000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10,  2, &
     2,  2,  2, &
     2,  2,  2, &
     2,  2,  2, &
     2,  5,  5, &
     5,  5,  5, &
     5,  5,  5, &
     5,  5,  5, &
     5,  5,  5, &
     5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.20D+00,  0.20D+00,  0.20D+00, &
    0.20D+00,  0.20D+00,  0.20D+00, &
    0.20D+00,  0.20D+00,  0.20D+00, &
    0.20D+00,  0.20D+00,  0.40D+00, &
    0.40D+00,  0.40D+00,  0.40D+00, &
    0.40D+00,  0.40D+00,  0.40D+00, &
    0.40D+00,  0.40D+00,  0.40D+00, &
    0.40D+00, -0.50D+00, -0.40D+00, &
   -0.30D+00, -0.20D+00, -0.10D+00, &
    0.00D+00,  0.10D+00,  0.20D+00, &
    0.30D+00,  0.40D+00,  0.50D+00, &
    0.60D+00,  0.70D+00,  0.80D+00, &
    0.90D+00,  1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0
    x = 0.0
    fx = 0.0
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gen_hermite_poly ( n, x, mu, p )

!*****************************************************************************80
!
!! GEN_HERMITE_POLY evaluates the generalized Hermite polynomials at X.
!
!  Discussion:
!
!    The generalized Hermite polynomials are orthogonal under the weight
!    function:
!
!      w(x) = |x|^(2*MU) * exp ( - x^2 )
!
!    over the interval (-oo,+oo).
!
!    When MU = 0, the generalized Hermite polynomial reduces to the standard
!    Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Theodore Chihara,
!    An Introduction to Orthogonal Polynomials,
!    Gordon and Breach, 1978,
!    ISBN: 0677041500,
!    LC: QA404.5 C44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials are 
!    to be evaluated.
!
!    Input, real ( kind = 8 ) MU, the parameter.
!    - 1 / 2 < MU.
!
!    Output, real ( kind = 8 ) P(0:N), the values of the first N+1
!    polynomials at the point X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) mu
  real ( kind = 8 ) p(0:n)
  real ( kind = 8 ) theta
  real ( kind = 8 ) x

  if ( n < 0 ) then
    return
  end if

  p(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1) = 2.0D+00 * x
 
  do i = 1, n - 1

    if ( mod ( i, 2 ) == 0 ) then
      theta = 0.0D+00
    else
      theta = 2.0D+00 * mu
    end if

    p(i+1) = 2.0D+00 * x * p(i) &
      - 2.0D+00 * ( real ( i, kind = 8 ) + theta ) * p(i-1)

  end do
 
  return
end
subroutine gen_laguerre_poly ( n, alpha, x, cx )

!*****************************************************************************80
!
!! GEN_LAGUERRE_POLY evaluates generalized Laguerre polynomials.
!
!  Differential equation:
!
!    X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0
!
!  Recursion:
!
!    L(0,ALPHA,X) = 1
!    L(1,ALPHA,X) = 1+ALPHA-X
!
!    L(N,ALPHA,X) = ( (2*N-1+ALPHA-X) * L(N-1,ALPHA,X) 
!                   - (N-1+ALPHA) * L(N-2,ALPHA,X) ) / N
!
!  Restrictions:
!
!    -1 < ALPHA
!
!  Special values:
!
!    For ALPHA = 0, the generalized Laguerre polynomial L(N,ALPHA,X)
!    is equal to the Laguerre polynomial L(N,X).
!
!    For ALPHA integral, the generalized Laguerre polynomial
!    L(N,ALPHA,X) equals the associated Laguerre polynomial L(N,ALPHA,X).
!
!  Norm:
!
!    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,ALPHA,X)^2 dX
!    = Gamma ( N + ALPHA + 1 ) / N!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2002
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
!    Input, integer ( kind = 4 ) N, the highest order function to compute.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.  -1 < ALPHA is required.
!
!    Input, real ( kind = 8 ) X, the point at which the functions are to be
!    evaluated.
!
!    Output, real ( kind = 8 ) CX(0:N), the functions of 
!    degrees 0 through N evaluated at the point X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) cx(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GEN_LAGUERRE_POLY - Fatal error!'
    write ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
    write ( *, '(a)' ) '  but ALPHA must be greater than -1.'
    stop
  end if
 
  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = 1.0D+00 + alpha - x

  do i = 2, n
    cx(i) = ( ( real ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( real (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
              / real (     i,     kind = 8 )
  end do

  return
end
function gud ( x )

!*****************************************************************************80
!
!! GUD evaluates the Gudermannian function.
!
!  Discussion:
!
!    The Gudermannian function relates the hyperbolic and trigonometric
!    functions.  For any argument X, there is a corresponding value
!    GAMMA so that
!
!      sinh(x) = tan(gamma).
!
!    The value GAMMA is called the Gudermannian of X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gudermannian.
!
!    Output, real ( kind = 8 ) GUD, the value of the Gudermannian.
!
  implicit none

  real ( kind = 8 ) gud
  real ( kind = 8 ) x

  gud = 2.0D+00 * atan ( tanh ( 0.5D+00 * x ) )

  return
end
subroutine gud_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GUD_VALUES returns some values of the Gudermannian function.
!
!  Discussion:
!
!    The Gudermannian function relates the hyperbolic and trigonometric
!    functions.  For any argument X, there is a corresponding value
!    gd(x) so that
!
!      SINH(x) = TAN(gd(x)).
!
!    This value GD is called the Gudermannian of X and symbolized
!    GD(X).  The inverse Gudermannian function is given as input a value 
!    GD and computes the corresponding value X.
!
!    GD(X) = 2 * arctan ( exp ( X ) ) - PI / 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
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
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -1.301760336D+00,  -0.8657694832D+00, 0.0000000000D+00, &
     0.09983374879D+00, 0.1986798470D+00, 0.4803810791D+00, &
     0.8657694832D+00,  1.131728345D+00,  1.301760336D+00,  &
     1.406993569D+00,   1.471304341D+00,  1.510419908D+00,  &
     1.534169144D+00 /)
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -2.0D+00, -1.0D+00,  0.0D+00, &
     0.1D+00,  0.2D+00,  0.5D+00, &
     1.0D+00,  1.5D+00,  2.0D+00, &
     2.5D+00,  3.0D+00,  3.5D+00, &
     4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
recursive function h_hofstadter ( n ) result ( value )

!*****************************************************************************80
!
!! H_HOFSTADTER computes the Hofstadter H sequence.
!
!  Discussion:
!
!    H(N) = 0                          if N = 0
!         = N - H ( H ( H ( N - 1 ) ), otherwise.
!
!    H(N) is defined for all nonnegative integers.
!
!  Example:
!
!     N  H(N)
!    --  ----
!
!     0     0
!     1     1
!     2     1
!     3     2
!     4     3
!     5     4
!     6     4
!     7     5
!     8     5
!     9     6
!    10     7
!    11     7
!    12     8
!    13     9
!    14    10
!    15    10
!    16    11
!    17    12
!    18    13
!    19    13
!    20    14
!    21    14
!    22    15
!    23    16
!    24    17
!    25    17
!    26    18
!    27    18
!    28    19
!    29    20
!    30    20
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Douglas Hofstadter,
!    Goedel, Escher, Bach,
!    Basic Books, 1979,
!    ISBN: 0465026567,
!    LC: QA9.8H63.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, integer ( kind = 4 ) H_HOFSTADTER, the value of the function.
!
  implicit none

  integer ( kind = 4 ) arg
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  if ( n <= 0 ) then
    value = 0
  else
    arg = n - 1
    value = n - h_hofstadter ( h_hofstadter ( h_hofstadter ( arg ) ) )
  end if

  return
end
subroutine hermite_poly ( n, x, cx )

!*****************************************************************************80
!
!! HERMITE_POLY evaluates the physicist's Hermite polynomials at X.
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X^2     -  2
!      8 X^3     - 12 X
!     16 X^4     - 48 X^2     + 12
!     32 X^5    - 160 X^3    + 120 X
!     64 X^6    - 480 X^4    + 720 X^2    - 120
!    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
!    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
!    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
!   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < oo ) exp ( - X^2 ) * H(N,X)^2 dX
!    = sqrt ( PI ) * 2^N * N!
!
!    H(N,X) = (-1)^N * exp ( X^2 ) * dn/dXn ( exp(-X^2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2002
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
!    Larry Andrews,
!    Special Functions of Mathematics for Engineers,
!    Second Edition, 
!    Oxford University Press, 1998.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials are 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 Hermite
!    polynomials at the point X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = 2.0D+00 * x
 
  do i = 2, n
    cx(i) = 2.0D+00 * x * cx(i-1) - 2.0D+00 * real ( i - 1, kind = 8 ) * cx(i-2)
  end do
 
  return
end
subroutine hermite_poly_coef ( n, c )

!*****************************************************************************80
!
!! HERMITE_POLY_COEF: coefficients of the physicist's Hermite polynomial H(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     2
!     2     -2     0      4
!     3      0   -12      0      8
!     4     12     0    -48      0      16
!     5      0   120      0   -160       0    32
!     6   -120     0    720      0    -480     0     64
!     7      0 -1680      0   3360       0 -1344      0   128
!     8   1680     0 -13440      0   13440     0  -3584     0    256
!     9      0 30240      0 -80640       0 48384      0 -9216      0 512
!    10 -30240     0 302400      0 -403200     0 161280     0 -23040   0 1024 
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2003
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the Hermite
!    polynomials.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 2.0D+00
 
  do i = 2, n
    c(i,0)     =  -2.0D+00 * real ( i - 1, kind = 8 ) * c(i-2,0)
    c(i,1:i-2) =   2.0D+00                            * c(i-1,0:i-3)  &
                  -2.0D+00 * real ( i - 1, kind = 8 ) * c(i-2,1:i-2)
    c(i,  i-1) =   2.0D+00                            * c(i-1,  i-2)
    c(i,  i  ) =   2.0D+00                            * c(i-1,  i-1)
  end do
 
  return
end
subroutine hermite_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! HERMITE_POLY_VALUES returns some values of the Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( nmax ) :: fx_vec = (/ &
     1.0D+00,            10.0D+00,           98.0D+00, &
     940.0D+00,          8812.0D+00,         80600.0D+00, &
     717880.0D+00,       6211600.0D+00,      520656800.0D+00, &
     421271200D+00,      3275529760.0D+00,   24329873600.0D+00, &
     171237081280.0D+00, 41.0D+00,          -8.0D+00, &
     3816.0D+00,         3041200.0D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( nmax ) :: x_vec = (/ &
    5.0D+00,  5.0D+00,  5.0D+00, &
    5.0D+00,  5.0D+00,  5.0D+00, &
    5.0D+00,  5.0D+00,  5.0D+00, &
    5.0D+00,  5.0D+00,  5.0D+00, &
    5.0D+00,  0.5D+00,  1.0D+00, &
    3.0D+00, 10.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hyper_2f1_values ( n_data, a, b, c, x, fx )

!*****************************************************************************80
!
!! HYPER_2F1_VALUES returns some values of the hypergeometric function 2F1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = Hypergeometric2F1 [ a, b, c, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 September 2007
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
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
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
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, C, X, the parameters.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 24

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00 /)
  real ( kind = 8 ) c
  real ( kind = 8 ), save, dimension ( n_max ) :: c_vec = (/ &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.72356129348997784913D+00, &
    0.97911109345277961340D+00, &
    1.0216578140088564160D+00, &
    1.4051563200112126405D+00, &
    0.46961431639821611095D+00, &
    0.95296194977446325454D+00, &
    1.0512814213947987916D+00, &
    2.3999062904777858999D+00, &
    0.29106095928414718320D+00, &
    0.92536967910373175753D+00, &
    1.0865504094806997287D+00, &
    5.7381565526189046578D+00, &
    15090.669748704606754D+00, &
   -104.31170067364349677D+00, &
    21.175050707768812938D+00, &
    4.1946915819031922850D+00, &
    1.0170777974048815592D+10, &
   -24708.635322489155868D+00, &
    1372.2304548384989560D+00, &
    58.092728706394652211D+00, &
    5.8682087615124176162D+18, &
   -4.4635010147295996680D+08, &
    5.3835057561295731310D+06, &
    20396.913776019659426D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if
 
  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    c = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    c = c_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K).
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
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
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end
subroutine i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

!*****************************************************************************80
!
!! I4_FACTOR factors an integer into prime factors.
!
!  Discussion:
!
!    This routine tries to decompose an integer into prime factors,
!    but the decomposition may be left incomplete either because
!    there are too many factors for the allowed space, or because
!    the factors are too large for the table of primes to handle.
!
!    The form of the factorization therefore includes an unresolved
!    factor "NLEFT":
!
!      N = NLEFT * Product ( I = 1 to NFACTOR ) FACTOR(I)^POWER(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be factored.  N may be 
!    positive, negative, or 0.
!
!    Input, integer ( kind = 4 ) MAXFACTOR, the maximum number of prime 
!    factors for which storage has been allocated.
!
!    Output, integer ( kind = 4 ) NFACTOR, the number of prime factors of 
!    N discovered by the routine.
!
!    Output, integer ( kind = 4 ) FACTOR(MAXFACTOR), the prime factors of N.
!
!    Output, integer ( kind = 4 ) POWER(MAXFACTOR).  POWER(I) is the power of
!    the FACTOR(I) in the representation of N.
!
!    Output, integer ( kind = 4 ) NLEFT, the factor of N that the routine 
!    could not divide out.  If NLEFT is 1, then N has been completely factored.
!    Otherwise, NLEFT represents factors of N involving large primes.
!
  implicit none

  integer ( kind = 4 ) maxfactor

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) maxprime
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) p
  integer ( kind = 4 ) power(maxfactor)
  integer ( kind = 4 ) prime

  nfactor = 0

  factor(1:maxfactor) = 0
  power(1:maxfactor) = 0

  nleft = n

  if ( n == 0 ) then
    return
  end if

  if ( abs ( n ) == 1 ) then
    nfactor = 1
    factor(1) = 1
    power(1) = 1
    return
  end if
!
!  Find out how many primes we stored.
!
  maxprime = prime ( -1 )
!
!  Try dividing the remainder by each prime.
!
  do i = 1, maxprime

    p = prime ( i )

    if ( mod ( abs ( nleft ), p ) == 0 ) then

      if ( nfactor < maxfactor ) then

        nfactor = nfactor + 1
        factor(nfactor) = p

        do

          power(nfactor) = power(nfactor) + 1
          nleft = nleft / p

          if ( mod ( abs ( nleft ), p ) /= 0 ) then
            exit
          end if

        end do

        if ( abs ( nleft ) == 1 ) then
          exit
        end if

      end if

    end if

  end do

  return
end
function i4_factorial ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!    0 <= N <= 13 is required.
!
!    Output, integer ( kind = 4 ) I4_FACTORIAL, the factorial of N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) n

  i4_factorial = 1

  if ( 13 < n ) then
    i4_factorial = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_FACTORIAL - Fatal error!'
    write ( *, '(a)' ) '  I4_FACTORIAL(N) cannot be computed as an integer'
    write ( *, '(a)' ) '  for 13 < N.'
    write ( *, '(a,i8)' ) '  Input value N = ', n
    stop
  end if

  do i = 1, n
    i4_factorial = i4_factorial * i
  end do

  return
end
subroutine i4_factorial_values ( n_data, n, fn )

!*****************************************************************************80
!
!! I4_FACTORIAL_VALUES returns values of the factorial function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 May 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA 
!    is set to the index of the test data.  On each subsequent call, N_DATA is
!    incremented and that test data is returned.  When there is no more
!    test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, integer ( kind = 4 ) FN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 13

  integer ( kind = 4 ), save, dimension ( nmax ) :: fnvec = (/ &
    1, 1, 2, 6, &
    24, 120, 720, 5040, &
    40320, 362880, 3628800, 39916800, &
    479001600 /)
  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: nvec = (/ &
     0,  1,  2,  3, &
     4,  5,  6,  7, &
     8,  9, 10, 11, &
    12 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    fn = 0
  else
    n = nvec(n_data)
    fn = fnvec(n_data)
  end if

  return
end
function i4_factorial2 ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    The formula is:
!
!      FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                      = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    Factorial2(N)
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
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial 
!    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
!
!    Output, integer ( kind = 4 ) I4_FACTORIAL2, the value of the function..
!
  implicit none

  integer ( kind = 4 ) i4_factorial2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_copy

  if ( n < 1 ) then
    i4_factorial2 = 1
    return
  end if

  n_copy = n
  i4_factorial2 = 1

  do while ( 1 < n_copy )
    i4_factorial2 = i4_factorial2 * n_copy
    n_copy = n_copy - 2
  end do

  return
end
subroutine i4_factorial2_values ( n_data, n, fn )

!*****************************************************************************80
!
!! I4_FACTORIAL2_VALUES returns values of the double factorial function.
!
!  Discussion:
!
!    The formula is:
!
!      FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                      = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    Factorial2(N)
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
!    14 August 2004
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
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, page 16.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, integer ( kind = 4 ) FN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  integer ( kind = 4 ), save, dimension ( n_max ) :: fn_vec = (/ &
        1, &
        1,     2,      3,      8,      15, &
       48,   105,    384,    945,    3840, &
    10395, 46080, 135135, 645120, 2027025 /)
  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, &
     1,  2,  3,  4,  5, &
     6,  7,  8,  9, 10, &
    11, 12, 13, 14, 15 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    fn = 0
  else
    n = n_vec(n_data)
    fn = fn_vec(n_data)
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
!    On an IEEE 32 bit machine, I4_HUGE should be 2^31 - 1, and its
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

  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
function i4_is_prime ( n )

!*****************************************************************************80
!
!! I4_IS_PRIME reports whether an integer is prime.
!
!  Discussion:
!
!    A simple, unoptimized sieve of Erasthosthenes is used to
!    check whether N can be divided by any integer between 2
!    and SQRT(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be tested.
!
!    Output, logical I4_IS_PRIME, is TRUE if N is prime, and FALSE
!    otherwise.  Note that negative numbers and 0 are not
!    considered prime.
!
  implicit none

  integer ( kind = 4 ) i
  logical i4_is_prime
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhi

  if ( n <= 0 ) then
    i4_is_prime = .false.
    return
  end if

  if ( n <= 3 ) then
    i4_is_prime = .true.
    return
  end if

  nhi = int ( sqrt ( real ( n, kind = 8 ) ) )

  do i = 2, nhi
    if ( mod ( n, i ) == 0 ) then
      i4_is_prime = .false.
      return
    end if
  end do

  i4_is_prime = .true.

  return
end
function i4_is_triangular ( i )

!*****************************************************************************80
!
!! I4_IS_TRIANGULAR determines whether an integer is triangular.
!
!  Discussion:
!
!    The N-th triangular number is equal to the sum of the first
!    N integers.
!
!  First Values:
!
!    Index  Value
!     0      0
!     1      1
!     2      3
!     3      6
!     4     10
!     5     15
!     6     21
!     7     28
!     8     36
!     9     45
!    10     55
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
!    Input, integer ( kind = 4 ) I, the integer to be checked.
!
!    Output, logical I4_IS_TRIANGULAR, is TRUE if I is triangular.
!
  implicit none

  integer ( kind = 4 ) i
  logical i4_is_triangular
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  if ( i < 0 ) then

    i4_is_triangular = .false.

  else if ( i == 0 ) then

    i4_is_triangular = .true.

  else

    call i4_to_triangle ( i, j, k )

    if ( j == k ) then
      i4_is_triangular = .true.
    else
      i4_is_triangular = .false.
    end if

  end if

  return
end
subroutine i4_partition_distinct_count ( n, q )

!*****************************************************************************80
!
!! I4_PARTITION_DISTINCT_COUNT returns any value of Q(N).
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  The number of partitions of N is symbolized
!    by P(N).  Thus, the number 5 has P(N) = 7, because it has the 
!    following partitions:
!
!    5 = 5
!      = 4 + 1 
!      = 3 + 2 
!      = 3 + 1 + 1 
!      = 2 + 2 + 1 
!      = 2 + 1 + 1 + 1 
!      = 1 + 1 + 1 + 1 + 1
!
!    However, if we require that each member of the partition
!    be distinct, we are computing something symbolized by Q(N).
!    The number 5 has Q(N) = 3, because it has the following partitions 
!    into distinct parts:
!
!    5 = 5
!      = 4 + 1 
!      = 3 + 2 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2003
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
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!
!    Output, integer ( kind = 4 ) Q, the number of partitions of the integer
!    into distinct parts.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(0:n)
  integer ( kind = 4 ) i
  logical i4_is_triangular
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k_sign
  integer ( kind = 4 ) q

  c(0) = 1

  do i = 1, n

    if ( i4_is_triangular ( i ) ) then
      c(i) = 1
    else
      c(i) = 0
    end if

    k = 0
    k_sign = -1

    do 

      k = k + 1
      k_sign = - k_sign
      k2 = k * ( 3 * k + 1 )

      if ( i < k2 ) then
        exit
      end if

      c(i) = c(i) + k_sign * c(i-k2)

    end do

    k = 0
    k_sign = -1

    do 

      k = k + 1
      k_sign = - k_sign
      k2 = k * ( 3 * k - 1 )

      if ( i < k2 ) then
        exit
      end if

      c(i) = c(i) + k_sign * c(i-k2)

    end do

  end do

  q = c(n)

  return
end
function i4_pochhammer ( i, j )

!*****************************************************************************80
!
!! I4_POCHHAMMER returns the value of ( I * (I+1) * ... * (J-1) * J ).
!
!  Discussion:
!
!    Pochhammer's symbol (A)_N is the value
!
!      (A)_N = Gamma ( A + N ) / Gamma ( A )
!
!    or, for integer arguments,
!
!      (I)_N = I * ( I + 1 ) * ... * ( I + N - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, values that define the product.
!
!    Output, integer ( kind = 4 ) I4_POCHHAMMER, the value of the product.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_pochhammer
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  i4_pochhammer = 1
  do k = i, j
    i4_pochhammer = i4_pochhammer * k
  end do

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4_to_triangle ( k, i, j )

!*****************************************************************************80
!
!! I4_TO_TRIANGLE converts an integer to triangular coordinates.
!
!  Discussion:
!
!    Triangular coordinates are handy when storing a naturally triangular
!    array (such as the lower half of a matrix) in a linear array.
!
!    Thus, for example, we might consider storing 
!
!    (1,1)
!    (2,1) (2,2)
!    (3,1) (3,2) (3,3)
!    (4,1) (4,2) (4,3) (4,4)
!
!    as the linear array
!
!    (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1) (4,2) (4,3) (4,4)    
!
!    Here, the quantities in parenthesis represent the natural row and
!    column indices of a single number when stored in a rectangular array.
!
!    In this routine, we are given the location K of an item in the 
!    linear array, and wish to determine the row I and column J
!    of the item when stored in the triangular array.
! 
!  Example:
!
!     K  I  J
!
!     0  0  0
!     1  1  1
!     2  2  1
!     3  2  2
!     4  3  1
!     5  3  2
!     6  3  3
!     7  4  1
!     8  4  2
!     9  4  3
!    10  4  4
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
!    Input, integer ( kind = 4 ) K, the linear index of the (I,J) element, 
!    which must be nonnegative.
!
!    Output, integer ( kind = 4 ) I, J, the row and column indices.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  if ( k < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_TRIANGLE - Fatal error!'
    write ( *, '(a)' ) '  K < 0.'
    write ( *, '(a,i8)' ) '  K = ', k
    stop

  else if ( k == 0 ) then

    i = 0
    j = 0
    return

  end if

  i = int ( sqrt ( real ( 2 * k, kind = 8 ) ) )

  if ( i * i + i < 2 * k ) then
    i = i + 1
  end if

  j = k - ( i * ( i - 1 ) ) / 2

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 8 )  ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8)' ) j
    end do

    write ( *, '(''  Col '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine jacobi_poly ( n, alpha, beta, x, cx )

!*****************************************************************************80
!
!! JACOBI_POLY evaluates the Jacobi polynomials at X.
!
!  Differential equation:
!
!    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
!
!  Recursion:
!
!    P(0,ALPHA,BETA,X) = 1,
!
!    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
!
!    P(N,ALPHA,BETA,X)  = 
!      ( 
!        (2*N+ALPHA+BETA-1) 
!        * ((ALPHA^2-BETA^2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X) 
!        * P(N-1,ALPHA,BETA,X)
!        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
!      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
!
!  Restrictions:
!
!    -1 < ALPHA
!    -1 < BETA
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA 
!      * P(N,ALPHA,BETA,X)^2 dX 
!    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
!      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
!
!  Special values:
!
!    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2002
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.  
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) ALPHA, one of the parameters defining the Jacobi
!    polynomials, ALPHA must be greater than -1.
!
!    Input, real ( kind = 8 ) BETA, the second parameter defining the Jacobi
!    polynomials, BETA must be greater than -1.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials are 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 Jacobi
!    polynomials at the point X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) cx(0:n)
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) c4
  integer ( kind = 4 ) i
  real ( kind = 8 ) r_i
  real ( kind = 8 ) x

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of ALPHA = ', alpha
    write ( *, '(a)' ) '  But ALPHA must be greater than -1.'
    stop
  end if
 
  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of BETA = ', beta
    write ( *, '(a)' ) '  But BETA must be greater than -1.'
    stop
  end if
  
  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = ( 1.0D+00 + 0.5D+00 * ( alpha + beta ) ) * x &
    + 0.5D+00 * ( alpha - beta )
 
  do i = 2, n

    r_i = real ( i, kind = 8 ) 

    c1 = 2.0D+00 * r_i * ( r_i + alpha + beta ) &
      * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

    c2 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) &
      * ( 2.0D+00 * r_i  + alpha + beta ) &
      * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

    c3 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) &
      * ( alpha + beta ) * ( alpha - beta )

    c4 = - 2.0D+00 * ( r_i - 1.0D+00 + alpha ) &
      * ( r_i - 1.0D+00 + beta )  * ( 2.0D+00 * r_i + alpha + beta )

    cx(i) = ( ( c3 + c2 * x ) * cx(i-1) + c4 * cx(i-2) ) / c1

  end do

  return
end
subroutine jacobi_poly_values ( n_data, n, a, b, x, fx )

!*****************************************************************************80
!
!! JACOBI_POLY_VALUES returns some values of the Jacobi polynomial.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      JacobiP[ n, a, b, x ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2012
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Output, real ( kind = 8 ) A, B, parameters of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 26

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, &
     3.0D+00, 4.0D+00, 5.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, &
    3.0D+00, 4.0D+00, 5.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.2500000000000000D+00, &
    -0.3750000000000000D+00, &
    -0.4843750000000000D+00, &
    -0.1328125000000000D+00, &
     0.2753906250000000D+00, &
    -0.1640625000000000D+00, &
    -0.1174804687500000D+01, &
    -0.2361328125000000D+01, &
    -0.2616210937500000D+01, &
     0.1171875000000000D+00, &
     0.4218750000000000D+00, &
     0.5048828125000000D+00, &
     0.5097656250000000D+00, &
     0.4306640625000000D+00, &
    -0.6000000000000000D+01, &
     0.3862000000000000D-01, &
     0.8118400000000000D+00, &
     0.3666000000000000D-01, &
    -0.4851200000000000D+00, &
    -0.3125000000000000D+00, &
     0.1891200000000000D+00, &
     0.4023400000000000D+00, &
     0.1216000000000000D-01, &
    -0.4396200000000000D+00, &
     0.1000000000000000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, 1, 2, 3, &
     4, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
    -1.0D+00, &
    -0.8D+00, &
    -0.6D+00, &
    -0.4D+00, &
    -0.2D+00, &
     0.0D+00, &
     0.2D+00, &
     0.4D+00, &
     0.6D+00, &
     0.8D+00, &
     1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine jacobi_symbol ( q, p, j )

!*****************************************************************************80
!
!! JACOBI_SYMBOL evaluates the Jacobi symbol (Q/P).
!
!  Discussion:
!
!    If P is prime, then
!
!      Jacobi Symbol (Q/P) = Legendre Symbol (Q/P)
!
!    Else 
!
!      let P have the prime factorization
!
!        P = Product ( 1 <= I <= N ) P(I)^E(I)
!
!      Jacobi Symbol (Q/P) =
!
!        Product ( 1 <= I <= N ) Legendre Symbol (Q/P(I))^E(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 86-87.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Q, an integer whose Jacobi symbol with
!    respect to P is desired.
!
!    Input, integer ( kind = 4 ) P, the number with respect to which the Jacobi
!    symbol of Q is desired.  P should be 2 or greater.
!
!    Output, integer ( kind = 4 ) J, the Jacobi symbol (Q/P).
!    Ordinarily, J will be -1, 0 or 1.
!    -2, not enough factorization space.
!    -3, an error during Legendre symbol calculation.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxfactor = 20

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) p
  integer ( kind = 4 ) power(maxfactor)
  integer ( kind = 4 ) pp
  integer ( kind = 4 ) q
  integer ( kind = 4 ) qq
!
!  P must be greater than 1.
!
  if ( p <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_SYMBOL - Fatal error!'
    write ( *, '(a)' ) '  P must be greater than 1.'
    l = -2
    stop
  end if
!
!  Decompose P into factors of prime powers.
!
  call i4_factor ( p, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_SYMBOL - Fatal error!'
    write ( *, '(a)' ) '  Not enough factorization space.'
    j = -2
    stop
  end if
!
!  Force Q to be nonnegative.
!
  qq = q

  do while ( qq < 0 )
    qq = qq + p
  end do
!
!  For each prime factor, compute the Legendre symbol, and
!  multiply the Jacobi symbol by the appropriate factor.
!
  j = 1
  do i = 1, nfactor
    pp = factor(i)
    call legendre_symbol ( qq, pp, l )
    if ( l < -1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'JACOBI_SYMBOL - Fatal error!'
      write ( *, '(a)' ) '  Error during Legendre symbol calculation.'
      j = -3
      stop
    end if
    j = j * l**power(i)
  end do

  return
end
subroutine krawtchouk ( n, p, x, m, v )

!*****************************************************************************80
!
!! KRAWTCHOUK evaluates the Krawtchouk polynomials at X.
!
!  Discussion:
!
!    The polynomial has a parameter P, which must be striclty between
!    0 and 1, and a parameter M which must be a nonnegative integer.
!
!    The Krawtchouk polynomial of order N, with parameters P and M,
!    evaluated at X, may be written K(N,P,X,M).
!
!    The first two terms are:
!
!      K(0,P,X,M) = 1
!      K(1,P,X,M) = X - P * M
!
!    and the recursion, for fixed P and M is
!
!                             ( N + 1 ) * K(N+1,P,X,M) =
!        ( X - ( N + P * ( M - 2 * N))) * K(N,  P,X,M)
!       - ( M - N + 1 ) * P * ( 1 - P ) * K(N-1,P,X,M)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Walter Gautschi,
!    Orthogonal Polynomials: Computation and Approximation,
!    Oxford, 2004,
!    ISBN: 0-19-850672-4,
!    LC: QA404.5 G3555.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    0 <= N.
!
!    Input, real ( kind = 8 ) P, the parameter.  0 < P < 1.
!
!    Input, real ( kind = 8 ) X, the evaluation parameter.
!
!    Input, integer ( kind = 4 ) M, the parameter.  0 <= M.
!
!    Output, real ( kind = 8 ) V(0:N), the values of the Krawtchouk polynomials
!    of orders 0 through N at X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) p
  real ( kind = 8 ) x
  real ( kind = 8 ) v(0:n)

  if ( n < 0 ) then
    write ( * , '(a)' ) ' '
    write ( * , '(a)' ) 'KRAWTCHOUK - Fatal error!'
    write ( * , '(a)' ) '  0 <= N is required.'
    stop
  end if

  if ( p <= 0.0D+00 .or. 1.0D+00 <= p ) then
    write ( * , '(a)' ) ' '
    write ( * , '(a)' ) 'KRAWTCHOUK - Fatal error!'
    write ( * , '(a)' ) '  0 < P < 1 is required.'
    stop
  end if

  if ( m < 0 ) then
    write ( * , '(a)' ) ' '
    write ( * , '(a)' ) 'KRAWTCHOUK - Fatal error!'
    write ( * , '(a)' ) '  0 <= M is required.'
    stop
  end if

  v(0) = 1.0D+00

  if ( 1 <= n ) then
    v(1) = x - p * real ( m, kind = 8 )
  end if

  do i = 1, n - 1
    v(i+1) = ( &
      ( x - ( real ( i, kind = 8 ) + p * real ( m - 2 * i, kind = 8 ) ) ) &
        * v(i) &
      - real ( m - i + 1, kind = 8 ) * p * ( 1.0D+00 - p ) * v(i-1)  &
      ) / real ( i + 1, kind = 8 )
  end do

  return
end
subroutine laguerre_associated ( n, m, x, cx )

!*****************************************************************************80
!
!! LAGUERRE_ASSOCIATED evaluates associated Laguerre polynomials L(N,M,X).
!
!  Differential equation:
!
!    X Y'' + (M+1-X) Y' + (N-M) Y = 0
!
!  First terms:
!
!    M = 0
!
!    L(0,0,X) =   1
!    L(1,0,X) =  -X   +  1
!    L(2,0,X) =   X^2 -  4 X   +  2
!    L(3,0,X) =  -X^3 +  9 X^2 -  18 X   +    6
!    L(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +     24
!    L(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x   +  120
!    L(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
!
!    M = 1
!
!    L(0,1,X) =    0
!    L(1,1,X) =   -1,
!    L(2,1,X) =    2 X - 4,
!    L(3,1,X) =   -3 X^2 + 18 X - 18,
!    L(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
!
!    M = 2
!
!    L(0,2,X) =    0
!    L(1,2,X) =    0,
!    L(2,2,X) =    2,
!    L(3,2,X) =   -6 X + 18,
!    L(4,2,X) =   12 X^2 - 96 X + 144
!
!    M = 3
!
!    L(0,3,X) =    0
!    L(1,3,X) =    0,
!    L(2,3,X) =    0,
!    L(3,3,X) =   -6,
!    L(4,3,X) =   24 X - 96
!
!    M = 4
!
!    L(0,4,X) =    0
!    L(1,4,X) =    0
!    L(2,4,X) =    0
!    L(3,4,X) =    0
!    L(4,4,X) =   24
!
!  Recursion:
!
!    if N = 0:
!
!      L(N,M,X)   = 0 
!
!    if N = 1:
!
!      L(N,M,X)   = (M+1-X)
!
!    if 2 <= N:
!
!      L(N,M,X)   = ( (M+2*N-1-X) * L(N-1,M,X) 
!                  +   (1-M-N)     * L(N-2,M,X) ) / N
!
!  Special values:
!
!    For M = 0, the associated Laguerre polynomials L(N,M,X) are equal 
!    to the Laguerre polynomials L(N,X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 February 2003
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer ( kind = 4 ) M, the parameter.  M must be nonnegative.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials are 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CX(0:N), the associated Laguerre polynomials of 
!    degrees 0 through N evaluated at the point X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) x

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop
  end if

  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = real ( m + 1, kind = 8 ) - x

  do i = 2, n
    cx(i) = ( ( real (   m + 2 * i - 1, kind = 8 ) - x ) * cx(i-1)   &
              + real ( - m     - i + 1, kind = 8       ) * cx(i-2) ) &
              / real (           i,     kind = 8 )
  end do

  return
end
subroutine laguerre_poly ( n, x, cx )

!*****************************************************************************80
!
!! LAGUERRE_POLY evaluates the Laguerre polynomials at X.
!
!  Differential equation:
!
!    X * Y'' + (1-X) * Y' + N * Y = 0
!
!  First terms:
!
!      1
!     -X    +  1
!   (  X^2 -  4 X     +  2 ) / 2
!   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
!   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
!   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
!   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
!   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3 
!      + 52920 X^2 - 35280 X + 5040 ) / 5040
!
!  Recursion:
!
!    L(0,X) = 1,
!    L(1,X) = 1-X,
!    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
!
!  Orthogonality:
!
!    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX
!    = 0 if N /= M
!    = 1 if N == M
!
!  Special values:
!
!    L(N,0) = 1.
!
!  Relations:
!
!    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2003
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials are 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CX(0:N), the Laguerre polynomials of 
!    degree 0 through N evaluated at the point X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = 1.0D+00 - x
 
  do i = 2, n

    cx(i) = ( ( real ( 2 * i - 1, kind = 8 ) - x ) * cx(i-1)   &
              - real (     i - 1, kind = 8 )       * cx(i-2) ) &
              / real (     i,     kind = 8 )

  end do

  return
end
subroutine laguerre_poly_coef ( n, c )

!*****************************************************************************80
!
!! LAGUERRE_POLY_COEF evaluates the Laguerre polynomial coefficients.
!
!  First terms:
!
!    0: 1
!    1: 1  -1
!    2: 1  -2  1/2
!    3: 1  -3  3/2  1/6
!    4: 1  -4  4   -2/3  1/24
!    5: 1  -5  5   -5/3  5/24  -1/120
!
!  Recursion:
!
!    L(0) = ( 1,  0, 0, ..., 0 )
!    L(1) = ( 1, -1, 0, ..., 0 )
!    L(N) = (2*N-1-X) * L(N-1) - (N-1) * L(N-2) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2003
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
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the
!    Laguerre polynomials of degree 0 through N.  Each polynomial 
!   is stored as a row.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0:n,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = -1.0D+00
 
  do i = 2, n

    c(i,1:n) = ( &
        real ( 2 * i - 1, kind = 8 ) * c(i-1,1:n)     &
      + real (   - i + 1, kind = 8 ) * c(i-2,1:n)     &
      -                                c(i-1,0:n-1) ) &
      / real (     i,     kind = 8 )

  end do

  return
end
subroutine laguerre_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_VALUES: some values of the Laguerre polynomial L(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LaguerreL[n,x]
!
!  Differential equation:
!
!    X * Y'' + (1-X) * Y' + N * Y = 0
!
!  First terms:
!
!      1
!     -X    +  1
!   (  X^2 -  4 X     +  2 ) / 2
!   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
!   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
!   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
!   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
!   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3 + 52920 X^2 - 35280 X 
!     + 5040 ) / 5040
!
!  Recursion:
!
!    L(0,X) = 1,
!    L(1,X) = 1-X,
!    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
!
!  Orthogonality:
!
!    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX
!    = 0 if N /= M
!    = 1 if N == M
!
!  Special values:
!
!    L(N,0) = 1.
!
!  Relations:
!
!    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2004
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 17

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.0000000000000000D+00, &
    -0.5000000000000000D+00, &
    -0.6666666666666667D+00, &
    -0.6250000000000000D+00, &
    -0.4666666666666667D+00, &
    -0.2569444444444444D+00, &
    -0.4047619047619048D-01, &
     0.1539930555555556D+00, &
     0.3097442680776014D+00, &
     0.4189459325396825D+00, &
     0.4801341790925124D+00, &
     0.4962122235082305D+00, &
    -0.4455729166666667D+00, &
     0.8500000000000000D+00, &
    -0.3166666666666667D+01, &
     0.3433333333333333D+02 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    0.5D+00, &
    3.0D+00, &
    5.0D+00, &
    1.0D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function lambert_w ( x )

!*****************************************************************************80
!
!! LAMBERT_W estimates the Lambert W function.
!
!  Discussion:
!
!    The function W(X) is defined implicitly by:
!
!      W(X) * e^W(X) = X
!
!    The function is also known as the "Omega" function.
!
!    In Mathematica, the function can be evaluated by:
!
!      W = ProductLog [ X ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Corless, Gaston Gonnet, David Hare, David Jeffrey, Donald Knuth,
!    On the Lambert W Function,
!    Advances in Computational Mathematics,
!    Volume 5, 1996, pages 329-359.
!
!    Brian Hayes,
!    "Why W?",
!    The American Scientist,
!    Volume 93, March-April 2005, pages 104-108.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) LAMBERT_W, an approximation to the 
!    Lambert W function.
!
  implicit none

  real ( kind = 8 ) lambert_w
  real ( kind = 8 ) lambert_w_crude
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  real ( kind = 8 ), parameter :: tol = 1.0D-10
  real ( kind = 8 ) w
  real ( kind = 8 ) x

  w = lambert_w_crude ( x )
  it = 0

  do

    if ( it_max < it ) then
      exit
    end if

    if ( abs ( ( x - w * exp ( w ) ) ) < &
      tol * abs ( ( w + 1.0D+00 ) * exp ( w ) ) ) then
      exit
    end if

    w = w - ( w * exp ( w ) - x ) &
      / ( ( w + 1.0D+00 ) * exp ( w ) &
      - ( w + 2.0D+00 ) * ( w * exp ( w ) - x ) &
      / ( 2.0D+00 * w + 2.0D+00 ) )

    it = it + 1

  end do

  lambert_w = w

  return
end
function lambert_w_crude ( x )

!*****************************************************************************80
!
!! LAMBERT_W_CRUDE is a crude estimate of the Lambert W function.
!
!  Discussion:
!
!    This crude approximation can be used as a good starting point
!    for an iterative process.
!
!    The function W(X) is defined implicitly by:
!
!      W(X) * e^W(X) = X
!
!    The function is also known as the "Omega" function.
!
!    In Mathematica, the function can be evaluated by:
!
!      W = ProductLog [ X ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Corless, Gaston Gonnet, David Hare, David Jeffrey, Donald Knuth,
!    On the Lambert W Function,
!    Advances in Computational Mathematics,
!    Volume 5, 1996, pages 329-359.
!
!    Brian Hayes,
!    "Why W?",
!    The American Scientist,
!    Volume 93, March-April 2005, pages 104-108.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) LAMBERT_W_CRUDE, a crude approximation 
!    to the Lambert W function.
!
  implicit none

  real ( kind = 8 ) lambert_w_crude
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( x <= 500.0D+00 ) then

    value = 0.04D+00 + 0.665D+00 &
      * ( 1.0D+00 + 0.0195D+00 * log ( x + 1.0D+00 ) ) * log ( x + 1.0D+00 )

  else

    value = log ( x - 4.0D+00 ) &
      - ( 1.0D+00 - 1.0D+00 / log ( x ) ) * log ( log ( x ) )

  end if

  lambert_w_crude = value

  return
end
subroutine lambert_w_values ( n_data, x, fx )

!*****************************************************************************80
!
!! LAMBERT_W_VALUES returns some values of the Lambert W function.
!
!  Discussion:
!
!    The function W(X) is defined implicitly by:
!
!      W(X) * e^W(X) = X
!
!    The function is also known as the "Omega" function.
!
!    In Mathematica, the function can be evaluated by:
!
!      W = ProductLog [ X ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Corless, Gaston Gonnet, David Hare, David Jeffrey, Donald Knuth,
!    On the Lambert W Function,
!    Advances in Computational Mathematics,
!    Volume 5, 1996, pages 329-359.
!
!    Brian Hayes,
!    "Why W?",
!    The American Scientist,
!    Volume 93, March-April 2005, pages 104-108.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.3517337112491958D+00, &
    0.5671432904097839D+00, &
    0.7258613577662263D+00, &
    0.8526055020137255D+00, &
    0.9585863567287029D+00, &
    0.1000000000000000D+01, &
    0.1049908894964040D+01, &
    0.1130289326974136D+01, &
    0.1202167873197043D+01, &
    0.1267237814307435D+01, &
    0.1326724665242200D+01, &
    0.1381545379445041D+01, &
    0.1432404775898300D+01, &
    0.1479856830173851D+01, &
    0.1524345204984144D+01, &
    0.1566230953782388D+01, &
    0.1605811996320178D+01, &
    0.1745528002740699D+01, &
    0.3385630140290050D+01, &
    0.5249602852401596D+01, &
    0.1138335808614005D+02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.2000000000000000D+01, &
    0.2500000000000000D+01, &
    0.2718281828459045D+01, &
    0.3000000000000000D+01, &
    0.3500000000000000D+01, &
    0.4000000000000000D+01, &
    0.4500000000000000D+01, &
    0.5000000000000000D+01, &
    0.5500000000000000D+01, &
    0.6000000000000000D+01, &
    0.6500000000000000D+01, &
    0.7000000000000000D+01, &
    0.7500000000000000D+01, &
    0.8000000000000000D+01, &
    0.1000000000000000D+02, &
    0.1000000000000000D+03, &
    0.1000000000000000D+04, &
    0.1000000000000000D+07 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_associated ( n, m, x, cx )

!*****************************************************************************80
!
!! LEGENDRE_ASSOCIATED evaluates the associated Legendre functions.
!
!  Differential equation:
!
!    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
!
!  First terms:
!
!    M = 0  ( = Legendre polynomials of first kind P(N,X) )
!
!    P00 =    1
!    P10 =    1 X
!    P20 = (  3 X^2 -     1)/2
!    P30 = (  5 X^3 -   3 X)/2
!    P40 = ( 35 X^4 -  30 X^2 +     3)/8
!    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
!    P60 = (231 X^6 - 315 X^4 + 105 X^2 -    5)/16
!    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
!
!    M = 1
!
!    P01 =   0
!    P11 =   1 * SQRT(1-X*X)
!    P21 =   3 * SQRT(1-X*X) * X
!    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
!    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
!
!    M = 2
!
!    P02 =   0
!    P12 =   0
!    P22 =   3 * (1-X*X)
!    P32 =  15 * (1-X*X) * X
!    P42 = 7.5 * (1-X*X) * (7*X*X-1)
!
!    M = 3
!
!    P03 =   0
!    P13 =   0
!    P23 =   0
!    P33 =  15 * (1-X*X)^1.5
!    P43 = 105 * (1-X*X)^1.5 * X
!
!    M = 4
!
!    P04 =   0
!    P14 =   0
!    P24 =   0
!    P34 =   0
!    P44 = 105 * (1-X*X)^2
!
!  Recursion:
!
!    if N < M:
!      P(N,M) = 0
!    if N = M:
!      P(N,M) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
!      all the odd integers less than or equal to N.
!    if N = M+1:
!      P(N,M) = X*(2*M+1)*P(M,M)
!    if M+1 < N:
!      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
!
!  Special values:
!
!    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
!    function of the first kind equals the Legendre polynomial of the
!    first kind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2004
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
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to be
!    evaluated.  X must satisfy -1 <= X <= 1.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 functions.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) somx2
  real ( kind = 8 ) x

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M is ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop
  end if
 
  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but M must be less than or equal to N.'
    stop
  end if
 
  if ( x < -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no less than -1.'
    stop
  end if
 
  if ( 1.0D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no more than 1.'
    stop
  end if
  
  cx(0:m-1) = 0.0D+00

  cx(m) = 1.0D+00
  somx2 = sqrt ( 1.0D+00 - x * x )
 
  fact = 1.0D+00
  do i = 1, m
    cx(m) = -cx(m) * fact * somx2
    fact = fact + 2.0D+00
  end do
 
  if ( m + 1 <= n ) then
    cx(m+1) = x * real ( 2 * m + 1, kind = 8 ) * cx(m)
  end if

  do i = m + 2, n
    cx(i) = ( real ( 2 * i     - 1, kind = 8 ) * x * cx(i-1) &
            + real (   - i - m + 1, kind = 8 ) *     cx(i-2) ) &  
            / real (     i - m,     kind = 8 )
  end do

  return
end
subroutine legendre_associated_normalized ( n, m, x, cx )

!*****************************************************************************80
!
!! LEGENDRE_ASSOCIATED_NORMALIZED: normalized associated Legendre functions.
!
!  Discussion:
!
!    The unnormalized associated Legendre functions P_N^M(X) have
!    the property that
!
!      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX 
!      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
!
!    By dividing the function by the square root of this term,
!    the normalized associated Legendre functions have norm 1.
!
!    However, we plan to use these functions to build spherical
!    harmonics, so we use a slightly different normalization factor of
!
!      sqrt ( ( ( 2 * N + 1 ) * ( N - M )! ) / ( 4 * pi * ( N + M )! ) ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 2005
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
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to be
!    evaluated.  X must satisfy -1 <= X <= 1.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 functions.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) somx2
  real ( kind = 8 ) x

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M is ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop
  end if
 
  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but M must be less than or equal to N.'
    stop
  end if
 
  if ( x < -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no less than -1.'
    stop
  end if
 
  if ( 1.0D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no more than 1.'
    stop
  end if
!
!  Entries 0 through M-1 are zero!
!
  cx(0:m-1) = 0.0D+00

  cx(m) = 1.0D+00
  somx2 = sqrt ( 1.0D+00 - x * x )
 
  factor = 1.0D+00
  do i = 1, m
    cx(m) = - cx(m) * factor * somx2
    factor = factor + 2.0D+00
  end do
 
  if ( m+1 <= n ) then
    cx(m+1) = x * real ( 2 * m + 1, kind = 8 ) * cx(m)
  end if

  do i = m+2, n
    cx(i) = ( real ( 2 * i     - 1, kind = 8 ) * x * cx(i-1) &
            + real (   - i - m + 1, kind = 8 ) *     cx(i-2) ) &  
            / real (     i - m,     kind = 8 )
  end do
!
!  Normalization.
!
  do mm = m, n
    factor = sqrt ( ( real ( 2 * mm + 1, kind = 8 ) &
      * r8_factorial ( mm - m ) ) &
      / ( 4.0D+00 * pi * r8_factorial ( mm + m ) ) )
    cx(mm) = cx(mm) * factor
  end do

  return
end
subroutine legendre_associated_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! LEGENDRE_ASSOCIATED_VALUES returns values of associated Legendre functions.
!
!  Discussion:
!
!    The function considered is the associated Legendre function P^M_N(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 June 2001
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and 
!    N_DATA is set to the index of the test data.  On each subsequent 
!    call, N_DATA is incremented and that test data is returned.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, integer M, real ( kind = 8 ) X, 
!    the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 19

  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     1,  1,  1,  1, &
     1,  2,  2,  2, &
     3,  3,  3,  3, &
     4,  5,  6,  7, &
     8,  9, 10 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( nmax ) :: m_vec = (/ &
     0,  0,  0,  0, &
     1,  0,  1,  2, &
     0,  1,  2,  3, &
     2,  2,  3,  3, &
     4,  4,  5 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( nmax ) :: fx_vec = (/ &
     0.000000D+00,  0.500000D+00,  0.707107D+00,  1.000000D+00, &
    -0.866025D+00, -0.125000D+00, -1.29904D+00,   2.25000D+00, &
    -0.437500D+00, -0.324759D+00,  5.62500D+00,  -9.74278D+00, &
     4.21875D+00,  -4.92187D+00,   12.7874D+00,   116.685D+00, &
    -1050.67D+00,  -2078.49D+00,   30086.2D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( nmax ) :: x_vec = (/ &
    0.0D+00,       0.5D+00,       0.7071067D+00, 1.0D+00, &
    0.5D+00,       0.5D+00,       0.5D+00,       0.5D+00, &
    0.5D+00,       0.5D+00,       0.5D+00,       0.5D+00, &
    0.5D+00,       0.5D+00,       0.5D+00,       0.5D+00, &
    0.5D+00,       0.5D+00,       0.5D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_associated_normalized_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! LEGENDRE_ASSOCIATED_NORMALIZED_VALUES: normalized associated Legendre.
!
!  Discussion:
!
!    The function considered is the associated Legendre polynomial P^M_N(X).
!
!    In Mathematica, the function can be evaluated by:
!
!      LegendreP [ n, m, x ]
!
!    The function is normalized by dividing by 
!
!      sqrt ( 4 * pi * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2010
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer M, real ( kind = 8 ) X, 
!    the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.2820947917738781D+00, &
     0.2443012559514600D+00, &
    -0.2992067103010745D+00, &
    -0.07884789131313000D+00, &
    -0.3345232717786446D+00, &
     0.2897056515173922D+00, &
    -0.3265292910163510D+00, &
    -0.06997056236064664D+00, &
     0.3832445536624809D+00, &
    -0.2709948227475519D+00, &
    -0.2446290772414100D+00, &
     0.2560660384200185D+00, &
     0.1881693403754876D+00, &
    -0.4064922341213279D+00, &
     0.2489246395003027D+00, &
     0.08405804426339821D+00, &
     0.3293793022891428D+00, &
    -0.1588847984307093D+00, &
    -0.2808712959945307D+00, &
     0.4127948151484925D+00, &
    -0.2260970318780046D+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 1, 0, &
    1, 2, 0, 1, &
    2, 3, 0, 1, &
    2, 3, 4, 0, &
    1, 2, 3, 4, &
    5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  1,  2, &
    2,  2,  3,  3, &
    3,  3,  4,  4, &
    4,  4,  4,  5, &
    5,  5,  5,  5, &
    5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_function_q ( n, x, cx )

!*****************************************************************************80
!
!! LEGENDRE_FUNCTION_Q evaluates the Legendre Q functions.
!
!  Differential equation:
!
!    (1-X*X) Y'' - 2 X Y' + N (N+1) = 0
!
!  First terms:
!
!    Q(0,X) = 0.5 * log((1+X)/(1-X))
!    Q(1,X) = Q(0,X)*    X - 1 
!    Q(2,X) = Q(0,X)*( 3*X^2-1)/4 - 1.5*X
!    Q(3,X) = Q(0,X)*( 5*X^3-3*X)/4 - 2.5*X^2 + 2/3
!    Q(4,X) = Q(0,X)*(35*X^4-30*X^2+3)/16 - 35/8 * X^3 + 55/24 * X
!    Q(5,X) = Q(0,X)*(63*X^5-70*X^3+15*X)/16 - 63/8*X^4 + 49/8*X^2 - 8/15
!
!  Recursion:
!
!    Q(0) = 0.5 * log ( (1+X) / (1-X) )
!    Q(1) = 0.5 * X * log ( (1+X) / (1-X) ) - 1.0
!
!    Q(N) = ( (2*N-1) * X * Q(N-1) - (N-1) * Q(N-2) ) / N
!
!  Restrictions:
!
!    -1 < X < 1
!
!  Special values:
!
!    Note that the Legendre function Q(N,X) is equal to the
!    associated Legendre function of the second kind,
!    Q(N,M,X) with M = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2003
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
!    Input, integer ( kind = 4 ) N, the highest order function to evaluate.
!
!    Input, real ( kind = 8 ) X, the point at which the functions are to be
!    evaluated.  X must satisfy -1 < X < 1.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 Legendre
!    functions at the point X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x
!
!  Check the value of X.
!
  if ( x <= -1.0D+00 .or. 1.0D+00 <= x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_FUNCTION_Q - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of X = ', x
    write ( *, '(a)' ) '  But X must be between -1 and 1.'
    stop
  end if
 
  if ( n < 0 ) then
    return
  end if

  cx(0) = 0.5D+00 * log ( ( 1.0D+00 + x ) / ( 1.0D+00 - x ) )

  if ( n == 0 ) then
    return
  end if

  cx(1) = x * cx(0) - 1.0D+00

  do i = 2, n
    cx(i) = ( real ( 2 * i - 1, kind = 8 ) * x * cx(i-1) &
            + real (   - i + 1, kind = 8 )     * cx(i-2) ) &
            / real (     i,     kind = 8 )
  end do 
 
  return
end
subroutine legendre_function_q_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! LEGENDRE_FUNCTION_Q_VALUES returns values of the Legendre Q function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( nmax ) :: fx_vec = (/ &
     0.00000000D+00, -1.00000000D+00,  0.00000000D+00, &
     0.66666667D+00, -0.40634921D+00,  0.00000000D+00, &
     0.54930614D+00, -0.72534693D+00, -0.81866327D+00, &
    -0.19865477D+00, -0.11616303D+00,  0.29165814D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  1,  2, &
     3,  9, 10, &
     0,  1,  2, &
     3,  9, 10 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( nmax ) :: x_vec = (/ &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.5D+00,  0.5D+00,  0.5D+00, &
    0.5D+00,  0.5D+00,  0.5D+00  /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_poly ( n, x, cx, cpx )

!*****************************************************************************80
!
!! LEGENDRE_POLY evaluates the Legendre polynomials P(N,X) at X.
!
!  Discussion:
!
!    P(N,1) = 1.
!    P(N,-1) = (-1)^N.
!    | P(N,X) | <= 1 in [-1,1].
!
!    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
!    function of the first kind and order N equals the Legendre polynomial
!    of the first kind and order N.
!
!    The N zeroes of P(N,X) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!    The Legendre polynomials are orthogonal under the inner product defined
!    as integration from -1 to 1:
!
!      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
!        = 0 if I =/= J
!        = 2 / ( 2*I+1 ) if I = J.
!
!    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!      C0*P(0,X) + C1*P(1,X) + ... + CN*P(N,X)
!    where
!      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,X) dx.
!
!    The formula is:
!
!      P(N,X) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
!
!  Differential equation:
!
!    (1-X*X) * P(N,X)'' - 2 * X * P(N,X)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0,X) =       1
!    P( 1,X) =       1 X
!    P( 2,X) =  (    3 X^2 -       1)/2
!    P( 3,X) =  (    5 X^3 -     3 X)/2
!    P( 4,X) =  (   35 X^4 -    30 X^2 +     3)/8
!    P( 5,X) =  (   63 X^5 -    70 X^3 +    15 X)/8
!    P( 6,X) =  (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
!    P( 7,X) =  (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
!    P( 8,X) =  ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
!    P( 9,X) =  (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
!    P(10,X) =  (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2
!                 -63 ) /256
!
!  Recursion:
!
!    P(0,X) = 1
!    P(1,X) = X
!    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!
!    P'(0,X) = 0
!    P'(1,X) = 1
!    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2004
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
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the Legendre polynomials 
!    of order 0 through N at the point X.
!
!    Output, real ( kind = 8 ) CPX(0:N), the values of the derivatives of the
!    Legendre polynomials of order 0 through N at the point X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  real ( kind = 8 ) cpx(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0D+00
  cpx(0) = 0.0D+00

  if ( n < 1 ) then
    return
  end if

  cx(1) = x
  cpx(1) = 1.0D+00
 
  do i = 2, n
 
    cx(i) = ( real ( 2 * i - 1, kind = 8 ) * x * cx(i-1)   &
            - real (     i - 1, kind = 8 ) *     cx(i-2) ) &
            / real (     i,     kind = 8 )
 
    cpx(i) = ( real ( 2 * i - 1, kind = 8 ) * ( cx(i-1) + x * cpx(i-1) ) &
             - real (     i - 1, kind = 8 ) *   cpx(i-2)               ) &
             / real (     i,     kind = 8 )
 
  end do
 
  return
end
subroutine legendre_poly_coef ( n, c )

!*****************************************************************************80
!
!! LEGENDRE_POLY_COEF evaluates the Legendre polynomial coefficients.
!
!  First terms:
!
!     1
!     0     1
!    -1/2   0      3/2
!     0    -3/2    0     5/2
!     3/8   0    -30/8   0     35/8
!     0    15/8    0   -70/8    0     63/8
!    -5/16  0    105/16  0   -315/16   0    231/16
!     0   -35/16   0   315/16   0   -693/16   0    429/16
!
!     1.00000
!     0.00000  1.00000
!    -0.50000  0.00000  1.50000
!     0.00000 -1.50000  0.00000  2.5000
!     0.37500  0.00000 -3.75000  0.00000  4.37500
!     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
!    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
!     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 February 2003
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
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the 
!    Legendre polynomials of degree 0 through N.  Each polynomial is 
!    stored as a row.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n <= 0 ) then
    return
  end if

  c(1,1) = 1.0D+00
 
  do i = 2, n
    c(i,0:i-2) =          real (   - i + 1, kind = 8 ) * c(i-2,0:i-2) &
                        / real (     i,     kind = 8 )
    c(i,1:i) = c(i,1:i) + real ( i + i - 1, kind = 8 ) * c(i-1,0:i-1) &
                        / real (     i,     kind = 8 )
  end do
 
  return
end
subroutine legendre_poly_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! LEGENDRE_POLY_VALUES returns values of the Legendre polynomials.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      LegendreP [ n, x ]
!
!  Differential equation:
!
!    (1-X*X) * P(N,X)'' - 2 * X * P(N,X)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0,X) =       1
!    P( 1,X) =       1 X
!    P( 2,X) =  (    3 X^2 -       1)/2
!    P( 3,X) =  (    5 X^3 -     3 X)/2
!    P( 4,X) =  (   35 X^4 -    30 X^2 +     3)/8
!    P( 5,X) =  (   63 X^5 -    70 X^3 +    15 X)/8
!    P( 6,X) =  (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
!    P( 7,X) =  (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
!    P( 8,X) =  ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
!    P( 9,X) =  (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
!    P(10,X) =  (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2
!                 -63 ) /256
!
!  Recursion:
!
!    P(0,X) = 1
!    P(1,X) = X
!    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!
!    P'(0,X) = 0
!    P'(1,X) = 1
!    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!
!  Formula:
!
!    P(N)(X) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
!
!  Orthogonality:
!
!    Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
!      = 0 if I =/= J
!      = 2 / ( 2*I+1 ) if I = J.
!
!  Approximation:
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!
!      C0*P(0,X) + C1*P(1,X) + ... + CN*P(N,X)
!
!    where
!
!      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,X) dx.
!
!  Special values:
!
!    P(N,1) = 1.
!    P(N,-1) = (-1)^N.
!    | P(N,X) | <= 1 in [-1,1].
!
!    Pm(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
!    function of the first kind and order N equals the Legendre polynomial
!    of the first kind and order N.
!
!    The N zeroes of P(N,X) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.2500000000000000D+00, &
    -0.4062500000000000D+00, &
    -0.3359375000000000D+00, &
     0.1577148437500000D+00, &
     0.3397216796875000D+00, &
     0.2427673339843750D-01, &
    -0.2799186706542969D+00, &
    -0.1524540185928345D+00, &
     0.1768244206905365D+00, &
     0.2212002165615559D+00, &
     0.0000000000000000D+00, &
    -0.1475000000000000D+00, &
    -0.2800000000000000D+00, &
    -0.3825000000000000D+00, &
    -0.4400000000000000D+00, &
    -0.4375000000000000D+00, &
    -0.3600000000000000D+00, &
    -0.1925000000000000D+00, &
     0.8000000000000000D-01, &
     0.4725000000000000D+00, &
     0.1000000000000000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine legendre_symbol ( q, p, l )

!*****************************************************************************80
!
!! LEGENDRE_SYMBOL evaluates the Legendre symbol (Q/P).
!
!  Discussion:
!
!    Let P be an odd prime.  Q is a QUADRATIC RESIDUE modulo P
!    if there is an integer R such that R^2 = Q ( mod P ).
!    The Legendre symbol ( Q / P ) is defined to be:
!
!      + 1 if Q ( mod P ) /= 0 and Q is a quadratic residue modulo P,
!      - 1 if Q ( mod P ) /= 0 and Q is not a quadratic residue modulo P,
!        0 if Q ( mod P ) == 0.
!
!    We can also define ( Q / P ) for P = 2 by:
!
!      + 1 if Q ( mod P ) /= 0
!        0 if Q ( mod P ) == 0
!
!    For any prime P, exactly half of the integers from 1 to P-1
!    are quadratic residues.
!
!    ( 0 / P ) = 0.
!
!    ( Q / P ) = ( mod ( Q, P ) / P ).
!
!    ( Q / P ) = ( Q1 / P ) * ( Q2 / P ) if Q = Q1 * Q2.
!
!    If Q is prime, and P is prime and greater than 2, then:
!
!      if ( Q == 1 ) then
!
!        ( Q / P ) = 1
!
!      else if ( Q == 2 ) then
!
!        ( Q / P ) = + 1 if mod ( P, 8 ) = 1 or mod ( P, 8 ) = 7,
!        ( Q / P ) = - 1 if mod ( P, 8 ) = 3 or mod ( P, 8 ) = 5.
!
!      else
!
!        ( Q / P ) = - ( P / Q ) if Q = 3 ( mod 4 ) and P = 3 ( mod 4 ),
!                  =   ( P / Q ) otherwise.
!
!  Example:
!
!    (0/7) =   0
!    (1/7) = + 1  ( 1^2 = 1 mod 7 )
!    (2/7) = + 1  ( 3^2 = 2 mod 7 )
!    (3/7) = - 1
!    (4/7) = + 1  ( 2^2 = 4 mod 7 )
!    (5/7) = - 1
!    (6/7) = - 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Pinter,
!    A Book of Abstract Algebra,
!    McGraw Hill, 1982, pages 236-237.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 86-87.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Q, an integer whose Legendre symbol with
!    respect to P is desired.
!
!    Input, integer ( kind = 4 ) P, a prime number, greater than 1, with respect
!    to which the Legendre symbol of Q is desired.
!
!    Output, integer ( kind = 4 ) L, the Legendre symbol (Q/P).
!    Ordinarily, L will be -1, 0 or 1.
!    L = -2, P is less than or equal to 1.
!    L = -3, P is not prime.
!    L = -4, the internal stack of factors overflowed.
!    L = -5, not enough factorization space.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxfactor = 20
  integer ( kind = 4 ), parameter :: maxstack = 50

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  logical i4_is_prime
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) nmore
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) p
  integer ( kind = 4 ) power(maxfactor)
  integer ( kind = 4 ) pp
  integer ( kind = 4 ) pstack(maxstack)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) qq
  integer ( kind = 4 ) qstack(maxstack)
!
!  P must be greater than 1.
!
  if ( p <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SYMBOL - Fatal error!'
    write ( *, '(a)' ) '  P must be greater than 1.'
    l = -2
    stop
  end if
!
!  P must be prime.
!
  if ( .not. i4_is_prime ( p ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SYMBOL - Fatal error!'
    write ( *, '(a)' ) '  P is not prime.'
    l = -3
    stop
  end if
!
!  ( k*P / P ) = 0.
!
  if ( mod ( q, p ) == 0 ) then
    l = 0
    return
  end if
!
!  For the special case P = 2, (Q/P) = 1 for all odd numbers.
!
  if ( p == 2 ) then
    l = 1
    return
  end if
!
!  Make a copy of Q, and force it to be nonnegative.
!
  qq = q

  do while ( qq < 0 )
    qq = qq + p
  end do

  nstack = 0
  pp = p
  l = 1

  do

    qq = mod ( qq, pp )
!
!  Decompose QQ into factors of prime powers.
!
    call i4_factor ( qq, maxfactor, nfactor, factor, power, nleft )

    if ( nleft /= 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEGENDRE_SYMBOL - Fatal error!'
      write ( *, '(a)' ) '  Not enough factorization space.'
      l = - 5
      stop
    end if
!
!  Each factor which is an odd power is added to the stack.
!
    nmore = 0

    do i = 1, nfactor

      if ( mod ( power(i), 2 ) == 1 ) then

        nmore = nmore + 1
        nstack = nstack + 1

        if ( maxstack < nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LEGENDRE_SYMBOL - Fatal error!'
          write ( *, '(a)' ) '  Stack overflow!'
          l = - 4
          stop
        end if

        pstack(nstack) = pp
        qstack(nstack) = factor(i)

      end if

    end do

    if ( nmore /= 0 ) then

      qq = qstack(nstack)
      nstack = nstack - 1
!
!  Check for a QQ of 1 or 2.
!
      if ( qq == 1 ) then

        l = + 1 * l

      else if ( qq == 2 .and. &
              ( mod ( pp, 8 ) == 1 .or. mod ( pp, 8 ) == 7 ) ) then

        l = + 1 * l

      else if ( qq == 2 .and. &
              ( mod ( pp, 8 ) == 3 .or. mod ( pp, 8 ) == 5 ) ) then

        l = - 1 * l

      else

        if ( mod ( pp, 4 ) == 3 .and. mod ( qq, 3 ) == 3 ) then
          l = - 1 * l
        end if

        call i4_swap ( pp, qq )

        cycle

      end if

    end if
!
!  If the stack is empty, we're done.
!
    if ( nstack == 0 ) then
      exit
    end if
!
!  Otherwise, get the last P and Q from the stack, and process them.
!
    pp = pstack(nstack)
    qq = qstack(nstack)
    nstack = nstack - 1

  end do

  return
end
function lerch ( z, s, a )

!*****************************************************************************80
!
!! LERCH estimates the Lerch transcendent function.
!
!  Discussion:
!
!    The Lerch transcendent function is defined as:
!
!      LERCH ( Z, S, A ) = Sum ( 0 <= K < +oo ) Z^K / ( A + K )^S
!
!    excluding any term with ( A + K ) = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      LerchPhi[z,s,a]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Thanks:
!
!    Oscar van Vlijmen
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Z, integer S, real ( kind = 8 ) A, 
!    the parameters of the function.
!
!    Output, real ( kind = 8 ) LERCH, an approximation to the Lerch
!    transcendent function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) eps
  integer ( kind = 4 ) k
  real ( kind = 8 ) lerch
  integer ( kind = 4 ) s
  real ( kind = 8 ) term
  real ( kind = 8 ) total
  real ( kind = 8 ) z
  real ( kind = 8 ) z_k

  if ( z <= 0.0D+00 ) then
    lerch = 0.0D+00
    return
  end if

  eps = 1.0D-10
  total = 0.0D+00
  k = 0
  z_k = 1.0D+00

  do

    if ( a + real ( k, kind = 8 ) /= 0.0D+00 ) then

      term = z_k / ( a + real ( k, kind = 8 ) )**s
      total = total + term

      if ( abs ( term ) <= eps * ( 1.0D+00 + abs ( total ) ) ) then
        exit
      end if

    end if

    k = k + 1
    z_k = z_k * z

  end do

  lerch = total

  return
end
subroutine lerch_values ( n_data, z, s, a, fx )

!*****************************************************************************80
!
!! LERCH_VALUES returns some values of the Lerch transcendent function.
!
!  Discussion:
!
!    The Lerch function is defined as
!
!      Phi(z,s,a) = Sum ( 0 <= k < +oo ) z^k / ( a + k )^s
!
!    omitting any terms with ( a + k ) = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      LerchPhi[z,s,a]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) Z, the parameters of the function.
!
!    Output, integer ( kind = 4 ) S, the parameters of the function.
!
!    Output, real ( kind = 8 ) A, the parameters of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    3.0D+00, &
    3.0D+00, &
    3.0D+00 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1644934066848226D+01, &
    0.1202056903159594D+01, &
    0.1000994575127818D+01, &
    0.1164481052930025D+01, &
    0.1074426387216080D+01, &
    0.1000492641212014D+01, &
    0.2959190697935714D+00, &
    0.1394507503935608D+00, &
    0.9823175058446061D-03, &
    0.1177910993911311D+00, &
    0.3868447922298962D-01, &
    0.1703149614186634D-04 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) s
  integer ( kind = 4 ), save, dimension ( n_max ) :: s_vec = (/ &
     2, 3, 10, &
     2, 3, 10, &
     2, 3, 10, &
     2, 3, 10 /)
  real ( kind = 8 ) z
  real ( kind = 8 ), save, dimension ( n_max ) :: z_vec = (/ &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.3333333333333333D+00, &
    0.3333333333333333D+00, &
    0.3333333333333333D+00, &
    0.1000000000000000D+00, &
    0.1000000000000000D+00, &
    0.1000000000000000D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    z = 0.0D+00
    s = 0
    a = 0.0D+00
    fx = 0.0D+00
  else
    z = z_vec(n_data)
    s = s_vec(n_data)
    a = a_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine lock ( n, a )

!*****************************************************************************80
!
!! LOCK returns the number of codes for a lock with N buttons.
!
!  Discussion:
!
!    A button lock has N numbered buttons.  To open the lock, groups
!    of buttons must be pressed in the correct order.  Each button
!    may be pushed no more than once.  Thus, a code for the lock is 
!    an ordered list of the groups of buttons to be pushed.
!
!    For this discussion, we will assume that EVERY button is pushed
!    at some time, as part of the code.  To count the total number
!    of codes, including those which don't use all the buttons, then
!    the number is 2 * A(N), or 2 * A(N) - 1 if we don't consider the
!    empty code to be valid.
!
!    If there are 3 buttons, then there are 13 possible "full button" codes:
!
!      (123)
!      (12) (3)
!      (13) (2)
!      (23) (1)
!      (1) (23)
!      (2) (13)
!      (3) (12)
!      (1) (2) (3)
!      (1) (3) (2)
!      (2) (1) (3)
!      (2) (3) (1)
!      (3) (1) (2)
!      (3) (2) (1)
!
!    and, if we don't need to push all the buttons, every "full button" code 
!    above yields a distinct "partial button" code by dropping the last set 
!    of buttons:
!
!      ()
!      (12)
!      (13)
!      (23)
!      (1)
!      (2)
!      (3)
!      (1) (2)
!      (1) (3)
!      (2) (1)
!      (2) (3)
!      (3) (1)
!      (3) (2)
!
!  Example:
!
!     N         A(N)
!     0           1
!     1           1
!     2           3
!     3          13
!     4          75
!     5         541
!     6        4683
!     7       47293
!     8      545835
!     9     7087261
!    10   102247563
!
!  Recursion:
!
!    A(I) = sum ( 0 <= J < I ) Binomial ( I, N-J ) * A(J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Velleman, Gregory Call,
!    Permutations and Combination Locks,
!    Mathematics Magazine,
!    Volume 68, Number 4, October 1995, pages 243-253.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum number of lock buttons.
!
!    Output, integer ( kind = 4 ) A(0:N), the number of lock codes.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    return
  end if

  a(0) = 1

  do i = 1, n
    a(i) = 0
    do j = 0, i - 1
      a(i) = a(i) + i4_choose ( i, i - j ) * a(j)
    end do
  end do

  return
end
subroutine meixner ( n, beta, c, x, v )

!*****************************************************************************80
!
!! MEIXNER evaluates Meixner polynomials at a point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Walter Gautschi,
!    Orthogonal Polynomials: Computation and Approximation,
!    Oxford, 2004,
!    ISBN: 0-19-850672-4,
!    LC: QA404.5 G3555.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum order of the polynomial.  
!    N must be at least 0.
!
!    Input, real ( kind = 8 ) BETA, the Beta parameter.  0 < BETA.
!
!    Input, real ( kind = 8 ) C, the C parameter.  0 < C < 1.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) V(0:N), the value of the polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  real ( kind = 8 ) v(0:n)
  real ( kind = 8 ) x

  if ( beta <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MEIXNER - Fatal error!'
    write ( *, '(a)' ) '  Parameter BETA must be positive.'
    stop
  end if

  if ( c <= 0.0D+00 .or. 1.0D+00 <= c ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MEIXNER - Fatal error!'
    write ( *, '(a)' ) '  Parameter C must be strictly between 0 and 1.'
    stop
  end if

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MEIXNER - Fatal error!'
    write ( *, '(a)' ) '  Parameter N must be nonnegative.'
    stop
  end if

  v(0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  v(1) = ( c - 1.0D+00 ) * x / beta / c + 1.0D+00

  if ( n == 1 ) then
    return
  end if

  do i = 1, n - 1
    v(i+1) = ( &
      ( ( c - 1.0D+00 ) * x + ( 1.0D+00 + c ) &
      * real ( i, kind = 8 ) + beta * c ) * v(i) &
      - real ( i, kind = 8 ) * v(i-1) &
      ) / ( real ( i, kind = 8 ) + beta )
  end do

  return
end
function mertens ( n )

!*****************************************************************************80
!
!! MERTENS evaluates the Mertens function.
!
!  Discussion:
!
!    The Mertens function M(N) is the sum from 1 to N of the Moebius
!    function MU.  That is,
!
!    M(N) = sum ( 1 <= I <= N ) MU(I)
!
!        N   M(N)
!        --  ----
!         1     1
!         2     0
!         3    -1
!         4    -1
!         5    -2
!         6    -1
!         7    -2
!         8    -2
!         9    -2
!        10    -1
!        11    -2
!        12    -2
!       100     1
!      1000     2
!     10000   -23
!    100000   -48
!
!    The determinant of the Redheffer matrix of order N is equal
!    to the Mertens function M(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    M Deleglise, J Rivat,
!    Computing the Summation of the Moebius Function,
!    Experimental Mathematics,
!    Volume 5, 1996, pages 291-295.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument.
!
!    Output, integer ( kind = 4 ) MERTENS, the value.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) mertens
  integer ( kind = 4 ) mu_i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  value = 0

  do i = 1, n
    call moebius ( i, mu_i )
    value = value + mu_i
  end do

  mertens = value

  return
end
subroutine mertens_values ( n_data, n, c )

!*****************************************************************************80
!
!! MERTENS_VALUES returns some values of the Mertens function.
!
!  Discussion:
!
!    The Mertens function M(N) is the sum from 1 to N of the Moebius
!    function MU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    M Deleglise, J Rivat,
!    Computing the Summation of the Moebius Function,
!    Experimental Mathematics,
!    Volume 5, 1996, pages 291-295.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the Mertens function.
!
!    Output, integer ( kind = 4 ) C, the value of the Mertens function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 15

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
      1,   0,  -1,   -1,  -2,  -1,  -2,  -2,   -2,  -1, &
     -2,  -2,   1,    2, -23 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     11,  12,  100, 1000, 10000 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine moebius ( n, mu )

!*****************************************************************************80
!
!! MOEBIUS returns the value of MU(N), the Moebius function of N.
!
!  Discussion:
!
!    MU(N) is defined as follows:
!
!      MU(N) = 1 if N = 1;
!              0 if N is divisible by the square of a prime;
!              (-1)^K, if N is the product of K distinct primes.
!
!    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
!    if N is a square, cube, etc.
!
!    The Moebius function is related to Euler's totient function:
!
!      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
!
!  Example:
!
!     N  MU(N)
!
!     1    1
!     2   -1
!     3   -1
!     4    0
!     5   -1
!     6    1
!     7   -1
!     8    0
!     9    0
!    10    1
!    11   -1
!    12    0
!    13   -1
!    14    1
!    15    1
!    16    0
!    17   -1
!    18    0
!    19   -1
!    20    0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value to be analyzed.
!
!    Output, integer ( kind = 4 ) MU, the value of MU(N).
!    If N is less than or equal to 0, MU will be returned as -2.
!    If there was not enough internal space for factoring, MU
!    is returned as -3.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxfactor = 20

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) power(maxfactor)

  if ( n <= 0 ) then
    mu = - 2
    return
  end if

  if ( n == 1 ) then
    mu = 1
    return
  end if
!
!  Factor N.
!
  call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOEBIUS - Fatal error!'
    write ( *, '(a)' ) '  Not enough factorization space.'
    mu = - 3
    stop
  end if

  mu = 1

  do i = 1, nfactor

    mu = - mu

    if ( 1 < power(i) ) then
      mu = 0
      return
    end if

  end do

  return
end
subroutine moebius_values ( n_data, n, c )

!*****************************************************************************80
!
!! MOEBIUS_VALUES returns some values of the Moebius function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the Moebius function.
!
!    Output, integer ( kind = 4 ) C, the value of the Moebius function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 20

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
      1,  -1,  -1,   0,  -1,   1,  -1,   0,   0,   1, &
     -1,   0,  -1,   1,   1,   0,  -1,   0,  -1,   0 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     11,  12,  13,  14,  15,  16,  17,  18,  19,  20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine motzkin ( n, a )

!*****************************************************************************80
!
!! MOTZKIN returns the Motzkin numbers up to order N.
!
!  Discussion:
!
!    The Motzkin number A(N) counts the number of distinct paths
!    from (0,0) to (0,N) in which the only steps used are
!    (1,1), (1,-1) and (1,0), and the path is never allowed to
!    go below the X axis.
!
!  Example:
!
!     N  A(N)
!
!     0    1
!     1    1
!     2    2
!     3    4
!     4    9
!     5   21
!     6   51
!     7  127
!     8  323
!     9  835
!    10 2188
!
!  Recursion:
!
!    A(N) = A(N-1) + sum ( 0 <= K <= N-2 ) A(K) * A(N-2-K)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order Motzkin number to compute.
!
!    Output, integer ( kind = 4 ) A(0:N), the Motzkin numbers.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    return
  end if

  a(0) = 1

  do i = 1, n
    a(i) = a(i-1)
    do j = 0, i - 2
      a(i) = a(i) + a(j) * a(i-2-j)
    end do
  end do

  return
end
subroutine normal_01_cdf_inv ( p, x )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10**16.
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
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an 
!    "infinite" value will be returned. 
!
!    Output, real ( kind = 8 ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = 8 ), parameter :: const1 = 0.180625D+00
  real ( kind = 8 ), parameter :: const2 = 1.6D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) r8poly_value
  real ( kind = 8 ), parameter :: split1 = 0.425D+00
  real ( kind = 8 ), parameter :: split2 = 5.0D+00
  real ( kind = 8 ) x

  if ( p <= 0.0D+00 ) then
    x = - huge ( x )
    return
  end if

  if ( 1.0D+00 <= p ) then
    x = huge ( x )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    x = q * r8poly_value ( 8, a, r ) / r8poly_value ( 8, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then

      x = huge ( x )

    else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value ( 8, c, r ) / r8poly_value ( 8, d, r )

      else

        r = r - split2
        x = r8poly_value ( 8, e, r ) / r8poly_value ( 8, f, r )

      end if

    end if

    if ( q < 0.0D+00 ) then
      x = -x
    end if

  end if

  return
end
subroutine omega ( n, ndiv )

!*****************************************************************************80
!
!! OMEGA returns OMEGA(N), the number of distinct prime divisors of N.
!
!  Discussion:
!
!    The formula is:
!
!      If N = 1, then
!
!        OMEGA(N) = 1
!
!      else if the prime factorization of N is
!
!        N = P1^E1 * P2^E2 * ... * PM^EM,
!
!      then
!
!        OMEGA(N) = M
!
!  Example:
!
!     N   OMEGA(N)
!
!     1    1
!     2    1
!     3    1
!     4    1
!     5    1
!     6    2
!     7    1
!     8    1
!     9    1
!    10    2
!    11    1
!    12    2
!    13    1
!    14    2
!    15    2
!    16    1
!    17    1
!    18    2
!    19    1
!    20    2
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value to be analyzed.  N must be 1 or
!    greater.
!
!    Output, integer ( kind = 4 ) NDIV, the value of OMEGA(N).  But if N is 0 or
!    less, NDIV is returned as 0, a nonsense value.  If there is
!    not enough room for factoring, NDIV is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxfactor = 20

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiv
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) power(maxfactor)

  if ( n <= 0 ) then
    ndiv = 0
    return
  end if

  if ( n == 1 ) then
    ndiv = 1
    return
  end if
!
!  Factor N.
!
  call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OMEGA - Fatal error!'
    write ( *, '(a)' ) '  Not enough factorization space.'
    ndiv = -1
    stop
  end if

  ndiv = nfactor

  return
end
subroutine omega_values ( n_data, n, c )

!*****************************************************************************80
!
!! OMEGA_VALUES returns some values of the OMEGA function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,s,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the OMEGA function.
!
!    Output, integer ( kind = 4 ) C, the value of the OMEGA function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 23

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
      1,   1,   1,   1,   1, &
      2,   1,   1,   1,   2, &
      3,   1,   4,   4,   3, &
      1,   5,   2,   2,   1, &
      6,   7,   8 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
           1,        2,        3,        4,        5, &
           6,        7,        8,        9,       10, &
          30,      101,      210,     1320,     1764, &
        2003,     2310,     2827,     8717,    12553, &
       30030,   510510,  9699690 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine partition_count_values ( n_data, n, c )

!*****************************************************************************80
!
!! PARTITION_COUNT_VALUES returns some values of the integer partition count.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  The number of partitions of N is symbolized
!    by P(N).  Thus, the number 5 has P(N) = 7, because it has the 
!    following partitions:
!
!    5 = 5
!      = 4 + 1 
!      = 3 + 2 
!      = 3 + 1 + 1 
!      = 2 + 2 + 1 
!      = 2 + 1 + 1 + 1 
!      = 1 + 1 + 1 + 1 + 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the integer.
!
!    Output, integer ( kind = 4 ) C, the number of partitions of the integer.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 21

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
      1, &
      1,   2,   3,   5,   7,  11,  15,  22,  30,  42, & 
     56,  77, 101, 135, 176, 231, 297, 385, 490, 627 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  &
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine partition_distinct_count_values ( n_data, n, c )

!*****************************************************************************80
!
!! PARTITION_DISTINCT_COUNT_VALUES returns some values of Q(N).
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  The number of partitions of N is symbolized
!    by P(N).  Thus, the number 5 has P(N) = 7, because it has the 
!    following partitions:
!
!    5 = 5
!      = 4 + 1 
!      = 3 + 2 
!      = 3 + 1 + 1 
!      = 2 + 2 + 1 
!      = 2 + 1 + 1 + 1 
!      = 1 + 1 + 1 + 1 + 1
!
!    However, if we require that each member of the partition
!    be distinct, so that no nonzero summand occurs more than once,
!    we are computing something symbolized by Q(N).
!    The number 5 has Q(N) = 3, because it has the following partitions 
!    into distinct parts:
!
!    5 = 5
!      = 4 + 1 
!      = 3 + 2 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the integer.
!
!    Output, integer ( kind = 4 ) C, the number of partitions of the integer
!    into distinct parts.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 21

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
      1, &
      1,   1,   2,   2,   3,   4,   5,   6,   8,  10, & 
     12,  15,  18,  22,  27,  32,  38,  46,  54,  64 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  &
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine pentagon_num ( n, p )

!*****************************************************************************80
!
!! PENTAGON_NUM computes the N-th pentagonal number.
!
!  Discussion:
!
!    The pentagonal number P(N) counts the number of dots in a figure of
!    N nested pentagons.  The pentagonal numbers are defined for both
!    positive and negative N.
!
!    The formula is:
!
!      P(N) = ( N * ( 3 * N - 1 ) ) / 2
!
!  First values:
!
!    N   P
!
!   -5  40
!   -4  26
!   -3  15
!   -2   7
!   -1   2
!    0   0
!    1   1
!    2   5
!    3  12
!    4  22
!    5  35
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the pentagonal number desired.
!
!    Output, integer ( kind = 4 ) P, the value of the N-th pentagonal number.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  p = ( n * ( 3 * n - 1 ) ) / 2

  return
end
subroutine phi ( n, phin )

!*****************************************************************************80
!
!! PHI computes the number of relatively prime predecessors of an integer.
!
!  Discussion:
!
!    PHI(N) is the number of integers between 1 and N which are
!    relatively prime to N.  I and J are relatively prime if they
!    have no common factors.  The function PHI(N) is known as
!    "Euler's totient function".
!
!    By convention, 1 and N are relatively prime.
!
!    The formula is:
!
!      PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
!
!      PHI(P^K) = P^(K-1) * ( P - 1 ) if P is prime.
!
!      PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
!
!      N = Sum ( D divides N ) PHI(D).
!
!  Example:
!
!     N  PHI(N)
!
!     1    1
!     2    1
!     3    2
!     4    2
!     5    4
!     6    2
!     7    6
!     8    4
!     9    6
!    10    4
!    11   10
!    12    4
!    13   12
!    14    6
!    15    8
!    16    8
!    17   16
!    18    6
!    19   18
!    20    8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value to be analyzed.
!
!    Output, integer ( kind = 4 ) PHIN, the value of PHI(N).  If N is less than
!    or equal to 0, PHI will be returned as 0.  If there is not enough
!    room for full factoring of N, PHI will be returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxfactor = 20

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) phin
  integer ( kind = 4 ) power(maxfactor)

  if ( n <= 0 ) then
    phin = 0
    return
  end if

  if ( n == 1 ) then
    phin = 1
    return
  end if
!
!  Factor N.
!
  call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PHI - Fatal error!'
    write ( *, '(a)' ) '  Not enough factorization space!'
    phin = -1
    stop
  end if

  phin = 1
  do i = 1, nfactor
    phin = phin * factor(i)**( power(i) - 1 ) * ( factor(i) - 1 )
  end do

  return
end
subroutine phi_values ( n_data, n, c )

!*****************************************************************************80
!
!! PHI_VALUES returns some values of the PHI function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the PHI function.
!
!    Output, integer ( kind = 4 ) C, the value of the PHI function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 20

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
      1,   1,   2,   2,   4,   2,   6,   4,   6,   4, &
      8,   8,  16,  20,  16,  40, 148, 200, 200, 648 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     20,  30,  40,  50,  60, 100, 149, 500, 750, 999 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine poly_bernoulli ( n, k, b )

!*****************************************************************************80
!
!! POLY_BERNOULLI evaluates the poly-Bernolli numbers with negative index.
!
!  Discussion:
!
!    The poly-Bernoulli numbers B_n^k were defined by M Kaneko
!    formally as the coefficients of X^n/n! in a particular power 
!    series.  He also showed that, when the super-index is negative,
!    we have
!
!      B_n^(-k) = Sum ( 0 <= j <= min ( n, k ) ) 
!        (j!)^2 * S(n+1,j+1) * S(k+1,j+1)
!
!    where S(n,k) is the Stirling number of the second kind, the number of
!    ways to partition a set of size n into k nonempty subset.
!
!    B_n^(-k) is also the number of "lonesum matrices", that is, 0-1
!    matrices of n rows and k columns which are uniquely reconstructable
!    from their row and column sums.
!
!    The poly-Bernoulli numbers get large very quickly.
!
!  Table:
!
!    \ K 0  1    2     3      4       5        6
!    N
!    0   1  1    1     1      1       1        1
!    1   1  2    4     8     16      32       64
!    2   1  4   14    46    146     454     1394
!    3   1  8   46   230   1066    4718    20266
!    4   1 16  146  1066   6902   41506   237686
!    5   1 32  454  4718  41506  329462  2441314
!    6   1 64 1394 20266 237686 2441314 22934774
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Chad Brewbaker,
!    Lonesum (0,1) Matrices and Poly-Bernoulli Numbers of Negative Index,
!    MS Thesis,
!    Iowa State University, 2005.
!
!    M Kaneko,
!    Poly-Bernoulli Numbers,
!    Journal Theorie des Nombres Bordeaux,
!    Volume 9, 1997, pages 221-228.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, the indices.  N and K should be 
!    nonnegative.
!
!    Output, integer ( kind = 4 ) B, the value of B_N^(-K).
!
  implicit none

  integer ( kind = 4 ) b
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jfact
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: s

  if ( n < 0 ) then
    b = 0
    return
  else if ( n == 0 ) then
    b = 1
    return
  end if

  if ( k < 0 ) then
    b = 0
    return
  else if ( k == 0 ) then
    b = 1
    return
  end if

  jhi = min ( n, k )
  m = max ( n, k ) + 1

  allocate ( s(1:m,1:m) )
  call stirling2 ( m, m, s )

  jfact = 1
  b = 0

  do j = 0, jhi

    b = b + jfact * jfact * s(n+1,j+1) * s(k+1,j+1)

    jfact = jfact * ( j + 1 )

  end do

  deallocate ( s )

  return
end
function poly_coef_count ( dim, degree )

!*****************************************************************************80
!
!! POLY_COEF_COUNT: polynomial coefficient count given dimension and degree.
!
!  Discussion:
!
!    To count all monomials of degree 5 or less in dimension 3,
!    we can count all monomials of degree 5 in dimension 4.
!
!    To count all monomials of degree 5 in dimension 4, we imagine
!    that each of the variables X, Y, Z and W is a "box" and that
!    we need to drop 5 pebbles into these boxes.  Every distinct
!    way of doing this represents a degree 5 monomial in dimension 4.
!    Ignoring W gives us monomials up to degree five in dimension 3.
!
!    To count them, we draw 3 lines as separators to indicate the
!    4 boxes, and then imagine all distinct sequences involving
!    the three lines and the 5 pebbles.  Indicate the lines by 1's
!    and the pebbles by 0's and we're asking for the number of
!    permutations of 3 1's and 5 0's, which is 8! / (3! 5!)
!
!    In other words, 56 = 8! / (3! 5!) is:
!    * the number of monomials of degree exactly 5 in dimension 4, 
!    * the number of monomials of degree 5 or less in dimension 3, 
!    * the number of polynomial coefficients of a polynomial of 
!      degree 5 in (X,Y,Z).
!
!    In general, the formula for the number of monomials of degree DEG
!    or less in dimension DIM is
!
!      (DEG+DIM)! / (DEG! * DIM!)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM, the dimension of the polynomial.
!    0 <= DIM.
!
!    Input, integer ( kind = 4 ) DEGREE, the degree of the polynomnial
!    0 <= DEGREE
!
!    Output, integer ( kind = 4 ) POLY_COEF_COUNT, the number of coefficients 
!    in the general polynomial of dimension DIM and degree DEGREE.
!
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) poly_coef_count

  if ( dim < 0 ) then
    poly_coef_count = -1
  else if ( degree < 0 ) then
    poly_coef_count = -1
  else
    poly_coef_count = i4_choose ( degree + dim, degree )
  end if

  return
end 
function prime ( n )

!*****************************************************************************80
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Discussion:
!
!    PRIME_MAX is 1600, and the largest prime stored is 13499.
!
!    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2005
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
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired prime number.
!    In general, is should be true that 0 <= N <= PRIME_MAX.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!
!    Output, integer ( kind = 4 ) PRIME, the N-th prime.  If N is out of range,
!    PRIME is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: prime_max = 1600

  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( prime_max ) :: npvec
  integer ( kind = 4 ) prime

  if ( icall == 0 ) then

    icall = 1

    npvec(1:100) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

    npvec(101:200) = (/ &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

    npvec(201:300) = (/ &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

    npvec(301:400) = (/ &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

    npvec(401:500) = (/ &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

    npvec(501:600) = (/ &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

    npvec(601:700) = (/ &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

    npvec(701:800) = (/ &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

    npvec(801:900) = (/ &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

    npvec(901:1000) = (/ &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

    npvec(1001:1100) = (/ &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

    npvec(1101:1200) = (/ &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

    npvec(1201:1300) = (/ &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

    npvec(1301:1400) = (/ &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

    npvec(1401:1500) = (/ &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

   npvec(1501:1600) = (/ &
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, &
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, &
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, &
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, &
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, &
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, &
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, &
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, &
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, &
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /)

  end if

  if ( n == -1 ) then
    prime = prime_max
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= prime_max ) then
    prime = npvec(n)
  else
    prime = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
  end if

  return
end
subroutine psi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! PSI_VALUES returns some values of the Psi or Digamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyGamma[x]
!
!    or
!
!      PolyGamma[0,x]
!
!    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
!
!    PSI(1) = -Euler's constant.
!
!    PSI(X+1) = PSI(X) + 1 / X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2004
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.5772156649015329D+00, &
    -0.4237549404110768D+00, &
    -0.2890398965921883D+00, &
    -0.1691908888667997D+00, &
    -0.6138454458511615D-01, &
     0.3648997397857652D-01, &
     0.1260474527734763D+00, &
     0.2085478748734940D+00, &
     0.2849914332938615D+00, &
     0.3561841611640597D+00, &
     0.4227843350984671D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function pyramid_num ( n )

!*****************************************************************************80
!
!! PYRAMID_NUM returns the N-th pyramidal number.
!
!  Discussion:
!
!    The N-th pyramidal number P(N) is formed by the sum of the first
!    N triangular numbers T(J):
!
!      T(J) = sum ( 1 <= J <= N ) J
!
!      P(N) = sum ( 1 <= I <= N ) T(I)
!
!    By convention, T(0) = 0.
!
!    The formula is:
!
!      P(N) = ( (N+1)^3 - (N+1) ) / 6
!
!  Example:
!
!    0    0
!    1    1
!    2    4
!    3   10
!    4   20
!    5   35
!    6   56
!    7   84
!    8  120
!    9  165
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired number, which 
!    must be at least 0.
!
!    Output, integer ( kind = 4 ) PYRAMID_NUM, the N-th pyramidal number.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) pyramid_num

  pyramid_num = ( ( n + 1 )**3 - ( n + 1 ) ) / 6

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r4_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
function r8_acosh ( x )

!*****************************************************************************80
!
!! R8_ACOSH returns the inverse hyperbolic cosine of a number.
!
!  Discussion:
!
!    One formula is:
!
!      R8_ACOSH = log ( X + SQRT ( X^2 - 1.0 ) )
!
!    but this formula suffers from roundoff and overflow problems.
!    The formula used here was recommended by W Kahan, as discussed
!    by Moler.
!
!    Applying the inverse function
!
!      Y = R8_ACOSH ( X )
!
!    implies that
!
!      X = COSH(Y) = 0.5 * ( EXP(Y) + EXP(-Y) ).
!
!    For every X greater than or equal to 1, there are two possible
!    choices Y such that X = COSH(Y), differing only in sign.  It
!    is usual to resolve this choice by taking the value of 
!    R8_ACOSH ( X ) to be nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cleve Moler,
!    Trigonometry is a Complex Subject,
!    MATLAB News and Notes,
!    Summer 1998.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose inverse hyperbolic 
!    cosine is desired.  X should be greater than or equal to 1.
!
!    Output, real ( kind = 8 ) R8_ACOSH, the inverse hyperbolic cosine of 
!    X.  The principal value (that is, the positive value of the two ) 
!    is returned.
!
  implicit none

  real ( kind = 8 ) r8_acosh
  real ( kind = 8 ) x

  if ( x < 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ACOSH - Fatal error!'
    write ( *, '(a)' ) '  Argument X must satisfy 1 <= X.'
    write ( *, '(a,g14.6)' ) '  The input X = ', x
    stop
  end if

  r8_acosh = 2.0D+00 * log ( &
      sqrt ( 0.5D+00 * ( x + 1.0D+00 ) ) &
    + sqrt ( 0.5D+00 * ( x - 1.0D+00 ) ) )

  return
end
function r8_asinh ( x )

!*****************************************************************************80
!
!! R8_ASINH returns the inverse hyperbolic sine of a number.
!
!  Discussion:
!
!    The assertion that:
!
!      Y = R8_ASINH ( X )
!
!    implies that
!
!      X = SINH(Y) = 0.5 * ( EXP(Y) - EXP(-Y) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose inverse hyperbolic 
!    sine is desired.
!
!    Output, real ( kind = 8 ) R8_ASINH, the inverse hyperbolic sine of X.
!
  implicit none

  real ( kind = 8 ) r8_asinh
  real ( kind = 8 ) x

  r8_asinh = log ( x + sqrt ( x * x + 1.0D+00 ) )

  return
end
function r8_atanh ( x )

!*****************************************************************************80
!
!! R8_ATANH returns the inverse hyperbolic tangent of a number.
!
!  Discussion:
!
!    Y = R8_ATANH ( X )
!
!    implies that
!
!    X = TANH(Y) = ( EXP(Y) - EXP(-Y) ) / ( EXP(Y) + EXP(-Y) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose inverse hyperbolic 
!    tangent is desired.  The absolute value of X should be less than 
!    or equal to 1.
!
!    Output, real ( kind = 8 ) R8_ATANH, the inverse hyperbolic tangent of X.
!
  implicit none

  real ( kind = 8 ) r8_atanh
  real ( kind = 8 ) x

  if ( 1.0D+00 <= abs ( x ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ATANH - Fatal error!'
    write ( *, '(a)' ) '  ABS(X) must be < 1.'
    write ( *, '(a,g14.6)' ) '  Your input is X = ', x
    stop
  end if

  r8_atanh = 0.5D+00 * log ( ( 1.0D+00 + x ) / ( 1.0D+00 - x ) )

  return
end
function r8_cas ( x )

!*****************************************************************************80
!
!! R8_CAS returns the "casine" of an R8.
!
!  Discussion:
!
!    The "casine", used in the discrete Hartley transform, is abbreviated
!    CAS(X), and defined by:
!
!      CAS(X) = cos ( X ) + sin( X )
!             = sqrt ( 2 ) * sin ( X + pi/4 )
!             = sqrt ( 2 ) * cos ( X - pi/4 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ralph Hartley,
!    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
!    Proceedings of the Institute of Radio Engineers,
!    Volume 30, pages 144-150, 1942.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose casine is desired.
!
!    Output, real ( kind = 8 ) R8_CAS, the casine of X, which will be between
!    plus or minus the square root of 2.
!
  implicit none

  real ( kind = 8 ) r8_cas
  real ( kind = 8 ) x

  r8_cas = cos ( x ) + sin ( x )

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the combinatorial coefficient C(N,K).
!
!  Discussion:
!
!    Real arithmetic is used, and C(N,K) is computed directly, via
!    Gamma functions, rather than recursively.
!
!    C(N,K) is the number of distinct combinations of K objects
!    chosen from a set of N distinct objects.  A combination is
!    like a set, in that order does not matter.
!
!    The formula is:
!
!      C(N,K) = N! / ( (N-K)! * K! )
!
!  Example:
!
!    The number of combinations of 2 things chosen from 5 is 10.
!
!    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
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
!    02 June 2007
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
 
  else if ( 1 < k .and. k < n-1 ) then
 
    arg = real ( n + 1, kind = 8 )
    facn = r8_gamma_log ( arg )
 
    arg = real ( k + 1, kind = 8 )
    fack = r8_gamma_log ( arg )
 
    arg = real ( n - k + 1, kind = 8 )
    facnmk = r8_gamma_log ( arg )
 
    value = anint ( exp ( facn - fack - facnmk ) )
 
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
function r8_cot ( angle )

!*****************************************************************************80
!
!! R8_COT returns the cotangent of an angle.
!
!  Discussion:
!
!    R8_COT ( THETA ) = COS ( THETA ) / SIN ( THETA )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in radians.
!
!    Output, real ( kind = 8 ) R8_COT, the cotangent of the angle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) r8_cot

  r8_cot = cos ( angle ) / sin ( angle )

  return
end
function r8_cot_deg ( angle )

!*****************************************************************************80
!
!! R8_COT_DEG returns the cotangent of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees.
!
!    Output, real ( kind = 8 ) R8_COT_DEG, the cotangent of the angle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) r8_cot_deg
  real ( kind = 8 ), parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00

  r8_cot_deg = cos ( degrees_to_radians * angle ) &
          / sin ( degrees_to_radians * angle )

  return
end
function r8_csc ( theta )

!*****************************************************************************80
!
!! R8_CSC returns the cosecant of X.
!
!  Discussion:
!
!    R8_CSC ( THETA ) = 1.0 / SIN ( THETA )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) THETA, the angle, in radians, whose 
!    cosecant is desired.  It must be the case that SIN ( THETA ) is not zero.
!
!    Output, real ( kind = 8 ) R8_CSC, the cosecant of THETA.
!
  implicit none

  real ( kind = 8 ) r8_csc
  real ( kind = 8 ) theta

  r8_csc = sin ( theta )

  if ( r8_csc == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_CSC - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Cosecant undefined for THETA = ', theta
    stop
  end if

  r8_csc = 1.0D+00 / r8_csc

  return
end
function r8_csc_deg ( angle )

!*****************************************************************************80
!
!! R8_CSC_DEG returns the cosecant of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees.
!
!    Output, real ( kind = 8 ) R8_CSC_DEG, the cosecant of the angle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) r8_csc_deg
  real ( kind = 8 ), parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00

  r8_csc_deg = 1.0D+00 / cos ( degrees_to_radians * angle )

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

  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) value

  value = 1.0D+00

  do while ( 1.0D+00 < 1.0D+00 + value )
    value = value / 2.0D+00
  end do

  r8_epsilon = 2.0D+00 * value

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    The formula is:
!
!      factorial ( N ) = product ( 1 <= I <= N ) I
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
function r8_factorial_log ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL_LOG computes the natural logarithm of the factorial of N.
!
!  Discussion:
!
!    The formula is:
!
!      log ( FACTORIAL ( N ) ) 
!        = log ( product ( 1 <= I <= N ) I )
!        = sum ( ( 1 <= I <= N ) log ( I ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the value is returned as 0.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL_LOG, the logarithm of
!    the factorial of N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial_log

  r8_factorial_log = 0.0D+00

  do i = 1, n
    r8_factorial_log = r8_factorial_log + log ( real ( i, kind = 8 ) )
  end do

  return
end
subroutine r8_factorial_log_values ( n_data, n, fn )

!*****************************************************************************80
!
!! R8_FACTORIAL_LOG_VALUES returns values of log(factorial(N)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2003
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
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA 
!    is set to the index of the test data.  On each subsequent call, N_DATA is
!    incremented and that test data is returned.  When there is no more
!    test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 27

  real ( kind = 8 ), save, dimension ( nmax ) :: fnvec = (/ &
      0.0D+00,         0.0D+00,       0.6931472D+00,  1.791757D+00, &
      3.178051D+00,    4.787489D+00,  6.579246D+00,   8.525160D+00, &
     10.60460D+00,    12.80182D+00,  15.10441D+00,   17.50232D+00, &
     19.98722D+00,    22.55216D+00,  25.19123D+00,   27.89927D+00, &
     30.67186D+00,    33.50508D+00,  36.39544D+00,   39.33987D+00, &
     42.33561D+00,    58.00362D+00, 148.4778D+00,   363.7394D+00, &
    605.0201D+00,   2611.331D+00,   5912.128D+00 /)
  real ( kind = 8 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: nvec = (/ &
     0,   1,    2,   3, &
     4,   5,    6,   7, &
     8,   9,   10,  11, &
    12,  13,   14,  15, &
    16,  17,   18,  19, &
    20,  25,   50, 100, &
   150, 500, 1000 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    fn = 0.0D+00
  else
    n = nvec(n_data)
    fn = fnvec(n_data)
  end if

  return
end
subroutine r8_factorial_values ( n_data, n, fn )

!*****************************************************************************80
!
!! R8_FACTORIAL_VALUES returns values of the real factorial function.
!
!  Discussion:
!
!    Although the factorial is an integer valued function, it quickly
!    becomes too large for an integer to hold.  This routine still accepts
!    an integer as the input argument, but returns the function value
!    as a real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 May 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA 
!    is set to the index of the test data.  On each subsequent call, N_DATA is
!    incremented and that test data is returned.  When there is no more
!    test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, real ( kind = 8 ) FN, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 23

  real ( kind = 8 ), save, dimension ( nmax ) :: fnvec = (/ &
    1.0D+00, 1.0D+00, 2.0D+00, 6.0D+00, &
    24.0D+00, 120.0D+00, 720.0D+00, 5040.0D+00, &
    40320.0D+00, 362880.0D+00, 3628800.0D+00, 39916800.0D+00, &
    479001600.0D+00, 6227020800.0D+00, 87178291200.0D+00, 1307674368000.0D+00, &
    2.0922789888D+13, 3.5568742810D+14, 6.4023737057D+15, 1.2164510041D+17, &
    2.4329020082D+18, 1.5511210043D+25, 2.6525285981D+32 /)
  real ( kind = 8 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: nvec = (/ &
     0,  1,  2,  3, &
     4,  5,  6,  7, &
     8,  9, 10, 11, &
    12, 13, 14, 15, &
    16, 17, 18, 19, &
    20, 25, 30 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    fn = 0.0D+00
  else
    n = nvec(n_data)
    fn = fnvec(n_data)
  end if

  return
end
function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    Factorial2(N)
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
!    Output, real ( kind = 8 ) R8_FACTORIAL2, the value of the function.
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
function r8_nint ( x )

!*****************************************************************************80
!
!! R8_NINT returns the nearest integer to an R8.
!
!  Example:
!
!        X        R8_NINT
!
!      1.3         1
!      1.4         1
!      1.5         1 or 2
!      1.6         2
!      0.0         0
!     -0.7        -1
!     -1.1        -1
!     -1.6        -2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value.
!
!    Output, integer ( kind = 4 ) R8_NINT, the nearest integer to X.
!
  implicit none

  integer ( kind = 4 ) r8_nint
  integer ( kind = 4 ) s
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    s = -1
  else
    s = 1
  end if

  r8_nint = s * int ( abs ( x ) + 0.5D+00 )

  return
end
function r8_pi ( )

!*****************************************************************************80
!
!! R8_PI returns the value of Pi.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_PI, the value of Pi.
!
  implicit none

  real ( kind = 8 ) r8_pi

  r8_pi = 3.141592653589793D+00

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
!      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
!             = d/dX LN ( GAMMA(X) )
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
!    23 January 2008
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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nq
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

  x = xx
  w = abs ( x )
  aug = 0.0D+00
!
!  Check for valid arguments, then branch to appropriate algorithm.
!
  if ( xmax1 <= - x .or. w < xmin1 ) then

    if ( 0.0D+00 < x ) then
      r8_psi = - xinf
    else
      r8_psi = xinf
    end if

    return
  end if

  if ( x < 0.5D+00 ) then
!
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!
    if ( w <= xsmall ) then

      aug = - 1.0D+00 / x
!
!  Argument reduction for cotangent.
!
    else

      if ( x < 0.0D+00 ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - real ( int ( w ), kind = 8 )
      nq = int ( w * 4.0D+00 )
      w = 4.0D+00 * ( w - real ( nq, kind = 8 ) * 0.25D+00 )
!
!  W is now related to the fractional part of 4.0 * X.
!  Adjust argument to correspond to values in the first
!  quadrant and determine the sign.
!
      n = nq / 2

      if ( n + n /= nq ) then
        w = 1.0D+00 - w
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
        if ( z == 0.0D+00 ) then

          if ( 0.0D+00 < x ) then
            r8_psi = - xinf
          else
            r8_psi = xinf
          end if

          return
        end if

        aug = sgn * ( 4.0D+00 / tan ( z ) )

      else

        aug = sgn * ( 4.0D+00 * tan ( z ) )

      end if

    end if

    x = 1.0D+00 - x

  end if
!
!  0.5 <= X <= 3.0.
!
  if ( x <= 3.0D+00 ) then

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
    w = 1.0D+00 / ( x * x )
    den = w
    upper = p2(1) * w
    do i = 1, 5
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do
    aug = ( upper + p2(7) ) / ( den + q2(6) ) - 0.5D+00 / x + aug
  end if

  r8_psi = aug + log ( x )

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8poly_degree ( na, a, degree )

!*****************************************************************************80
!
!! R8POLY_DEGREE returns the degree of a polynomial in power sum form.
!
!  Discussion:
!
!    The power sum form of a polynomial is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x^(n-1) + a(n) * x^(n)
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) DEGREE, the degree of A.
!
  implicit none

  integer ( kind = 4 ) na

  real ( kind = 8 ) a(0:na)
  integer ( kind = 4 ) degree

  degree = na

  do while ( 0 < degree )

    if ( a(degree) /= 0.0D+00 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X^N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mag
  integer ( kind = 4 ) n2
  character plus_minus
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  call r8poly_degree ( n, a, n2 )

  if ( a(n2) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( 2 <= n2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine r8poly_val_horner ( n, c, x, cx )

!*****************************************************************************80
!
!! R8POLY_VAL_HORNER evaluates a polynomial using Horner's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Input, real ( kind = 8 ) C(0:N), the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) CX, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) cx
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  cx = c(n)
  do i = n - 1, 0, -1
    cx = cx * x + c(i)
  end do

  return
end
function r8poly_value ( n, a, x )

!*****************************************************************************80
!
!! R8POLY_VALUE evaluates an R8POLY
!
!  Discussion:
!
!    For sanity's sake, the value of N indicates the NUMBER of 
!    coefficients, or more precisely, the ORDER of the polynomial,
!    rather than the DEGREE of the polynomial.  The two quantities
!    differ by 1, but cause a great deal of confusion.
!
!    Given N and A, the form of the polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8poly_value
  real ( kind = 8 ) x

  r8poly_value = a(n)
  do i = n - 1, 1, -1
    r8poly_value = r8poly_value * x + a(i)
  end do

  return
end
subroutine r8vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do
    write ( *, '(a)' ) '  ......  ..............'
    i = n
    write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,2x,g14.6,2x,a)' ) i, a(i), '...more entries...'

  end if

  return
end
subroutine r8vec_uniform ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
function sec_deg ( angle )

!*****************************************************************************80
!
!! SEC_DEG returns the secant of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees.
!
!    Output, real ( kind = 8 ) SEC_DEG, the secant of the angle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00
  real ( kind = 8 ) sec_deg

  sec_deg = 1.0D+00 / sin ( degrees_to_radians * angle )

  return
end
subroutine sigma ( n, sigma_n )

!*****************************************************************************80
!
!! SIGMA returns the value of SIGMA(N), the divisor sum.
!
!  Discussion:
!
!    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
!
!    The formula is:
!
!      SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
!
!      SIGMA(P**K) = ( P**(K+1) - 1 ) / ( P - 1 ) if P is prime.
!
!  First values:
!
!     N  SIGMA(N)
!
!     1    1
!     2    3
!     3    4
!     4    7
!     5    6
!     6   12
!     7    8
!     8   15
!     9   13
!    10   18
!    11   12
!    12   28
!    13   14
!    14   24
!    15   24
!    16   31
!    17   18
!    18   39
!    19   20
!    20   42
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value to be analyzed.
!
!    Output, integer ( kind = 4 ) SIGMA_N, the value of SIGMA(N).  If N is
!    less than or equal to 0, SIGMA_N will be returned as 0.  If there is not
!    enough room for factoring N, SIGMA_N is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxfactor = 20

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) power(maxfactor)
  integer ( kind = 4 ) sigma_n

  if ( n <= 0 ) then
    sigma_n = 0
    return
  end if

  if ( n == 1 ) then
    sigma_n = 1
    return
  end if
!
!  Factor N.
!
  call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIGMA - Fatal error!'
    write ( *, '(a)' ) '  Not enough factorization space.'
    sigma_n = -1
    stop
  end if

  sigma_n = 1
  do i = 1, nfactor
    sigma_n = ( sigma_n * ( factor(i)**( power(i) + 1 ) - 1 ) ) &
      / ( factor(i) - 1 )
  end do

  return
end
subroutine sigma_values ( n_data, n, c )

!*****************************************************************************80
!
!! SIGMA_VALUES returns some values of the Sigma function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the Sigma function.
!
!    Output, integer ( kind = 4 ) C, the value of the Sigma function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 20

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
     1,    3,    4,    7,    6,   12,    8,   15,   13,   18, &
    72,  128,  255,  176,  576, 1170,  618,  984, 2232, 2340 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     30, 127, 128, 129, 210, 360, 617, 815, 816,1000 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
function sin_deg ( angle )

!*****************************************************************************80
!
!! SIN_DEG returns the sine of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees.
!
!    Output, real ( kind = 8 ) SIN_DEG, the sine of the angle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00
  real ( kind = 8 ) sin_deg

  sin_deg = sin ( degrees_to_radians * angle ) 

  return
end
function sin_power_int ( a, b, n )

!*****************************************************************************80
!
!! SIN_POWER_INT evaluates the sine power integral.
!
!  Discussion:
!
!    The function is defined by
!
!      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
!
!    The algorithm uses the following fact:
!
!      Integral sin^n ( t ) = (1/n) * (
!        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, integer ( kind = 4 ) N, the power of the sine function.
!
!    Output, real ( kind = 8 ) SIN_POWER_INT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mlo
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) sin_power_int
  real ( kind = 8 ) value

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIN_POWER_INT - Fatal error!'
    write ( *, '(a)' ) '  Power N < 0.'
    value = 0.0D+00
    stop
  end if

  sa = sin ( a )
  sb = sin ( b )
  ca = cos ( a )
  cb = cos ( b )

  if ( mod ( n, 2 ) == 0 ) then
    value = b - a
    mlo = 2
  else
    value = ca - cb
    mlo = 3
  end if

  do m = mlo, n, 2
    value = ( real ( m - 1, kind = 8 ) * value &
              + sa**(m-1) * ca - sb**(m-1) * cb ) &
      / real ( m, kind = 8 )
  end do

  sin_power_int = value

  return
end
subroutine sin_power_int_values ( n_data, a, b, n, fx )

!*****************************************************************************80
!
!! SIN_POWER_INT_VALUES returns some values of the sine power integral.
!
!  Discussion:
!
!    The function has the form
!
!      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin(T) )^N dt
!
!    In Mathematica, the function can be evaluated by:
!
!      Integrate [ ( Sin[x] )^n, { x, a, b } ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, integer ( kind = 4 ) N, the power.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
     0.10D+02, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.00D+00, &
     0.10D+01, &
     0.00D+00, &
     0.00D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
     0.20D+02, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.10D+01, &
     0.20D+01, &
     0.20D+01, &
     0.10D+01, &
     0.10D+01 /)
  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.10000000000000000000D+02, &
    0.45969769413186028260D+00, &
    0.27267564329357957615D+00, &
    0.17894056254885809051D+00, &
    0.12402556531520681830D+00, &
    0.88974396451575946519D-01, &
    0.90393123848149944133D+00, &
    0.81495684202992349481D+00, &
    0.21887522421729849008D-01, &
    0.17023439374069324596D-01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0, &
     1, &
     2, &
     3, &
     4, &
     5, &
     5, &
     5, &
    10, &
    11 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    n = 0
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine slice ( dim_num, slice_num, piece_num )

!*****************************************************************************80
!
!! SLICE: maximum number of pieces created by a given number of slices.
!
!  Discussion:
!
!    If we imagine slicing a pizza, each slice produce more pieces.  
!    The position of the slice affects the number of pieces created, but there
!    is a maximum.  
!
!    This function determines the maximum number of pieces created by a given
!    number of slices, applied to a space of a given dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Banks,
!    Slicing Pizzas, Racing Turtles, and Further Adventures in 
!    Applied Mathematics,
!    Princeton, 1999,
!    ISBN13: 9780691059471,
!    LC: QA93.B358.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SLICE_NUM, the number of slices.
!
!    Input, real ( kind = 8 ) PIECE_NUM, the maximum number of pieces that can
!    be created by the given number of slices applied in the given dimension.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) j
  integer ( kind = 4 ) piece_num
  integer ( kind = 4 ) slice_num

  piece_num = 0
  do j = 0, min ( dim_num, slice_num )
    piece_num = piece_num + i4_choose ( slice_num, j )
  end do

  return
end
subroutine spherical_harmonic ( l, m, theta, phi, c, s )

!*****************************************************************************80
!
!! SPHERICAL_HARMONIC evaluates spherical harmonic functions.
!
!  Discussion:
!
!    The spherical harmonic function Y(L,M,THETA,PHI,X) is the
!    angular part of the solution to Laplace's equation in spherical
!    coordinates.
!
!    Y(L,M,THETA,PHI,X) is related to the associated Legendre
!    function as follows:
!
!      Y(L,M,THETA,PHI,X) = FACTOR * P(L,M,cos(THETA)) * exp ( i * M * PHI )
!
!    Here, FACTOR is a normalization factor:
!
!      FACTOR = sqrt ( ( 2 * L + 1 ) * ( L - M )! / ( 4 * PI * ( L + M )! ) )
!
!    In Mathematica, a spherical harmonic function can be evaluated by
!
!      SphericalHarmonicY [ l, m, theta, phi ]
!
!    Note that notational tradition in physics requires that THETA
!    and PHI represent the reverse of what they would normally mean
!    in mathematical notation; that is, THETA goes up and down, and
!    PHI goes around.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 2005
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
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the first index of the spherical harmonic
!    function.  Normally, 0 <= L.
!
!    Input, integer ( kind = 4 ) M, the second index of the spherical harmonic 
!    function.  Normally, -L <= M <= L.
!
!    Input, real ( kind = 8 ) THETA, the polar or latitudinal angle, for which
!    0 <= THETA <= PI.
!
!    Input, real ( kind = 8 ) PHI, the longitudinal angle, for which
!    0 <= PHI <= 2*PI.
!
!    Output, real ( kind = 8 ) C(0:L), S(0:L), the real and imaginary
!    parts of the functions Y(L,0:L,THETA,PHI).
!
  implicit none

  integer ( kind = 4 ) l

  real ( kind = 8 ) c(0:l)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_abs
  real ( kind = 8 ) phi
  real ( kind = 8 ) plm(0:l)
  real ( kind = 8 ) s(0:l)
  real ( kind = 8 ) theta

  m_abs = abs ( m )

  call legendre_associated_normalized ( l, m_abs, cos ( theta ), plm )

  c(0:l) = plm(0:l) * cos ( real ( m, kind = 8 ) * phi )
  s(0:l) = plm(0:l) * sin ( real ( m, kind = 8 ) * phi )

  if ( m < 0 ) then
    c(0:l) = - c(0:l)
    s(0:l) = - s(0:l)
  end if

  return
end
subroutine spherical_harmonic_values ( n_data, l, m, theta, phi, yr, yi )

!*****************************************************************************80
!
!! SPHERICAL_HARMONIC_VALUES returns values of spherical harmonic functions.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by
!
!      SphericalHarmonicY [ l, m, theta, phi ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 2005
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
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and
!    N_DATA is set to the index of the test data.  On each subsequent
!    call, N_DATA is incremented and that test data is returned.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) L, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) THETA, PHI, the arguments of the function.
!
!    Output, real ( kind = 8 ) YR, YI, the real and imaginary parts of
!    the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 20

  integer ( kind = 4 ) l
  integer ( kind = 4 ), save, dimension ( nmax ) :: l_vec = (/ &
     0,  1,  2,  &
     3,  4,  5,  &
     5,  5,  5,  &
     5,  4,  4,  &
     4,  4,  4,  &
     3,  3,  3,  &
     3,  3 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( nmax ) :: m_vec = (/ &
     0,  0,  1,  &
     2,  3,  5,  &
     4,  3,  2,  &
     1,  2,  2,  &
     2,  2,  2,  &
    -1, -1, -1,  &
    -1, -1 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) phi
  real ( kind = 8 ), save, dimension ( nmax ) :: phi_vec = (/ &
    0.1047197551196598D+01, 0.1047197551196598D+01, 0.1047197551196598D+01, &
    0.1047197551196598D+01, 0.1047197551196598D+01, 0.6283185307179586D+00, &
    0.6283185307179586D+00, 0.6283185307179586D+00, 0.6283185307179586D+00, &
    0.6283185307179586D+00, 0.7853981633974483D+00, 0.7853981633974483D+00, &
    0.7853981633974483D+00, 0.7853981633974483D+00, 0.7853981633974483D+00, &
    0.4487989505128276D+00, 0.8975979010256552D+00, 0.1346396851538483D+01, &
    0.1795195802051310D+01, 0.2243994752564138D+01 /)
  real ( kind = 8 ) theta
  real ( kind = 8 ), save, dimension ( nmax ) :: theta_vec = (/ &
    0.5235987755982989D+00, 0.5235987755982989D+00, 0.5235987755982989D+00, &
    0.5235987755982989D+00, 0.5235987755982989D+00, 0.2617993877991494D+00, &
    0.2617993877991494D+00, 0.2617993877991494D+00, 0.2617993877991494D+00, &
    0.2617993877991494D+00, 0.6283185307179586D+00, 0.1884955592153876D+01, &
    0.3141592653589793D+01, 0.4398229715025711D+01, 0.5654866776461628D+01, &
    0.3926990816987242D+00, 0.3926990816987242D+00, 0.3926990816987242D+00, &
    0.3926990816987242D+00, 0.3926990816987242D+00 /)
  real ( kind = 8 ) yi
  real ( kind = 8 ), save, dimension ( nmax ) :: yi_vec = (/ &
    0.0000000000000000D+00,  0.0000000000000000D+00, -0.2897056515173922D+00, &
    0.1916222768312404D+00,  0.0000000000000000D+00,  0.0000000000000000D+00, &
    0.3739289485283311D-02, -0.4219517552320796D-01,  0.1876264225575173D+00, &
   -0.3029973424491321D+00,  0.4139385503112256D+00, -0.1003229830187463D+00, &
    0.0000000000000000D+00, -0.1003229830187463D+00,  0.4139385503112256D+00, &
   -0.1753512375142586D+00, -0.3159720118970196D+00, -0.3940106541811563D+00, &
   -0.3940106541811563D+00, -0.3159720118970196D+00 /)
  real ( kind = 8 ) yr
  real ( kind = 8 ), save, dimension ( nmax ) :: yr_vec = (/ &
   0.2820947917738781D+00,  0.4231421876608172D+00, -0.1672616358893223D+00, &
  -0.1106331731112457D+00,  0.1354974113737760D+00,  0.5390423109043568D-03, &
  -0.5146690442951909D-02,  0.1371004361349490D-01,  0.6096352022265540D-01, &
  -0.4170400640977983D+00,  0.0000000000000000D+00,  0.0000000000000000D+00, &
   0.0000000000000000D+00,  0.0000000000000000D+00,  0.0000000000000000D+00, &
   0.3641205966137958D+00,  0.2519792711195075D+00,  0.8993036065704300D-01, &
  -0.8993036065704300D-01, -0.2519792711195075D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    l = 0
    m = 0
    theta = 0.0D+00
    phi = 0.0D+00
    yr = 0.0D+00
    yi = 0.0D+00
  else
    l = l_vec(n_data)
    m = m_vec(n_data)
    theta = theta_vec(n_data)
    phi = phi_vec(n_data)
    yr = yr_vec(n_data)
    yi = yi_vec(n_data)
  end if

  return
end
subroutine stirling1 ( n, m, s1 )

!*****************************************************************************80
!
!! STIRLING1 computes the Stirling numbers of the first kind.
!
!  Discussion:
!
!    The absolute value of the Stirling number S1(N,M) gives the number
!    of permutations on N objects having exactly M cycles, while the
!    sign of the Stirling number records the sign (odd or even) of
!    the permutations.  For example, there are six permutations on 3 objects:
!
!      A B C   3 cycles (A) (B) (C)
!      A C B   2 cycles (A) (BC)
!      B A C   2 cycles (AB) (C)
!      B C A   1 cycle  (ABC)
!      C A B   1 cycle  (ABC)
!      C B A   2 cycles (AC) (B)
!
!    There are 
!
!      2 permutations with 1 cycle, and S1(3,1) = 2
!      3 permutations with 2 cycles, and S1(3,2) = -3,
!      1 permutation with 3 cycles, and S1(3,3) = 1.
!
!    Since there are N! permutations of N objects, the sum of the absolute 
!    values of the Stirling numbers in a given row, 
!
!      sum ( 1 <= I <= N ) abs ( S1(N,I) ) = N!
!
!  First terms:
!
!    N/M:  1     2      3     4     5    6    7    8
!
!    1     1     0      0     0     0    0    0    0
!    2    -1     1      0     0     0    0    0    0
!    3     2    -3      1     0     0    0    0    0
!    4    -6    11     -6     1     0    0    0    0
!    5    24   -50     35   -10     1    0    0    0
!    6  -120   274   -225    85   -15    1    0    0
!    7   720 -1764   1624  -735   175  -21    1    0
!    8 -5040 13068 -13132  6769 -1960  322  -28    1
!
!  Recursion:
!
!    S1(N,1) = (-1)^(N-1) * (N-1)! for all N.
!    S1(I,I) = 1 for all I.
!    S1(I,J) = 0 if I < J.
!
!    S1(N,M) = S1(N-1,M-1) - (N-1) * S1(N-1,M)
!
!  Properties:
!
!    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
!
!    X_N = sum ( 0 <= K <= N ) S1(N,K) X^K
!    where X_N is the falling factorial function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows of the table.
!
!    Input, integer ( kind = 4 ) M, the number of columns of the table.
!
!    Output, integer ( kind = 4 ) S1(N,M), the Stirling numbers of the 
!    first kind.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s1(n,m)

  if ( n <= 0 ) then
    return
  end if

  if ( m <= 0 ) then
    return
  end if

  s1(1,1) = 1
  s1(1,2:m) = 0
 
  do i = 2, n

    s1(i,1) = - ( i - 1 ) * s1(i-1,1)

    do j = 2, m
      s1(i,j) = s1(i-1,j-1) - ( i - 1 ) * s1(i-1,j)
    end do

  end do
 
  return
end
subroutine stirling2 ( n, m, s2 )

!*****************************************************************************80
!
!! STIRLING2 computes the Stirling numbers of the second kind.
!
!  Discussion:
!
!    S2(N,M) represents the number of distinct partitions of N elements
!    into M nonempty sets.  For a fixed N, the sum of the Stirling
!    numbers S2(N,M) is represented by B(N), called "Bell's number",
!    and represents the number of distinct partitions of N elements.
!
!    For example, with 4 objects, there are:
!
!    1 partition into 1 set:
!
!      (A,B,C,D)
!
!    7 partitions into 2 sets:
!
!      (A,B,C) (D)
!      (A,B,D) (C)
!      (A,C,D) (B)
!      (A) (B,C,D)
!      (A,B) (C,D)
!      (A,C) (B,D)
!      (A,D) (B,C)
!
!    6 partitions into 3 sets:
!
!      (A,B) (C) (D)
!      (A) (B,C) (D)
!      (A) (B) (C,D)
!      (A,C) (B) (D)
!      (A,D) (B) (C)
!      (A) (B,D) (C)
!
!    1 partition into 4 sets:
!
!      (A) (B) (C) (D)
!
!    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
!
!    The Stirling numbers of the second kind S(N,1:N) are the coefficients of
!    the Bell polynomial B(N,X):
!
!      B(0,X) = 1
!      B(N,X) = sum ( 1 <= M <= N ) S(N,M) * X^M
!
!  First terms:
!
!    N/M: 1    2    3    4    5    6    7    8
!
!    1    1    0    0    0    0    0    0    0
!    2    1    1    0    0    0    0    0    0
!    3    1    3    1    0    0    0    0    0
!    4    1    7    6    1    0    0    0    0
!    5    1   15   25   10    1    0    0    0
!    6    1   31   90   65   15    1    0    0
!    7    1   63  301  350  140   21    1    0
!    8    1  127  966 1701 1050  266   28    1
!
!  Recursion:
!
!    S2(N,1) = 1 for all N.
!    S2(I,I) = 1 for all I.
!    S2(I,J) = 0 if I < J.
!
!    S2(N,M) = M * S2(N-1,M) + S2(N-1,M-1)
!
!  Properties:
!
!    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
!
!    X^N = sum ( 0 <= K <= N ) S2(N,K) X_K
!    where X_K is the falling factorial function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows of the table.
!
!    Input, integer ( kind = 4 ) M, the number of columns of the table.
!
!    Output, integer ( kind = 4 ) S2(N,M), the Stirling numbers of the 
!    second kind.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s2(n,m)

  if ( n <= 0 ) then
    return
  end if

  if ( m <= 0 ) then
    return
  end if

  s2(1,1) = 1
  s2(1,2:m) = 0
 
  do i = 2, n

    s2(i,1) = 1

    do j = 2, m
      s2(i,j) = j * s2(i-1,j) + s2(i-1,j-1)
    end do

  end do
 
  return
end
function tan_deg ( angle )

!*****************************************************************************80
!
!! TAN_DEG returns the tangent of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in degrees.
!
!    Output, real ( kind = 8 ) TAN_DEG, the tangent of the angle.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00
  real ( kind = 8 ) tan_deg

  tan_deg = sin ( degrees_to_radians * angle ) &
          / cos ( degrees_to_radians * angle )

  return
end
subroutine tau ( n, taun )

!*****************************************************************************80
!
!! TAU returns the value of TAU(N), the number of distinct divisors of N.
!
!  Discussion:
!
!    TAU(N) is the number of distinct divisors of N, including 1 and N.
!
!    If the prime factorization of N is
!
!      N = P1**E1 * P2**E2 * ... * PM**EM,
!
!    then
!
!      TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
!
!    One consequence of this fact is that TAU is odd if and only
!    if N is a perfect square.
!
!  First values:
!
!     N   TAU(N)
!
!     1    1
!     2    2
!     3    2
!     4    3
!     5    2
!     6    4
!     7    2
!     8    4
!     9    3
!    10    4
!    11    2
!    12    6
!    13    2
!    14    4
!    15    4
!    16    5
!    17    2
!    18    6
!    19    2
!    20    6
!    21    4
!    22    4
!    23    2
!    24    8
!    25    3
!    26    4
!    27    4
!    28    6
!    29    2
!    30    8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value to be analyzed.  N must be 1 or
!    greater.
!
!    Output, integer ( kind = 4 ) TAUN, the value of TAU(N).  But if N is 0 or
!    less, TAUN is returned as 0, a nonsense value.  If there is
!    not enough room for factoring, TAUN is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxfactor = 20

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) power(maxfactor)
  integer ( kind = 4 ) taun

  if ( n <= 0 ) then
    taun = 0
    return
  end if

  if ( n == 1 ) then
    taun = 1
    return
  end if
!
!  Factor N.
!
  call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TAU - Fatal error!'
    write ( *, '(a)' ) '  Not enough factorization space.'
    taun = -1
    stop
  end if

  taun = 1
  do i = 1, nfactor
    taun = taun * ( power(i) + 1 )
  end do

  return
end
subroutine tau_values ( n_data, n, c )

!*****************************************************************************80
!
!! TAU_VALUES returns some values of the Tau function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 February 2003
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
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When 
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the argument of the Tau function.
!
!    Output, integer ( kind = 4 ) C, the value of the Tau function.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 20

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
    1,  2,  2,  3,  2,  4,  2,  4,  3,  4, &
    2, 12, 12,  4, 18, 24,  2,  8, 14, 28 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
     23,  72, 126, 226, 300, 480, 521, 610, 832, 960 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
function tetrahedron_num ( n )

!*****************************************************************************80
!
!! TETRAHEDRON_NUM returns the N-th tetrahedral number.
!
!  Discussion:
!
!    The N-th tetrahedral number T3(N) is formed by the sum of the first
!    N triangular numbers:
!
!      T3(N) = sum ( 1 <= I <= N ) T2(I)
!            = sum ( 1 <= I <= N ) sum ( 1 <= J < I ) J
!
!    By convention, T3(0) = 0.
!
!    The formula is:
!
!      T3(N) = ( N * ( N + 1 ) * ( N + 2 ) ) / 6
!
!  First Values:
!
!     0
!     1
!     4
!    10
!    20
!    35
!    56
!    84
!   120
!   165
!   220
!   275
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired number, which 
!    must be at least 0.
!
!    Output, integer ( kind = 4 ) TETRAHEDRON_NUM, the N-th tetrahedron number.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) tetrahedron_num

  tetrahedron_num = ( n * ( n + 1 ) * ( n + 2 ) ) / 6

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
function triangle_num ( n )

!*****************************************************************************80
!
!! TRIANGLE_NUM returns the N-th triangular number.
!
!  Discussion:
!
!    The N-th triangular number T(N) is formed by the sum of the first
!    N integers:
!
!      T(N) = sum ( 1 <= I <= N ) I
!
!    By convention, T(0) = 0.
!
!    T(N) can be computed quickly by the formula:
!
!      T(N) = ( N * ( N + 1 ) ) / 2
!
!  First Values:
!
!     0
!     1
!     3
!     6
!    10
!    15
!    21
!    28
!    36
!    45
!    55
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired number, 
!    which must be at least 0.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NUM, the N-th triangular number.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) triangle_num

  triangle_num = ( n * ( n + 1 ) ) / 2

  return
end
subroutine triangle_to_i4 ( i, j, k )

!*****************************************************************************80
!
!! TRIANGLE_TO_I4 converts a triangular coordinate to an integer.
!
!  Discussion:
!
!    Triangular coordinates are handy when storing a naturally triangular
!    array (such as the lower half of a matrix) in a linear array.
!
!    Thus, for example, we might consider storing 
!
!    (1,1)
!    (2,1) (2,2)
!    (3,1) (3,2) (3,3)
!    (4,1) (4,2) (4,3) (4,4)
!
!    as the linear array
!
!    (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1) (4,2) (4,3) (4,4)    
!
!    Here, the quantities in parenthesis represent the natural row and
!    column indices of a single number when stored in a rectangular array.
!
!    Thus, our goal is, given the row I and column J of the data,
!    to produce the value K which indicates its position in the linear
!    array.
!
!    The triangular numbers are the indices associated with the
!    diagonal elements of the original array, T(1,1), T(2,2), T(3,3)
!    and so on.
!
!    The formula is:
!
!      K = J + ( (I-1) * I ) / 2
!
!  First Values:
!
!     I  J  K
!
!     0  0  0
!     1  1  1
!     2  1  2
!     2  2  3
!     3  1  4
!     3  2  5
!     3  3  6
!     4  1  7
!     4  2  8
!     4  3  9
!     4  4 10
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
!    Input, integer ( kind = 4 ) I, J, the row and column indices.  I and J must
!    be nonnegative, and J must not be greater than I.
!
!    Output, integer ( kind = 4 ) K, the linear index of the (I,J) element.

  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_TO_I4 - Fatal error!'
    write ( *, '(a)' ) '  I < 0.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  else if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_TO_I4 - Fatal error!'
    write ( *, '(a)' ) '  J < 0.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  else if ( i < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_TO_I4 - Fatal error!'
    write ( *, '(a)' ) '  I < J.'
    write ( *, '(a,i8)' ) '  I = ', i
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  end if

  k = j + ( ( i - 1 ) * i ) / 2

  return
end
recursive function v_hofstadter ( n ) result ( value )

!*****************************************************************************80
!
!! V_HOFSTADTER computes the Hofstadter V sequence.
!
!  Discussion:
!
!    V(N) = 0                          if N = 0
!         = 1                          if 1 <= N <= 4
!         = V(N-V(N-1)) + V(N-V(N-4)), otherwise.
!
!    V(N) is defined for all nonnegative integers.
!
!  Table:
!
!     N  V(N)
!    --  ----
!
!     0     0
!     1     1
!     2     1
!     3     1
!     4     1
!     5     2
!     6     3
!     7     4
!     8     5
!     9     5
!    10     6
!    11     6
!    12     7
!    13     8
!    14     8
!    15     9
!    16     9
!    17    10
!    18    11
!    19    11
!    20    11
!    21    12
!    22    12
!    23    13
!    24    14
!    25    14
!    26    15
!    27    15
!    28    16
!    29    17
!    30    17
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the function.
!
!    Output, integer ( kind = 4 ) V_HOFSTADTER, the value of the function.
!
  implicit none

  integer ( kind = 4 ) arg1
  integer ( kind = 4 ) arg2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  if ( n <= 0 ) then

    value = 0

  else if ( n <= 4 ) then

    value = 1

  else

    arg1 = n - 1
    arg2 = n - 4

    value = v_hofstadter ( n - v_hofstadter ( arg1 ) ) &
          + v_hofstadter ( n - v_hofstadter ( arg2 ) )

  end if

  return
end
subroutine vibonacci ( n, seed, v )

!*****************************************************************************80
!
!! VIBONACCI computes the first N Vibonacci numbers.
!
!  Discussion:
!
!    The "Vibonacci numbers" are a generalization of the Fibonacci numbers:
!      V(N+1) = +/- V(N) +/- V(N-1)
!    where the signs are chosen randomly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Brian Hayes,
!    The Vibonacci Numbers,
!    American Scientist,
!    July-August 1999, Volume 87, Number 4.
!
!    Divakar Viswanath,
!    Random Fibonacci sequences and the number 1.13198824,
!    Mathematics of Computation,
!    1998.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest number to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, integer ( kind = 4 ) V(N), the first N Vibonacci numbers.  By 
!    convention, V(1) and V(2) are taken to be 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) v(n)

  if ( n <= 0 ) then
    return
  end if

  v(1) = 1

  if ( n <= 1 ) then
    return
  end if

  v(2) = 1

  do i = 3, n
    
    j = i4_uniform ( 0, 1, seed )

    if ( j == 0 ) then
      s1 = -1
    else
      s1 = +1
    end if

    j = i4_uniform ( 0, 1, seed )

    if ( j == 0 ) then
      s2 = -1
    else
      s2 = +1
    end if

    v(i) = s1 * v(i-1) + s2 * v(i-2)

  end do
 
  return
end
subroutine zeckendorf ( n, m_max, m, i_list, f_list )

!*****************************************************************************80
!
!! ZECKENDORF produces the Zeckendorf decomposition of a positive integer.
!
!  Discussion:
!
!    Zeckendorf proved that every positive integer can be represented
!    uniquely as the sum of non-consecutive Fibonacci numbers.
!
!    N = sum ( 1 <= I <= M ) F_LIST(I)
!
!  Example:
!
!     N    Decomposition
!
!    50    34 + 13 + 3
!    51    34 + 13 + 3 + 1
!    52    34 + 13 + 5
!    53    34 + 13 + 5 + 1
!    54    34 + 13 + 5 + 2
!    55    55
!    56    55 + 1
!    57    55 + 2
!    58    55 + 3
!    59    55 + 3 + 1
!    60    55 + 5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the positive integer to be decomposed.
!
!    Input, integer ( kind = 4 ) M_MAX, the maximum dimension of I_LIST 
!    and F_LIST.
!
!    Output, integer ( kind = 4 ) M, the number of parts in the decomposition.
!
!    Output, integer ( kind = 4 ) I_LIST(M_MAX), contains in entries 1 
!    through M the index of the Fibonacci numbers in the decomposition.
!
!    Output, integer ( kind = 4 ) F_LIST(M_MAX), contains in entries 1 
!    through M the value of the Fibonacci numbers in the decomposition.
!
  implicit none

  integer ( kind = 4 ) m_max

  integer ( kind = 4 ) f
  integer ( kind = 4 ) f_list(m_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_list(m_max)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_copy

  m = 0

  n_copy = n
!
!  Extract a sequence of Fibonacci numbers.
!
  do while ( 0 < n_copy .and. m < m_max ) 
    call fibonacci_floor ( n_copy, f, i )
    m = m + 1
    i_list(m) = i
    n_copy = n_copy - f
  end do
!
!  Replace any pair of consecutive indices ( I, I-1 ) by I+1.
!
  do i = m, 2, -1

    if ( i_list(i-1) == i_list(i) + 1 ) then
      i_list(i-1) = i_list(i-1) + 1
      i_list(i:m-1) = i_list(i+1:m)
      i_list(m) = 0
      m = m - 1
    end if

  end do
!
!  Fill in the actual values of the Fibonacci numbers.
!
  do i = 1, m
    call fibonacci_direct ( i_list(i), f_list(i) )
  end do

  return
end
subroutine zernike_poly ( m, n, rho, z )

!*****************************************************************************80
!
!! ZERNIKE_POLY evaluates a Zernike polynomial at RHO.
!
!  Discussion:
!
!    This routine uses the facts that:
!
!    *) R^M_N = 0 if M < 0, or N < 0, or N < M.
!    *) R^M_M = RHO^M
!    *) R^M_N = 0 if mod ( N - M ) = 1.
!
!    and the recursion:
!
!    R^M_(N+2) = A * [ ( B * RHO^2 - C ) * R^M_N - D * R^M_(N-2) ]
!
!    where 
!
!    A = ( N + 2 ) / ( ( N + 2 )^2 - M^2 )
!    B = 4 * ( N + 1 )
!    C = ( N + M )^2 / N + ( N - M + 2 )^2 / ( N + 2 )
!    D = ( N^2 - M^2 ) / N
!
!    I wish I could clean up the recursion in the code, but for
!    now, I have to treat the case M = 0 specially.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the upper index.
!
!    Input, integer ( kind = 4 ) N, the lower index.
!
!    Input, real ( kind = 8 ) RHO, the radial coordinate.
!
!    Output, real ( kind = 8 ) Z, the value of the Zernike 
!    polynomial R^M_N at the point RHO.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nn
  real ( kind = 8 ) rho
  real ( kind = 8 ) z
  real ( kind = 8 ) zm2
  real ( kind = 8 ) zp2
!
!  Do checks.
!
  if ( m < 0 ) then
    z = 0.0D+00
    return
  end if

  if ( n < 0 ) then
    z = 0.0D+00
    return
  end if

  if ( n < m ) then
    z = 0.0D+00
    return
  end if

  if ( mod ( n - m, 2 ) == 1 ) then
    z = 0.0D+00
    return
  end if

  zm2 = 0.0D+00
  z = rho**m

  if ( m == 0 ) then

    if ( n == 0 ) then
      return
    end if

    zm2 = z
    z = 2.0D+00 * rho * rho - 1.0D+00

    do nn = m+2, n-2, 2

      a = real ( nn + 2, kind = 8 ) / real ( ( nn + 2 )**2 - m**2, kind = 8 )
      b = real ( 4 * ( nn + 1 ), kind = 8 )
      c = real ( ( nn + m )**2, kind = 8 ) / real ( nn, kind = 8 ) &
        + real ( ( nn - m + 2 )**2, kind = 8 ) / real ( nn + 2, kind = 8 )
      d = real ( nn**2 - m**2, kind = 8 ) / real ( nn, kind = 8 )

      zp2 = a * ( ( b * rho**2 - c ) * z - d * zm2 ) 
      zm2 = z
      z = zp2

    end do

  else

    do nn = m, n-2, 2

      a = real ( nn + 2, kind = 8 ) / real ( ( nn + 2 )**2 - m**2, kind = 8 )
      b = real ( 4 * ( nn + 1 ), kind = 8 )
      c = real ( ( nn + m )**2, kind = 8 ) / real ( nn, kind = 8 ) &
        + real ( ( nn - m + 2 )**2, kind = 8 ) / real ( nn + 2, kind = 8 )
      d = real ( nn**2 - m**2, kind = 8 ) / real ( nn, kind = 8 )

      zp2 = a * ( ( b * rho**2 - c ) * z - d * zm2 ) 
      zm2 = z
      z = zp2

    end do

  end if

  return
end
subroutine zernike_poly_coef ( m, n, c )

!*****************************************************************************80
!
!! ZERNIKE_POLY_COEF: coefficients of a Zernike polynomial.
!
!  Discussion:
!
!    With our coefficients stored in C(0:N), the
!    radial function R^M_N(RHO) is given by
!
!      R^M_N(RHO) = C(0) 
!                 + C(1) * RHO
!                 + C(2) * RHO^2
!                 + ...
!                 + C(N) * RHO^N
!
!    and the odd and even Zernike polynomials are
!
!      Z^M_N(RHO,PHI,odd)  = R^M_N(RHO) * sin(PHI)
!      Z^M_N(RHO,PHI,even) = R^M_N(RHO) * cos(PHI)
!
!    The first few "interesting" values of R are:
!
!    R^0_0 = 1
!
!    R^1_1 = RHO
!
!    R^0_2 = 2 * RHO^2 - 1
!    R^2_2 =     RHO^2
!
!    R^1_3 = 3 * RHO^3 - 2 * RHO
!    R^3_3 =     RHO^3
!
!    R^0_4 = 6 * RHO^4 - 6 * RHO^2 + 1
!    R^2_4 = 4 * RHO^4 - 3 * RHO^2
!    R^4_4 =     RHO^4
!
!    R^1_5 = 10 * RHO^5 - 12 * RHO^3 + 3 * RHO
!    R^3_5 =  5 * RHO^5 -  4 * RHO^3
!    R^5_5 =      RHO^5
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
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the parameters of the polynomial.
!    Normally, 0 <= M <= N and 0 <= N.
!  
!    Output, real ( kind = 8 ) C(0:N), the coefficients of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nm_minus
  integer ( kind = 4 ) nm_plus
  real ( kind = 8 ) r8_choose

  c(0:n) = 0.0D+00

  if ( n < 0 ) then
    return
  end if

  if ( m < 0 ) then
    return
  end if
      
  if ( n < m ) then
    return
  end if

  if ( mod ( n - m, 2 ) == 1 ) then
    return
  end if

  nm_plus = ( m + n ) / 2
  nm_minus = ( n - m ) / 2

  c(n) = r8_choose ( n, nm_plus )

  do l = 0, nm_minus - 1

    c(n-2*l-2) = - real ( ( nm_plus - l ) * ( nm_minus - l ), kind = 8 ) &
      * c(n-2*l) / real ( ( n - l ) * ( l + 1 ), kind = 8 )

  end do

  return
end
function zeta ( p )

!*****************************************************************************80
!
!! ZETA estimates the Riemann Zeta function.
!
!  Discussion:
!
!    For 1 < P, the Riemann Zeta function is defined as:
!
!      ZETA ( P ) = Sum ( 1 <= N < +oo ) 1 / N^P
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2004
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
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the power to which the integers are raised.
!    P must be greater than 1.
!
!    Output, real ( kind = 8 ) ZETA, an approximation to the Riemann 
!    Zeta function.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) p
  real ( kind = 8 ) total
  real ( kind = 8 ) total_old
  real ( kind = 8 ) zeta

  if ( p <= 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ZETA - Fatal error!'
    write ( *, '(a)' ) '  Exponent P <= 1.0.'
    zeta = -1.0D+00
    stop
  end if

  total = 0.0D+00
  n = 0

  do

    n = n + 1
    total_old = total
    total = total + 1.0D+00 / ( real ( n, kind = 8 ) )**p

    if ( total <= total_old .or. 1000 <= n ) then
      exit
    end if

  end do

  zeta = total

  return
end
subroutine zeta_values ( n_data, n, z )

!*****************************************************************************80
!
!! ZETA_VALUES returns some values of the Riemann Zeta function.
!
!  Discussion:
!
!    ZETA(N) = sum ( 1 <= I < +oo ) 1 / I^N
!
!    In Mathematica, the function can be evaluated by:
!
!      Zeta[n]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 2004
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1, 
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the argument of the Zeta function.
!
!    Output, real ( kind = 8 ) Z, the value of the Zeta function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     2, &
     3, &
     4, &
     5, &
     6, &
     7, &
     8, &
     9, &
    10, &
    11, &
    12, &
    16, &
    20, &
    30, &
    40 /)
  real ( kind = 8 ) z
  real ( kind = 8 ), save, dimension ( n_max ) :: zeta_vec = (/ &
    0.164493406684822643647D+01, &
    0.120205690315959428540D+01, &
    0.108232323371113819152D+01, &
    0.103692775514336992633D+01, &
    0.101734306198444913971D+01, &
    0.100834927738192282684D+01, &
    0.100407735619794433939D+01, &
    0.100200839292608221442D+01, &
    0.100099457512781808534D+01, &
    0.100049418860411946456D+01, &
    0.100024608655330804830D+01, &
    0.100001528225940865187D+01, &
    0.100000095396203387280D+01, &
    0.100000000093132743242D+01, &
    0.100000000000090949478D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    z = 0.0D+00
  else
    n = n_vec(n_data)
    z = zeta_vec(n_data)
  end if

  return
end
