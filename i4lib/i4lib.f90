function i4_bit_hi1 ( n )

!*****************************************************************************80
!
!! I4_BIT_HI1 returns the position of the high 1 bit base 2 in an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!       N    Binary    Hi 1
!    ----    --------  ----
!       0           0     0
!       1           1     1
!       2          10     2
!       3          11     2
!       4         100     3
!       5         101     3
!       6         110     3
!       7         111     3
!       8        1000     4
!       9        1001     4
!      10        1010     4
!      11        1011     4
!      12        1100     4
!      13        1101     4
!      14        1110     4
!      15        1111     4
!      16       10000     5
!      17       10001     5
!    1023  1111111111    10
!    1024 10000000000    11
!    1025 10000000001    11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be measured.
!    N should be nonnegative.  If N is nonpositive, the function
!    will always be 0.
!
!    Output, integer ( kind = 4 ) I4_BIT_HI1, the position of the highest bit.
!
  implicit none

  integer ( kind = 4 ) bit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_bit_hi1
  integer ( kind = 4 ) n

  i = n
  bit = 0

  do

    if ( i <= 0 ) then
      exit
    end if

    bit = bit + 1
    i = i / 2

  end do

  i4_bit_hi1 = bit

  return
end
function i4_bit_lo0 ( n )

!*****************************************************************************80
!
!! I4_BIT_LO0 returns the position of the low 0 bit base 2 in an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!       N    Binary    Lo 0
!    ----    --------  ----
!       0           0     1
!       1           1     2
!       2          10     1
!       3          11     3
!       4         100     1
!       5         101     2
!       6         110     1
!       7         111     4
!       8        1000     1
!       9        1001     2
!      10        1010     1
!      11        1011     3
!      12        1100     1
!      13        1101     2
!      14        1110     1
!      15        1111     5
!      16       10000     1
!      17       10001     2
!    1023  1111111111     1
!    1024 10000000000     1
!    1025 10000000001     1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be measured.
!    N should be nonnegative.
!
!    Output, integer ( kind = 4 ) I4_BIT_LO0, the position of the low 1 bit.
!
  implicit none

  integer ( kind = 4 ) bit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_bit_lo0
  integer ( kind = 4 ) n

  bit = 0
  i = n

  do

    bit = bit + 1
    i2 = i / 2

    if ( i == 2 * i2 ) then
      exit
    end if

    i = i2

  end do

  i4_bit_lo0 = bit

  return
end
function i4_bit_lo1 ( n )

!*****************************************************************************80
!
!! I4_BIT_LO1 returns the position of the low 1 bit base 2 in an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!       N    Binary    Lo 1
!    ----    --------  ----
!       0           0     0
!       1           1     1
!       2          10     2
!       3          11     1
!       4         100     3
!       5         101     1
!       6         110     2
!       7         111     1
!       8        1000     4
!       9        1001     1
!      10        1010     2
!      11        1011     1
!      12        1100     3
!      13        1101     1
!      14        1110     2
!      15        1111     1
!      16       10000     5
!      17       10001     1
!    1023  1111111111     1
!    1024 10000000000    11
!    1025 10000000001     1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be measured.
!    N should be nonnegative.
!
!    Output, integer ( kind = 4 ) I4_BIT_LO1, the position of the low 1 bit.
!
  implicit none

  integer ( kind = 4 ) bit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_bit_lo1
  integer ( kind = 4 ) n

  bit = 0
  i = n

  do

    bit = bit + 1
    i2 = i / 2

    if ( i /= 2 * i2 ) then
      exit
    end if

    i = i2

  end do

  i4_bit_lo1 = bit

  return
end
function i4_bit_reverse ( i, n )

!*****************************************************************************80
!
!! I4_BIT_REVERSE reverses the bits in an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!       I      N  2^N     I4_BIT_REVERSE ( I, N )
!    ----    --------  -----------------------
!       0      0    1     0
!       1      0    1     1
!
!       0      3    8     0
!       1      3    8     4
!       2      3    8     2
!       3      3    8     6
!       4      3    8     1
!       5      3    8     5
!       6      3    8     3
!       7      3    8     7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the integer to be bit reversed.
!    I should be nonnegative.  Normally I < 2^N.
!
!    Input, integer ( kind = 4 ) N, indicates the number of bits to
!    be reverse (N+1) or the base with respect to which the integer is to
!    be reversed (2^N).  N should be nonnegative.
!
!    Output, integer ( kind = 4 ) I4_BIT_REVERSE, the bit reversed value.
!
  implicit none

  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_bit_reverse
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  if ( i < 0 ) then

    value = -1

  else if ( n < 0 ) then

    value = -1

  else

    b = 2**n
    j = mod ( i, b )

    value = 0

    do

      if ( b == 1 ) then

        value = value + j
        j = 0
        exit

      else

        if ( mod ( j, 2 ) == 1 ) then
          value = value + b / 2
          j = j - 1
        end if

        j = j / 2
        b = b / 2

      end if

    end do

  end if

  i4_bit_reverse = value

  return
end
function i4_ceiling ( r )

!*****************************************************************************80
!
!! I4_CEILING rounds an R8 "up" (towards +oo) to the next I4.
!
!  Example:
!
!    R     Value
!
!   -1.1  -1
!   -1.0  -1
!   -0.9   0
!    0.0   0
!    5.0   5
!    5.1   6
!    5.9   6
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the value to be rounded up.
!
!    Output, integer ( kind = 4 ) I4_CEILING, the rounded value.
!
  implicit none

  integer ( kind = 4 ) i4_ceiling
  real ( kind = 8 ) r
  integer ( kind = 4 ) value

  value = int ( r )
  if ( real ( value, kind = 8 ) < r ) then
    value = value + 1
  end if

  i4_ceiling = value

  return
end
function i4_characteristic ( q )

!*****************************************************************************80
!
!! I4_CHARACTERISTIC gives the characteristic for an I4.
!
!  Discussion:
!
!    For any positive integer Q, the characteristic is:
!
!    Q, if Q is a prime;
!    P, if Q = P^N for some prime P and some integer N;
!    0, otherwise, that is, if Q is negative, 0, 1, or the product
!       of more than one distinct prime.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738:
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, 1994, pages 494-495.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Q, the value to be tested.
!
!    Output, integer ( kind = 4 ) I4_CHARACTERISTIC, the characteristic of Q.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_characteristic
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) q
  integer ( kind = 4 ) q_copy

  if ( q <= 1 ) then
    i4_characteristic = 0
    return
  end if
!
!  If Q is not prime, then there is at least one prime factor
!  of Q no greater than SQRT(Q)+1.
!
!  A faster code would only consider prime values of I,
!  but that entails storing a table of primes and limiting the
!  size of Q.  Simplicity and flexibility for now!
!
  i_max = int ( sqrt ( real ( q ) ) ) + 1
  q_copy = q

  do i = 2, i_max

    if ( mod ( q_copy, i ) == 0 ) then

      do while ( mod ( q_copy, i ) == 0 )
        q_copy = q_copy / i
      end do

      if ( q_copy == 1 ) then
        i4_characteristic = i
      else
        i4_characteristic = 0
      end if

      return

    end if

  end do
!
!  If no factor was found, then Q is prime.
!
  i4_characteristic = q

  return
end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K) as an I4.
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
function i4_div_rounded ( a, b )

!*****************************************************************************80
!
!! I4_DIV_ROUNDED computes the rounded result of I4 division.
!
!  Discussion:
!
!    This routine computes C = A / B, where A, B and C are integers
!    and C is the closest integer value to the exact real result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the number to be divided.
!
!    Input, integer ( kind = 4 ) B, the divisor, or the number of parts.
!
!    Output, integer ( kind = 4 ) I4_DIV_ROUNDED, the rounded result
!    of the division.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) a_abs
  integer ( kind = 4 ) b
  integer ( kind = 4 ) b_abs
  integer ( kind = 4 ) i4_div_rounded
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) value

  if ( a == 0 .and. b == 0 ) then

    value = i4_huge
 
  else if ( a == 0 ) then

    value = 0

  else if ( b == 0 ) then

    if ( a < 0 ) then
      value = - i4_huge
    else
      value = + i4_huge
    end if

  else

    a_abs = abs ( a )
    b_abs = abs ( b )

    value = a_abs / b_abs
!
!  Round the value.
!
    if ( ( 2 * value + 1 ) * b_abs < 2 * a_abs ) then
      value = value + 1
    end if
!
!  Set the sign.
!
    if ( ( a < 0 .and. 0 < b ) .or. ( 0 < a .and. b < 0 ) ) then
      value = - value
    end if

  end if

  i4_div_rounded = value

  return
end
function i4_divp ( i, j )

!*****************************************************************************80
!
!! I4_DIVP returns the smallest multiple of J greater than or equal to I.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!    I  J  I4_DIVP(I,J)
!
!    0  4    0
!    1  4    1
!    2  4    1
!    3  4    1
!    4  4    1
!    5  4    2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be analyzed.
!
!    Input, integer ( kind = 4 ) J, the number, multiples of which will
!    be compared against I.  J may not be zero.
!
!    Output, integer ( kind = 4 ) I4_DIVP, the smallest multiple of J that
!    is greater than or equal to I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_divp
  integer ( kind = 4 ) j

  if ( j /= 0 ) then
    i4_divp = 1 + ( i - 1 ) / j
  else
    i4_divp = 0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_DIVP - Fatal error!'
    write ( *, '(a)' ) '  The input value of J was zero!'
    stop
  end if

  return
end
function i4_even ( i )

!*****************************************************************************80
!
!! I4_EVEN returns TRUE if an I4 is even.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the integer to be tested.
!
!    Output, logical I4_EVEN, is TRUE if I is even.
!
  implicit none

  integer ( kind = 4 ) i
  logical i4_even

  i4_even = ( mod ( i, 2 ) == 0 )

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
function i4_factorial2 ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
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
!    Output, integer ( kind = 4 ) I4_FACTORIAL2, the value of N!!.
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
function i4_floor ( r )

!*****************************************************************************80
!
!! I4_FLOOR rounds an R8 "down" (towards -oo) to the nearest I4.
!
!  Example:
!
!    R     Value
!
!   -1.1  -2
!   -1.0  -1
!   -0.9  -1
!    0.0   0
!    5.0   5
!    5.1   5
!    5.9   5
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the value to be rounded down.
!
!    Output, integer ( kind = 4 ) I4_FLOOR, the rounded value.
!
  implicit none

  integer ( kind = 4 ) i4_floor
  real ( kind = 8 ) r
  integer ( kind = 4 ) value

  value = int ( r )
  if ( r < real ( value, kind = 8 ) ) then
    value = value - 1
  end if

  i4_floor = value

  return
end
subroutine i4_fraction ( i, j, k )

!*****************************************************************************80
!
!! I4_FRACTION computes a ratio and returns an integer result.
!
!  Discussion:
!
!    Given integer variables I and J, FORTRAN will evaluate the expression
!    "I/J" using integer arithmetic.  This routine, which carries out the
!    same operation, is thus not needed in FORTRAN.  It is provided simply
!    to match the corresponding function in MATLAB, where the default
!    result of "I/J" is a real number.
!
!  Example:
!
!       I     J     Real     K = I4_FRACTION ( I, J)
!
!       1     2     0.5      0
!       8     4     2.00     2
!       9     4     2.25     2
!       7     4     1.75     1
!      -7     4    -1.75    -1
!       7    -4    -1.75    -1
!      -7    -4     1.75     1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the arguments.
!
!    Output, integer ( kind = 4 ) K, the value of the ratio.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i / j

  return
end
function i4_gcd ( i, j )

!*****************************************************************************80
!
!! I4_GCD finds the greatest common divisor of two I4's.
!
!  Discussion:
!
!    Note that only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I4_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I4_GCD is the
!    greatest common divisor of I and J.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, two numbers whose GCD is desired.
!
!    Output, integer ( kind = 4 ) I4_GCD, the greatest common divisor
!    of I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r

  i4_gcd = 1
!
!  Return immediately if either I or J is zero.
!
  if ( i == 0 ) then
    i4_gcd = max ( 1, abs ( j ) )
    return
  else if ( j == 0 ) then
    i4_gcd = max ( 1, abs ( i ) )
    return
  end if
!
!  Set P to the larger of I and J, Q to the smaller.
!  This way, we can alter P and Q as we go.
!
  p = max ( abs ( i ), abs ( j ) )
  q = min ( abs ( i ), abs ( j ) )
!
!  Carry out the Euclidean algorithm.
!
  do

    r = mod ( p, q )

    if ( r == 0 ) then
      exit
    end if

    p = q
    q = r

  end do

  i4_gcd = q

  return
end
function i4_gcdb ( i, j, k )

!*****************************************************************************80
!
!! I4_GCDB finds the greatest common divisor of the form K**N of two I4's.
!
!  Discussion:
!
!    Note that if J is negative, I4_GCDB will also be negative.
!    This is because it is likely that the caller is forming
!    the fraction I/J, and so any minus sign should be
!    factored out of J.
!
!    If I and J are both zero, I4_GCDB is returned as 1.
!
!    If I is zero and J is not, I4_GCDB is returned as J,
!    and vice versa.
!
!    If I and J are nonzero, and have no common divisor of the
!    form K**N, I4_GCDB is returned as 1.
!
!    Otherwise, I4_GCDB is returned as the largest common divisor
!    of the form K**N shared by I and J.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, two numbers whose greatest common
!    divisor K**N is desired.
!
!    Input, integer ( kind = 4 ) K, the possible divisor of I and J.
!
!    Output, integer ( kind = 4 ) I4_GCDB, the greatest common divisor of
!    the form K**N shared by I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icopy
  integer ( kind = 4 ) i4_gcdb
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcopy
  integer ( kind = 4 ) k

  i4_gcdb = 1
!
!  If both I and J are zero, I4_GCDB is 1.
!
  if ( i == 0 .and. j == 0 ) then
    i4_gcdb = 1
    return
  end if
!
!  If just one of I and J is zero, I4_GCDB is the other one.
!
  if ( i == 0 ) then
    i4_gcdb = j
    return
  else if ( j == 0 ) then
    i4_gcdb = i
    return
  end if
!
!  Divide out K as long as you can.
!
  if ( 0 < j ) then
    i4_gcdb = 1
  else
    i4_gcdb = -1
  end if

  icopy = i
  jcopy = j

  do

    if ( mod ( icopy, k ) /= 0 .or. mod ( jcopy, k ) /= 0 ) then
      exit
    end if

    i4_gcdb = i4_gcdb * k
    icopy = icopy / k
    jcopy = jcopy / k

  end do

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
function i4_huge_normalizer ( )

!*****************************************************************************80
!
!! I4_HUGE_NORMALIZER returns the "normalizer" for I4_HUGE.
!
!  Discussion:
!
!    The value returned is 1 / ( I4_HUGE + 1 ).
!
!    For any I4, it should be the case that
!
!     -1 < I4 * I4_HUGE_NORMALIZER < 1.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) I4_HUGE_NORMALIZER, the "normalizer"
!    for I4_HUGE.
!
  implicit none

  real ( kind = 8 ) i4_huge_normalizer

  i4_huge_normalizer = 4.656612873077392578125D-10

  return
end
function i4_is_power_of_2 ( n )

!*****************************************************************************80
!
!! I4_IS_POWER_OF_2 reports whether an I4 is a power of 2.
!
!  Discussion:
!
!    The powers of 2 are 1, 2, 4, 8, 16, and so on.
!
!    An I4 is an integer ( kind = 4 ) value.
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
!    Input, integer ( kind = 4 ) N, the integer to be tested.
!
!    Output, logical I4_IS_POWER_OF_2, is TRUE if N is a power of 2.
!
  implicit none

  logical i4_is_power_of_2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_copy

  n_copy = n
  i4_is_power_of_2 = .false.

  if ( n_copy <= 0 ) then
    return
  end if

  do while ( n_copy /= 1 )

    if ( mod ( n_copy, 2 ) == 1 ) then
      return
    end if

    n_copy = n_copy / 2

  end do

  i4_is_power_of_2 = .true.

  return
end
function i4_is_prime ( n )

!*****************************************************************************80
!
!! I4_IS_PRIME reports whether an I4 is prime.
!
!  Discussion:
!
!    A simple, unoptimized sieve of Erasthosthenes is used to
!    check whether N can be divided by any integer between 2
!    and SQRT(N).
!
!    Note that negative numbers, 0 and 1 are not considered prime.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer to be tested.
!
!    Output, logical I4_IS_PRIME, is TRUE if N is prime, and FALSE
!    otherwise.
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

  if ( n == 1 ) then
    i4_is_prime = .false.
    return
  end if

  if ( n <= 3 ) then
    i4_is_prime = .true.
    return
  end if

  nhi = int ( sqrt ( real ( n ) ) )

  do i = 2, nhi
    if ( mod ( n, i ) == 0 ) then
      i4_is_prime = .false.
      return
    end if
  end do

  i4_is_prime = .true.

  return
end
function i4_lcm ( i, j )

!*****************************************************************************80
!
!! I4_LCM computes the least common multiple of two I4's.
!
!  Discussion:
!
!    The least common multiple may be defined as
!
!      LCM(I,J) = ABS ( I * J ) / GCD(I,J)
!
!    where GCD(I,J) is the greatest common divisor of I and J.
!
!    An I4 is an integer ( kind = 4 ) value.
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
!    Input, integer ( kind = 4 ) I, J, the integers whose I4_LCM is desired.
!
!    Output, integer ( kind = 4 ) I4_LCM, the least common multiple of I and J.
!    I4_LCM is never negative.  I4_LCM is 0 if either I or J is zero.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) j
  integer ( kind = 4 ) i4_lcm

  i4_lcm = abs ( i * ( j / i4_gcd ( i, j ) ) )

  return
end
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 10
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
!
!  Discussion:
!
!    For positive I4_LOG_2(I), it should be true that
!      2^I4_LOG_2(X) <= |I| < 2^(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 2
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_2, the integer part of the
!    logarithm base 2 of the absolute value of I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647

  if ( i == 0 ) then

    i4_log_2 = - i4_huge

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end
function i4_log_i4 ( i4, j4 )

!*****************************************************************************80
!
!! I4_LOG_I4 returns the logarithm of an I4 to an I4 base.
!
!  Discussion:
!
!    Only the integer part of the logarithm is returned.
!
!    If
!
!      K4 = I4_LOG_J4 ( I4, J4 ),
!
!    then we ordinarily have
!
!      J4^(K4-1) < I4 <= J4^K4.
!
!    The base J4 should be positive, and at least 2.  If J4 is negative,
!    a computation is made using the absolute value of J4.  If J4 is
!    -1, 0, or 1, the logarithm is returned as 0.
!
!    The number I4 should be positive and at least 2.  If I4 is negative,
!    a computation is made using the absolute value of I4.  If I4 is
!    -1, 0, or 1, then the logarithm is returned as 0.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!    I4  J4  K4
!
!     0   3   0
!     1   3   0
!     2   3   0
!     3   3   1
!     4   3   1
!     8   3   1
!     9   3   2
!    10   3   2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, the number whose logarithm is desired.
!
!    Input, integer ( kind = 4 ) J4, the base of the logarithms.
!
!    Output, integer ( kind = 4 ) I4_LOG_I4, the integer part of the logarithm
!    base abs(J4) of abs(I4).
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_abs
  integer ( kind = 4 ) i4_log_i4
  integer ( kind = 4 ) j4
  integer ( kind = 4 ) j4_abs
  integer ( kind = 4 ) value

  value = 0

  i4_abs = abs ( i4 )

  if ( 2 <= i4_abs ) then

    j4_abs = abs ( j4 )

    if ( 2 <= j4_abs ) then

      do while ( j4_abs <= i4_abs )
        i4_abs = i4_abs / j4_abs
        value = value + 1
      end do

    end if

  end if

  i4_log_i4 = value

  return
end
function i4_log_r8 ( x, b )

!*****************************************************************************80
!
!! I4_LOG_R8 returns the logarithm of an I4 to an R8 base.
!
!  Discussion:
!
!    The base B should be positive, but in any case only the absolute
!    value of B is considered.
!
!    The number X whose logarithm is desired should be positive, but
!    in any case only the absolute value of X is considered.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    An R8 is a real ( kind = 8 ) value.
!
!  Example:
!
!    If B is greater than 1, and X is positive:
!
!    if 1/B^2  <  X <= 1/B   I4_LOG_R8(X) = -1,
!    if 1/B    <  X <= 1     I4_LOG_R8(X) = 0,
!    if 1      <= X <  B,    I4_LOG_R8(X) = 0,
!    if B      <= X <  B^2  I4_LOG_R8(X) = 1,
!    if B^2    <= X <  B^3  I4_LOG_R8(X) = 2.
!
!    For positive I4_LOG_R8(X), it should be true that
!
!      ABS(B)^I4_LOG_R8(X) <= ABS(X) < ABS(B)^(I4_LOG_R8(X)+1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number whose logarithm base B is
!    desired.  If X is 0, then I4_LOG_B is returned as -I4_HUGE().
!
!    Input, real ( kind = 8 ) B, the absolute value of the base of the
!    logarithms.  B must not be -1, 0, or 1.
!
!    Output, integer ( kind = 4 ) I4_LOG_R8, the integer part of the logarithm
!    base abs(B) of abs(X).
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_abs
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_log_r8
  integer ( kind = 4 ) value_sign
  integer ( kind = 4 ) x
  real ( kind = 8 ) x_abs

  if ( x == 0 ) then
    i4_log_r8 = - i4_huge
    return
  end if

  b_abs = abs ( b )
  i4_log_r8 = 0

  if ( b_abs == 1.0D+00 ) then
    return
  end if

  if ( b == 0.0D+00 ) then
    return
  end if

  x_abs = abs ( real ( x ) )

  if ( b_abs < 1.0D+00 ) then
    value_sign = -1
    b_abs = 1.0D+00 / b_abs
  else
    value_sign = +1
  end if

  if ( 1.0D+00 <= x_abs .and. x_abs < b_abs ) then
    i4_log_r8 = value_sign * i4_log_r8
    return
  end if

  do while ( b_abs < x_abs )
    x_abs = x_abs / b_abs
    i4_log_r8 = i4_log_r8 + 1
  end do

  do while ( x_abs * b_abs <= 1.0D+00 )
    x_abs = x_abs * b_abs
    i4_log_r8 = i4_log_r8 - 1
  end do
!
!  If the absolute value of the base was less than 1, we inverted
!  earlier.  Now negate the logarithm to account for that.
!
  i4_log_r8 = value_sign * i4_log_r8

  return
end
subroutine i4_mant ( x, s, j, k, l )

!*****************************************************************************80
!
!! I4_MANT computes the "mantissa" of a double precision number.
!
!  Discussion:
!
!    I4_MANT computes the "mantissa" or "fraction part" of a real
!    number X, which it stores as a pair of integers, (J/K).
!
!    It also computes the sign, and the integer part of the logarithm
!    (base 2) of X.
!
!    On return:
!
!      X = S * (J/K) * 2^L
!
!    where
!
!      S is +1 or -1,
!      K is a power of 2,
!      1 <= (J/K) < 2,
!      L is an integer.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, integer ( kind = 4 ) S, the "sign" of the number.
!    S will be -1 if X is less than 0, and +1 if X is greater
!    than or equal to zero.
!
!    Output, integer ( kind = 4 ) J, the top part of the mantissa fraction.
!
!    Output, integer ( kind = 4 ) K, the bottom part of the mantissa
!    fraction.  K is a power of 2.
!
!    Output, integer ( kind = 4 ) L, the integer part of the logarithm
!    (base 2) of X.
!
  implicit none

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) xtemp
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0D+00 ) then
    s = 1
    j = 0
    k = 1
    l = 0
    return
  end if
!
!  2: Determine the sign S.
!
  if ( 0.0D+00 < x ) then
    s = + 1
    xtemp = + x
  else
    s = - 1
    xtemp = - x
  end if
!
!  3: Force XTEMP to lie between 1 and 2, and compute the logarithm L.
!
  l = 0

  do while ( 2.0D+00 <= xtemp )
    xtemp = xtemp / 2.0D+00
    l = l + 1
  end do

  do while ( xtemp < 1.0D+00 )
    xtemp = xtemp * 2.0D+00
    l = l - 1
  end do
!
!  4: Now strip out the mantissa as J/K.
!
  j = 0
  k = 1

  do

    j = 2 * j

    if ( 1.0D+00 <= xtemp ) then
      j = j + 1
      xtemp = xtemp - 1.0D+00
    end if

    if ( xtemp == 0.0D+00 ) then
      exit
    end if

    k = 2 * k
    xtemp = xtemp * 2.0D+00

  end do

  return
end
subroutine i4_mod_inv ( b, n, y )

!*****************************************************************************80
!
!! I4_MOD_INV calculates the inverse of B mod N.
!
!  Discussion:
!
!    This function uses the extended Euclidean algorithm.
!
!    Unless the algorithm fails, the output value Y will satisfy
!
!      ( B * Y ) mod N = 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2011
!
!  Author:
!
!    Original MATLAB version by Wade Trappe, Lawrence Washington.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wade Trappe, Lawrence Washington,
!    Introduction to Cryptography with Coding Theory,
!    Prentice Hall, 2005,
!    ISBN13: 978-0131862395,
!    LC: QA268.T73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) B, the value whose inverse is desired.
!    B must not be 0, or a multiple of N.  However, B can be negative.
!
!    Input, integer ( kind = 4 ) N, the value with respect to which the inverse
!    is desired.  N must be 2 or greater.
!
!    Output, integer ( kind = 4 ) Y, the inverse of B mod N.  However, if the
!    inverse does not exist, Y is returned as 0.
!
  implicit none

  integer ( kind = 4 ) b
  integer ( kind = 4 ) b0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) t
  integer ( kind = 4 ) t0
  integer ( kind = 4 ) temp
  integer ( kind = 4 ) y

  n0 = n
  b0 = abs ( b )
  t0 = 0
  t = 1

  q = ( n0 / b0 )
  r = n0 - q * b0

  do while ( 0 < r )

    temp = t0 - q * t

    if ( 0 <= temp ) then
      temp =     mod (   temp, n )
    else
      temp = n - mod ( - temp, n )
    end if

    n0 = b0
    b0 = r
    t0 = t
    t = temp

    q = ( n0 / b0 )
    r = n0 - q * b0

  end do

  if ( b0 /= 1 ) then
    y = 0
  else
    y = mod ( t, n )
    if ( b < 0 ) then
      y = - y
    end if
  end if

  return
end
subroutine i4_moddiv ( n, d, m, r )

!*****************************************************************************80
!
!! I4_MODDIV breaks an I4 into a multiple of a divisor and remainder.
!
!  Discussion:
!
!    The formula used is:
!
!      N = M * D + R
!
!      0 <= || R || < || D ||
!
!    and R has the sign of N.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!      N         D       M      R
!
!     107       50      2      7
!     107      -50     -2      7
!    -107       50     -2     -7
!    -107      -50      2     -7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number to be decomposed.
!
!    Input, integer ( kind = 4 ) D, the divisor.  D may not be zero.
!
!    Output, integer ( kind = 4 ) M, the number of times N
!    is evenly divided by D.
!
!    Output, integer ( kind = 4 ) R, a remainder.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) r

  if ( d == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODDIV - Fatal error!'
    write ( *, '(a)' ) '  Input divisor D = 0'
    stop
  end if

  m = n / d
  r = n - d * m

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
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
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_mop ( i )

!*****************************************************************************80
!
!! I4_MOP returns the I-th power of -1 as an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, integer ( kind = 4 ) I4_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_mop

  if ( mod ( i, 2 ) == 0 ) then
    i4_mop = 1
  else
    i4_mop = -1
  end if

  return
end
function i4_odd ( i )

!*****************************************************************************80
!
!! I4_ODD returns TRUE if an I4 is odd.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the integer to be tested.
!
!    Output, logical I4_ODD, is TRUE if I is odd.
!
  implicit none

  integer ( kind = 4 ) i
  logical i4_odd

  i4_odd = ( mod ( i + 1, 2 ) == 0 )

  return
end
function i4_power ( i, j )

!*****************************************************************************80
!
!! I4_POWER returns the integer power of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
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
!    Input, integer ( kind = 4 ) I, J, the base and the power.
!    J should be nonnegative.
!
!    Output, integer ( kind = 4 ) I4_POWER, the value of I^J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_power
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  if ( j < 0 ) then

    if ( i == 1 ) then
      i4_power = 1
    else if ( i == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4_POWER - Fatal error!'
      write ( *, '(a)' ) '  I^J requested, with I = 0 and J negative.'
      stop
    else
      i4_power = 0
    end if

  else if ( j == 0 ) then

    if ( i == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4_POWER - Fatal error!'
      write ( *, '(a)' ) '  I^J requested, with I = 0 and J = 0.'
      stop
    else
      i4_power = 1
    end if

  else if ( j == 1 ) then

    i4_power = i

  else

    i4_power = 1
    do k = 1, j
      i4_power = i4_power * i
    end do

  end if

  return
end
function i4_sign ( x )

!*****************************************************************************80
!
!! I4_SIGN evaluates the sign of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number whose sign is desired.
!
!    Output, integer ( kind = 4 ) I4_SIGN, the sign of X:
!
  implicit none

  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    i4_sign = -1
  else
    i4_sign = +1
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
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
subroutine i4_swap3 ( i, j, k )

!*****************************************************************************80
!
!! I4_SWAP3 swaps three I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On output, the
!    values of I, J, and K have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  l = i
  i = j
  j = k
  k = l

  return
end
subroutine i4_to_angle ( i, angle )

!*****************************************************************************80
!
!! I4_TO_ANGLE maps I4's to points on a circle.
!
!  Discussion:
!
!    The angles are intended to be used to select colors on a color
!    hexagon whose 6 vertices are red, yellow, green, cyan, blue,
!    magenta.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!     I   X      ANGLE
!
!     0   0/3      0
!     1   1/3    120
!     2   2/3    240
!
!     3   1/6     60
!     4   3/6    180
!     5   5/6    300
!
!     6   1/12    30
!     7   3/12    90
!     8   5/12   150
!     9   7/12   210
!    10   9/12   270
!    11  11/12   330
!
!    12   1/24    15
!    13   3/24    45
!    14   5/24    75
!    etc
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired color.
!
!    Output, real ( kind = 8 ) ANGLE, an angle, measured in degrees,
!    between 0 and 360.
!
  implicit none

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4

  if ( 0 <= abs ( i ) .and. abs ( i ) <= 2 ) then

    angle = 120.0D+00 * real ( abs ( i ), kind = 8 )

  else

    i1 = i4_log_2 ( abs ( i ) / 3 )
    i2 = abs ( i ) + 1 - 3 * 2**i1
    i3 = 2 * ( i2 - 1 ) + 1
    i4 = 3 * 2**( i1 + 1 )

    angle = 360.0D+00 * real ( i3, kind = 8 ) / real ( i4, kind = 8 )

  end if

  return
end
subroutine i4_to_digits_binary ( i, n, c )

!*****************************************************************************80
!
!! I4_TO_DIGITS_BINARY produces the binary digits of an I4.
!
!  Discussion:
!
!    An I4 is an integer.
!
!  Example:
!
!     I    N     C               Binary
!    --  ---   ---         ------------
!     0    1   0                      0
!     0    2   0, 0                  00
!     1    3   1, 0, 0              100
!     2    3   0, 1, 0              010
!     3    3   1, 1, 0              011
!     4    3   0, 0, 1              100
!     8    3   0, 0, 0           (1)000
!     8    5   0, 0, 0, 1, 0      01000
!    -8    5   0, 0, 0, 1, 0  (-) 01000
!
!     0    3   0, 0, 0
!     1    3   1, 0, 0
!     2    3   0, 1, 0
!     3    3   1, 1, 0
!     4    3   0, 0, 1
!     5    3   1, 0, 1
!     6    3   0, 1, 1
!     7    3   1, 1, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, an integer to be represented.
!
!    Input, integer ( kind = 4 ) N, the number of binary digits to produce.
!
!    Output, integer ( kind = 4 ) C(N), the first N binary digits of I,
!    with C(1) being the units digit.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_copy
  integer ( kind = 4 ) j

  i_copy = abs ( i )

  do j = 1, n

    c(j) = mod ( i_copy, 2 )
    i_copy = i_copy / 2

  end do

  return
end
subroutine i4_to_digits_decimal ( i, n, digit )

!*****************************************************************************80
!
!! I4_TO_DIGITS_DECIMAL determines the last N decimal digits of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the integer to be analyzed.
!
!    Input, integer ( kind = 4 ) N, the number of digits to determine.
!
!    Output, integer ( kind = 4 ) DIGIT(N), the last N decimal digits of I.
!    DIGIT(I) is the "coefficient" of 10**(I-1).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) digit(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_copy
  integer ( kind = 4 ) j

  i_copy = abs ( i )

  do j = 1, n
    digit(j) = mod ( i_copy, 10 )
    i_copy = ( i_copy - digit(j) ) / 10
  end do

  return
end
subroutine i4_to_fac ( intval, prime_num, npower )

!*****************************************************************************80
!
!! I4_TO_FAC converts an I4 into a product of prime factors.
!
!  Discussion:
!
!    This routine will fail if the input integer is not positive,
!    or if PRIME_NUM is too small to account for the factors of the integer.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The formula is:
!
!      INTVAL = Product ( 1 <= I <= PRIME_NUM ) PRIME(I)**NPOWER(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INTVAL, the integer to be factored.
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the number of prime factors for
!    which storage has been allocated.
!
!    Output, integer ( kind = 4 ) NPOWER(PRIME_NUM), the powers of the primes.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) intcopy
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) npower(prime_num)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) prime

  if ( intval <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_FAC - Fatal error!'
    write ( *, '(a)' ) '  Input integer is not positive.'
    stop
  end if
!
!  Try dividing the remainder by each prime.
!
  intcopy = intval

  do i = 1, prime_num

    npower(i) = 0

    p = prime ( i )

    do while ( mod ( intcopy, p ) == 0 )
      npower(i) = npower(i) + 1
      intcopy = intcopy / p
    end do

  end do

  return
end
subroutine i4_to_halton ( dim_num, step, seed, leap, base, r )

!*****************************************************************************80
!
!! I4_TO_HALTON computes one element of a leaped Halton subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Halton sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Halton sequence.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 84-90
!
!    John Halton, GB Smith,
!    Algorithm 247:
!    Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, Number 12, December 1964, pages 701-702.
!
!    Ladislav Kocis, William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, June 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the Halton sequence index
!    corresponding to STEP = 0.
!    0 <= SEED(1:DIM_NUM) is required.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in the
!    Halton sequence.  1 <= LEAP(1:DIM_NUM) is required.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the Halton bases.
!    1 < BASE(1:DIM_NUM) is required.
!
!    Output, real ( kind = 8 ) R(DIM_NUM), the STEP-th element of the leaped
!    Halton subsequence.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  real ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) leap(dim_num)
  real ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) seed2
  integer ( kind = 4 ) step
!
!  Check the input.
!
  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( step < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) ' STEP < 0.'
    stop
  end if

  if ( any ( seed(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  Some SEED(*) < 0.'
    stop
  end if

  if ( any ( leap(1:dim_num) < 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  Some LEAP < 1.'
    stop
  end if

  if ( any ( base(1:dim_num) <= 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  Some BASE <= 1.'
    stop
  end if
!
!  Calculate the data.
!
  do i = 1, dim_num

    seed2 = seed(i) + step * leap(i)

    r(i) = 0.0D+00

    base_inv = real ( 1.0D+00, kind = 8 ) / real ( base(i), kind = 8 )

    do while ( seed2 /= 0 )
      digit = mod ( seed2, base(i) )
      r(i) = r(i) + real ( digit, kind = 8 ) * base_inv
      base_inv = base_inv / real ( base(i), kind = 8 )
      seed2 = seed2 / base(i)
    end do

  end do

  return
end
function i4_to_isbn ( i )

!*****************************************************************************80
!
!! I4_TO_ISBN converts an I4 to an ISBN digit.
!
!  Discussion:
!
!    Only the integers 0 through 10 can be input.  The representation
!    of 10 is 'X'.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Book Industry Study Group,
!    The Evolution in Product Identification:
!    Sunrise 2005 and the ISBN-13,
!    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, an integer between 0 and 10.
!
!    Output, character I4_TO_ISBN, the ISBN character code of the integer.
!    If I is illegal, then I4_TO_ISBN is set to '?'.
!
  implicit none

  integer ( kind = 4 ) i
  character i4_to_isbn

       if ( i == 0 ) then
    i4_to_isbn = '0'
  else if ( i == 1 ) then
    i4_to_isbn = '1'
  else if ( i == 2 ) then
    i4_to_isbn = '2'
  else if ( i == 3 ) then
    i4_to_isbn = '3'
  else if ( i == 4 ) then
    i4_to_isbn = '4'
  else if ( i == 5 ) then
    i4_to_isbn = '5'
  else if ( i == 6 ) then
    i4_to_isbn = '6'
  else if ( i == 7 ) then
    i4_to_isbn = '7'
  else if ( i == 8 ) then
    i4_to_isbn = '8'
  else if ( i == 9 ) then
    i4_to_isbn = '9'
  else if ( i == 10 ) then
    i4_to_isbn = 'X'
  else
    i4_to_isbn = '?'
  end if

  return
end
function i4_to_l ( i4 )

!*****************************************************************************80
!
!! I4_TO_L converts an I4 to a logical value.
!
!  Discussion:
!
!    0 is FALSE, and anything else if TRUE.
!
!    An I4 is an integer value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, an integer.
!
!    Output, logical I4_TO_L, the logical value of I4.
!
  implicit none

  integer ( kind = 4 ) i4
  logical i4_to_l
  logical value

  value = ( i4 /= 0 )

  i4_to_l = value

  return
end
function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
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
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
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

  i4_uniform_ab = value

  return
end
subroutine i4_unswap3 ( i, j, k )

!*****************************************************************************80
!
!! I4_UNSWAP3 unswaps three I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On output, the values
!    of I, J, and K have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  l = k
  k = j
  j = i
  i = l

  return
end
function i4_walsh_1d ( x, digit )

!*****************************************************************************80
!
!! I4_WALSH_1D evaluates the Walsh function.
!
!  Discussion:
!
!    Consider the binary representation of X, and number the digits
!    in descending order, from leading to lowest, with the units digit
!    being numbered 0.
!
!    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
!
!    An I4 is an integer ( kind = 4 ) value.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Walsh function.
!
!    Input, integer ( kind = 4 ) DIGIT, the index of the Walsh function.
!
!    Output, integer ( kind = 4 ) I4_WALSH_1D, the value of the Walsh function.
!
  implicit none

  integer ( kind = 4 ) digit
  integer ( kind = 4 ) i4_walsh_1d
  integer ( kind = 4 ) n
  real ( kind = 8 ) x
  real ( kind = 8 ) x_copy
!
!  Hide the effect of the sign of X.
!
  x_copy = abs ( x )
!
!  If DIGIT is positive, divide by 2 DIGIT times.
!  If DIGIT is negative, multiply by 2 (-DIGIT) times.
!
  x_copy = x_copy / 2.0D+00**digit
!
!  Make it an integer.
!  Because it's positive, and we're using INT, we don't change the
!  units digit.
!
  n = int ( x_copy )
!
!  Is the units digit odd or even?
!
  if ( mod ( n, 2 ) == 0 ) then
    i4_walsh_1d = 0
  else
    i4_walsh_1d = 1
  end if

  return
end
function i4_width ( i )

!*****************************************************************************80
!
!! I4_WIDTH returns the "width" of an I4.
!
!  Discussion:
!
!    The width of an integer is the number of characters necessary to print it.
!
!    The width of an integer can be useful when setting the appropriate output
!    format for a vector or array of values.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I  I4_WIDTH
!    -----  -------
!    -1234    5
!     -123    4
!      -12    3
!       -1    2
!        0    1
!        1    1
!       12    2
!      123    3
!     1234    4
!    12345    5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose width is desired.
!
!    Output, integer ( kind = 4 ) I4_WIDTH, the number of characters
!    necessary to represent the integer in base 10, including a negative
!    sign if necessary.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) i4_width

  if ( 0 < i ) then
    i4_width = i4_log_10 ( i ) + 1
  else if ( i == 0 ) then
    i4_width = 1
  else if ( i < 0 ) then
    i4_width = i4_log_10 ( i ) + 2
  end if

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
function i4_xor ( i, j )

!*****************************************************************************80
!
!! I4_XOR calculates the exclusive OR of two I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    FORTRAN offers the function IEOR ( I, J ) which should be used instead.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, two values whose exclusive OR is needed.
!
!    Output, integer ( kind = 4 ) I4_XOR, the exclusive OR of I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_xor
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  i1 = i
  j1 = j
  k = 0
  l = 1

  do while ( i1 /= 0 .or. j1 /= 0 )

    i2 = i1 / 2
    j2 = j1 / 2

    if ( &
      ( ( i1 == 2 * i2 ) .and. ( j1 /= 2 * j2 ) ) .or. &
      ( ( i1 /= 2 * i2 ) .and. ( j1 == 2 * j2 ) ) ) then
      k = k + l
    end if

    i1 = i2
    j1 = j2
    l = 2 * l

  end do

  i4_xor = k

  return
end
subroutine i43mat_flip_cols ( m, n, a )

!*****************************************************************************80
!
!! I43MAT_FLIP_COLS swaps the columns of an I43MAT.
!
!  Discussion:
!
!    An I43MAT is a matrix, each of whose entries is an I43,
!    a triple of I4's.
!
!    An I43MAT can be stored as a 3 x M x N array, where M counts the "columns"
!    and N counts the "rows".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(3,M,N), the matrix whose columns
!    are to be flipped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(3,m,n)
  integer ( kind = 4 ) b(3,m,1)
  integer ( kind = 4 ) j

  do j = 1, n / 2
    b(1:3,1:m,    1) = a(1:3,1:m,    j)
    a(1:3,1:m,    j) = a(1:3,1:m,n+1-j)
    a(1:3,1:m,n+1-j) = b(1:3,1:m,    1)
  end do

  return
end
subroutine i43mat_flip_rows ( m, n, a )

!*****************************************************************************80
!
!! I43MAT_FLIP_ROWS swaps the rows of an I43MAT.
!
!  Discussion:
!
!    An I43MAT is a matrix, each of whose entries is an I43,
!    a triple of I4's.
!
!    An I43MAT can be stored as a 3 x M x N array, where M counts the "columns"
!    and N counts the "rows".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(3,M,N), the matrix whose rows
!    are to be flipped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(3,m,n)
  integer ( kind = 4 ) b(3,1,n)
  integer ( kind = 4 ) i

  do i = 1, m / 2
    b(1:3,    1,1:n) = a(1:3,    i,1:n)
    a(1:3,    i,1:n) = a(1:3,m+1-i,1:n)
    a(1:3,m+1-i,1:n) = b(1:3,    1,1:n)
  end do

  return
end
subroutine i4block_print ( l, m, n, a, title )

!*****************************************************************************80
!
!! I4BLOCK_PRINT prints an I4BLOCK.
!
!  Discussion:
!
!    An I4BLOCK is a 3D array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, N, the dimensions of the block.
!
!    Input, integer ( kind = 4 ) A(L,M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(l,m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do k = 1, n

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  K = ', k

    do jlo = 1, m, 10
      jhi = min ( jlo + 10 - 1, m )
      write ( *, '(a)' ) ' '
      write ( *, '(8x,a2,10(2x,i6))' ) 'J:', ( j, j = jlo, jhi )
      write ( *, '(7x,a2)' ) 'I:'
      do i = 1, l
        write ( *, '(2x,i6,a1,1x,10(2x,i6))' ) i, ':', a(i,jlo:jhi,k)
      end do
    end do

  end do

  return
end
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of
!    vectors of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  Column index I = ', i, ' is less than 1.'
    stop
  end if

  if ( n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index I = ', i
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  Column index J = ', j, ' is less than 1.'
    stop
  end if

  if ( n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index J = ', j
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_find ( m, n, a, ivec, col )

!*****************************************************************************80
!
!! I4COL_FIND searches an I4COL for a particular column value.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    M = 3, N = 4,
!
!    A = (
!      1  2  3  4
!      5  6  7  8
!      9 10 11 12 )
!
!    IVEC = ( 3, 7, 11 )
!
!    COL = 3
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the table.  M is also the length of IVEC.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) IVEC(M), a vector to be matched with the data
!    in the array.
!
!    Output, integer ( kind = 4 ) COL, the index of the first column of
!    the table which exactly matches every entry of IVEC, or -1 if no match
!    could be found.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) col
  integer ( kind = 4 ) ivec(m)
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    col = -1
    return
  end if

  do j = 1, n

    i = 1

    do while ( ivec(i) == a(i,j) )

      if ( i == m ) then
        col = j
        return
      end if

      i = i + 1

    end do

  end do

  col = -1

  return
end
subroutine i4col_find_item ( m, n, a, item, row, col )

!*****************************************************************************80
!
!! I4COL_FIND_ITEM searches an I4COL for a given scalar value.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the table.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) ITEM, the value to search for.
!
!    Output, integer ( kind = 4 ) ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM.  The search
!    is conducted by columns.  If the item is not found, then
!    ROW = COL = -1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  do j = 1, n
    do i = 1, m
      if ( a(i,j) == item ) then
        row = i
        col = j
        return
      end if
    end do
  end do

  row = -1
  col = -1

  return
end
subroutine i4col_find_pair_wrap ( m, n, a, item1, item2, row, col )

!*****************************************************************************80
!
!! I4COL_FIND_PAIR_WRAP searches an I4COL for a pair of items.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!    The items (ITEM1, ITEM2) must occur consecutively.
!    However, wrapping is allowed, that is, if ITEM1 occurs
!    in the last row, and ITEM2 "follows" it in the first row
!    of the same column, a match is declared.
!
!    If the pair of items is not found, then ROW = COL = -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input, integer ( kind = 4 ) A(M,N), the array to search.
!
!    Input, integer ( kind = 4 ) ITEM1, ITEM2, the values to search for.
!
!    Output, integer ( kind = 4 ) ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM1 followed immediately
!    by ITEM2.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) item1
  integer ( kind = 4 ) item2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  do j = 1, n
    do i = 1, m

      if ( a(i,j) == item1 ) then

        i2 = i + 1

        if ( m < i2 ) then
          i2 = 1
        end if

        if ( a(i2,j) == item2 ) then
          row = i
          col = j
          return
        end if

      end if

    end do
  end do

  row = -1
  col = -1

  return
end
subroutine i4col_first_index ( m, n, a, first_index )

!*****************************************************************************80
!
!! I4COL_FIRST_INDEX indexes the first occurrence of values in an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's.
!    It is regarded as an array of N columns of length M.
!
!    For element A(1:M,J) of the matrix, FIRST_INDEX(J) is the index in A of
!    the first column whose entries are equal to A(1:M,J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!    The length of an "element" of A, and the number of "elements".
!
!    Input, integer ( kind = 4 ) A(M,N), the array.
!
!    Output, integer ( kind = 4 ) FIRST_INDEX(N), the first occurrence index.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) first_index(n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2

  first_index(1:n) = -1

  do j1 = 1, n

    if ( first_index(j1) == -1 ) then

      first_index(j1) = j1

      do j2 = j1 + 1, n
        if ( all ( a(1:m,j1) == a(1:m,j2) ) ) then
          first_index(j2) = j1
        end if
      end do

    end if

  end do

  return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_sort_d ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_D descending sorts an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in descending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )
      isgn = -isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_sort2_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M vectors.
!    On output, the elements of each column of A have been sorted in ascending
!    order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  do col = 1, n

    i = 0
    indx = 0
    isgn = 0
    j = 0
!
!  Call the external heap sorter.
!
    do

      call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
      if ( 0 < indx ) then

        t        = a(i,col)
        a(i,col) = a(j,col)
        a(j,col) = t
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(j,col) < a(i,col) ) then
          isgn = +1
        else
          isgn = -1
        end if

      else if ( indx == 0 ) then

        exit

      end if

    end do

  end do

  return
end
subroutine i4col_sort2_d ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT2_D descending sorts elements of each column of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M vectors.
!    On output, the elements of each column of A have been sorted in descending
!    order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  do col = 1, n

    i = 0
    indx = 0
    isgn = 0
    j = 0
!
!  Call the external heap sorter.
!
    do

      call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
      if ( 0 < indx ) then

        t        = a(i,col)
        a(i,col) = a(j,col)
        a(j,col) = t
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(i,col) < a(j,col) ) then
          isgn = +1
        else
          isgn = -1
        end if

      else if ( indx == 0 ) then

        exit

      end if

    end do

  end do

  return
end
subroutine i4col_sorted_singleton_count ( m, n, a, singleton_num )

!*****************************************************************************80
!
!! I4COL_SORTED_SINGLETON_COUNT counts singletons in an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!    A "singleton" is an item that occurs exactly once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) SINGLETON_NUM, the number of singletons.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  logical differ_from_next
  logical differ_from_previous
  integer ( kind = 4 ) j
  integer ( kind = 4 ) singleton_num

  singleton_num = 0

  if ( n <= 0 ) then
    return
  end if

  differ_from_next = .true.

  do j = 1, n

    differ_from_previous = differ_from_next

    if ( j < n ) then
      differ_from_next = any ( a(1:m,j) /= a(1:m,j+1) )
    else
      differ_from_next = .true.
    end if

    if ( differ_from_previous .and. differ_from_next ) then
      singleton_num = singleton_num + 1
    end if

  end do

  return
end
subroutine i4col_sorted_unique ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE keeps unique elements in a sorted I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!    The array can be sorted into ascending or descending order.
!    The important point is that identical elements must be stored
!    in adjacent positions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the sorted array of N columns of M-vectors.
!    On output, a sorted array of columns of M-vectors.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      j1 = j1 + 1
      a(1:m,j1) = a(1:m,j2)
    end if

  end do

  unique_num = j1

  return
end
subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine i4col_swap ( m, n, a, j1, j2 )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns J1 and J2 of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, J1 = 2, J2 = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns
!    of length M.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2

  if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
    write ( *, '(a,i8)' ) '  J1 =    ', j1
    write ( *, '(a,i8)' ) '  J2 =    ', j2
    write ( *, '(a,i8)' ) '  N =     ', n
    stop

  end if

  if ( j1 == j2 ) then
    return
  end if

  col(1:m)  = a(1:m,j1)
  a(1:m,j1) = a(1:m,j2)
  a(1:m,j2) = col(1:m)

  return
end
subroutine i4col_unique_index ( m, n, a, unique_index )

!*****************************************************************************80
!
!! I4COL_UNIQUE_INDEX indexes the unique occurrence of values in an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of I4's.
!    It is regarded as an array of N columns of length M.
!
!    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
!    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
!    gathered in order, then
!
!      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!    The length of an "element" of A, and the number of "elements".
!
!    Input, integer ( kind = 4 ) A(M,N), the array.
!
!    Output, integer ( kind = 4 ) UNIQUE_INDEX(N), the unique index.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_index(n)
  integer ( kind = 4 ) unique_num

  unique_index(1:n) = -1
  unique_num = 0

  do j1 = 1, n

    if ( unique_index(j1) == -1 ) then

      unique_num = unique_num + 1
      unique_index(j1) = unique_num

      do j2 = j1 + 1, n
        if ( all ( a(1:m,j1) == a(1:m,j2) ) ) then
          unique_index(j2) = unique_num
        end if
      end do

    end if

  end do

  return
end
subroutine i4i4_sort_a ( i1, i2, j1, j2 )

!*****************************************************************************80
!
!! I4I4_SORT_A ascending sorts a pair of integers.
!
!  Discussion:
!
!    An I4I4 is a pair of integers, regarded as a single data item.
!
!    The program allows the reasonable call:
!
!      call i4i4_sort_a ( i1, i2, i1, i2 )
!
!    and this will return the reasonable result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, I2, the values to sort.
!
!    Output, integer ( kind = 4 ) J1, J2, the sorted values.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4_sort_a ( i1, i2, i1, i2 )
!
  k1 = i1
  k2 = i2

  j1 = min ( k1, k2 )
  j2 = max ( k1, k2 )

  return
end
subroutine i4i4i4_sort_a ( i1, i2, i3, j1, j2, j3 )

!*****************************************************************************80
!
!! I4I4I4_SORT_A ascending sorts a triple of integers.
!
!  Discussion:
!
!    An I4I4I4 is a triple of integers, regarded as a single data item.
!
!    The program allows the reasonable call:
!
!      call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
!
!    and this will return the reasonable result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, I2, I3, the values to sort.
!
!    Output, integer ( kind = 4 ) J1, J2, J3, the sorted values.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
!
  k1 = i1
  k2 = i2
  k3 = i3

  j1 = min ( min ( k1, k2 ), min ( k2, k3 ) )
  j2 = min ( max ( k1, k2 ), &
       min ( max ( k2, k3 ), max ( k3, k1 ) ) )
  j3 = max ( max ( k1, k2 ), max ( k2, k3 ) )

  return
end
subroutine i4int_to_r4int ( imin, imax, i, rmin, rmax, r )

!*****************************************************************************80
!
!! I4INT_TO_R4INT maps an I4INT to an R4INT.
!
!  Discussion:
!
!    The formula used is:
!
!      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IMIN, IMAX, the range.
!
!    Input, integer ( kind = 4 ) I, the integer to be converted.
!
!    Input, real ( kind = 4 ) RMIN, RMAX, the range.
!
!    Output, real ( kind = 4 ) R, the corresponding value in [RMIN,RMAX].
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 4 ) r
  real ( kind = 4 ) rmax
  real ( kind = 4 ) rmin

  if ( imax == imin ) then

    r = 0.5E+00 * ( rmin + rmax )

  else

    r = ( real ( imax - i,        kind = 4 ) * rmin   &
        + real (        i - imin, kind = 4 ) * rmax ) &
        / real ( imax     - imin, kind = 4 )

  end if

  return
end
subroutine i4int_to_r8int ( imin, imax, i, rmin, rmax, r )

!*****************************************************************************80
!
!! I4INT_TO_R8INT maps an I4INT to an R8INT.
!
!  Discussion:
!
!    The formula used is:
!
!      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IMIN, IMAX, the range.
!
!    Input, integer ( kind = 4 ) I, the integer to be converted.
!
!    Input, real ( kind = 8 ) RMIN, RMAX, the range.
!
!    Output, real ( kind = 8 ) R, the corresponding value in [RMIN,RMAX].
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) r
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmin

  if ( imax == imin ) then

    r = 0.5D+00 * ( rmin + rmax )

  else

    r = ( real ( imax - i,        kind = 8 ) * rmin   &
        + real (        i - imin, kind = 8 ) * rmax ) &
        / real ( imax     - imin, kind = 8 )

  end if

  return
end
subroutine i4list_print ( n, first, list_num, list, title )

!*****************************************************************************80
!
!! I4LIST_PRINT prints an I4LIST.
!
!  Discussion:
!
!    An I4LIST is a list of integers grouped into N segments.
!    An index vector locates the first entry of each segment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of segments.
!
!    Input, integer ( kind = 4 ) FIRST(N+1), indexes the first entry
!    of each segment.
!
!    Input, integer ( kind = 4 ) LIST_NUM, the number of entries.
!
!    Input, integer ( kind = 4 ) LIST(LIST_NUM), the data.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) list_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) first(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) list(list_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n

    do jlo = first(i), first(i+1) - 1, 5
      jhi = min ( jlo + 4, first(i+1) - 1 )
      if ( jlo == first(i) ) then
        write ( *, '(i5,a,5(2x,i8))' ) i, ':', list(jlo:jhi)
      else
        write ( *, '(6x,  5(2x,i8))' )         list(jlo:jhi)
      end if
    end do

  end do

  return
end
subroutine i4mat_border_add ( m, n, table, table2 )

!*****************************************************************************80
!
!! I4MAT_BORDER_ADD adds a "border" to an I4MAT.
!
!  Discussion:
!
!    We suppose the input data gives values of a quantity on nodes
!    in the interior of a 2D grid, and we wish to create a new table
!    with additional positions for the nodes that would be on the
!    border of the 2D grid.
!
!                  0 0 0 0 0 0
!      * * * *     0 * * * * 0
!      * * * * --> 0 * * * * 0
!      * * * *     0 * * * * 0
!                  0 0 0 0 0 0
!
!    The illustration suggests the situation in which a 3 by 4 array
!    is input, and a 5 by 6 array is to be output.
!
!    The old data is shifted to its correct positions in the new array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input,integer TABLE(M,N), the table data.
!
!    Output, integer ( kind = 4 ) TABLE2(M+2,N+2), the augmented table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) table(m,n)
  integer ( kind = 4 ) table2(m+2,n+2)

  table2(1,1:n+2) = 0
  table2(m+2,1:n+2) = 0
  table2(2:m+1,1) = 0
  table2(2:m+1,n+2) = 0

  table2(2:m+1,2:n+1) = table(1:m,1:n)

  return
end
subroutine i4mat_border_cut ( m, n, table, table2 )

!*****************************************************************************80
!
!! I4MAT_BORDER_CUT cuts the "border" of an I4MAT.
!
!  Discussion:
!
!    We suppose the input data gives values of a quantity on nodes
!    on a 2D grid, and we wish to create a new table corresponding only
!    to those nodes in the interior of the 2D grid.
!
!      0 0 0 0 0 0
!      0 * * * * 0    * * * *
!      0 * * * * 0 -> * * * *
!      0 * * * * 0    * * * *
!      0 0 0 0 0 0
!
!    The illustration suggests the situation in which a 5 by 6 array
!    is input, and a 3 by 4 array is to be output.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
!
!    Output, integer ( kind = 4 ) TABLE2(M-2,N-2), the new table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) table(m,n)
  integer ( kind = 4 ) table2(m-2,n-2)

  if ( m <= 2 .or. n <= 2 ) then
    return
  end if

  table2(1:m-2,1:n-2) = table(2:m-1,2:n-1)

  return
end
subroutine i4mat_copy ( m, n, a1, a2 )

!*****************************************************************************80
!
!! I4MAT_COPY copies an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A1(M,N), the matrix to be copied.
!
!    Output, integer ( kind = 4 ) A2(M,N), the copied matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(m,n)
  integer ( kind = 4 ) a2(m,n)

  a2(1:m,1:n) = a1(1:m,1:n)

  return
end
subroutine i4mat_elim ( m, n, a )

!*****************************************************************************80
!
!! I4MAT_ELIM carries out exact Gauss elimination on an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 1999
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
!    Input/output, integer ( kind = 4 ) A(M,N).  On input, the M by N matrix to
!    be Gauss eliminated.  On output, the Gauss-eliminated matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) amax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol(n)
  integer ( kind = 4 ) ifact
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imult
  integer ( kind = 4 ) irow(m)
  integer ( kind = 4 ) iswap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jmult
!
!  Initialize the swap parity counter.
!
  iswap = 1
!
!  For each column JCOL...
!
  do jcol = 1, min ( m, n )
!
!  Find the maximum element in rows JCOL through M.
!
    amax = abs ( a(jcol,jcol) )
    imax = jcol

    do i = jcol + 1, m
      if ( amax < abs ( a(i,jcol) ) ) then
        amax = abs ( a(i,jcol) )
        imax = i
      end if
    end do
!
!  If the maximum entry is nonzero, then...
!
    if ( amax /= 0 ) then
!
!  If the maximum entry does not occur in row JCOL, then swap rows.
!
      if ( imax /= jcol ) then
        iswap = - iswap
        call i4vec_swap ( n, a(jcol,1:n), a(imax,1:n) )
      end if
!
!  Eliminate all nonzero entries in column JCOL, below the diagonal entry.
!
      do i = jcol + 1, m

        if ( a(i,jcol) /= 0 ) then

          jmult = a(i,jcol)
          imult = a(jcol,jcol)
          ifact = i4_gcd ( imult, jmult )
          imult = imult / ifact
          jmult = jmult / ifact

          do j = jcol, n
            a(i,j) = jmult * a(jcol,j) - imult * a(i,j)
          end do

        end if

      end do
!
!  Remove any row or column factors.
!
      call i4mat_red ( m, n, a, irow, icol )

    end if

  end do

  return
end
subroutine i4mat_flip_cols ( m, n, a )

!*****************************************************************************80
!
!! I4MAT_FLIP_COLS swaps the columns of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an integer matrix.
!
!    To "flip" the columns of an I4MAT is to start with something like
!
!      11 12 13 14 15
!      21 22 23 24 25
!      31 32 33 34 35
!      41 42 43 44 45
!      51 52 53 54 55
!
!    and return
!
!      15 14 13 12 11
!      25 24 23 22 21
!      35 34 33 32 31
!      45 44 43 42 41
!      55 54 53 52 51
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), the matrix whose columns
!    are to be flipped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) b(m)
  integer ( kind = 4 ) j

  do j = 1, n / 2
    b(1:m      ) = a(1:m,    j)
    a(1:m,    j) = a(1:m,n+1-j)
    a(1:m,n+1-j) = b(1:m)
  end do

  return
end
subroutine i4mat_flip_rows ( m, n, a )

!*****************************************************************************80
!
!! I4MAT_FLIP_ROWS swaps the rows of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an integer matrix.
!
!    To "flip" the rows of an I4MAT is to start with something like
!
!      11 12 13 14 15
!      21 22 23 24 25
!      31 32 33 34 35
!      41 42 43 44 45
!      51 52 53 54 55
!
!    and return
!
!      51 52 53 54 55
!      41 42 43 44 45
!      31 32 33 34 35
!      21 22 23 24 25
!      11 12 13 14 15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), the matrix whose rows
!    are to be flipped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) i

  do i = 1, m / 2
    b(      1:n) = a(    i,1:n)
    a(    i,1:n) = a(m+1-i,1:n)
    a(m+1-i,1:n) = b(      1:n)
  end do

  return
end
subroutine i4mat_histogram ( m, n, a, histo_num, histo_gram )

!*****************************************************************************80
!
!! I4MAT_HISTOGRAM computes a histogram of the elements of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!    It is assumed that the entries in the vector A are nonnegative.
!    Only values between 0 and HISTO_NUM will be histogrammed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of A.
!
!    Input, integer ( kind = 4 ) A(M,N), the array to examine.
!
!    Input, integer ( kind = 4 ) HISTO_NUM, the maximum value for which a
!    histogram entry will be computed.
!
!    Output, integer ( kind = 4 ) HISTO_GRAM(0:HISTO_NUM), contains the
!    number of entries of A with the values of 0 through HISTO_NUM.
!
  implicit none

  integer ( kind = 4 ) histo_num
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) histo_gram(0:histo_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  histo_gram(0:histo_num) = 0

  do j = 1, n
    do i = 1, m

      if ( 0 <= a(i,j) .and. a(i,j) <= histo_num ) then
        histo_gram(a(i,j)) = histo_gram(a(i,j)) + 1
      end if

    end do
  end do

  return
end
subroutine i4mat_indicator ( m, n, table )

!*****************************************************************************80
!
!! I4MAT_INDICATOR sets up an "indicator" I4MAT.
!
!  Discussion:
!
!    The value of each entry suggests its location, as in:
!
!      11  12  13  14
!      21  22  23  24
!      31  32  33  34
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) TABLE(M,N), the table.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) table(m,n)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, m
    do j = 1, n
      table(i,j) = fac * i + j
    end do
  end do

  return
end
subroutine i4mat_l1_inverse ( n, a, b )

!*****************************************************************************80
!
!! I4MAT_L1_INVERSE inverts a unit lower triangular I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!    A unit lower triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's above the main diagonal.
!
!    The inverse of an integer unit lower triangular matrix is also
!    an integer unit lower triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call i4mat_l1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, integer ( kind = 4 ) A(N,N), the unit lower triangular matrix.
!
!    Output, integer ( kind = 4 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n

    do j = 1, i - 1
      b(i,j) = - dot_product ( a(i,1:i-1), b(1:i-1,j) )
    end do

    b(i,i) = 1
    b(i,i+1:n) = 0

  end do

  return
end
function i4mat_max ( m, n, a )

!*****************************************************************************80
!
!! I4MAT_MAX returns the maximum of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 2010
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
!    Input, integer ( kind = 4 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I4MAT_MAX, the maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i4mat_max

  i4mat_max = maxval ( a )

  return
end
subroutine i4mat_max_index ( m, n, a, i_max, j_max )

!*****************************************************************************80
!
!! I4MAT_MAX_INDEX returns the location of the maximum of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
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
!    Input, integer ( kind = 4 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I_MAX, J_MAX, the indices of the
!    maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max

  i_max = -1;
  j_max = -1;

  do j = 1, n
    do i = 1, m
      if ( i == 1 .and. j == 1 ) then
        i_max = i
        j_max = j
      else if ( a(i_max,j_max) < a(i,j) ) then
        i_max = i
        j_max = j
      end if
    end do
  end do

  return
end
function i4mat_min ( m, n, a )

!*****************************************************************************80
!
!! I4MAT_MIN returns the minimum of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 2010
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
!    Input, integer ( kind = 4 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I4MAT_MIN, the minimum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i4mat_min

  i4mat_min = minval ( a )

  return
end
subroutine i4mat_min_index ( m, n, a, i_min, j_min )

!*****************************************************************************80
!
!! I4MAT_MIN_INDEX returns the location of the minimum of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2005
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
!    Input, integer ( kind = 4 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I_MIN, J_MIN, the indices of the
!    minimum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_min

  i_min = -1
  j_min = -1

  do j = 1, n
    do i = 1, m
      if ( i == 1 .and. j == 1 ) then
        i_min = i
        j_min = j
      else if ( a(i,j) < a(i_min,j_min) ) then
        i_min = i
        j_min = j
      end if
    end do
  end do

  return
end
subroutine i4mat_mm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! I4MAT_MM multiplies two I4MAT's.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!    In FORTRAN90, this operation is more efficiently done by the
!    command:
!
!      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the order of the matrices.
!
!    Input, integer ( kind = 4 ) A(N1,N2), B(N2,N3), the matrices to multiply.
!
!    Output, integer ( kind = 4 ) C(N1,N3), the product matrix C = A * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  integer ( kind = 4 ) a(n1,n2)
  integer ( kind = 4 ) b(n2,n3)
  integer ( kind = 4 ) c(n1,n3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  do i = 1, n1
    do j = 1, n3
      c(i,j) = 0
      do k = 1, n2
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
  end do

  return
end
subroutine i4mat_perm ( n, a, p )

!*****************************************************************************80
!
!! I4MAT_PERM permutes the rows and columns of a square I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   30 September 2009
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf,
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, integer ( kind = 4 ) A(N,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) is the new
!    number of row and column I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) t

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_PERM - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if

  call perm_cycle ( n, p, is, nc, 1 )

  do i = 1, n

    i1 = - p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              t        = a(i1,j1)
              a(i1,j1) = it
              it       = t

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

  return
end
subroutine i4mat_perm_uniform ( n, a, seed )

!*****************************************************************************80
!
!! I4MAT_PERM_UNIFORM selects a random permutation of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!    The matrix is assumed to be square.  A single permutation is
!    applied to both rows and columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns
!    in the array.
!
!    Input/output, integer ( kind = 4 ) A(N,N), the array to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) seed
!
!  Permute the rows and columns together.
!
  do i = 1, n

    i2 = i4_uniform_ab ( i, n, seed )

    call i4vec_swap ( n, a(i2,1:n), a(i,1:n) )
    call i4vec_swap ( n, a(1:n,i2), a(1:n,i) )

  end do

  return
end
subroutine i4mat_perm2 ( m, n, a, p, q )

!*****************************************************************************80
!
!! I4MAT_PERM2 permutes the rows and columns of a rectangular I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, number of columns in the matrix.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) P(M), the row permutation.  P(I) is the
!    new number of row I.
!
!    Input, integer ( kind = 4 ) Q(N), the column permutation.  Q(I) is the
!    new number of column I.  Note that this routine allows you to pass a
!    single array as both P and Q.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(n)
  integer ( kind = 4 ) t

  call perm_check ( m, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_PERM2 - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if

  call perm_cycle ( m, p, is, nc, 1 )

  call perm_check ( n, q, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_PERM2 - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if

  if ( 0 < q(1) ) then
    call perm_cycle ( n, q, is, nc, 1 )
  end if

  do i = 1, m

    i1 = - p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( q(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( q(j1) )

              t        = a(i1,j1)
              a(i1,j1) = it
              it       = t

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                cycle
              end if

            end do

            j2 = abs ( q(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:m) = abs ( p(1:m) )

  if ( q(1) <= 0 ) then
    q(1:n) = abs ( q(1:n) )
  end if

  return
end
subroutine i4mat_perm2_uniform ( m, n, a, seed )

!*****************************************************************************80
!
!! I4MAT_PERM2_UNIFORM selects a random permutation of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!    The matrix may be rectangular.  Separate permutations are
!    applied to the rows and columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), the array to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) seed
!
!  Permute the rows.
!
  do i = 1, m
    i2 = i4_uniform_ab ( i, m, seed )
    call i4vec_swap ( n, a(i2,1:n), a(i,1:n) )
  end do
!
!  Permute the columns.
!
  do j = 1, n
    j2 = i4_uniform_ab ( j, n, seed )
    call i4vec_swap ( m, a(1:m,j2), a(1:m,j) )
  end do

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
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
  character ( len = * ) title

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
!    An I4MAT is a rectangular array of I4's.
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
  character ( len = * ) title

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
subroutine i4mat_red ( m, n, a, row, col )

!*****************************************************************************80
!
!! I4MAT_RED divides out common factors in a row or column of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix.
!
!    Input/output, integer ( kind = 4 ) A(M,N), on input, the M by N matrix
!    to be reduced.  On output, A has been reduced.  The greatest common
!    factor in any row or column is 1.
!
!    Output, integer ( kind = 4 ) ROW(M), the row factors that were divided out.
!
!    Output, integer ( kind = 4 ) COL(N), the column factors that were divided
!    out.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(n)
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row(m)

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IMAT_RED - Fatal error!'
    write ( *, '(a)' ) '  M must be greater than 0.'
    write ( *, '(a,i8)' ) '  Input M = ', m
    stop
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IMAT_RED - Fatal error!'
    write ( *, '(a)' ) '  N must be greater than 0.'
    write ( *, '(a,i8)' ) '  Input N = ', n
    stop
  end if
!
!  Remove factors common to a column.
!
  do j = 1, n
    call i4vec_red ( m, a(1:m,j), factor )
    col(j) = factor
  end do
!
!  Remove factors common to a row.
!
  do i = 1, m
    call i4vec_red ( n, a(i,1:n), factor )
    row(i) = factor
  end do

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 2004
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
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
  character ( len = 8 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4mat_u1_inverse ( n, a, b )

!*****************************************************************************80
!
!! I4MAT_U1_INVERSE inverts a unit upper triangular I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!    A unit upper triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's below the main diagonal.
!
!    The inverse of an integer unit upper triangular matrix is also
!    an integer unit upper triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call i4mat_u1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, integer ( kind = 4 ) A(N,N), the unit upper triangular matrix.
!
!    Output, integer ( kind = 4 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = n, 1, -1

    b(j+1:n,j) = 0
    b(j,j) = 1

    do i = j - 1, 1, -1
      b(i,j) = - dot_product ( a(i,i+1:j), b(i+1:j,j) )
    end do

  end do

  return
end
subroutine i4mat_uniform_ab ( m, n, a, b, seed, x )

!*****************************************************************************80
!
!! I4MAT_UNIFORM_AB returns a scaled pseudorandom I4MAT.
!
!  Discussion:
!
!    An I4MAT is a matrix of integer values.
!
!    The pseudorandom numbers will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the row and column dimensions
!    of the matrix.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(M,N), a matrix of values between A and B.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
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

      x(i,j) = value

    end do
  end do

  return
end
subroutine i4mat_zero ( m, n, a )

!*****************************************************************************80
!
!! I4MAT_ZERO zeroes out an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the row and column dimensions of the matrix.
!
!    Output, integer A(M,N), a matrix of zeroes.
!
  implicit none

  integer m
  integer n

  integer a(m,n)

  a(1:m,1:n) = 0

  return
end
subroutine i4row_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4ROW_COMPARE compares two rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    Input:
!
!    M = 3, N = 4, I = 2, J = 3
!
!    A = (
!    1  2  3  4
!    5  6  7  8
!    9 10 11 12 )
!
!    Output:
!
!    ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of M rows of vectors
!    of length N.
!
!    Input, integer ( kind = 4 ) I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row J < row I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check that I and J are legal.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is less than 1.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  else if ( m < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is less than 1.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  else if ( m < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      isgn = -1
      return
    else if ( a(j,k) < a(i,k) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4row_find_item ( m, n, a, item, row, col )

!*****************************************************************************80
!
!! I4ROW_FIND_ITEM searches the rows of an I4ROW for a given value.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the table to search.
!
!    Input, integer ( kind = 4 ) ITEM, the value to search for.
!
!    Output, integer ( kind = 4 ) ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM.  The search
!    is conducted by rows.  If the item is not found, then
!    ROW = COL = -1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  row = -1
  col = -1

  do i = 1, m
    do j = 1, n
      if ( a(i,j) == item ) then
        row = i
        col = j
        return
      end if
    end do
  end do

  return
end
subroutine i4row_find_pair_wrap ( m, n, a, item1, item2, row, col )

!*****************************************************************************80
!
!! I4ROW_FIND_PAIR_WRAP searches rows of an I4ROW for a pair of items.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!    The items must occur consecutively, with ITEM1 occurring
!    first.  However, wrapping is allowed.  That is, if ITEM1
!    occurs in the last column, and ITEM2 in the first, this
!    is also regarded as a match.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the table to search.
!
!    Input, integer ( kind = 4 ) ITEM1, ITEM2, the values to search for.
!
!    Output, integer ( kind = 4 ) ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM1 followed immediately
!    by ITEM2.  The search is conducted by rows.  If the pair of
!    items is not found, then ROW = COL = -1.  If COL = N,
!    the ITEM1 occurs in column N and ITEM2 occurs in column 1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item1
  integer ( kind = 4 ) item2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) row

  row = -1
  col = -1

  do i = 1, m
    do j = 1, n

      if ( a(i,j) == item1 ) then

        if ( j < n ) then
          jp1 = j + 1
        else
          jp1 = 1
        end if

        if ( a(i,jp1) == item2 ) then
          row = i
          col = j
          return
        end if

      end if

    end do
  end do

  return
end
subroutine i4row_max ( m, n, a, amax )

!*****************************************************************************80
!
!! I4ROW_MAX returns the maximums of the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) AMAX(M), the maximums of the rows
!    of the array.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) amax(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m

    amax(i) = a(i,1)
    do j = 2, n
      if ( amax(i) < a(i,j) ) then
        amax(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine i4row_mean ( m, n, a, mean )

!*****************************************************************************80
!
!! I4ROW_MEAN returns the means of the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of data.
!
!    Output, real ( kind = 8 ) MEAN(M), the mean of each row.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mean(m)

  do i = 1, m
    mean(i) = sum ( a(i,1:n) ) / real ( n, kind = 8 )
  end do

  return
end
subroutine i4row_min ( m, n, a, amin )

!*****************************************************************************80
!
!! I4ROW_MIN returns the minimums of the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) AMIN(M), the minimums of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) amin(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m

    amin(i) = a(i,1)
    do j = 2, n
      if ( a(i,j) < amin(i) ) then
        amin(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine i4row_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT_A ascending sorts the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!    In lexicographic order, the statement "X < Y", applied to two
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, X is less than Y if, at the first index where they
!    differ, the X value is less than the Y value.
!
!  Example:
!
!    Input:
!
!      M = 5, N = 3
!
!      A =
!        3  2  1
!        2  4  3
!        3  1  8
!        2  4  2
!        1  9  9
!
!    Output:
!
!      A =
!        1  9  9
!        2  4  2
!        2  4  3
!        3  1  8
!        3  2  1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4row_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4row_sort_d ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT_D descending sorts the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows and columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in descending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4row_compare ( m, n, a, i, j, isgn )
      isgn = -isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4row_sort2_d ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT2_D descending sorts the elements of each row of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the elements of each row of A have been sorted in descending
!    order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row
  integer ( kind = 4 ) t

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  do row = 1, m

    i = 0
    indx = 0
    isgn = 0
    j = 0
!
!  Call the external heap sorter.
!
    do

      call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
      if ( 0 < indx ) then

        t        = a(row,i)
        a(row,i) = a(row,j)
        a(row,j) = t
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(row,i) < a(row,j) ) then
          isgn = +1
        else
          isgn = -1
        end if

      else if ( indx == 0 ) then

        exit

      end if

    end do

  end do

  return
end
subroutine i4row_sorted_unique ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4ROW_SORTED_UNIQUE keeps unique elements in an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!    The rows of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    M rows of data.  On output, the first UNIQUE_NUM rows
!    contain the unique rows.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  i1 = 1

  do i2 = 2, m

    if ( any ( a(i1,1:n) /= a(i2,1:n) ) ) then
      i1 = i1 + 1
      a(i1,1:n) = a(i2,1:n)
    end if

  end do

  unique_num = i1

  return
end
subroutine i4row_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4ROW_SORTED_UNIQUE_COUNT counts unique elements in an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!    The rows of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    M rows of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  i1 = 1

  do i2 = 2, m

    if ( any ( a(i1,1:n) /= a(i2,1:n) ) ) then
      unique_num = unique_num + 1
      i1 = i2
    end if

  end do

  return
end
subroutine i4row_sum ( m, n, a, rsum )

!*****************************************************************************80
!
!! I4ROW_SUM returns the sums of the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of data.
!
!    Output, integer ( kind = 4 ) RSUM(M), the sum of the entries
!    of each row of TABLE.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) rsum(m)

  do i = 1, m
    rsum(i) = sum ( a(i,1:n) )
  end do

  return
end
subroutine i4row_swap ( m, n, a, i1, i2 )

!*****************************************************************************80
!
!! I4ROW_SWAP swaps two rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of data.
!
!    Input, integer ( kind = 4 ) I1, I2, the two rows to swap.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) row(n)
!
!  Check.
!
  if ( i1 < 1 .or. m < i1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I1 is out of range.'
    stop
  end if

  if ( i2 < 1 .or. m < i2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I2 is out of range.'
    stop
  end if

  if ( i1 == i2 ) then
    return
  end if

  row(1:n)  = a(i1,1:n)
  a(i1,1:n) = a(i2,1:n)
  a(i2,1:n) = row(1:n)

  return
end
subroutine i4row_variance ( m, n, a, variance )

!*****************************************************************************80
!
!! I4ROW_VARIANCE returns the variances of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of I4's, regarded
!    as an array of M rows of length N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the array of data.
!
!    Output, real ( kind = 8 ) VARIANCE(M), the variance of each row.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance(m)

  if ( n < 2 ) then

    variance(1:m) = 0.0D+00

  else

    do i = 1, m

      mean = sum ( a(i,1:n) ) / real ( n, kind = 8 )

      variance(i) = 0.0D+00
      do j = 1, n
        variance(i) = variance(i) + ( real ( a(i,j), kind = 8 ) - mean )**2
      end do

      variance(i) = variance(i) / real ( n - 1, kind = 8 )

    end do

  end if

  return
end
subroutine i4vec_add ( n, a, b, c )

!*****************************************************************************80
!
!! I4VEC_ADD computes C = A + B for I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) A(N), the first vector.
!
!    Input, integer ( kind = 4 ) B(N), the second vector.
!
!    Output, integer ( kind = 4 ) C(N), the sum of the vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) c(n)

  c(1:n) = a(1:n) + b(1:n)

  return
end
function i4vec_all_nonpositive ( n, a )

!*****************************************************************************80
!
!! I4VEC_ALL_NONPOSITIVE: ( all ( A <= 0 ) ) for I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, logical I4VEC_ALL_NONPOSITIVE is TRUE if all entries
!    of A are less than or equal to zero.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical i4vec_all_nonpositive

  i4vec_all_nonpositive = all ( a(1:n) <= 0 )

  return
end
subroutine i4vec_amax ( n, a, aamax )

!*****************************************************************************80
!
!! I4VEC_AMAX returns the largest magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be searched.
!
!    Output, integer ( kind = 4 ) AAMAX, the value of the entry of
!    largest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aamax
  integer ( kind = 4 ) i

  if ( n <= 0 ) then

    aamax = 0

  else

    aamax = abs ( a(1) )

    do i = 2, n
      aamax = max ( aamax, abs ( a(i) ) )
    end do

  end if

  return
end
subroutine i4vec_amax_index ( n, a, amax_index )

!*****************************************************************************80
!
!! I4VEC_AMAX_INDEX returns the index of the largest magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be searched.
!
!    Output, integer ( kind = 4 ) AMAX_INDEX, the index of the entry
!    of largest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aamax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) amax_index

  if ( n <= 0 ) then

    amax_index = 0

  else

    aamax = abs ( a(1) )
    amax_index = 1

    do i = 2, n

      if ( aamax < abs ( a(i) ) ) then
        aamax = abs ( a(i) )
        amax_index = i
      end if

    end do

  end if

  return
end
subroutine i4vec_amin ( n, a, aamin )

!*****************************************************************************80
!
!! I4VEC_AMIN returns the smallest magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries to be checked.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be checked.
!
!    Output, integer AAMIN, the value of the smallest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aamin
  integer ( kind = 4 ) i

  if ( n <= 0 ) then

    aamin = 0

  else

    aamin = abs ( a(1) )

    do i = 2, n
      aamin = min ( aamin, abs ( a(i) ) )
    end do

  end if

  return
end
subroutine i4vec_amin_index ( n, a, amin_index )

!*****************************************************************************80
!
!! I4VEC_AMIN_INDEX returns the index of the smallest magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries to be checked.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be checked.
!
!    Output, integer ( kind = 4 ) AMIN_INDEX, the entry of the smallest
!    magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aamin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) amin_index

  if ( n <= 0 ) then

    amin_index = 0

  else

    aamin = a(1)
    amin_index = 1

    do i = 2, n

      if ( abs ( a(i) ) < aamin ) then
        aamin = abs ( a(i) )
        amin_index = i
      end if

    end do

  end if

  return
end
subroutine i4vec_aminz ( n, a, aminz )

!*****************************************************************************80
!
!! I4VEC_AMINZ returns the smallest nonzero magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries to be checked.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be checked.
!
!    Output, integer ( kind = 4 ) AMINZ, the value of the smallest nonzero
!    magnitude.  If all entries are zero, AMINZ is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aminz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iset

  aminz = 0
  iset = 0

  do i = 1, n

    if ( a(i) /= 0 ) then

      if ( iset == 0 ) then
        aminz = abs ( a(i) )
        iset = 1
      else
        aminz = min ( aminz, abs ( a(i) ) )
      end if

    end if

  end do

  return
end
subroutine i4vec_aminz_index ( n, a, aminz_index )

!*****************************************************************************80
!
!! I4VEC_AMINZ_INDEX returns the smallest nonzero magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries to be checked.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be checked.
!
!    Output, integer ( kind = 4 ) AMINZ_INDEX, the entry of the smallest
!    nonzero magnitude.  If all entries are zero, AMINZ_INDEX is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aminz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) aminz_index

  aminz = 0
  aminz_index = 0

  do i = 1, n

    if ( a(i) /= 0 ) then

      if ( aminz_index == 0 .or. abs ( a(i) ) < aminz ) then
        aminz = abs ( a(i) )
        aminz_index = i
      end if

    end if

  end do

  return
end
function i4vec_any_lt ( n, a, b )

!*****************************************************************************80
!
!! I4VEC_ANY_LT: ( any ( A < B ) ) for I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) A(N), the first vector.
!
!    Input, integer ( kind = 4 ) B(N), the second vector.
!
!    Output, logical I4VEC_ANY_LT is TRUE if any entry
!    of A is less than the corresponding entry of B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  logical i4vec_any_lt

  i4vec_any_lt = any ( a(1:n) < b(1:n) )

  return
end
function i4vec_any_negative ( n, a )

!*****************************************************************************80
!
!! I4VEC_ANY_NEGATIVE: ( any A < 0 ) for I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, logical I4VEC_ANY_NEGATIVE is TRUE if any entry is negative.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical i4vec_any_negative

  i4vec_any_negative = any ( a(1:n) < 0 )

  return
end
function i4vec_any_nonzero ( n, a )

!*****************************************************************************80
!
!! I4VEC_ANY_NONZERO: ( any A nonzero ) for I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, logical I4VEC_ANY_NONZERO is TRUE if any entry is nonzero.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical i4vec_any_nonzero

  i4vec_any_nonzero = any ( a(1:n) /= 0 )

  return
end
subroutine i4vec_ascend_sub ( n, a, length, sub )

!*****************************************************************************80
!
!! I4VEC_ASCEND_SUB computes the longest ascending subsequence of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The subsequence is required to be strictly increasing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be examined.
!
!    Output, integer ( kind = 4 ) LENGTH, the length of the longest
!    increasing subsequence.
!
!    Output, integer ( kind = 4 ) SUB(N), contains in entries 1 through LENGTH
!    a longest increasing subsequence of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) length
  integer ( kind = 4 ) sub(n)
  integer ( kind = 4 ) top(n)
  integer ( kind = 4 ) top_prev(n)

  top(1:n) = 0
  top_prev(1:n) = 0
  sub(1:n) = 0

  if ( n <= 0 ) then
    length = 0
    return
  end if

  length = 0

  do i = 1, n

    k = -1

    do j = 1, length
      if ( a(i) <= a(top(j)) ) then
        k = j
        exit
      end if
    end do

    if ( k == -1 ) then
      length = length + 1
      k = length
    end if

    top(k) = i

    if ( 1 < k ) then
      top_prev(i) = top(k-1)
    else
      top_prev(i) = 0
    end if

  end do
!
!  Extract the subsequence.
!
  j = top(length)
  sub(length) = a(j)

  do i = length - 1, 1, -1
    j = top_prev(j)
    sub(i) = a(j)
  end do

  return
end
function i4vec_ascends ( n, x )

!*****************************************************************************80
!
!! I4VEC_ASCENDS determines if an I4VEC is (weakly) ascending.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    X = ( -8, 1, 2, 3, 7, 7, 9 )
!
!    I4VEC_ASCENDS = TRUE
!
!    The sequence is not required to be strictly ascending.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), the array to be examined.
!
!    Output, logical I4VEC_ASCENDS, is TRUE if the entries of X ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical i4vec_ascends
  integer ( kind = 4 ) x(n)

  i4vec_ascends = .false.

  do i = 1, n - 1
    if ( x(i+1) < x(i) ) then
      return
    end if
  end do

  i4vec_ascends = .true.

  return
end
subroutine i4vec_axpy ( n, ia, x, incx, y, incy )

!*****************************************************************************80
!
!! I4VEC_AXPY adds a scaled multiple of one I4VEC to another.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    If X and Y are simple vectors, then IAXPY is equivalent to:
!
!      DO I = 1, N
!        Y(I) = Y(I) + IA * X(I)
!      END DO
!
!    However, by using the increments correctly, IAXPY can also be used
!    to manipulate rows or columns of matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of X and Y.
!
!    Input, integer ( kind = 4 ) IA, the scalar value by which each entry
!    of X is multiplied before being added to Y.
!
!    Input, integer ( kind = 4 ) X(*), the vector, a multiple of which is to be
!    added to Y.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
!    Input/output, integer ( kind = 4 ) Y(*).
!    On output, each entry of Y has been increased by
!    IA times the corresponding entry of X.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    entries of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) indy
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(*)
  integer ( kind = 4 ) y(*)

  indx = 1
  indy = 1

  do i = 1, n

    y(indy) = y(indy) + ia * x(indx)

    indx = indx + incx
    indy = indy + incy

  end do

  return
end
subroutine i4vec_bracket ( n, a, xval, left, right )

!*****************************************************************************80
!
!! I4VEC_BRACKET searches a sorted I4VEC for successive brackets of a value.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    If the values in the vector are thought of as defining intervals
!    on the number line, then this routine searches for the interval
!    containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, integer ( kind = 4 ) A(N), an array that has been sorted
!    into ascending order.
!
!    Input, integer ( kind = 4 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    In the most common case, 1 <= LEFT < LEFT + 1 = RIGHT <= N,
!    and A(LEFT) <= XVAL <= A(RIGHT).
!
!    Special cases:
!      Value is less than all data values:
!        LEFT = -1, RIGHT = 1, and XVAL < A(RIGHT).
!      Value is greater than all data values:
!        LEFT = N, RIGHT = -1, and A(LEFT) < XVAL.
!      Value is equal to a data value:
!        LEFT = RIGHT, and A(LEFT) = A(RIGHT) = XVAL.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) high
  integer ( kind = 4 ) left
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) right
  integer ( kind = 4 ) xval
!
!  XVAL < A(1).
!
  if ( xval < a(1) ) then
    left = -1
    right = 1
!
!  A(N) < XVAL.
!
  else if ( a(n) < xval ) then
    left = n
    right = -1
!
!  N = 1
!
  else if ( n == 1 ) then
    left = 1
    right = 1
!
!  A(1) <= XVAL <= A(N).
!
  else

    low = 1
    high = n - 1

    do

      mid = ( low + high ) / 2

      if ( high < low ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_BRACKET - Fatal error!'
        write ( *, '(a)' ) '  Algorithm or data failure.'
        stop
      end if

      if ( a(mid) == xval ) then
        left = mid
        right = mid
        exit
      else if ( a(mid+1) == xval ) then
        left = mid + 1
        right = mid + 1
        exit
      else if ( a(mid) < xval .and. xval < a(mid+1) ) then
        left = mid
        right = mid + 1
        exit
      else if ( a(mid+1) < xval ) then
        low = mid + 1
      else if ( xval < a(mid) ) then
        high = mid - 1
      end if

    end do

  end if

  return
end
subroutine i4vec_compare ( n, a1, a2, isgn )

!*****************************************************************************80
!
!! I4VEC_COMPARE compares two I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2, 6, 2 )
!      A2 = ( 2, 8, 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A2 < A1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) k

  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = -1
      return
    else if ( a2(k) < a1(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_COPY copies an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), the vector to be copied.
!
!    Output, integer ( kind = 4 ) A2(N), a copy of A1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)

  a2(1:n) = a1(1:n)

  return
end
subroutine i4vec_cum ( n, a, a_cum )

!*****************************************************************************80
!
!! I4VEC_CUM computes the cumulutive sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    Input:
!
!      A = (/ 1, 2, 3, 4 /)
!
!    Output:
!
!      A_CUM = (/ 1, 3, 6, 10 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be summed.
!
!    Output, integer ( kind = 4 ) A_CUM(N), the cumulative sum of the
!    entries of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_cum(n)
  integer ( kind = 4 ) i

  a_cum(1) = a(1)

  do i = 2, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
subroutine i4vec_cum0 ( n, a, a_cum )

!*****************************************************************************80
!
!! I4VEC_CUM0 computes the cumulutive sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    This routine returns a vector of length N+1, with the first value
!    being 0.
!
!  Example:
!
!    Input:
!
!      A = (/ 1, 2, 3, 4 /)
!
!    Output:
!
!      A_CUM = (/ 0, 1, 3, 6, 10 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be summed.
!
!    Output, integer ( kind = 4 ) A_CUM(0:N), the cumulative sum of the
!    entries of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_cum(0:n)
  integer ( kind = 4 ) i

  a_cum(0) = 0

  do i = 1, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
function i4vec_descends ( n, x )

!*****************************************************************************80
!
!! I4VEC_DESCENDS determines if an I4VEC is decreasing.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    X = ( 9, 7, 7, 3, 2, 1, -8 )
!
!    I4VEC_DESCENDS = TRUE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), the array to be examined.
!
!    Output, logical I4VEC_DESCENDS, is TRUE if the entries of X descend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical i4vec_descends
  integer ( kind = 4 ) x(n)

  i4vec_descends = .false.

  do i = 1, n - 1
    if ( x(i) < x(i+1) ) then
      return
    end if
  end do

  i4vec_descends = .true.

  return
end
subroutine i4vec_direct_product ( factor_index, factor_order, factor_value, &
  factor_num, point_num, x )

!*****************************************************************************80
!
!! I4VEC_DIRECT_PRODUCT creates a direct product of I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out that task for the abscissas X.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      X(1:4) = ( 1, 2, 3, 4 )
!
!    Rule 2:
!      Order = 3
!      X(1:3) = ( 10, 20, 30 )
!
!    Rule 3:
!      Order = 2
!      X(1:2) = ( 100, 200 )
!
!    Product Rule:
!      Order = 24
!      X(1:24) =
!        ( 1, 10, 100 )
!        ( 2, 10, 100 )
!        ( 3, 10, 100 )
!        ( 4, 10, 100 )
!        ( 1, 20, 100 )
!        ( 2, 20, 100 )
!        ( 3, 20, 100 )
!        ( 4, 20, 100 )
!        ( 1, 30, 100 )
!        ( 2, 30, 100 )
!        ( 3, 30, 100 )
!        ( 4, 30, 100 )
!        ( 1, 10, 200 )
!        ( 2, 10, 200 )
!        ( 3, 10, 200 )
!        ( 4, 10, 200 )
!        ( 1, 20, 200 )
!        ( 2, 20, 200 )
!        ( 3, 20, 200 )
!        ( 4, 20, 200 )
!        ( 1, 30, 200 )
!        ( 2, 30, 200 )
!        ( 3, 30, 200 )
!        ( 4, 30, 200 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, integer ( kind = 4 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, integer X(FACTOR_NUM,POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of
!    values to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values
!    to set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  integer ( kind = 4 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  integer ( kind = 4 ) x(factor_num,point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    x(1:factor_num,1:point_num) = 0
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      x(factor_index,start:start+contig-1) = factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine i4vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! I4VEC_DIRECT_PRODUCT2 creates a direct product of I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, integer ( kind = 4 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, integer ( kind = 4 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of
!    values to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values to
!    set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  integer ( kind = 4 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  integer ( kind = 4 ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = 1
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
function i4vec_dot_product ( n, x, y )

!*****************************************************************************80
!
!! I4VEC_DOT_PRODUCT computes the dot product of two I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), Y(N), the arrays.
!
!    Output, integer ( kind = 4 ) I4VEC_DOT_PRODUCT, the dot product of X and Y.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i4vec_dot_product
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) y(n)

  i4vec_dot_product = dot_product ( x(1:n), y(1:n) )

  return
end
function i4vec_eq ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_EQ is true if two I4VECs are equal.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), two vectors to compare.
!
!    Output, logical I4VEC_EQ, is TRUE if every pair of elements A1(I)
!    and A2(I) are equal, and FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  logical i4vec_eq

  i4vec_eq = ( all ( a1(1:n) == a2(1:n) ) )

  return
end
function i4vec_even_all ( n, a )

!*****************************************************************************80
!
!! I4VEC_EVEN_ALL is TRUE if all entries of an I4VEC are even.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, logical I4VEC_EVEN_ALL, TRUE if all entries are even.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical i4vec_even_all

  i4vec_even_all = all ( mod ( a(1:n), 2 ) == 0 )

  return
end
function i4vec_even_any ( n, a )

!*****************************************************************************80
!
!! I4VEC_EVEN_ANY is TRUE if any entry of an I4VEC is even.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, logical I4VEC_EVEN_ANY, TRUE if any entry is even.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical i4vec_even_any

  i4vec_even_any = any ( mod ( a(1:n), 2 ) == 0 )

  return
end
subroutine i4vec_find ( n, a, value, location )

!*****************************************************************************80
!
!! I4VEC_FIND finds the first occurrence of a value in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Input, integer ( kind = 4 ) VALUE, the value being sought.
!
!    Output, integer ( kind = 4 ) LOCATION, the first location in A where 
!    VALUE occurs, or -1 if VALUE never occurs.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) location
  integer ( kind = 4 ) value

  location = -1

  do i = 1, n

    if ( a(i) == value ) then
      location = i
      return
    end if

  end do

  return
end
subroutine i4vec_first_index ( n, a, first_index )

!*****************************************************************************80
!
!! I4VEC_FIRST_INDEX indexes the first occurrence of values in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    For element A(I) of the vector, FIRST_INDEX(I) is the index in A of
!    the first occurrence of the value A(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) FIRST_INDEX(N), the first occurrence index.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) first_index(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  first_index(1:n) = -1

  do i = 1, n

    if ( first_index(i) == -1 ) then

      first_index(i) = i

      do j = i + 1, n
        if ( a(i) == a(j) ) then
          first_index(j) = i
        end if
      end do

    end if

  end do

  return
end
subroutine i4vec_frac ( n, a, k, frac )

!*****************************************************************************80
!
!! I4VEC_FRAC searches for the K-th smallest element in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    Hoare's algorithm is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), array to search.  On output,
!    the elements of A have been somewhat rearranged.
!
!    Input, integer ( kind = 4 ) K, the fractile to be sought.  If K = 1, the
!    minimum entry is sought.  If K = N, the maximum is sought.
!    Other values of K search for the entry which is K-th in size.
!    K must be at least 1, and no greater than N.
!
!    Output, integer ( kind = 4 ) FRAC, the value of the K-th fractile of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) frac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iryt
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left
  integer ( kind = 4 ) t

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal N < K, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( iryt <= left ) then
      frac = a(k)
      exit
    end if

    ix = a(k)
    i = left
    j = iryt

    do

      if ( j < i ) then

        if ( j < k ) then
          left = i
        end if

        if ( k < i ) then
          iryt = j
        end if

        exit

      end if
!
!  Find I so that IX <= A(I).
!
      do while ( a(i) < ix )
        i = i + 1
      end do
!
!  Find J so that A(J) <= IX.
!
      do while ( ix < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then

        t    = a(i)
        a(i) = a(j)
        a(j) = t

        i = i + 1
        j = j - 1

      end if

    end do

  end do

  return
end
subroutine i4vec_gcd ( n, v, gcd )

!*****************************************************************************80
!
!! I4VEC_GCD returns the greatest common divisor of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The value GCD returned has the property that it is the greatest integer
!    which evenly divides every entry of V.
!
!    The entries in V may be negative.
!
!    Any zero entries in V are ignored.  If all entries of V are zero,
!    GCD is returned as 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of V.
!
!    Input, integer ( kind = 4 ) V(N), the vector.
!
!    Output, integer ( kind = 4 ) GCD, the greatest common divisor of V.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) gcd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) v(n)

  gcd = 0

  do i = 1, n

    if ( v(i) /= 0 ) then
      if ( gcd == 0 ) then
        gcd = abs ( v(i) )
      else
        gcd = i4_gcd ( gcd, v(i) )
      end if
    end if

  end do
!
!  If GCD is 0, that can only happen because all entries of V are zero.
!
  if ( gcd == 0 ) then
    gcd = 1
  end if

  return
end
subroutine i4vec_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_A reorders an I4VEC into an ascending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
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
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( key <= a(m) ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an I4VEC into an descending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_heap_d_extract ( n, a, value )

!*****************************************************************************80
!
!! I4VEC_HEAP_D_EXTRACT extracts the maximum value from a heap descending I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In other words, the routine finds the maximum value in the
!    heap, returns that value to the user, deletes that value from
!    the heap, and restores the heap to its proper form.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input/output, integer ( kind = 4 ) A(N), the heap.
!
!    Output, integer ( kind = 4 ) VALUE, the item of maximum value, which has
!    been removed from the heap.
!
  implicit none

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_HEAP_D_EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  The heap is empty.'
    stop
  end if
!
!  Get the maximum value.
!
  value = a(1)

  if ( n == 1 ) then
    n = 0
    return
  end if
!
!  Shift the last value down.
!
  a(1) = a(n)
!
!  Restore the heap structure.
!
  n = n - 1
  call i4vec_sort_heap_d ( n, a )

  return
end
subroutine i4vec_heap_d_insert ( n, a, value )

!*****************************************************************************80
!
!! I4VEC_HEAP_D_INSERT inserts a new I4 into a heap descending I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input/output, integer ( kind = 4 ) A(N), the heap.
!
!    Input, integer ( kind = 4 ) VALUE, the value to be inserted.
!
  implicit none

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) parent
  integer ( kind = 4 ) value

  n = n + 1
  i = n

  do while ( 1 < i )

    parent = i / 2

    if ( value <= a(parent) ) then
      exit
    end if

    a(i) = a(parent)
    i = parent

  end do

  a(i) = value

  return
end
subroutine i4vec_heap_d_max ( n, a, value )

!*****************************************************************************80
!
!! I4VEC_HEAP_D_MAX returns the maximum value in a heap descending I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input, integer ( kind = 4 ) A(N), the heap.
!
!    Output, integer ( kind = 4 ) VALUE, the maximum value in the heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) value

  value = a(1)

  return
end
subroutine i4vec_histogram ( n, a, histo_num, histo_gram )

!*****************************************************************************80
!
!! I4VEC_HISTOGRAM computes a histogram of the elements of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    It is assumed that the entries in the vector A are nonnegative.
!    Only values between 0 and HISTO_NUM will be histogrammed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the array to examine.
!
!    Input, integer ( kind = 4 ) HISTO_NUM, the maximum value for which a
!    histogram entry will be computed.
!
!    Output, integer ( kind = 4 ) HISTO_GRAM(0:HISTO_NUM), contains the
!    number of entries of A with the values of 0 through HISTO_NUM.
!
  implicit none

  integer ( kind = 4 ) histo_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) histo_gram(0:histo_num)
  integer ( kind = 4 ) i

  histo_gram(0:histo_num) = 0

  do i = 1, n

    if ( 0 <= a(i) .and. a(i) <= histo_num ) then
      histo_gram(a(i)) = histo_gram(a(i)) + 1
    end if

  end do

  return
end
function i4vec_index ( n, a, aval )

!*****************************************************************************80
!
!! I4VEC_INDEX returns the first location of a given value in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be searched.
!
!    Input, integer ( kind = 4 ) AVAL, the value to be indexed.
!
!    Output, integer ( kind = 4 ) I4VEC_INDEX, the first location in A which
!    has the value AVAL, or -1 if no such index exists.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_index

  do i = 1, n
    if ( a(i) == aval ) then
      i4vec_index = i
      return
    end if
  end do

  i4vec_index = -1

  return
end
subroutine i4vec_index_delete_all ( n, x, indx, xval )

!*****************************************************************************80
!
!! I4VEC_INDEX_DELETE_ALL deletes a value in an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, integer ( kind = 4 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) equal1
  integer ( kind = 4 ) equal2
  integer ( kind = 4 ) get
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) put
  integer ( kind = 4 ) x(*)
  integer ( kind = 4 ) xval

  if ( n < 1 ) then
    n = 0
    return
  end if

  call i4vec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    return
  end if

  equal1 = equal

  do

    if ( equal1 <= 1 ) then
      exit
    end if

    if ( x(indx(equal1-1)) /= xval ) then
      exit
    end if

    equal1 = equal1 - 1

  end do

  equal2 = equal

  do

    if ( n <= equal2 ) then
      exit
    end if

    if ( x(indx(equal2+1)) /= xval ) then
      exit
    end if

    equal2 = equal2 + 1

  end do
!
!  Discard certain X values.
!
  put = 0

  do get = 1, n

    if ( x(get) /= xval ) then
      put = put + 1
      x(put) = x(get)
    end if

  end do

  x(put+1:n) = 0
!
!  Adjust the INDX values.
!
  do equal = equal1, equal2
    do i = 1, n
      if ( indx(equal) < indx(i) ) then
        indx(i) = indx(i) - 1
      end if
    end do
  end do
!
!  Discard certain INDX values.
!
  indx(equal1:n+equal1-equal2-1) = indx(equal2+1:n)
  indx(n+equal1-equal2:n) = 0
!
!  Adjust N.
!
  n = put

  return
end
subroutine i4vec_index_delete_dupes ( n, x, indx, n2, x2, indx2 )

!*****************************************************************************80
!
!! I4VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The output quantities N2, X2, and INDX2 are computed from the
!    input quantities by sorting, and eliminating duplicates.
!
!    The output arrays should be dimensioned of size N, unless the user
!    knows in advance what the value of N2 will be.
!
!    The output arrays may be identified with the input arrays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 )N, the size of the input list.
!
!    Input, integer ( kind = 4 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Output, integer ( kind = 4 ) N2, the number of unique entries in X.
!
!    Output, integer ( kind = 4 ) X2(N2), a copy of the list which has
!    been sorted, and made unique.
!
!    Output, integer  ( kind = 4 ) INDX2(N2), the sort index of the new list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x2(n)
  integer ( kind = 4 ) x3(n)

  i = 0
  n3 = 0

  do

    i = i + 1

    if ( n < i ) then
      exit
    end if

    if ( 1 < i ) then
      if ( x(indx(i)) == x3(n3) ) then
        cycle
      end if
    end if

    n3 = n3 + 1
    x3(n3) = x(indx(i))

  end do
!
!  Copy data into output arrays.
!
  n2 = n3
  x2(1:n2) = x3(1:n3)
  call i4vec_indicator ( n2, indx2 )

  return
end
subroutine i4vec_index_delete_one ( n, x, indx, xval, n2, x2, indx2 )

!*****************************************************************************80
!
!! I4VEC_INDEX_DELETE_ONE deletes one copy of I4 from an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    If the value occurs in the list more than once, only one copy is deleted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, integer ( kind = 4 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) N2, the size of the current list.
!
!    Output, integer ( kind = 4 ) X2(N2), the list.
!
!    Output, integer ( kind = 4 ) INDX2(N2), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x2(n)
  integer ( kind = 4 ) xval

  if ( n < 1 ) then
    n2 = 0
    return
  end if

  n2 = n
  indx2(1:n2) = indx(1:n2)
  x2(1:n2) = x(1:n2)

  call i4vec_index_search ( n2, x2, indx2, xval, less, equal, more )

  if ( equal /= 0 ) then
    j = indx2(equal)
    x2(j:n2-1) = x2(j+1:n2)
    indx2(equal:n2-1) = indx2(equal+1:n2)
    do i = 1, n2 - 1
      if ( j < indx2(i) ) then
        indx2(i) = indx2(i) - 1
      end if
    end do
    n2 = n2 - 1
  end if

  return
end
subroutine i4vec_index_insert ( n, x, indx, xval )

!*****************************************************************************80
!
!! I4VEC_INDEX_INSERT inserts an I4 into an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, integer ( kind = 4 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) x(*)
  integer ( kind = 4 ) xval

  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if

  call i4vec_index_search ( n, x, indx, xval, less, equal, more )

  x(n+1) = xval
  indx(n+1:more+1:-1) = indx(n:more:-1)
  indx(more) = n + 1
  n = n + 1

  return
end
subroutine i4vec_index_insert_unique ( n, x, indx, xval )

!*****************************************************************************80
!
!! I4VEC_INDEX_INSERT_UNIQUE inserts a unique I4 into an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!    If the input value XVAL does not already occur in X, then N is increased.
!
!    Input/output, integer ( kind = 4 ) X(N), the list.
!    If the input value XVAL does not already occur in X, then it is added
!    to X.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!    If the input value XVAL does not already occur in X, then INDX is updated.
!
!    Input, integer ( kind = 4 ) XVAL, the value which will be inserted into
!    the X vector if it is not there already.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) x(*)
  integer ( kind = 4 ) xval

  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if
!
!  Does XVAL already occur in X?
!
  call i4vec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    x(n+1) = xval
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1
  end if

  return
end
subroutine i4vec_index_order ( n, x, indx )

!*****************************************************************************80
!
!! I4VEC_INDEX_ORDER sorts an I4VEC using an index vector.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The index vector itself is not modified.  Therefore, the pair
!    (X,INDX) no longer represents an index sorted vector.  If this
!    relationship is to be preserved, then simply set INDX(1:N)=(1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, integer ( kind = 4 ) X(N), the list.  On output, the list
!    has been sorted.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) y(n)

  y(1:n) = x(indx(1:n))
  x(1:n) = y(1:n)

  return
end
subroutine i4vec_index_search ( n, x, indx, xval, less, equal, more )

!*****************************************************************************80
!
!! I4VEC_INDEX_SEARCH searches for an I4 in an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, integer ( kind = 4 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, the value to be sought.
!
!    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) xhi
  integer ( kind = 4 ) xlo
  integer ( kind = 4 ) xmid
  integer ( kind = 4 ) xval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n
  xlo = x(indx(lo))
  xhi = x(indx(hi))

  if ( xval < xlo ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( xval == xlo ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  if ( xhi < xval ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( xval == xhi ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))

    if ( xval == xmid ) then
      equal = mid
      less = mid - 1
      more = mid + 1
      return
    else if ( xval < xmid ) then
      hi = mid
    else if ( xmid < xval ) then
      lo = mid
    end if

  end do

  return
end
subroutine i4vec_index_sort_unique ( n, x, n2, x2, indx2 )

!*****************************************************************************80
!
!! I4VEC_INDEX_SORT_UNIQUE creates a sorted unique index for an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, integer ( kind = 4 ) X(N), the list.
!
!    Output, integer ( kind = 4 ) N2, the number of unique elements in X.
!
!    Output, integer ( kind = 4 ) X2(N2), a list of the unique elements of X.
!
!    Output, integer ( kind = 4 ) INDX2(N2), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x2(n)

  n2 = 0

  do i = 1, n
    call i4vec_index_insert_unique ( n2, x2, indx2, x(i) )
  end do

  x2(n2+1:n) = -1
  indx2(n2+1:n) = -1

  return
end
subroutine i4vec_indexed_heap_d ( n, a, indx )

!*****************************************************************************80
!
!! I4VEC_INDEXED_HEAP_D creates a descending heap from an indexed I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
!    each referencing an entry of the data vector.
!
!    The function adjusts the index vector INDX so that, for 1 <= J <= N/2,
!    we have:
!      A(INDX(2*J))   <= A(INDX(J))
!    and
!      A(INDX(2*J+1)) <= A(INDX(J))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the index array.
!
!    Input, integer ( kind = 4 ) A(*), the data vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the index array.
!    Each entry of INDX must be a valid index for the array A.
!    On output, the indices have been reordered into a descending heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = indx(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(indx(m)) < a(indx(m+1)) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(indx(m)) <= a(key) ) then
        exit
      end if

      indx(ifree) = indx(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    indx(ifree) = key

  end do

  return
end
subroutine i4vec_indexed_heap_d_extract ( n, a, indx, indx_extract )

!*****************************************************************************80
!
!! I4VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
!    each referencing an entry of the data vector.
!
!    The routine finds the maximum value in the heap, returns that value to the
!    user, deletes that value from the heap, and restores the heap to its
!    proper form.
!
!    Note that the argument N must be a variable, which will be decremented
!    before return, and that INDX will hold one less value on output than it
!    held on input.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the
!    index vector.
!
!    Input, integer ( kind = 4 ) A(*), the data vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the index vector.
!
!    Output, integer ( kind = 4 ) INDX_EXTRACT, the index in A of the item of
!    maximum value, which has now been removed from the heap.
!
  implicit none

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) indx_extract
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_INDEXED_HEAP_D_EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  The heap is empty.'
    stop
  end if
!
!  Get the index of the maximum value.
!
  indx_extract = indx(1)

  if ( n == 1 ) then
    n = 0
    return
  end if
!
!  Shift the last index down.
!
  indx(1) = indx(n)
!
!  Restore the heap structure.
!
  n = n - 1
  call i4vec_indexed_heap_d ( n, a, indx )

  return
end
subroutine i4vec_indexed_heap_d_insert ( n, a, indx, indx_insert )

!*****************************************************************************80
!
!! I4VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
!    each referencing an entry of the data vector.
!
!    Note that the argument N must be a variable, and will be incremented before
!    return, and that INDX must be able to hold one more entry on output than
!    it held on input.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the
!    index vector.
!
!    Input, integer ( kind = 4 ) A(*), the data vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the index vector.
!
!    Input, integer ( kind = 4 ) INDX_INSERT, the index in A of the value
!    to be inserted into the heap.
!
  implicit none

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) indx_insert
  integer ( kind = 4 ) n
  integer ( kind = 4 ) parent

  n = n + 1
  i = n

  do while ( 1 < i )

    parent = i / 2

    if ( a(indx_insert) <= a(indx(parent)) ) then
      exit
    end if

    indx(i) = indx(parent)
    i = parent

  end do

  indx(i) = indx_insert

  return
end
subroutine i4vec_indexed_heap_d_max ( n, a, indx, indx_max )

!*****************************************************************************80
!
!! I4VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
!    each referencing an entry of the data vector.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the index vector.
!
!    Input, integer ( kind = 4 ) A(*), the data vector.
!
!    Input, integer ( kind = 4 ) INDX(N), the index vector.
!
!    Output, integer ( kind = 4 ) INDX_MAX, the index in A of the maximum value
!    in the heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx_max

  indx_max = indx(1)

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_insert ( n, a, pos, value )

!*****************************************************************************80
!
!! I4VEC_INSERT inserts a value into an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the array on input.
!
!    Input/output, integer ( kind = 4 ) A(N+1), the array.  On input, A is
!    assumed to contain N entries.  On output, A actually contains N+1 entries.
!
!    Input, integer ( kind = 4 ) POS, the position to be assigned the new entry.
!    1 <= POS <= N+1.
!
!    Input, integer ( kind = 4 ) VALUE, the value to be inserted at the given
!    position.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) value

  if ( pos < 1 .or. n + 1 < pos ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_INSERT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal insertion position = ', pos
    stop

  else

    do i = n + 1, pos + 1, -1
      a(i) = a(i-1)
    end do

    a(pos) = value

  end if

  return
end
function i4vec_lcm ( n, v )

!*****************************************************************************80
!
!! I4VEC_LCM returns the least common multiple of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The value LCM returned has the property that it is the smallest integer
!    which is evenly divisible by every element of V.
!
!    The entries in V may be negative.
!
!    If any entry of V is 0, then LCM is 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of V.
!
!    Input, integer ( kind = 4 ) V(N), the vector.
!
!    Output, integer ( kind = 4 ) I4VEC_LCM, the least common multiple of V.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_lcm
  integer ( kind = 4 ) i4vec_lcm
  integer ( kind = 4 ) lcm
  integer ( kind = 4 ) v(n)

  lcm = 1

  do i = 1, n

    if ( v(i) == 0 ) then
      lcm = 0
      return
    end if

    lcm = i4_lcm ( lcm, v(i) )

  end do

  i4vec_lcm = lcm

  return
end
subroutine i4vec_mask_print ( n, a, mask_num, mask, title )

!*****************************************************************************80
!
!! I4VEC_MASK_PRINT prints a masked I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MASK_NUM, the number of masked elements.
!
!    Input, integer ( kind = 4 ) MASK(MASK_NUM), the indices of the vector
!    to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) mask_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mask(mask_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Masked vector printout:'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, mask_num
    write ( *, '(2x,i8,a,1x,i8,2x,i10)' ) i, ':', mask(i), a(mask(i))
  end do

  return
end
subroutine i4vec_max ( n, a, amax )

!*****************************************************************************80
!
!! I4VEC_MAX computes the maximum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMAX, the value of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amax

  amax = maxval ( a(1:n) )

  return
end
subroutine i4vec_max_index ( n, a, max_index )

!*****************************************************************************80
!
!! I4VEC_MAX_INDEX computes the index of a maximum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    If more than one element has the maximum value, this routine returns
!    the index of the first such element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) MAX_INDEX, the index of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_index

  if ( n <= 0 ) then

    max_index = 0

  else

    amax = a(1)
    max_index = 1

    do i = 2, n

      if ( amax < a(i) ) then
        amax = a(i)
        max_index = i
      end if

    end do

  end if

  return
end
function i4vec_max_index_last ( n, x )

!*****************************************************************************80
!
!! I4VEC_MAX_INDEX_LAST returns the last maximal element location in an I4VEC
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    X = ( 5, 1, 2, 5, 0, 5, 3 )
!
!    I4VEC_MAX_INDEX_LAST = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), the array to be examined.
!
!    Output, integer ( kind = 4 ) I4VEC_MAX_INDEX_LAST, the index of the
!    last element of X of maximal value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_max_index_last
  integer ( kind = 4 ) max_last
  integer ( kind = 4 ) x(n)

  i4vec_max_index_last = 0

  do i = 1, n
    if ( i == 1 ) then
      i4vec_max_index_last = 1
      max_last = x(1)
    else if ( max_last <= x(i) ) then
      i4vec_max_index_last = i
      max_last = x(i)
    end if
  end do

  return
end
subroutine i4vec_mean ( n, a, mean )

!*****************************************************************************80
!
!! I4VEC_MEAN returns the mean of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector whose mean is desired.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  real ( kind = 8 ) mean

  mean = real ( sum ( a(1:n) ), kind = 8 ) &
       / real ( n, kind = 8 )

  return
end
subroutine i4vec_median ( n, a, median )

!*****************************************************************************80
!
!! I4VEC_MEDIAN returns the median of an unsorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    Hoare's algorithm is used.  The values of the vector are
!    rearranged by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to search.  On output,
!    the order of the elements of A has been somewhat changed.
!
!    Output, integer ( kind = 4 ) MEDIAN, the value of the median of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) median

  k = ( n + 1 ) / 2

  call i4vec_frac ( n, a, k, median )

  return
end
subroutine i4vec_merge_a ( na, a, nb, b, nc, c )

!*****************************************************************************80
!
!! I4VEC_MERGE_A merges two ascending sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The elements of A and B should be sorted in ascending order.
!
!    The elements in the output array C will also be in ascending order,
!    and unique.
!
!    The output vector C may share storage with A or B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, integer ( kind = 4 ) A(NA), the first sorted array.
!
!    Input, integer ( kind = 4 ) NB, the dimension of B.
!
!    Input, integer ( kind = 4 ) B(NB), the second sorted array.
!
!    Output, integer ( kind = 4 ) NC, the number of elements in the output
!    array.  Note that C should usually be dimensioned at least NA+NB in the
!    calling routine.
!
!    Output, integer ( kind = 4 ) C(NC), the merged unique sorted array.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb

  integer ( kind = 4 ) a(na)
  integer ( kind = 4 ) b(nb)
  integer ( kind = 4 ) c(na+nb)
  integer ( kind = 4 ) d(na+nb)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) na2
  integer ( kind = 4 ) nb2
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) order

  na2 = na
  nb2 = nb

  ja = 0
  jb = 0
  nc = 0

  call i4vec_order_type ( na2, a, order )

  if ( order < 0 .or. 2 < order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_MERGE - Fatal error!'
    write ( *, '(a)') '  The input array A is not ascending sorted.'
    stop
  end if

  call i4vec_order_type ( nb2, b, order )

  if ( order < 0 .or. 2 < order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_MERGE - Fatal error!'
    write ( *, '(a)' ) '  The input array B is not ascending sorted.'
    stop
  end if

  do
!
!  If we've used up all the entries of A, stick the rest of B on the end.
!
    if ( na2 <= ja ) then

      do j = 1, nb2 - jb
        jb = jb + 1
        if ( nc == 0 .or. d(nc) < b(jb) ) then
          nc = nc + 1
          d(nc) = b(jb)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  If we've used up all the entries of B, stick the rest of A on the end.
!
    else if ( nb2 <= jb ) then

      do j = 1, na2 - ja
        ja = ja + 1
        if ( nc == 0 .or. d(nc) < a(ja) ) then
          nc = nc + 1
          d(nc) = a(ja)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  Otherwise, if the next entry of A is smaller, that's our candidate.
!
    else if ( a(ja+1) <= b(jb+1) ) then

      ja = ja + 1
      if ( nc == 0 .or. d(nc) < a(ja) ) then
        nc = nc + 1
        d(nc) = a(ja)
      end if
!
!  ...or if the next entry of B is the smaller, consider that.
!
    else

      jb = jb + 1
      if ( nc == 0 .or. d(nc) < b(jb) ) then
        nc = nc + 1
        d(nc) = b(jb)
      end if
    end if

  end do

  return
end
subroutine i4vec_min ( n, a, amin )

!*****************************************************************************80
!
!! I4VEC_MIN computes the minimum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMIN, the value of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amin

  amin = minval ( a(1:n) )

  return
end
subroutine i4vec_min_index ( n, a, imin )

!*****************************************************************************80
!
!! I4VEC_MIN_INDEX computes the index of the minimum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) IMIN, the index of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imin

  if ( n <= 0 ) then

    imin = 0

  else

    amin = a(1)
    imin = 1

    do i = 2, n

      if ( a(i) < amin ) then
        amin = a(i)
        imin = i
      end if

    end do

  end if

  return
end
subroutine i4vec_min_mv ( m, n, u, v, w )

!*****************************************************************************80
!
!! I4VEC_MIN_MV determines U(1:N) /\ V for vectors U and a single vector V.
!
!  Discussion:
!
!    For two vectors U and V, each of length M, we define
!
!      ( U /\ V ) (I) = min ( U(I), V(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the vectors.
!
!    Input, integer ( kind = 4 ) N, the number of vectors in U.
!
!    Input, integer ( kind = 4 ) U(M,N), N vectors, each of length M.
!
!    Input, integer ( kind = 4 ) V(M), a vector of length M.
!
!    Output, integer ( kind = 4 ) W(M,N), the value of U /\ W.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) u(m,n)
  integer ( kind = 4 ) v(m)
  integer ( kind = 4 ) w(m,n)

  do j = 1, n
    do i = 1, m
      w(i,j) = min ( u(i,j), v(i) )
    end do
  end do

  return
end
function i4vec_nonzero_count ( n, a )

!*****************************************************************************80
!
!! I4VEC_NONZERO_COUNT counts the nonzero entries in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input, integer ( kind = 4 ) A(N), an array.
!
!    Output, integer ( kind = 4 ) I4VEC_NONZERO_COUNT, the number of
!    nonzero entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_nonzero_count

  i4vec_nonzero_count = 0

  do i = 1, n
    if ( a(i) /= 0 ) then
      i4vec_nonzero_count = i4vec_nonzero_count + 1
    end if
  end do

  return
end
subroutine i4vec_nonzero_first ( n, x, nz, indx )

!*****************************************************************************80
!
!! I4VEC_NONZERO_FIRST left-shifts all nonzeros in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The routine preserves the ordering of the nonzero entries.  It counts
!    the nonzeros, and returns an index vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) X(N), the vector to be shifted.
!
!    Output, integer ( kind = 4 ) NZ, the number of nonzero entries in
!    the vector.
!
!    Output, integer ( kind = 4 ) INDX(N), contains the original location
!    of each entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) x(n)

  nz = 0

  do j = 1, n
    indx(j) = j
  end do

  j = 0

  do while ( j < n )

    j = j + 1

    if ( x(j) /= 0 ) then

      nz = nz + 1

      if ( nz /= j ) then

        x(nz) = x(j)
        x(j) = 0

        k = indx(nz)
        indx(nz) = j
        indx(j) = k

      end if
    end if
  end do

  return
end
function i4vec_norm_l0 ( n, a )

!*****************************************************************************80
!
!! I4VEC_NORM_L0 returns the l0 "norm" of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The l0 "norm" simply counts the number of nonzero entries in the vector.
!    It is not a true norm, but has some similarities to one.  It is useful
!    in the study of compressive sensing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, integer ( kind = 4 ) I4VEC_NORM_L0, the value of the norm.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_norm_l0
  integer ( kind = 4 ) value

  value = 0
  do i = 1, n
    if ( a(i) /= 0 ) then
      value = value + 1
    end if
  end do

  i4vec_norm_l0 = value

  return
end
function i4vec_odd_all ( n, a )

!*****************************************************************************80
!
!! I4VEC_ODD_ALL is TRUE if all entries of an I4VEC are odd.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, logical I4VEC_ODD_ALL, TRUE if all entries are odd.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical i4vec_odd_all

  i4vec_odd_all = all ( mod ( a(1:n), 2 ) == 1 )

  return
end
function i4vec_odd_any ( n, a )

!*****************************************************************************80
!
!! I4VEC_ODD_ANY is TRUE if any entry of an I4VEC is odd.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, logical I4VEC_ODD_ANY, TRUE if any entry is odd.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical i4vec_odd_any

  i4vec_odd_any = any ( mod ( a(1:n), 2 ) == 1 )

  return
end
subroutine i4vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! I4VEC_ORDER_TYPE determines if I4VEC is (non)strictly ascending/descending.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the array.
!
!    Input, integer ( kind = 4 ) A(N), the array to be checked.
!
!    Output, integer ( kind = 4 ) ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( a(1) < a(i) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do while ( i < n )

    i = i + 1

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
function i4vec_pairwise_prime ( n, a )

!*****************************************************************************80
!
!! I4VEC_PAIRWISE_PRIME checks whether an I4VEC's entries are pairwise prime.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    Two positive integers I and J are pairwise prime if they have no common
!    factor greater than 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to check.
!
!    Input, integer ( kind = 4 ) A(N), the vector of integers.
!
!    Output, logical I4VEC_PAIRWISE_PRIME, is TRUE if the vector of integers
!    is pairwise prime.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  logical i4vec_pairwise_prime
  integer ( kind = 4 ) j

  i4vec_pairwise_prime = .false.

  do i = 1, n
    do j = i + 1, n
      if ( i4_gcd ( a(i), a(j) ) /= 1 ) then
        return
      end if
    end do
  end do

  i4vec_pairwise_prime = .true.

  return
end
subroutine i4vec_part ( n, nval, a )

!*****************************************************************************80
!
!! I4VEC_PART partitions an integer NVAL into N nearly equal parts.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    Input:
!
!      N = 5, NVAL = 17
!
!    Output:
!
!      A = ( 4, 4, 3, 3, 3 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) NVAL, the integer to be partitioned.
!    NVAL may be positive, zero, or negative.
!
!    Output, integer ( kind = 4 ) A(N), the partition of NVAL.  The entries of
!    A add up to NVAL.  The entries of A are either all equal, or
!    differ by at most 1.  The entries of A all have the same sign
!    as NVAL, and the "largest" entries occur first.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nval

  a(1:n) = 0

  if ( 0 < nval ) then

    j = 1
    do i = 1, nval
      a(j) = a(j) + 1
      j = j + 1
      if ( n < j ) then
        j = 1
      end if
    end do

  else if ( nval < 0 ) then

    j = 1
    do i = nval, -1
      a(j) = a(j) - 1
      j = j + 1
      if ( n < j ) then
        j = 1
      end if
    end do

  end if

  return
end
subroutine i4vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! I4VEC_PART_QUICK_A reorders an I4VEC as part of a quick sort.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The routine reorders the entries of A.  Using A(1) as a key,
!    all entries of A that are less than or equal to the key will
!    precede the key which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
!
!    Output:
!
!      L = 3, R = 6
!
!      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
!            -------        -------
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, integer ( kind = 4 ) A(N).  On input, the array to be
!    checked.  On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the
!    three segments.
!    Let KEY = the input value of A(1).  Then
!    I <= L                 A(I) < KEY;
!         L < I < R         A(I) = KEY;
!                 R <= I    KEY < A(I).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) key
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  integer ( kind = 4 ) t

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key = a(1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( key < a(l+1) ) then
      r = r - 1
      t      = a(r)
      a(r)   = a(l+1)
      a(l+2) = t
    else if ( a(l+1) == key ) then
      m = m + 1
      t      = a(m)
      a(m)   = a(l+1)
      a(l+1) = t
      l = l + 1
    else if ( a(l+1) < key ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(i) = a(i+m)
  end do
!
!  Out of bounds here, occasionally.
!
  l = l - m

  a(l+1:l+m) = key

  return
end
subroutine i4vec_permute ( n, p, a )

!*****************************************************************************80
!
!! I4VEC_PERMUTE permutes an I4VEC in place.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    This routine permutes an array of integer "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = (   1,   2,   3,   4,   5 )
!
!    Output:
!
!      A    = (   2,   4,   5,   1,   3 ).
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
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_temp
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine i4vec_permute_uniform ( n, a, seed )

!*****************************************************************************80
!
!! I4VEC_PERMUTE_UNIFORM randomly permutes an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  call perm_uniform ( n, base, seed, p )

  call i4vec_permute ( n, p, a )

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine i4vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! I4VEC_PRINT_PART prints "part" of an I4VEC.
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
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
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
      write ( *, '(2x,i8,a,1x,i8)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,i8)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ........'
    i = n
    write ( *, '(2x,i8,a,1x,i8)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,i8)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,i8,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end
subroutine i4vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! I4VEC_PRINT_SOME prints "some" of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
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
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices
!    to print.  The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,a,2x,i8)' ) i, ':', a(i)
  end do

  return
end
function i4vec_product ( n, a )

!*****************************************************************************80
!
!! I4VEC_PRODUCT returns the product of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In FORTRAN90, this facility is offered by the built in
!    PRODUCT function:
!
!      I4VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
!
!    In MATLAB, this facility is offered by the built in
!    PROD function:
!
!      I4VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) I4VEC_PRODUCT, the product of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4vec_product

  i4vec_product = product ( a(1:n) )

  return
end
subroutine i4vec_red ( n, a, factor )

!*****************************************************************************80
!
!! I4VEC_RED divides out common factors in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    On output, the entries of A have no common factor
!    greater than 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) A(N), the vector to be reduced.
!
!    Output, integer ( kind = 4 ) FACTOR, the common factor that was divided
!    out.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
!
!  Find the smallest nonzero value.
!
  factor = 0

  do i = 1, n

    if ( a(i) /= 0 ) then

      if ( factor == 0 ) then
        factor = abs ( a(i) )
      else
        factor = min ( factor, abs ( a(i) ) )
      end if

    end if

  end do

  if ( factor == 0 ) then
    return
  end if
!
!  Find the greatest common factor of the entire vector.
!
  do i = 1, n
    factor = i4_gcd ( a(i), factor )
  end do

  if ( factor == 1 ) then
    return
  end if
!
!  Divide out the common factor.
!
  do i = 1, n
    a(i) = a(i) / factor
  end do

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
!
!      A(1:N) = A(N:1:-1)
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)

  a(1:n) = a(n:1:-1)

  return
end
subroutine i4vec_rotate ( n, m, a )

!*****************************************************************************80
!
!! I4VEC_ROTATE rotates an I4VEC in place.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      A = ( 1, 2, 3, 4, 5 )
!
!    Output:
!
!      A = ( 4, 5, 1, 2, 3 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, integer A(N), the array to be rotated.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mcopy
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) temp
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i4_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

  do

    istart = istart + 1

    if ( n < istart ) then
      exit
    end if

    temp = a(istart)
    iget = istart
!
!  Copy the new value into the vacated entry.
!
    do

      iput = iget

      iget = iget - mcopy

      if ( iget < 1 ) then
        iget = iget + n
      end if

      if ( iget == istart ) then
        exit
      end if

      a(iput) = a(iget)
      nset = nset + 1

    end do

    a(iput) = temp
    nset = nset + 1

    if ( n <= nset ) then
      exit
    end if

  end do

  return
end
subroutine i4vec_run_count ( n, a, run_count )

!*****************************************************************************80
!
!! I4VEC_RUN_COUNT counts runs of equal values in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    A run is a sequence of equal values.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be examined.
!
!    Output, integer ( kind = 4 ) RUN_COUNT, the number of runs.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) run_count
  integer ( kind = 4 ) test

  run_count = 0

  if ( n < 1 ) then
    return
  end if

  test = 0

  do i = 1, n

    if ( i == 1 .or. a(i) /= test ) then
      run_count = run_count + 1
      test = a(i)
    end if

  end do

  return
end
subroutine i4vec_search_binary_a ( n, a, b, indx )

!*****************************************************************************80
!
!! I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    Binary search is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.9,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the array to be searched.  A must
!    be sorted in ascending order.
!
!    Input, integer ( kind = 4 ) B, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    -1, B does not occur in A.
!    I, A(I) = B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid

  indx = - 1

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( a(mid) < b ) then
      low = mid + 1
    else if ( b < a(mid) ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine i4vec_search_binary_d ( n, a, b, indx )

!*****************************************************************************80
!
!! I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC for a value.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    Binary search is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.9,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the array to be searched.  A must
!    be sorted in descending order.
!
!    Input, integer ( kind = 4 ) B, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    -1, B does not occur in A.
!    I, A(I) = B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid

  indx = - 1

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( b < a(mid) ) then
      low = mid + 1
    else if ( a(mid) < b ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine i4vec_sort_bubble_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_BUBBLE_A ascending sorts an I4VEC using bubble sort.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  do i = 1, n - 1
    do j = i + 1, n
      if ( a(j) < a(i) ) then
        k    = a(i)
        a(i) = a(j)
        a(j) = k
      end if
    end do
  end do

  return
end
subroutine i4vec_sort_bubble_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_BUBBLE_D descending sorts an I4VEC using bubble sort.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  do i = 1, n - 1
    do j = i + 1, n
      if ( a(i) < a(j) ) then
        k    = a(i)
        a(i) = a(j)
        a(j) = k
      end if
    end do
  end do

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) t

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  t    = a(1)
  a(1) = a(n)
  a(n) = t
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    t    = a(1)
    a(1) = a(n1)
    a(n1) = t

  end do

  return
end
subroutine i4vec_sort_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_D descending sorts an I4VEC using heap sort.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) t

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into ascending heap form.
!
  call i4vec_heap_a ( n, a )
!
!  2: Sort A.
!
!  The smallest object in the heap is in A(1).
!  Move it to position A(N).
!
  t    = a(1)
  a(1) = a(n)
  a(n) = t
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_a ( n1, a )
!
!  Take the smallest object from A(1) and move it to A(N1).
!
    t     = a(1)
    a(1)  = a(n1)
    a(n1) = t

  end do

  return
end
subroutine i4vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call i4vec_permute ( n, indx, a )
!
!    after which A(1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) value

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      value = a(indxt)

    else

      indxt = indx(ir)
      value = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( value < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine i4vec_sort_heap_index_d ( n, a, indx )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call i4vec_permute ( n, indx, a )
!
!    after which A(1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) value

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      value = a(indxt)

    else

      indxt = indx(ir)
      value = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j+1)) < a(indx(j)) ) then
          j = j + 1
        end if
      end if

      if ( a(indx(j)) < value ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine i4vec_sort_insert_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_INSERT_A uses an ascending insertion sort on an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, A contains data to be sorted.
!    On output, the entries of A have been sorted in ascending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( a(j) <= x ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine i4vec_sort_insert_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_INSERT_D uses a descending insertion sort on an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, A contains data to be sorted.
!    On output, the entries of A have been sorted in descending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( x <= a(j) ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine i4vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_QUICK_A ascending sorts an I4VEC using quick sort.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    Input:
!
!      N = 7
!
!      A = (/ 6, 7, 3, 2, 9, 1, 8 /)
!
!    Output:
!
!      A = (/ 1, 2, 3, 6, 7, 8, 9 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 30
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( n <= 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call i4vec_part_quick_a ( n_segment, a(base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine i4vec_sort_shell_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_SHELL_A ascending sorts an I4VEC using Shell's sort.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an array to be sorted.
!    On output, the sorted array.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) asave
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) ipow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) maxpow

  if ( n <= 1 ) then
    return
  end if
!
!  Determine the smallest MAXPOW so that
!    N <= ( 3^MAXPOW - 1 ) / 2
!
  maxpow = 1

  do while ( 3**maxpow < 2 * n + 1 )
    maxpow = maxpow + 1
  end do

  if ( 1 < maxpow ) then
    maxpow = maxpow - 1
  end if
!
!  Now sort groups of size ( 3^IPOW - 1 ) / 2.
!
  do ipow = maxpow, 1, -1

    inc = ( 3**ipow - 1 ) / 2
!
!  Sort the values with indices equal to K mod INC.
!
    do k = 1, inc
!
!  Insertion sort of the items with index
!  INC+K, 2*INC+K, 3*INC+K, ...
!
      do i = inc + k, n, inc

        asave = a(i)
        ifree = i
        j = i - inc

        do

          if ( j < 1 ) then
            exit
          end if

          if ( a(j) <= asave ) then
            exit
          end if

          ifree = j
          a(j+inc) = a(j)
          j = j - inc

        end do

        a(ifree) = asave

      end do

    end do

  end do

  return
end
subroutine i4vec_sorted_undex ( x_num, x_val, x_unique_num, undx, xdnu )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNDEX returns unique sorted indexes for a sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of X, in sorted order,
!    and a vector XDNU, which identifies, for each entry of X, the index of
!    the unique sorted element of X.
!
!    This is all done with index vectors, so that the elements of
!    X are never moved.
!
!    Assuming X is already sorted, we examine the entries of X in order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the vector X could be
!    replaced by a compressed vector XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(I) = X(UNDX(I)).
!
!    We could then, if we wished, reconstruct the entire vector X, or
!    any element of it, by index, as follows:
!
!      X(I) = XU(XDNU(I)).
!
!    We could then replace X by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of X, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector X, the unique sort and
!    inverse unique sort vectors and the compressed unique sorted vector.
!
!      I    X    XU  Undx  Xdnu
!    ----+----+----+-----+-----+
!      1 | 11 |  11    1     1
!      2 | 11 |  22    5     1
!      3 | 11 |  33    8     1
!      4 | 11 |  55    9     1
!      5 | 22 |              2
!      6 | 22 |              2
!      7 | 22 |              2
!      8 | 33 |              3
!      9 | 55 |              4
!
!    UNDX(3) = 8 means that unique sorted item(3) is at X(8).
!    XDNU(6) = 2 means that X(6) is at unique sorted item(2).
!
!    XU(XDNU(I))) = X(I).
!    XU(I)        = X(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of data values.
!
!    Input, integer ( kind = 4 ) X_VAL(X_NUM), the data values.
!
!    Input, integer ( kind = 4 ) X_UNIQUE_NUM, the number of unique values i
!    n X_VAL.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Output, integer ( kind = 4 ) UNDX(X_UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(X_NUM), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) x_num
  integer ( kind = 4 ) x_unique_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) undx(x_unique_num)
  integer ( kind = 4 ) x_val(x_num)
  integer ( kind = 4 ) xdnu(x_num)
!
!  Walk through the sorted array.
!
  i = 1

  j = 1
  undx(j) = i

  xdnu(i) = j

  do i = 2, x_num

    if ( x_val(i) /= x_val(undx(j)) ) then
      j = j + 1
      undx(j) = i
    end if

    xdnu(i) = j

  end do

  return
end
subroutine i4vec_sorted_unique ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in A.
!
!    Input/output, integer ( kind = 4 ) A(N).  On input, the sorted
!    integer array.  On output, the unique elements in A.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a(itest) /= a(unique_num) ) then
      unique_num = unique_num + 1
      a(unique_num) = a(itest)
    end if

  end do

  return
end
subroutine i4vec_sorted_unique_count ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE_COUNT counts the unique elements in a sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    Because the array is sorted, this algorithm is O(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the sorted array to examine.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num

  if ( n < 1 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do i = 2, n

    if ( a(i-1) /= a(i) ) then
      unique_num = unique_num + 1
    end if

  end do

  return
end
subroutine i4vec_sorted_unique_hist ( n, a, maxuniq, unique_num, auniq, acount )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE_HIST histograms the unique elements of a sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the array to examine.  The elements of A
!    should have been sorted.
!
!    Input, integer ( kind = 4 ) MAXUNIQ, the maximum number of unique elements
!    that can be handled.  If there are more than MAXUNIQ unique
!    elements in A, the excess will be ignored.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements.
!
!    Output, integer ( kind = 4 ) AUNIQ(UNIQUE_NUM), the unique elements of A.
!
!    Output, integer ( kind = 4 ) ACOUNT(UNIQUE_NUM), the number of times
!    each element of AUNIQ occurs in A.
!
  implicit none

  integer ( kind = 4 ) maxuniq
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) acount(maxuniq)
  integer ( kind = 4 ) auniq(maxuniq)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
!
!  Start taking statistics.
!
  unique_num = 0

  do i = 1, n

    if ( i == 1 ) then

      unique_num = 1
      auniq(unique_num) = a(1)
      acount(unique_num) = 1

    else if ( a(i) == auniq(unique_num) ) then

      acount(unique_num) = acount(unique_num) + 1

    else if ( unique_num < maxuniq ) then

      unique_num = unique_num + 1
      auniq(unique_num) = a(i)
      acount(unique_num) = 1

    end if

  end do

  return
end
subroutine i4vec_split ( n, a, split, split_index )

!*****************************************************************************80
!
!! I4VEC_SPLIT "splits" an unsorted I4VEC based on a splitting value.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    If the vector is already sorted, it is simpler to do a binary search
!    on the data than to call this routine.
!
!    The vector is not assumed to be sorted before input, and is not
!    sorted during processing.  If sorting is not needed, then it is
!    more efficient to use this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to split.  On output,
!    all the entries of A that are less than or equal to SPLIT
!    are in A(1:SPLIT_INDEX).
!
!    Input, integer ( kind = 4 ) SPLIT, the value used to split the vector.
!    It is not necessary that any value of A actually equal SPLIT.
!
!    Output, integer ( kind = 4 ) SPLIT_INDEX, indicates the position of the
!    last entry of the split vector that is less than or equal to SPLIT.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) split
  integer ( kind = 4 ) split_index
  integer ( kind = 4 ) t
!
!  Partition the vector into A1, A2, A3, where
!    A1 = A(I1:J1) holds values <= SPLIT,
!    A2 = A(I2:J2) holds untested values,
!    A3 = A(I3:J3) holds values > SPLIT.
!
  i1 = 1
  j1 = 0

  i2 = 1
  j2 = n

  i3 = n + 1
  j3 = n
!
!  Pick the next item from A2, and move it into A1 or A3.
!  Adjust indices appropriately.
!
  do i = 1, n

    if ( a(i2) <= split ) then

      i2 = i2 + 1
      j1 = j1 + 1

    else

      t       = a(i2)
      a(i2)   = a(i3-1)
      a(i3-1) = t

      i3 = i3 - 1
      j2 = j2 - 1

    end if

  end do

  split_index = j1

  return
end
subroutine i4vec_std ( n, a, std )

!*****************************************************************************80
!
!! I4VEC_STD returns the standard deviation of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The standard deviation of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      std ( X(1:n) ) = sqrt ( sum ( ( X(1:n) - mean )^2 ) / ( n - 1 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector whose variance is desired.
!
!    Output, real ( kind = 8 ) STD, the standard deviation of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  real ( kind = 8 ) mean
  real ( kind = 8 ) std

  if ( n < 2 ) then

    std = 0.0D+00

  else

    mean = real ( sum ( a(1:n) ), kind = 8 ) / real ( n, kind = 8 )

    std = sum ( ( real ( a(1:n), kind = 8 ) - mean )**2 )

    std = sqrt ( std / real ( n - 1, kind = 8 ) )

  end if

  return
end
function i4vec_sum ( n, a )

!*****************************************************************************80
!
!! I4VEC_SUM returns the sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In FORTRAN90, this facility is offered by the built in
!    SUM function:
!
!      I4VEC_SUM ( N, A ) = SUM ( A(1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) I4VEC_SUM, the sum of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4vec_sum

  i4vec_sum = sum ( a(1:n) )

  return
end
subroutine i4vec_swap ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_SWAP swaps the entries of two I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the arrays.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the vectors to swap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) a3(n)

  a3(1:n) = a1(1:n)
  a1(1:n) = a2(1:n)
  a2(1:n) = a3(1:n)

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!
!        1    2    3    4    5
!        6    7    8    9   10
!       11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5i12)' ) a(ilo:ihi)
  end do

  return
end
subroutine i4vec_undex ( x_num, x_val, x_unique_num, undx, xdnu )

!*****************************************************************************80
!
!! I4VEC_UNDEX returns unique sorted indexes for an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of X, in sorted order,
!    and a vector XDNU, which identifies, for each entry of X, the index of
!    the unique sorted element of X.
!
!    This is all done with index vectors, so that the elements of
!    X are never moved.
!
!    The first step of the algorithm requires the indexed sorting
!    of X, which creates arrays INDX and XDNI.  (If all the entries
!    of X are unique, then these arrays are the same as UNDX and XDNU.)
!
!    We then use INDX to examine the entries of X in sorted order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the vector X could be
!    replaced by a compressed vector XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(1:X_UNIQUE_NUM) = X(UNDX(1:X_UNIQUE_NUM)).
!
!    We could then, if we wished, reconstruct the entire vector X, or
!    any element of it, by index, as follows:
!
!      X(I) = XU(XDNU(I)).
!
!    We could then replace X by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of X, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector X, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I    X  Indx  Xdni      XU  Undx  Xdnu
!    ----+----+-----+-----+-------+-----+-----+
!      1 | 11     1     1 |    11     1     1
!      2 | 22     3     5 |    22     2     2
!      3 | 11     6     2 |    33     4     1
!      4 | 33     9     8 |    55     5     3
!      5 | 55     2     9 |                 4
!      6 | 11     7     3 |                 1
!      7 | 22     8     6 |                 2
!      8 | 22     4     7 |                 2
!      9 | 11     5     4 |                 1
!
!    INDX(2) = 3 means that sorted item(2) is X(3).
!    XDNI(2) = 5 means that X(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
!    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = X(I).
!    XU(I)        = X(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of data values.
!
!    Input, integer ( kind = 4 ) X_VAL(X_NUM), the data values.
!
!    Input, integer ( kind = 4 ) X_UNIQUE_NUM, the number of unique values
!    in X_VAL.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Output, integer ( kind = 4 ) UNDX(X_UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(X_NUM), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) x_num
  integer ( kind = 4 ) x_unique_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(x_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) undx(x_unique_num)
  integer ( kind = 4 ) x_val(x_num)
  integer ( kind = 4 ) xdnu(x_num)
!
!  Implicitly sort the array.
!
  call i4vec_sort_heap_index_a ( x_num, x_val, indx )
!
!  Walk through the implicitly sorted array.
!
  i = 1

  j = 1
  undx(j) = indx(i)

  xdnu(indx(i)) = j

  do i = 2, x_num

    if ( x_val(indx(i)) /= x_val(undx(j)) ) then
      j = j + 1
      undx(j) = indx(i)
    end if

    xdnu(indx(i)) = j

  end do

  return
end
subroutine i4vec_uniform_ab ( n, a, b, seed, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM_AB returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(N), a vector of numbers between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
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

    x(i) = value

  end do

  return
end
subroutine i4vec_unique_count ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_UNIQUE_COUNT counts the unique elements in an unsorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    Because the array is unsorted, this algorithm is O(N^2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the unsorted array to examine.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
!    of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) unique_num

  unique_num = 0

  do i = 1, n

    unique_num = unique_num + 1

    do j = 1, i - 1

      if ( a(i) == a(j) ) then
        unique_num = unique_num - 1
        exit
      end if

    end do

  end do

  return
end
subroutine i4vec_unique_index ( n, a, unique_index )

!*****************************************************************************80
!
!! I4VEC_UNIQUE_INDEX indexes the unique occurrence of values in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    For element A(I) of the vector, UNIQUE_INDEX(I) is the uniqueness index
!    of A(I).  That is, if A_UNIQUE contains the unique elements of A,
!    gathered in order, then
!
!      A_UNIQUE ( UNIQUE_INDEX(I) ) = A(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) UNIQUE_INDEX(N), the unique index.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) unique_index(n)
  integer ( kind = 4 ) unique_num

  unique_index(1:n) = -1
  unique_num = 0

  do i = 1, n

    if ( unique_index(i) == -1 ) then

      unique_num = unique_num + 1
      unique_index(i) = unique_num

      do j = i + 1, n
        if ( a(i) == a(j) ) then
          unique_index(j) = unique_num
        end if
      end do

    end if

  end do

  return
end
subroutine i4vec_value_index ( n, a, value, max_index, n_index, value_index )

!*****************************************************************************80
!
!! I4VEC_VALUE_INDEX indexes entries equal to a given value in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    Input:
!
!      N = 10
!      A = (  2, 3, 1, 3, 2, 4, 2, 3, 5, 3 )
!      X_VALUE = 3
!
!    Output:
!
!      N_INDEX = 4
!      VALUE_INDEX = ( 2, 4, 8, 10 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) A(N), the array to be indexed.
!
!    Input, integer ( kind = 4 ) VALUE, a value to be searched for.
!
!    Input, integer ( kind = 4 ) MAX_INDEX, the maximum number of indices
!    to find.
!
!    Output, integer ( kind = 4 ) N_INDEX, the number of entries equal to VALUE.
!
!    Output, integer ( kind = 4 ) VALUE_INDEX(MAX_INDEX), the indices of entries
!    equal to VALUE.
!
  implicit none

  integer ( kind = 4 ) max_index
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n_index
  integer ( kind = 4 ) value
  integer ( kind = 4 ) value_index(max_index)

  n_index = 0

  do i = 1, n

    if ( a(i) == value ) then

      if ( max_index <= n_index ) then
        return
      end if

      n_index = n_index + 1
      value_index(n_index) = i

    end if

  end do

  return
end
subroutine i4vec_value_num ( n, a, value, value_num )

!*****************************************************************************80
!
!! I4VEC_VALUE_NUM counts entries equal to a given value in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) A(N), the array to be indexed.
!
!    Input, integer ( kind = 4 ) VALUE, a value to be searched for.
!
!    Input, integer ( kind = 4 ) VALUE_NUM, the number of times the
!    value occurs.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) value
  integer ( kind = 4 ) value_num

  value_num = 0

  do i = 1, n

    if ( a(i) == value ) then
      value_num = value_num + 1
    end if

  end do

  return
end
subroutine i4vec_variance ( n, a, variance )

!*****************************************************************************80
!
!! I4VEC_VARIANCE returns the variance of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector whose variance is desired.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance

  if ( n < 2 ) then

    variance = 0.0D+00

  else

    mean = real ( sum ( a(1:n) ), kind = 8 ) / real ( n, kind = 8 )

    variance = sum ( ( real ( a(1:n), kind = 8 ) - mean )**2 )

    variance = variance / real ( n - 1, kind = 8 )

  end if

  return
end
function i4vec_width ( n, a )

!*****************************************************************************80
!
!! I4VEC_WIDTH returns the "width" of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The width of an integer vector is simply the maximum of the widths of
!    its entries.
!
!    The width of a single integer is the number of characters
!    necessary to print it.
!
!    The width of an integer vector can be useful when the vector is
!    to be printed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, integer ( kind = 4 ) I4VEC_WIDTH, the width of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_width
  integer ( kind = 4 ) i4vec_width

  i4vec_width = -1

  do i = 1, n
    i4vec_width = max ( i4vec_width, i4_width ( a(i) ) )
  end do

  return
end
subroutine i4vec_zero ( n, a )

!*****************************************************************************80
!
!! I4VEC_ZERO sets the entries of an I4VEC to 0.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Output, integer ( kind = 4 ) A(N), the vector, which has been set to zero.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)

  a(1:n) = 0

  return
end
subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares entries of an I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), contain the two components
!    of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item J < item I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end
subroutine i4vec2_print ( n, a, b, title )

!*****************************************************************************80
!
!! I4VEC2_PRINT prints a pair of integer vectors.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), B(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,1x,i10,2x,i10)' ) i, ':', a(i), b(i)
  end do

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      t     = a1(i)
      a1(i) = a1(j)
      a1(j) = t

      t     = a2(i)
      a2(i) = a2(j)
      a2(j) = t
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_sort_d ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_D descending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      t     = a1(i)
      a1(i) = a1(j)
      a1(j) = t

      t     = a2(i)
      a2(i) = a2(j)
      a2(j) = t
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )
      isgn = -isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
function l_to_i4 ( l )

!*****************************************************************************80
!
!! L_TO_I4 converts an L to an I4.
!
!  Discussion:
!
!    0 is FALSE, and anything else if TRUE.
!
!    An I4 is an integer value.
!    An L is a logical value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical L, a logical value.
!
!    Output, integer ( kind = 4 ) L_TO_I4, the integer value of L.
!
  implicit none

  logical l
  integer ( kind = 4 ) l_to_i4
  integer ( kind = 4 ) value

  if ( l ) then
    value = 1
  else
    value = 0
  end if

  l_to_i4 = value

  return
end
subroutine i4vec2_sorted_unique ( n, a1, a2, unique_num )

!*****************************************************************************80
!
!! I4VEC2_SORTED_UNIQUE keeps the unique elements in a sorted I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of unique items.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a1(itest) /= a1(unique_num) .or. a2(itest) /= a2(unique_num) ) then

      unique_num = unique_num + 1

      a1(unique_num) = a1(itest)
      a2(unique_num) = a2(itest)

    end if

  end do

  return
end
subroutine perm_check ( n, p, base, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from BASE to
!    to BASE+N-1 occurs among the N entries of the permutation.
!
!    Set the input quantity BASE to 0, if P is a 0-based permutation,
!    or to 1 if P is a 1-based permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Input, integer ( kind = 4 ) BASE, the index base.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) find
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seek

  ierror = 0

  do seek = base, base + n - 1

    ierror = 1

    do find = 1, n
      if ( p(find) == seek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.'
      stop
    end if

  end do

  return
end
subroutine perm_cycle ( n, iopt, p, isgn, ncycle )

!*****************************************************************************80
!
!! PERM_CYCLE analyzes a permutation.
!
!  Discussion:
!
!    The routine will count cycles, find the sign of a permutation,
!    and tag a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      IOPT = 1
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input, integer ( kind = 4 ) IOPT, requests tagging.
!    0, the permutation will not be tagged.
!    1, the permutation will be tagged.
!
!    Input/output, integer ( kind = 4 ) P(N).  On input, P describes a
!    permutation, in the sense that entry I is to be moved to P(I).
!    If IOPT = 0, then P will not be changed by this routine.
!    If IOPT = 1, then on output, P will be "tagged".  That is,
!    one element of every cycle in P will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of P
!    which is negative, moving to I2 = ABS(P(I1)), then to
!    P(I2), and so on, until returning to I1.
!
!    Output, integer ( kind = 4 ) ISGN, the "sign" of the permutation, which is
!    +1 if the permutation is even, -1 if odd.  Every permutation
!    may be produced by a certain number of pairwise switches.
!    If the number of switches is even, the permutation itself is
!    called even.
!
!    Output, integer ( kind = 4 ) NCYCLE, the number of cycles in the
!    permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) is
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ncycle
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_CYCLE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if

  is = 1
  ncycle = n

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      ncycle = ncycle - 1
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = - i4_sign ( p(i) )
    end if

    p(i) = is * abs ( p(i) )

  end do

  isgn = 1 - 2 * mod ( n - ncycle, 2 )

  return
end
subroutine perm_uniform ( n, base, seed, p )

!*****************************************************************************80
!
!! PERM_UNIFORM selects a random permutation of N objects.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input, integer ( kind = 4 ) BASE, is 0 for a 0-based permutation and 1 for
!    a 1-based permutation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) P(N), the permutation.  P(I) is the "new"
!    location of the object originally at I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  do i = 1, n
    p(i) = ( i - 1 ) + base
  end do

  do i = 1, n
    j = i4_uniform_ab ( i, n, seed )
    k    = p(i)
    p(i) = p(j)
    p(j) = k
  end do

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
!    US Department of Commerce, 1964, pages 870-873.
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
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

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
