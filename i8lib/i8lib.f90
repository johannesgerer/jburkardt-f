function i8_bit_hi1 ( n )

!*****************************************************************************80
!
!! I8_BIT_HI1 returns the position of the high 1 bit base 2 in an I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
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
!    28 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) N, the integer to be measured.
!    N should be nonnegative.  If N is nonpositive, I8_BIT_HI1
!    will always be 0.
!
!    Output, integer ( kind = 8 ) I8_BIT_HI1, the number of bits base 2.
!
  implicit none

  integer ( kind = 8 ) :: bit
  integer ( kind = 8 ) :: i8_bit_hi1
  integer ( kind = 8 ) :: i
  integer ( kind = 8 ) :: n

  i = n
  bit = 0

  do

    if ( i <= 0 ) then
      exit
    end if

    bit = bit + 1
    i = i / 2

  end do

  i8_bit_hi1 = bit

  return
end
function i8_bit_lo0 ( n )

!*****************************************************************************80
!
!! I8_BIT_LO0 returns the position of the low 0 bit base 2 in an I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
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
!    28 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) N, the integer to be measured.
!    N should be nonnegative.
!
!    Output, integer ( kind = 8 ) I8_BIT_LO0, the position of the low 1 bit.
!
  implicit none

  integer ( kind = 8 ) bit
  integer ( kind = 8 ) i
  integer ( kind = 8 ) i2
  integer ( kind = 8 ) i8_bit_lo0
  integer ( kind = 8 ) n

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

  i8_bit_lo0 = bit

  return
end
function i8_choose ( n, k )

!*****************************************************************************80
!
!! I8_CHOOSE computes the binomial coefficient C(N,K) as an I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
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
!    27 June 2010
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
!    Input, integer ( kind = 8 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 8 ) I8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 8 ) i
  integer ( kind = 8 ) i8_choose
  integer ( kind = 8 ) k
  integer ( kind = 8 ) mn
  integer ( kind = 8 ) mx
  integer ( kind = 8 ) n
  integer ( kind = 8 ) value

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

  i8_choose = value

  return
end
function i8_huge ( )

!*****************************************************************************80
!
!! I8_HUGE returns a "huge" I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
!
!    On an IEEE machine, I8_HUGE should be 2**63 - 1, and its
!    bit pattern should be
!
!     0111111111111111111111111111111111111111111111111111111111111111
!
!    In this case, its numerical value is 9223372036854775807.
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
!    Output, integer ( kind = 8 ) I8_HUGE, a "huge" I8.
!
  implicit none

  integer ( kind = 8 ) i8_huge

  i8_huge = 9223372036854775807_8

  return
end
function i8_huge_normalizer ( )

!*****************************************************************************80
!
!! I8_HUGE_NORMALIZER returns the "normalizer" for I8_HUGE.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
!
!    The value returned is 1 / ( I8_HUGE + 1 ).
!
!    For any I8, it should be the case that
!
!     -1 < I8 * I8_HUGE_NORMALIZER < 1.
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
!    Output, real ( kind = 8 ) I8_HUGE_NORMALIZER, the "normalizer"
!    for I8_HUGE.
!
  implicit none

  real    ( kind = 8 ) i8_huge_normalizer

  i8_huge_normalizer = 1.084202172485504434007D-19

  return
end
function i8_power ( i, j )

!*****************************************************************************80
!
!! I8_POWER returns the integer power of an I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) I, J, the base and the power.  
!    J should be nonnegative.
!
!    Output, integer ( kind = 8 ) I8_POWER, the value of I^J.
!
  implicit none

  integer ( kind = 8 ) i
  integer ( kind = 8 ) i8_power
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k

  if ( j < 0 ) then

    if ( i == 1 ) then
      i8_power = 1
    else if ( i == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I8_POWER - Fatal error!'
      write ( *, '(a)' ) '  I^J requested, with I = 0 and J negative.'
      stop
    else
      i8_power = 0
    end if

  else if ( j == 0 ) then

    if ( i == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I8_POWER - Fatal error!'
      write ( *, '(a)' ) '  I^J requested, with I = 0 and J = 0.'
      stop
    else
      i8_power = 1
    end if

  else if ( j == 1 ) then

    i8_power = i

  else

    i8_power = 1
    do k = 1, j
      i8_power = i8_power * i
    end do

  end if

  return
end
function i8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I8_UNIFORM returns a scaled pseudorandom I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
!
!    Note that ALL integer variables in this routine are
!    of type integer ( kind = 8 )!
!
!    The input arguments to this function should NOT be constants; they should
!    be variables of type integer ( kind = 8 )!
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
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
!    Input, integer ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 8 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 8 ) I8_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 8 ) a
  integer ( kind = 8 ) b
  integer ( kind = 8 ) i8_uniform
  real    ( kind = 8 ) r
  real    ( kind = 8 ) r8i8_uniform_01
  integer ( kind = 8 ) seed
  integer ( kind = 8 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  r = r8i8_uniform_01 ( seed )
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0D+00 - r ) * ( real ( min ( a, b ), kind = 8 ) - 0.5D+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 8 ) + 0.5D+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 8 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i8_uniform = value

  return
end
subroutine i8_uniform2 ( seed )

!*****************************************************************************80
!
!! I8_UNIFORM2 returns a pseudorandom I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 8 ) SEED, the pseudorandom integer.
!
  implicit none

  integer ( kind = 8 ), parameter :: c = 3037000493_8
  integer ( kind = 8 ), parameter :: m = 2862933555777941757_8
  integer ( kind = 8 ) seed

  seed = m * seed + c

  return
end
subroutine i8_uniform3 ( seed )

!*****************************************************************************80
!
!! I8_UNIFORM3 returns a pseudorandom I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 8 ) SEED, the pseudorandom integer.
!
  implicit none

  integer ( kind = 8 ), parameter :: c = 3037000493_8
  integer ( kind = 8 ), parameter :: i8_huge = 9223372036854775807_8 
  integer ( kind = 8 ), parameter :: m = 2862933555777941757_8
  integer ( kind = 8 ) seed

  seed = m * seed + c

  if ( seed < 0 ) then
    seed = seed + i8_huge
  end if

  return
end
function i8_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I8_UNIFORM_AB returns a pseudorandom I8 between two limits.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 8 ) SEED, the seed for the
!    random number generator.
!
!    Output, integer ( kind = 8 ) I8_UNIFORM_AB, the pseudorandom
!    value, between A and B.
!
  implicit none

  integer ( kind = 8 ) a
  integer ( kind = 8 ) b
  integer ( kind = 8 ), parameter :: c = 3037000493_8
  integer ( kind = 8 ), parameter :: i8_huge = 9223372036854775807_8 
  integer ( kind = 8 ) i8_uniform_ab
  integer ( kind = 8 ), parameter :: m = 2862933555777941757_8

  integer ( kind = 8 ) seed
  integer ( kind = 8 ) value

  seed = m * seed + c

  if ( seed < 0 ) then
    seed = seed + i8_huge
  end if

  i8_uniform_ab = a + mod ( seed, b - a + 1 )

  return
end
function i8_xor ( i, j )

!*****************************************************************************80
!
!! I8_XOR calculates the exclusive OR of two I8's.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
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
!    Input, integer ( kind = 8 ) I, J, two values whose exclusive OR is needed.
!
!    Output, integer ( kind = 8 ) I8_XOR, the exclusive OR of I and J.
!
  implicit none

  integer ( kind = 8 ) i
  integer ( kind = 8 ) i1
  integer ( kind = 8 ) i2
  integer ( kind = 8 ) i8_xor
  integer ( kind = 8 ) j
  integer ( kind = 8 ) j1
  integer ( kind = 8 ) j2
  integer ( kind = 8 ) k
  integer ( kind = 8 ) l

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

  i8_xor = k

  return
end
function r8_uniform2 ( seed )

!*****************************************************************************80
!
!! I8_UNIFORM2 returns a pseudorandom I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 8 ) SEED, the pseudorandom integer.
!
  implicit none

  integer ( kind = 8 ), parameter :: a = 2862933555777941757_8
  integer ( kind = 8 ), parameter :: b = 3037000493_8
! integer ( kind = 8 ), parameter :: i4_t32 = 2147483648_8
  integer ( kind = 8 ), parameter :: i4_t32 = 4294967296_8
  real    ( kind = 8 ) r8_uniform2
  integer ( kind = 8 ) s1
  integer ( kind = 8 ) s2
  integer ( kind = 8 ) seed
  real    ( kind = 8 ), parameter :: t32 = 4.656612873077392578125E-10_8

  seed = a * seed + b

  s2 = mod ( seed, i4_t32 )
  s1 = seed / i4_t32

  r8_uniform2 = ( real ( s1, kind = 8 ) + real ( s2, kind = 8 ) * t32 ) * t32

  return
end
function r8_uniform3 ( seed )

!*****************************************************************************80
!
!! I8_UNIFORM3 returns a pseudorandom I8.
!
!  Discussion:
!
!    An I8 is an integer ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 8 ) SEED, the pseudorandom integer.
!
  implicit none

  integer ( kind = 8 ), parameter :: a = 2862933555777941757_8
  integer ( kind = 8 ), parameter :: b = 3037000493_8
! integer ( kind = 8 ), parameter :: i4_t32 = 2147483648_8
  integer ( kind = 8 ), parameter :: i4_t32 = 4294967296_8
  integer ( kind = 8 ), parameter :: i8_huge = 9223372036854775807_8 
  real    ( kind = 8 ) r8_uniform3
  integer ( kind = 8 ) s1
  integer ( kind = 8 ) s2
  integer ( kind = 8 ) seed
  real    ( kind = 8 ), parameter :: t32 = 4.656612873077392578125E-10_8

  seed = a * seed + b

  if ( seed < 0 ) then
    seed = seed + i8_huge
  end if

  s2 = mod ( seed, i4_t32 )
  s1 = seed / i4_t32

  r8_uniform3 = ( real ( s1, kind = 8 ) + real ( s2, kind = 8 ) * t32 ) * t32

  return
end
function r8i8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8I8_UNIFORM returns a scaled pseudorandom R8 using an I8 seed.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    An I8 is an integer ( kind = 8 ) value.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 8 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8I8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  integer ( kind = 8 ) k
  real    ( kind = 8 ) r8i8_uniform
  integer ( kind = 8 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8I8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + huge ( seed )
  end if

  r8i8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8i8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8I8_UNIFORM_01 returns a unit pseudorandom R8 using an I8 seed.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    An I8 is an integer ( kind = 8 ) value.
!
!    This routine implements the recursion
!
!      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8I8_UNIFORM_01
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
!    31 May 2007
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
!    Input/output, integer ( kind = 8 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8I8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 8 ) k
  real    ( kind = 8 ) r8i8_uniform_01
  integer ( kind = 8 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8I8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + huge ( seed )
  end if

  r8i8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
