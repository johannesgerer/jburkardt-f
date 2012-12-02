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
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

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
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
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
  integer ( kind = 4 ) i4_huge

  if ( i == 0 ) then

    i4_log_2 = - i4_huge ( )

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
function i4_seed_advance ( seed )

!*****************************************************************************80
!
!! I4_SEED_ADVANCE "advances" the seed.
!
!  Discussion:
!
!    This routine implements one step of the recursion
!
!      SEED = 16807 * SEED mod ( 2**31 - 1 )
!
!    This version of the routine does not check whether the input value of
!    SEED is zero.  If the input value is zero, the output value will be zero.
!
!    If we repeatedly use the output of SEED_ADVANCE as the next input,
!    and we start with SEED = 12345, then the first few iterates are:
!
!         Input      Output
!          SEED        SEED
!
!         12345   207482415
!     207482415  1790989824
!    1790989824  2035175616
!    2035175616    77048696
!      77048696    24794531
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
!    Input, integer ( kind = 4 ) SEED, the seed value.
!
!    Output, integer ( kind = 4 ) I4_SEED_ADVANCE, the "next" seed.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) i4_seed_advance
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_new

  k = seed / 127773

  seed_new = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed_new < 0 ) then
    seed_new = seed_new + i4_huge ( )
  end if

  i4_seed_advance = seed_new

  return
end
function r4_ieee_uniform ( seed )

!*****************************************************************************80
!
!! R4_IEEE_UNIFORM computes an "IEEE uniform" pseudorandom real number.
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
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 4 ) R4_IEEE_RANDOM, a pseudorandom value,
!    uniformly distributed in the IEEE distribution.
!
  implicit none

  integer ( kind = 4 ) e
  integer ( kind = 4 ) e2
  integer ( kind = 4 ) f
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_seed_advance
  integer ( kind = 4 ) len
  integer ( kind = 4 ) pos
  real ( kind = 4 ) r4
  real ( kind = 4 ) r4_ieee_uniform
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) seed
!
!  Use bits of SEED to define the value.
!
  pos = 0
  len = 1
  s = ibits ( seed, pos, len )

  if ( s == 0 ) then
    s2 = 1
  else
    s2 = -1
  end if

  pos = 1
  len = 8
  e = ibits ( seed, pos, len )

  e2 = e - 127

  pos = 9
  len = 23
  f = ibits ( seed, pos, len )

! write ( *, '(2x,i2,2x,i8,2x,i12)' ) s, e, f

  f2 = f + 2**23
! write ( *, '(2x,i2,2x,i8,2x,i12)' ) s2, e2, f2

  r4 = real ( f2, kind = 4 ) / 2.0E+00**23
! write ( *, '(2x,g14.6)' ) r4

  if ( .true. )then

    r4 = s2 * real ( f2, kind = 4 ) * 2.0E+00 ** ( e2 - 23 )

  else

    r4 = s2 * real ( f2, kind = 4 ) 

    if ( 0 < e2 - 23 ) then
      do i = 1, e2 - 23
        r4 = r4 * 2.0E+00
      end do
    else if ( e2 - 23 < 0 ) then
      do i = 1, 23 - e2
        r4 = r4 / 2.0E+00
      end do
    end if

  end if
!
!  Advance the seed.
!
  seed = i4_seed_advance ( seed )

  r4_ieee_uniform = r4

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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
