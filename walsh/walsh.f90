subroutine ffwt ( n, x )

!*****************************************************************************80
!
!! FFWT performs an in-place fast Walsh transform.
!
!  Discussion:
!
!    This routine performs a fast Walsh transform on an input series X
!    leaving the transformed results in X. 
!    X is dimensioned N, which must be a power of 2.
!    The results of this Walsh transform are in sequency order.
!
!    The output sequence could be normalized by dividing by N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = 8 ) X(N), the data to be transformed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) hold
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) js
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mw
  integer ( kind = 4 ) mw1
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) nz2
  integer ( kind = 4 ) nzi
  integer ( kind = 4 ) nzn
  integer ( kind = 4 ) two_power(24)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) z

  m = i4_log_2 ( n )

  do i = 1, m
    two_power(i) = 2**( m - i )
  end do

  do l = 1, m

    nz = 2**( l - 1 )
    nzi = 2 * nz
    nzn = n / nzi
    nz2 = nz / 2
    if ( nz2 == 0 ) then
      nz2 = 1
    end if

    do i = 1, nzn

      js = ( i - 1 ) * nzi
      z = 1.0D+00

      do ii = 1, 2

        do j = 1, nz2
          js = js + 1
          j2 = js + nz
          hold = x(js) + z * x(j2)
          z = - z
          x(j2) = x(js) + z * x(j2)
          x(js) = hold
          z = - z
        end do

        if ( l == 1 ) then
          exit
        end if

        z = - 1.0D+00

      end do

    end do
  end do
!
!  Bit reversal section.
!
  nw = 0
  do k = 1, n
!
!  Choose correct index and switch elements if not already switched.
!
    if ( k < nw + 1 ) then
      hold = x(nw+1)
      x(nw+1) = x(k)
      x(k) = hold
    end if
!
!  Bump up series by 1.
!
    do i = 1, m

      ii = i
      if ( nw < two_power(i) ) then
        exit
      end if
      mw = nw / two_power(i)
      mw1 = mw / 2
      if ( mw <= 2 * mw1 ) then
        exit
      end if

      nw = nw - two_power(i)

    end do

    nw = nw + two_power(ii)

  end do

  return
end
subroutine fwt ( n, x, y )

!*****************************************************************************80
!
!! FWT performs a fast Walsh transform.
!
!  Discussion:
!
!    This routine performs a fast Walsh transform on an input series X
!    leaving the transformed results in X. 
!    X is dimensioned N, which must be a power of 2.
!    The results of this Walsh transform are in sequency order.
!
!    The output sequence could be normalized by dividing by N.
!
!    Note that the program text in the reference included the line
!      y(jd) = abs ( x(j) - x(j2) )
!    which has been corrected to:
!      y(jd) = x(j) - x(j2)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = 8 ) X(N), the data to be transformed.
!
!    Workspace, real Y(N).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jd
  integer ( kind = 4 ) js
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) nzi
  integer ( kind = 4 ) nzn
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  n2 = n / 2
  m = i4_log_2 ( n )

  do l = 1, m

    ny = 0
    nz = 2**(l-1)
    nzi = 2 * nz
    nzn = n / nzi

    do i = 1, nzn

      nx = ny + 1
      ny = ny + nz
      js = ( i - 1 ) * nzi
      jd = js + nzi + 1

      do j = nx, ny
        js = js + 1
        j2 = j + n2
        y(js) = x(j) + x(j2)
        jd = jd - 1
        y(jd) = x(j) - x(j2)
      end do

    end do

    x(1:n) = y(1:n)

  end do

  return
end
subroutine haar ( n, x, y )

!*****************************************************************************80
!
!! HAAR performs a Haar transform.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = 8 ) X(N), the data to be transformed.
!
!    Workspace, real Y(N).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) l3
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  k = i4_log_2 ( n )

  do i = 1, k

    l = k + 1 - i
    l2 = 2**( l - 1 )

    y(1:2*l2) = x(1:2*l2)

    do j = 1, l2
       l3 = l2 + j
       jj = 2 * j - 1
       x(j) = y(jj) + y(jj+1)
       x(l3) = y(jj) - y(jj+1)
    end do

  end do

  return
end
subroutine haarin ( n, x, y )

!*****************************************************************************80
!
!! HAARIN inverts a Haar transform.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = 8 ) X(N), the data to be transformed.
!
!    Workspace, real Y(N).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jj1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lj
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  k = i4_log_2 ( n )

  do i = 1, k

    l = 2**( i - 1 )

    y(1:2*l) = x(1:2*l)

    do j = 1, l
      lj = l + j
      jj = 2 * j
      jj1 = jj - 1
      x(jj) = y(j) - y(lj)
      x(jj1) = y(j) + y(lj)
    end do

  end do

  return
end
subroutine hnorm ( n, x )

!*****************************************************************************80
!
!! HNORM computes normalization factors for a forward or inverse Haar transform.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = 8 ) X(N), the data to be transformed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) k
  real ( kind = 8 ) wlk
  real ( kind = 8 ) x(n)

  k = i4_log_2 ( n )

  x(1) = x(1) / 2.0D+00**k

  if ( 1 <= k ) then
    x(2) = x(2) / 2.0D+00**k
  end if

  do ii = 2, k

    i = ii - 1
    wlk = 1.0D+00 / 2.0D+00**( k - i )
    jmin = 2**i + 1
    jmax = 2**ii

    x(jmin:jmax) = x(jmin:jmax) * wlk

  end do

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
subroutine r8vec_shift_circular ( shift, n, x )

!*****************************************************************************80
!
!! R8VEC_SHIFT_CIRCULAR performs a circular shift on an R8VEC.
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
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SHIFT, the amount by which each entry is to
!    be shifted.
!
!    Input, integer ( kind = 4 ) N, the length of the vector.
!
!    Input/output, real ( kind = 8 ) X(N), the vector to be shifted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) shift
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = x(1:n)

  do i = 1, n
    j = i4_wrap ( i - shift, 1, n )
    x(i) = y(j)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

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
subroutine walsh ( n, x, y )

!*****************************************************************************80
!
!! WALSH performs a fast Walsh transform.
!
!  Discussion:
!
!    This routine performs a fast Wash transform on an input series X
!    leaving the transformed results in X.  The array Y is working space.
!    X and Y are dimensioned N, which must be a power of 2.
!    The results of this Walsh transform are in sequency order.
!
!    The output sequence could be normalized by dividing by N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = 8 ) X(N), the data to be transformed.
!
!    Workspace, real Y(N).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n/2)
  real ( kind = 8 ) z

  n2 = n / 2
  m = i4_log_2 ( n )
  z = - 1.0D+00

  do j = 1, m

    n1 = 2**( m - j + 1 )
    j1 = 2**( j - 1 )

    do l = 1, j1

      is = ( l - 1 ) * n1 + 1
      i1 = 0
      w = z

      do i = is, is + n1 - 1, 2
        a = x(i)
        x(is+i1) = a + x(i+1)
        i1 = i1 + 1
        y(i1) = ( x(i+1) - a ) * w
        w = w * z
      end do

      x(n1/2+is:n1+is-1) = y(1:n1/2)

    end do

  end do

  return
end
