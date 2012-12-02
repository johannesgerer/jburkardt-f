subroutine f_alpha ( n, q_d, alpha, seed, x )

!*****************************************************************************80
!
!! F_ALPHA generates a 1/F^ALPHA noise sequence.
!
!  Discussion:
!
!    Thanks to Miro Stoyanov for pointing out that the second half of
!    the data returned by the inverse Fourier transform should be
!    discarded, 24 August 2010.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    Original C version by Todd Walter.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of samples to generate.
!
!    Input, real ( kind = 8 ) Q_D, the variance of the noise.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent for the noise.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sequence sampled with the given
!    power law.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) h_a(n)
  real ( kind = 8 ) h_azero
  real ( kind = 8 ) h_b(n)
  real ( kind = 8 ) hfa(2*n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) q_d
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w_a(n)
  real ( kind = 8 ) w_azero
  real ( kind = 8 ) w_b(n)
  real ( kind = 8 ) wfa(2*n)
  real ( kind = 8 ) wi
  real ( kind = 8 ) wr
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x2(2*n)
!
!  Set the deviation of the noise.
!
  q_d = sqrt ( q_d )
!
!  Generate the coefficients Hk.
!
  hfa(1) = 1.0D+00
  do i = 2, n
    hfa(i) = hfa(i-1) * ( 0.5D+00 * alpha + real ( i - 2, kind = 8 ) ) &
      / ( real ( i - 1, kind = 8 ) )
  end do
  hfa(n+1:2*n) = 0.0D+00
!
!  Fill Wk with white noise.
!
  do i = 1, n
    wfa(i) = q_d * r8_normal_01 ( seed )
  end do
  wfa(n+1:2*n) = 0.0D+00
!
!  Perform the discrete Fourier transforms of Hk and Wk.
!
  call r8vec_sftf ( 2 * n, hfa, h_azero, h_a, h_b )

  call r8vec_sftf ( 2 * n, wfa, w_azero, w_a, w_b )
!
!  Multiply the two complex vectors.
!
  w_azero = w_azero * h_azero

  do i = 1, n
    wr = w_a(i)
    wi = w_b(i)
    w_a(i) = wr * h_a(i) - wi * h_b(i)
    w_b(i) = wi * h_a(i) + wr * h_b(i)
  end do
!
!  This scaling is introduced only to match the behavior
!  of the Numerical Recipes code...
!
  w_azero = w_azero * real ( 2 * n, kind = 8 )

  w_a(1:n-1) = w_a(1:n-1) * real ( n, kind = 8 )
  w_b(1:n-1) = w_b(1:n-1) * real ( n, kind = 8 )

  w_a(n) = w_a(n) * real ( 2 * n, kind = 8 )
  w_b(n) = w_b(n) * real ( 2 * n, kind = 8 )
!
!  Take the inverse Fourier transform of the result.
!
  call r8vec_sftb ( 2 * n, w_azero, w_a, w_b, x2 )
!
!  Only return the first N inverse Fourier transform values.
!
  x(1:n) = x2(1:n)

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen
 
  iunit = 0
 
  do i = 1, 99
 
    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )
 
      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if
 
  end do

  return
end
function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!    This means that on every other call, the code can use the second
!    value that it calculated.
!
!    However, if the user has changed the SEED value between calls,
!    the routine automatically resets itself and discards the saved data.
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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a normally distributed
!    random value.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed1 = 0
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: seed3 = 0
  integer ( kind = 4 ), parameter :: two = 2
  integer ( kind = 4 ), save :: used = 0
  real ( kind = 8 ) v1
  real ( kind = 8 ), save :: v2 = 0.0D+00
!
!  If USED is odd, but the input SEED does not match
!  the output SEED on the previous call, then the user has changed
!  the seed.  Wipe out internal memory.
!
  if ( mod ( used, two ) == 1 ) then

    if ( seed /= seed2 ) then
      used = 0
      seed1 = 0
      seed2 = 0
      seed3 = 0
      v2 = 0.0D+00
    end if

  end if
!
!  If USED is even, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, two ) == 0 ) then

    seed1 = seed

    r1 = r8_uniform_01 ( seed )

    if ( r1 == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed

    r2 = r8_uniform_01 ( seed )

    seed3 = seed

    v1 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    v2 = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )

    r8_normal_01 = v1
    seed = seed2
!
!  If USED is odd (and the input SEED matched the output value from
!  the previous call), return the second normal and its corresponding seed.
!
  else

    r8_normal_01 = v2
    seed = seed3

  end if

  used = used + 1

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
!    For now, the input quantity SEED is an integer variable.
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_PART prints "part" of an R8VEC.
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
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ......  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end
subroutine r8vec_sftb ( n, azero, a, b, r )

!*****************************************************************************80
!
!! R8VEC_SFTB computes a "slow" backward Fourier transform of real data.
!
!  Discussion:
!
!    SFTB and SFTF are inverses of each other.  If we begin with data
!    AZERO, A, and B, and apply SFTB to it, and then apply SFTF to the
!    resulting R vector, we should get back the original AZERO, A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) AZERO, the constant Fourier coefficient.
!
!    Input, real ( kind = 8 ) A(N/2), B(N/2), the Fourier coefficients.
!
!    Output, real ( kind = 8 ) R(N), the reconstructed data sequence.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n/2)
  real ( kind = 8 ) azero
  real ( kind = 8 ) b(n/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) r(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta

  r(1:n) = azero
  do i = 1, n
    do k = 1, n / 2
      theta = real ( k * ( i - 1 ) * 2, kind = 8 ) * pi &
        / real ( n, kind = 8 )
      r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
    end do
  end do

  return
end
subroutine r8vec_sftf ( n, r, azero, a, b )

!*****************************************************************************80
!
!! R8VEC_SFTF computes a "slow" forward Fourier transform of real data.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    R and apply SFTB to it, and then apply SFTB to the resulting AZERO, 
!    A, and B, we should get back the original R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) R(N), the data to be transformed.
!
!    Output, real ( kind = 8 ) AZERO, = sum ( 1 <= I <= N ) R(I) / N.
!
!    Output, real ( kind = 8 ) A(N/2), B(N/2), the Fourier coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(1:n/2)
  real ( kind = 8 ) azero
  real ( kind = 8 ) b(1:n/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) theta

  azero = sum ( r(1:n) ) / real ( n, kind = 8 )

  do i = 1, n / 2

    a(i) = 0.0D+00
    b(i) = 0.0D+00

    do j = 1, n
      theta = real ( 2 * i * ( j - 1 ), kind = 8 ) * pi &
        / real ( n, kind = 8 )
      a(i) = a(i) + r(j) * cos ( theta )
      b(i) = b(i) + r(j) * sin ( theta )
    end do

    a(i) = a(i) / real ( n, kind = 8 )
    b(i) = b(i) / real ( n, kind = 8 )

    if ( i /= ( n / 2 ) ) then
      a(i) = 2.0D+00 * a(i)
      b(i) = 2.0D+00 * b(i)
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
