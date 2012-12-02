function randlc ( x )

!*****************************************************************************80
!
!! RANDLC returns a uniform pseudorandom value.
!
!  Discussion:
!
!    The number returned is in the range (0, 1).  
!
!    The algorithm uses the linear congruential generator:
!
!      X(K+1) = A * X(K)  mod 2^46
!
!    This scheme generates 2^44 numbers before repeating.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Eric Barszcz, John Barton, D Browning, Robert Carter, 
!    Leonardo Dagum, Rod Fatoohi,
!    Samuel Fineberg, Paul Frederickson, Thomas Lasinski, Robert Schreiber, 
!    Horst Simon, V Venkatakrishnan, Sisira Weeratunga,
!    The NAS Parallel Benchmarks,
!    RNR Technical Report RNR-94-007,
!    March 1994.
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 2, Seminumerical Algorithms,
!    Third Edition,
!    Addison Wesley, 1997,
!    ISBN: 0201896842,
!    LC: QA76.6.K64.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, the seed.  X should be an 
!    odd integer such that 1 <= X <= 2^46.
!
!    Output, real ( kind = 8 ) RANDLC, the next pseudorandom value.
!
  implicit none

  real ( kind = 8 ), parameter :: a = 1220703125.0D+00
  real ( kind = 8 ), save :: a1
  real ( kind = 8 ), save :: a2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: ks = 0
  real ( kind = 8 ), save :: r23
  real ( kind = 8 ), save :: r46
  real ( kind = 8 ) randlc
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ), save :: t23
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ), save :: t46
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) z
!
!  If this is the first call, compute 
!
!    R23 = 2 ^ -23, 
!    R46 = 2 ^ -46,
!    T23 = 2 ^ 23, 
!    T46 = 2 ^ 46.  
!
!  These are computed in loops, rather than by merely using the power operator, 
!  in order to insure that the results are exact on all systems.  
!
  if ( ks == 0 ) then

    r23 = 1.0D+00
    r46 = 1.0D+00
    t23 = 1.0D+00
    t46 = 1.0D+00

    do i = 1, 23
      r23 = 0.5D+00 * r23
      t23 = 2.0D+00 * t23
    end do

    do i = 1, 46
      r46 = 0.50D+00 * r46
      t46 = 2.0D+00 * t46
    end do
!
!  Break A into two parts such that A = 2^23 * A1 + A2.
!
    t1 = r23 * a
    a1 = real ( int ( t1 ), kind = 8 )
    a2 = a - t23 * a1

    ks = 1

  end if
!
!  Deal with a 0 input value of X.
!
  if ( x == 0.0D+00 ) then
    x = 314159265.0D+00
  end if
!
!  Deal somewhat arbitrarily with negative input X.
!
  if ( x < 0.0D+00 ) then
    x = - x
  end if
!
!  Break X into two parts X1 and X2 such that:
!
!    X = 2^23 * X1 + X2, 
!
!  then compute
!
!    Z = A1 * X2 + A2 * X1  (mod 2^23)
!    X = 2^23 * Z + A2 * X2  (mod 2^46).
!
  t1 = r23 * x
  x1 = real ( int ( t1 ), kind = 8 )
  x2 = x - t23 * x1

  t1 = a1 * x2 + a2 * x1
  t2 = real ( int ( r23 * t1 ), kind = 8 )
  z = t1 - t23 * t2

  t3 = t23 * z + a2 * x2
  t4 = real ( int ( r46 * t3 ), kind = 8 )
  x = t3 - t46 * t4

  randlc = r46 * x

  return
end
function randlc_jump ( x, k )

!*****************************************************************************80
!
!! RANDLC_JUMP returns the K-th element of a uniform pseudorandom sequence.
!
!  Discussion:
!
!    The sequence uses the linear congruential generator:
!
!      X(K+1) = A * X(K)  mod 2^46
!
!    The K-th element, which can be represented as
!
!      X(K) = A^K * X(0)  mod 2^46
!
!    is computed directly using the binary algorithm for exponentiation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Eric Barszcz, John Barton, D Browning, Robert Carter, 
!    Leonardo Dagum, Rod Fatoohi,
!    Samuel Fineberg, Paul Frederickson, Thomas Lasinski, Robert Schreiber, 
!    Horst Simon, V Venkatakrishnan, Sisira Weeratunga,
!    The NAS Parallel Benchmarks,
!    RNR Technical Report RNR-94-007,
!    March 1994.
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 2, Seminumerical Algorithms,
!    Third Edition,
!    Addison Wesley, 1997,
!    ISBN: 0201896842,
!    LC: QA76.6.K64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the initial seed (with index 0).  
!
!    Input, integer ( kind = 4 ) K, the index of the desired value.
!
!    Output, real ( kind = 8 ) RANDLC_JUMP, the K-th value in the sequence.
!
  implicit none

  real ( kind = 8 ), parameter :: a = 1220703125.0D+00
  real ( kind = 8 ), save :: a1
  real ( kind = 8 ), save :: a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k2
  integer ( kind = 4 ), save :: ks = 0
  integer ( kind = 4 ) m
  real ( kind = 8 ), save :: r23
  real ( kind = 8 ), save :: r46
  real ( kind = 8 ) randlc_jump
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ), save :: t23
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ), save :: t46
  integer ( kind = 4 ) twom
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xk
  real ( kind = 8 ) z
!
!  If this is the first call, compute 
!
!    R23 = 2 ^ -23, 
!    R46 = 2 ^ -46,
!    T23 = 2 ^ 23, 
!    T46 = 2 ^ 46.  
!
!  These are computed in loops, rather than by merely using the power operator, 
!  in order to insure that the results are exact on all systems.  
!
  if ( ks == 0 ) then

    r23 = 1.0D+00
    r46 = 1.0D+00
    t23 = 1.0D+00
    t46 = 1.0D+00

    do i = 1, 23
      r23 = 0.5D+00 * r23
      t23 = 2.0D+00 * t23
    end do

    do i = 1, 46
      r46 = 0.50D+00 * r46
      t46 = 2.0D+00 * t46
    end do
!
!  Break A into two parts such that A = 2^23 * A1 + A2.
!
    t1 = r23 * a
    a1 = real ( int ( t1 ), kind = 8 )
    a2 = a - t23 * a1

    ks = 1

  end if

  if ( k < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDLC_JUMP - Fatal error!'
    write ( *, '(a)' ) '  K < 0.'
    stop

  else if ( k == 0 ) then

    xk = x
!
!  Find M so that K < 2^M.
!
  else

    k2 = k
    xk = x

    m = 1
    twom = 2
    do while ( twom <= k )
      twom = twom * 2
      m = m + 1
    end do

    b = a
    b1 = a1
    b2 = a2

    do i = 1, m

      j = k2 / 2
!
!  Replace X by A * X, if appropriate.
!
      if ( 2 * j /= k2 ) then

        t1 = r23 * xk
        x1 = real ( int ( t1 ), kind = 8 )
        x2 = xk - t23 * x1

        t1 = b1 * x2 + b2 * x1
        t2 = real ( int ( r23 * t1 ), kind = 8 )
        z = t1 - t23 * t2

        t3 = t23 * z + b2 * x2
        t4 = real ( int ( r46 * t3 ), kind = 8 )
        xk = t3 - t46 * t4

      end if
!
!  Replace A by A * A mod 2^46.
!
      t1 = r23 * b
      x1 = real ( int ( t1 ), kind = 8 )
      x2 = b - t23 * x1

      t1 = b1 * x2 + b2 * x1
      t2 = real ( int ( r23 * t1 ), kind = 8 )
      z = t1 - t23 * t2

      t3 = t23 * z + b2 * x2
      t4 = real ( int ( r46 * t3 ), kind = 8 )
      b = t3 - t46 * t4
!
!  Update A1, A2.
!
      t1 = r23 * b
      b1 = real ( int ( t1 ), kind = 8 )
      b2 = b - t23 * b1

      k2 = j

    end do

  end if

  randlc_jump = r46 * xk

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
