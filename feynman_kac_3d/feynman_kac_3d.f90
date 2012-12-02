program main

!*****************************************************************************80
!
!! MAIN is the main program for FEYNMAN_KAC_3D.
!
!  Discussion:
!
!    This program is derived from section 2.5, exercise 2.2 of Petersen 
!    and Arbenz.
!
!    The problem is to determine the solution U(X,Y,Z) of the following 
!    partial differential equation:
!
!      (1/2) Laplacian U - V(X,Y,Z) * U = 0,
!
!    inside the elliptic domain D:
! 
!      D = { (X,Y,Z) | (X/A)^2+(Y/B)^2+(Z/C)^2 <= 1 }
!   
!    with the boundary condition U(boundary(D)) = 1.
!
!    The V(X,Y,Z) is the potential function:
!
!      V = 2 * ( (X/A^2)^2 + (Y/B^2)^2 + (Z/C^2)^2 ) + 1/A^2 + 1/B^2 + 1/C^2.
!
!    The analytic solution of this problem is already known:
!
!      U(X,Y,Z) = exp ( (X/A)^2 + (Y/B)^2 + (Z/C)^2 - 1 ).
!
!    Our method is via the Feynman-Kac Formula.
!
!    The idea is to start from any (x,y,z) in D, and
!    compute (x+Wx(t),y+Wy(t),z+Wz(t)) where 3-D Brownian motion
!    (Wx,Wy,Wz) is updated each step by sqrt(h)*(z1,z2,z3),
!    each z1,z2,z3 are independent approximately Gaussian 
!    random variables with zero mean and variance 1. 
!
!    Each (x1(t),x2(t),x3(t)) is advanced until (x1,x2,x3) exits 
!    the domain D.  
!
!    Upon its first exit from D, the sample path (x1,x2,x3) is stopped and a 
!    new sample path at (x,y,z) is started until N such paths are completed.
! 
!    The Feynman-Kac formula gives the solution here as
!
!      U(X,Y,Z) = (1/N) sum(1 <= I <= N) Y(tau_i),
!
!    where
!
!      Y(tau) = exp( -int(s=0..tau) v(x1(s),x2(s),x3(s)) ds),
!
!    and tau = first exit time for path (x1,x2,x3). 
!
!    The integration procedure is a second order weak accurate method:
!
!      X(t+h)  = [ x1(t) + sqrt ( h ) * z1 ]
!                [ x2(t) + sqrt ( h ) * z2 ]
!                [ x3(t) + sqrt ( h ) * z3 ]
!
!    Here Z1, Z2, and Z3 are approximately normal univariate Gaussians. 
!
!    An Euler predictor approximates Y at the end of the step
!
!      Y_e     = (1 - h*v(X(t)) * Y(t), 
!
!    A trapezoidal rule completes the step:
!
!      Y(t+h)  = Y(t) - (h/2)*[v(X(t+h))*Y_e + v(X(t))*Y(t)].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2012
!
!  Author:
!
!    Original C version by Wesley Petersen.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Peter Arbenz, Wesley Petersen,
!    Introduction to Parallel Computing:
!    A Practical Guide with Examples in C,
!    Oxford, 2004,
!    ISBN: 0-19-851577-4,
!    LC: QA76.59.P47.
!
  implicit none

  real ( kind = 8 ) :: a = 3.0D+00
  real ( kind = 8 ) :: b = 2.0D+00
  real ( kind = 8 ) :: c = 1.0D+00
  real ( kind = 8 ) chk
  integer ( kind = 4 ) :: dim = 3
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) dz
  real ( kind = 8 ) err
  real ( kind = 8 ) :: h = 0.001D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_ceiling
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) :: n = 10000
  integer ( kind = 4 ) n_inside
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nk
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) rth
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) steps
  integer ( kind = 4 ) steps_ave
  integer ( kind = 4 ) trial
  real ( kind = 8 ) us
  real ( kind = 8 ) ut
  real ( kind = 8 ) vh
  real ( kind = 8 ) vs
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y
  real ( kind = 8 ) w
  real ( kind = 8 ) w_exact
  real ( kind = 8 ) we
  real ( kind = 8 ) wt
  real ( kind = 8 ) z
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2
  real ( kind = 8 ) z3

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEYNMAN_KAC_3D:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Program parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The calculation takes place inside a 3D ellipsoid.'
  write ( *, '(a)' ) '  A rectangular grid of points will be defined.'
  write ( *, '(a)' ) '  The solution will be estimated for those grid points'
  write ( *, '(a)' ) '  that lie inside the ellipsoid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) &
    '  Each solution will be estimated by computing ', n, &
    ' trajectories'
  write ( *, '(a)' ) '  from the point to the boundary.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    (X/A)^2 + (Y/B)^2 + (Z/C)^2 = 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The ellipsoid parameters A, B, C are set to:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    A = ', a
  write ( *, '(a,g14.6)' ) '    B = ', b
  write ( *, '(a,g14.6)' ) '    C = ', c
  write ( *, '(a,g14.6)' ) '  Stepsize H = ', h
!
!  RTH is the scaled stepsize.
!
  rth = sqrt ( real ( dim, kind = 8 ) * h )
!
!  Choose the spacing so we have about 10 points in the shortest direction.
!
  if ( a == min ( min ( a, b ), c ) ) then
    ni = 11
    nj = 1 + i4_ceiling ( b / a ) * ( ni - 1 )
    nk = 1 + i4_ceiling ( c / a ) * ( ni - 1 )
  else if ( b == min ( min ( a, b ), c ) ) then
    nj = 11
    ni = 1 + i4_ceiling ( a / b ) * ( nj - 1 )
    nk = 1 + i4_ceiling ( c / b ) * ( nj - 1 )
  else
    nk = 11
    ni = 1 + i4_ceiling ( a / c ) * ( nk - 1 )
    nj = 1 + i4_ceiling ( b / c ) * ( nk - 1 )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a)' ) '  X coordinate marked by ', ni, ' points'
  write ( *, '(a,i4,a)' ) '  Y coordinate marked by ', nj, ' points'
  write ( *, '(a,i4,a)' ) '  Z coordinate marked by ', nk, ' points'
!
!  Loop over the grid points.
!
  err = 0.0D+00
  n_inside = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     X            Y             Z          W Approx' // &
    '      W Exact          Error      Ave Steps'
  write ( *, '(a)' ) ' '

  do i = 1, ni

    x = ( real ( ni - i,     kind = 8 ) * ( - a ) &
        + real (      i - 1, kind = 8 ) *     a ) &
        / real ( ni     - 1, kind = 8 )

    do j = 1, nj

      y = ( real ( nj - j,     kind = 8 ) * ( - b ) &
          + real (      j - 1, kind = 8 ) *     b ) & 
          / real ( nj     - 1, kind = 8 )

      do k = 1, nk

        z = ( real ( nk - k,     kind = 8 ) * ( - c ) &
            + real (      k - 1, kind = 8 ) *     c ) & 
            / real ( nk     - 1, kind = 8 )

        chk = ( x / a )**2 + ( y / b )**2 + ( z / c )**2

        if ( 1.0D+00 < chk ) then
          w_exact = 1.0D+00
          wt = 1.0D+00
          steps_ave = 0
          write ( *, &
            '(6(2x,g12.4),2x,i8)' ) &
            x, y, z, wt, w_exact, abs ( w_exact - wt ), steps_ave
          cycle
        end if

        n_inside = n_inside + 1
!
!  Compute the exact solution at this point (x,y,z).
!
        w_exact = exp ( ( x / a )**2 + ( y / b )**2 + ( z / c )**2 - 1.0D+00 )
!
!  Now try to estimate the solution at this point.
!
        wt = 0.0D+00
        steps = 0

        do trial = 1, n

          x1 = x
          x2 = y
          x3 = z
! 
!  W = exp(-int(s=0..t) v(X)ds) 
!
          w = 1.0D+00
!
!  CHK is < 1.0 while the point is inside the ellipsoid.
!
          chk = 0.0D+00
 
          do while ( chk < 1.0D+00 )
!
!  Determine DX, DY, DZ.
!
            ut = r8_uniform_01 ( seed )
            if ( ut < 1.0D+00 / 3.0D+00 ) then
              us = r8_uniform_01 ( seed ) - 0.5D+00
              if ( us < 0.0D+00 ) then
                dx = - rth
              else
                dx = + rth
              end if
            else
              dx = 0.0D+00
            end if

            ut = r8_uniform_01 ( seed )
            if ( ut < 1.0D+00 / 3.0D+00 ) then
              us = r8_uniform_01 ( seed ) - 0.5D+00
              if ( us < 0.0D+00 ) then
                dy = - rth
              else
                dy = + rth
              end if
            else
              dy = 0.0D+00
            end if

            ut = r8_uniform_01 ( seed )
            if ( ut < 1.0D+00 / 3.0D+00 ) then
              us = r8_uniform_01 ( seed ) - 0.5D+00
              if ( us < 0.0D+00 ) then
                dz = - rth
              else
                dz = + rth
              end if
            else
              dz = 0.0D+00
            end if

            call potential ( a, b, c, x1, x2, x3, vs )
! 
!  Move to the new point.
!
            x1 = x1 + dx
            x2 = x2 + dy
            x3 = x3 + dz

            steps = steps + 1
  
            call potential ( a, b, c, x1, x2, x3, vh )

            we = ( 1.0D+00 - h * vs ) * w
            w = w - 0.5D+00 * h * ( vh * we + vs * w ) 

            chk = ( x1 / a )**2 + ( x2 / b )**2 + ( x3 / c )**2

          end do

          wt = wt + w

        end do
!
!  WT is the average of the sum of the different trials.
!
        wt = wt / real ( n, kind = 8 )
        steps_ave = steps / n
!
!  Add error in WT to the running L2 error in the solution.
!
        err = err + ( w_exact - wt )**2

        write ( *, &
          '(6(2x,g12.4),2x,i8)' ) &
          x, y, z, wt, w_exact, abs ( w_exact - wt ), steps_ave

      end do
    end do
  end do
!
!  Compute the RMS error for all the points.
!
  err = sqrt ( err / real ( n_inside, kind = 8 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  RMS absolute error in solution = ', err
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEYNMAN_KAC_3D:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '

  call timestamp ( )

  stop
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
subroutine potential ( a, b, c, x, y, z, v )

!*****************************************************************************80
!
!! POTENTIAL evaluates the potential function V(X,Y,Z).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters that define the ellipse.
!
!    Input, real ( kind = 8 ) X, Y, Z, the coordinates of the point.
!
!    Output, real ( kind = 8 ) V, the value of the potential function at (X,Y,Z).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) v
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  v = 2.0D+00 * ( ( x / a**2 )**2 + ( y / b**2 )**2 + ( z / c**2 )**2 ) &
    + 1.0D+00 / a**2 + 1.0D+00 / b**2 + 1.0D+00 / c**2

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
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
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

