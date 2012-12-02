program main

!*****************************************************************************80
!
!! MAIN is the main program for SPIRAL_DATA.
!
!  Discussion:
!
!    Note that the continuous velocity field (U,V)(X,Y) that is discretely
!    sampled here satisfies the homogeneous continuity equation, that is,
!    it has zero divergence.  In other words:
!
!      dU/dX + dV/dY = 0.
!
!    This is by construction, since we have
!
!      U(X,Y) =  10 * d/dY ( PHI(X) * PHI(Y) )
!      V(X,Y) = -10 * d/dX ( PHI(X) * PHI(Y) )
!
!    which guarantees zero divergence.
!
!    The underlying function PHI is defined by
!
!      PHI(Z) = ( 1 - cos ( C * pi * Z ) ) * ( 1 - Z )**2
!
!    where C is a parameter.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  character ( len = 255 ) :: uv_file = 'spiral_uv.txt'
  real ( kind = 8 ) v
  real ( kind = 8 ) x
  character ( len = 255 ) :: xy_file = 'spiral_xy.txt'
  real ( kind = 8 ) y

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPIRAL_DATA'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Write data describing a "spiral" vector flow'
  write ( *, '(a)' ) '  which satisfies the continuity equation.'
  write ( *, '(a)' ) '  (Actually, the flow does not spiral, '
  write ( *, '(a)' ) '  it''s simply concentric loops.)'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The flow is defined on an N+1 x N+1 grid.'
  write ( *, '(a)' ) '  A typical value might be 10 or 20.'
  write ( *, '(a)' ) '  Enter your desired value of N'
  read ( *, * ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The flow depends on a parameter C'
  write ( *, '(a)' ) '  A typical value might be 0.75 or 1.0'
  write ( *, '(a)' ) '  Enter your desired value of C.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Generating data on N+1 x N+1 grid, N = ', n
  write ( *, '(a,g14.6)' ) '  Using parameter C = ', c

  open ( unit = 1, file = xy_file, status = 'replace' )
  open ( unit = 2, file = uv_file, status = 'replace' )

  do i = 0, n

    x = real ( i, kind = 8 ) / real ( n, kind = 8 )

    do j = 0, n

      y = real ( j, kind = 8 ) / real ( n, kind = 8 )

      u =   10 * ( 1.0D+00 - cos ( c * pi * x ) )             &
               * ( 1.0D+00 - x )**2                           &
               * (                                            &
                   c * pi * sin ( c * pi * y ) * ( 1 - y )**2 &
                 - ( 1 - cos ( c * pi * y ) ) * 2 * ( 1 - y ) &
                 )

      v = - 10 * ( 1.0D+00 - cos ( c * pi * y ) )             &
               * ( 1.0D+00 - y )**2                           &
               * (                                            &
                   c * pi * sin ( c * pi * x ) * ( 1 - x )**2 &
                 - ( 1 - cos ( c * pi * x ) ) * 2 * ( 1 - x ) &
                 )

      write ( 1, '(2x,f10.4,2x,f10.4)' ) x, y
      write ( 2, '(2x,f10.4,2x,f10.4)' ) u, v

    end do
  end do

  close ( unit = 1 )
  close ( unit = 2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XY data written to "' // trim ( xy_file ) // '".'
  write ( *, '(a)' ) '  UV data written to "' // trim ( uv_file ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This data can be plotted by:'
  write ( *, '(a)' ) '    the FORTRAN90 program VECTOR_PLOT or'
  write ( *, '(a)' ) '    the MATLAB program DIRECTION_ARROWS or'
  write ( *, '(a)' ) '    the MATLAB program VELOCITY_ARROWS.'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPIRAL_DATA'
  write ( *, '(a)' ) '  Normal end of execution'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
