program main

!*****************************************************************************80
!
!! MAIN is the main program for ANYPLT_PRB.
!
!  Discussion:
!
!    ANYPLT_PRB tests the ANYPLT routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) carray
  character ( len = 8 ) date
  integer ( kind = 4 ) icom
  integer ( kind = 4 ) iplt1
  integer ( kind = 4 ) iplt2
  integer ( kind = 4 ) ixplt1
  integer ( kind = 4 ) ixplt2
  integer ( kind = 4 ) iyplt1
  integer ( kind = 4 ) iyplt2
  integer ( kind = 4 ) marray
  character ( len = 10 ) time
  real xplt1
  real xplt2
  real yplt1
  real yplt2

  common /anycom/ iplt1, iplt2, ixplt1, ixplt2, iyplt1, &
                  iyplt2, marray, xplt1, xplt2, yplt1, yplt2
  common /anychr/ carray

  call date_and_time ( date, time )

  write ( *, * ) ' '
  write ( *, * ) 'ANYPLT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, * ) '  Test the ANYPLT graphics interface.'
  write ( *, * ) ' '
  write ( *, * ) '  Today''s date: ', date
  write ( *, * ) '  Today''s time: ', time
!
!  Enable graphics, and select window portion.
!  If we're using the ANYCGM interface, request PS output.
!
  icom = 0

  xplt1 = 0.0E+00
  xplt2 = 1.0E+00
  yplt1 = 0.0E+00
  yplt2 = 1.0E+00
  iplt1 = 1

  carray = 'PS'

  call anyplt ( icom )
!
!  Get the version number
!
  icom = 13
  call anyplt ( icom )

  write ( *, * ) ' '
  write ( *, * ) trim ( carray )

  call test01
  call test02
  call test03
!
!  End graphics
!
  icom = 1
  call anyplt ( icom )

  write ( *, * ) ' '
  write ( *, * ) 'ANYPLT_PRB'
  write ( *, * ) '  Normal end of ANYPLT tests.'

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 draws a sine curve.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
  implicit none

  character ( len = 80 ) carray
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icom
  integer ( kind = 4 ) iplt1
  integer ( kind = 4 ) iplt2
  integer ( kind = 4 ) ixplt1
  integer ( kind = 4 ) ixplt2
  integer ( kind = 4 ) iyplt1
  integer ( kind = 4 ) iyplt2
  integer ( kind = 4 ) marray
  real xlen
  real xmax
  real xmin
  real xplt1
  real xplt2
  real ylen
  real ymax
  real ymin
  real yplt1
  real yplt2

  common /anycom/ iplt1, iplt2, ixplt1, ixplt2, iyplt1, &
                  iyplt2, marray, xplt1, xplt2, yplt1, yplt2
  common /anychr/ carray

  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  Draw a sine curve.'
!
!  Begin plot
!
  icom = 2
  call anyplt ( icom )
!
!  Define user coordinate window
!
  xmin = -0.5E+00
  xlen = 11.0E+00
  xmax = xmin + xlen

  ymin = -1.2E+00
  ylen = 2.4E+00
  ymax = ymin + ylen

  icom = 3
  xplt1 = xmin
  xplt2 = xlen
  yplt1 = ymin
  yplt2 = ylen
  call anyplt ( icom )
!
!  Move to first point, and draw
!
  icom = 4
  do i = 1, 101
    xplt1 = real ( i - 1 ) / 10.0E+00
    yplt1 = sin ( xplt1 )
    call anyplt ( icom )
    icom = 5
  end do
!
!  Label every ten-th point
!
  carray = '.'
  do i = 1, 101, 10
    xplt1 = real ( i - 1 ) / 10.0E+00
    yplt1 = sin(xplt1)
    icom = 11
    call anyplt ( icom )
  end do
!
!  Print label
!
  icom = 7
  xplt1 = 0.10E+00
  xplt2 = 0.05E+00
  yplt1 = 0.85E+00
  yplt2 = 0.0
  carray = 'Sine curve'
  marray = 10
  call anyplt ( icom )
!
!  Label Y axis
!
  icom = 7
  xplt1 = 0.05E+00
  xplt2 = 0.025E+00
  yplt1 = 0.5E+00
  yplt2 = 90.0E+00
  carray = 'Y axis'
  marray = 6
  call anyplt ( icom )
!
!  Draw Y axis
!
  icom = 4
  xplt1 = 0.0E+00
  yplt1 = ymin
  call anyplt ( icom )
  icom = 5
  xplt1 = 0.0E+00
  yplt1 = ymax
  call anyplt ( icom )
!
!  Label X axis
!
  icom = 7
  xplt1 = 0.3E+00
  xplt2 = 0.025E+00
  yplt1 = 0.05E+00
  yplt2 = 0.0E+00
  carray = 'X axis'
  marray = 6
  call anyplt ( icom )
!
!  Draw X axis
!
  icom = 4
  xplt1 = xmin
  yplt1 = 0.0E+00
  call anyplt ( icom )
  icom = 5
  xplt1 = xmax
  yplt1 = 0.0E+00
  call anyplt ( icom )
!
!  End plot
!
  icom = 9
  call anyplt ( icom )

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
  implicit none

  real alpha
  character ( len = 80 ) carray
  real del
  real delt
  real dely
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icom
  integer ( kind = 4 ) iplt1
  integer ( kind = 4 ) iplt2
  integer ( kind = 4 ) ixplt1
  integer ( kind = 4 ) ixplt2
  integer ( kind = 4 ) iyplt1
  integer ( kind = 4 ) iyplt2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) marray
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real, parameter :: pi = 3.141592653589793E+00
  real theta
  real tnorm
  real u1
  real u2
  real uu
  real v1
  real v2
  real vv
  real x1
  real x2
  real xlen
  real xmax
  real xmin
  real xplt1
  real xplt2
  real y1
  real y2
  real ylen
  real ymax
  real ymin
  real yplt1
  real yplt2

  common /anycom/ iplt1, iplt2, ixplt1, ixplt2, iyplt1, &
                  iyplt2, marray, xplt1, xplt2, yplt1, yplt2
  common /anychr/ carray

  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  ???.'
!
!  Begin plot
!
  icom = 2
  call anyplt ( icom )
!
!  Define user coordinate window.
!  For this plot, make sure axes are of quite different lengths!
!
  xmin = 0.0E+00
  xlen = 10.0E+00
  xmax = xmin + xlen
  ymin = 0.0E+00
  ylen = 1.0E+00
  ymax = ymin + ylen

  icom = 3
  xplt1 = xmin
  xplt2 = xlen
  yplt1 = ymin
  yplt2 = ylen
  call anyplt ( icom )
!
!  Move to first point, and draw
!
  ny = 5
  nx = 5
  delt = 0.5E+00 * ( xmax - xmin ) / real ( nx - 1 )
  dely = 0.5E+00 * ( ymax - ymin ) / real ( ny - 1 )

  carray = '.'

  do i = 1, ny
    do j = 1, nx

      x1 = xmin + real ( j - 1 ) * ( xmax - xmin ) / real ( nx - 1 )
      y1 = ymin + real ( i - 1 ) * ( ymax - ymin ) / real ( ny - 1 )

      xplt1 = x1
      yplt1 = y1
      icom = 4
      call anyplt ( icom )

      x2 = x1 + delt
      y2 = y1 + dely * exp ( x1 * y1 ) / exp ( xmax * ymax )
      xplt1 = x2
      yplt1 = y2
      icom = 5
      call anyplt ( icom )

      uu = x2 - x1
      vv = y2 - y1
      tnorm = sqrt ( uu * uu + vv * vv )

      if ( 0.0E+00 < tnorm )then
        theta = pi - 0.5 * atan ( 2.0E+00 )
        alpha = atan2 ( vv, uu )
        del = sqrt ( 5.0E+00 ) * tnorm / 3.0E+00
        u1 = x1 + del * cos ( alpha - theta )
        v1 = y1 + del * sin ( alpha - theta )
        u2 = x1 + del * cos ( alpha + theta )
        v2 = y1 + del * sin ( alpha + theta )

        xplt1 = u1
        yplt1 = v1
        icom = 4
        call anyplt ( icom )

        xplt1 = x2
        yplt1 = y2
        icom = 5
        call anyplt ( icom )

        xplt1 = u2
        yplt1 = v2
        icom = 5
        call anyplt ( icom )

      end if

      write ( *, * ) i,j

    end do
  end do
!
!  Print label
!
  icom = 7
  xplt1 = 0.10E+00
  xplt2 = 0.05E+00
  yplt1 = 0.85E+00
  yplt2 = 0.0E+00
  carray = 'Comparison of arrows'
  marray = 20
  call anyplt ( icom )
!
!  Label Y axis
!
  icom = 7
  xplt1 = 0.05E+00
  xplt2 = 0.025E+00
  yplt1 = 0.5E+00
  yplt2 = 90.0E+00
  carray = 'y axis'
  marray = 6
  call anyplt ( icom )
!
!  Draw Y axis
!
  icom = 4
  xplt1 = 0.0E+00
  yplt1 = ymin
  call anyplt ( icom )
  icom = 5
  xplt1 = 0.0E+00
  yplt1 = ymax
  call anyplt ( icom )
!
!  Label X axis
!
  icom = 7
  xplt1 = 0.3E+00
  xplt2 = 0.025E+00
  yplt1 = 0.05E+00
  yplt2 = 0.0E+00
  carray = 'X axis'
  marray = 6
  call anyplt ( icom )
!
!  Draw X axis
!
  icom = 4
  xplt1 = xmin
  yplt1 = 0.0E+00
  call anyplt ( icom )
  icom = 5
  xplt1 = xmax
  yplt1 = 0.0E+00
  call anyplt ( icom )
!
!  End plot
!
  icom = 9
  call anyplt ( icom )

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
  implicit none

  character ( len = 80 ) carray
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icom
  integer ( kind = 4 ) iplt1
  integer ( kind = 4 ) iplt2
  integer ( kind = 4 ) ixplt1
  integer ( kind = 4 ) ixplt2
  integer ( kind = 4 ) iyplt1
  integer ( kind = 4 ) iyplt2
  integer ( kind = 4 ) marray
  real, parameter :: pi = 3.14159265E+00
  real xlen
  real xmax
  real xmin
  real xplt1
  real xplt2
  real ylen
  real ymax
  real ymin
  real yplt1
  real yplt2

  common /anycom/ iplt1, iplt2, ixplt1, ixplt2, iyplt1, &
                  iyplt2, marray, xplt1, xplt2, yplt1, yplt2
  common /anychr/ carray

  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  ???.'
!
!  Begin plot
!
  icom = 2
  call anyplt ( icom )
!
!  Define user coordinate window.
!  For this plot, make sure the axes are of quite different lengths!
!
  xmin = 0.0E+00
  xlen = 1.0E+00
  xmax = xmin + xlen

  ymin = 0.0E+00
  ylen = 1.0E+00
  ymax = ymin + ylen

  icom = 3
  xplt1 = xmin
  xplt2 = xlen
  yplt1 = ymin
  yplt2 = ylen
  call anyplt ( icom )
!
!  Draw the arrows.
!
  icom = 14
  do i = 1, 12
    xplt1 = 0.5E+00
    yplt1 = 0.5E+00
    xplt2 = real ( i - 1 ) * 2.0E+00 * pi / 12.0E+00
    yplt2 = 0.5E+00 * real ( i ) / 12.0E+00
    call anyplt ( icom )
  end do
!
!  Print label
!
  write ( *, * ) 'DEBUG: TEST03, ICOM = 7'
  icom = 7
  xplt1 = 0.10E+00
  xplt2 = 0.03E+00
  yplt1 = 0.85E+00
  yplt2 = 0.0E+00
  carray = 'Equal angles despite scale?'
  marray = 27
  call anyplt ( icom )
!
!  Label Y axis
!
  icom = 7
  xplt1 = 0.05E+00
  xplt2 = 0.025E+00
  yplt1 = 0.5E+00
  yplt2 = 90.0E+00
  carray = 'Y axis'
  marray = 6
  call anyplt ( icom )
!
!  Draw Y axis
!
  icom = 4
  xplt1 = xmin
  yplt1 = ymin
  call anyplt ( icom )
  icom = 5
  xplt1 = xmin
  yplt1 = ymax
  call anyplt ( icom )
!
!  Label X axis
!
  icom = 7
  xplt1 = 0.3E+00
  xplt2 = 0.025E+00
  yplt1 = 0.05E+00
  yplt2 = 0.0E+00
  carray = 'X axis'
  marray = 6
  call anyplt ( icom )
!
!  Draw X axis
!
  icom = 4
  xplt1 = xmin
  yplt1 = ymin
  call anyplt ( icom )
  icom = 5
  xplt1 = xmax
  yplt1 = ymin
  call anyplt ( icom )
!
!  End the plot.
!
  icom = 9
  call anyplt ( icom )

  return
end
