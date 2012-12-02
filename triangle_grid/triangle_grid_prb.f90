program main

!*****************************************************************************80
!
!! MAIN is the main program for TRIANGLE_GRID_PRB.
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
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_GRID_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TRIANGLE_GRID library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_GRID_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests TRIANGLE_GRID.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ng = ((n+1)*(n+2))/2

  character ( len = 255 ) filename
  integer ( kind = 4 ) j
  real ( kind = 8 ) t(2,3)
  real ( kind = 8 ) tg(2,ng)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  TRIANGLE_GRID can define a triangular grid of points'
  write ( *, '(a)' ) '  with N+1 points on a side, based on any triangle.'

  t = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.5D+00, 0.86602540378443860D+00 /), (/ 2, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Defining triangle:'
  write ( *, '(a)' ) '     J      X      Y'
  write ( *, '(a)' ) ' '
  do j = 1, 3
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) j, t(1:2,j)
  end do
  call triangle_grid ( n, t, tg )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     J      X      Y'
  write ( *, '(a)' ) ' '
  do j = 1, ng
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) j, tg(1:2,j)
  end do

  filename = 'triangle_grid_test01.xy'

  open ( unit = 1, file = filename, status = 'replace' )
  do j = 1, ng
    write ( 1, '(2x,g14.6,2x,g14.6)' ) tg(1:2,j)
  end do
  close ( unit = 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to "' // trim ( filename ) // '".'

  return
end

