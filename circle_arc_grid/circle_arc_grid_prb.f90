program main

!*****************************************************************************80
!
!! MAIN is the main program for CIRCLE_ARC_GRID_PRB.
!
!  Discussion:
!
!    CIRCLE_ARC_GRID_PRB tests CIRCLE_ARC_GRID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_ARC_GRID_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test CIRCLE_ARC_GRID.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_ARC_GRID_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of CIRCLE_ARC_GRID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) c(2)
  character ( len = 80 ) filename
  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  real ( kind = 8 ), allocatable :: xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_POINTS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Write points on a circle to a file.'

  r = 2.0D+00
  c(1) = 5.0D+00
  c(2) = 5.0D+00
  a(1) = 0.0D+00
  a(2) = 90.0D+00
  n = 10
!
!  Echo the input.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Radius =           ', r
  write ( *, '(a,g14.6,2x,g14.6)' ) '  Center =           ', c(1), c(2)
  write ( *, '(a,g14.6)' ) '  Angle 1 =          ', a(1)
  write ( *, '(a,g14.6)' ) '  Angle 2 =          ', a(2)
  write ( *, '(a,i8)' ) '  Number of points = ', n

  allocate ( xy(1:2,1:n) )
!
!  Compute the data.
!
  call circle_arc_grid ( r, c, a, n, xy )
!
!  Print a little of the data.
!
  call r82vec_print_part ( n, xy, 5, '  A few of the points:' )
!
!  Write the data.
!
  filename = 'arc.txt'
  call r8mat_write ( filename, 2, n, xy )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to "' // trim ( filename ) // '".'
!
!  Free memory.
!
  deallocate ( xy )

  return
end
