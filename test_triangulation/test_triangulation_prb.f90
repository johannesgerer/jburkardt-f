program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_TRIANGULATION_PRB.
!
!  Discussion:
!
!    TEST_TRIANGULATION_PRB runs the TEST_TRIANGULATION tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TRIANGULATION_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_TRIANGULATION library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TRIANGULATION_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests P00_TEST_NUM, P00_TITLE and P00_HEADER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  P00_TEST_NUM reports the number of problems.'
  write ( *, '(a)' ) '  P00_TITLE returns a title for each problem.'
  write ( *, '(a)' ) '  P00_HEADER prints some information about each problem.'

  call p00_test_num ( test_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of tests available = ', test_num

  do test = 1, test_num

    call p00_title ( test, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Test number:        ', test
    write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

    call p00_header ( test )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests P00_TEST_NUM and P00_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_start = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  P00_TEST_NUM reports the number of problems.'
  write ( *, '(a)' ) '  P00_SAMPLE returns sample points from the region.'

  call p00_test_num ( test_num )

  do test = 1, test_num

    seed = seed_start

    call p00_title ( test, title )

    n = 20

    write ( *, '(a)'     ) ' '
    write ( *, '(a,i6)'  ) '  Test number         =  ', test
    write ( *, '(a, a)'  ) '  Title:              = "', trim ( title ) // '"'
    write ( *, '(a,i6)'  ) '  Spatial dimension M =  ', m
    write ( *, '(a,i6)'  ) '  Number of samples N =  ', n
    write ( *, '(a,i12)' ) '  Initial SEED:       =  ', seed
    write ( *, '(a)'     ) ' '

    allocate ( point(1:m,1:n) )

    call p00_sample ( test, m, n, seed, point )

    call r8mat_transpose_print ( m, n, point, '  The sample points:' )

    deallocate ( point )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests P00_BOUNDARY_NEAREST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: boundary
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) :: seed_start = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  P00_BOUNDARY_NEAREST returns the nearest boundary'
  write ( *, '(a)' ) '  point for a set of points.'
  write ( *, '(a)' ) ' '

  call p00_test_num ( test_num )

  do test = 1, test_num

    seed = seed_start

    call p00_title ( test, title )

    n = 20

    write ( *, '(a)'     ) ' '
    write ( *, '(a,i6)'  ) '  Test number         =  ', test
    write ( *, '(a, a)'  ) '  Title:              = "', trim ( title ) // '"'
    write ( *, '(a,i6)'  ) '  Spatial dimension M =  ', m
    write ( *, '(a,i6)'  ) '  Number of samples N =  ', n
    write ( *, '(a,i12)' ) '  Initial SEED:       =  ', seed
    write ( *, '(a)'     ) ' '

    allocate ( point(1:m,1:n) )
    allocate ( boundary(1:m,1:n) )

    call p00_sample ( test, m, n, seed, point )

    call p00_boundary_nearest ( test, m, n, point, boundary )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Pairs of Point and Nearest Boundary Point:'
    write ( *, '(a)' ) '  (This data can be passed to the PLOT_POINTS'
    write ( *, '(a)' ) '  routine and plotted with the DASH option.)'
    write ( *, '(a)' ) ' '

    do j = 1, n
      write ( *, '(2x,g14.6,2x,g14.6)' ) point(1:2,j)
      write ( *, '(2x,g14.6,2x,g14.6)' ) boundary(1:2,j)
      write ( *, '(a)' ) ' '
    end do

    deallocate ( point )
    deallocate ( boundary )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests P00_BOUNDARY_EPS.
!
!  Discussion:
!
!    This routine creates an EPS file containing an image of the boundary.
!
!    The boundary is broken up into segments of about H = 0.05 in length.
!
!    If SHOW_NODES is true, these points are shown in the image.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  real ( kind = 8 ) :: box(m,2)
  character ( len = 80 ) eps_file_name
  real ( kind = 8 ) h
  real ( kind = 8 ) hi(m)
  real ( kind = 8 ) :: h_relative = 0.05D+00
  real ( kind = 8 ) lo(m)
  logical, parameter :: show_nodes = .false.
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  P00_BOUNDARY_EPS makes an EPS file containing'
  write ( *, '(a)' ) '  an image of the boundary.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The boundary will be drawn to an absolute fineness'
  write ( *, '(a,g14.6)' ) '  of H = ', h
  write ( *, '(a)' ) ' '

  if ( show_nodes ) then
    write ( *, '(a)' ) '  The boundary nodes will be shown.'
  else
    write ( *, '(a)' ) '  The boundary nodes will NOT be shown.'
  end if

  eps_file_name = 'p00_boundary.eps'

  call p00_test_num ( test_num )

  do test = 1, test_num

    call p00_box ( test, m, lo, hi )

    h = h_relative * min ( hi(1) - lo(1), hi(2) - lo(2) )

    call p00_title ( test, title )

    write ( *, '(a)'       ) ' '
    write ( *, '(a,i6)'    ) '  Test number         = ', test
    write ( *, '(a, a)'    ) '  Title:              = "', trim ( title ) // '"'
    write ( *, '(a,g14.6)' ) '  H                   = ', h

    call file_name_inc ( eps_file_name )

    call p00_boundary_eps ( test, h, show_nodes, eps_file_name )

    write ( *, '(a)' ) '  Boundary image file is "' &
      // trim ( eps_file_name ) // '".'

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests P00_POINTS_EPS.
!
!  Discussion:
!
!    This routine creates an EPS file containing an image of
!    a set of points, and the boundary of the region.
!
!    The boundary is broken up into segments of about H = 0.05 in
!    relative length, and these points are shown in the image.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) eps_file_name
  real ( kind = 8 ) :: h
  real ( kind = 8 ) :: h_relative = 0.025D+00
  real ( kind = 8 ) hi(2)
  real ( kind = 8 ) lo(2)
  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 400
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  P00_POINTS_EPS makes an EPS file containing an'
  write ( *, '(a)' ) '  image of a set of points, and the boundary.'

  eps_file_name = 'p00_points.eps'

  call p00_test_num ( test_num )

  do test = 1, test_num

    call p00_box ( test, m, lo, hi )

    h = h_relative * max ( hi(1) - lo(1), hi(2) - lo(2) )

    seed = 123456789

    call p00_title ( test, title )

    write ( *, '(a)'       ) ' '
    write ( *, '(a,i6)'    ) '  Test number         =  ', test
    write ( *, '(a,i6)'    ) '  Spatial dimension   =  ', m
    write ( *, '(a,g14.6)' ) '  Relative spacing    =  ', h_relative
    write ( *, '(a,g14.6)' ) '  Actual spacing H    =  ', h
    write ( *, '(a, a)'    ) '  Title:              = "', trim ( title ) // '"'

    call file_name_inc ( eps_file_name )
    write ( *, * ) 'DEBUG :about to allocate'
    allocate ( point(1:m,1:n) )
    write ( *, * ) 'DEBUG :about to sample'
    call p00_sample ( test, m, n, seed, point )
    write ( *, * ) 'DEBUG :about to EPS'
    call p00_points_eps ( test, h, m, n, point, eps_file_name )

    write ( *, '(a)' ) '  Boundary image file is "' &
      // trim ( eps_file_name ) // '".'

    deallocate ( point )

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests P00_POLY_WRITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) file_name
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num

  file_name = 'p00.poly'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  P00_POLY_WRITE creates a POLY file out of the'
  write ( *, '(a)' ) '  boundary data for each problem.'

  call p00_test_num ( test_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of tests available = ', test_num

  do test = 1, test_num

    call file_name_inc ( file_name )

    call p00_poly_write ( test, file_name )

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests P00_HEX_GRID.
!
!  Discussion:
!
!    P00_HEX_GRID_COUNT and DTABLE_WRITE are also tested.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: file_name = 'p08_hex_grid_points.txt'
  real ( kind = 8 ) h
  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  P00_HEX_GRID returns points inside a region that'
  write ( *, '(a)' ) '    lie on a hexagonal grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  P00_HEX_GRID_COUNT counts the number of hex'
  write ( *, '(a) ') '    grid points inside a region.'
  write ( *, '(a)' ) '  DTABLE_WRITE writes sets of points to a file.'

  test = 8
  h = 0.0250D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  We use a hex grid spacing of H = ', h

  call p00_hex_grid_count ( test, m, h, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  P00_HEX_GRID_COUNT reports that the number of'
  write ( *, '(a,i6)' ) '  hex grid points will be ', n

  allocate ( point(1:m,1:n) )

  call p00_hex_grid ( test, m, h, n, point )

  call r8mat_write ( file_name, m, n, point )

  call r8mat_transpose_print_some ( m, n, point, 1, 1, m, 5, &
    '  A few of the points:' )

  deallocate ( point )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The hex grid points were written to the file "' // &
    trim ( file_name ) // '".'

  return
end
