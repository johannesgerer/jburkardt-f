program main

!*****************************************************************************80
!
!! MAIN is the main program for CVT_TET_MESH.
!
!  Discussion:
!
!    CVT_TET_MESH uses CVT sampling on problems from TEST_TET_MESH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_TET_MESH:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Apply simple CVT sampling routines to produce'
  write ( *, '(a)' ) '  a set of sample points in regions from'
  write ( *, '(a)' ) '  the TEST_TET_MESH package.'
!
!  How many test cases are available?
!
  call p00_test_num ( test_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of test cases = ', test_num
!
!  Run test 1 on the test cases.  This simply computes a CVT in the
!  interior of the test region.
!
  if ( .false. ) then
    do test = 1, test_num
      call test01 ( test )
    end do
  end if
!
!  Run test 2 on the test cases.  This computes a nonstandard CVT in the
!  interior and on the boundary of the test region.
!
  if ( .false. ) then
    do test = 1, test_num
      call test02 ( test )
    end do
  end if
!
!  Run test 3 on the test cases.  This computes a nonstandard CVT in the
!  interior and on the boundary of the test region, plus fixed points.
!
  if ( .true. ) then
    do test = 2, test_num
      call test03 ( test )
    end do
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_TET_MESH:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( test )

!*****************************************************************************80
!
!! TEST01 creates a standard CVT in a given test region.
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
!    Input, integer ( kind = 4 ) TEST, the index of the problem to be treated.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  logical :: comment = .false.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: count
  real ( kind = 8 ) energy
  character ( len = 80 ) file_eps_name
  character ( len = 80 ) file_txt_name
  integer ( kind = 4 ) iteration
  integer ( kind = 4 ), parameter :: iteration_max = 40
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: n = 400
  integer ( kind = 4 ), allocatable, dimension ( : ) :: nearest
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_new
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ), parameter :: sample_num = 100000
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) :: seed_start = 123456789
  integer ( kind = 4 ) test
  character ( len = 80 ) title

  seed = seed_start
  energy = -1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Compute a standard CVT in the interior of'
  write ( *, '(a)' ) '  a test region from TEST_TET_MESH.'

  call p00_title ( test, title )

  write ( *, '(a)' ) ' '
  write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

  call p00_header ( test )
!
!  Initialize the sampling points.
!
  allocate ( count(1:n) )
  allocate ( nearest(1:sample_num) )
  allocate ( point(1:dim_num,1:n) )
  allocate ( point_new(1:dim_num,1:n) )
  allocate ( sample(1:dim_num,1:sample_num) )

  call p00_sample ( test, n, seed, point )

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated Voronoi energy:'
  write ( *, '(a)' ) ' '
!
!  Carry out an iteration
!
  do iteration = 1, iteration_max

    call p00_sample ( test, sample_num, seed, sample )
!
!  Find the nearest cell generator.
!
    call find_closest ( dim_num, n, sample_num, point, sample, nearest )
!
!  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
!
    point_new(1:dim_num,1:n) = 0.0D+00
    count(1:n) = 0
    energy = 0.0D+00

    do j = 1, sample_num

      energy = energy &
        + sum ( ( point(1:dim_num,nearest(j)) - sample(1:dim_num,j) )**2 )

      point_new(1:dim_num,nearest(j)) = point_new(1:dim_num,nearest(j)) &
        + sample(1:dim_num,j)

      count(nearest(j)) = count(nearest(j)) + 1

    end do
!
!  Compute the new generators.
!
    do j = 1, n
      if ( count(j) /= 0 ) then
        point_new(1:dim_num,j) = point_new(1:dim_num,j) &
          / real ( count(j), kind = 8 )
      end if
    end do

    energy = energy / real ( sample_num, kind = 8 )
    write ( *, '(2x,i6,2x,g14.6)' ) iteration, energy
!
!  Update.
!
    point(1:dim_num,1:n) = point_new(1:dim_num,1:n)

  end do

  if ( test < 10 ) then
    write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '_test01.txt'
  else
    write ( file_txt_name, '(a,i2,a)' ) 'cvt_p', test, '_test01.txt'
  end if

  call cvt_write ( dim_num, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST01: wrote data to "' &
    // trim ( file_txt_name ) // '".'

  deallocate ( count )
  deallocate ( nearest )
  deallocate ( point )
  deallocate ( point_new )
  deallocate ( sample )

  return
end
subroutine test02 ( test )

!*****************************************************************************80
!
!! TEST02 creates a CVT dataset, modified so some points are on the boundary.
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
!    Input, integer ( kind = 4 ) TEST, the index of the problem to be treated.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  logical :: comment = .false.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: count
  real ( kind = 8 ) energy
  character ( len = 80 ) file_txt_name
  real ( kind = 8 ) :: h = 0.05D+00
  real ( kind = 8 ) hi(dim_num)
  integer ( kind = 4 ) iteration
  integer ( kind = 4 ), parameter :: iteration_max = 40
  integer ( kind = 4 ) j
  real ( kind = 8 ) lo(dim_num)
  integer ( kind = 4 ), parameter :: n = 400
  integer ( kind = 4 ), allocatable, dimension ( : ) :: nearest
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_new
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ), parameter :: sample_num = 100000
  real ( kind = 8 ) scale
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) :: seed_start = 123456789
  integer ( kind = 4 ) test
  character ( len = 80 ) title

  seed = seed_start
  energy = -1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Try to compute a nonstandard CVT in the interior and'
  write ( *, '(a)' ) '  on the boundary of a test region from TEST_TET_MESH.'

  call p00_title ( test, title )

  write ( *, '(a)' ) ' '
  write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

  call p00_header ( test )
!
!  Allocate space for some arrays.
!
  allocate ( count(1:n) )
  allocate ( nearest(1:sample_num) )
  allocate ( point(1:dim_num,1:n) )
  allocate ( point_new(1:dim_num,1:n) )
  allocate ( sample(1:dim_num,1:sample_num) )
!
!  Determine the amount by which the region should be expanded.
!
  call p00_box ( test, lo, hi )
  scale = maxval ( hi(1:3) - lo(1:3) )
  h = 0.05D+00 * scale
!
!  Initialize the sampling points.
!
  call p00_sample ( test, n, seed, point )

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated Voronoi energy:'
  write ( *, '(a)' ) ' '
!
!  Carry out an iteration
!
  do iteration = 1, iteration_max
!
!  Get sample points from a "slightly enlarged" region.
!
    call p00_sample_h1 ( test, sample_num, h, seed, sample )
!
!  Find the nearest cell generator.
!
    call find_closest ( dim_num, n, sample_num, point, sample, nearest )
!
!  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
!
    point_new(1:dim_num,1:n) = 0.0D+00
    count(1:n) = 0
    energy = 0.0D+00

    do j = 1, sample_num

      energy = energy &
        + sum ( ( point(1:dim_num,nearest(j)) - sample(1:dim_num,j) )**2 )

      point_new(1:dim_num,nearest(j)) = point_new(1:dim_num,nearest(j)) &
        + sample(1:dim_num,j)

      count(nearest(j)) = count(nearest(j)) + 1

    end do
!
!  Compute the new generators.
!
    do j = 1, n
      if ( count(j) /= 0 ) then
        point_new(1:dim_num,j) = point_new(1:dim_num,j) &
          / real ( count(j), kind = 8 )
      end if
    end do
!
!  Project generators back into region.
!
    call p00_boundary_project ( test, n, point_new )

    energy = energy / real ( sample_num, kind = 8 )
    write ( *, '(2x,i6,2x,g14.6)' ) iteration, energy
!
!  Update.
!
    point(1:dim_num,1:n) = point_new(1:dim_num,1:n)

  end do

  if ( test < 10 ) then
    write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '_test02.txt'
  else
    write ( file_txt_name, '(a,i2,a)' ) 'cvt_p', test, '_test02.txt'
  end if

  call cvt_write ( dim_num, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST02: wrote data to "' &
    // trim ( file_txt_name ) // '".'

  deallocate ( count )
  deallocate ( nearest )
  deallocate ( point )
  deallocate ( point_new )
  deallocate ( sample )

  return
end
subroutine test03 ( test )

!*****************************************************************************80
!
!! TEST03 creates a CVT dataset with points on the boundary or set by user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the index of the problem to be treated.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) box_volume
  logical :: comment = .false.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: count
  real ( kind = 8 ) energy
  character ( len = 80 ) file_txt_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: fixed
  integer ( kind = 4 ) fixed_num
  real ( kind = 8 ) h
  real ( kind = 8 ) hi(dim_num)
  integer ( kind = 4 ) iteration
  integer ( kind = 4 ), parameter :: iteration_max = 20
  integer ( kind = 4 ) j
  real ( kind = 8 ) lo(dim_num)
  integer ( kind = 4 ), parameter :: n = 400
  integer ( kind = 4 ), allocatable, dimension ( : ) :: nearest
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_new
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer ( kind = 4 ), parameter :: sample_num = 1000000
  real ( kind = 8 ) scale
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) :: seed_start = 123456789
  integer ( kind = 4 ) test
  character ( len = 80 ) title

  call p00_fixed_num ( test, fixed_num )

  seed = seed_start
  energy = -1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Try to compute a nonstandard CVT in the interior and'
  write ( *, '(a)' ) '  on the boundary of a test region from TEST_TET_MESH.'

  call p00_title ( test, title )

  write ( *, '(a)' ) ' '
  write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

  call p00_header ( test )
!
!  Allocate space for some arrays.
!
  allocate ( count(1:n) )

  if ( 0 < fixed_num ) then
    allocate ( fixed(1:dim_num,1:fixed_num) )
  else
  end if

  allocate ( nearest(1:sample_num) )
  allocate ( point(1:dim_num,1:n) )
  allocate ( point_new(1:dim_num,1:n) )
  allocate ( sample(1:dim_num,1:sample_num) )
!
!  Determine the amount by which the region should be expanded.
!
  call p00_box ( test, lo, hi )
  box_volume = product ( hi(1:3) - lo(1:3) )
  h = ( box_volume / real ( n, kind = 8 ) )**( 1.0D+00 / 3.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Number of points requested = ', n
  write ( *, '(a,i12)' ) '  Number of CVT iterations = ', iteration_max
  write ( *, '(a,i12)' ) '  Number of CVT sample points = ', sample_num
  write ( *, '(a,g14.6)' ) '  Using expansion increment H = ', h
!
!  Get the fixed points.
!
  if ( 0 < fixed_num ) then
    call p00_fixed_points ( test, fixed_num, fixed(1:dim_num,1:fixed_num) )
  end if
!
!  Initialize the sampling points.
!
  if ( 0 < fixed_num ) then
    point(1:dim_num,1:fixed_num) = fixed(1:dim_num,1:fixed_num)
  end if

  call p00_sample ( test, n-fixed_num, seed, point(1:dim_num,fixed_num+1:n) )

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated Voronoi energy:'
  write ( *, '(a)' ) ' '
!
!  Carry out an iteration
!
  do iteration = 1, iteration_max
!
!  Get sample points from a "slightly enlarged" region.
!
    write ( *, * ) 'DEBUG1'
    call p00_sample_h1 ( test, sample_num, h, seed, sample )
!
!  Find the nearest cell generator.
!
    write ( *, * ) 'DEBUG2'
    call find_closest ( dim_num, n, sample_num, point, sample, nearest )
!
!  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
!
    write ( *, * ) 'DEBUG3'
    point_new(1:dim_num,1:n) = 0.0D+00
    count(1:n) = 0
    energy = 0.0D+00

    do j = 1, sample_num

      energy = energy &
        + sum ( ( point(1:dim_num,nearest(j)) - sample(1:dim_num,j) )**2 )

      point_new(1:dim_num,nearest(j)) = point_new(1:dim_num,nearest(j)) &
        + sample(1:dim_num,j)

      count(nearest(j)) = count(nearest(j)) + 1

    end do
!
!  Compute the new generators.
!  But the fixed points don't move.
!
    write ( *, * ) 'DEBUG4'
    do j = fixed_num+1, n
      if ( count(j) /= 0 ) then
        point_new(1:dim_num,j) = point_new(1:dim_num,j) &
          / real ( count(j), kind = 8 )
      end if
    end do
!
!  Project generators back into region.
!
    write ( *, * ) 'DEBUG5'
    call p00_boundary_project ( test, n, point_new )

    energy = energy / real ( sample_num, kind = 8 )
    write ( *, '(2x,i6,2x,g14.6)' ) iteration, energy
!
!  Update.
!
    point(1:dim_num,fixed_num+1:n) = point_new(1:dim_num,fixed_num+1:n)

  end do

  if ( test < 10 ) then
    write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '_test03.txt'
  else
    write ( file_txt_name, '(a,i2,a)' ) 'cvt_p', test, '_test03.txt'
  end if

  call cvt_write ( dim_num, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST03: wrote data to "' &
    // trim ( file_txt_name ) // '".'

  deallocate ( count )
  if ( 0 < fixed_num ) then
    deallocate ( fixed )
  end if
  deallocate ( nearest )
  deallocate ( point )
  deallocate ( point_new )
  deallocate ( sample )

  return
end
subroutine cvt_write ( dim_num, n, seed_start, seed, sample_num, &
  iteration_max, energy, point, file_out_name, comment )

!*****************************************************************************80
!
!! CVT_WRITE writes a CVT dataset to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the M-dimensional
!    components of a CVT generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) SEED_START, the initial random number seed.
!
!    Input, integer ( kind = 4 ) SEED, the current random number seed.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of sampling points used on
!    each CVT iteration.
!
!    Input, integer ( kind = 4 ) STEP, the number of iterations used in the CVT
!    calculation.
!
!    Input, real ( kind = 8 ) ENERGY, the estimated Voronoi energy.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of
!    the output file.
!
!    Input, logical COMMENT, is true if comments may be included.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  logical comment
  real ( kind = 8 ) energy
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iteration_max
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(dim_num,n)
  integer ( kind = 4 ) sample_num
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_start
  character ( len = 40 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, &
    status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(a)' ) '  "' // trim ( file_out_name ) //'".'
    stop
  end if

  if ( comment ) then

    call timestring ( string )

    write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
    write ( file_out_unit, '(a)'       ) &
      '#  created by routine CVT_WRITE in CVT_TETRA.F90'
    write ( file_out_unit, '(a)'       ) '#  at ' // trim ( string )
    write ( file_out_unit, '(a)'       ) '#'
    write ( file_out_unit, '(a,i12)'   ) &
      '#  Spatial dimension DIM_NUM = ', dim_num
    write ( file_out_unit, '(a,i12)'   ) '#  Number of points N =  ', n
    write ( file_out_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) = ', &
      epsilon ( point(1,1) )
    write ( file_out_unit, '(a,i12)'   ) '#  Initial SEED = ', seed_start
    write ( file_out_unit, '(a,i12)'   ) '#  Current SEED = ', seed
    write ( file_out_unit, '(a,i12)'   ) &
      '#  Number of sample points =       ', sample_num
    write ( file_out_unit, '(a,i12)'   ) &
      '#  Number of sampling iterations = ', iteration_max
    write ( file_out_unit, '(a,g14.6)' ) &
      '#  Estimated Voronoi energy =      ', energy
    write ( file_out_unit, '(a)'       )  '#'

  end if

  write ( string, '(a,i3,a)' ) '(', dim_num, '(f12.6,2x))'
  do j = 1, n
    write ( file_out_unit, string ) point(1:dim_num,j)
  end do

  close ( unit = file_out_unit )

  return
end
subroutine find_closest ( dim_num, n, sample_num, point, sample, nearest )

!*****************************************************************************80
!
!! FIND_CLOSEST finds the closest point to each sample.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) SAMPLE_NUM, the number of samples.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,N), the point coordinates.
!
!    Input, real ( kind = 8 ) SAMPLE(DIM_NUM,SAMPLE_NUM), the sample
!    coordinates.
!
!    Output, integer ( kind = 4 ) NEAREST(SAMPLE_NUM), the index of the nearest
!    point to each sample.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n
  integer ( kind = 4 ) sample_num

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest(sample_num)
  real ( kind = 8 ) point(dim_num,n)
  real ( kind = 8 ) sample(dim_num,sample_num)

  do i = 1, sample_num

    dist_min = huge ( dist_min )
    nearest(i) = -1

    do j = 1, n

      dist = sum ( ( point(1:dim_num,j) - sample(1:dim_num,i) )**2 )

      if ( dist < dist_min ) then
        dist_min = dist
        nearest(i) = j
      end if

    end do

  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = '31 May 2001   9:45:54.872 AM'
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
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
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
  character ( len = * ) string
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

  write ( string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
