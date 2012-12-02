program main

!*****************************************************************************80
!
!! MAIN creates data for a CVT "movie".
!
!  Discussion:
!
!    The computational region is TEST_TRIANGULATION region #8, the
!    "holey pie".
!
!    We initialize with hex grid points.
!
!    We include fixed points at the vertices of the region.
!
!    We sample on a "fattened" region, but force wayward centroids
!    back to the boundary.
!
!    We perform 100 CVT iterations, and save each set of generators
!    in a file, "p08_hbf_000.txt" contains the initial data, followed
!    by "p08_hbf_001.txt" through "p08_hbf_100.txt".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2005
!
!  Parameters:
!
!    Input, integer TEST, the index of the problem to be treated.
!
  implicit none

  real ( kind = 8 ) area_box
  real ( kind = 8 ) area_hex_unit
  real ( kind = 8 ) :: boundary_h = 0.05D+00
  logical comment
  integer, allocatable, dimension ( : ) :: count
  real ( kind = 8 ) energy
  character ( len = 80 ) file_eps_name
  character ( len = 80 ) file_txt_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: fixed
  integer fixed_num
  real ( kind = 8 ) h
  real ( kind = 8 ), allocatable, dimension ( : ) :: hi
  integer iteration
  integer, parameter :: iteration_max = 100
  integer j
  real ( kind = 8 ), allocatable, dimension ( : ) :: lo
  integer m
  integer n
  integer n_box
  integer, allocatable, dimension ( : ) :: nearest
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_new
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sample
  integer, parameter :: sample_num = 100000
  real ( kind = 8 ) scale
  integer seed
  integer :: seed_start = 123456789
  integer test
  character ( len = 80 ) title

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_MOVIE5_DATA:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Carry out a CVT iteration, saving the generators'
  write ( *, '(a)' ) '  at each step, to be used as data for a movie.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Apply CVT techniques to a test from TEST_TRIANGULATION'
  write ( *, '(a)' ) '  using:'
  write ( *, '(a)' ) '  * hex initialization;'
  write ( *, '(a)' ) '  * fixed points;'
  write ( *, '(a)' ) '  * sampling on an enlarged region;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Instead of specifying N, the number of points,'
  write ( *, '(a)' ) '  specify H, the desired average cell size.'

  test = 8

  call p00_fixed_num ( test, fixed_num )

  seed = seed_start
  energy = -1.0D+00

  call p00_title ( test, title )

  write ( *, '(a)' ) ' '
  write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

  call p00_header ( test )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of fixed points = ', fixed_num

  call p00_dimension ( test, m )

  allocate ( hi(1:m) )
  allocate ( lo(1:m) )

  call p00_box ( test, m, lo, hi )
!
!  Determine the value of N.
!
  area_box = ( hi(1) - lo(1) ) * ( hi(2) - lo(2) )
  area_hex_unit = 3.0D+00 * sqrt ( 3.0D+00 ) / 2.0D+00

  n_box = 400
  h = ( 2.0D+00 / sqrt ( 3.0D+00 ) ) &
    * sqrt ( area_box / ( real ( n_box, kind = 8 ) * area_hex_unit ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) &
    '  To get about N_BOX = ', n_box, '  points in the box,'
  write ( *, '(a,g14.6)' ) '  we request a hex grid spacing H = ', h

  call p00_hex_grid_count ( test, m, h, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  This actually gives us N = ', n
  write ( *, '(a)' ) '  points in the region, which may be less than N_BOX,'
  write ( *, '(a)' ) '  since the region may be a subset of the box.'

  n = n + fixed_num

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) '  We just throw in ', fixed_num, ' fixed points,'
  write ( *, '(a,i6)' ) '  Now this actually gives us N = ', n

  allocate ( point(1:m,1:n) )
  allocate ( fixed(1:m,1:fixed_num) )

  call p00_fixed_points ( test, m, fixed_num, fixed(1:m,1:fixed_num) )

  point(1:m,1:fixed_num) = fixed(1:m,1:fixed_num)

  call p00_hex_grid ( test, m, h, n, point(1:m,fixed_num+1:n) )

  call r8mat_transpose_print_some ( m, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )
!
!  Write the initial points to a file.
!
  file_txt_name = 'p08_hbf_000.txt'

  comment = .false.

  call cvt_write ( m, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )
!
!  Prepare for iteration.
!
  allocate ( count(1:n) )
  allocate ( nearest(1:sample_num) )
  allocate ( point_new(1:m,1:n) )
  allocate ( sample(1:m,1:sample_num) )

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
    h = 0.05D+00
    call p00_sample_h1 ( test, m, sample_num, h, seed, sample )
!
!  Find the nearest cell generator.
!
    call find_closest ( m, n, sample_num, point, sample, nearest )
!
!  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
!
    point_new(1:m,1:fixed_num) = fixed(1:m,1:fixed_num)
    point_new(1:m,fixed_num+1:n) = 0.0D+00
    count(1:n) = 0
    energy = 0.0D+00

    do j = 1, sample_num

      energy = energy + sum ( ( point(1:m,nearest(j)) - sample(1:m,j) )**2 )

      point_new(1:m,nearest(j)) = point_new(1:m,nearest(j)) + sample(1:m,j)

      count(nearest(j)) = count(nearest(j)) + 1

    end do
!
!  Compute the new generators.
!  But the fixed points don't move.
!
    do j = fixed_num+1, n
      if ( count(j) /= 0 ) then
        point_new(1:m,j) = point_new(1:m,j) / real ( count(j), kind = 8 )
      end if
    end do
!
!  Project generators back into region.
!
    call p00_boundary_project ( test, m, n-fixed_num, &
      point_new(1:m,fixed_num+1:n) )

    energy = energy / real ( sample_num, kind = 8 )
    write ( *, '(2x,g14.6)' ) energy
!
!  Update.
!
    point(1:m,fixed_num+1:n) = point_new(1:m,fixed_num+1:n)
!
!  Write to file.
!
    call file_name_inc ( file_txt_name )

    comment = .false.

    call cvt_write ( m, n, seed_start, seed, sample_num, iteration_max, &
      energy, point, file_txt_name, comment )

  end do

! if ( test < 10 ) then
!   write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '_hbf.txt'
! else
!   write ( file_txt_name, '(a,i2,a)' ) 'cvt_p', test, '_hbf.txt'
! end if

! call cvt_write ( m, n, seed_start, seed, sample_num, iteration_max, &
!   energy, point, file_txt_name, comment )

  if ( test < 10 ) then
    write ( file_eps_name, '(a,i1,a)' ) 'cvt_p0', test, '_hbf.eps'
  else
    write ( file_eps_name, '(a,i2,a)' ) 'cvt_p', test, '_hbf.eps'
  end if

  call p00_points_eps ( test, boundary_h, m, n, point, file_eps_name )

  deallocate ( count )
  deallocate ( fixed )
  deallocate ( hi )
  deallocate ( lo )
  deallocate ( nearest )
  deallocate ( point )
  deallocate ( point_new )
  deallocate ( sample )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_MOVIE5_DATA:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine find_closest ( m, n, sample_num, point, sample, nearest )

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
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer SAMPLE_NUM, the number of samples.
!
!    Input, real ( kind = 8 ) POINT(M,N), the point coordinates.
!
!    Input, real ( kind = 8 ) SAMPLE(M,SAMPLE_NUM), the sample coordinates.
!
!    Output, integer NEAREST(SAMPLE_NUM), the index of the nearest
!    point to each sample.
!
  implicit none

  integer m
  integer n
  integer sample_num

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  integer i
  integer j
  integer nearest(sample_num)
  real ( kind = 8 ) point(m,n)
  real ( kind = 8 ) sample(m,sample_num)

  do i = 1, sample_num

    dist_min = huge ( dist_min )
    nearest(i) = -1

    do j = 1, n

      dist = sum ( ( point(1:m,j) - sample(1:m,i) )**2 )

      if ( dist < dist_min ) then
        dist_min = dist
        nearest(i) = j
      end if

    end do

  end do

  return
end
subroutine cvt_write ( m, n, seed_start, seed, sample_num, iteration_max, &
  energy, point, file_out_name, comment )

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
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer SEED_START, the initial random number seed.
!
!    Input, integer SEED, the current random number seed.
!
!    Input, integer SAMPLE_NUM, the number of sampling points used on
!    each CVT iteration.
!
!    Input, integer STEP, the number of iterations used in the CVT
!    calculation.
!
!    Input, real ( kind = 8 ) ENERGY, the estimated Voronoi energy.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(M,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of
!    the output file.
!
!    Input, logical COMMENT, is true if comments may be included.
!
  implicit none

  integer m
  integer n

  logical comment
  real ( kind = 8 ) energy
  character ( len = * ) file_out_name
  integer file_out_unit
  integer ios
  integer iteration_max
  integer j
  real ( kind = 8 ) point(m,n)
  integer sample_num
  integer seed
  integer seed_start
  character ( len = 80 ) string

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

    write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
    write ( file_out_unit, '(a)'       ) &
      '#  created by routine CVT_WRITE in CVT_TRIANGULATION.F90'
    write ( file_out_unit, '(a)'       ) '#'
    write ( file_out_unit, '(a,i12)'   ) '#  Spatial dimension M = ', m
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

  write ( string, '(a,i3,a)' ) '(', m, '(f12.6,2x))'
  do j = 1, n
    write ( file_out_unit, string ) point(1:m,j)
  end do

  close ( unit = file_out_unit )

  return
end
