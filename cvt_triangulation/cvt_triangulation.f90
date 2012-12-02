program main

!*****************************************************************************80
!
!! MAIN is the main program for CVT_TRIANGULATION.
!
!  Discussion:
!
!    CVT_TRIANGULATION uses CVT sampling on problems from TEST_TRIANGULATION.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
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
  write ( *, '(a)' ) 'CVT_TRIANGULATION:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Apply simple CVT sampling routines to produce'
  write ( *, '(a)' ) '  a set of sample points in regions from'
  write ( *, '(a)' ) '  the TEST_TRIANGULATION package.'

  call p00_test_num ( test_num )

! do test = 1, test_num
!   call test01 ( test )
! end do

! do test = 1, test_num
!   call test02 ( test )
! end do

! do test = 1, test_num
!   call test03 ( test )
! end do

  test = 10

! do test = 1, test_num
    call test04 ( test )
! end do

! test = 10
! do test = 1, test_num
!   call test05 ( test )
! end do

! call test06
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_TRIANGULATION:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( test )

!*****************************************************************************80
!
!! TEST01 uses CVT techniques on a given test problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the problem to be treated.
!
  implicit none

  logical :: comment = .true.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: count
  integer ( kind = 4 ), parameter :: dim_num = 2
  real ( kind = 8 ) energy
  character ( len = 80 ) file_eps_name
  character ( len = 80 ) file_txt_name
  real ( kind = 8 ) :: boundary_h = 0.05D+00
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
  write ( *, '(a)' ) '  Apply CVT techniques to a test problem'
  write ( *, '(a)' ) '  from TEST_TRIANGULATION.'

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

  call p00_sample ( test, dim_num, n, seed, point )

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated Voronoi energy:'
  write ( *, '(a)' ) ' '
!
!  Carry out an iteration
!
  do iteration = 1, iteration_max

    call p00_sample ( test, dim_num, sample_num, seed, sample )
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
    write ( *, '(2x,g14.6)' ) energy
!
!  Update.
!
    point(1:dim_num,1:n) = point_new(1:dim_num,1:n)

  end do

  if ( test < 10 ) then
    write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '.txt'
  else
    write ( file_txt_name, '(a,i2,a)' ) 'cvt_p', test, '.txt'
  end if

  call cvt_write ( dim_num, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )

  if ( test < 10 ) then
    write ( file_eps_name, '(a,i1,a)' ) 'cvt_p0', test, '.eps'
  else
    write ( file_eps_name, '(a,i2,a)' ) 'cvt_p', test, '.eps'
  end if

  call p00_points_eps ( test, boundary_h, dim_num, n, point, file_eps_name )

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
!! TEST02 reruns TEST01 cases with fixed points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the problem to be treated.
!
  implicit none

  logical :: comment = .true.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: count
  integer ( kind = 4 ), parameter :: dim_num = 2
  real ( kind = 8 ) energy
  character ( len = 80 ) file_txt_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: fixed
  integer ( kind = 4 ) fixed_num
  integer ( kind = 4 ) iteration
  integer ( kind = 4 ), parameter :: iteration_max = 20
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

  call p00_fixed_num ( test, fixed_num )

  if ( fixed_num <= 0 ) then
    return
  end if

  seed = seed_start
  energy = -1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Apply CVT techniques to a test problem'
  write ( *, '(a)' ) '  from TEST_TRIANGULATION.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we only rerun cases which involve'
  write ( *, '(a)' ) '  fixed points, and show how to handle them.'
 
  call p00_title ( test, title )

  write ( *, '(a)' ) ' '
  write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

  call p00_header ( test )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of fixed points = ', fixed_num
!
!  Initialize the sampling points.
!
  allocate ( count(1:n) )
  allocate ( fixed(1:dim_num,1:fixed_num) )
  allocate ( nearest(1:sample_num) )
  allocate ( point(1:dim_num,1:n) )
  allocate ( point_new(1:dim_num,1:n) )
  allocate ( sample(1:dim_num,1:sample_num) )

  call p00_fixed_points ( test, dim_num, fixed_num, &
    fixed(1:dim_num,1:fixed_num) )

  point(1:dim_num,1:fixed_num) = fixed(1:dim_num,1:fixed_num)

  call p00_sample ( test, dim_num, n-fixed_num, seed, &
    point(1:dim_num,fixed_num+1:n) )

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated Voronoi energy:'
  write ( *, '(a)' ) ' '
!
!  Carry out an iteration
!
  do iteration = 1, iteration_max

    call p00_sample ( test, dim_num, sample_num, seed, sample )
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
!  But the fixed points don't move.
!
    do j = fixed_num+1, n
      if ( count(j) /= 0 ) then
        point_new(1:dim_num,j) = point_new(1:dim_num,j) &
          / real ( count(j), kind = 8 )
      end if
    end do

    energy = energy / real ( sample_num, kind = 8 )
    write ( *, '(2x,g14.6)' ) energy
!
!  Update.
!
    point(1:dim_num,fixed_num+1:n) = point_new(1:dim_num,fixed_num+1:n)

  end do

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Final points (first 10 only)' )

  if ( test < 10 ) then
    write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '_fixed.txt'
  else
    write ( file_txt_name, '(a,i2,a)' ) 'cvt_p', test, '_fixed.txt'
  end if

  call cvt_write ( dim_num, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )

  deallocate ( count )
  deallocate ( fixed )
  deallocate ( nearest )
  deallocate ( point )
  deallocate ( point_new ) 
  deallocate ( sample )

  return
end
subroutine test03 ( test )

!*****************************************************************************80
!
!! TEST03 uses constrained CVT techniques on a given test problem.
!
!  Discussion:
!
!    Using sampling on the region after it has been enlarged by H,
!    and projection of generators that are exterior to the original
!    region back onto the boundary, this routine tries to get a
!    mesh that is CVT in the interior, and that has many generators
!    on the boundary (after they have been pushed there!).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the problem to be treated.
!
  implicit none

  logical :: comment = .true.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: count
  integer ( kind = 4 ), parameter :: dim_num = 2
  real ( kind = 8 ) energy
  character ( len = 80 ) file_eps_name
  character ( len = 80 ) file_txt_name
  real ( kind = 8 ) :: h
  real ( kind = 8 ), allocatable, dimension ( : ) :: hi
  integer ( kind = 4 ) iteration
  integer ( kind = 4 ), parameter :: iteration_max = 40
  integer ( kind = 4 ) j
  real ( kind = 8 ), allocatable, dimension ( : ) :: lo
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
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Try to get an approximate CVT mesh in a region'
  write ( *, '(a)' ) '  that also has many points ON the boundary.'

  call p00_title ( test, title )

  write ( *, '(a)' ) ' '
  write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

  call p00_header ( test )
!
!  Initialize the sampling points.
!
  allocate ( count(1:n) )
  allocate ( hi(1:dim_num) )
  allocate ( lo(1:dim_num) )
  allocate ( nearest(1:sample_num) )
  allocate ( point(1:dim_num,1:n) )
  allocate ( point_new(1:dim_num,1:n) )
  allocate ( sample(1:dim_num,1:sample_num) )
!
!  Determine the amount by which the region should be expanded.
!
  call p00_box ( test, dim_num, lo, hi )
  scale = maxval ( hi(1:dim_num) - lo(1:dim_num) )
  h = 0.05D+00 * scale
!
!  Get starting values.
!
  call p00_sample ( test, dim_num, n, seed, point )

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated Voronoi energy (before projection):'
  write ( *, '(a)' ) ' '
!
!  Carry out an iteration
!
  do iteration = 1, iteration_max
!
!  Get sample points from a "slightly enlarged" region.
!
    call p00_sample_h1 ( test, dim_num, sample_num, h, seed, sample )
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
!  TEMPORARY:
!
    if ( iteration == 1 ) then
      file_txt_name = 'centroids1.txt'
      call cvt_write ( dim_num, n, seed_start, seed, sample_num, &
        iteration_max, energy, point_new, file_txt_name, comment )
    end if
!
!  Project generators back into region.
!
    call p00_boundary_project ( test, dim_num, n, point_new )
!
!  TEMPORARY:
!
    if ( iteration == 1 ) then
      file_txt_name = 'centroids2.txt'
      call cvt_write ( dim_num, n, seed_start, seed, sample_num, &
        iteration_max, energy, point_new, file_txt_name, comment )
    end if

    energy = energy / real ( sample_num, kind = 8 )
    write ( *, '(2x,i6,2x,g14.6)' ) iteration, energy
!
!  Update.
!
    point(1:dim_num,1:n) = point_new(1:dim_num,1:n)

  end do

  if ( test < 10 ) then
    write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '_boundary.txt'
  else
    write ( file_txt_name, '(a,i2,a)' ) 'cvt_p', test, '_boundary.txt'
  end if

  call cvt_write ( dim_num, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )

  if ( test < 10 ) then
    write ( file_eps_name, '(a,i1,a)' ) 'cvt_p0', test, '_boundary.eps'
  else
    write ( file_eps_name, '(a,i2,a)' ) 'cvt_p', test, '_boundary.eps'
  end if

  call p00_points_eps ( test, h, dim_num, n, point, file_eps_name )

  deallocate ( count )
  deallocate ( hi )
  deallocate ( lo )
  deallocate ( nearest )
  deallocate ( point )
  deallocate ( point_new ) 
  deallocate ( sample )

  return
end
subroutine test04 ( test )

!*****************************************************************************80
!
!! TEST04 repeats test 3 with fixed points.
!
!  Discussion:
!
!    Using sampling on the region after it has been enlarged by H,
!    and projection of generators that are exterior to the original
!    region back onto the boundary, this routine tries to get a
!    mesh that is CVT in the interior, and that has many generators
!    on the boundary (after they have been pushed there!).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the problem to be treated.
!
  implicit none

  logical :: comment = .true.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: count
  integer ( kind = 4 ), parameter :: dim_num = 2
  real ( kind = 8 ) energy
  character ( len = 80 ) file_eps_name
  character ( len = 80 ) file_txt_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: fixed
  integer ( kind = 4 ) fixed_num
  real ( kind = 8 ) :: h
  real ( kind = 8 ), allocatable, dimension ( : ) :: hi
  integer ( kind = 4 ) iteration
  integer ( kind = 4 ), parameter :: iteration_max = 20
  integer ( kind = 4 ) j
  real ( kind = 8 ), allocatable, dimension ( : ) :: lo
  integer ( kind = 4 ), parameter :: n = 100
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

  if ( fixed_num <= 0 ) then
    return
  end if

  seed = seed_start
  energy = -1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Try to get an approximate CVT mesh in a region'
  write ( *, '(a)' ) '  that also has many points ON the boundary.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we only rerun cases which involve'
  write ( *, '(a)' ) '  fixed points, and show how to handle them.'

  call p00_title ( test, title )

  write ( *, '(a)' ) ' '
  write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

  call p00_header ( test )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of fixed points = ', fixed_num
!
!  Initialize the sampling points.
!
  allocate ( count(1:n) )
  allocate ( fixed(1:dim_num,1:fixed_num) )
  allocate ( hi(1:dim_num) )
  allocate ( lo(1:dim_num) )
  allocate ( nearest(1:sample_num) )
  allocate ( point(1:dim_num,1:n) )
  allocate ( point_new(1:dim_num,1:n) )
  allocate ( sample(1:dim_num,1:sample_num) )
!
!  Determine the amount by which the region should be expanded.
!
  call p00_box ( test, dim_num, lo, hi )
  scale = maxval ( hi(1:dim_num) - lo(1:dim_num) )
  h = 0.10D+00 * scale
!
!  Get starting values.
!
  call p00_fixed_points ( test, dim_num, fixed_num, &
    fixed(1:dim_num,1:fixed_num) )

  point(1:dim_num,1:fixed_num) = fixed(1:dim_num,1:fixed_num)

  call p00_sample ( test, dim_num, n-fixed_num, seed, &
    point(1:dim_num,fixed_num+1:n) )

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated Voronoi energy (before projection):'
  write ( *, '(a)' ) ' '
!
!  Carry out an iteration
!
  do iteration = 1, iteration_max
!
!  Get sample points from a "slightly enlarged" region.
!
    call p00_sample_h1 ( test, dim_num, sample_num, h, seed, sample )
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
!  But the fixed points don't move.
!
    do j = fixed_num+1, n
      if ( count(j) /= 0 ) then
        point_new(1:dim_num,j) = point_new(1:dim_num,j) &
          / real ( count(j), kind = 8 )
      end if
    end do
!
!  Project generators back into region.
!
    call p00_boundary_project ( test, dim_num, n-fixed_num, &
      point_new(1:dim_num,fixed_num+1:n) )

    energy = energy / real ( sample_num, kind = 8 )
    write ( *, '(2x,i6,2x,g14.6)' ) iteration, energy
!
!  Update.
!
    point(1:dim_num,fixed_num+1:n) = point_new(1:dim_num,fixed_num+1:n)

  end do

  if ( test < 10 ) then
    write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '_boundary_fixed.txt'
  else
    write ( file_txt_name, '(a,i2,a)' ) 'cvt_p', test, '_boundary_fixed.txt'
  end if

  write ( *, '(a)' ) '  Creating data file "' // trim ( file_txt_name ) // '".'

  call cvt_write ( dim_num, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )

  if ( test < 10 ) then
    write ( file_eps_name, '(a,i1,a)' ) 'cvt_p0', test, '_boundary_fixed.eps'
  else
    write ( file_eps_name, '(a,i2,a)' ) 'cvt_p', test, '_boundary_fixed.eps'
  end if

  write ( *, '(a)' ) '  Creating graphics file "' // &
    trim ( file_eps_name ) // '".'

  call p00_points_eps ( test, h, dim_num, n, point, file_eps_name )

  deallocate ( count )
  deallocate ( fixed )
  deallocate ( hi )
  deallocate ( lo )
  deallocate ( nearest )
  deallocate ( point )
  deallocate ( point_new ) 
  deallocate ( sample )

  return
end
subroutine test05 ( test )

!*****************************************************************************80
!
!! TEST05 uses CVT techniques with hex initialization.
!
!  Discussion:
!
!    This test is similar to test #1, but instead of specifying
!    a fixed number of points, we specify a desired average cell size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
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
  integer ( kind = 4 ), allocatable, dimension ( : ) :: count
  integer ( kind = 4 ), parameter :: dim_num = 2
  real ( kind = 8 ) energy
  character ( len = 80 ) file_eps_name
  character ( len = 80 ) file_txt_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: fixed
  integer ( kind = 4 ) fixed_num
  real ( kind = 8 ) h
  real ( kind = 8 ), allocatable, dimension ( : ) :: hi
  integer ( kind = 4 ) iteration
  integer ( kind = 4 ), parameter :: iteration_max = 100
  integer ( kind = 4 ) j
  real ( kind = 8 ), allocatable, dimension ( : ) :: lo
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_box
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

  call p00_fixed_num ( test, fixed_num )

  seed = seed_start
  energy = -1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Apply CVT techniques to a test from TEST_TRIANGULATION'
  write ( *, '(a)' ) '  using:'
  write ( *, '(a)' ) '  * hex initialization;'
  write ( *, '(a)' ) '  * fixed points;'
  write ( *, '(a)' ) '  * sampling on an enlarged region;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Instead of specifying N, the number of points,'
  write ( *, '(a)' ) '  specify H, the desired average cell size.'
  
  call p00_title ( test, title )

  write ( *, '(a)' ) ' '
  write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

  call p00_header ( test )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of fixed points = ', fixed_num

  allocate ( hi(1:dim_num) )
  allocate ( lo(1:dim_num) )

  call p00_box ( test, dim_num, lo, hi )
!
!  Determine the value of N.
!
  area_box = product ( hi(1:dim_num) - lo(1:dim_num) )
  area_hex_unit = 3.0D+00 * sqrt ( 3.0D+00 ) / 2.0D+00

  n_box = 400
  h = ( 2.0D+00 / sqrt ( 3.0D+00 ) ) &
    * sqrt ( area_box / ( real ( n_box, kind = 8 ) * area_hex_unit ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) &
    '  To get about N_BOX = ', n_box, '  points in the box,'
  write ( *, '(a,g14.6)' ) '  we request a hex grid spacing H = ', h
  
  call p00_hex_grid_count ( test, dim_num, h, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  This actually gives us N = ', n
  write ( *, '(a)' ) '  points in the region, which may be less than N_BOX,'
  write ( *, '(a)' ) '  since the region may be a subset of the box.'

  n = n + fixed_num

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) '  We just throw in ', fixed_num, ' fixed points,'
  write ( *, '(a,i6)' ) '  Now this actually gives us N = ', n

  allocate ( point(1:dim_num,1:n) )
  allocate ( fixed(1:dim_num,1:fixed_num) )

  call p00_fixed_points ( test, dim_num, fixed_num, fixed(1:dim_num,1:fixed_num) )

  point(1:dim_num,1:fixed_num) = fixed(1:dim_num,1:fixed_num)

  call p00_hex_grid ( test, dim_num, h, n, point(1:dim_num,fixed_num+1:n) )

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )
!
!  Prepare for iteration.
!
  allocate ( count(1:n) )
  allocate ( nearest(1:sample_num) )
  allocate ( point_new(1:dim_num,1:n) )
  allocate ( sample(1:dim_num,1:sample_num) )

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
    call p00_sample_h1 ( test, dim_num, sample_num, h, seed, sample )
!
!  Find the nearest cell generator.
!
    call find_closest ( dim_num, n, sample_num, point, sample, nearest )
!
!  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
!
    point_new(1:dim_num,1:fixed_num) = fixed(1:dim_num,1:fixed_num)
    point_new(1:dim_num,fixed_num+1:n) = 0.0D+00
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
    do j = fixed_num+1, n
      if ( count(j) /= 0 ) then
        point_new(1:dim_num,j) = point_new(1:dim_num,j) &
          / real ( count(j), kind = 8 )
      end if
    end do
!
!  Project generators back into region.
!
    call p00_boundary_project ( test, dim_num, n-fixed_num, &
      point_new(1:dim_num,fixed_num+1:n) )

    energy = energy / real ( sample_num, kind = 8 )
    write ( *, '(2x,g14.6)' ) energy
!
!  Update.
!
    point(1:dim_num,fixed_num+1:n) = point_new(1:dim_num,fixed_num+1:n)

  end do
!
!  Write data to file.
!
  if ( test < 10 ) then
    write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '_hbf.txt'
  else
    write ( file_txt_name, '(a,i2,a)' ) 'cvt_p', test, '_hbf.txt'
  end if

  call cvt_write ( dim_num, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )

  write ( *, '(a)' ) '  TEST05: Final data written to "' &
    // trim ( file_txt_name ) // '".'
!
!  Write EPS image to file.
!
  if ( test < 10 ) then
    write ( file_eps_name, '(a,i1,a)' ) 'cvt_p0', test, '_hbf.eps'
  else
    write ( file_eps_name, '(a,i2,a)' ) 'cvt_p', test, '_hbf.eps'
  end if

  call p00_points_eps ( test, boundary_h, dim_num, n, point, file_eps_name )

  write ( *, '(a)' ) '  TEST05: EPS image written to "' &
    // trim ( file_eps_name ) // '".'

  deallocate ( count )
  deallocate ( fixed )
  deallocate ( hi )
  deallocate ( lo )
  deallocate ( nearest )
  deallocate ( point )
  deallocate ( point_new ) 
  deallocate ( sample )

  return
end
subroutine test06

!*****************************************************************************80
!
!! TEST06 uses CVT techniques with a nonuniform density.
!
!  Discussion:
!
!    This test is similar to test #5, but a user-specified sampling
!    routine is called for sampling test region 10, the unit square.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area_box
  real ( kind = 8 ) :: boundary_h = 0.05D+00
  logical :: comment = .true.
  integer ( kind = 4 ), allocatable, dimension ( : ) :: count
  integer ( kind = 4 ), parameter :: dim_num = 2
  real ( kind = 8 ) energy
  character ( len = 80 ) file_eps_name
  character ( len = 80 ) file_txt_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: fixed
  integer ( kind = 4 ) fixed_num
  real ( kind = 8 ) h
  integer ( kind = 4 ) iteration
  integer ( kind = 4 ), parameter :: iteration_max = 40
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: n = 400
  integer ( kind = 4 ) n_box
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
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
!
!  We will work on test 10, the unit square.
!
  test = 10

  call p00_fixed_num ( test, fixed_num )

  seed = seed_start
  energy = -1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Apply CVT techniques to a test from TEST_TRIANGULATION'
  write ( *, '(a)' ) '  using:'
  write ( *, '(a)' ) '  * hex initialization;'
  write ( *, '(a)' ) '  * fixed points;'
  write ( *, '(a)' ) '  * sampling on an enlarged region;'
  write ( *, '(a)' ) '  * a nonuniform density is applied.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Instead of specifying N, the number of points,'
  write ( *, '(a)' ) '  specify H, the desired average cell size.'
  
  call p00_title ( test, title )

  write ( *, '(a)' ) ' '
  write ( *, '(a, a)' ) '  Title:             "', trim ( title ) // '"'

  call p00_header ( test )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of fixed points = ', fixed_num

  write ( *, '(a,i6)' ) '  Number of generators N = ', n

  allocate ( point(1:dim_num,1:n) )
  allocate ( fixed(1:dim_num,1:fixed_num) )
  allocate ( count(1:n) )
  allocate ( nearest(1:sample_num) )
  allocate ( point_new(1:dim_num,1:n) )
  allocate ( sample(1:dim_num,1:sample_num) )
  allocate ( weight(1:sample_num) )

  if ( .false. ) then

    call p00_fixed_points ( test, dim_num, fixed_num, &
      fixed(1:dim_num,1:fixed_num) )

    point(1:dim_num,1:fixed_num) = fixed(1:dim_num,1:fixed_num)

  else

    fixed_num = 0

  end if

  h = 0.0D+00
  call p99_sample_h1 ( dim_num, n-fixed_num, h, seed, &
    point(1:dim_num,fixed_num+1:n) )

  call r8mat_transpose_print_some ( dim_num, n, point, &
    1, 1, 2, 10, '  Initial points (first 10 only)' )
!
!  Prepare for iteration.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated Weighted Voronoi energy:'
  write ( *, '(a)' ) ' '
!
!  Carry out an iteration
!
  do iteration = 1, iteration_max
!
!  Get sample points from a "slightly enlarged" region.
!
    if ( .false. ) then
      h = 0.10D+00
      call p99_sample_h1 ( dim_num, sample_num, h, seed, sample )
    else
      h = 0.0D+00
      call p99_sample_h1 ( dim_num, sample_num, h, seed, sample )
    end if
!
!  Find the nearest cell generator.
!
    call find_closest ( dim_num, n, sample_num, point, sample, nearest )
!
!  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
!
    point_new(1:dim_num,1:fixed_num) = fixed(1:dim_num,1:fixed_num)
    point_new(1:dim_num,fixed_num+1:n) = 0.0D+00

    count(1:n) = 0

    call p99_weight ( dim_num, n, sample, weight )

    energy = 0.0D+00

    do j = 1, sample_num

      energy = energy &
        + weight(j) &
        * sum ( ( point(1:dim_num,nearest(j)) - sample(1:dim_num,j) )**2 )

      point_new(1:dim_num,nearest(j)) = point_new(1:dim_num,nearest(j)) &
        + sample(1:dim_num,j)

      count(nearest(j)) = count(nearest(j)) + 1

    end do
!
!  Compute the new generators.
!  But the fixed points don't move.
!
    do j = fixed_num+1, n
      if ( count(j) /= 0 ) then
        point_new(1:dim_num,j) = point_new(1:dim_num,j) &
          / real ( count(j), kind = 8 )
      end if
    end do
!
!  Project generators back into region.
!
    if ( .false. ) then
      call p00_boundary_project ( test, dim_num, n-fixed_num, &
        point_new(1:dim_num,fixed_num+1:n) )
    end if

    energy = energy / real ( sample_num, kind = 8 )
    write ( *, '(2x,g14.6)' ) energy
!
!  Update.
!
    point(1:dim_num,fixed_num+1:n) = point_new(1:dim_num,fixed_num+1:n)

  end do

! if ( test < 10 ) then
!   write ( file_txt_name, '(a,i1,a)' ) 'cvt_p0', test, '_weighted.txt'
! else
!   write ( file_txt_name, '(a,i2,a)' ) 'cvt_p99_weighted.txt'
! end if

  file_txt_name = 'cvt_p98_weighted.txt'

  call cvt_write ( dim_num, n, seed_start, seed, sample_num, iteration_max, &
    energy, point, file_txt_name, comment )

! if ( test < 10 ) then
!   write ( file_eps_name, '(a,i1,a)' ) 'cvt_p0', test, '_weighted.eps'
! else
!   write ( file_eps_name, '(a,i2,a)' ) 'cvt_p99_weighted.eps'
! end if

  file_eps_name = 'cvt_p98_weighted.eps'

  call p00_points_eps ( test, boundary_h, dim_num, n, point, file_eps_name )

  deallocate ( count )
  deallocate ( fixed )
  deallocate ( nearest )
  deallocate ( point )
  deallocate ( point_new ) 
  deallocate ( sample )
  deallocate ( weight )

  return
end
subroutine p99_sample_h1 ( dim_num, n, h, seed, point )

!*****************************************************************************80
!
!! P99_SAMPLE_H1 samples points from the enlarged region in problem 99.
!
!  Discussion:
!
!    The region is enlarged by an amount H in each direction,
!    and a nonuniform density is applied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real H, the enlargement amount.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) POINT(DIM_NUM,N), the coordinates
!    of the points.
!
  implicit none

  integer dim_num
  integer n
  integer, parameter :: batch = 100

  integer ( kind = 4 ) found
  real ( kind = 8 ) h
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(dim_num,n)
  real ( kind = 8 ) point2(dim_num,batch)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) test(batch)
  real ( kind = 8 ) weight(batch)
  real ( kind = 8 ), parameter :: x1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 = +1.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 = +1.0D+00

  found = 0

  do while ( found < n )
!
!  Generate a batch of points in [0,1]x[0,1].
!
    call random_number ( harvest = point2(1:dim_num,1:batch) )
!
!  Remap the points to the box [X1,X2] x [Y1,Y2].
!
    point2(1,1:batch) = ( x1 - h ) &
      + point2(1,1:batch) * ( x2 - x1 + 2.0D+00 * h )

    point2(2,1:batch) = ( y1 - h ) &
      + point2(2,1:batch) * ( y2 - y1 + 2.0D+00 * h )

    call random_number ( harvest = test(1:batch) )

    call p99_weight ( dim_num, batch, point2, weight )

    do j = 1, batch

      if ( test(j) <= weight(j) ) then
        found = found + 1
        if ( found <= n ) then
          point(1:2,found) = point2(1:2,j)
        end if
      end if

    end do

  end do

  return
end
subroutine p99_weight ( dim_num, n, p, w )

!*****************************************************************************80
!
!! P99_WEIGHT returns the weight for problem 99.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer P(DIM_NUM,N), the point coordinates.
!
!    Output, real ( kind = 8 ) W(N), the weight associated with each point.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) w(n)

  do j = 1, n

    w(j) = ( 0.5D+00 &
      + 15.5D+00 * ( 1.0D+00 + p(1,j) ) * ( 1.0D+00 - p(1,j) ) &
                 * ( 1.0D+00 + p(2,j) ) * ( 1.0D+00 - p(2,j) ) ) / 16.0D+00

    w(j) = max ( w(j), 0.03125D+00 )
    w(j) = min ( w(j), 1.0D+00 )

  end do

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer SAMPLE_NUM, the number of samples.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,N), the point coordinates.
!
!    Input, real ( kind = 8 ) SAMPLE(DIM_NUM,SAMPLE_NUM), the sample
!    coordinates.
!
!    Output, integer NEAREST(SAMPLE_NUM), the index of the nearest
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
!    Input, integer DIM_NUM, the spatial dimension.
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
      '#  created by routine CVT_WRITE in CVT_TRIANGULATION.F90'
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
