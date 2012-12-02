program main

!*****************************************************************************80
!
!! MAIN is the main program for QUALITY_PRB.
!
!  Discussion:
!
!    QUALITY_PRB calls the QUALITY routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUALITY_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the QUALITY library.'

  call test_cvt ( )

  call test_halton ( )

  call test_sphere ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUALITY_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test_cvt ( )

!*****************************************************************************80
!
!! TEST_CVT tests a dataset in the unit hypercube.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  character ( len = 80 ) input_filename
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns
  external sample_hypercube_uniform
  integer ( kind = 4 ) seed_init
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node
  integer ( kind = 4 ) triangle_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: z

  ns = 100000
  seed_init = 123456789
  input_filename = 'cvt_02_00100.txt'

  call r8mat_header_read ( input_filename, dim_num, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_CVT:'
  write ( *, '(a)' ) '  Measures of uniform point dispersion.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The pointset was read from "' &
    // trim ( input_filename ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sampling routine is "SAMPLE_HYPERCUBE_UNIFORM".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Spatial dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)' ) '  The number of points N =         ', n
  write ( *, '(a,i12)' ) '  The number of sample points NS = ', ns
  write ( *, '(a,i12)' ) '  The random number SEED_INIT =    ', seed_init
  write ( *, '(a)' ) ' '

  allocate ( z(1:dim_num,1:n) )

  call r8mat_data_read ( input_filename, dim_num, n, z )

  call r8mat_transpose_print_some ( dim_num, n, z, 1, 1, 5, 5, &
    '  5x5 portion of data read from file:' )
!
!  For 2D datasets, determine the Delaunay triangulation.
!
  if ( dim_num == 2 )then

    allocate ( triangle_node(3,3*n) )
    allocate ( triangle_neighbor(3,3*n) )

    call dtris2 ( n, z, triangle_num, triangle_node, triangle_neighbor )

  else

    triangle_num = 0
    allocate ( triangle_node(0,0) )
    allocate ( triangle_neighbor(0,0) )

  end if

  if ( dim_num == 2 ) then
    call test005 ( n, z, triangle_num, triangle_node )
    call test006 ( n, z, triangle_num, triangle_node )
  end if

  call test007 ( dim_num, n, z )
  call test01 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test02 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test03 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test04 ( dim_num, n, z )
  call test05 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test06 ( dim_num, n, z )
  call test07 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test08 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )

  if ( dim_num == 2 ) then
    call test083 ( n, z, triangle_num, triangle_node )
  end if

  call test085 ( dim_num, n, z )
  call test09 ( dim_num, n, z )
  call test10 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test11 ( dim_num, n, z )

  deallocate ( triangle_node )
  deallocate ( triangle_neighbor )
  deallocate ( z )

  return
end
subroutine test_halton ( )

!*****************************************************************************80
!
!! TEST_HALTON tests a dataset in the unit hypercube.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  character ( len = 80 ) input_filename
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns
  external sample_hypercube_uniform
  integer ( kind = 4 ) seed_init
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node
  integer ( kind = 4 ) triangle_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: z

  ns = 100000
  seed_init = 123456789
  input_filename = 'halton_02_00100.txt'

  call r8mat_header_read ( input_filename, dim_num, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_HALTON:'
  write ( *, '(a)' ) '  Measures of uniform point dispersion.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The pointset was read from "' &
    // trim ( input_filename ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sampling routine is "SAMPLE_HYPERCUBE_UNIFORM".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Spatial dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)' ) '  The number of points N =         ', n
  write ( *, '(a,i12)' ) '  The number of sample points NS = ', ns
  write ( *, '(a,i12)' ) '  The random number SEED_INIT =    ', seed_init
  write ( *, '(a)' ) ' '

  allocate ( z(1:dim_num,1:n) )

  call r8mat_data_read ( input_filename, dim_num, n, z )

  call r8mat_transpose_print_some ( dim_num, n, z, 1, 1, 5, 5, &
    '  5x5 portion of data read from file:' )
!
!  For 2D datasets, determine the Delaunay triangulation.
!
  if ( dim_num == 2 )then

    allocate ( triangle_node(3,3*n) )
    allocate ( triangle_neighbor(3,3*n) )

    call dtris2 ( n, z, triangle_num, triangle_node, triangle_neighbor )

  else

    triangle_num = 0
    allocate ( triangle_node(0,0) )
    allocate ( triangle_neighbor(0,0) )

  end if

  if ( dim_num == 2 ) then
    call test005 ( n, z, triangle_num, triangle_node )
    call test006 ( n, z, triangle_num, triangle_node )
  end if

  call test007 ( dim_num, n, z )
  call test01 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test02 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test03 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test04 ( dim_num, n, z )
  call test05 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test06 ( dim_num, n, z )
  call test07 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test08 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )

  if ( dim_num == 2 ) then
    call test083 ( n, z, triangle_num, triangle_node )
  end if

  call test085 ( dim_num, n, z )
  call test09 ( dim_num, n, z )
  call test10 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init )
  call test11 ( dim_num, n, z )

  deallocate ( triangle_node )
  deallocate ( triangle_neighbor )
  deallocate ( z )

  return
end
subroutine test_sphere ( )

!*****************************************************************************80
!
!! TEST_SPHERE tests a dataset that is not in the unit hypercube.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) input_filename
  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) ns
  external sample_sphere_uniform
  integer ( kind = 4 ) seed_init
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node
  integer ( kind = 4 ) triangle_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: z

  ns = 100000
  seed_init = 123456789
  input_filename = 'sphere_02_00100.txt'

  call r8mat_header_read ( input_filename, dim_num, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SPHERE:'
  write ( *, '(a)' ) '  Measures of uniform point dispersion.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The pointset was read from "' &
    // trim ( input_filename ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sampling routine is "SAMPLE_SPHERE_UNIFORM".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Spatial dimension DIM_NUM =      ', dim_num
  write ( *, '(a,i12)' ) '  The number of points N =         ', n
  write ( *, '(a,i12)' ) '  The number of sample points NS = ', ns
  write ( *, '(a,i12)' ) '  The random number SEED_INIT =    ', seed_init
  write ( *, '(a)' ) ' '

  allocate ( z(1:dim_num,1:n) )

  call r8mat_data_read ( input_filename, dim_num, n, z )

  call r8mat_transpose_print_some ( dim_num, n, z, 1, 1, 5, 5, &
    '  5x5 portion of data read from file:' )
!
!  For 2D datasets, determine the Delaunay triangulation.
!
  if ( dim_num == 2 )then

    allocate ( triangle_node(3,3*n) )
    allocate ( triangle_neighbor(3,3*n) )

    call dtris2 ( n, z, triangle_num, triangle_node, triangle_neighbor )

  else

    triangle_num = 0
    allocate ( triangle_node(0,0) )
    allocate ( triangle_neighbor(0,0) )

  end if

  if ( dim_num == 2 ) then
    call test005 ( n, z, triangle_num, triangle_node )
    call test006 ( n, z, triangle_num, triangle_node )
  end if

  call test007 ( dim_num, n, z )
  call test01 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init )
  call test02 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init )
  call test03 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init )
  call test04 ( dim_num, n, z )
  call test05 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init )
  call test06 ( dim_num, n, z )
  call test07 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init )
  call test08 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init )

  if ( dim_num == 2 ) then
    call test083 ( n, z, triangle_num, triangle_node )
  end if

  call test085 ( dim_num, n, z )
!
!  We don't call TEST09, because the test only works in the unit hypercube.
!
! call test09 ( dim_num, n, z )
  call test10 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init )
  call test11 ( dim_num, n, z )

  deallocate ( triangle_node )
  deallocate ( triangle_neighbor )
  deallocate ( z )

  return
end
subroutine test005 ( n, z, triangle_num, triangle_node )

!*****************************************************************************80
!
!! TEST005 tests ALPHA_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), parameter :: triangle_order = 3

  real ( kind = 8 ) alpha_measure
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  real ( kind = 8 ) z(2,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  ALPHA_MEASURE computes the ALPHA measure of quality.'

  write ( *, '(a,f14.6)' ) '  The minimal angle measure    ALPHA = ', &
    alpha_measure ( n, z, triangle_order, triangle_num, triangle_node )

  return
end
subroutine test006 ( n, z, triangle_num, triangle_node )

!*****************************************************************************80
!
!! TEST006 tests AREA_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), parameter :: triangle_order = 3

  real ( kind = 8 ) area_measure
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  real ( kind = 8 ) z(2,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  AREA_MEASURE computes the AREA measure of quality.'
  write ( *, '(a,f14.6)' ) '  The area ratio measure        AREA = ', &
    area_measure ( n, z, triangle_order, triangle_num, triangle_node )

  return
end
subroutine test007 ( dim_num, n, z )

!*****************************************************************************80
!
!! TEST007 tests BETA_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) beta_measure
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  BETA_MEASURE computes the BETA measure of quality.'

  write ( *, '(a,f14.6)' ) '  Relative spacing deviation BETA =    ', &
    beta_measure ( dim_num, n, z )

  return
end
subroutine test01 ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! TEST01 tests CHI_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) chi_measure
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CHI_MEASURE computes the CHI measure of quality.'

  write ( *, '(a,f14.6)' ) '  The regularity measure         Chi = ', &
    chi_measure ( dim_num, n, z, ns, sample_routine, seed_init )

  return
end
subroutine test02 ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! TEST02 tests D_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) d_measure
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  D_MEASURE computes the D measure of quality.'

  write ( *, '(a,g14.6)' ) '  2nd moment determinant measure   D = ', &
    d_measure ( dim_num, n, z, ns, sample_routine, seed_init )

  return
end
subroutine test03 ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! TEST03 tests E_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) e_measure
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  E_MEASURE computes the E measure of quality.'

  write ( *, '(a,g14.6)' ) '  Voronoi energy measure           E = ', &
    e_measure ( dim_num, n, z, ns, sample_routine, seed_init )

  return
end
subroutine test04 ( dim_num, n, z )

!*****************************************************************************80
!
!! TEST04 tests GAMMA_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) gamma_measure
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  GAMMA_MEASURE computes the Gamma measure of quality.'

  write ( *, '(a,f14.6)' ) '  The mesh ratio               Gamma = ', &
    gamma_measure ( dim_num, n, z )

  return
end
subroutine test05 ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! TEST05 tests H_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) h_measure
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  H_MEASURE computes the H measure of quality.'

  write ( *, '(a,f14.6)' ) '  The point distribution norm      H = ', &
    h_measure ( dim_num, n, z, ns, sample_routine, seed_init )

  return
end
subroutine test06 ( dim_num, n, z )

!*****************************************************************************80
!
!! TEST06 tests LAMBDA_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) lambda_measure
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  LAMBDA_MEASURE computes the Lambda measure of quality.'

  write ( *, '(a,f14.6)' ) '  The COV measure             Lambda = ', &
    lambda_measure ( dim_num, n, z )

  return
end
subroutine test07 ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! TEST07 tests MU_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) mu_measure
  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  MU_MEASURE computes the Mu measure of quality.'

  write ( *, '(a,f14.6)' ) '  The point distribution ratio    Mu = ', &
    mu_measure ( dim_num, n, z, ns, sample_routine, seed_init )

  return
end
subroutine test08 ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! TEST08 tests NU_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) ns
  real ( kind = 8 ) nu_measure
  external sample_routine
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  NU_MEASURE computes the Nu measure of quality.'

  write ( *, '(a,f14.6)' ) '  The cell volume deviation       Nu = ', &
    nu_measure ( dim_num, n, z, ns, sample_routine, seed_init )

  return
end
subroutine test083 ( n, z, triangle_num, triangle_node )

!*****************************************************************************80
!
!! TEST083 tests Q_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), parameter :: triangle_order = 3

  real ( kind = 8 ) q_measure
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  real ( kind = 8 ) z(2,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST083'
  write ( *, '(a)' ) '  Q_MEASURE computes the Q measure of quality.'

  write ( *, '(a,f14.6)' ) '  The triangle shape measure       Q = ', &
    q_measure ( n, z, triangle_order, triangle_num, triangle_node )

  return
end
subroutine test085 ( dim_num, n, z )

!*****************************************************************************80
!
!! TEST085 tests R0_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) r0_measure
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  R0_MEASURE computes the R0 measure of quality.'

  write ( *, '(a,f14.6)' ) '  The Riesz energy with S = 0,    R0 = ', &
    r0_measure ( dim_num, n, z )

  return
end
subroutine test09 ( dim_num, n, z )

!*****************************************************************************80
!
!! TEST09 tests SPHERE_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) sphere_measure
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  SPHERE_MEASURE computes the sphere measure of quality.'

  write ( *, '(a,f14.6)' ) '  Nonintersecting sphere volume    S = ', &
    sphere_measure ( dim_num, n, z )

  return
end
subroutine test10 ( dim_num, n, z, ns, sample_routine, seed_init )

!*****************************************************************************80
!
!! TEST10 tests TAU_MEASURE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) ns
  external sample_routine
  integer ( kind = 4 ) seed_init
  real ( kind = 8 ) tau_measure
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  TAU_MEASURE computes the Tau measure of quality.'

  write ( *, '(a,f14.6)' ) '  2nd moment trace measure       Tau = ', &
    tau_measure ( dim_num, n, z, ns, sample_routine, seed_init )

  return
end
subroutine test11 ( dim_num, n, z )

!*****************************************************************************80
!
!! TEST11 tests POINTSET_SPACING.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) gamma(n)
  real ( kind = 8 ) gamma_ave
  real ( kind = 8 ) gamma_max
  real ( kind = 8 ) gamma_min
  real ( kind = 8 ) gamma_std
  real ( kind = 8 ) z(dim_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  POINTSET_SPACING computes pointset spacing parameters.'

  call pointset_spacing ( dim_num, n, z, gamma )

  gamma_min = minval ( gamma(1:n) )
  gamma_max = maxval ( gamma(1:n) )

  gamma_ave = sum ( gamma(1:n) ) / real ( n, kind = 8 )

  if ( 1 < n ) then
    gamma_std = sqrt ( sum ( ( gamma(1:n) - gamma_ave )**2 ) &
      / real ( n - 1, kind = 8 ) )
  else
    gamma_std = 0.0D+00
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  Minimum spacing          GAMMA_MIN = ', &
    gamma_min
  write ( *, '(a,f14.6)' ) '  Average spacing          GAMMA_AVE = ', &
    gamma_ave
  write ( *, '(a,f14.6)' ) '  Maximum spacing          GAMMA_MAX = ', &
    gamma_max
  write ( *, '(a,f14.6)' ) '  Spacing standard dev     GAMMA_STD = ', &
    gamma_std

  return
end
