program main

!*****************************************************************************80
!
!! MAIN is the main program for CVT_MOD_PRB.
!
!  Discussion:
!
!    CVT_MOD_PRB calls a set of problems for CVT_MOD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_MOD_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CVT_MOD library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_MOD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CVT_MOD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 400

  real ( kind = 8 ) change_l2
  integer ( kind = 4 ), parameter :: cvt_steps = 50
  character ( len = 80 ) :: file_in_name = 'none'
  character ( len = 80 ) :: file_out_name = 'cvt_mod_02_00400.txt'
  real ( kind = 8 ) generator(m,n)
  integer ( kind = 4 ) i
  logical reset
  integer ( kind = 4 ) sample_function_cvt
  integer ( kind = 4 ) sample_function_init
  integer ( kind = 4 ), parameter :: sample_num_cvt = 100000
  integer ( kind = 4 ), parameter :: sample_num_steps = 50
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789
  real ( kind = 8 ), dimension ( m ) :: width = (/ 1.0D+00, 2.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CVT_MOD computes a Centroidal Voronoi Tessellation'
  write ( *, '(a)' ) '    modulo 1.'
  write ( *, '(a)' ) ' '

  sample_function_init = -1
  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Spatial dimension M =        ', m
  write ( *, '(a,i12)' ) '  Number of generators =       ', n
  write ( *, '(a,i12)' ) '  Initial random number seed = ', seed
  write ( *, '(a)' ) ' '

  if ( sample_function_init == -1 ) then
    write ( *, '(a)' ) '  Initialize using RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_init == 0 ) then
    write ( *, '(a)' ) '  Initialize using UNIFORM.'
  else if ( sample_function_init == 1 ) then
    write ( *, '(a)' ) '  Initialize using HALTON.'
  else if ( sample_function_init == 2 ) then
    write ( *, '(a)' ) '  Initialize using GRID.'
  else if ( sample_function_init == 3 ) then
    write ( *, '(a)' ) '  Initialize from file "' &
      // trim ( file_in_name ) // '".'
  end if

  if ( sample_function_cvt == -1 ) then
    write ( *, '(a)' ) '  Sample using RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_cvt == 0 ) then
    write ( *, '(a)' ) '  Sample using UNIFORM.'
  else if ( sample_function_cvt == 1 ) then
    write ( *, '(a)' ) '  Sample using HALTON.'
  else if ( sample_function_cvt == 2 ) then
    write ( *, '(a)' ) '  Sample using GRID.'
  end if

  write ( *, '(a,i6)' ) '  Number of sample points = ', sample_num_cvt
  write ( *, '(a,i6)' ) '  Number of sample steps =  ', sample_num_steps
!
!  Initialize the generators.
!
  reset = .true.

  call region_sampler_mod ( m, n, sample_num_cvt, generator, &
    sample_function_init, reset, seed, width )

  do i = 1, cvt_steps

    call cvt_iteration_mod ( m, n, generator, width, sample_num_cvt, &
      sample_function_cvt, seed, change_l2 )

  end do

  call r8mat_write ( file_out_name, m, n, generator )

  return
end
