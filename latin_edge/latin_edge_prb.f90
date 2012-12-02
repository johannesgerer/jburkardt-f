program main

!*****************************************************************************80
!
!! MAIN is the main program for LATIN_EDGE_PRB.
!
!  Discussion:
!
!    LATIN_EDGE_PRB tests the Latin hypercube routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer seed
  integer seed_save

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_EDGE_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LATIN_EDGE library.'

  call test00 ( seed )

  seed_save = seed
  call test01 ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_EDGE_PRB:'
  write ( *, '(a)' ) '  Repeat TEST01 with a different seed from the first run.'

  call test01 ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_EDGE_PRB:'
  write ( *, '(a)' ) '  Repeat TEST01 with the same seed as the first run.'

  seed = seed_save
  call test01 ( seed )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_EDGE_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test00 ( seed )

!*****************************************************************************80
!
!! TEST00 tests GET_SEED, RANDOM_INITIALIZE.
!
!  Discussion:
!
!    If the LATIN routines are set to use the system random number
!    generator, then setting SEED here is enough.
!
!    If they are set to use the portable UNIFORM routine, then
!    SEED must be passed to those routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST00'
  write ( *, '(a)' ) '  GET_SEED returns a seed for the random number'
  write ( *, '(a)' ) '  generator, based on the current time.'
  write ( *, '(a)' ) '  RANDOM_INITIALIZE uses that seed to initialize'
  write ( *, '(a)' ) '  the FORTRAN90 random number generator.'

  call get_seed ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  GET_SEED returns SEED = ', seed

  call random_initialize ( seed )

  return
end
subroutine test01 ( seed )

!*****************************************************************************80
!
!! TEST01 tests LATIN_EDGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: dim_num = 2
  integer, parameter :: point_num = 11

  integer j
  integer seed
  real ( kind = 8 ) x(dim_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  LATIN_EDGE chooses a Latin cell arrangement,'
  write ( *, '(a)' ) '  which includes the edge points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)'  ) '  Spatial dimension =  ', dim_num
  write ( *, '(a,i6)'  ) '  Number of points =   ', point_num
  write ( *, '(a,i12)' ) '  Random number SEED = ', seed

  call latin_edge ( dim_num, point_num, seed, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Latin Edge Square points:'
  write ( *, '(a)' ) ' '

  do j = 1, point_num
    write ( *, '(2x,2f10.4)' ) x(1:dim_num,j)
  end do

  return
end
