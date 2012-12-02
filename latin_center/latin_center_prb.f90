program main

!*****************************************************************************80
!
!! MAIN is the main program for LATIN_CENTER_PRB.
!
!  Discussion:
!
!    LATIN_CENTER_PRB tests the Latin Center Square routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_save

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_CENTER_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LATIN_CENTER library.'

  call test00 ( seed )

  seed_save = seed
  call test01 ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_CENTER_PRB:'
  write ( *, '(a)' ) '  Repeat TEST01 with a different seed from the first run.'

  call test01 ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_CENTER_PRB:'
  write ( *, '(a)' ) '  Repeat TEST01 with the same seed as the first run.'

  seed = seed_save
  call test01 ( seed )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_CENTER_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test00 ( seed )

!*****************************************************************************80
!
!! TEST00 tests GET_SEED and RANDOM_INITIALIZE.
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
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) seed

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
!! TEST01 tests LATIN_CENTER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: point_num = 10

  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(dim_num,point_num)

  write ( *, '(a)'     ) ' '
  write ( *, '(a)'     ) 'TEST01'
  write ( *, '(a)'     ) '  LATIN_CENTER chooses a Latin cell arrangement,'
  write ( *, '(a)'     ) '  and returns the centers of those cells.'
  write ( *, '(a)'     ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension =  ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points =   ', point_num
  write ( *, '(a,i12)' ) '  Random number SEED = ', seed

  call latin_center ( dim_num, point_num, seed, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Latin Center Square points:'
  write ( *, '(a)' ) ' '

  do j = 1, point_num
    write ( *, '(2x,f10.4,2x,f10.4)' ) x(1:dim_num,j)
  end do

  return
end
