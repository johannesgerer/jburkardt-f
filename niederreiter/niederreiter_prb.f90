program main

!*****************************************************************************80
!
!! MAIN is the main program for NIEDERREITER_PRB.
!
!  Discussion:
!
!    NIEDERREITER_PRB calls a set of problems for NIEDERREITER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) dim_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NIEDERREITER_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NIEDERREITER library.'

  base = 2
  call test01 ( base )

  base = 3
  call test01 ( base )

  base = 13
  call test01 ( base )

  base = 2
  call test02 ( base )

  base = 3
  call test02 ( base )

  base = 2
  dim_num = 20
  call test03 ( base, dim_num )

  base = 2
  dim_num = 29
  call test03 ( base, dim_num )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NIEDERREITER_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( base )

!*****************************************************************************80
!
!! TEST01 tests NIEDERREITER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE, the base to use in the computation.
!    BASE should be a prime, or a power of a prime.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_max = 4

  integer ( kind = 4 ) base
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(dim_max)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_out

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  NIEDERREITER computes the next element of '
  write ( *, '(a)' ) '  a Niederreiter quasirandom sequence using base BASE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we call NIEDERREITER repeatedly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Using base BASE =      ', base

  do dim_num = 2, dim_max

    seed = 0

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Seed    Seed     Niederreiter'
    write ( *, '(a)' ) '      In     Out'
    write ( *, '(a)' ) ' '
    do i = 0, 110
      seed_in = seed
      call niederreiter ( dim_num, base, seed, r )
      seed_out = seed
      if ( i <= 11 .or. 95 <= i ) then
        write ( *, '(i8,i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
      else if ( i == 12 ) then
        write ( *, '(a)' ) '......................'
      end if
    end do

  end do

  return
end
subroutine test02 ( base )

!*****************************************************************************80
!
!! TEST02 tests NIEDERREITER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE, the base to use in the computation.
!    BASE should be a prime, or a power of a prime.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_out

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  NIEDERREITER computes the next element of '
  write ( *, '(a)' ) '  a Niederreiter quasirandom sequence using base BASE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we demonstrate how the SEED can be'
  write ( *, '(a)' ) '  manipulated to skip ahead in the sequence, or'
  write ( *, '(a)' ) '  to come back to any part of the sequence.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Using base BASE =           ', base
  write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num

  seed = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed     Niederreiter'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 0, 10
    seed_in = seed
    call niederreiter ( dim_num, base, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Jump ahead by increasing SEED:'
  write ( *, '(a)' ) ' '

  seed = 100

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed     Niederreiter'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    seed_in = seed
    call niederreiter ( dim_num, base, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Jump back by decreasing SEED:'
  write ( *, '(a)' ) ' '

  seed = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed     Niederreiter'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 0, 10
    seed_in = seed
    call niederreiter ( dim_num, base, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Jump ahead by increasing SEED:'
  write ( *, '(a)' ) ' '

  seed = 98

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed     Niederreiter'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    seed_in = seed
    call niederreiter ( dim_num, base, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
  end do

  return
end
subroutine test03 ( base, dim_num )

!*****************************************************************************80
!
!! TEST03 tests NIEDERREITER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE, the base to use in the computation.
!    BASE should be a prime, or a power of a prime.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  real ( kind = 8 ), allocatable :: r(:)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_out

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  NIEDERREITER computes the next element of '
  write ( *, '(a)' ) '  a Niederreiter quasirandom sequence using base BASE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Using base BASE =           ', base
  write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num

  seed = 0
  allocate ( r(1:dim_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed     Niederreiter'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    seed_in = seed
    call niederreiter ( dim_num, base, seed, r )
    seed_out = seed
    write ( *, '(2i8)', ADVANCE = 'NO' ) seed_in, seed_out
    do dim = 1, dim_num
      write ( *, '(f10.4)', ADVANCE = 'NO' ) r(dim)
      if ( mod ( dim, 5 ) == 0 .and. dim /= dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)', ADVANCE = 'NO' ) '                '
      end if
    end do
    write ( *, * ) ' '
  end do

  deallocate ( r )

  return
end
