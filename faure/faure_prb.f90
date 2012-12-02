program main

!*****************************************************************************80
!
!! MAIN is the main program for FAURE_PRB.
!
!  Discussion:
!
!    FAURE_PRB calls a set of problems for FAURE.
!
!  Modified:
!
!    04 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FAURE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FAURE library.'

  call test005 ( )
  call test006 ( )
  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FAURE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests BINOMIAL_TABLE.
!
!  Modified:
!
!    08 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ), dimension ( 0:m, 0:n ) :: coef
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: qs = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  BINOMIAL_TABLE computes a table of binomial.'
  write ( *, '(a)' ) '  coefficients mod QS.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here, QS = ', qs

  call binomial_table ( qs, m, n, coef )

  write ( *, '(a)' ) ' '
  write ( *, '(a,8i8)' ) '   I/J', ( j, j = 0, n )
  write ( *, '(a)' ) ' '

  do i = 0, m
    write ( *, '(2x,i2,2x,8i8)' ) i, coef(i,0:n)
  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests I4_LOG_I4.
!
!  Modified:
!
!    09 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_log_i4
  integer ( kind = 4 ) j4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  I4_LOG_I4: logarith of I4 base J4,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        I4        J4 I4_LOG_I4'
  write ( *, '(a)' ) ' '

  do j4 = 2, 5
    do i4 = 0, 10
      write ( *, '(2x, i8, 2x, i8, 2x, i8 )' ) i4, j4, i4_log_i4 ( i4, j4 )
    end do
    write ( *, '(a)' ) ' '
  end do

  return
end

subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests FAURE.
!
!  Modified:
!
!    31 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_max = 4

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) qs
  integer ( kind = 4 ) prime_ge
  real ( kind = 8 ) r(dim_max)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_out

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  FAURE computes the next element of a Faure sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we call FAURE repeatedly.'

  do dim_num = 2, dim_max

    seed = -1
    qs = prime_ge ( dim_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num
    write ( *, '(a,i8)' ) '  The prime base QS       =   ', qs
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Seed    Seed    Faure'
    write ( *, '(a)' ) '      In     Out'
    write ( *, '(a)' ) ' '
    do i = 1, 10
      seed_in = seed
      call faure ( dim_num, seed, r )
      seed_out = seed
      write ( *, '(2i8,4f10.6)' ) seed_in, seed_out, r(1:dim_num)
    end do

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests FAURE.
!
!  Modified:
!
!    31 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) qs
  integer ( kind = 4 ) prime_ge
  real ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_out

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  FAURE computes the next element of a Faure sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we demonstrate how the SEED can be'
  write ( *, '(a)' ) '  manipulated to skip ahead in the sequence, or'
  write ( *, '(a)' ) '  to come back to any part of the sequence.'

  qs = prime_ge ( dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Using dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  The prime base QS       = ', qs

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Note that on the first call to FAURE, if'
  write ( *, '(a)' ) '  SEED is negative, it is reset to a value that'
  write ( *, '(a)' ) '  is the recommended starting point:'

  seed = -1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed    Faure'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    seed_in = seed
    call faure ( dim_num, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.6)' ) seed_in, seed_out, r(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  However, if the input value of SEED is 0,'
  write ( *, '(a)' ) '  then no initial skipping is done.'

  seed = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed    Faure'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    seed_in = seed
    call faure ( dim_num, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.6)' ) seed_in, seed_out, r(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Jump ahead by increasing SEED:'
  write ( *, '(a)' ) ' '

  seed = 100

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed    Faure'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    seed_in = seed
    call faure ( dim_num, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.6)' ) seed_in, seed_out, r(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Jump back by decreasing SEED:'
  write ( *, '(a)' ) ' '

  seed = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed    Faure'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    seed_in = seed
    call faure ( dim_num, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.6)' ) seed_in, seed_out, r(1:dim_num)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests FAURE.
!
!  Modified:
!
!    04 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_base = 10
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) qs
  integer ( kind = 4 ) prime_ge
  real ( kind = 8 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_out

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  FAURE computes the next element of a Faure sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we try some large dimensions.'

  do dim_num = dim_base, 6 * dim_base, dim_base

    allocate ( r(1:dim_num) )

    seed = -1
    qs = prime_ge ( dim_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Using dimension DIM_NUM = ', dim_num
    write ( *, '(a,i8)' ) '  The prime base QS       = ', qs
    do i = 1, 2
      seed_in = seed
      call faure ( dim_num, seed, r )
      seed_out = seed
      write ( *, '(a)' ) ' '
      write ( *, '(a,i10)' ) '  Seed in =  ', seed_in
      write ( *, '(a,i10)' ) '  Seed out = ', seed_out
      write ( *, '(a,i2,a)' ) '  R(1:', dim_num,') = '
      write ( *, '(5(2x,f10.6))' ) r(1:dim_num)
    end do

    deallocate ( r )

  end do

  return
end
