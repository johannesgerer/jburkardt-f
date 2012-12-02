program main

!*****************************************************************************80
!
!! MAIN is the main program for UNIFORM_PRB.
!
!  Discussion:
!
!    UNIFORM_PRB calls sample problems for the UNIFORM library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNIFORM_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the UNIFORM library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test065 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test111 ( )
  call test112 ( )
  call test118 ( )
  call test119 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test21 ( )
  call test22 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNIFORM_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests C4_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ) c4_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  C4_UNIFORM_01 computes pseudorandom complex values '
  write ( *, '(a)' ) '  uniformly distributed in the unit circle.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8,2x,f14.8)' ) i, c4_uniform_01 ( seed )
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests C4VEC_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  C4VEC_UNIFORM_01 computes pseudorandom complex values '
  write ( *, '(a)' ) '  uniformly distributed in the unit circle.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

  call c4vec_uniform_01 ( n, seed, c )

  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8,2x,f14.8)' ) i, c(i)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests C8_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  C8_UNIFORM_01 computes pseudorandom double precision'
  write ( *, '(a)' ) '  complex values uniformly distributed in the unit'
  write ( *, '(a)' ) '  circle.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8,2x,f14.8)' ) i, c8_uniform_01 ( seed )
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests C8VEC_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  C8VEC_UNIFORM_01 computes pseudorandom '
  write ( *, '(a)' ) '  double precision complex values uniformly distributed '
  write ( *, '(a)' ) '  in the unit circle.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

  call c8vec_uniform_01 ( n, seed, c )

  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8,2x,f14.8)' ) i, c(i)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests CH_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ch_uniform_ab
  character chi
  character clo
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  CH_UNIFORM_AB computes pseudorandom characters '
  write ( *, '(a)' ) '  in an interval [CLO,CHI].'

  clo = 'A'
  chi = 'J'
  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The lower endpoint CLO = "' // clo // '".'
  write ( *, '(a)' ) '  The upper endpoint CHI = "' // chi // '".'
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,a1)' ) i, ch_uniform_ab ( clo, chi, seed )
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests GET_SEED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_old

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  GET_SEED picks an initial seed value for UNIFORM.'
  write ( *, '(a)' ) '  The value chosen should vary over time, because'
  write ( *, '(a)' ) '  the seed is based on reading the clock.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This is just the "calendar" clock, which does'
  write ( *, '(a)' ) '  not change very fast, so calling GET_SEED several'
  write ( *, '(a)' ) '  times in a row may result in the same value.'

  seed = 12345678
  seed_old = seed

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Initial SEED is ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Next 3 values of R8_UNIFORM:'
  write ( *, '(a)' ) ' '

  do j = 1, 3
    write ( *, '(2x,f14.8)' ) r8_uniform_01 ( seed )
  end do

  do i = 1, 4

    do 

      call get_seed ( seed )

      if ( seed /= seed_old ) then
        seed_old = seed
        exit
      end if

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  New seed from GET_SEED is ', seed
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Next 3 values of R8_UNIFORM_01:'
    write ( *, '(a)' ) ' '

    do j = 1, 3
      write ( *, '(2x,f14.8)' ) r8_uniform_01 ( seed )
    end do

  end do

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests I4_SEED_ADVANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i4_seed_advance
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_new
  integer ( kind = 4 ) step

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  I4_SEED_ADVANCE advances the seed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step        SEED input       SEED output'
  write ( *, '(a)' ) ' '

  seed_new = 12345

  do step = 1, 10

    seed = seed_new
    seed_new = i4_seed_advance ( seed )

    write ( *, '(2x,i4,2x,i16,2x,i16)' ) step, seed, seed_new

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests I4_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: a = 6
  integer ( kind = 4 ), parameter :: b = 10

  integer ( kind = 4 ) freq(a:b)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10000

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  I4_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  freq(a:b) = 0

  do test = 1, test_num

    j = i4_uniform_ab ( a, b, seed )

    if ( j < a ) then
      write ( *, '(a,i8)' ) '  Illegal value J = ', j
    else if ( j <= b ) then
      freq(j) = freq(j) + 1
    else
      write ( *, '(a,i8)' ) '  Illegal value J = ', j
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I    Frequency'
  write ( *, '(a)' ) ' '
  do i = a, b
    write ( *, '(2x,i8,2x,i8)' ) i, freq(i)
  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests I4_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: a = -100
  integer ( kind = 4 ), parameter :: b = 200
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  I4_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 20

    j = i4_uniform_ab ( a, b, seed )

    write ( *, '(2x,i8,2x,i8)' ) i, j

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests I4_UNIFORM_0I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_0i
  real ( kind = 4 ) mean
  integer ( kind = 4 ), parameter :: n = 1000
  integer ( kind = 4 ) seed
  real ( kind = 4 ) variance
  integer ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  I4_UNIFORM_0I samples a uniform random'
  write ( *, '(a)' ) '  integer distribution in [0,2**31-1].'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Starting with seed = ', seed

  do i = 1, n
    x(i) = i4_uniform_0i ( seed )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First few values:'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,i12)' ) i, x(i)
  end do

  mean = sum ( real ( x(1:n), kind = 4 ) / real ( n, kind = 4 ) )

  variance = sum ( ( real ( x(1:n), kind = 4 ) - mean )**2 ) &
                / real ( n - 1, kind = 4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of values computed was N = ', n
  write ( *, '(a,g14.6)' ) '  Average value was ', mean
  write ( *, '(a,i12)' ) '  Minimum value was ', minval ( x(1:n) )
  write ( *, '(a,i12)' ) '  Maximum value was ', maxval ( x(1:n) )
  write ( *, '(a,g14.6)' ) '  Variance was ', variance

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests I4VEC_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10000
  integer ( kind = 4 ), parameter :: a = 6
  integer ( kind = 4 ), parameter :: b = 10

  integer ( kind = 4 ) freq(a:b)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  I4VEC_UNIFORM_AB computes a vector of pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  freq(a:b) = 0

  call i4vec_uniform_ab ( n, a, b, seed, i4vec )

  do i = 1, n

    if ( i4vec(i) < a ) then
      write ( *, '(a,i8)' ) '  Illegal value J = ', i4vec(i)
    else if ( i4vec(i) <= b ) then
      freq(i4vec(i)) = freq(i4vec(i)) + 1
    else
      write ( *, '(a,i8)' ) '  Illegal value J = ', i4vec(i)
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I    Frequency'
  write ( *, '(a)' ) ' '
  do i = a, b
    write ( *, '(2x,i8,2x,i8)' ) i, freq(i)
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests I8_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 8 ) a
  integer ( kind = 8 ) b
  integer ( kind = 8 ) i
  integer ( kind = 8 ) i8_uniform_ab
  integer ( kind = 8 ) seed
! integer ( kind = 8 ), parameter :: seed_init = 123456789
  integer ( kind = 8 ), parameter :: seed_init = 123456789987654321_8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  I8_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'

! a = 100000
  a = 1000000000_8
! b = 800000
  b = 8000000000_8

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i24)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i24)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i24)' ) '  The initial seed is    ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,i24)' ) i, i8_uniform_ab ( a, b, seed )
  end do

  return
end
subroutine test111 ( )

!*****************************************************************************80
!
!! TEST111 tests L_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  logical l_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST111'
  write ( *, '(a)' ) '  L_UNIFORM computes pseudorandom logical values.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,l1)' ) i, l_uniform ( seed )
  end do

  return
end
subroutine test112 ( )

!*****************************************************************************80
!
!! TEST112 tests LVEC_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  logical lvec(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST112'
  write ( *, '(a)' ) '  LVEC_UNIFORM_01 computes a vector of'
  write ( *, '(a)' ) '  pseudorandom logical values.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

  call lvec_uniform ( n, seed, lvec )

  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,l1)' ) i, lvec(i)
  end do

  return
end
subroutine test118 ( )

!*****************************************************************************80
!
!! TEST118 tests LCRG_ANBN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) an
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bn
  integer ( kind = 4 ) c
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST118'
  write ( *, '(a)' ) '  LCRG_ANBN determines a linear congruential random'
  write ( *, '(a)' ) '  number generator equivalent to N steps of a given one.'
!
!  These parameters define the old (1969) IBM 360 random number generator:
!
  a = 16807
  b = 0
  c = 2147483647

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LCRG parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  A = ', a
  write ( *, '(a,i12)' ) '  B = ', b
  write ( *, '(a,i12)' ) '  C = ', c
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             N             A             B'
  write ( *, '(a)' ) ' '

  do n = 0, 10
    call lcrg_anbn ( a, b, c, n, an, bn )
    write ( *, '(2x,i12,2x,i12,2x,i12)' ) n, an, bn
  end do

  return
end
subroutine test119 ( )

!*****************************************************************************80
!
!! TEST119 tests LCRG_ANBN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) an
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bn
  integer ( kind = 4 ) c
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) u
  integer ( kind = 4 ) v
  integer ( kind = 4 ), allocatable, dimension ( : ) :: x
  integer ( kind = 4 ), allocatable, dimension ( : ) :: y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST119'
  write ( *, '(a)' ) '  LCRG_ANBN determines a linear congruential random'
  write ( *, '(a)' ) '  number generator equivalent to N steps of a given one.'
!
!  These parameters define the old (1969) IBM 360 random number generator:
!
  a = 16807
  b = 0
  c = 2147483647

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LCRG parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  A  = ', a
  write ( *, '(a,i12)' ) '  B  = ', b
  write ( *, '(a,i12)' ) '  C  = ', c
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                           N            In           Out'
  write ( *, '(a)' ) ' '

  k = 0
  u = 12345
  write ( *, '(2x,12x,2x,i12,2x,12x,2x,i12)' ) k, u
  do k = 1, 11
    call lcrg_evaluate ( a, b, c, u, v )
    write ( *, '(2x,12x,2x,i12,2x,i12,2x,i12)' ) k, u, v
    u = v
  end do
!
!  Now try to replicate these results using N procesors.
!
  n = 4
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  call lcrg_anbn ( a, b, c, n, an, bn )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LCRG parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  AN = ', an
  write ( *, '(a,i12)' ) '  BN = ', bn
  write ( *, '(a,i12)' ) '  C  = ', c
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             J             N            In           Out'
  write ( *, '(a)' ) ' '

  x(1) = 12345
  do j = 2, n
    call lcrg_evaluate ( a, b, c, x(j-1), x(j) )
  end do

  do j = 1, n
    write ( *, '(2x,i12,2x,i12,2x,12x,2x,i12)' ) j, j-1, x(j)
  end do

  do k = n + 1, 12, n
    do j = 1, n
      call lcrg_evaluate ( an, bn, c, x(j), y(j) )
      write ( *, '(2x,i12,2x,i12,2x,i12,2x,i12)' ) j, k+j-2, x(j), y(j)
      x(j) = y(j)
    end do
  end do

  deallocate ( x )
  deallocate ( y )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests LCRG_SEED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_lcrg
  integer ( kind = 4 ) seed_out
  integer ( kind = 4 ) seed_start
  real ( kind = 4 ) u

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  LCRG_SEED directly computes the updated value of a'
  write ( *, '(a)' ) '  seed used by an linear congruential random number'
  write ( *, '(a)' ) '  generator.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I          SEED          SEED          SEED    U'
  write ( *, '(a)' ) '                 Input        Output          LCRG'
  write ( *, '(a)' ) ' '
!
!  These parameters define the old (1969) IBM 360 random number generator:
!
  a = 16807
  b = 0
  c = 2147483647
!
!  This seed value was used in Pierre L'Ecuyer's article.
!
  seed_start = 12345

  seed = seed_start
!
!  Compute 1000 random numbers "the hard way", that is, sequentially.
!  Every now and then, call LCRG_SEED to compute SEED directly.
!
  do i = 1, 1000

    seed_in = seed
    u = r4_uniform_01 ( seed )
    seed_out = seed

    if ( i <= 10 .or. i == 100 .or. i == 1000 ) then

      call lcrg_seed ( a, b, c, i, seed_start, seed_lcrg )

      write ( *, '(2x,i8,2x,i12,2x,i12,2x,i12,2x,g14.6)' ) &
        i, seed_in, seed_out, seed_lcrg, u

    end if

  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests R4_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_uniform_ab
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  R4_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'

  a = 5.0E+00
  b = 10.0E+00
  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
  write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6)' ) i, r4_uniform_ab ( a, b, seed )
  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests R4_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  R4_UNIFORM_01 computes pseudorandom values '
  write ( *, '(a)' ) '  in the interval [0,1].'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6)' ) i, r4_uniform_01 ( seed )
  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests R8_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  R8_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'

  a = 5.0D+00
  b = 10.0D+00
  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
  write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6)' ) i, r8_uniform_ab ( a, b, seed )
  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests R8_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  R8_UNIFORM_01 computes pseudorandom values '
  write ( *, '(a)' ) '  in the interval [0,1].'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6)' ) i, r8_uniform_01 ( seed )
  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests R8_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 1000

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_out
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) u_avg
  real ( kind = 8 ) u_var

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  R8_UNIFORM_01 computes a sequence of '
  write ( *, '(a)' ) '  uniformly distributed pseudorandom numbers.'
!
!  Start with a known seed.
!
  seed = 12345

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Initial SEED = ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 10 values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I         Input        Output   R8_UNIFORM_01'
  write ( *, '(a)' ) '                SEED          SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    seed_in = seed
    u(i) = r8_uniform_01 ( seed )
    seed_out = seed

    write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) i, seed_in, seed_out, u(i)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i10,a)' ) '  Now compute ', n, ' elements.'
  write ( *, '(a)' ) ' '

  u_avg = 0.0D+00
  do i = 1, n
    u(i) = r8_uniform_01 ( seed )
    u_avg = u_avg + u(i)
  end do

  u_avg = u_avg / real ( n, kind = 8 )

  u_var = sum ( ( u(1:n) - u_avg )**2 ) / real ( n - 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.6)' ) '  Average value = ', u_avg
  write ( *, '(a,f10.6)' ) '  Expecting       ', 0.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,f10.6)' ) '  Variance =      ', u_var
  write ( *, '(a,f10.6)' ) '  Expecting       ', 1.0D+00 / 12.0D+00

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests R8_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in
  integer ( kind = 4 ) seed_out
  integer ( kind = 4 ) seed_save
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  R8_UNIFORM_01 computes a sequence of pseudorandom'
  write ( *, '(a)' ) '  numbers but all computations depend on the seed value.'
  write ( *, '(a)' ) '  In this test, we show how a sequence of "random"'
  write ( *, '(a)' ) '  values can be manipulated by accessing the seed.'

  seed = 1066

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Initial SEED is ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call R8_UNIFORM_01 10 times, and watch SEED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I         Input        Output   R8_UNIFORM_01'
  write ( *, '(a)' ) '                SEED          SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    seed_in = seed

    if ( i == 5 ) then
      seed_save = seed
    end if

    x = r8_uniform_01 ( seed )
    seed_out = seed
    write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) i, seed_in, seed_out,x

  end do

  seed = seed_save

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Reset SEED to its value at step 5, = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call R8_UNIFORM_01 10 times, and watch how SEED'
  write ( *, '(a)' ) '  and R8_UNIFORM_01 restart themselves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I         Input        Output    R8_UNIFORM_01'
  write ( *, '(a)' ) '                SEED          SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    x = r8_uniform_01 ( seed )
    seed_out = seed
    write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) i, seed_in, seed_out,x
  end do

  seed = -12345678

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  What happens with a negative SEED?'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I         Input        Output   R8_UNIFORM_01'
  write ( *, '(a)' ) '                SEED          SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    x = r8_uniform_01 ( seed )
    seed_out = seed
    write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) i, seed_in, seed_out,x
  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests R8_UNIFORM_01 and R8MAT_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 100
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  R8_UNIFORM_01 computes pseudorandom values '
  write ( *, '(a)' ) '    one at a time.'
  write ( *, '(a)' ) '  R8MAT_UNIFORM_01 computes a matrix of values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For the same initial seed, the results should '
  write ( *, '(a)' ) '  be identical, but R8MAT_UNIFORM_01 might be faster.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

  seed = seed_init

  do j = 1, n
    do i = 1, m
      a(i,j) = r8_uniform_01 ( seed )
    end do
  end do

  seed = seed_init;
  call r8mat_uniform_01 ( m, n, seed, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I       J      A(I,J)        B(I,J)'
  write ( *, '(a)' ) '                  (R8_UNIFORM_01)  (R8MAT_UNIFORM_01)'
  write ( *, '(a)' ) ' '

  do k = 0, 10
    i = ( k * m + ( 10 - k ) * 1 ) / 10
    j = ( k * n + ( 10 - k ) * 1 ) / 10
    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) i, j, a(i,j), b(i,j)
  end do
  
  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests R8_UNIFORM_01 and R8VEC_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  R8_UNIFORM_01 computes pseudorandom values '
  write ( *, '(a)' ) '  one at a time.'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_01 computes a vector of values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For the same initial seed, the results should '
  write ( *, '(a)' ) '  be identical, but R8VEC_UNIFORM_01 might be faster.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

  seed = seed_init

  do i = 1, n
    a(i) = r8_uniform_01 ( seed )
  end do

  seed = seed_init;
  call r8vec_uniform_01 ( n, seed, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I      A(I)            B(I)'
  write ( *, '(a)' ) '          (R8_UNIFORM_01)  (R8VEC_UNIFORM_01)'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, a(i), b(i)
  end do
  
  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests R8I8_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r8i8_uniform_ab
  integer ( kind = 8 ) i
  integer ( kind = 8 ) seed
  integer ( kind = 8 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  R8I8_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B] using an I8 seed.'

  a = 5.0D+00
  b = 10.0D+00
  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
  write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6)' ) i, r8i8_uniform_ab ( a, b, seed )
  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests R8I8_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8i8_uniform_01
  integer ( kind = 8 ) i
  integer ( kind = 8 ) seed
  integer ( kind = 8 ), parameter :: seed_init = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  R8I8_UNIFORM_01 computes pseudorandom values '
  write ( *, '(a)' ) '  in the interval [0,1] using an I8 seed.'

  seed = seed_init

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6)' ) i, r8i8_uniform_01 ( seed )
  end do

  return
end
