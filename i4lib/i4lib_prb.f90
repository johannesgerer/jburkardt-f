program main

!*****************************************************************************80
!
!! MAIN is the main program for I4LIB_PRB.
!
!  Discussion:
!
!    I4LIB_PRB tests routines from the I4LIB library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4LIB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the I4LIB library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
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
  call test23 ( )
  call test24 ( ) 
  call test245 ( )
  call test25 ( )
  call test26 ( )
  call test27 ( )
  call test28 ( )
  call test29 ( )

  call test30 ( )
  call test31 ( )
  call test32 ( )
  call test33 ( )
  call test34 ( )
  call test35 ( )
  call test36 ( )
  call test37 ( )
  call test38 ( )
  call test39 ( )

  call test40 ( )
  call test41 ( )
  call test42 ( )
  call test43 ( )
  call test44 ( )
  call test45 ( )
  call test46 ( )
  call test47 ( )
  call test48 ( )
  call test49 ( )

  call test50 ( )
  call test51 ( )
  call test52 ( )
  call test53 ( )
  call test54 ( )
  call test55 ( )
  call test56 ( )
  call test57 ( )
  call test58 ( )
  call test59 ( )

  call test60 ( )
  call test602 ( )
  call test605 ( )
  call test61 ( )
  call test62 ( )
  call test63 ( )
  call test64 ( )
  call test65 ( )
  call test66 ( )
  call test67 ( )
  call test68 ( )
  call test69 ( )

  call test70 ( )
  call test71 ( )
  call test72 ( )
  call test73 ( )
  call test74 ( )
  call test75 ( )
  call test76 ( )
  call test77 ( )
  call test78 ( )
  call test79 ( )

  call test80 ( )
  call test81 ( )
  call test82 ( )
  call test83 ( )
  call test84 ( )
  call test85 ( )
  call test86 ( )
  call test87 ( )
  call test88 ( )
  call test89 ( )

  call test90 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4LIB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests I4_BIT_HI1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) i4_bit_hi1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  I4_BIT_HI1 returns the location of the high 1 bit.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I  I4_BIT_HI1(I)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    i = i4_uniform_ab ( 0, 100, seed )
    j = i4_bit_hi1 ( i )
    write ( *, '(2x,i8,2x,i8)' ) i, j
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests I4_BIT_LO0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) i4_bit_lo0
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  I4_BIT_LO0 returns the location of the low 0 bit.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I  I4_BIT_LO0(I)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    i = i4_uniform_ab ( 0, 100, seed )
    j = i4_bit_lo0 ( i )
    write ( *, '(2x,i8,2x,i8)' ) i, j
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests I4_BIT_LO1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) i4_bit_lo1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  I4_BIT_LO1 returns the location of the low 1 bit.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I  I4_BIT_LO1(I)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    i = i4_uniform_ab ( 0, 100, seed )
    j = i4_bit_lo1 ( i )
    write ( *, '(2x,i8,2x,i8)' ) i, j
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests I4_BIT_REVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i4_bit_reverse
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  I4_BIT_REVERSE bit reverses I with respect to 2^J'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J  I4_BIT_REVERSE(I,J)'
  write ( *, '(a)' ) ' '

  do j = 0, 4
    i_hi = 2**j - 1
    do i = 0, i_hi
      k = i4_bit_reverse ( i, j )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
    end do
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests I4_CHARACTERISTIC.
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_characteristic

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  I4_CHARACTERISTIC computes the characteristic'
  write ( *, '(a)' ) '  of an integer Q, which is  '
  write ( *, '(a)' ) '    Q if Q is prime;'
  write ( *, '(a)' ) '    P, if Q = P**N for some prime P;'
  write ( *, '(a)' ) '    0, if Q is negative, 0, 1, or the product of '
  write ( *, '(a)' ) '      more than 1 distinct prime.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  I4_CHARACTERISTIC'
  write ( *, '(a)' ) ' '

  do i = 1, 50
    write ( *, '(2x,i2,13x,i4)' ) i, i4_characteristic ( i )
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests I4_DIV_ROUNDED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ), parameter :: a_hi =  100
  integer ( kind = 4 ), parameter :: a_lo = -100
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: b_hi =  10
  integer ( kind = 4 ), parameter :: b_lo = -10
  real ( kind = 8 ) c0
  integer ( kind = 4 ) c1
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) c3
  integer ( kind = 4 ) c4
  integer ( kind = 4 ) i4_div_rounded
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  I4_DIV_ROUNDED performs rounded integer division.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C0 = real ( a ) / real ( b )'
  write ( *, '(a)' ) '  C1 = I4_DIV_ROUNDED ( A, B )'
  write ( *, '(a)' ) '  C2 = nint ( real ( a ) / real ( b ) )'
  write ( *, '(a)' ) '  C3 = A / B'
  write ( *, '(a)' ) '  C4 = int ( real ( a ) / real ( b ) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C1 and C2 should be equal;'
  write ( *, '(a)' ) '  C3 and C4 should be equal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     A     B           C0         C1    C2      C3    C4'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num
    a = i4_uniform_ab ( a_lo, a_hi, seed )
    b = i4_uniform_ab ( b_lo, b_hi, seed )
    if ( b == 0 ) then
      b = 7
    end if
    c0 = real ( a, kind = 8 ) / real ( b, kind = 8 )
    c1 = i4_div_rounded ( a, b )
    c2 = nint ( real ( a, kind = 8 ) / real ( b, kind = 8 ) )
    c3 = a / b
    c4 = int ( real ( a, kind = 8 ) / real ( b, kind = 8 ) )
    write ( *, '(2x,i4,2x,i4,4x,f14.6,2x,i4,2x,i4,4x,i4,2x,i4)' ) &
      a, b, c0, c1, c2, c3, c4
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests I4_DIVP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ), parameter :: a_hi =  100
  integer ( kind = 4 ), parameter :: a_lo = -100
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: b_hi =  10
  integer ( kind = 4 ), parameter :: b_lo = -10
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i4_divp
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  I4_DIVP(A,B) returns the smallest multiplier of J'
  write ( *, '(a)' ) '  that is less than I'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     A     B     C     D'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num
    a = i4_uniform_ab ( a_lo, a_hi, seed )
    b = i4_uniform_ab ( b_lo, b_hi, seed )
    if ( b == 0 ) then
      b = 7
    end if
    c = i4_divp ( a, b )
    d = c * b
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) a, b, c, d
  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests I4_GCD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ), dimension(test_num) :: i_test = (/ &
    36, 49, 0, 12, 36, 1, 91 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension(test_num) :: j_test = (/ &
    30, -7, 71, 12, 49, 42, 28 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  I4_GCD computes the greatest common factor,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I       J  I4_GCD'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    i = i_test(test)
    j = j_test(test)
    write ( *, '(2x,3i8)') i, j, i4_gcd ( i, j )
  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests I4_HUGE.
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

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) dummy

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  I4_HUGE returns a huge integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  I4_HUGE() = ', i4_huge ( )
  write ( *, '(a,i12)' ) '  HUGE(1) =   ', huge ( dummy )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests I4_HUGE_NORMALIZER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge
  real ( kind = 8 ) i4_huge_normalizer
  real ( kind = 8 ) r8
  real ( kind = 8 ) value

  i4 = i4_huge ( )
  r8 = i4_huge_normalizer ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  I4_HUGE_NORMALIZER returns 1/(I4_HUGE+1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  I4_HUGE() = ', i4
  write ( *, '(a,g14.6)' ) '  I4_HUGE_NORMALIZER() = ', r8

  value = real ( i4, kind = 8 ) * r8

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  I4_HUGE * I4_HUGE_NORMALIZER = ', value

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests I4_IS_PRIME.
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

  integer ( kind = 4 ) i
  logical i4_is_prime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  I4_IS_PRIME reports whether an integer is prime.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I     I4_IS_PRIME(I)'
  write ( *, '(a)' ) ' '

  do i = -2, 25
    write ( *, '(2x,i8,2x,l1)' ) i, i4_is_prime ( i )
  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests I4_LCM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_lcm
  integer ( kind = 4 ), dimension(test_num) :: i_test = (/ &
    36, 49, 0, 12, 36, 1, 91 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension(test_num) :: j_test = (/ &
    30, -7, 71, 12, 49, 42, 28 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  I4_LCM computes the least common multiple.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I     J   I4_LCM'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    i = i_test(test)
    j = j_test(test)
    write ( *, '(2x,4i8)') i, j, i4_lcm ( i, j )
  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests I4_LOG_10.
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

  integer ( kind = 4 ), parameter :: test_num = 13

  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x
  integer ( kind = 4 ), dimension ( test_num ) :: x_test = (/ &
    0, 1, 2, 3, 9, 10, 11, 99, 101, -1, -2, -3, -9 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  I4_LOG_10: whole part of log base 10,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, I4_LOG_10'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = x_test(test)
    write ( *, '( 2x, i8, i12 )' ) x, i4_log_10 ( x )
  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests I4_LOG_2.
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

  integer ( kind = 4 ), parameter :: test_num = 17

  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x
  integer ( kind = 4 ), dimension ( test_num ) :: x_test = (/ &
      0,    1,    2,    3,    9, &
     10,   11,   99,  101,   -1, &
     -2,   -3,   -9, 1000, 1023, &
   1024, 1025 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  I4_LOG_2: whole part of log base 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X     I4_LOG_2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = x_test(test)
    write ( *, '( 2x, i8, i12 )' ) x, i4_log_2 ( x )
  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests I4_LOG_I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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

  integer i4
  integer i4_log_i4
  integer j4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
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
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests I4_LOG_R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 10

  real ( kind = 8 ) b
  real ( kind = 8 ), dimension(test_num) :: b_test = (/ &
    2.0D+00, 3.0D+00,  4.0D+00,  5.0D+00,   6.0D+00, &
    7.0D+00, 8.0D+00, 16.0D+00, 32.0D+00, 256.0D+00 /)
  integer ( kind = 4 ) i4_log_r8
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x

  x = 16

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  I4_LOG_R8: whole part of log base B,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, B, I4_LOG_R8'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    b = b_test(test)

    write ( *, '(2x, i8, g14.6, i12 )' ) x, b, i4_log_r8 ( x, b )

  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests I4_MANT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) s
  real ( kind = 8 ) x

  x = -314.159D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  I4_MANT decomposes an integer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Number to be decomposed is X = ', x

  call i4_mant ( x, s, j, k, l )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i24,a,i24,a,i24,a,i8)' ) &
    '  I4_MANT: X = ', s, ' * (', j, '/', k, ') * 2**', l

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests I4_MODDIV;
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

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ), dimension ( test_num ) :: ndivid = (/ &
    50, -50, 50, -50 /)
  integer ( kind = 4 ) nmult
  integer ( kind = 4 ) nrem
  integer ( kind = 4 ), dimension ( test_num ) :: number = (/ &
    107, 107, -107, -107 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  I4_MODDIV factors a number'
  write ( *, '(a)' ) '  into a multiple and a remainder.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Number   Divisor  Multiple Remainder'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call i4_moddiv ( number(test), ndivid(test), nmult, nrem )
    write ( *, '(2x,4i10)' ) number(test), ndivid(test), nmult, nrem
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat using FORTRAN MOD:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    nrem = mod ( number(test), ndivid(test) )
    nmult = number(test) / ndivid(test)
    write ( *, '(2x,4i10)' ) number(test), ndivid(test), nmult, nrem
  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests I4_MODP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ), dimension ( test_num ) :: ndivid = (/ &
    50, -50, 50, -50 /)
  integer ( kind = 4 ) nmult
  integer ( kind = 4 ) nrem
  integer ( kind = 4 ), dimension ( test_num ) :: number = (/ &
    107, 107, -107, -107 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  I4_MODP factors a number'
  write ( *, '(a)' ) '  into a multiple and a remainder.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Number   Divisor  Multiple Remainder'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    nrem = i4_modp ( number(test), ndivid(test) )
    nmult = ( number(test) - nrem ) / ndivid(test)
    write ( *, '(2x,4i10)' ) number(test), ndivid(test), nmult, nrem
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat using FORTRAN MOD:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    nrem = mod ( number(test), ndivid(test) )
    nmult = number(test) / ndivid(test)
    write ( *, '(2x,4i10)' ) number(test), ndivid(test), nmult, nrem
  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests I4_SIGN.
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

  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x
  integer ( kind = 4 ), dimension(test_num) :: x_test = (/ -10, -7, 0, 5, 9 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  I4_SIGN returns the sign of a number.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = x_test(test)
    write ( *, '(2i8)' ) x, i4_sign ( x )
  end do

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests I4_SWAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  I4_SWAP swaps two integers.'

  i = 1
  j = 202

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Before swapping: '
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  I = ', i
  write ( *, '(a,i8)' ) '  J = ', j

  call i4_swap ( i, j )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After swapping: '
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  I = ', i
  write ( *, '(a,i8)' ) '  J = ', j

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests I4_WALSH_1D;
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_walsh_1d
  integer ( kind = 4 ) w0
  integer ( kind = 4 ) wm1
  integer ( kind = 4 ) wm2
  integer ( kind = 4 ) wm3
  integer ( kind = 4 ) wp1
  integer ( kind = 4 ) wp2
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  I4_WALSH_1D evaluates 1D Walsh functions:'
  write ( *, '(a)' ) ' '
  write ( *, * ) 'X  W(+2) W(+1) W(0) W(-1) W(-2) W(-3)'
  write ( *, '(a)' ) ' '

  do i = 0, 32

    x = real ( i, kind = 8 ) / 4.0D+00

    wp2 = i4_walsh_1d ( x,  2 )
    wp1 = i4_walsh_1d ( x,  1 )
    w0  = i4_walsh_1d ( x,  0 )
    wm1 = i4_walsh_1d ( x, -1 )
    wm2 = i4_walsh_1d ( x, -2 )
    wm3 = i4_walsh_1d ( x, -3 )

    write ( *, '(2x,f10.6,6i2)' ) x, wp2, wp1, w0, wm1, wm2, wm3

  end do

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests I4_WRAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo

  ilo = 4
  ihi = 8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  I4_WRAP forces an integer to lie within given limits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ILO = ', ilo
  write ( *, '(a,i8)' ) '  IHI = ', ihi
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  I4_WRAP(I)'
  write ( *, '(a)' ) ' '

  do i = -10, 20
    write ( *, '(2x,2i8)' ) i, i4_wrap ( i, ilo, ihi )
  end do

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests I4_XOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: i4_hi = 100
  integer ( kind = 4 ) :: i4_lo = 0
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) i4_xor
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  I4_XOR returns the bitwise exclusive OR of'
  write ( *, '(a)' ) '  two integers.'
  write ( *, '(a)' ) '  Compare with FORTRAN90 intrinsic IEOR.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J    I4_XOR      IEOR'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    i = i4_uniform_ab ( i4_lo, i4_hi, seed )
    j = i4_uniform_ab ( i4_lo, i4_hi, seed )
    k = i4_xor ( i, j )
    l = ieor ( i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i, j, k, l
  end do

  return
end
subroutine test245 ( )

!*****************************************************************************80
!
!! TEST245 tests I4BLOCK_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 4
  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), dimension(l,m,n) :: x = reshape ( (/ &
        1,  2,  3,   4,  1, &
        4,  9, 16,   1,  8, &
       27, 64,  2,   4,  6, &
        8,  2,  8,  18, 32, &
        2, 16, 54, 128 /), (/ l, m, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST245'
  write ( *, '(a)' ) '  I4BLOCK_PRINT prints an I4BLOCK.'

  call i4block_print ( l, m, n, x, '  The 3D array:' )

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests I4COL_FIND_ITEM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item
  integer ( kind = 4 ), dimension ( test_num ) :: item_test = (/ &
    34, 12, 90 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  I4COL_FIND_ITEM finds the first occurrence of'
  write ( *, '(a)' ) '  an item in an integer array of columns.'

  do i = 1, m
    do j = 1, n
      a(i,j) = 10 * i + j
    end do
  end do

  call i4mat_print ( m, n, a, '  The matrix of columns:' )

  do test = 1, test_num

    item = item_test(test)

    call i4col_find_item ( m, n, a, item, row, col )

    write ( *, '(a,i8,a,i8,a,i8)' ) '  Item ', item, '  occurs in row ', &
      row, ' and column ', col

  end do

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests I4COL_FIND_PAIR_WRAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item1
  integer ( kind = 4 ), dimension ( test_num ) :: item1_test = (/ &
    22, 32, 22, 54, 54 /)
  integer ( kind = 4 ) item2
  integer ( kind = 4 ), dimension ( test_num ) :: item2_test = (/ &
    32, 22, 23, 14, 11 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  I4COL_FIND_PAIR_WRAP finds the first occurrence of'
  write ( *, '(a)' ) '  a pair of item in an integer array of columns.'
  write ( *, '(a)' ) '  Items in the array are ordered by column, and'
  write ( *, '(a)' ) '  wraparound is allowed.'

  do i = 1, m
    do j = 1, n
      a(i,j) = 10 * i + j
    end do
  end do

  call i4mat_print ( m, n, a, '  The matrix of columns:' )

  do test = 1, test_num

    item1 = item1_test(test)
    item2 = item2_test(test)

    call i4col_find_pair_wrap ( m, n, a, item1, item2, row, col )

    write ( *, '(a,i8,a,i8,a,i8,a,i8)' ) '  Item ', item1, &
      ' followed by item ', item2, ' occurs in row ', &
      row, ' and column ', col

  end do

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests I4COL_SORT_A and I4COL_SORT_D.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) :: b = 1
  integer ( kind = 4 ) :: c = 10
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  I4COL_SORT_A ascending sorts an integer array'
  write ( *, '(a)' ) '  as a table of columns.'
  write ( *, '(a)' ) '  I4COL_SORT_D descending sorts an integer array'
  write ( *, '(a)' ) '  as a table of columns.'

  seed = 123456789

  call i4mat_uniform_ab ( m, n, b, c, seed, a )

  call i4mat_print ( m, n, a, '  The original matrix:' )

  call i4col_sort_a ( m, n, a )

  call i4mat_print ( m, n, a, '  Ascending sorted:' )

  call i4col_sort_d ( m, n, a )

  call i4mat_print ( m, n, a, '  Descending sorted:' )

  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests I4COL_SORT2_A;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = 20
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  For a rectangular integer matrix:'
  write ( *, '(a)' ) '  I4COL_SORT2_D sorts the elements of the columns.'

  call i4mat_uniform_ab ( m, n, b, c, seed, a )

  call i4mat_print ( m, n, a, '  The matrix:' )

  call i4col_sort2_a ( m, n, a )

  call i4mat_print ( m, n, a, '  The element-sorted column matrix:' )

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests I4COL_SORTED_SINGLETON_COUNT;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) singleton_num
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  I4COL_SORTED_SINGLETON_COUNT counts singletons'
  write ( *, '(a)' ) '  in a sorted I4COL;'

  seed = 123456789

  do test = 1, test_num

    b = 0
    c = 3

    call i4mat_uniform_ab ( m, n, b, c, seed, a )

    call i4col_sort_a ( m, n, a )

    call i4mat_print ( m, n, a, '  Ascending sorted ICOL:' )

    call i4col_sorted_singleton_count ( m, n, a, singleton_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) '  Number of singletons = ', singleton_num

  end do

  return
end
subroutine test30 ( )

!*****************************************************************************80
!
!! TEST30 tests I4COL_SORTED_UNIQUE_COUNT;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  I4COL_SORTED_UNIQUE_COUNT counts the unique entries'
  write ( *, '(a)' ) '  of a sorted I4COL;'

  seed = 123456789

  do test = 1, test_num

    b = 0
    c = 3

    call i4mat_uniform_ab ( m, n, b, c, seed, a )

    call i4col_sort_a ( m, n, a )

    call i4mat_print ( m, n, a, '  Ascending sorted I4COL:' )

    call i4col_sorted_unique_count ( m, n, a, unique_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) '  Number of unique entries = ', unique_num

  end do

  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests I4MAT_ELIM and I4MAT_RED.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(n)
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(m)
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  I4MAT_ELIM does exact Gauss elimination.'
  write ( *, '(a)' ) '  I4MAT_RED divides common factors in a matrix;'

  do test = 1, test_num

    if ( test == 1 ) then

      k = 0
      do i = 1, m
        do j = 1, n
          k = k + 1
          a(i,j) = k
        end do
      end do

    else if ( test == 2 ) then

      factor = 8 * 7 * 6 * 5 * 4 * 3 * 2

      do i = 1, m
        do j = 1, n
          a(i,j) = factor / ( i + j - 1 )
        end do
      end do

    else if ( test == 3 ) then

      do i = 1, m
        do j = 1, n
          a(i,j) = i * j
        end do
      end do

    end if

    call i4mat_print ( m, n, a, '  The original matrix:' )

    call i4mat_red ( m, n, a, row, col )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The matrix, as returned by I4MAT_RED:'
    write ( *, '(a)' ) '  (Factors are displayed in an extra row and column.'
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(2x,6i8)' )  ( a(i,j), j = 1, n ), row(i)
    end do
    write ( *, '(2x,5i8)' )  ( col(j), j = 1, n )

    call i4mat_elim ( m, n, a )

    call i4mat_print ( m, n, a, '  The matrix returned by I4MAT_ELIM:' )

  end do

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests I4MAT_MAX_INDEX and I4MAT_MIN_INDEX.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = 10
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  I4MAT_MAX_INDEX locates the maximum;'
  write ( *, '(a)' ) '  I4MAT_MIN_INDEX locates the minimum;'

  seed = 123456789

  call i4mat_uniform_ab ( m, n, b, c, seed, a )

  call i4mat_print ( m, n, a, '  Random array:' )

  write ( *, '(a)' ) ' '
  call i4mat_max_index ( m, n, a, i, j )
  write ( *, '(a,2i8)' ) '  Maximum I,J indices            ', i, j
  call i4mat_min_index ( m, n, a, i, j )
  write ( *, '(a,2i8)' ) '  Minimum I,J indices            ', i, j

  return
end
subroutine test33 ( )

!*****************************************************************************80
!
!! TEST33 tests I4MAT_L1_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6
!
!  Each row of this definition is a COLUMN of the matrix.
!
  integer ( kind = 4 ), dimension ( n, n ) :: a = reshape ( (/ &
     1,  2,  0,  5,  0, 75, &
     0,  1,  0,  0,  0,  0, &
     0,  0,  1,  3,  0,  0, &
     0,  0,  0,  1,  0,  6, &
     0,  0,  0,  0,  1,  4, &
     0,  0,  0,  0,  0,  1 /), (/ n, n /) )
  integer ( kind = 4 ) b(n,n)
  integer ( kind = 4 ) c(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  I4MAT_L1_INVERSE inverts a unit lower triangular matrix.'

  call i4mat_print ( n, n, a, '  The original matrix:' )

  call i4mat_l1_inverse ( n, a, b )

  call i4mat_print ( n, n, b, '  The inverse matrix:' )

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call i4mat_print ( n, n, c, '  The product:' )

  return
end
subroutine test34 ( )

!*****************************************************************************80
!
!! TEST34 tests I4MAT_PERM_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  I4MAT_PERM_UNIFORM applies a random permutation'
  write ( *, '(a)' ) '  to a square integer matrix.'

  seed = 123456789

  do i = 1, n
    do j = 1, n
      a(i,j) = 10 * i + j
    end do
  end do

  call i4mat_print ( n, n, a, '  The original matrix:' )

  call i4mat_perm_uniform ( n, a, seed )

  call i4mat_print ( n, n, a, '  The permuted matrix:' )

  return
end
subroutine test35 ( )

!*****************************************************************************80
!
!! TEST35 tests I4MAT_U1_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6
!
!  Each row of this definition is a COLUMN of the matrix.
!
  integer ( kind = 4 ), dimension ( n, n ) :: a = reshape ( (/ &
    1,  0,  0,  0,  0,  0, &
    2,  1,  0,  0,  0,  0, &
    0,  0,  1,  0,  0,  0, &
    5,  0,  3,  1,  0,  0, &
    0,  0,  0,  0,  1,  0, &
   75,  0,  0,  6,  4,  1 /), (/ n, n /) )
  integer ( kind = 4 ) b(n,n)
  integer ( kind = 4 ) c(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  I4MAT_U1_INVERSE inverts a unit upper triangular matrix.'

  call i4mat_print ( n, n, a, '  The original matrix:' )

  call i4mat_u1_inverse ( n, a, b )

  call i4mat_print ( n, n, b, '  The inverse matrix:' )

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call i4mat_print ( n, n, c, '  The product:' )

  return
end
subroutine test36 ( )

!*****************************************************************************80
!
!! TEST36 tests I4ROW_MAX and I4ROW_MIN;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) amax(m)
  integer ( kind = 4 ) amin(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36'
  write ( *, '(a)' ) '  I4ROW_MAX computes row maximums;'
  write ( *, '(a)' ) '  I4ROW_MIN computes row minimums;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = k
    end do
  end do

  call i4mat_print ( m, n, a, '  The matrix:' )

  call i4row_max ( m, n, a, amax )

  call i4row_min ( m, n, a, amin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Maximum, minimum:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,i3,3x,2i8)' ) &
      i, amax(i), amin(i)
  end do

  return
end
subroutine test37 ( )

!*****************************************************************************80
!
!! TEST37 tests I4ROW_MEAN and I4ROW_SUM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) mean(m)
  integer ( kind = 4 ) rowsum(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  I4ROW_MEAN computes row means;'
  write ( *, '(a)' ) '  I4ROW_SUM computes row sums;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = k
    end do
  end do

  call i4mat_print ( m, n, a, '  The matrix:' )

  call i4row_sum ( m, n, a, rowsum )

  call i4row_mean ( m, n, a, mean )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sum, mean:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,i3,3x,i8,3x,f10.4)' ) &
      i, rowsum(i), mean(i)
  end do

  return
end
subroutine test38 ( )

!*****************************************************************************80
!
!! TEST38 tests I4ROW_SORT_A;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = 10
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST38'
  write ( *, '(a)' ) '  For a rectangular integer matrix:'
  write ( *, '(a)' ) '  I4ROW_SORT_A sorts the rows;'

  seed = 123456789

  call i4mat_uniform_ab ( m, n, b, c, seed, a )

  call i4mat_print ( m, n, a, '  The original matrix:' )

  call i4row_sort_a ( m, n, a )

  call i4mat_print ( m, n, a, '  The row-sorted matrix:' )

  return
end
subroutine test39 ( )

!*****************************************************************************80
!
!! TEST39 tests I4ROW_SORT_D and I4ROW_SORT2_D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST39'
  write ( *, '(a)' ) '  For a rectangular integer matrix:'
  write ( *, '(a)' ) '  I4ROW_SORT_D sorts the rows;'
  write ( *, '(a)' ) '  I4ROW_SORT2_D sorts the elements of the rows.'

  seed = 123456789

  do i = 1, m
    do j = 1, n
      a(i,j) = 10 * i + j
    end do
  end do

  call i4mat_print ( m, n, a, '  The original matrix:' )

  call i4mat_perm2_uniform ( m, n, a, seed )

  call i4mat_print ( m, n, a, '  The matrix, permuted by I4MAT_PERM2_UNIFORM:' )

  call i4row_sort_d ( m, n, a )

  call i4mat_print ( m, n, a, '  The row-sorted matrix:' )

  call i4row_sort2_d ( m, n, a )

  call i4mat_print ( m, n, a, '  The element-sorted row-sorted matrix:' )

  return
end
subroutine test40 ( )

!*****************************************************************************80
!
!! TEST40 tests I4ROW_SWAP;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST40'
  write ( *, '(a)' ) '  For an integer matrix of rows,'
  write ( *, '(a)' ) '  I4ROW_SWAP swaps two rows;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = k
    end do
  end do

  call i4mat_print ( m, n, a, '  The matrix:' )

  row1 = 1
  row2 = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Swap rows ', row1, ' and ', row2
  write ( *, '(a)' ) ' '

  call i4row_swap ( m, n, a, row1, row2 )

  call i4mat_print ( m, n, a, '  The new matrix:' )

  return
end
subroutine test41 ( )

!*****************************************************************************80
!
!! TEST41 tests I4ROW_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) variance(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST41'
  write ( *, '(a)' ) '  I4ROW_VARIANCE computes row variances;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = k
    end do
  end do

  call i4mat_print ( m, n, a, '  The matrix:' )

  call i4row_variance ( m, n, a, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Row variances:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,i3,2x,f10.4)' ) i, variance(i)
  end do

  return
end
subroutine test42 ( )

!*****************************************************************************80
!
!! TEST42 tests I4VEC_AMAX;
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

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST42'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_AMAX:   maximum absolute entry;'

  seed = 123456789

  b = -n
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  call i4vec_amax ( n, a, aval )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Maximum absolute value: ', aval

  return
end
subroutine test43 ( )

!*****************************************************************************80
!
!! TEST43 tests I4VEC_AMIN.
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

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST43'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_AMIN:   minimum absolute entry;'
  write ( *, '(a)' ) ' '

  seed = 123456789

  b = -n
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  write ( *, '(a)' ) ' '

  call i4vec_amin ( n, a, aval )

  write ( *, '(a,i8)' ) '  Minimum absolute value: ', aval

  return
end
subroutine test44 ( )

!*****************************************************************************80
!
!! TEST44 tests I4VEC_AMINZ and I4VEC_AMINZ_INDEX;
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

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST44'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_AMINZ:  minimum nonzero absolute entry;'
  write ( *, '(a)' ) '  I4VEC_AMINZ_INDEX: index of minimum nonzero absolute entry;'

  seed = 123456789

  b = -n
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  write ( *, '(a)' ) ' '

  call i4vec_aminz ( n, a, aval )
  call i4vec_aminz_index ( n, a, ival )

  write ( *, '(a,i8)' ) '  Minimum abs nonzero:      ', aval
  write ( *, '(a,i8)' ) '  Minimum abs nonzero index:', ival

  return
end
subroutine test45 ( )

!*****************************************************************************80
!
!! TEST45 tests I4VEC_AMAX_INDEX and I4VEC_AMIN_INDEX;
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

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST45'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_AMAX_INDEX:  index of maximum absolute entry;'
  write ( *, '(a)' ) '  I4VEC_AMIN_INDEX:  index minimum absolute entry;'

  seed = 123456789

  b = -n
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  write ( *, '(a)' ) ' '

  call i4vec_amax_index ( n, a, ival )

  write ( *, '(a,i8)' ) '  Maximum abs index:        ', ival

  call i4vec_amin_index ( n, a, ival )

  write ( *, '(a,i8)' ) '  Minimum abs index:      ', ival

  return
end
subroutine test46 ( )

!*****************************************************************************80
!
!! TEST46 tests I4VEC_MAX_INDEX, I4VEC_MAX_INDEX_LAST and I4VEC_INDEX;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) i4vec_max_index_last
  integer ( kind = 4 ) i4vec_index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST46'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_MAX_INDEX:          a maximal index;'
  write ( *, '(a)' ) '  I4VEC_MAX_INDEX_LAST:     last maximal index;'
  write ( *, '(a)' ) '  I4VEC_INDEX:              first index of given value;'

  seed = 123456789

  b = -n
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  write ( *, '(a)' ) ' '

  call i4vec_max_index ( n, a, ival )
  write ( *, '(a,i8)' ) '  Maximum index:            ', ival
  ival = i4vec_max_index_last ( n, a )
  write ( *, '(a,i8)' ) '  Last maximum index:       ', ival

  call i4vec_min_index ( n, a, ival )
  write ( *, '(a,i8)' ) '  Minimum index:            ', ival

  aval = a(n/2)
  write ( *, '(a)' ) ' '
  j = i4vec_index ( n, a, aval )
  write ( *, '(a,i8,a,i8)' ) '  Index of first occurrence of ', aval, ' is ', j

  aval = aval + 1
  j = i4vec_index ( n, a, aval )
  write ( *, '(a,i8,a,i8)' ) '  Index of first occurrence of ', aval, ' is ', j

  return
end
subroutine test47 ( )

!*****************************************************************************80
!
!! TEST47 tests I4VEC_ASCEND_SUB
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

  integer ( kind = 4 ), parameter :: n = 14

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) :: b = 1
  integer ( kind = 4 ) :: c = 10
  integer ( kind = 4 ) length
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) sub(n)
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST47'
  write ( *, '(a)' ) '  I4VEC_ASCEND_SUB computes a longest ascending'
  write ( *, '(a)' ) '  subsequence of an integer vector.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  do test = 1, test_num
    call i4vec_uniform_ab ( n, b, c, seed, a )
    call i4vec_print ( n, a, '  The vector to be tested:' )
    call i4vec_ascend_sub ( n, a, length, sub )
    call i4vec_print ( length, sub, '  A longest ascending subsequence:' )
  end do

  return
end
subroutine test48 ( )

!*****************************************************************************80
!
!! TEST48 tests I4VEC_ASCENDS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 6

  logical i4vec_ascends
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x(n)
!
!  Each ROW of this definition is a COLUMN of the matrix.
!
  integer ( kind = 4 ), dimension(n,test_num) :: x_test = reshape ( (/ &
    1, 3, 2, 4, &
    2, 2, 2, 2, &
    1, 2, 2, 4, &
    1, 2, 3, 4, &
    4, 4, 3, 1, &
    9, 7, 3, 0 /), (/ n, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST48'
  write ( *, '(a)' ) '  I4VEC_ASCENDS determines if an integer vector ascends.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    x(1:n) = x_test(1:n,test)

    call i4vec_print ( n, x, '  Test vector:' )

    write ( *, '(a,l1)' ) '  I4VEC_ASCENDS =  ', i4vec_ascends ( n, x )

  end do

  return
end
subroutine test49 ( )

!*****************************************************************************80
!
!! TEST49 tests I4VEC_BRACKET and I4VEC_INSERT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) a(n_max)
  integer ( kind = 4 ), dimension (test_num) :: atest = (/ &
    -10, 2, 9, 10, 20, 24 /)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) n
  integer ( kind = 4 ) right
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST49'
  write ( *, '(a)' ) '  I4VEC_BRACKET finds a pair of entries in a'
  write ( *, '(a)' ) '  sorted integer array which bracket a value.'
  write ( *, '(a)' ) '  I4VEC_INSERT inserts a value into a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use these two routines to bracket a value,'
  write ( *, '(a)' ) '  and then insert it.'

  n = 10
  do i = 1, n
    a(i) = 2 * i
  end do
  a(6) = a(5)

  call i4vec_print ( n, a, '  Sorted array:' )

  do test = 1, test_num

    aval = atest(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Search for AVAL = ', aval

    call i4vec_bracket ( n, a, aval, left, right )

    write ( *, '(a,i8)' ) '  Left = ', left
    write ( *, '(a,i8)' ) '  Right = ', right

    if ( 1 <= left ) then
      write ( *, '(a,i8)' ) '  A(LEFT)=', a(left)
    end if

    if ( 1 <= right ) then
      write ( *, '(a,i8)' ) '  A(RIGHT) = ', a(right)
    end if
!
!  Insert the value.
!
    if ( left == -1 ) then
      left = 0
    end if

    if ( left == right ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No insertion necessary.'

    else

      call i4vec_insert ( n, a, left+1, aval )

      n = n + 1

      call i4vec_print ( n, a, '  Sorted, augmented array:' )

    end if

  end do

  return
end
subroutine test50 ( )

!*****************************************************************************80
!
!! TEST50 tests I4VEC_CUM and I4VEC_CUM0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_cum(n)
  integer ( kind = 4 ) a_cum0(0:n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST50'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_CUM:   cumulative sum;'
  write ( *, '(a)' ) '  I4VEC_CUM0:  cumulative sum, zero based;'

  seed = 123456789

  b = -n
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  call i4vec_cum ( n, a, a_cum )

  call i4vec_print ( n, a_cum, '  Cumulative sums:' )

  call i4vec_cum0 ( n, a, a_cum0 )

  call i4vec_print ( n + 1, a_cum0, '  0-based Cumulative sums:' )

  return
end
subroutine test51 ( )

!*****************************************************************************80
!
!! TEST51 tests I4VEC_DESCENDS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 6

  logical i4vec_descends
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x(n)
!
!  Each ROW of this definition is a COLUMN of the matrix.
!
  integer ( kind = 4 ), dimension(n,test_num) :: x_test = reshape ( (/ &
    1, 3, 2, 4, &
    2, 2, 2, 2, &
    1, 2, 2, 4, &
    1, 2, 3, 4, &
    4, 4, 3, 1, &
    9, 7, 3, 0 /), (/ n, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST51'
  write ( *, '(a)' ) '  I4VEC_DESCENDS determines if an integer vector descends.'

  do test = 1, test_num

    x(1:n) = x_test(1:n,test)

    call i4vec_print ( n, x, '  Test vector:' )

    write ( *, '(a,l1)' ) '  I4VEC_DESCENDS = ', i4vec_descends ( n, x )

  end do

  return
end
subroutine test52 ( )

!*****************************************************************************80
!
!! TEST52 tests I4VEC_DIRECT_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: factor_num = 3
  integer ( kind = 4 ), parameter :: point_num = 24

  integer ( kind = 4 ) factor_index
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ), allocatable, dimension ( : ) :: factor_value
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x(factor_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST52'
  write ( *, '(a)' ) '  I4VEC_DIRECT_PRODUCT forms the entries of a'
  write ( *, '(a)' ) '  direct product of a given number of I4VEC factors.'

  x(1:factor_num,1:point_num) = 0

  do factor_index = 1, factor_num

    if ( factor_index == 1 ) then
      factor_order = 4
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 1, 2, 3, 4 /)
    else if ( factor_index == 2 ) then
      factor_order = 3
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 50, 60, 70 /)
    else if ( factor_index == 3 ) then
      factor_order = 2
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 800, 900 /)
    end if

    call i4vec_direct_product ( factor_index, factor_order, factor_value,  &
      factor_num, point_num, x )

    deallocate ( factor_value )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     J     X(1)  X(2)  X(3)'
  write ( *, '(a)' ) ' '

  do j = 1, point_num
    write ( *, '(2x,i4,4x,i4,2x,i4,2x,i4)' ) j, x(1:factor_num,j)
  end do

  return
end
subroutine test53 ( )

!*****************************************************************************80
!
!! TEST53 tests I4VEC_DIRECT_PRODUCT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: factor_num = 3
  integer ( kind = 4 ), parameter :: point_num = 24

  integer ( kind = 4 ) factor_index
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ), allocatable, dimension ( : ) :: factor_value
  integer ( kind = 4 ) j
  integer ( kind = 4 ) w(point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST53'
  write ( *, '(a)' ) '  I4VEC_DIRECT_PRODUCT2 forms the entries of a'
  write ( *, '(a)' ) '  direct product of a given number of I4VEC factors.'

  w(1:point_num) = 1

  do factor_index = 1, factor_num

    if ( factor_index == 1 ) then
      factor_order = 4
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 2, 3, 5, 7 /)
    else if ( factor_index == 2 ) then
      factor_order = 3
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 11, 13, 17 /)
    else if ( factor_index == 3 ) then
      factor_order = 2
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 19, 21 /)
    end if

    call i4vec_direct_product2 ( factor_index, factor_order, factor_value,  &
      factor_num, point_num, w )

    deallocate ( factor_value )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     J        W(J)'
  write ( *, '(a)' ) ' '

  do j = 1, point_num
    write ( *, '(2x,i4,4x,i8)' ) j, w(j)
  end do

  return
end
subroutine test54 ( )

!*****************************************************************************80
!
!! TEST54 tests I4VEC_FRAC;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) afrac
  integer ( kind = 4 ) :: b = 1
  integer ( kind = 4 ) :: c = 2 * n
  integer ( kind = 4 ) k
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST54'
  write ( *, '(a)' ) '  I4VEC_FRAC: K-th smallest integer vector entry.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  The array to search:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Fractile    Value'
  write ( *, '(a)' ) ' '

  do k = 1, n, n/2

    call i4vec_frac ( n, a, k, afrac )

    write ( *, '(2x,2i8)' ) k, afrac

  end do

  return
end
subroutine test55 ( )

!*****************************************************************************80
!
!! TEST55 tests I4VEC_HEAP_A and I4VEC_HEAP_D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = n
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST55'
  write ( *, '(a)' ) '  For an integer vector,'
  write ( *, '(a)' ) '  I4VEC_HEAP_A puts into ascending heap form.'
  write ( *, '(a)' ) '  I4VEC_HEAP_D puts into descending heap form.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted array:' )

  call i4vec_heap_a ( n, a )

  call i4vec_print ( n, a, '  Ascending heap form:' )

  call i4vec_heap_d ( n, a )

  call i4vec_print ( n, a, '  Descending heap form:' )

  return
end
subroutine test56 ( )

!*****************************************************************************80
!
!! TEST56 tests I4VEC_HEAP_D_EXTRACT, I4VEC_HEAP_D_INSERT and I4VEC_HEAP_D_MAX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  integer ( kind = 4 ) a(n_max)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) val

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST56'
  write ( *, '(a)' ) '  For a descending heap of integers,'
  write ( *, '(a)' ) '  I4VEC_HEAP_D_INSERT inserts a value into the heap.'
  write ( *, '(a)' ) '  I4VEC_HEAP_D_EXTRACT extracts the maximum value;'
  write ( *, '(a)' ) '  I4VEC_HEAP_D_MAX reports the maximum value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  These 3 operations are enough to model a priority queue.'

  n = 0

  seed = 123456789

  do i = 1, n_max

    b = 0
    c = 10

    val = i4_uniform_ab ( b, c, seed )

    call i4vec_heap_d_insert ( n, a, val )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Inserting value          ', val

    call i4vec_heap_d_max ( n, a, val )

    write ( *, '(a,i8)' ) '  Current maximum value is ', val

  end do

  call i4vec_print ( n, a, '  Current heap as a vector:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now extract the maximum several times.'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call i4vec_heap_d_extract ( n, a, val )
    write ( *, '(a,i8)' ) '  Extracting maximum element = ', val
  end do

  call i4vec_print ( n, a, '  Current heap as a vector:' )

  return
end
subroutine test57 ( )

!*****************************************************************************80
!
!! TEST57 tests I4VEC_HISTOGRAM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 1000

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: histo_gram
  integer ( kind = 4 ) histo_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST57'
  write ( *, '(a)' ) '  I4VEC_HISTOGRAM histograms an integer vector.'

  call i4vec_uniform_ab ( n, 0, 25, seed, a )

  histo_num = 20
  allocate ( histo_gram(0:histo_num) )

  call i4vec_histogram ( n, a, histo_num, histo_gram )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Histogram of data from 0 to ', histo_num
  write ( *, '(a)' ) ' '

  do i = 0, histo_num
    if ( 0 < histo_gram(i) ) then
      write ( *, '(2x,i8,2x,i8)' ) i, histo_gram(i)
    end if
  end do

  return
end
subroutine test58 ( )

!*****************************************************************************80
!
!! TEST58 tests I4VEC_INDEX_INSERT_UNIQUE and I4VEC_INDEX_SEARCH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ) xval

  n = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST58'
  write ( *, '(a)' ) '  I4VEC_INDEX_INSERT_UNIQUE inserts unique values into an'
  write ( *, '(a)' ) '  index sorted array.'
  write ( *, '(a)' ) '  I4VEC_INDEX_SEARCH searches for an entry with '
  write ( *, '(a)' ) '  a given value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate some random values:'

  b = 0
  c = n_max
  seed = 123456789

  do i = 1, n_max
    xval = i4_uniform_ab ( b, c, seed )
    call i4vec_index_insert_unique ( n, x, indx, xval )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I   INDX(I)  X(I)  X(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,i3,9x,i3)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Results of search for given XVAL:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XVAL  Less Equal More'
  write ( *, '(a)' ) ' '

  do xval = 0, 20
    call i4vec_index_search ( n, x, indx, xval, less, equal, more )
    write ( *, '(2x,i3,3x,i3,3x,i3,3x,i3)' ) xval, less, equal, more
  end do

  return
end
subroutine test59 ( )

!*****************************************************************************80
!
!! TEST59 tests I4VEC_INDEX_INSERT, I4VEC_INDEX_DELETE_DUPES, I4VEC_INDEX_DELETE_ALL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ) xval

  n = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST59'
  write ( *, '(a)' ) '  I4VEC_INDEX_INSERT inserts values into an'
  write ( *, '(a)' ) '  index sorted array of integers.'
  write ( *, '(a)' ) '  I4VEC_INDEX_DELETE_ALL deletes all copies of a'
  write ( *, '(a)' ) '  particular value.'
  write ( *, '(a)' ) '  I4VEC_INDEX_DELETE_ONE deletes one copies of a'
  write ( *, '(a)' ) '  particular value.'
  write ( *, '(a)' ) '  I4VEC_INDEX_DELETE_DUPES deletes duplicates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate some random values:'
  write ( *, '(a)' ) ' '

  xval = 8
  call i4vec_index_insert ( n, x, indx, xval )

  xval = 7
  call i4vec_index_insert ( n, x, indx, xval )

  b = 0
  c = 20
  seed = 123456789

  do i = 1, 20
    xval = i4_uniform_ab ( b, c, seed )
    write ( *, '(4x,i3)' ) xval
    call i4vec_index_insert ( n, x, indx, xval )
  end do

  xval = 7
  call i4vec_index_insert ( n, x, indx, xval )

  xval = 8
  call i4vec_index_insert ( n, x, indx, xval )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,i3,9x,i3)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call I4VEC_INDEX_DELETE_ONE to delete a value of 8:'

  xval = 8
  call i4vec_index_delete_one ( n, x, indx, xval, n, x, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call I4VEC_INDEX_DELETE_ALL to delete values of 7:'

  xval = 7
  call i4vec_index_delete_all ( n, x, indx, xval )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,i3,9x,i3)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call I4VEC_INDEX_DELETE_DUPES to delete duplicates:'

  call i4vec_index_delete_dupes ( n, x, indx, n, x, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of unique entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,i3)' ) i, indx(i), x(i)
  end do

  return
end
subroutine test60 ( )

!*****************************************************************************80
!
!! TEST60 tests I4VEC_INDEX_INSERT_UNIQUE and I4VEC_INDEX_ORDER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ) xval

  n = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST60'
  write ( *, '(a)' ) '  I4VEC_INDEX_INSERT_UNIQUE inserts unique values into'
  write ( *, '(a)' ) '  an index sorted array.'
  write ( *, '(a)' ) '  I4VEC_INDEX_ORDER sorts an index sorted array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate some random values:'
  write ( *, '(a)' ) ' '

  b = 0
  c = 20
  seed = 123456789

  do i = 1, 20
    xval = i4_uniform_ab ( b, c, seed )
    write ( *, '(4x,i3)' ) xval
    call i4vec_index_insert_unique ( n, x, indx, xval )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of unique entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,i3,9x,i3)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call I4VEC_INDEX_ORDER to carry out the sorting:'

  call i4vec_index_order ( n, x, indx )

  call i4vec_print ( n, x, '  X:' )

  return
end
subroutine test602 ( )

!*****************************************************************************80
!
!! TEST602 tests I4VEC_INDEXED_HEAP_D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 20
  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) :: a(m) = (/ &
    101, 102, 103, 104, 105, 106, 107, 108, 109, 110, &
    111, 112, 113, 114, 115, 116, 117, 118, 119, 120 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: indx(n) = (/ &
    1, 11, 17, 5, 7, 13, 15, 3, 19, 9 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST602'
  write ( *, '(a)' ) '  I4VEC_INDEXED_HEAP_D creates a descending heap'
  write ( *, '(a)' ) '  from an indexed vector.'
!
!  Print before.
!
  call i4vec_print ( m, a, '  The data vector:' )
  call i4vec_print ( n, indx, '  The index vector:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX):'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
  end do
!
!  Heap the data.
!
  call i4vec_indexed_heap_d ( n, a, indx )
!
!  Print afterwards.  Only INDX should change.
!
  call i4vec_print ( m, a, '  The data vector (should NOT change):' )
  call i4vec_print ( n, indx, '  The index vector (may change):' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX) is now a descending heap:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
  end do

  return
end
subroutine test605 ( )

!*****************************************************************************80
!
!! TEST605 tests I4VEC_INDEXED_HEAP_D_EXTRACT and related routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 20
  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) a(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) indx_extract
  integer ( kind = 4 ) indx_insert
  integer ( kind = 4 ) indx_max
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST605'
  write ( *, '(a)' ) '  For an indexed I4VEC,'
  write ( *, '(a)' ) '  I4VEC_INDEXED_HEAP_D_INSERT inserts a value into the heap.'
  write ( *, '(a)' ) '  I4VEC_INDEXED_HEAP_D_EXTRACT extracts the maximum value;'
  write ( *, '(a)' ) '  I4VEC_INDEXED_HEAP_D_MAX reports the maximum value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  These 3 operations are enough to model a priority queue.'
!
!  Set the data array.  To keep things easy, we will use the indicator vector.
!
  call i4vec_indicator ( m, a )
!
!  The index array will initially be a random subset of the numbers 1 to M,
!  in random order.
!
  n = 5
  indx(1:11) = (/ 9, 2, 8, 14, 5, 7, 15, 1, 19, 20, 3 /)

  call i4vec_print ( m, a, '  The data vector:' )
  call i4vec_print ( n, indx, '  The index vector:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX):'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
  end do
!
!  Create a descending heap from the indexed array.
!
  call i4vec_indexed_heap_d ( n, a, indx )

  call i4vec_print ( n, indx, '  The index vector after heaping:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX) after heaping:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
  end do
!
!  Insert five entries, and monitor the maximum.
!
  do i = 1, 5

    indx_insert = indx(n+1)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Inserting value          ', a(indx_insert)

    call i4vec_indexed_heap_d_insert ( n, a, indx, indx_insert )

    call i4vec_indexed_heap_d_max ( n, a, indx, indx_max )

    write ( *, '(a,i8)' ) '  Current maximum is ', a(indx_max)

  end do
  call i4vec_print ( m, a, '  The data vector after insertions:' )
  call i4vec_print ( n, indx, '  The index vector after insertions:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX) after insertions:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
  end do
!
!  Extract the first 5 largest elements.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now extract the maximum several times.'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call i4vec_indexed_heap_d_extract ( n, a, indx, indx_extract )
    write ( *, '(a,i8,a,i8)' ) '  Extracting maximum element A(', &
      indx_extract,') = ', a(indx_extract)
  end do

  call i4vec_print ( m, a, '  The data vector after extractions:' )
  call i4vec_print ( n, indx, '  The index vector after extractions:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX) after extractions:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
  end do

  return
end
subroutine test61 ( )

!*****************************************************************************80
!
!! TEST61 tests I4VEC_INDICATOR;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST61'
  write ( *, '(a)' ) '  I4VEC_INDICATOR sets A = (1,2,3...,N)'

  call i4vec_indicator ( n, a )

  call i4vec_print ( n, a, '  The "indicator" vector:' )

  return
end
subroutine test62 ( )

!*****************************************************************************80
!
!! TEST62 tests I4VEC_MAX and I4VEC_MIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_max
  integer ( kind = 4 ) a_min
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST62'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_MAX:           maximum entry;'
  write ( *, '(a)' ) '  I4VEC_MIN:           minimum entry;'

  b = 1
  c = 30
  seed = 123456789

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  call i4vec_max ( n, a, a_max )
  call i4vec_min ( n, a, a_min )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Maximum:                  ', a_max
  write ( *, '(a,i8)' ) '  Minimum:                  ', a_min

  return
end
subroutine test63 ( )

!*****************************************************************************80
!
!! TEST63 tests I4VEC_MEAN and I4VEC_MEDIAN;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  integer ( kind = 4 ) median
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST63'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_MEAN:          mean value;'
  write ( *, '(a)' ) '  I4VEC_MEDIAN:        median value;'

  seed = 123456789

  b = -n
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  call i4vec_mean ( n, a, mean )
  call i4vec_median ( n, a, median )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Mean:                     ', mean
  write ( *, '(a,i8)' ) '  Median:                   ', median

  return
end
subroutine test64 ( )

!*****************************************************************************80
!
!! TEST64 tests I4VEC_MERGE_A;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 10
  integer ( kind = 4 ), parameter :: n2 = 10

  integer ( kind = 4 ) a1(n1)
  integer ( kind = 4 ) a2(n2)
  integer ( kind = 4 ) a3(n1+n2)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) index
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) search_val
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST64'
  write ( *, '(a)' ) '  For ascending order:'
  write ( *, '(a)' ) '  I4VEC_MERGE_A merges two sorted integer arrays;'

  seed = 123456789

  b = 0
  c = n1

  call i4vec_uniform_ab ( n1, b, c, seed, a1 )

  search_val = a1(1)

  call i4vec_sort_heap_a ( n1, a1 )

  b = 0
  c = n2

  call i4vec_uniform_ab ( n2, b, c, seed, a2 )

  call i4vec_sort_heap_a ( n2, a2 )

  call i4vec_print ( n1, a1, '  Input vector A1:' )

  call i4vec_print ( n2, a2, '  Input vector A2:' )

  call i4vec_merge_a ( n1, a1, n2, a2, n3, a3 )

  call i4vec_print ( n3, a3, '  Merged vector A3:' )

  return
end
subroutine test65 ( )

!*****************************************************************************80
!
!! TEST65 tests I4VEC_NONZERO_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 15

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i4vec_nonzero_count
  integer ( kind = 4 ) nonzero
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST65'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_NONZERO_COUNT: number of nonzeroes;'

  seed = 123456789

  b = -3
  c = 4

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  nonzero = i4vec_nonzero_count ( n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nonzeroes :     ', nonzero

  return
end
subroutine test66 ( )

!*****************************************************************************80
!
!! TEST66 tests I4VEC_NONZERO_FIRST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_save(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST66'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_NONZERO_FIRST left shifts the nonzero entries'
  write ( *, '(a)' ) '  of an I4VEC so they appear first.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  ----------Before--------------    ----------After---------------'
  write ( *, '(a)' ) ' '
  seed = 123456789

  ilo = -1
  ihi = +2

  do test = 1, test_num

    call i4vec_uniform_ab ( n, ilo, ihi, seed, a )
    a_save(1:n) = a(1:n)
    call i4vec_nonzero_first ( n, a, nz, indx )
    write ( *, '(2x,10i3,4x,10i3)' ) a_save(1:n), a(1:n)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The value NZ counts the nonzeros, and'
  write ( *, '(a)' ) '  the vector INDX indicates the original positions:'
  write ( *, '(a)' ) ' '

  call i4vec_uniform_ab ( n, ilo, ihi, seed, a )
  a_save(1:n) = a(1:n)
  call i4vec_nonzero_first ( n, a, nz, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,10i3)' ) a_save(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i2)' ) '  Number of nonzeros NZ = ', nz
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Shifted vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,10i3)' ) a(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,10i3)' ) indx(1:n)

  return
end
subroutine test67 ( )

!*****************************************************************************80
!
!! TEST67 tests I4VEC_ORDER_TYPE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) j
  integer ( kind = 4 ) order
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x(n)
!
!  Each ROW of the definition is a COLUMN of the matrix.
!
  integer ( kind = 4 ), dimension ( n, test_num ) :: x_test = reshape ( (/ &
    1, 3, 2, 4, &
    2, 2, 2, 2, &
    1, 2, 2, 4, &
    1, 2, 3, 4, &
    4, 4, 3, 1, &
    9, 7, 3, 0 /), (/ n, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST67'
  write ( *, '(a)' ) '  I4VEC_ORDER_TYPE classifies an integer vector as'
  write ( *, '(a)' ) '  -1: no order'
  write ( *, '(a)' ) '   0: all equal;'
  write ( *, '(a)' ) '   1: ascending;'
  write ( *, '(a)' ) '   2: strictly ascending;'
  write ( *, '(a)' ) '   3: descending;'
  write ( *, '(a)' ) '   4: strictly descending.'

  do test = 1, test_num

    x(1:n) = x_test(1:n,test)

    call i4vec_order_type ( n, x, order )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The following vector has order type ', order
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(2x,i8,i8)' ) j, x(j)
    end do

  end do

  return
end
subroutine test68 ( )

!*****************************************************************************80
!
!! TEST68 tests I4VEC_PAIRWISE_PRIME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 6

  logical i4vec_pairwise_prime
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x(n)
!
!  Each ROW of the definition is a COLUMN of the matrix.
!
  integer ( kind = 4 ), dimension ( n, test_num ) :: x_test = reshape ( (/ &
     1,  3,  2,  4, &
     2,  2,  2,  2, &
     5,  7, 12, 29, &
     1, 13,  1, 11, &
     1,  4,  9, 16, &
     6, 35, 13, 77 /), (/ n, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST68'
  write ( *, '(a)' ) '  I4VEC_PAIRWISE_PRIME determines if a vector of'
  write ( *, '(a)' ) '  integers is pairwise prime.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '              Pairwise'
  write ( *, '(a)' ) '  Row Vector     Prime?'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    x(1:n) = x_test(1:n,test)

    write ( *, '(2x,4i3,3x,l1)' ) x(1:n), i4vec_pairwise_prime ( n, x )

  end do

  return
end
subroutine test69 ( )

!*****************************************************************************80
!
!! TEST69 tests I4VEC_PART.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) nval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST69'
  write ( *, '(a)' ) '  I4VEC_PART partitions an integer.'

  nval = 17
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NVAL = ', nval

  call i4vec_part ( n, nval, a )

  call i4vec_print ( n, a, '  Partitioned:' )

  nval = -49
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NVAL = ', nval

  call i4vec_part ( n, nval, a )

  call i4vec_print ( n, a, '  Partitioned:' )

  return
end
subroutine test70 ( )

!*****************************************************************************80
!
!! TEST70 tests I4VEC_PART_QUICK_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = n
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST70'
  write ( *, '(a)' ) '  I4VEC_PART_QUICK_A reorders an integer vector'
  write ( *, '(a)' ) '  as part of a quick sort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Before rearrangement:' )

  call i4vec_part_quick_a ( n, a, l, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rearranged array'
  write ( *, '(a,i8)' ) '  Left index =  ', l
  write ( *, '(a,i8)' ) '  Key index =   ', l+1
  write ( *, '(a,i8)' ) '  Right index = ', r

  call i4vec_print ( l,     a(1:l),   '  Left half:' )
  call i4vec_print ( 1,     a(l+1),   '  Key:' )
  call i4vec_print ( n-l-1, a(l+2:n), '  Right half:' )

  return
end
subroutine test71 ( )

!*****************************************************************************80
!
!! TEST71 tests I4VEC_PERMUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) :: c = n
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST71'
  write ( *, '(a)' ) '  I4VEC_PERMUTE reorders an integer vector'
  write ( *, '(a)' ) '  according to a given permutation.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  A, before rearrangement:' )

  call perm_uniform ( n, base, seed, p )

  call i4vec_print ( n, p, '  Permutation vector P:' )

  call i4vec_permute ( n, p, a )

  call i4vec_print ( n, a, '  A, after rearrangement:' )

  return
end
subroutine test72 ( )

!*****************************************************************************80
!
!! TEST72 tests I4VEC_REVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = 3 * n
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST72'
  write ( *, '(a)' ) '  I4VEC_REVERSE reverses a list of integers.'

  seed = 123456789

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Original vector:' )

  call i4vec_reverse ( n, a )

  call i4vec_print ( n, a, '  Reversed:' )

  a(1:n) = a(n:1:-1)

  call i4vec_print ( n, a, '  Re-reversed array using a(1:n) = a(n:1:-1):' )

  return
end
subroutine test73 ( )

!*****************************************************************************80
!
!! TEST73 tests I4VEC_RUN_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) run_count
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST73'
  write ( *, '(a)' ) '  I4VEC_RUN_COUNT counts runs in an I4VEC'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Run Count        Sequence'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    call i4vec_uniform_ab ( n, 0, 1, seed, a )

    call i4vec_run_count ( n, a, run_count )

    write ( *, '(2x,i8,8x,20i2)' ) run_count, a(1:n)

  end do

  return
end
subroutine test74 ( )

!*****************************************************************************80
!
!! TEST74 tests I4VEC_SEARCH_BINARY_A;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) index
  integer ( kind = 4 ) search_val
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST74'
  write ( *, '(a)' ) '  For ascending order:'
  write ( *, '(a)' ) '  I4VEC_SEARCH_BINARY_A searchs an array for a value;'

  seed = 123456789

  b = 0
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  search_val = a(1)

  call i4vec_sort_heap_a ( n, a )

  call i4vec_print ( n, a, '  Input vector A:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Search the array A for the value ', search_val

  call i4vec_search_binary_a ( n, a, search_val, index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEARCH RESULT:'
  if ( 0 < index ) then
    write ( *, '(a,i8)' ) '    The value occurs in index ', index
  else
    write ( *, '(a)' ) '    The value does not occur in the array.'
  end if

  return
end
subroutine test75 ( )

!*****************************************************************************80
!
!! TEST75 tests I4VEC_SORT_BUBBLE_A and I4VEC_SORTED_UNIQUE_HIST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: unique_max = 20
  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) acount(unique_max)
  integer ( kind = 4 ) auniq(unique_max)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = 3 * n
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST75'
  write ( *, '(a)' ) '  For a vector of integers,'
  write ( *, '(a)' ) '  I4VEC_SORT_BUBBLE_A ascending sorts,'
  write ( *, '(a)' ) '  I4VEC_SORTED_UNIQUE_HIST makes a histogram '
  write ( *, '(a)' ) '  of unique entries.'

  seed = 123456789

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted:' )

  call i4vec_sort_bubble_a ( n, a )

  call i4vec_print ( n, a, '  Ascending sorted:' )

  call i4vec_sorted_unique_hist ( n, a, unique_max, unique_num, auniq, acount );

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  I4VEC_UNIQ3 counts ',  unique_num, ' unique entries.'

  call i4vec2_print ( unique_num, auniq, acount, '  Value and Multiplicity' )

  return
end
subroutine test76 ( )

!*****************************************************************************80
!
!! TEST76 tests I4VEC_SORT_HEAP_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = 3 * n
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST76'
  write ( *, '(a)' ) '  For a vector of integers,'
  write ( *, '(a)' ) '  I4VEC_SORT_HEAP_A ascending sorts,'

  seed = 123456789

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted:' )

  call i4vec_sort_heap_a ( n, a )

  call i4vec_print ( n, a, '  Ascending sorted:' )

  return
end
subroutine test77 ( )

!*****************************************************************************80
!
!! TEST77 tests I4VEC_SORT_HEAP_D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = 3 * n
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST77'
  write ( *, '(a)' ) '  For a vector of integers,'
  write ( *, '(a)' ) '  I4VEC_SORT_HEAP_D descending sorts.'

  seed = 123456789

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted:' )

  call i4vec_sort_heap_d ( n, a )

  call i4vec_print ( n, a, '  Descending sorted:' )

  return
end
subroutine test78 ( )

!*****************************************************************************80
!
!! TEST78 tests I4VEC_SORT_HEAP_INDEX_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST78'
  write ( *, '(a)' ) '  I4VEC_SORT_HEAP_INDEX_A creates an ascending'
  write ( *, '(a)' ) '  sort index for an integer array.'

  seed = 123456789

  b = 0
  c = 3 * n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted array:' )

  call i4vec_sort_heap_index_a ( n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After indexed ascending sort:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), A(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(3i8)' ) i, indx(i), a(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now use the index array to carry out the'
  write ( *, '(a)' ) '  permutation implicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), A(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(3i8)' ) i, indx(i), a(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call I4VEC_PERMUTE to carry out the permutation'
  write ( *, '(a)' ) '  explicitly.'

  call i4vec_permute ( n, indx, a )

  call i4vec_print ( n, a, '  I, A(I)' )

  return
end
subroutine test79 ( )

!*****************************************************************************80
!
!! TEST79 tests I4VEC_SORT_HEAP_INDEX_D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST79'
  write ( *, '(a)' ) '  I4VEC_SORT_HEAP_INDEX_D creates a descending'
  write ( *, '(a)' ) '  sort index for an I4VEC.'

  seed = 123456789

  b = 0
  c = 3 * n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted array:' )

  call i4vec_sort_heap_index_d ( n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After indexed descending sort:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), A(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,3i8)' ) i, indx(i), a(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now use the index array to carry out the'
  write ( *, '(a)' ) '  permutation implicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), A(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,3i8)' ) i, indx(i), a(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call I4VEC_PERMUTE to carry out the permutation'
  write ( *, '(a)' ) '  explicitly.'

  call i4vec_permute ( n, indx, a )

  call i4vec_print ( n, a, '  I, A(I)' )

  return
end
subroutine test80 ( )

!*****************************************************************************80
!
!! TEST80 tests I4VEC_SORT_INSERT_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST80'
  write ( *, '(a)' ) '  I4VEC_SORT_INSERT_A sorts an integer array.'

  seed = 123456789

  b = 0
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted array:' )

  call i4vec_sort_insert_a ( n, a )

  call i4vec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test81 ( )

!*****************************************************************************80
!
!! TEST81 tests I4VEC_SORT_QUICK_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST81'
  write ( *, '(a)' ) '  I4VEC_SORT_QUICK_A sorts an integer vector'
  write ( *, '(a)' ) '  using quick sort.'

  seed = 123456789

  b = 0
  c = 3 * n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted array:' )

  call i4vec_sort_quick_a ( n, a )

  call i4vec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test82 ( )

!*****************************************************************************80
!
!! TEST82 tests I4VEC_SORT_SHELL_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST82'
  write ( *, '(a)' ) '  I4VEC_SORT_SHELL_A sorts an integer vector'
  write ( *, '(a)' ) '  using Shell''s sort.'

  seed = 123456789

  b = 0
  c = 3 * n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted array:' )

  call i4vec_sort_shell_a ( n, a )

  call i4vec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test83 ( )

!*****************************************************************************80
!
!! TEST83 tests I4VEC_SORTED_UNDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: x_num = 9

  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
  integer ( kind = 4 ) x_unique_num
  integer ( kind = 4 ), dimension ( x_num ) :: x_val = (/ &
    11, 11, 11, 22, 22, 33, 33, 55, 55 /)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: xu_val
  integer ( kind = 4 ) xdnu(x_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST83'
  write ( *, '(a)' ) '  I4VEC_SORTED_UNDEX produces index vectors which create a sorted'
  write ( *, '(a)' ) '  list of the unique elements of a sorted I4VEC,'
  write ( *, '(a)' ) '  and a map from the original vector to the (implicit)'
  write ( *, '(a)' ) '  vector of sorted unique elements.'

  call i4vec_print ( x_num, x_val, '  The vector X:' )

  call i4vec_sorted_unique_count ( x_num, x_val, x_unique_num )

  allocate ( undx(1:x_unique_num) )
  allocate ( xu_val(1:x_unique_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique entries in X is ', x_unique_num

  call i4vec_sorted_undex ( x_num, x_val, x_unique_num, undx, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to list the unique elements of X'
  write ( *, '(a)' ) '  in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX   X(UNDX)'
  write ( *, '(a)' ) ' '

  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,i8)' ) i, undx(i), x_val(undx(i))
  end do

  xu_val(1:x_unique_num) = x_val(undx(1:x_unique_num))

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to created XU, a copy of X'
  write ( *, '(a)' ) '  containing only the unique elements, in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX XU(I)'
  write ( *, '(a)' ) ' '
  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,i4)' ) i, undx(i), xu_val(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XDNU can be used to match each element of X with one of the'
  write ( *, '(a)' ) '  unique elements'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  XDNU  X(I)   XU(XDNU(I))'
  write ( *, '(a)' ) ' '

  do i = 1, x_num
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,i12)' ) i, xdnu(i), x_val(i), xu_val(xdnu(i))
  end do

  deallocate ( undx )
  deallocate ( xu_val )

  return
end
subroutine test84 ( )

!*****************************************************************************80
!
!! TEST84 tests I4VEC_SORTED_UNIQUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) :: b = 0
  integer ( kind = 4 ) :: c = n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) unique_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST84'
  write ( *, '(a)' ) '  I4VEC_SORTED_UNIQUE finds unique entries in a sorted array.'

  seed = 123456789

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_sort_heap_a ( n, a )

  call i4vec_print ( n, a, '  Input vector:' )

  call i4vec_sorted_unique ( n, a, unique_num )

  call i4vec_print ( unique_num, a, '  Unique entries:' )

  return
end
subroutine test85 ( )

!*****************************************************************************80
!
!! TEST85 tests I4VEC_TRANSPOSE_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  integer ( kind = 4 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST85'
  write ( *, '(a)' ) '  I4VEC_TRANSPOSE_PRINT prints an integer vector'
  write ( *, '(a)' ) '  with 5 entries to a row, and an optional title.'

  call i4vec_indicator ( n, a )

  call i4vec_print ( n, a, '  Output from I4VEC_PRINT:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call I4VEC_TRANSPOSE_PRINT with a short title:'

  call i4vec_transpose_print ( n, a, '  My array:  ' )

  return
end
subroutine test86 ( )

!*****************************************************************************80
!
!! TEST86 tests I4VEC_UNDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: x_num = 9

  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
  integer ( kind = 4 ) x_unique_num
  integer ( kind = 4 ), dimension ( x_num ) :: x_val = (/ &
    33, 55, 11, 11, 55, 33, 22, 22, 11 /)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: xu_val
  integer ( kind = 4 ) xdnu(x_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST86'
  write ( *, '(a)' ) '  I4VEC_UNDEX produces index vectors which create a sorted'
  write ( *, '(a)' ) '  list of the unique elements of an (unsorted) I4VEC,'
  write ( *, '(a)' ) '  and a map from the original vector to the (implicit)'
  write ( *, '(a)' ) '  vector of sorted unique elements.'

  call i4vec_print ( x_num, x_val, '  The vector X:' )

  call i4vec_unique_count ( x_num, x_val, x_unique_num )

  allocate ( undx(1:x_unique_num) )
  allocate ( xu_val(1:x_unique_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique entries in X is ', x_unique_num

  call i4vec_undex ( x_num, x_val, x_unique_num, undx, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to list the unique elements of X'
  write ( *, '(a)' ) '  in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX   X(UNDX)'
  write ( *, '(a)' ) ' '

  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,i8)' ) i, undx(i), x_val(undx(i))
  end do

  xu_val(1:x_unique_num) = x_val(undx(1:x_unique_num))

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to created XU, a copy of X'
  write ( *, '(a)' ) '  containing only the unique elements, in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX XU(I)'
  write ( *, '(a)' ) ' '
  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,i4)' ) i, undx(i), xu_val(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XDNU can be used to match each element of X with one of the'
  write ( *, '(a)' ) '  unique elements'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  XDNU  X(I)   XU(XDNU(I))'
  write ( *, '(a)' ) ' '

  do i = 1, x_num
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,i12)' ) i, xdnu(i), x_val(i), xu_val(xdnu(i))
  end do

  deallocate ( undx )
  deallocate ( xu_val )

  return
end
subroutine test87 ( )

!*****************************************************************************80
!
!! TEST87 tests I4VEC_UNIQUE_INDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) unique_index(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST87'
  write ( *, '(a)' ) '  I4VEC_UNIQUE_INDEX, for each entry in an I4VEC'
  write ( *, '(a)' ) '  indexes the unique elements.'

  b = 1
  c = 5

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_unique_index ( n, a, unique_index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I      A(I)    UNIQUE'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, a(i), unique_index(i)
  end do

  return
end
subroutine test88 ( )

!*****************************************************************************80
!
!! TEST88 tests I4VEC_VALUE_INDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_index = 3
  integer ( kind = 4 ), parameter :: n = 25

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) n_index
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) value_index(max_index)

  seed = 123456789

  value = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST88'
  write ( *, '(a)' ) '  I4VEC_VALUE_INDEX indexes entries equal to'
  write ( *, '(a)' ) '  a given value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The desired value is ', value
  write ( *, '(a,i8)' ) '  Maximum number of indices to find is ', max_index

  b = 1
  c = 5

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector A:' )

  call i4vec_value_index ( n, a, value, max_index, n_index, value_index )

  call i4vec_print ( n_index, value_index, &
    '  Indices of entries equal to given value: ' )

  return
end
subroutine test89 ( )

!*****************************************************************************80
!
!! TEST89 tests I4VEC_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) seed
  real ( kind = 8 ) variance

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST89'
  write ( *, '(a)' ) '  For an integer vector:'
  write ( *, '(a)' ) '  I4VEC_VARIANCE:      variance.'

  seed = 123456789

  b = -n
  c = n

  call i4vec_uniform_ab ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Input vector:' )

  call i4vec_variance ( n, a, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Variance:                 ', variance

  return
end
subroutine test90 ( )

!*****************************************************************************80
!
!! TEST90 tests I4VEC2_SORT_A, I4VEC2_SORT_D, and I4VEC2_SORTED_UNIQUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) jvec(n)
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST90'
  write ( *, '(a)' ) '  For a pair of integer vectors:'
  write ( *, '(a)' ) '  I4VEC2_SORT_A ascending sorts;'
  write ( *, '(a)' ) '  I4VEC2_SORT_D descending sorts;'
  write ( *, '(a)' ) '  I4VEC2_SORTED_UNIQUE counts unique entries.'

  b = 1
  c = 3

  call i4vec_uniform_ab ( n, b, c, seed, ivec )

  call i4vec_uniform_ab ( n, b, c, seed, jvec )

  ivec(3) = ivec(1)
  jvec(3) = jvec(1)

  ivec(5) = ivec(2)
  jvec(5) = jvec(2)

  ivec(9) = ivec(1)
  jvec(9) = jvec(1)

  call i4vec2_print ( n, ivec, jvec, '  The array:' )

  call i4vec2_sort_a ( n, ivec, jvec )

  call i4vec2_print ( n, ivec, jvec, '  After ascending sort:' )

  call i4vec2_sort_d ( n, ivec, jvec )

  call i4vec2_print ( n, ivec, jvec, '  After descending sort:' )

  call i4vec2_sorted_unique ( n, ivec, jvec, unique_num )

  call i4vec2_print ( unique_num, ivec, jvec, '  After UNIQ:' )

  return
end
