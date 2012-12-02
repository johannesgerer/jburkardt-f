program main

!*****************************************************************************80
!
!! MAIN is the main program for R16LIB_PRB.
!
!  Discussion:
!
!    R16LIB_PRB tests routines from the R16LIB library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R16LIB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the R16LIB library.'

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test004 ( )
  call test005 ( )
  call test006 ( )
  call test007 ( )
  call test008 ( )
  call test009 ( )

  call test010 ( )
  call test011 ( )
  call test013 ( )

  call test0183 ( )
  call test1515 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R16LIB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests R16_ABS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 16 ) r16
  real ( kind = 16 ) r16_abs
  real ( kind = 16 ) r16_absolute
  real ( kind = 16 ) r16_uniform
  real ( kind = 16 ) :: r16_hi = 5.0_16
  real ( kind = 16 ) :: r16_lo = -3.0_16
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  R16_ABS returns the absolute value of an R16.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    r16 = r16_uniform ( r16_lo, r16_hi, seed )
    r16_absolute = r16_abs ( r16 )
    write ( *, '(2x,f10.6,2x,f10.6)' ) r16, r16_absolute
  end do

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests R16_ATAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  real ( kind = 16 ) r16_atan
  integer ( kind = 4 ) test
  real ( kind = 16 ) x
  real ( kind = 16 ), dimension ( test_num ) :: xtest = (/ &
     1.0_16,  1.0_16,  0.0_16, -1.0_16, &
    -1.0_16, -1.0_16,  0.0_16,  1.0_16 /)
  real ( kind = 16 ) y
  real ( kind = 16 ), dimension ( test_num ) :: ytest = (/ &
     0.0_16,  1.0_16,  1.0_16,  1.0_16, &
     0.0_16, -1.0_16, -1.0_16, -1.0_16 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  R16_ATAN computes the arc-tangent given Y and X;'
  write ( *, '(a)' ) '  ATAN2 is the system version of this routine.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             Y          ATAN2(Y,X)    R16_ATAN(Y,X)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = xtest(test)
    y = ytest(test)
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      x, y, atan2 ( y, x ), r16_atan ( y, x )
  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests R16_CAS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 16 ) r16_cas
  real ( kind = 16 ) r16_pi
  integer ( kind = 4 ), parameter :: test_num = 12
  integer ( kind = 4 ) test
  real ( kind = 16 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  R16_CAS evaluates the casine of a number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        X          R16_CAS ( X )'
  write ( *, '(a)' ) ' '
  do test = 0, test_num
    x = r16_pi ( ) * real ( test, kind = 16 ) / real ( test_num, kind = 16 )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, r16_cas ( x )
  end do

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests R16_CEILING.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) r16_ceiling
  integer ( kind = 4 ) ival
  real ( kind = 16 ) rval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  R16_CEILING rounds a value up.'
  write ( *, '(a)' ) ' '

  do i = -6, 6
    rval = real ( i, kind = 16 ) / 5.0_16
    ival = r16_ceiling ( rval )
    write ( *, '(2x,g14.6,2x,i8)' ) rval, ival
  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests R16_DIFF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 15

  integer ( kind = 4 ) ndig
  real ( kind = 16 ) r16_diff
  integer ( kind = 4 ) test
  real ( kind = 16 ) x
  real ( kind = 16 ), dimension ( test_num ) :: y_test = (/ &
    0.0625_16, 0.125_16, 0.25_16, 0.50_16,  0.874_16, &
    0.876_16,  0.90_16,  0.95_16, 0.99_16,  1.0_16, &
    1.01_16,   1.05_16,  1.10_16, 3.0_16,  10.0_16 /)
  real ( kind = 16 ) y

  ndig = 3
  x = 1.0_16

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  R16_DIFF computes a difference X-Y to a given'
  write ( *, '(a)' ) '    number of binary places.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  For this test, we use ', ndig, ' binary places.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       Y       X-Y     R16_DIFF(X,Y)'
  write ( *, '(a)' ) ' '
  do test = 1, test_num
    y = y_test(test)
    write ( *, '(4f10.4)' ) x, y, x-y, r16_diff ( x, y, ndig )
  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests R16_DIGIT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxdig = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) digit(-2:maxdig)
  integer ( kind = 4 ) idigit
  real ( kind = 16 ) r16_pi
  real ( kind = 16 ) x

  x = r16_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  R16_DIGIT extracts decimal digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  Here, we get digits of ', x
  write ( *, '(a)' ) ' '

  do idigit = -2, maxdig
    call r16_digit ( x, idigit, digit(idigit) )
  end do

  write ( *, '(2x,25i3)' ) ( i, i = -2, maxdig )
  write ( *, '(2x,25i3)' ) digit(-2:maxdig)

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests R16_EPSILON
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 16 ) r16_epsilon
  real ( kind = 16 ) r
  real ( kind = 16 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  R16_EPSILON produces the R16 machine precision.'
  write ( *, '(a)' ) ' '

  r = r16_epsilon ( )
  write ( *, '(a,g14.6)' ) '  R = R16_EPSILON()   = ', r
  s = ( 1.0_16 + r ) - 1.0_16
  write ( *, '(a,g14.6)' ) '  ( 1 + R ) - 1       = ', s
  s = ( 1.0_16 + ( r / 2.0_16 ) ) - 1.0_16
  write ( *, '(a,g14.6)' ) '  ( 1 + (R/2) ) - 1   = ', s

  write ( *, '(a)' ) ' '
  r = epsilon ( r )
  write ( *, '(a,g14.6)' ) '  R = EPSILON ( R )   = ', r
  s = ( 1.0_16 + r ) - 1.0_16
  write ( *, '(a,g14.6)' ) '  ( 1 + R ) - 1       = ', s
  s = ( 1.0_16 + ( r / 2.0_16 ) ) - 1.0_16
  write ( *, '(a,g14.6)' ) '  ( 1 + (R/2) ) - 1   = ', s

  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests R16_FRACTION.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 16 ) fraction
  real ( kind = 16 ) r16
  real ( kind = 16 ) r16_fraction
  real ( kind = 16 ) r16_uniform
  real ( kind = 16 ) :: r16_hi = 5.0_16
  real ( kind = 16 ) :: r16_lo = -3.0_16
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  R16_FRACTION returns the fraction part of an R16.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    r16 = r16_uniform ( r16_lo, r16_hi, seed )
    fraction = r16_fraction ( r16 )
    write ( *, '(2x,f10.6,2x,f10.6)' ) r16, fraction
  end do

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests R16_HUGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 16 ) r16_huge

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  R16_HUGE returns a "huge" R16;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '    R16_HUGE ( ) =      ', r16_huge ( )
  write ( *, '(a,g24.16)' ) '    HUGE ( 1.0_16 ) = ', huge ( 1.0_16 )

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests R16_LOG_2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 18

  real ( kind = 16 ) r16_log_2
  integer ( kind = 4 ) test
  real ( kind = 16 ) x
  real ( kind = 16 ), dimension(test_num) :: x_test = (/ &
    0.0_16,  1.0_16,  2.0_16,   3.0_16,  9.0_16, &
   10.0_16, 11.0_16, 99.0_16, 101.0_16, -1.0_16, &
   -2.0_16, -3.0_16, -9.0_16,   0.5_16,  0.33_16, &
    0.25_16, 0.20_16, 0.01_16 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  R16_LOG_2 computes the logarithm base 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, R16_LOG_2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = x_test(test)
    write ( *, '( 2g14.6 )' ) x, r16_log_2 ( x )
  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests R16_LOG_B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 10

  real ( kind = 16 ) b
  real ( kind = 16 ), dimension(test_num) :: b_test = (/ &
    2.0_16, 3.0_16, 4.0_16, 5.0_16, 6.0_16, &
    7.0_16, 8.0_16, 16.0_16, 32.0_16, 256.0_16 /)
  real ( kind = 16 ) r16_log_b
  integer ( kind = 4 ) test
  real ( kind = 16 ) x

  x = 16.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  R16_LOG_B computes the logarithm base B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, B, R16_LOG_B'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    b = b_test(test)

    write ( *, '( 2x,3g14.6, i12 )' ) x, b, r16_log_b ( x, b )

  end do

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests R16_MOD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 16 ) r16_mod
  real    ( kind = 16 ) r16_uniform
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ) seed
  real    ( kind = 16 ) x
  real    ( kind = 16 ) :: x_hi =  10.0D+00
  real    ( kind = 16 ) :: x_lo = -10.0D+00
  real    ( kind = 16 ) y
  real    ( kind = 16 ) z1
  real    ( kind = 16 ) z2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  R8_MOD returns the remainder after division.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X         Y     MOD(X,Y)    R16_MOD(X,Y)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    x = r16_uniform ( x_lo, x_hi, seed )
    y = r16_uniform ( x_lo, x_hi, seed )

    z1 =     mod ( x, y )
    z2 = r16_mod ( x, y )

    write ( * , '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) x, y, z1, z2

  end do

  return
end
subroutine test0183 ( )

!*****************************************************************************80
!
!! TEST0183 tests R16_PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 16 ) r16_pi
  real ( kind = 16 ) v1
  real ( kind = 16 ) v2
  real ( kind = 16 ) v3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0183'
  write ( *, '(a)' ) '  R16_PI returns the value of PI.'
  write ( *, '(a)' ) ' '
  v1 = r16_pi ( )
  write ( *, '(a,g42.32)' ) '  R16_PI =    ', v1
  v2 = 4.0_16 * atan ( 1.0_16 )
  write ( *, '(a,g42.32)' ) '  4*atan(1) = ', v2

  return
end
subroutine test1515 ( )

!*****************************************************************************80
!
!! TEST1515 tests R16VEC_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 16 ) r(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1515'
  write ( *, '(a)' ) '  R16VEC_UNIFORM_01 returns a random R8VEC '
  write ( *, '(a)' ) '  with entries in [0,1].'

  seed = 123456789

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Input SEED = ', seed

    call r16vec_uniform_01 ( n, seed, r )

    call r16vec_print_some ( n, r, 1, 10, '  Random vector:' )

  end do

  return
end
