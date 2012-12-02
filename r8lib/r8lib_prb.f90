program main

!*****************************************************************************80
!
!! MAIN is the main program for R8LIB_PRB.
!
!  Discussion:
!
!    R8LIB_PRB tests routines from the R8LIB library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8LIB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the R8LIB library.'

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
  call test012 ( )
  call test013 ( )
  call test014 ( )
  call test015 ( )
  call test016 ( )
  call test017 ( )
  call test018 ( )
  call test019 ( )

  call test020 ( )
  call test021 ( )
  call test022 ( )
  call test023 ( )
  call test0235 ( )
  call test024 ( )
  call test025 ( )
  call test026 ( )
  call test027 ( )
  call test028 ( )
  call test029 ( )
  call test0295 ( )

  call test031 ( )
  call test032 ( )
  call test033 ( )
  call test034 ( )
  call test035 ( )
  call test036 ( )
  call test0365 ( )
  call test037 ( )
  call test038 ( )
  call test0383 ( )
  call test0385 ( )
  call test039 ( )
  call test0393 ( )
  call test0395 ( )
  call test0397 ( )

  call test040 ( )
  call test041 ( )
  call test0415 ( )
  call test042 ( )
  call test043 ( )
  call test044 ( )
  call test0442 ( )
  call test0443 ( )
  call test0445 ( )
  call test045 ( )
  call test046 ( )
  call test047 ( )
  call test048 ( )
  call test049 ( )

  call test050 ( )
  call test051 ( )
  call test052 ( )
  call test053 ( )
  call test054 ( )
  call test055 ( )
  call test0555 ( )
  call test056 ( )
  call test057 ( )
  call test058 ( )
  call test059 ( )

  call test060 ( )
  call test061 ( )
  call test062 ( )
  call test063 ( )
  call test064 ( )
  call test065 ( )
  call test066 ( )
  call test067 ( )
  call test068 ( )
  call test069 ( )

  call test070 ( )
  call test071 ( )
  call test072 ( )
  call test073 ( )
  call test0737 ( )
  call test074 ( )
  call test075 ( )
  call test076 ( )
  call test0764 ( )
  call test0766 ( )
  call test077 ( )
  call test0775 ( )
  call test0776 ( )
  call test078 ( )
  call test079 ( )

  call test080 ( )
  call test081 ( )
  call test082 ( )
  call test083 ( )
  call test084 ( )
  call test085 ( )
  call test086 ( )
  call test087 ( )
  call test088 ( )
  call test089 ( )

  call test090 ( )
  call test091 ( )
  call test092 ( )
  call test093 ( )
  call test094 ( )
  call test095 ( )
  call test098 ( )
  call test099 ( )

  call test100 ( )
  call test101 ( )
  call test102 ( )
  call test103 ( )
  call test104 ( )
  call test105 ( )
  call test106 ( )
  call test1064 ( )
  call test1065 ( )
  call test1066 ( )
  call test1067 ( )
  call test107 ( )
  call test108 ( )
  call test109 ( )

  call test110 ( )
  call test111 ( )
  call test112 ( )
  call test113 ( )
  call test114 ( )
  call test1143 ( )
  call test1145 ( )
  call test1147 ( )
  call test115 ( )
  call test116 ( )
  call test1165 ( )
  call test1166 ( )
  call test117 ( )
  call test118 ( )
  call test119 ( )

  call test120 ( )
  call test121 ( )
  call test1215 ( )
  call test122 ( )
  call test123 ( )
  call test124 ( )
  call test125 ( )
  call test1251 ( )
  call test1252 ( )
  call test1255 ( )
  call test1256 ( )
  call test1258 ( )
  call test126 ( )
  call test127 ( )
  call test128 ( )
  call test129 ( )

  call test130 ( )
  call test152 ( )
  call test131 ( )
  call test132 ( )
  call test133 ( )
  call test134 ( )
  call test135 ( )
  call test136 ( )
  call test137 ( )
  call test138 ( )
  call test139 ( )

  call test140 ( )
  call test141 ( )
  call test142 ( )
  call test143 ( )
  call test144 ( )
  call test145 ( )
  call test146 ( )
  call test1465 ( )
  call test147 ( )
  call test1475 ( )
  call test148 ( )
  call test149 ( )

  call test150 ( )
  call test1504 ( )
  call test1505 ( )
  call test151 ( )
  call test1515 ( )
  call test153 ( )
  call test154 ( )
  call test155 ( )
  call test156 ( )
  call test157 ( )
  call test158 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8LIB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests R8_ABS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8
  real ( kind = 8 ) r8_abs
  real ( kind = 8 ) r8_absolute
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) :: r8_hi = 5.0D+00
  real ( kind = 8 ) :: r8_lo = -3.0D+00
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  R8_ABS returns the absolute value of an R8.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X         R8_ABS(X)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    r8 = r8_uniform_ab ( r8_lo, r8_hi, seed )
    r8_absolute = r8_abs ( r8 )
    write ( *, '(2x,f10.6,2x,f10.6)' ) r8, r8_absolute
  end do

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests R8_ATAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  real ( kind = 8 ) r8_atan
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension ( test_num ) :: xtest = (/ &
     1.0D+00,  1.0D+00,  0.0D+00, -1.0D+00, &
    -1.0D+00, -1.0D+00,  0.0D+00,  1.0D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), dimension ( test_num ) :: ytest = (/ &
     0.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
     0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  R8_ATAN computes the arc-tangent given Y and X;'
  write ( *, '(a)' ) '  ATAN2 is the system version of this routine.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             Y          ATAN2(Y,X)    R8_ATAN(Y,X)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = xtest(test)
    y = ytest(test)
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      x, y, atan2 ( y, x ), r8_atan ( y, x )
  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests R8_CAS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_cas
  real ( kind = 8 ) r8_pi
  integer ( kind = 4 ), parameter :: test_num = 12
  integer ( kind = 4 ) test
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  R8_CAS evaluates the casine of a number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        X           R8_CAS ( X )'
  write ( *, '(a)' ) ' '
  do test = 0, test_num
    x = r8_pi ( ) * real ( test, kind = 8 ) / real ( test_num, kind = 8 )
    write ( *, '(2x,g14.6,2x,g14.6)' ) x, r8_cas ( x )
  end do

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests R8_CEILING.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_ceiling
  real ( kind = 8 ) rval
  real ( kind = 8 ) rval2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  R8_CEILING rounds a value up.'
  write ( *, '(a)' ) ' '

  do i = -6, 6
    rval = real ( i, kind = 8 ) / 5.0D+00
    rval2 = r8_ceiling ( rval )
    write ( *, '(2x,g14.6,2x,g14.6)' ) rval, rval2
  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests R8_DIFF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 15

  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r8_diff
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension ( test_num ) :: y_test = (/ &
    0.0625D+00, 0.125D+00, 0.25D+00, 0.50D+00,  0.874D+00, &
    0.876D+00,  0.90D+00,  0.95D+00, 0.99D+00,  1.0D+00, &
    1.01D+00,   1.05D+00,  1.10D+00, 3.0D+00,  10.0D+00 /)
  real ( kind = 8 ) y

  ndig = 3
  x = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  R8_DIFF computes a difference X-Y to a given'
  write ( *, '(a)' ) '  number of binary places.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  For this test, we use ', ndig, ' binary places.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       Y       X-Y     R8_DIFF(X,Y)'
  write ( *, '(a)' ) ' '
  do test = 1, test_num
    y = y_test(test)
    write ( *, '(4f10.4)' ) x, y, x-y, r8_diff ( x, y, ndig )
  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests R8_DIGIT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2006
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
  real ( kind = 8 ) r8_pi
  real ( kind = 8 ) x

  x = r8_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  R8_DIGIT extracts decimal digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  Here, we get digits of ', x
  write ( *, '(a)' ) ' '

  do idigit = -2, maxdig
    call r8_digit ( x, idigit, digit(idigit) )
  end do

  write ( *, '(2x,25i3)' ) ( i, i = -2, maxdig )
  write ( *, '(2x,25i3)' ) digit(-2:maxdig)

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests R8_EPSILON and R8_EPSILON_COMPUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) r8_epsilon_compute
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  R8_EPSILON returns the R8 machine precision.'
  write ( *, '(a)' ) '  R8_EPSILON_COMPUTE computes the R8 machine precision.'
  write ( *, '(a)' ) ' '

  r1 = r8_epsilon ( )
  write ( *, '(a,g24.16)' ) '  R1 = R8_EPSILON()         = ', r1

  r2 = r8_epsilon_compute ( )
  write ( *, '(a,g24.16)' ) '  R2 = R8_EPSILON_COMPUTE() = ', r2

  s = ( 1.0D+00 + r2 ) - 1.0D+00
  write ( *, '(a,g24.16)' ) '  ( 1 + R2 ) - 1            = ', s

  s = ( 1.0D+00 + ( r2 / 2.0D+00 ) ) - 1.0D+00
  write ( *, '(a,g14.6)' ) '  ( 1 + (R2/2) ) - 1        = ', s

  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests R8_FRACTIONAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fractional
  real ( kind = 8 ) r8
  real ( kind = 8 ) r8_fractional
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) :: r8_hi = 5.0D+00
  real ( kind = 8 ) :: r8_lo = -3.0D+00
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 10

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  R8_FRACTIONAL returns the fractional part of an R8.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    r8 = r8_uniform_ab ( r8_lo, r8_hi, seed )
    fractional = r8_fractional ( r8 )
    write ( *, '(2x,f10.6,2x,f10.6)' ) r8, fractional
  end do

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests R8_HUGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_huge

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  R8_HUGE returns a "huge" R8;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '    R8_HUGE ( ) =      ', r8_huge ( )
  write ( *, '(a,g24.16)' ) '    HUGE ( 1.0D+00 ) = ', huge ( 1.0D+00 )

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests R8_LOG_2.
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

  integer ( kind = 4 ), parameter :: test_num = 18

  real ( kind = 8 ) r8_log_2
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension(test_num) :: x_test = (/ &
    0.0D+00,  1.0D+00,  2.0D+00,   3.0D+00,  9.0D+00, &
   10.0D+00, 11.0D+00, 99.0D+00, 101.0D+00, -1.0D+00, &
   -2.0D+00, -3.0D+00, -9.0D+00,   0.5D+00,  0.33D+00, &
    0.25D+00, 0.20D+00, 0.01D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  R8_LOG_2 computes the logarithm base 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, R8_LOG_2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = x_test(test)
    write ( *, '( 2g14.6 )' ) x, r8_log_2 ( x )
  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests R8_LOG_B.
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

  integer ( kind = 4 ), parameter :: test_num = 10

  real ( kind = 8 ) b
  real ( kind = 8 ), dimension(test_num) :: b_test = (/ &
    2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00, 6.0D+00, &
    7.0D+00, 8.0D+00, 16.0D+00, 32.0D+00, 256.0D+00 /)
  real ( kind = 8 ) r8_log_b
  integer ( kind = 4 ) test
  real ( kind = 8 ) x

  x = 16.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  R8_LOG_B computes the logarithm base B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, B, R8_LOG_B'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    b = b_test(test)

    write ( *, '( 2x,3g14.6, i12 )' ) x, b, r8_log_b ( x, b )

  end do

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests R8_MANT.
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

  integer ( kind = 4 ) is
  integer ( kind = 4 ) l
  real ( kind = 8 ) r
  real ( kind = 8 ) x

  x = -314.159D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  R8_MANT decomposes a value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number to be decomposed:'
  write ( *, '(2x,g14.6)' ) x

  call r8_mant ( x, is, r, l )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,g14.6,a,i8)' ) &
    '  R8_MANT: X = ', is, ' * ', r, ' * 2**', l

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests R8_MOD.
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

  real ( kind = 8 ) r8_mod
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi =  10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  R8_MOD returns the remainder after division.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X         Y     MOD(X,Y)    R8_MOD(X,Y)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    x = r8_uniform_ab ( x_lo, x_hi, seed )
    y = r8_uniform_ab ( x_lo, x_hi, seed )

    z1 =    mod ( x, y )
    z2 = r8_mod ( x, y )

    write ( * , '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) x, y, z1, z2

  end do

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests R8_MODP.
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

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_hi =  10.0D+00
  real ( kind = 8 ) :: x_lo = -10.0D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  R8_MODP returns the remainder after division.'
  write ( *, '(a)' ) '  Unlike the FORTRAN MOD, R8_MODP ( X, Y ) is positive.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X       Y      MOD(X,Y)  R8_MODP(X,Y)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    x = r8_uniform_ab ( x_lo, x_hi, seed )
    y = r8_uniform_ab ( x_lo, x_hi, seed )

    z1 =    mod  ( x, y )
    z2 = r8_modp ( x, y )

    write ( * , '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) x, y, z1, z2

  end do

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests R8_NINT
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

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) r8_nint
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  R8_NINT produces the nearest integer to an R8.'
  write ( *, '(a)' ) ' '

  b = -10.0D+00
  c = +10.0D+00

  do test = 1, test_num
    x = r8_uniform_ab ( b, c, seed )
    write ( *, '(2x,f10.4,2x,i8)' ) x, r8_nint ( x )
  end do

  return;
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests R8_NORMAL_01.
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

  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  R8_NORMAL_01 generates normally distributed'
  write ( *, '(a)' ) '  random values.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    x = r8_normal_01 ( seed )
    write ( *, '(2x,g14.6)' ) x

  end do

  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests R8_PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) four
  real ( kind = 8 ) one
  real ( kind = 8 ) r8_pi
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2

  four = real ( 4, kind = 8 )
  one = real ( 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  R8_PI returns the value of PI.'
  write ( *, '(a)' ) ' '
  v1 = r8_pi ( )
  write ( *, '(a,g24.16)' ) '  R8_PI =     ', v1
  v2 = four * atan ( one )
  write ( *, '(a,g24.16)' ) '  4*atan(1) = ', v2

  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests R8_POWER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_power
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p
  real ( kind = 8 ) r
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018'
  write ( *, '(a)' ) '  R8_POWER computes R**P.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      R          P       R**P'
  write ( *, '(a)' ) ' '

  do i = -5, 5

    r = 2.0D+00
    p = i
    value = r8_power ( r, p )
    write ( *, '(2x,g14.6,i5,g14.6,i5)' ) r, p, value

  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests R8_POWER_FAST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) mults
  integer ( kind = 4 ) p
  real ( kind = 8 ) r
  real ( kind = 8 ) rp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  R8_POWER_FAST computes R**P, economizing on'
  write ( *, '(a)' ) '  multiplications.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      R          P       R**P       Mults'
  write ( *, '(a)' ) ' '

  do i = -10, 40

    r = 2.0D+00
    p = i
    call r8_power_fast ( r, p, rp, mults )
    write ( *, '(2x,g14.6,i5,g14.6,i5)' ) r, p, rp, mults

  end do

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests R8_ROUND2.
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nplace
  real ( kind = 8 ) r8_pi
  real ( kind = 8 ) x
  real ( kind = 8 ) xround

  x = r8_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020'
  write ( *, '(a)' ) '  R8_ROUND2 rounds a number to a'
  write ( *, '(a)' ) '  specified number of base 2 digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test effect on PI:'
  write ( *, '(a,g24.16)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPLACE  XROUND'
  write ( *, '(a)' ) ' '

  do i = 0, 20
    nplace = i
    call r8_round2 ( nplace, x, xround )
    write ( *, '(2x,i8,g24.16)' ) i, xround
  end do

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests R8_ROUNDB.
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

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nplace
  real ( kind = 8 ) r8_pi
  real ( kind = 8 ) x
  real ( kind = 8 ) xround

  base = 3
  x = r8_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  R8_ROUNDB rounds a number to a '
  write ( *, '(a)' ) '  specified number of base BASE digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here, we will use BASE = ',base
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test effect on PI:'
  write ( *, '(a,g24.16)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPLACE  XROUND'
  write ( *, '(a)' ) ' '

  do i = 0, 20
    nplace = i
    call r8_roundb ( base, nplace, x, xround )
    write ( *, '(2x,i8,g24.16)' ) i, xround
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Try with a negative base:'
  x = 121.0D+00
  base = -3
  nplace = 3
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  Input quantity is X = ', x
  write ( *, '(a,i8)' ) '  to be rounded in base ', base

  do nplace = 1, 5

    call r8_roundb ( base, nplace, x, xround )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,g24.16)' ) '  Output value to ', nplace, &
      ' places is ', xround

  end do

  return
end
subroutine test022 ( )

!*****************************************************************************80
!
!! TEST022 tests R8_ROUNDX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r8_pi
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nplace
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) xround

  seed = 123456789
  x = r8_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022'
  write ( *, '(a)' ) '  R8_ROUNDX rounds a number to a '
  write ( *, '(a)' ) '  specified number of decimal digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test effect on PI:'
  write ( *, '(a,g24.16)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPLACE  XROUND'
  write ( *, '(a)' ) ' '

  do i = 0, 10
    nplace = i
    call r8_roundx ( nplace, x, xround )
    write ( *, '(2x,i8,g24.16)' ) i, xround
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test effect on random values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NPLACE  X     XROUND'
  write ( *, '(a)' ) ' '

  do i = 1, 5

    x = r8_uniform_01 ( seed )

    write ( *, '(a)' ) ' '

    do nplace = 0, 10, 2
      call r8_roundx ( nplace, x, xround )
      write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) nplace, x, xround
    end do

  end do

  return
end
subroutine test023 ( )

!*****************************************************************************80
!
!! TEST023 tests R8_SIGN.
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

  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) r8_sign
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter, dimension ( test_num ) :: x_test = (/ &
    -1.25D+00, -0.25D+00, 0.0D+00, +0.5D+00, +9.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023'
  write ( *, '(a)' ) '  R8_SIGN returns the sign of a number.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = x_test(test)
    write ( *, '(2x,f8.4,2x,f8.4)' ) x, r8_sign ( x )
  end do

  return
end
subroutine test0235 ( )

!*****************************************************************************80
!
!! TEST0235 tests R8_SWAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0235'
  write ( *, '(a)' ) '  R8_SWAP swaps two reals.'

  x = 1.0D+00
  y = 3.141592653589793D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Before swapping:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    X = ', x
  write ( *, '(a,g14.6)' ) '    Y = ', y

  call r8_swap ( x, y )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After swapping:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    X = ', x
  write ( *, '(a,g14.6)' ) '    Y = ', y

  return
end
subroutine test024 ( )

!*****************************************************************************80
!
!! TEST024 tests R8_TO_R8_DISCRETE.
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

  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) :: ndx = 19
  real ( kind = 8 ) r
  real ( kind = 8 ) rd
  real ( kind = 8 ) :: rhi = 10.0D+00
  real ( kind = 8 ) rhi2
  real ( kind = 8 ) :: rlo = 1.0D+00
  real ( kind = 8 ) rlo2
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 15

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024'
  write ( *, '(a)' ) '  R8_TO_R8_DISCRETE maps numbers to a discrete set'
  write ( *, '(a)' ) '  of equally spaced numbers in an interval.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of discrete values = ', ndx
  write ( *, '(a,2g14.6)' ) '  Real interval: ', rlo, rhi
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R   RD'
  write ( *, '(a)' ) ' '

  seed = 123456789

  rlo2 = rlo - 2.0D+00
  rhi2 = rhi + 2.0D+00

  do test = 1, test_num
    r = r8_uniform_ab ( rlo2, rhi2, seed )
    call r8_to_r8_discrete ( r, rlo, rhi, ndx, rd )
    write ( *, '(2x,g14.6,g14.6)' ) r, rd
  end do

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests R8_TO_I4.
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

  integer ( kind = 4 ) ix
  integer ( kind = 4 ) ixmax
  integer ( kind = 4 ) ixmin
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  R8_TO_I4 finds an integer IX in [IXMIN,IXMAX]'
  write ( *, '(a)' ) '  corresponding to X in [XMIN,XMAX].'

  xmin = 2.5D+00
  x = 3.5D+00
  xmax = 5.5D+00

  ixmin = 10
  ixmax = 40

  call r8_to_i4 ( x, xmin, xmax, ixmin, ixmax, ix )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6,a,g14.6)' ) &
    '  XMIN ',  xmin, '   X = ',  x, '  XMAX = ', xmax
  write ( *, '(a,i14,a,i14,a,i14)' ) &
    ' IXMIN ', ixmin, '  IX = ', ix, ' IXMAX = ', ixmax

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests R8_UNIFORM.
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

  real ( kind = 8 ), parameter :: b = 10.0D+00
  real ( kind = 8 ), parameter :: c = 20.0D+00
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  R8_UNIFORM returns random values in a given range:'
  write ( *, '(a)' ) '  [ B, C ]'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this problem:'
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a,g14.6)' ) '  C = ', c
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 10
    r = r8_uniform_ab ( b, c, seed )
    write ( *, '(2x,g14.6)' ) r
  end do

  return
end
subroutine test027 ( )

!*****************************************************************************80
!
!! TEST027 tests R8_UNIFORM_01.
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

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_old
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  R8_UNIFORM_01 produces a sequence of random values.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Using random seed ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED   R8_UNIFORM_01(SEED)'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    seed_old = seed
    x = r8_uniform_01 ( seed )
    write ( *, '(2x,i12,2x,g14.6)' ) seed, x
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the sequence can be restarted.'
  write ( *, '(a)' ) '  Set the seed back to its original value, and see that'
  write ( *, '(a)' ) '  we generate the same sequence.'

  seed = 123456789
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED   R8_UNIFORM_01(SEED)'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_old = seed
    x = r8_uniform_01 ( seed )
    write ( *, '(2x,i12,2x,g14.6)' ) seed, x
  end do

  return
end
subroutine test028 ( )

!*****************************************************************************80
!
!! TEST028 tests R8_UNIFORM_01
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

  integer ( kind = 4 ) i
  real ( kind = 8 ) mean
  integer ( kind = 4 ), parameter :: n = 1000
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028'
  write ( *, '(a)' ) '  R8_UNIFORM_01 samples a uniform random'
  write ( *, '(a)' ) '  distribution in [0,1].'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Starting with seed = ', seed

  do i = 1, n
    x(i) = r8_uniform_01 ( seed )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First few values:'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
  end do

  call r8vec_mean ( n, x, mean )

  call r8vec_variance ( n, x, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of values computed was N = ', n
  write ( *, '(a,g14.6)' ) '  Average value was ', mean
  write ( *, '(a,g14.6)' ) '  Minimum value was ', minval ( x(1:n) )
  write ( *, '(a,g14.6)' ) '  Maximum value was ', maxval ( x(1:n) )
  write ( *, '(a,g14.6)' ) '  Variance was ', variance

  return
end
subroutine test029 ( )

!*****************************************************************************80
!
!! TEST029 tests R8_WALSH_1D;
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

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_walsh_1d
  real ( kind = 8 ) w0
  real ( kind = 8 ) wm1
  real ( kind = 8 ) wm2
  real ( kind = 8 ) wm3
  real ( kind = 8 ) wp1
  real ( kind = 8 ) wp2
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029'
  write ( *, '(a)' ) '  R8_WALSH_1D evaluates 1D Walsh functions:'
  write ( *, '(a)' ) ' '
  write ( *, * ) 'X  W(+2) W(+1) W(0) W(-1) W(-2) W(-3)'
  write ( *, '(a)' ) ' '

  do i = 0, 32

    x = real ( i, kind = 8 ) / 4.0D+00

    wp2 = r8_walsh_1d ( x,  2 )
    wp1 = r8_walsh_1d ( x,  1 )
    w0  = r8_walsh_1d ( x,  0 )
    wm1 = r8_walsh_1d ( x, -1 )
    wm2 = r8_walsh_1d ( x, -2 )
    wm3 = r8_walsh_1d ( x, -3 )

    write ( *, '(2x,f10.6,6f4.1)' ) x, wp2, wp1, w0, wm1, wm2, wm3

  end do

  return
end
subroutine test0295 ( )

!*****************************************************************************80
!
!! TEST0295 tests R8_WRAP;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = - 2.0D+00
  real ( kind = 8 ), parameter :: b = 12.0D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) r8_wrap
  real ( kind = 8 ), parameter :: rhi = 6.5D+00
  real ( kind = 8 ), parameter :: rlo = 3.0D+00
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0295'
  write ( *, '(a)' ) '  R8_WRAP "wraps" an R8 to lie within an interval:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Wrapping interval is ', rlo, ', ', rhi
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      R      R8_WRAP ( R )'
  write ( *, '(a)' ) ' '
  seed = 123456789

  do test = 1, test_num

    r = r8_uniform_ab ( a, b, seed )
    r2 = r8_wrap ( r, rlo, rhi )
    write ( *, '(2x,g14.6,2x,g14.6)' ) r, r2

  end do

  return
end
subroutine test031 ( )

!*****************************************************************************80
!
!! TEST031 tests R82POLY2_TYPE.
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

  integer ( kind = 4 ), parameter :: test_num = 12

  real ( kind = 8 ) a
  real ( kind = 8 ), save, dimension ( test_num ) :: a_test = (/  &
    9.0D+00, 4.0D+00, 9.0D+00,  1.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00,  0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), save, dimension ( test_num ) :: b_test = (/ &
    -4.0D+00,   1.0D+00,  16.0D+00,   1.0D+00, 0.0D+00,  &
     2.0D+00, 1.0D+00,   1.0D+00,  1.0D+00,  0.0D+00, &
     0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) c
  real ( kind = 8 ), save, dimension ( test_num ) :: c_test = (/  &
     0.0D+00,  -4.0D+00,   0.0D+00,   0.0D+00, 1.0D+00,  &
     0.0D+00, 0.0D+00,   0.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) d
  real ( kind = 8 ), save, dimension ( test_num ) :: d_test = (/ &
    -36.0D+00,  3.0D+00,  36.0D+00,  -6.0D+00, 3.0D+00, &
    -2.0D+00, 0.0D+00,   0.0D+00,  0.0D+00,  2.0D+00, &
     0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) e
  real ( kind = 8 ), save, dimension ( test_num ) :: e_test = (/ &
    -24.0D+00, -4.0D+00, -32.0D+00, -10.0D+00, -1.0D+00, &
     16.0D+00, -6.0D+00, -6.0D+00, -2.0D+00, -1.0D+00, &
     0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) f
  real ( kind = 8 ), save, dimension ( test_num ) :: f_test = (/ &
    -36.0D+00,  1.0D+00, -92.0D+00, 115.0D+00, -3.0D+00, &
     33.0D+00, +8.0D+00, 10.0D+00,  +1.0D+00,  1.0D+00, &
      0.0D+00, 1.0D+00 /)
  integer ( kind = 4 ) test
  integer ( kind = 4 ) type

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST031'
  write ( *, '(a)' ) '  R82POLY2_TYPE determines the type of a second order'
  write ( *, '(a)' ) '  equation in two variables.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a = a_test(test)
    b = b_test(test)
    c = c_test(test)
    d = d_test(test)
    e = e_test(test)
    f = f_test(test)

    write ( *, '(a)' ) ' '

    call r82poly2_print ( a, b, c, d, e, f )

    call r82poly2_type ( a, b, c, d, e, f, type )

    write ( *, '(a,i8)' ) '  Type = ', type

    call r82poly2_type_print ( type )

  end do

  return
end
subroutine test032 ( )

!*****************************************************************************80
!
!! TEST032 tests R82VEC_ORDER_TYPE.
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

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 10

  integer ( kind = 4 ) order
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  character ( len = 40 ) title
  real ( kind = 8 ) x(2,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032'
  write ( *, '(a)' ) '  R82VEC_ORDER_TYPE classifies a R8VEC as'
  write ( *, '(a)' ) '  -1: no order'
  write ( *, '(a)' ) '   0: all equal;'
  write ( *, '(a)' ) '   1: ascending;'
  write ( *, '(a)' ) '   2: strictly ascending;'
  write ( *, '(a)' ) '   3: descending;'
  write ( *, '(a)' ) '   4: strictly descending.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call r8mat_uniform_01 ( 2, n, seed, x )

    x(1:2,1:n) = real ( nint ( 3.0D+00 * x(1:2,1:n) ), kind = 8 )

    call r82vec_order_type ( n, x, order )

    write ( title, '(a,i8)' ) '  Order type = ', order

    call r82vec_print ( n, x, title )

  end do

  return
end
subroutine test033 ( )

!*****************************************************************************80
!
!! TEST033 tests R82VEC_PART_QUICK_A.
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

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) :: b = 0.0D+00
  real ( kind = 8 ) :: c = 10.0D+00
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST033'
  write ( *, '(a)' ) '  R82VEC_PART_QUICK_A reorders an R82VEC'
  write ( *, '(a)' ) '  as part of a quick sort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r8mat_uniform_ab ( 2, n, b, c, seed, a )

  call r82vec_print ( n, a, '  Before rearrangement:' )

  call r82vec_part_quick_a ( n, a, l, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rearranged array'
  write ( *, '(a,i8)' ) '  Left index =  ', l
  write ( *, '(a,i8)' ) '  Key index =   ', l+1
  write ( *, '(a,i8)' ) '  Right index = ', r
  write ( *, '(a)' ) ' '

  call r82vec_print ( l,     a(1:2,1:l),   '  Left half:' )
  call r82vec_print ( 1,     a(1:2,l+1),   '  Key:' )
  call r82vec_print ( n-l-1, a(1:2,l+2:n), '  Right half:' )

  return
end
subroutine test034 ( )

!*****************************************************************************80
!
!! TEST034 tests R82VEC_SORT_HEAP_INDEX_A.
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

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) :: b = 0.0D+00
  real ( kind = 8 ) :: c = 10.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST034'
  write ( *, '(a)' ) '  R82VEC_SORT_HEAP_INDEX_A index sorts an R82VEC'
  write ( *, '(a)' ) '  using heapsort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r8mat_uniform_ab ( 2, n, b, c, seed, a )
!
!  Give a few elements the same first component.
!
  a(1,3) = a(1,5)
  a(1,4) = a(1,12)
!
!  Give a few elements the same second component.
!
  a(2,6) = a(2,1)
  a(2,2) = a(2,9)
!
!  Make two entries equal.
!
  a(1:2,7) = a(1:2,11)

  call r82vec_print ( n, a, '  Before rearrangement:' )

  call r82vec_sort_heap_index_a ( n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I Index A(Index)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,i8,g14.6,g14.6)' ) i, indx(i), a(1:2,indx(i))
  end do

  call r82vec_permute ( n, indx, a )

  call r82vec_print ( n, a, '  After rearrangement by R82VEC_PERMUTE:' )

  return
end
subroutine test035 ( )

!*****************************************************************************80
!
!! TEST035 tests R82VEC_SORT_QUICK_A.
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

  integer ( kind = 4 ), parameter :: n = 12

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) :: b = 0.0D+00
  real ( kind = 8 ) :: c = 10.0D+00
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  R82VEC_SORT_QUICK_A sorts an R82VEC'
  write ( *, '(a)' ) '  using quick sort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r8mat_uniform_ab ( 2, n, b, c, seed, a )
!
!  Give a few elements the same first component.
!
  a(1,3) = a(1,5)
  a(1,4) = a(1,12)
!
!  Give a few elements the same second component.
!
  a(2,6) = a(2,1)
  a(2,2) = a(2,9)
!
!  Make two entries equal.
!
  a(1:2,7) = a(1:2,11)

  call r82vec_print ( n, a, '  Before rearrangement:' )

  call r82vec_sort_quick_a ( n, a )

  call r82vec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test036 ( )

!*****************************************************************************80
!
!! TEST036 tests R8BLOCK_EXPAND_LINEAR.
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

  integer ( kind = 4 ), parameter :: l = 4
  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: lfat = 1
  integer ( kind = 4 ), parameter :: mfat = 2
  integer ( kind = 4 ), parameter :: nfat = 1
  integer ( kind = 4 ), parameter :: l2 = ( l - 1 ) * ( lfat + 1 ) + 1
  integer ( kind = 4 ), parameter :: m2 = ( m - 1 ) * ( mfat + 1 ) + 1
  integer ( kind = 4 ), parameter :: n2 = ( n - 1 ) * ( nfat + 1 ) + 1
  real ( kind = 8 ), dimension(l,m,n) :: x = reshape ( (/ &
        1.0D+00,  2.0D+00,  3.0D+00,   4.0D+00,  1.0D+00, &
        4.0D+00,  9.0D+00, 16.0D+00,   1.0D+00,  8.0D+00, &
       27.0D+00, 64.0D+00,  2.0D+00,   4.0D+00,  6.0D+00, &
        8.0D+00,  2.0D+00,  8.0D+00,  18.0D+00, 32.0D+00, &
        2.0D+00, 16.0D+00, 54.0D+00, 128.0D+00 /), (/ l, m, n /) )
  real ( kind = 8 ) xfat(l2,m2,n2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036'
  write ( *, '(a)' ) '  R8BLOCK_EXPAND_LINEAR linearly interpolates new data'
  write ( *, '(a)' ) '  between old values in a 3D block.'

  call r8block_print ( l, m, n, x, '  Original block:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LFAT = ', lfat
  write ( *, '(a,i8)' ) '  MFAT = ', mfat
  write ( *, '(a,i8)' ) '  NFAT = ', nfat

  call r8block_expand_linear ( l, m, n, x, lfat, mfat, nfat, xfat )

  call r8block_print ( l2, m2, n2, xfat, '  Fattened block:' )

  return
end
subroutine test0365 ( )

!*****************************************************************************80
!
!! TEST0365 tests R8BLOCK_PRINT.
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
  real ( kind = 8 ), dimension(l,m,n) :: x = reshape ( (/ &
        1.0D+00,  2.0D+00,  3.0D+00,   4.0D+00,  1.0D+00, &
        4.0D+00,  9.0D+00, 16.0D+00,   1.0D+00,  8.0D+00, &
       27.0D+00, 64.0D+00,  2.0D+00,   4.0D+00,  6.0D+00, &
        8.0D+00,  2.0D+00,  8.0D+00,  18.0D+00, 32.0D+00, &
        2.0D+00, 16.0D+00, 54.0D+00, 128.0D+00 /), (/ l, m, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0365'
  write ( *, '(a)' ) '  R8BLOCK_PRINT prints an R8BLOCK.'

  call r8block_print ( l, m, n, x, '  The 3D array:' )

  return
end
subroutine test037 ( )

!*****************************************************************************80
!
!! TEST037 tests R8COL_FIND.
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

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: m = 3

  real ( kind = 8 ) dtab(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8vec(m)

  k = 1
  do i = 1, m
    do j = 1, n

      dtab(i,j) = real ( k, kind = 8 )

      if ( j == 3 ) then
        r8vec(i) = real ( k, kind = 8 )
      end if

      k = k + 1

    end do
  end do

  call r8col_find ( m, n, dtab, r8vec, icol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037'
  write ( *, '(a)' ) '  For an R8COL;'
  write ( *, '(a)' ) '  R8COL_FIND seeks a column matching given data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  R8COL_FIND returns ICOL = ', icol

  return
end
subroutine test038 ( )

!*****************************************************************************80
!
!! TEST038 tests R8COL_INSERT and R8COL_SORT_HEAP_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n_max = 10

  real ( kind = 8 ), dimension (m,n_max) :: a = reshape ( (/ &
    2.0D+00,  6.0D+00, 10.0D+00, &
    4.0D+00,  8.0D+00, 12.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00 /), (/ m, n_max /) )
  integer ( kind = 4 ) col
  real ( kind = 8 ), dimension(m) :: r8vec1 = (/ 3.0D+00, 7.0D+00, 11.0D+00 /)
  real ( kind = 8 ), dimension(m) :: r8vec2 = (/ 3.0D+00, 4.0D+00, 18.0D+00 /)
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST038'
  write ( *, '(a)' ) '  For an R8COL;'
  write ( *, '(a)' ) '  R8COL_SORT_HEAP_A does an ascending heap sort'
  write ( *, '(a)' ) '  R8COL_INSERT inserts new columns.'

  n = 4

  call r8mat_print ( m, n, a, '  The unsorted matrix:' )

  call r8col_sort_heap_a ( m, n, a )

  call r8mat_print ( m, n, a, '  The sorted matrix:' )

  call r8vec_print ( m, r8vec1, '  New column:' )

  call r8col_insert ( n_max, m, n, a, r8vec1, col )

  if ( col < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The data was already in column ', abs ( col )
  else
    call r8mat_print ( m, n, a, '  The updated matrix:' )
  end if

  call r8vec_print ( m, r8vec2, '  New column:' )

  call r8col_insert ( n_max, m, n, a, r8vec2, col )

  if ( col < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The data was already in column ', abs ( col )
  else
    call r8mat_print ( m, n, a, '  The updated matrix:' )
  end if

  return
end
subroutine test0383 ( )

!*****************************************************************************80
!
!! TEST0383 tests R8COL_PART_QUICK_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n = 8

  real ( kind = 8 ), dimension ( m, n ) :: a = reshape ( (/ &
     2.0D+00, 4.0D+00, &
     8.0D+00, 8.0D+00, &
     6.0D+00, 2.0D+00, &
     0.0D+00, 2.0D+00, &
    10.0D+00, 6.0D+00, &
    10.0D+00, 0.0D+00, &
     0.0D+00, 6.0D+00, &
     5.0D+00, 8.0D+00 /), (/ m, n /) )
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0383'
  write ( *, '(a)' ) '  For an R8COL;'
  write ( *, '(a)' ) '  R8COL_PART_QUICK_A partitions the matrix.'

  call r8mat_print ( m, n, a, '  The matrix:' )

  l = 2
  r = 4
  call r8col_part_quick_a ( m, n, a, l, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  L = ', l
  write ( *, '(a,i4)' ) '  R = ', r

  call r8mat_print ( m, n, a, '  The partitioned matrix:' )

  return
end
subroutine test0385 ( )

!*****************************************************************************80
!
!! TEST0385 tests R8COL_SORT_HEAP_INDEX_A.
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 15

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
    2.0D+00,  6.0D+00, 10.0D+00, &
    4.0D+00,  8.0D+00, 12.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  6.0D+00,  0.0D+00, &
    3.0D+00,  4.0D+00, 18.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  6.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    1.0D+00,  5.0D+00,  9.1D+00, &
    1.0D+00,  5.1D+00,  9.0D+00 /), (/ m, n /) )
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0385'
  write ( *, '(a)' ) '  R8COL_SORT_HEAP_INDEX_A computes an index vector which'
  write ( *, '(a)' ) '  ascending sorts an R8COL.'

  call r8mat_transpose_print ( m, n, a, '  The unsorted R8COL (transposed):' )

  call r8col_sort_heap_index_a ( m, n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The implicitly sorted R8COL (transposed)'
  write ( *, '(a)' ) ' '

  do j = 1, n
    j2 = indx(j)
    write ( *, '(2x,i4,a,2x,f10.1,2x,f10.1,2x,f10.1)' ) j2, ':', a(1:m,j2)
  end do

  return
end
subroutine test039 ( )

!*****************************************************************************80
!
!! TEST039 tests R8COL_SORT_QUICK_A.
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ), parameter :: b = 0.0D+00
  real ( kind = 8 ), parameter :: c = 10.0D+00
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST039'
  write ( *, '(a)' ) '  For an R8COL;'
  write ( *, '(a)' ) '  R8COL_SORT_QUICK_A does a quicksort.'

  seed = 123456789

  call r8mat_uniform_ab ( m, n, b, c, seed, a )

  call r8mat_print ( m, n, a, '  The unsorted matrix:' )

  call r8col_sort_quick_a ( m, n, a )

  call r8mat_print ( m, n, a, '  The sorted matrix:' )

  return
end
subroutine test0393 ( )

!*****************************************************************************80
!
!! TEST0393 tests R8COL_SORTED_TOL_UNIQUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 22

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
    1.9D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    4.0D+00,  8.0D+00, 12.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  6.0D+00,  0.0D+00, &
    2.0D+00,  0.0D+00, 10.1D+00, &
    2.0D+00,  0.1D+00, 10.0D+00, &
    3.0D+00,  4.0D+00, 18.0D+00, &
    1.9D+00,  8.0D+00, 10.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  6.0D+00, 10.0D+00, &
    2.1D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    2.0D+00,  0.0D+00, 10.1D+00, &
    1.0D+00,  5.0D+00,  9.1D+00, &
    1.0D+00,  5.1D+00,  9.0D+00 /), (/ m, n /) )
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) tol
  integer ( kind = 4 ) unique_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0393'
  write ( *, '(a)' ) '  R8COL_SORTED_TOL_UNIQUE finds tolerably unique columns'
  write ( *, '(a)' ) '  in a sorted R8COL.'

  call r8mat_transpose_print ( m, n, a, '  The unsorted R8COL (transposed):' )

  call r8col_sort_heap_a ( m, n, a )

  call r8mat_transpose_print ( m, n, a, '  The sorted R8COL (transposed):' )

  tol = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Using tolerance = ', tol

  call r8col_sorted_tol_unique ( m, n, a, tol, unique_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of tolerably unique columns is ', unique_num

  call r8mat_transpose_print ( m, unique_num, a, &
    '  The sorted tolerably unique R8COL (transposed):' )

  return
end
subroutine test0395 ( )

!*****************************************************************************80
!
!! TEST0395 tests R8COL_SORTED_TOL_UNIQUE_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 22

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
    1.9D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    4.0D+00,  8.0D+00, 12.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  6.0D+00,  0.0D+00, &
    2.0D+00,  0.0D+00, 10.1D+00, &
    2.0D+00,  0.1D+00, 10.0D+00, &
    3.0D+00,  4.0D+00, 18.0D+00, &
    1.9D+00,  8.0D+00, 10.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  6.0D+00, 10.0D+00, &
    2.1D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    2.0D+00,  0.0D+00, 10.1D+00, &
    1.0D+00,  5.0D+00,  9.1D+00, &
    1.0D+00,  5.1D+00,  9.0D+00 /), (/ m, n /) )
  integer ( kind = 4 ) n2
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) tol
  integer ( kind = 4 ) unique_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0395'
  write ( *, '(a)' ) '  R8COL_SORTED_TOL_UNIQUE_COUNT counts tolerably '
  write ( *, '(a)' ) '  unique columns in a sorted R8COL.'

  call r8mat_transpose_print ( m, n, a, '  The unsorted R8COL (transposed):' )

  call r8col_sort_heap_a ( m, n, a )

  call r8mat_transpose_print ( m, n, a, '  The sorted R8COL (transposed):' )

  tol = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Using tolerance = ', tol

  call r8col_sorted_tol_unique_count ( m, n, a, tol, unique_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of tolerably unique columns is ', unique_num

  return
end
subroutine test0397 ( )

!*****************************************************************************80
!
!! TEST0397 tests R8COL_SORTED_TOL_UNDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 22

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
    1.9D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    4.0D+00,  8.0D+00, 12.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  6.0D+00,  0.0D+00, &
    2.0D+00,  0.0D+00, 10.1D+00, &
    2.0D+00,  0.1D+00, 10.0D+00, &
    3.0D+00,  4.0D+00, 18.0D+00, &
    1.9D+00,  8.0D+00, 10.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  6.0D+00, 10.0D+00, &
    2.1D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    2.0D+00,  0.0D+00, 10.1D+00, &
    1.0D+00,  5.0D+00,  9.1D+00, &
    1.0D+00,  5.1D+00,  9.0D+00 /), (/ m, n /) )
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: au
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n_unique
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) tol
  integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
  integer ( kind = 4 ) xdnu(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0397'
  write ( *, '(a)' ) '  R8COL_SORTED_TOL_UNDEX produces index vectors which '
  write ( *, '(a)' ) '  create a sorted list of the tolerably unique columns'
  write ( *, '(a)' ) '  of a sorted R8COL,'
  write ( *, '(a)' ) '  and a map from the original R8COL to the (implicit)'
  write ( *, '(a)' ) '  R8COL of sorted tolerably unique elements.'

  call r8mat_transpose_print ( m, n, a, '  The unsorted R8COL (transposed):' )

  call r8col_sort_heap_a ( m, n, a )

  call r8mat_transpose_print ( m, n, a, '  The sorted R8COL (transposed):' )

  tol = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Tolerance for equality = ', tol

  call r8col_sorted_tol_unique_count ( m, n, a, tol, n_unique )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of tolerably unique columns is ', n_unique

  allocate ( au(1:m,1:n_unique) )
  allocate ( undx(1:n_unique) )

  call r8col_sorted_tol_undex ( m, n, a, n_unique, tol, undx, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XDNU points to the representative for each item.'
  write ( *, '(a)' ) '  UNDX selects the representatives.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  XDNU  UNDX'
  write ( *, '(a)' ) ' '
  do i = 1, n_unique
    write ( *, '(2x,i4,2x,i4,2x,i4)' ) i, xdnu(i), undx(i)
  end do
  do i = n_unique + 1, n
    write ( *, '(2x,i4,2x,i4)'       ) i, xdnu(i)
  end do

  do j = 1, n_unique
    au(1:m,j) = a(1:m,undx(j))
  end do

  call r8mat_transpose_print ( m, n_unique, au, &
    '  The tolerably unique R8COL (transposed):' )

  deallocate ( au )
  deallocate ( undx )

  return
end
subroutine test040 ( )

!*****************************************************************************80
!
!! TEST040 tests R8COL_MAX and R8COL_MIN;
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) amax(n)
  real ( kind = 8 ) amin(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST040'
  write ( *, '(a)' ) '  For an R8COL;'
  write ( *, '(a)' ) '  R8COL_MAX computes maximums;'
  write ( *, '(a)' ) '  R8COL_MIN computes minimums;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k, kind = 8 )
    end do
  end do

  call r8mat_print ( m, n, a, '  The array:' )

  call r8col_max ( m, n, a, amax )

  call r8col_min ( m, n, a, amin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Column, maximum, minimum:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    write ( *, '(2x,i3,3x,f10.4,2x,f10.4)' ) j, amax(j), amin(j)
  end do

  return
end
subroutine test041 ( )

!*****************************************************************************80
!
!! TEST041 tests R8COL_MEAN and R8COL_SUM;
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) colsum(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) mean(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST041'
  write ( *, '(a)' ) '  For an R8COL;'
  write ( *, '(a)' ) '  R8COL_MEAN computes means;'
  write ( *, '(a)' ) '  R8COL_SUM computes sums;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k, kind = 8 )
    end do
  end do

  call r8mat_print ( m, n, a, '  The array:' )

  call r8col_sum ( m, n, a, colsum )

  call r8col_mean ( m, n, a, mean )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Column  sum, mean:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    write ( *, '(2x,i3,3x,f10.4,2x,f10.4)' ) j, colsum(j), mean(j)
  end do

  return
end
subroutine test0415 ( )

!*****************************************************************************80
!
!! TEST0415 tests R8COL_PERMUTE.
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
    11.0D+00, 21.0D+00, 31.0D+00, &
    12.0D+00, 22.0D+00, 32.0D+00, &
    13.0D+00, 23.0D+00, 33.0D+00, &
    14.0D+00, 24.0D+00, 34.0D+00, &
    15.0D+00, 25.0D+00, 35.0D+00 /), (/ m, n /) )
  integer ( kind = 4 ), dimension ( n ) :: perm = (/ 2, 4, 5, 1, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0415'
  write ( *, '(a)' ) '  R8COL_PERMUTE permutes an R8COL in place.'

  call r8mat_print ( m, n, a, '  A (unpermuted)' )

  call i4vec_print ( n, perm, '  The (column) permutation vector:' )

  call r8col_permute ( m, n, perm, a )

  call r8mat_print ( m, n, a, '  A (permuted)' )

  return
end
subroutine test042 ( )

!*****************************************************************************80
!
!! TEST042 tests R8COL_SORTR_A.
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

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ), parameter :: b = 0.0D+00
  real ( kind = 8 ), parameter :: c = 10.0D+00
  integer ( kind = 4 ) key
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST042'
  write ( *, '(a)' ) '  R8COL_SORTR_A is given an array, and reorders'
  write ( *, '(a)' ) '  it so that a particular column is sorted.'
  write ( *, '(a)' ) ' '

  key = 2
  write ( *, '(a,i8)' ) '  Here, the special column is ', key

  seed = 123456789

  call r8mat_uniform_ab ( m, n, b, c, seed, a )

  call r8mat_print ( m, n, a, '  Unsorted array:' )

  call r8col_sortr_a ( m, n, a, key )

  call r8mat_print ( m, n, a, '  Sorted array:' )

  return
end
subroutine test043 ( )

!*****************************************************************************80
!
!! TEST043 tests R8COL_SWAP;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol1
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST043'
  write ( *, '(a)' ) '  R8COL_SWAP swaps two columns of an R8COL;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k, kind = 8 )
    end do
  end do

  call r8mat_print ( m, n, a, '  The array:' )

  icol1 = 1
  icol2 = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Swap columns ', icol1, ' and ', icol2

  call r8col_swap ( m, n, a, icol1, icol2 )

  call r8mat_print ( m, n, a, '  The updated matrix:' )

  return
end
subroutine test044 ( )

!*****************************************************************************80
!
!! TEST044 tests R8COL_TO_R8VEC.
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST044'
  write ( *, '(a)' ) '  R8COL_TO_R8VEC converts an array of columns to a vector.'
  write ( *, '(a)' ) ' '

  do i = 1, m
    do j = 1, n
      a(i,j) = real ( 10 * i + j, kind = 8 )
    end do
  end do

  call r8mat_print ( m, n, a, '  The array of columns:' )

  call r8col_to_r8vec ( m, n, a, x )

  call r8vec_print ( m*n, x, '  The resulting vector of columns:' )

  return
end
subroutine test0442 ( )

!*****************************************************************************80
!
!! TEST0442 tests R8COL_TOL_UNDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 22

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
    1.9D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    4.0D+00,  8.0D+00, 12.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  6.0D+00,  0.0D+00, &
    2.0D+00,  0.0D+00, 10.1D+00, &
    2.0D+00,  0.1D+00, 10.0D+00, &
    3.0D+00,  4.0D+00, 18.0D+00, &
    1.9D+00,  8.0D+00, 10.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  6.0D+00, 10.0D+00, &
    2.1D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    2.0D+00,  0.0D+00, 10.1D+00, &
    1.0D+00,  5.0D+00,  9.1D+00, &
    1.0D+00,  5.1D+00,  9.0D+00 /), (/ m, n /) )
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: au
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n_unique
  real ( kind = 8 ) tol
  integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
  integer ( kind = 4 ), allocatable, dimension ( : ) :: xdnu

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0442'
  write ( *, '(a)' ) '  R8COL_TOL_UNDEX produces index vectors which '
  write ( *, '(a)' ) '  create a sorted list of the tolerably unique columns'
  write ( *, '(a)' ) '  of an R8COL,'
  write ( *, '(a)' ) '  and a map from the original R8COL to the (implicit)'
  write ( *, '(a)' ) '  R8COL of sorted tolerably unique elements.'

  call r8mat_transpose_print ( m, n, a, '  The unsorted R8COL (transposed):' )

  tol = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Tolerance for equality = ', tol

  call r8col_tol_unique_count ( m, n, a, tol, n_unique )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of tolerably unique columns is ', n_unique

  allocate ( au(1:m,1:n_unique) )
  allocate ( undx(1:n_unique) )
  allocate ( xdnu(1:n) )

  call r8col_tol_undex ( m, n, a, n_unique, tol, undx, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XDNU points to the representative for each item.'
  write ( *, '(a)' ) '  UNDX selects the representatives.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  XDNU  UNDX'
  write ( *, '(a)' ) ' '
  do i = 1, n_unique
    write ( *, '(2x,i4,2x,i4,2x,i4)' ) i, xdnu(i), undx(i)
  end do
  do i = n_unique + 1, n
    write ( *, '(2x,i4,2x,i4)'       ) i, xdnu(i)
  end do

  do j = 1, n_unique
    au(1:m,j) = a(1:m,undx(j))
  end do

  call r8mat_transpose_print ( m, n_unique, au, &
    '  The tolerably unique R8COL (transposed):' )

  deallocate ( au )
  deallocate ( undx )
  deallocate ( xdnu )

  return
end
subroutine test0443 ( )

!*****************************************************************************80
!
!! TEST0443 tests R8COL_UNDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 15

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
    2.0D+00,  6.0D+00, 10.0D+00, &
    4.0D+00,  8.0D+00, 12.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  6.0D+00,  0.0D+00, &
    3.0D+00,  4.0D+00, 18.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  6.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    1.0D+00,  5.0D+00,  9.1D+00, &
    1.0D+00,  5.1D+00,  9.0D+00 /), (/ m, n /) )
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: au
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n_unique
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
  integer ( kind = 4 ) xdnu(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0443'
  write ( *, '(a)' ) '  R8COL_UNDEX produces index vectors which create a sorted'
  write ( *, '(a)' ) '  list of the unique columns of an (unsorted) R8COL,'
  write ( *, '(a)' ) '  and a map from the original R8COL to the (implicit)'
  write ( *, '(a)' ) '  R8COL of sorted unique elements.'

  call r8mat_transpose_print ( m, n, a, '  The R8COL (transposed):' )

  call r8col_unique_count ( m, n, a, n_unique )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique columns is ', n_unique

  allocate ( au(1:m,1:n_unique) )
  allocate ( undx(1:n_unique) )

  call r8col_undex ( m, n, a, n_unique, undx, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XDNU points to the representative for each item.'
  write ( *, '(a)' ) '  UNDX selects the representatives.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  XDNU  UNDX'
  write ( *, '(a)' ) ' '
  do i = 1, n_unique
    write ( *, '(2x,i4,2x,i4,2x,i4)' ) i, xdnu(i), undx(i)
  end do
  do i = n_unique + 1, n
    write ( *, '(2x,i4,2x,i4)'       ) i, xdnu(i)
  end do

  do j = 1, n_unique
    au(1:m,j) = a(1:m,undx(j))
  end do

  call r8mat_transpose_print ( m, n_unique, au, &
    '  The Unique R8COL (transposed):' )

  deallocate ( au )
  deallocate ( undx )

  return
end
subroutine test0445 ( )

!*****************************************************************************80
!
!! TEST0445 tests R8COL_UNIQUE_COUNT.
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 15

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
    2.0D+00,  6.0D+00, 10.0D+00, &
    4.0D+00,  8.0D+00, 12.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  6.0D+00,  0.0D+00, &
    3.0D+00,  4.0D+00, 18.0D+00, &
    0.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  6.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    3.0D+00,  7.0D+00, 11.0D+00, &
    2.0D+00,  0.0D+00, 10.0D+00, &
    2.0D+00,  6.0D+00, 10.0D+00, &
    1.0D+00,  5.0D+00,  9.0D+00, &
    1.0D+00,  5.0D+00,  9.1D+00, &
    1.0D+00,  5.1D+00,  9.0D+00 /), (/ m, n /) )
  integer ( kind = 4 ) unique_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0445'
  write ( *, '(a)' ) '  R8COL_UNIQUE_COUNT counts unique columns'
  write ( *, '(a)' ) '  in an unsorted R8COL.'

  call r8mat_transpose_print ( m, n, a, '  The R8COL (transposed):' )

  call r8col_unique_count ( m, n, a, unique_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique columns is ', unique_num

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests R8COL_VARIANCE.
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) variance(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045'
  write ( *, '(a)' ) '  R8COL_VARIANCE computes variances of an R8COL;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k, kind = 8 )
    end do
  end do

  call r8mat_print ( m, n, a, '  The array:' )

  call r8col_variance ( m, n, a, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Column       variance:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    write ( *, '(2x,i8,2x,f10.4)' ) j, variance(j)
  end do

  return
end
subroutine test046 ( )

!*****************************************************************************80
!
!! TEST046 tests R8R8VEC_INDEX_INSERT_UNIQUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 30

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) :: x_max = 4.0D+00
  real ( kind = 8 ) :: x_min = 1.0D+00
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(n_max)
  real ( kind = 8 ) :: y_max = 3.0D+00
  real ( kind = 8 ) :: y_min = 1.0D+00
  real ( kind = 8 ) yval

  n = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST046'
  write ( *, '(a)' ) '  R8R8VEC_INDEX_INSERT_UNIQUE inserts unique values into'
  write ( *, '(a)' ) '  an index sorted array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Generate ', n_max, ' random values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index    XVAL    YVAL'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, n_max

    xval = r8_uniform_ab ( x_min, x_max, seed )
    xval = real ( nint ( xval ), kind = 8 )
    yval = r8_uniform_ab ( y_min, y_max, seed )
    yval = real ( nint ( yval ), kind = 8 )

    call r8r8vec_index_insert_unique ( n_max, n, x, y, indx, xval, yval, &
      ival, ierror )

    write ( *, '(2x,i3,6x,f6.2,9x,f6.2)' ) ival, xval, yval

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector of unique X Y values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  X(I)   Y(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,f6.2,9x,f6.2)' ) i, x(i), y(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, Y sorted by index'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(INDX(I))  Y(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,f6.2,9x,f6.2,9x,f6.2)' ) i, indx(i), &
      x(indx(i)), y(indx(i))
  end do

  return
end
subroutine test047 ( )

!*****************************************************************************80
!
!! TEST047 tests R8R8R8VEC_INDEX_INSERT_UNIQUE.
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

  integer ( kind = 4 ), parameter :: n_max = 30

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(n_max)
  real ( kind = 8 ) yval
  real ( kind = 8 ) z(n_max)
  real ( kind = 8 ) zval

  n = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST047'
  write ( *, '(a)' ) '  R8R8R8VEC_INDEX_INSERT_UNIQUE inserts unique values'
  write ( *, '(a)' ) '  into an index sorted array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of random values to generate = ', n_max
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    XVAL    YVAL  ZVAL  Index'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, n_max

    xval = r8_uniform_ab ( 1.0D+00, 4.0D+00, seed )
    xval = real ( nint ( xval ), kind = 8 )
    yval = r8_uniform_ab ( 1.0D+00, 3.0D+00, seed )
    yval = real ( nint ( yval ), kind = 8 )
    zval = r8_uniform_ab ( 1.0D+00, 4.0D+00, seed )
    zval = real ( nint ( zval ), kind = 8 )

    call r8r8r8vec_index_insert_unique ( n_max, n, x, y, z, indx, &
      xval, yval, zval, ival, ierror )

    write ( *, '(2x,i3,6x,f6.2,9x,f6.2,9x,f6.2)' ) ival, xval, yval, zval

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector of unique X Y Z values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  X(I)   Y(I)    Z(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,f6.2,9x,f6.2,9x,f6.2)' ) i, x(i), y(i), z(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X Y Z sorted by index:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,f6.2,9x,f6.2,9x,f6.2)' ) i, indx(i), &
      x(indx(i)), y(indx(i)), z(indx(i))
  end do

  return
end
subroutine test048 ( )

!*****************************************************************************80
!
!! TEST048 tests R8INT_TO_I4INT and I4INT_TO_R8INT;
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ir
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_uniform_ab
  real ( kind = 8 ) :: rhi = 200.0D+00
  real ( kind = 8 ) rhi2
  real ( kind = 8 ) :: rlo = 100.0D+00
  real ( kind = 8 ) rlo2
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 10

  ilo = 1
  ihi = 11

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST048'
  write ( *, '(a)' ) '  For data in an interval,'
  write ( *, '(a)' ) '  I4INT_TO_R8INT converts an integer to a real;'
  write ( *, '(a)' ) '  R8INT_TO_I4INT converts a real to an integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2i8)' ) '  Integer interval: ', ilo, ihi
  write ( *, '(a,2g14.6)' ) '  Real interval: ', rlo, rhi
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R   I(R)  R(I(R))'
  write ( *, '(a)' ) ' '

  seed = 123456789

  rlo2 = rlo - 15.0D+00
  rhi2 = rhi + 15.0D+00

  do test = 1, test_num
    r = r8_uniform_ab ( rlo2, rhi2, seed )
    call r8int_to_i4int ( rlo, rhi, r, ilo, ihi, ir )
    call i4int_to_r8int ( ilo, ihi, ir, rlo, rhi, r2 )
    write ( *, '(2x,g14.6,i8,g14.6)' ) r, ir, r2
  end do

  return
end
subroutine test049 ( )

!*****************************************************************************80
!
!! TEST049 tests R8MAT_CHOLESKY_FACTOR, R8MAT_CHORESKY_FACTOR and R8MAT_CHOLESKY_SOLVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) d(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  real ( kind = 8 ) l(n,n)
  real ( kind = 8 ) r(n,n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST049'
  write ( *, '(a)' ) '  For a positive definite symmetric matrix,'
  write ( *, '(a)' ) '  R8MAT_CHOLESKY_FACTOR computes the lower'
  write ( *, '(a)' ) '  triangular Cholesky factor.'
  write ( *, '(a)' ) '  R8MAT_CHORESKY_FACTOR computes the upper'
  write ( *, '(a)' ) '  triangular Cholesky factor.'
  write ( *, '(a)' ) '  R8MAT_CHOLESKY_SOLVE solves a linear system'
  write ( *, '(a)' ) '  using the Cholesky factorization.'

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a, '  Matrix to be factored:' )
!
!  Compute the lower Cholesky factor.
!
  call r8mat_cholesky_factor ( n, a, l, ierror )

  call r8mat_print ( n, n, l, '  Cholesky factor L:' )

  d(1:n,1:n) = matmul ( l(1:n,1:n), transpose ( l(1:n,1:n) ) )

  call r8mat_print ( n, n, d, '  Product L * L'':' )
!
!  Compute the upper Cholesky factor.
!
  call r8mat_choresky_factor ( n, a, r, ierror )

  call r8mat_print ( n, n, r, '  Cholesky factor R:' )

  d(1:n,1:n) = matmul ( r(1:n,1:n), transpose ( r(1:n,1:n) ) )

  call r8mat_print ( n, n, d, '  Product R * R'':' )
!
!  Solve a linear system.
!
  b(1:n-1) = 0.0D+00
  b(n) = real ( n + 1, kind = 8 )

  call r8vec_print ( n, b, '  Right hand side:' )

  call r8mat_cholesky_solve ( n, l, b, x )

  call r8vec_print ( n, x, '  Computed solution:' )

  return
end
subroutine test050 ( )

!*****************************************************************************80
!
!! TEST050 tests R8MAT_DET_2D;
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

  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_det_2d
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    1.0D+00, 10.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST050'
  write ( *, '(a)' ) '  R8MAT_DET_2D: determinant of a 2 by 2 matrix;'

  call r8mat_vand2 ( n, x, a )

  det = r8mat_det_2d ( a )

  call r8mat_print ( n, n, a, '  Matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8MAT_DET_2D computes determinant:', det
!
!  Special formula for the determinant of a Vandermonde matrix:
!
  det = 1.0D+00
  do i = 1, n
    do j = 1, i-1
      det = det * ( x(i) - x(j) )
    end do
  end do
  write ( *, '(a,g14.6)' ) '  Exact determinant is ', det

  return
end
subroutine test051 ( )

!*****************************************************************************80
!
!! TEST051 tests R8MAT_DET_3D;
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

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_det_3d
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    1.0D+00, 10.0D+00, 4.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST051'
  write ( *, '(a)' ) '  R8MAT_DET_3D: determinant of a 3 by 3 matrix;'

  call r8mat_vand2 ( n, x, a )
  det = r8mat_det_3d ( a )

  call r8mat_print ( n, n, a, '  Matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8MAT_DET_3D computes determinant:', det
!
!  Special formula for the determinant of a Vandermonde matrix:
!
  det = 1.0D+00
  do i = 1, n
    do j = 1, i-1
      det = det * ( x(i) - x(j) )
    end do
  end do
  write ( *, '(a,g14.6)' ) '  Exact determinant is ', det

  return
end
subroutine test052 ( )

!*****************************************************************************80
!
!! TEST052 tests R8MAT_DET_4D;
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

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_det_4d
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    1.0D+00, 10.0D+00, 4.0D+00, 2.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST052'
  write ( *, '(a)' ) '  R8MAT_DET_4D determinant of a 4 by 4 matrix;'

  call r8mat_vand2 ( n, x, a )
  det = r8mat_det_4d ( a )

  call r8mat_print ( n, n, a, '  Matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8MAT_DET_4D computes determinant:', det
!
!  Special formula for the determinant of a Vandermonde matrix:
!
  det = 1.0D+00
  do i = 1, n
    do j = 1, i-1
      det = det * ( x(i) - x(j) )
    end do
  end do
  write ( *, '(a,g14.6)' ) '  Exact determinant is ', det

  return
end
subroutine test053 ( )

!*****************************************************************************80
!
!! TEST053 tests R8MAT_DET_5D;
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_det_5d
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    1.0D+00, 10.0D+00, 4.0D+00, 2.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST053'
  write ( *, '(a)' ) '  R8MAT_DET_5D: determinant of 5 by 5 matrix.'

  call r8mat_vand2 ( n, x, a )
  det = r8mat_det_5d ( a )

  call r8mat_print ( n, n, a, '  Matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8MAT_DET_5D computes determinant: ', det
!
!  Special formula for the determinant of a Vandermonde matrix:
!
  det = 1.0D+00
  do i = 1, n
    do j = 1, i-1
      det = det * ( x(i) - x(j) )
    end do
  end do
  write ( *, '(a,g14.6)' ) '  Exact determinant is ', det

  return
end
subroutine test054 ( )

!*****************************************************************************80
!
!! TEST054 tests R8MAT_EXPAND_LINEAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: mfat = 2
  integer ( kind = 4 ), parameter :: nfat = 1
  integer ( kind = 4 ), parameter :: m2 = ( m - 1 ) * ( mfat + 1 ) + 1
  integer ( kind = 4 ), parameter :: n2 = ( n - 1 ) * ( nfat + 1 ) + 1
!
!  Each row of this definition is a COLUMN of the matrix.
!
  real ( kind = 8 ), dimension(m,n) :: x = reshape ( (/ &
    1.0D+00, 2.0D+00,  3.0D+00,  4.0D+00, &
    1.0D+00, 4.0D+00,  9.0D+00, 16.0D+00, &
    1.0D+00, 8.0D+00, 27.0D+00, 64.0D+00 /), (/ m, n /) )
  real ( kind = 8 ) xfat(m2,n2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST054'
  write ( *, '(a)' ) '  R8MAT_EXPAND_LINEAR linearly interpolates new data'
  write ( *, '(a)' ) '  between old values in a matrix.'

  call r8mat_print ( m, n, x, '  Original matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  MFAT = ', mfat
  write ( *, '(a,i8)' ) '  NFAT = ', nfat

  call r8mat_expand_linear ( m, n, x, mfat, nfat, xfat )

  call r8mat_print ( m2, n2, xfat, '  Fattened matrix:' )

  return
end
subroutine test055 ( )

!*****************************************************************************80
!
!! TEST055 tests R8MAT_EXPAND_LINEAR2.
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: m2 = 10
  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: n2 = 5

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a2(m2,n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST055'
  write ( *, '(a)' ) '  R8MAT_EXPAND_LINEAR2 fills in a large array by'
  write ( *, '(a)' ) '  interpolating data from a small array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original matrix has dimensions:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,2i8)' ) m, n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expanded matrix has dimensions:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,2i8)' ) m2, n2

  do i = 1, m
    do j = 1, n
      a(i,j) = 10.0D+00 * real ( i, kind = 8 ) + real ( j, kind = 8 )
    end do
  end do

  call r8mat_print ( m, n, a, '  The little matrix A:' )

  call r8mat_expand_linear2 ( m, n, a, m2, n2, a2 )

  call r8mat_print ( m2, n2, a2, '  Expanded array A2:' )

  return
end
subroutine test0555 ( )

!*****************************************************************************80
!
!! TEST0555 tests R8MAT_FSS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: nb = 3

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0555'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  R8MAT_FSS factors and solves multiple linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8mat_uniform_01 ( n, n, seed, a )
!
!  Set the desired solutions.
!
  x(1:n) = 1.0D+00
  call r8mat_mv ( n, n, a, x, b(1:n,1) )

  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do
  call r8mat_mv ( n, n, a, x, b(1:n,2) )

  do i = 1, n
    x(i) = real ( 1 + mod ( i - 1, 3 ), kind = 8 )
  end do
  call r8mat_mv ( n, n, a, x, b(1:n,3) )
!
!  Factor and solve the system.
!
  call r8mat_fss ( n, a, nb, b, info )
  
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST0555 - Fatal error!'
    write ( *, '(a)' ) '  R8MAT_FSS reports the matrix is singular.'
    return
  end if

  call r8mat_print ( n, nb, b, '  Solutions:' )

  return
end
subroutine test056 ( )

!*****************************************************************************80
!
!! TEST056 tests R8MAT_GIVENS_POST.
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

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) ag(n,n)
  integer ( kind = 4 ) col
  real ( kind = 8 ) g(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST056'
  write ( *, '(a)' ) '  R8MAT_GIVENS_POST computes a Givens ' // &
    'postmultiplier rotation matrix.'

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( i**(j-1), kind = 8 )
    end do
  end do

  call r8mat_print ( n, n, a, '  Matrix A:' )

  row = 3
  col = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,2i8)' ) '  I, J=', row, col

  call r8mat_givens_post ( n, a, row, col, g )

  call r8mat_print ( n, n, g, '  G' )

  ag(1:n,1:n) = matmul ( a(1:n,1:n), g(1:n,1:n) )

  call r8mat_print ( n, n, ag, '  A*G' )

  return
end
subroutine test057 ( )

!*****************************************************************************80
!
!! TEST057 tests R8MAT_GIVENS_PRE.
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

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) col
  real ( kind = 8 ) g(n,n)
  real ( kind = 8 ) ga(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST057'
  write ( *, '(a)' ) '  R8MAT_GIVENS_PRE computes a Givens ' // &
    'premultiplier rotation matrix.'

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( i**(j-1), kind = 8 )
    end do
  end do

  call r8mat_print ( n, n, a, '  Matrix A:' )

  row = 3
  col = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,2i8)' ) '  I, J=', row, col

  call r8mat_givens_pre ( n, a, row, col, g )

  call r8mat_print ( n, n, g, '  G' )

  ga(1:n,1:n) = matmul ( g(1:n,1:n), a(1:n,1:n) )

  call r8mat_print ( n, n, ga, '  G*A' )

  return
end
subroutine test058 ( )

!*****************************************************************************80
!
!! TEST058 tests R8MAT_HESS.
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

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) h(n,n)
  external test058_f
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    1.0D+00, 2.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST058'
  write ( *, '(a)' ) '  R8MAT_HESS estimates the Hessian matrix'
  write ( *, '(a)' ) '  of a scalar function.'

  call r8mat_hess ( test058_f, n, x, h )

  call r8mat_print ( n, n, h, '  Estimated jacobian:' )

  call test058_hess ( n, x, h )

  call r8mat_print ( n, n, h, '  Exact jacobian:' )

  return
end
subroutine test058_f ( n, x, f )

!*****************************************************************************80
!
!! TEST058_F is a sample nonlinear function for treatment by R8MAT_HESS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of parameters.
!
!    Input, real ( kind = 8 ) X(N), the parameter values.
!
!    Output, real ( kind = 8 ) F, the function value.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  real ( kind = 8 ) x(n)

  f = x(1)**2 + x(1) * x(2) + x(2) * cos ( 10.0D+00 * x(3) )

  return
end
subroutine test058_hess ( n, x, h )

!*****************************************************************************80
!
!! TEST058_HESS is the exact Hessian of TEST058_F.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of parameters.
!
!    Input, real ( kind = 8 ) X(N), the parameter values.
!
!    Output, real ( kind = 8 ) H(N,N), the Hessian values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) x(n)

  h(1,1) = 2.0D+00
  h(1,2) = 1.0D+00
  h(1,3) = 0.0D+00

  h(2,1) = 1.0D+00
  h(2,2) = 0.0D+00
  h(2,3) = - 10.0D+00 * sin ( 10.0D+00 * x(3) )

  h(3,1) = 0.0D+00
  h(3,2) = - 10.0D+00 * sin ( 10.0D+00 * x(3) )
  h(3,3) = - 100.0D+00 * x(2) * cos ( 10.0D+00 * x(3) )

  return
end
subroutine test059 ( )

!*****************************************************************************80
!
!! TEST059 tests R8MAT_HOUSE_FORM and R8VEC_HOUSE_COLUMN.
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

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), parameter :: b = 0.0D+00
  real ( kind = 8 ), parameter :: c = 5.0D+00
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) ha(n,n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST059'
  write ( *, '(a)' ) '  R8VEC_HOUSE_COLUMN returns the compact form of'
  write ( *, '(a)' ) '  a Householder matrix that "packs" a column'
  write ( *, '(a)' ) '  of a matrix.'
!
!  Get a random matrix.
!
  seed = 123456789

  call r8mat_uniform_ab ( n, n, b, c, seed, a )

  call r8mat_print ( n, n, a, '  Matrix A:' )

  do k = 1, n-1

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Working on column K = ', k

    call r8vec_house_column ( n, a(1,k), k, v )

    call r8mat_house_form ( n, v, h )

    call r8mat_print ( n, n, h, '  Householder matrix H:' )

    ha(1:n,1:n) = matmul ( h(1:n,1:n), a(1:n,1:n) )

    call r8mat_print ( n, n, ha, '  Product H*A:' )
!
!  If we set A := HA, then we can successively convert A to upper
!  triangular form.
!
    a(1:n,1:n) = ha(1:n,1:n)

  end do

  return
end
subroutine test060 ( )

!*****************************************************************************80
!
!! TEST060 tests R8MAT_HOUSE_FORM.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ), dimension ( n ) :: v = (/ &
    0.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST060'
  write ( *, '(a)' ) '  R8MAT_HOUSE_FORM forms a Householder'
  write ( *, '(a)' ) '  matrix from its compact form.'

  call r8vec_print ( n, v, '  Compact vector form V:' )

  call r8mat_house_form ( n, v, h )

  call r8mat_print ( n, n, h, '  Householder matrix H:' )

  return
end
subroutine test061 ( )

!*****************************************************************************80
!
!! TEST061 tests R8MAT_HOUSE_POST.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) ah(n,n)
  real ( kind = 8 ), parameter :: b = 0.0D+00
  real ( kind = 8 ), parameter :: c = 5.0D+00
  real ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) row
  integer ( kind = 4 ) col
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST061'
  write ( *, '(a)' ) '  R8MAT_HOUSE_POST computes a Householder'
  write ( *, '(a)' ) '  postmultiplier;'

  seed = 123456789

  call r8mat_uniform_ab ( n, n, b, c, seed, a )

  call r8mat_print ( n, n, a, '  Matrix A:' )

  row = 2
  col = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,2i8)' ) '  I, J=', row, col

  call r8mat_house_post ( n, a, row, col, h )

  call r8mat_print ( n, n, h, '  Householder matrix H:' )

  ah(1:n,1:n) = matmul ( a(1:n,1:n), h(1:n,1:n) )

  call r8mat_print ( n, n, ah, '  Product A*H:' )

  return
end
subroutine test062 ( )

!*****************************************************************************80
!
!! TEST062 tests R8MAT_HOUSE_PRE.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), parameter :: b = 0.0D+00
  real ( kind = 8 ), parameter :: c = 5.0D+00
  integer ( kind = 4 ) col
  real ( kind = 8 ) h(n,n)
  real ( kind = 8 ) ha(n,n)
  integer ( kind = 4 ) row
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST062'
  write ( *, '(a)' ) '  R8MAT_HOUSE_PRE computes a Householder premultiplier;'

  seed = 123456789

  call r8mat_uniform_ab ( n, n, b, c, seed, a )

  call r8mat_print ( n, n, a, '  Matrix A:' )

  row = 2
  col = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,2i8)' ) '  I, J=', row, col

  call r8mat_house_pre ( n, a, row, col, h )

  call r8mat_print ( n, n, h, '  Householder matrix H:' )

  ha(1:n,1:n) = matmul ( h(1:n,1:n), a(1:n,1:n) )

  call r8mat_print ( n, n, ha, '  Product H*A:' )

  return
end
subroutine test063 ( )

!*****************************************************************************80
!
!! TEST063 tests R8MAT_MAX_INDEX and R8MAT_MIN_INDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ), parameter :: b = 0.0D+00
  real ( kind = 8 ), parameter :: c = 10.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST063'
  write ( *, '(a)' ) '  R8MAT_MAX_INDEX locates the maximum entry of an R8COL;'
  write ( *, '(a)' ) '  R8MAT_MIN_INDEX locates the minimum entry of an R8COL;'

  seed = 123456789

  call r8mat_uniform_ab ( m, n, b, c, seed, a )

  call r8mat_print ( m, n, a, '  Random array:' )

  call r8mat_max_index ( m, n, a, i, j )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2i8)' ) '  Maximum I,J indices            ', i, j
  call r8mat_min_index ( m, n, a, i, j )
  write ( *, '(a,2i8)' ) '  Minimum I,J indices            ', i, j

  return
end
subroutine test064 ( )

!*****************************************************************************80
!
!! TEST064 tests R8MAT_INVERSE_2D.
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

  integer ( kind = 4 ), parameter :: n = 2
!
!  Each ROW of this definion is a COLUMN of the matrix.
!
  real ( kind = 8 ), dimension(n,n) :: a = reshape ( (/ &
    1.0D+00, 3.0D+00, &
    2.0D+00, 4.0D+00 /), (/ 2, 2 /) )
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  real ( kind = 8 ) det

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST064'
  write ( *, '(a)' ) '  R8MAT_INVERSE_2D inverts a 2 by 2 matrix.'

  call r8mat_print ( n, n, a, '  Matrix A to invert:' )
!
!  Compute the inverse matrix.
!
  call r8mat_inverse_2d ( a, b, det )

  if ( det == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input matrix was singular, no inverse'
    write ( *, '(a)' ) '  could be computed.'
    return
  end if

  call r8mat_print ( n, n, b, '  Inverse matrix B:' )

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call r8mat_print ( n, n, c, '  Product C = A * B:' )

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests R8MAT_INVERSE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
!
!  Each ROW of this definion is a COLUMN of the matrix.
!
  real ( kind = 8 ), dimension ( n, n ) :: a = reshape ( (/ &
    1.0D+00, 4.0D+00, 7.0D+00, &
    2.0D+00, 5.0D+00, 8.0D+00, &
    3.0D+00, 6.0D+00, 0.0D+00 /), (/ n, n /) )
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  real ( kind = 8 ) det

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  R8MAT_INVERSE_3D inverts a 3 by 3 matrix.'

  call r8mat_print ( n, n, a, '  Matrix A to be inverted:' )
!
!  Compute the inverse matrix.
!
  call r8mat_inverse_3d ( a, b, det )

  if ( det == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input matrix was singular, no inverse'
    write ( *, '(a)' ) '  could be computed.'
    return
  end if

  call r8mat_print ( n, n, b, '  Inverse matrix B:' )

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call r8mat_print ( n, n, c, '  Product C = A * B:' )

  return
end
subroutine test066 ( )

!*****************************************************************************80
!
!! TEST066 tests R8MAT_INVERSE_4D.
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

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST066'
  write ( *, '(a)' ) '  R8MAT_INVERSE_4D inverts a 4 x 4 matrix.'

  do i = 1, n
    do j = 1, n

      if ( i <= j ) then
        a(i,j) = real ( n + 1 - j, kind = 8 )
      else if ( j == i - 1 ) then
        a(i,j) = n - j
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  call r8mat_print ( n, n, a, '  Matrix A to be inverted:' )

  call r8mat_inverse_4d ( a, b, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Determinant is ', det

  call r8mat_print ( n, n, b, '  Inverse B:' )

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call r8mat_print ( n, n, c, '  Product C = A * B:' )

  return
end
subroutine test067 ( )

!*****************************************************************************80
!
!! TEST067 tests R8MAT_JAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) :: eps = 0.00001D+00
  real ( kind = 8 ) fprime(m,n)
  external test067_f
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST067'
  write ( *, '(a)' ) '  R8MAT_JAC estimates the M by N jacobian matrix'
  write ( *, '(a)' ) '  of a nonlinear function.'

  call r8mat_jac ( m, n, eps, test067_f, x, fprime )

  call r8mat_print ( m, n, fprime, '  Estimated jacobian:' )

  call test067_jac ( m, n, x, fprime )

  call r8mat_print (  m, n, fprime, '  Exact jacobian:' )

  return
end
subroutine test067_f ( m, n, x, f )

!*****************************************************************************80
!
!! TEST067_F is a sample nonlinear function for treatment by R8MAT_JAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of functions.
!
!    Input, integer N, the number of parameters.
!
!    Input, real ( kind = 8 ) X(N), the parameter values.
!
!    Output, real ( kind = 8 ) F(M), the function values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(m)
  real ( kind = 8 ) x(n)

  f(1) = sin ( x(1) * x(2) )
  f(2) = sqrt ( 1.0D+00 + x(1)**2 ) + x(3)
  f(3) = x(1) + 2.0D+00 * x(2) + 3.0D+00 * x(3) + 4.0D+00 * x(4)

  return
end
subroutine test067_jac ( m, n, x, fprime )

!*****************************************************************************80
!
!! TEST067_JAC is the exact jacobian of TEST067_F.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of functions.
!
!    Input, integer N, the number of parameters.
!
!    Input, real ( kind = 8 ) X(N), the parameter values.
!
!    Output, real ( kind = 8 ) FPRIME(M,N), the jacobian values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fprime(m,n)
  real ( kind = 8 ) x(n)

  fprime(1,1) = cos ( x(1) * x(2) ) * x(2)
  fprime(1,2) = cos ( x(1) * x(2) ) * x(1)
  fprime(1,3) = 0.0D+00
  fprime(1,4) = 0.0D+00

  fprime(2,1) = x(1) / sqrt ( 1.0D+00 + x(1)**2 )
  fprime(2,2) = 0.0D+00
  fprime(2,3) = 1.0D+00
  fprime(2,4) = 0.0D+00

  fprime(3,1) = 1.0D+00
  fprime(3,2) = 2.0D+00
  fprime(3,3) = 3.0D+00
  fprime(3,4) = 4.0D+00

  return
end
subroutine test068 ( )

!*****************************************************************************80
!
!! TEST068 tests R8MAT_L_INVERSE.
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

  integer ( kind = 4 ), parameter :: n = 4
!
!  Each row of this definition is a COLUMN of the matrix.
!
  real ( kind = 8 ), dimension ( n, n ) :: a = reshape ( (/ &
    1.0D+00, 2.0D+00, 4.0D+00,  7.0D+00, &
    0.0D+00, 3.0D+00, 5.0D+00,  8.0D+00, &
    0.0D+00, 0.0D+00, 6.0D+00,  9.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00, 10.0D+00 /), (/ n, n /) )
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST068'
  write ( *, '(a)' ) '  R8MAT_L_INVERSE inverts a lower triangular matrix.'

  call r8mat_print ( n, n, a, '  Matrix A to be inverted:' )

  call r8mat_l_inverse ( n, a, b )

  call r8mat_print ( n, n, b, '  Inverse matrix B:' )

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call r8mat_print ( n, n, c, '  Product C = A * B:' )

  return
end
subroutine test069 ( )

!*****************************************************************************80
!
!! TEST069 tests R8MAT_L_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), dimension(28) :: a1 = (/ &
    11.0D+00, 21.0D+00, 31.0D+00, 41.0D+00, 51.0D+00, 61.0D+00, 71.0D+00, &
              22.0D+00, 32.0D+00, 42.0D+00, 52.0D+00, 62.0D+00, 72.0D+00, &
                        33.0D+00, 43.0D+00, 53.0D+00, 63.0D+00, 73.0D+00, &
                                  44.0D+00, 54.0D+00, 64.0D+00, 74.0D+00, &
                                            55.0D+00, 65.0D+00, 75.0D+00, &
                                                      66.0D+00, 76.0D+00, &
                                                                77.0D+00 /)
  real ( kind = 8 ), dimension(18) :: a2 = (/ &
    11.0D+00, 21.0D+00, 31.0D+00, 41.0D+00, 51.0D+00, 61.0D+00, 71.0D+00, &
              22.0D+00, 32.0D+00, 42.0D+00, 52.0D+00, 62.0D+00, 72.0D+00, &
                        33.0D+00, 43.0D+00, 53.0D+00, 63.0D+00, 73.0D+00 /)
  real ( kind = 8 ), dimension(10) :: a3 = (/ &
    11.0D+00, 21.0D+00, 31.0D+00, 41.0D+00, &
              22.0D+00, 32.0D+00, 42.0D+00, &
                        33.0D+00, 43.0D+00, &
                                  44.0D+00 /)
  integer ( kind = 4 ), parameter :: m1 = 7
  integer ( kind = 4 ), parameter :: m2 = 7
  integer ( kind = 4 ), parameter :: m3 = 4
  integer ( kind = 4 ), parameter :: n1 = 7
  integer ( kind = 4 ), parameter :: n2 = 3
  integer ( kind = 4 ), parameter :: n3 = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST069'
  write ( *, '(a)' ) '  R8MAT_L_PRINT prints a lower triangular matrix'
  write ( *, '(a)' ) '  stored compactly.  Only the (possibly) nonzero '
  write ( *, '(a)' ) '  elements are printed.'

  call r8mat_l_print ( m1, n1, a1, '  A 7 by 7 matrix.' )

  call r8mat_l_print ( m2, n2, a2, '  A 7 by 3 matrix.' )

  call r8mat_l_print ( m3, n3, a3, '  A 4 by 7 matrix.' )

  return
end
subroutine test070 ( )

!*****************************************************************************80
!
!! TEST070 tests R8MAT_L1_INVERSE.
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
  real ( kind = 8 ), dimension ( n, n ) :: a = reshape ( (/ &
     1.0D+00, 2.0D+00, 0.0D+00, 5.0D+00, 0.0D+00, 75.0D+00, &
     0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,  0.0D+00, &
     0.0D+00, 0.0D+00, 1.0D+00, 3.0D+00, 0.0D+00,  0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00,  6.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00,  4.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,  1.0D+00 /), (/ n, n /) )
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST070'
  write ( *, '(a)' ) '  R8MAT_L1_INVERSE inverts a unit lower triangular matrix.'

  call r8mat_print ( n, n, a, '  Matrix A to be inverted:' )

  call r8mat_l1_inverse ( n, a, b )

  call r8mat_print ( n, n, b, '  Inverse matrix B:' )

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call r8mat_print ( n, n, c, '  Product C = A * B:' )

  return
end
subroutine test071 ( )

!*****************************************************************************80
!
!! TEST071 tests R8MAT_LU.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) l(m,m)
  real ( kind = 8 ) p(m,m)
  real ( kind = 8 ) plu(m,n)
  real ( kind = 8 ) u(m,n)
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    1.0D+00, 10.0D+00, 4.0D+00, 2.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST071'
  write ( *, '(a)' ) '  R8MAT_LU computes the LU factors of a matrix.'

  call r8mat_vand2 ( n, x, a )

  call r8mat_print ( m, n, a, '  Matrix to be factored:' )

  call r8mat_lu ( m, n, a, l, p, u )

  call r8mat_print ( m, m, p, '  P factor:' )

  call r8mat_print ( m, m, l, '  L factor:' )

  call r8mat_print ( m, n, u, '  U factor:' )

  plu(1:m,1:n) = matmul ( transpose ( p(1:m,1:m) ), &
    matmul ( l(1:m,1:m), u(1:m,1:n) ) )

  call r8mat_print ( m, n, plu, '  P*L*U:' )

  return
end
subroutine test072 ( )

!*****************************************************************************80
!
!! TEST072 tests R8MAT_MAX and R8MAT_MIN.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ), parameter :: b = 0.0D+00
  real ( kind = 8 ), parameter :: c = 10.0D+00
  real ( kind = 8 ) r8mat_max
  real ( kind = 8 ) r8mat_min
  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST072'
  write ( *, '(a)' ) '  For a real matrix,'
  write ( *, '(a)' ) '  R8MAT_MAX computes the maximum value;'
  write ( *, '(a)' ) '  R8MAT_MIN computes the minimum value;'

  seed = 123456789

  call r8mat_uniform_ab ( m, n, b, c, seed, a )

  call r8mat_print ( m, n, a, '  Random array:' )

  temp1 = r8mat_min ( m, n, a )
  temp2 = r8mat_max ( m, n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Minimum value = ', temp1
  write ( *, '(a,g14.6)' ) '  Maximum value = ', temp2

  return
end
subroutine test073 ( )

!*****************************************************************************80
!
!! TEST073 tests R8MAT_MAXCOL_MINROW, R8MAT_MAXROW_MINCOL, R8MAT_MINCOL_MAXROW and R8MAT_MINROW_MAXCOL.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ), parameter :: b = 0.0D+00
  real ( kind = 8 ), parameter :: c = 10.0D+00
  real ( kind = 8 ) r8mat_maxcol_minrow
  real ( kind = 8 ) r8mat_maxrow_mincol
  real ( kind = 8 ) r8mat_mincol_maxrow
  real ( kind = 8 ) r8mat_minrow_maxcol
  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST073'
  write ( *, '(a)' ) '  R8MAT_MAXCOL_MINROW computes the maximum over'
  write ( *, '(a)' ) '  columns of the mininum over rows;'
  write ( *, '(a)' ) '  R8MAT_MAXROW_MINCOL computes the maximum over'
  write ( *, '(a)' ) '  rows of the mininum over columns;'
  write ( *, '(a)' ) '  R8MAT_MINCOL_MAXROW computes the minimum over'
  write ( *, '(a)' ) '  columns of the maxinum over rows;'
  write ( *, '(a)' ) '  R8MAT_MINROW_MAXCOL computes the minimum over'
  write ( *, '(a)' ) '  rows of the maxinum over columns;'
  write ( *, '(a)' ) ' '

  seed = 123456789

  call r8mat_uniform_ab ( m, n, b, c, seed, a )

  call r8mat_print ( m, n, a, '  Random array:' )

  temp1 = r8mat_maxcol_minrow ( m, n, a )
  temp2 = r8mat_minrow_maxcol ( m, n, a )

  write ( *, '(a,2g14.6)' ) '  MAXCOL_MINROW, MINROW_MAXCOL = ', temp1, temp2

  temp1 = r8mat_maxrow_mincol ( m, n, a )
  temp2 = r8mat_mincol_maxrow ( m, n, a )

  write ( *, '(a,2g14.6)' ) '  MAXROW_MINCOL, MINCOL_MAXROW = ', temp1, temp2

  return
end
subroutine test0737 ( )

!*****************************************************************************80
!
!! TEST0737 tests R8MAT_NULLSPACE_SIZE and R8MAT_NULLSPACE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 7

  real ( kind = 8 ), dimension ( m, n ) :: a = reshape ( (/ &
    1.0, -2.0, 3.0, -1.0, &
    3.0, -6.0, 9.0, -3.0, &
    0.0,  0.0, 0.0,  0.0, &
    2.0, -2.0, 0.0,  1.0, &
    6.0, -8.0, 6.0,  0.0, &
    3.0,  3.0, 6.0,  9.0, &
    1.0,  1.0, 2.0,  3.0 /), (/ m, n /) )
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: ax
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: nullspace
  integer ( kind = 4 ) nullspace_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0737'
  write ( *, '(a)' ) &
    '  R8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.'
  write ( *, '(a)' ) '  R8MAT_NULLSPACE computes the nullspace of a matrix.'

  call r8mat_print ( m, n, a, '  Input A:' )

  call r8mat_nullspace_size ( m, n, a, nullspace_size )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Nullspace size is ', nullspace_size

  allocate ( nullspace(1:n,1:nullspace_size) )

  call r8mat_nullspace ( m, n, a, nullspace_size, nullspace )

  call r8mat_print ( n, nullspace_size, nullspace, '  Nullspace vectors:' )

  allocate ( ax(1:m,1:nullspace_size) )

  ax = matmul ( a, nullspace )

  call r8mat_print ( m, nullspace_size, ax, &
    '  Product A * Nullspace vectors:' )

  deallocate ( ax )
  deallocate ( nullspace )

  return
end
subroutine test074 ( )

!*****************************************************************************80
!
!! TEST074 tests R8MAT_ORTH_UNIFORM.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) ata(n,n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST074'
  write ( *, '(a)' ) '  R8MAT_ORTH_UNIFORM computes a random orthogonal matrix.'

  seed = 123456789

  call r8mat_orth_uniform ( n, seed, a )

  call r8mat_print ( n, n, a, '  Random orthogonal matrix A' )

  ata(1:n,1:n) = matmul ( transpose ( a(1:n,1:n) ), a(1:n,1:n) )

  call r8mat_print ( n, n, ata, '  AT*A' )

  return
end
subroutine test075 ( )

!*****************************************************************************80
!
!! TEST075 tests R8MAT_PLOT.
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

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 100

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) ip1

  a(1:m,1:n) = 0.0D+00

  do i = 1, m

    a(i,i) = -2.0D+00

    if ( i+1 <= n ) then
      ip1 = i+1
    else
      ip1 = 1
    end if

    a(i,ip1) = 1.0D+00

    if ( 1 <= i-1 ) then
      im1 = i-1
    else
      im1 = n
    end if

    a(i,im1) = 1.0D+00

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST075'
  write ( *, '(a)' ) '  R8MAT_PLOT prints a symbolic picture of a matrix.'
  write ( *, '(a)' ) '  Typically, '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    "-" for negative, '
  write ( *, '(a)' ) '    " " for zero, and'
  write ( *, '(a)' ) '    "+" for positive entries'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  or'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    "X" for nonzero and, '
  write ( *, '(a)' ) '    " " for zero.'
  write ( *, '(a)' ) ' '

  call r8mat_plot ( m, n, a, '  A plot of the matrix:' )

  return
end
subroutine test076 ( )

!*****************************************************************************80
!
!! TEST076 tests R8MAT_POWER_METHOD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) av(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST076'
  write ( *, '(a)' ) '  R8MAT_POWER_METHOD applies the power method'
  write ( *, '(a)' ) '  to a matrix.'

  do i = 1, n
    do j = 1, n
      if ( j == i - 1 .or. j == i + 1 ) then
        a(i,j) = -1.0D+00
      else if ( j == i ) then
        a(i,j) = 2.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  v(1:n) = 0.0D+00

  call r8mat_power_method ( n, a, r, v )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated eigenvalue = ', r

  call r8vec_print ( n, v, '  Estimated eigenvector V:' )

  av(1:n) = matmul ( a(1:n,1:n), v(1:n) )

  call r8vec_print ( n, av, '  Value of A*V:' )

  return
end
subroutine test0764 ( )

!*****************************************************************************80
!
!! TEST076 tests R8MAT_REF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 7

  real ( kind = 8 ), dimension ( m, n ) :: a = reshape ( (/ &
    1.0, -2.0, 3.0, -1.0, &
    3.0, -6.0, 9.0, -3.0, &
    0.0,  0.0, 0.0,  0.0, &
    2.0, -2.0, 0.0,  1.0, &
    6.0, -8.0, 6.0,  0.0, &
    3.0,  3.0, 6.0,  9.0, &
    1.0,  1.0, 2.0,  3.0 /), (/ m, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0764'
  write ( *, '(a)' ) '  R8MAT_REF computes the row echelon form of a matrix.'

  call r8mat_print ( m, n, a, '  Input A:' )

  call r8mat_ref ( m, n, a )

  call r8mat_print ( m, n, a, '  REF form:' )

  return
end
subroutine test0766 ( )

!*****************************************************************************80
!
!! TEST0766 tests R8MAT_RREF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 7

  real ( kind = 8 ), dimension ( m, n ) :: a = reshape ( (/ &
    1.0, -2.0, 3.0, -1.0, &
    3.0, -6.0, 9.0, -3.0, &
    0.0,  0.0, 0.0,  0.0, &
    2.0, -2.0, 0.0,  1.0, &
    6.0, -8.0, 6.0,  0.0, &
    3.0,  3.0, 6.0,  9.0, &
    1.0,  1.0, 2.0,  3.0 /), (/ m, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0766'
  write ( *, '(a)' ) '  R8MAT_RREF computes the reduced row echelon form of a matrix.'

  call r8mat_print ( m, n, a, '  Input A:' )

  call r8mat_rref ( m, n, a )

  call r8mat_print ( m, n, a, '  RREF form:' )

  return
end
subroutine test077 ( )

!*****************************************************************************80
!
!! TEST077 tests R8MAT_SOLVE.
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

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: rhs_num = 2
!
!  Each row of this definition is a COLUMN of the matrix.
!
  real ( kind = 8 ), dimension (n,n+rhs_num) :: a = reshape ( &
    (/ 1.0D+00,  4.0D+00,  7.0D+00, &
       2.0D+00,  5.0D+00,  8.0D+00, &
       3.0D+00,  6.0D+00,  0.0D+00, &
      14.0D+00, 32.0D+00, 23.0D+00, &
       7.0D+00, 16.0D+00,  7.0D+00 /), &
    (/ n, n+rhs_num /) )
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST077'
  write ( *, '(a)' ) '  R8MAT_SOLVE solves linear systems.'
!
!  Print out the matrix to be inverted.
!
  call r8mat_print ( n, n+rhs_num, a, '  The linear system:' )
!
!  Solve the systems.
!
  call r8mat_solve ( n, rhs_num, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input matrix was singular.'
    write ( *, '(a)' ) '  The solutions could not be computed.'
    write ( *, '(a)' ) ' '
    return
  end if

  call r8mat_print ( n, rhs_num, a(1:n,n+1:n+rhs_num), &
    '  The computed solutions' )

  return
end
subroutine test0775 ( )

!*****************************************************************************80
!
!! TEST0775 tests R8MAT_SOLVE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2005
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ), dimension (n,n) :: a
  real ( kind = 8 ), dimension ( n ) :: b
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  real ( kind = 8 ), dimension ( n ) :: x
  real ( kind = 8 ), dimension ( n ) :: x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0775'
  write ( *, '(a)' ) '  R8MAT_SOLVE_2D solves 2D linear systems.'

  seed = 123456789

  do test = 1, test_num

    call r8mat_uniform_01 ( n, n, seed, a )
    call r8vec_uniform_01 ( n, seed, x )
    b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

    call r8mat_solve_2d ( a, b, det, x2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Solution / Computed:'
    write ( *, '(a)' ) ' '

    do i = 1, n
      write ( *, '(2x,g14.6,2x,g14.6)' ) x(i), x2(i)
    end do

  end do

  return
end
subroutine test0776 ( )

!*****************************************************************************80
!
!! TEST0776 tests R8MAT_SOLVE_3D.
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
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension (n,n) :: a
  real ( kind = 8 ), dimension ( n ) :: b
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  real ( kind = 8 ), dimension ( n ) :: x
  real ( kind = 8 ), dimension ( n ) :: x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0776'
  write ( *, '(a)' ) '  R8MAT_SOLVE_3D solves 3D linear systems.'

  seed = 123456789

  do test = 1, test_num

    call r8mat_uniform_01 ( n, n, seed, a )
    call r8vec_uniform_01 ( n, seed, x )
    b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

    call r8mat_solve_3d ( a, b, det, x2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Solution / Computed:'
    write ( *, '(a)' ) ' '

    do i = 1, n
      write ( *, '(2x,g14.6,2x,g14.6)' ) x(i), x2(i)
    end do

  end do

  return
end
subroutine test078 ( )

!*****************************************************************************80
!
!! TEST078 tests R8MAT_SOLVE2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ), dimension ( 2, 2 ) :: a1 = reshape ( (/ &
    1.0D+00, 3.0D+00, &
    2.0D+00, 4.0D+00 /), (/ 2, 2 /) )
  real ( kind = 8 ), dimension ( 3, 3 ) :: a2 = reshape ( (/ &
    2.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 1.0D+00 /), (/ 3, 3 /) )
  real ( kind = 8 ), dimension ( 4, 4 ) :: a3 = reshape ( (/ &
    1.0D+00, 2.0D+00, 1.0D+00, 3.0D+00, &
    0.0D+00, 1.0D+00, 2.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, 3.0D+00, 2.0D+00, &
    1.0D+00, 3.0D+00, 0.0D+00, 1.0D+00 /), (/ 4, 4 /) )
  real ( kind = 8 ), dimension ( 3, 3 ) :: a4 = reshape ( (/ &
    2.0D+00, 1.0D+00, 3.0D+00, &
    4.0D+00, 2.0D+00, 6.0D+00, &
    1.0D+00, 4.0D+00, 5.0D+00 /), (/ 3, 3 /) )
  real ( kind = 8 ), allocatable, dimension ( : ) :: b
  real ( kind = 8 ), dimension ( 2 ) :: b1 = (/ &
    5.0D+00, 11.0D+00 /)
  real ( kind = 8 ), dimension ( 3 ) :: b2 = (/ &
    4.0D+00, 2.0D+00, 2.0D+00 /)
  real ( kind = 8 ), dimension ( 4 ) :: b3 = (/ &
    5.0D+00, 11.0D+00, 16.0D+00, 15.0D+00 /)
  real ( kind = 8 ), dimension ( 3 ) :: b4 = (/ &
    13.0D+00, 17.0D+00, 20.0D+00 /)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ), dimension ( test_num ) :: n_test = (/ 2, 3, 4, 3 /)
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST078'
  write ( *, '(a)' ) '  R8MAT_SOLVE2 is a linear solver.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    n = n_test ( test )

    allocate ( a(1:n,1:n) )
    allocate ( b(1:n) )
    allocate ( x(1:n) )

    if ( test == 1 ) then
      a(1:n,1:n) = a1(1:n,1:n)
      b(1:n) = b1(1:n)
    else if ( test == 2 ) then
      a(1:n,1:n) = a2(1:n,1:n)
      b(1:n) = b2(1:n)
    else if ( test == 3 ) then
      a(1:n,1:n) = a3(1:n,1:n)
      b(1:n) = b3(1:n)
    else if ( test == 4 ) then
      a(1:n,1:n) = a4(1:n,1:n)
      b(1:n) = b4(1:n)
    end if

    call r8vec_print ( n, b, '  Right hand side:' )

    call r8mat_solve2 ( n, a, b, x, ierror )

    write ( *, '(a)' ) ' '
    if ( ierror == 0 ) then
      write ( *, '(a)' ) '  The system is nonsingular.'
    else if ( ierror == 1 ) then
      write ( *, '(a)' ) '  The system is singular, but consistent.'
    else if ( ierror == 2 ) then
      write ( *, '(a)' ) '  The system is singular and inconsistent.'
    end if

    call r8vec_print ( n, x, '  Computed solution:' )

    deallocate ( a )
    deallocate ( b )
    deallocate ( x )

  end do

  return
end
subroutine test079 ( )

!*****************************************************************************80
!
!! TEST079 tests R8MAT_SYMM_JACOBI;
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) q(n,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST079'
  write ( *, '(a)' ) '  For a symmetric R8MAT:'
  write ( *, '(a)' ) '  R8MAT_SYMM_JACOBI diagonalizes;'
!
!  Choose the eigenvalues.
!
  call r8vec_indicator ( n, x )
!
!  Choose the eigenvectors.
!
  seed = 123456789
  call r8mat_orth_uniform ( n, seed, q )
!
!  Compute A = Q*X*Q.
!
  call r8mat_symm_eigen ( n, x, q, a )

  call r8mat_print ( n, n, a, '  Matrix to diagonalize:' )

  call r8mat_symm_jacobi ( n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed Eigenvalues:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) a(i,i)
  end do

  return
end
subroutine test080 ( )

!*****************************************************************************80
!
!! TEST080 tests R8MAT_TO_R8PLU and R8PLU_TO_R8MAT;
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) :: b = 0.0D+00
  real ( kind = 8 ) :: c = 1.0D+00
  integer ( kind = 4 ) info
  real ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST080'
  write ( *, '(a)' ) '  R8MAT_TO_R8PLU determines the compressed PLU factors'
  write ( *, '(a)' ) '  of a R8MAT.'
  write ( *, '(a)' ) '  R8PLU_TO_R8MAT determines the original matrix from'
  write ( *, '(a)' ) '  the compressed PLU factors.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8mat_uniform_ab ( n, n, b, c, seed, a )

  call r8mat_print ( n, n, a, '  The matrix A:' )
!
!  Factor the matrix.
!
  call r8mat_to_r8plu ( n, a, pivot, lu, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Warning!'
    write ( *, '(a)' ) '  R8MAT_TO_R8PLU declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
  end if
!
!  Display the gory details.
!
  call i4vec_print ( n, pivot, '  The pivot vector P:' )

  call r8mat_print ( n, n, lu, '  The compressed LU factors:' )
!
!  Recover the matrix from the PLU factors.
!
  call r8plu_to_r8mat ( n, pivot, lu, a2 )

  call r8mat_print ( n, n, a2, '  The recovered matrix A2:' )

  return
end
subroutine test081 ( )

!*****************************************************************************80
!
!! TEST081 tests R8MAT_TRACE.
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

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_trace
  real ( kind = 8 ) trace

  do i = 1, n
    do j = 1, n

      if ( i <= j ) then
        a(i,j) = real ( n + 1 - j, kind = 8 )
      else if ( j == i - 1 ) then
        a(i,j) = real ( n - j, kind = 8 )
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST081'
  write ( *, '(a)' ) '  R8MAT_TRACE computes the trace of a matrix'

  call r8mat_print ( n, n, a, '  Matrix:' )

  trace = r8mat_trace ( n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Trace is ', trace

  return
end
subroutine test082 ( )

!*****************************************************************************80
!
!! TEST082 tests R8MAT_TRANSPOSE_PRINT;
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

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ), parameter :: n = 12

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST082'
  write ( *, '(a)' ) '  R8MAT_TRANSPOSE_PRINT prints a R8MAT,'
  write ( *, '(a)' ) '  transposed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix row order M =    ', m
  write ( *, '(a,i8)' ) '  Matrix column order N = ', n
!
!  Set the matrix.
!
  do i = 1, m
    do j = 1, n
      a(i,j) = real ( i * 100 + j, kind = 8 )
    end do
  end do

  call r8mat_transpose_print ( m, n, a, '  The transposed matrix A:' )

  return
end
subroutine test083 ( )

!*****************************************************************************80
!
!! TEST083 tests R8MAT_U_INVERSE.
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

  integer ( kind = 4 ), parameter :: n = 4
!
!  Each row of this definition is a COLUMN of the matrix.
!
  real ( kind = 8 ), dimension ( n, n ) :: a = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00,  0.0D+00, &
    2.0D+00, 3.0D+00, 0.0D+00,  0.0D+00, &
    4.0D+00, 5.0D+00, 6.0D+00,  0.0D+00, &
    7.0D+00, 8.0D+00, 9.0D+00, 10.0D+00 /), (/ n, n /) )
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST083'
  write ( *, '(a)' ) '  R8MAT_U_INVERSE inverts an upper triangular matrix.'

  call r8mat_print ( n, n, a, '  Input matrix A' )

  call r8mat_u_inverse ( n, a, b )

  call r8mat_print ( n, n, b, '  Inverse matrix B:' )

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call r8mat_print ( n, n, c, '  Product C = A * B:' )

  return
end
subroutine test084 ( )

!*****************************************************************************80
!
!! TEST084 tests R8MAT_U1_INVERSE.
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
  real ( kind = 8 ), dimension ( n, n ) :: a = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,  0.0D+00, &
    2.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,  0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00,  0.0D+00, &
    5.0D+00, 0.0D+00, 3.0D+00, 1.0D+00, 0.0D+00,  0.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00,  0.0D+00, &
   75.0D+00, 0.0D+00, 0.0D+00, 6.0D+00, 4.0D+00,  1.0D+00 /), (/ n, n /) )
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST084'
  write ( *, '(a)' ) '  R8MAT_U1_INVERSE inverts a unit upper triangular matrix.'

  call r8mat_print ( n, n, a, '  Input matrix A' )

  call r8mat_u1_inverse (  n, a, b )

  call r8mat_print ( n, n, b, '  Inverse matrix B:' )

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call r8mat_print ( n, n, c, '  Product C = A * B:' )

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests R8MAT_UNIFORM.
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

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ), dimension (m,n) :: a
  real ( kind = 8 ), parameter :: b = 2.0D+00
  real ( kind = 8 ), parameter :: c = 10.0D+00
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  R8MAT_UNIFORM sets a matrix to random values.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r8mat_uniform_ab ( m, n, b, c, seed, a )
!
!  Print out the matrix to be inverted.
!
  call r8mat_print ( m, n, a, '  The random matrix:' )

  return
end
subroutine test086 ( )

!*****************************************************************************80
!
!! TEST086 tests R8PLU_DET;
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) info
  real ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST086'
  write ( *, '(a)' ) '  R8PLU_DET determines the determinant of a matrix'
  write ( *, '(a)' ) '  from its compressed PLU factors.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  call r8mat_print ( n, n, a, '  The matrix A:' )
!
!  Factor the matrix.
!
  call r8mat_to_r8plu ( n, a, pivot, lu, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  R8MAT_TO_R8PLU declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Compute the determinant.
!
  call r8plu_det ( n, pivot, lu, det )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The determinant = ', det

  return
end
subroutine test087 ( )

!*****************************************************************************80
!
!! TEST087 tests R8PLU_INVERSE;
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) info
  real ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST087'
  write ( *, '(a)' ) '  R8PLU_INVERSE determines the inverse of a matrix'
  write ( *, '(a)' ) '  from its compressed PLU factors.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  call r8mat_print ( n, n, a, '  The matrix A:' )
!
!  Factor the matrix.
!
  call r8mat_to_r8plu ( n, a, pivot, lu, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST0852 - Fatal error!'
    write ( *, '(a)' ) '  R8MAT_TO_R8PLU declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Compute the inverse.
!
  call r8plu_inverse ( n, pivot, lu, b )

  call r8mat_print ( n, n, b, '  The inverse matrix B:' )
!
!  Compute the product C = A * B.
!
  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  call r8mat_print ( n, n, c, '  The product matrix C = A * B:' )

  return
end
subroutine test088 ( )

!*****************************************************************************80
!
!! TEST088 tests R8PLU_MUL;
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST088'
  write ( *, '(a)' ) '  R8PLU_MUL computes the product A*x=b'
  write ( *, '(a)' ) '  using the compressed PLU factors of A.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  call r8mat_print ( n, n, a, '  The matrix A:' )
!
!  Set the right hand side B1.
!
  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

  call r8vec_print ( n, b, '  The right hand side B (computed from A):' )
!
!  Factor the matrix.
!
  call r8mat_to_r8plu ( n, a, pivot, lu, info )
!
!  Compute the matrix-vector product.
!
  call r8plu_mul ( n, pivot, lu, x, b )

  call r8vec_print ( n, b, '  The right hand side B (computed from PLU):' )

  return
end
subroutine test089 ( )

!*****************************************************************************80
!
!! TEST089 tests R8PLU_SOL;
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST089'
  write ( *, '(a)' ) '  R8PLU_SOL solves a linear system A*x=b'
  write ( *, '(a)' ) '  using the compressed PLU factors of A.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8mat_uniform_01 ( n, n, seed, a )

  call r8mat_print ( n, n, a, '  The matrix A:' )
!
!  Set the right hand side.
!
  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do
  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )
  x(1:n) = 0.0D+00

  call r8vec_print ( n, b, '  The right hand side B:' )
!
!  Factor the matrix.
!
  call r8mat_to_r8plu ( n, a, pivot, lu, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  R8MAT_TO_R8PLU declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call r8plu_sol ( n, pivot, lu, b, x )

  call r8vec_print ( n, x, '  The computed solution X:' )

  return
end
subroutine test090 ( )

!*****************************************************************************80
!
!! TEST090 tests R8POLY_DERIV.
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

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) cp(0:n)
  integer ( kind = 4 ) d
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST090'
  write ( *, '(a)' ) '  R8POLY_DERIV computes the coefficients of'
  write ( *, '(a)' ) '  the derivative of a polynomial.'

  call r8vec_indicator ( n, x )

  call roots_to_r8poly ( n, x, c )

  call r8poly_print ( n, c, '  The initial polynomial' )

  do d = 0, n
    call r8poly_deriv ( n, c, d, cp )
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The derivative of order ', d
    write ( *, '(a)' ) ' '
    call r8poly_print ( n-d, cp, ' ' )
  end do

  return
end
subroutine test091 ( )

!*****************************************************************************80
!
!! TEST091 tests R8POLY_LAGRANGE_COEF.
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

  integer ( kind = 4 ), parameter :: npol = 5

  integer ( kind = 4 ) ipol
  real ( kind = 8 ) pcof(0:npol-1)
  real ( kind = 8 ) xpol(npol)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST091'
  write ( *, '(a)' ) '  R8POLY_LAGRANGE_COEF returns the coefficients'
  write ( *, '(a)' ) '  for a Lagrange basis polynomial.'
  write ( *, '(a)' ) '  R8POLY_PRINT prints a polynomial.'

  call r8vec_indicator ( npol, xpol )

  call r8vec_print ( npol, xpol, '  Abscissas:' )

  do ipol = 1, npol

    call r8poly_lagrange_coef ( npol, ipol, xpol, pcof )

    call r8poly_print ( npol-1, pcof, '  The Lagrange basis polynomial:' )

  end do

  return
end
subroutine test092 ( )

!*****************************************************************************80
!
!! TEST092 tests R8POLY_LAGRANGE_COEF and R8POLY_DERIV.
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

  integer ( kind = 4 ), parameter :: npol = 5

  integer ( kind = 4 ) d
  integer ( kind = 4 ) ipol
  real ( kind = 8 ) pcof(0:npol-1)
  real ( kind = 8 ) pprime(0:npol-1)
  real ( kind = 8 ) xpol(npol)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST092'
  write ( *, '(a)' ) '  R8POLY_LAGRANGE_COEF returns the coefficients'
  write ( *, '(a)' ) '  for a Lagrange basis polynomial.'
  write ( *, '(a)' ) '  R8POLY_DERIV computes derivatives of a polynomial.'

  call r8vec_indicator ( npol, xpol )

  call r8vec_print ( npol, xpol, '  Abscissas:' )

  do ipol = 1, npol

    call r8poly_lagrange_coef ( npol, ipol, xpol, pcof )

    call r8poly_print ( npol-1, pcof, '  The Lagrange basis polynomial:' )

    do d = 1, npol-1
      call r8poly_deriv ( npol-1, pcof, d, pprime )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The derivative of order ', d
      call r8poly_print ( npol-1-d, pprime, ' ' )
    end do

  end do

  return
end
subroutine test093 ( )

!*****************************************************************************80
!
!! TEST093 tests R8POLY_LAGRANGE_0, R8POLY_LAGRANGE_1, and R8POLY_LAGRANGE_2.
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

  integer ( kind = 4 ), parameter :: npol = 5

  real ( kind = 8 ) dw2dx2
  real ( kind = 8 ) dwdx
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) nx
  real ( kind = 8 ) wval
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST093'
  write ( *, '(a)' ) '  R8POLY_LAGRANGE_0 evaluates the Lagrange'
  write ( *, '(a)' ) '  factor W(X) at a point.'
  write ( *, '(a)' ) '  R8POLY_LAGRANGE_1 evaluates the Lagrange'
  write ( *, '(a)' ) '  factor W''(X) at a point.'
  write ( *, '(a)' ) '  R8POLY_LAGRANGE_2 evaluates the Lagrange'
  write ( *, '(a)' ) '  factor W"(X) at a point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data points is ', npol
!
!  Set the abscissas of the polynomials.
!
  xlo = 0.0D+00
  xhi = real ( npol - 1, kind = 8 )

  call r8vec_even ( npol, xlo, xhi, xpol )

  call r8vec_print ( npol, xpol, '  Abscissas:' )
!
!  Evaluate W(X), W'(X), W''.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          W(X)          W''(X)        W"(X)'
  write ( *, '(a)' ) ' '

  nx = 4 * npol - 1

  do ival = 1, nx

    call r8vec_even_select ( nx, xlo, xhi, ival, xval )

    call r8poly_lagrange_0 ( npol, xpol, xval, wval )
    call r8poly_lagrange_1 ( npol, xpol, xval, dwdx )
    call r8poly_lagrange_2 ( npol, xpol, xval, dw2dx2 )

    write ( *, '(6g12.4)' ) xval, wval, dwdx, dw2dx2

  end do

  return
end
subroutine test094 ( )

!*****************************************************************************80
!
!! TEST094 tests R8POLY_LAGRANGE_FACTOR.
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

  integer ( kind = 4 ), parameter :: npol = 5

  real ( kind = 8 ) dwdx
  integer ( kind = 4 ) i
  real ( kind = 8 ) wval
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST094'
  write ( *, '(a)' ) '  R8POLY_LAGRANGE_FACTOR evaluates the Lagrange'
  write ( *, '(a)' ) '  factor W(X) at a point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  For this test, we use ', npol, ' functions.'
!
!  Set the abscissas of the polynomials.
!
  xlo = 0.0D+00
  xhi = real ( npol - 1, kind = 8 )

  do i = 1, npol
    xpol(i) = ( real ( npol - ( i - 1 ), kind = 8 ) * xlo &
              + real (        ( i - 1 ), kind = 8 ) * xhi ) &
              / real ( npol,             kind = 8 )
  end do

  call r8vec_print ( npol, xpol, '  Abscissas:' )
!
!  Evaluate W(X) and W'(X).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          W(X)          W''(X)'
  write ( *, '(a)' ) ' '

  do i = 0, 2 * npol - 1

    call r8vec_even_select ( 2 * npol - 1, xhi, xlo, i, xval )

    call r8poly_lagrange_factor ( npol, xpol, xval, wval, dwdx )

    write ( *, '(2x,f10.4,2x,f10.4,2x,f10.4)' ) xval, wval, dwdx

  end do

  return
end
subroutine test095 ( )

!*****************************************************************************80
!
!! TEST095 tests R8POLY_LAGRANGE_VAL.
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

  integer ( kind = 4 ), parameter :: npol = 5

  real ( kind = 8 ) dpdx(npol)
  integer ( kind = 4 ) ipol
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) nx
  real ( kind = 8 ) pval(npol)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST095'
  write ( *, '(a)' ) '  R8POLY_LAGRANGE_VAL evaluates a Lagrange'
  write ( *, '(a)' ) '  interpolating polynomial at a point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of data points = ', npol
!
!  Set the abscissas of the polynomials.
!
  xlo = 0.0D+00
  xhi = real ( npol - 1, kind = 8 )
  call r8vec_even ( npol, xlo, xhi, xpol )

  call r8vec_print ( npol, xpol, '  Abscissas:' )
!
!  Evaluate the polynomials.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here are the values of the functions at '
  write ( *, '(a)' ) '  several points:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          L1          L2          L3      L4' // &
    '          L5'
  write ( *, '(a)' ) ' '

  nx = 2 * npol - 1

  do ival = 1, nx

    call r8vec_even_select ( nx, xlo, xhi, ival, xval )

    do ipol = 1, npol
      call r8poly_lagrange_val ( npol, ipol, xpol, xval, pval(ipol), dpdx(ipol) )
    end do

    write ( *, '(6g12.4)' ) xval, ( pval(ipol), ipol = 1, npol )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  And the derivatives:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          L''1         L''2         L''3' // &
    '     L''4         L''5'
  write ( *, '(a)' ) ' '

  nx = 2 * npol - 1

  do ival = 1, nx

    call r8vec_even_select ( nx, xlo, xhi, ival, xval )

    do ipol = 1, npol
      call r8poly_lagrange_val ( npol, ipol, xpol, xval, pval(ipol), dpdx(ipol) )
    end do

    write ( *, '(6g12.4)' ) xval, ( dpdx(ipol), ipol = 1, npol )

  end do

  return
end
subroutine test098 ( )

!*****************************************************************************80
!
!! TEST098 tests R8POLY_VALUE_HORNER.
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

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension (0:n) :: p = (/ &
    24.0D+00, -50.0D+00, +35.0D+00, -10.0D+00, 1.0D+00 /)
  real ( kind = 8 ) pval
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST098'
  write ( *, '(a)' ) '  R8POLY_VALUE_HORNER evaluates a polynomial at a'
  write ( *, '(a)' ) '  point, using Horner''s method.'

  call r8poly_print ( n, p, '  The polynomial:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        X            P(X)'
  write ( *, '(a)' ) ' '

  do i = 0, 15
    x = real ( i, kind = 8 ) / 3.0D+00
    call r8poly_value_horner ( n, p, x, pval )
    write ( *, '(2x,g14.6,g14.6)' ) x, pval
  end do

  return
end
subroutine test099 ( )

!*****************************************************************************80
!
!! TEST099 tests R8POLY2_EX and R8POLY2_EX2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST099'
  write ( *, '(a)' ) '  R8POLY2_EX finds the extreme value'
  write ( *, '(a)' ) '  of a parabola determined by three points.'
  write ( *, '(a)' ) '  R8POLY2_EX2 finds the extreme value'
  write ( *, '(a)' ) '  of a parabola determined by three points.'
  a =  2.0D+00
  b = -4.0D+00
  c = 10.0D+00

  x1 = 1.0D+00
  y1 = a * x1**2 + b * x1 + c
  x2 = 2.0D+00
  y2 = a * x2**2 + b * x2 + c
  x3 = 3.0D+00
  y3 = a * x3**2 + b * x3 + c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Parabolic coefficients A, B, C ='
  write ( *, '(2x,3g14.6)' ) a, b, c
  write ( *, '(a)' ) ' '

  call r8r8_print ( x1, y1, '  Point 1' )
  call r8r8_print ( x2, y2, '  Point 2' )
  call r8r8_print ( x3, y3, '  Point 3' )

  a = 0.0D+00
  b = 0.0D+00
  c = 0.0D+00

  call r8poly2_ex ( x1, y1, x2, y2, x3, y3, xmin, ymin, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  R8POLY2_EX returns XMIN, YMIN = ', xmin, ymin

  call r8poly2_ex2 ( x1, y1, x2, y2, x3, y3, xmin, ymin, a, b, c, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  R8POLY2_EX2 returns XMIN, YMIN = ', xmin, ymin
  write ( *, '(a,3g14.6)' ) '  and A, B, C = ', a, b, c

  return
end
subroutine test100 ( )

!*****************************************************************************80
!
!! TEST100 tests R8POLY2_VAL.
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
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yp
  real ( kind = 8 ) ypp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST100'
  write ( *, '(a)' ) '  R8POLY2_VAL evaluates a parabola given'
  write ( *, '(a)' ) '  3 data points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our parabola will be 2*x*x + 3 * x + 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Case 1: 3 distinct data points:'
  write ( *, '(a)' ) ' '

  x1 = -1.0D+00
  x2 = 1.0D+00
  x3 = 3.0D+00

  call test100_f ( x1, y1, yp, ypp )
  call test100_f ( x2, y2, yp, ypp )
  call test100_f ( x3, y3, yp, ypp )

  write ( *, '(2x,2g14.6)' ) x1, y1
  write ( *, '(2x,2g14.6)' ) x2, y2
  write ( *, '(2x,2g14.6)' ) x3, y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sampled data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, Y, Y'', Y"'
  write ( *, '(a)' ) ' '
  do i = 0, 3
    x = real ( i, kind = 8 )
    call r8poly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )
    write ( *, '(2x,4g14.6)' ) x, y, yp, ypp
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Case 2: X1=X2, X3 distinct:'
  write ( *, '(a)' ) ' '

  x1 = - 1.0D+00
  x2 = - 1.0D+00
  x3 = 3.0D+00

  call test100_f ( x1, y1, y2, ypp )
  call test100_f ( x3, y3, yp, ypp )
  write ( *, '(2x,2g14.6)' ) x1, y1
  write ( *, '(2x,2g14.6)' ) x2, y2
  write ( *, '(2x,2g14.6)' ) x3, y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sampled data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, Y, Y'', Y"'
  write ( *, '(a)' ) ' '
  do i = 0, 3
    x = real ( i, kind = 8 )
    call r8poly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )
    write ( *, '(2x,4g14.6)' ) x, y, yp, ypp
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Case 3: X1=X2=X3:'
  write ( *, '(a)' ) ' '

  x1 = - 1.0D+00
  x2 = - 1.0D+00
  x3 = - 1.0D+00

  call test100_f ( x1, y1, y2, y3 )

  write ( *, '(2x,2g14.6)' ) x1, y1
  write ( *, '(2x,2g14.6)' ) x2, y2
  write ( *, '(2x,2g14.6)' ) x3, y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sampled data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, Y, Y'', Y"'
  write ( *, '(a)' ) ' '
  do i = 0, 3
    x = real ( i, kind = 8 )
    call r8poly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )
    write ( *, '(2x,4g14.6)' ) x, y, yp, ypp
  end do

  return
end
subroutine test100_f ( x, y, yp, ypp )

!*****************************************************************************80
!
!! TEST100_F evaluates a parabola for us.
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

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) yp
  real ( kind = 8 ) ypp

  y = 2.0D+00 * x**2 + 3.0D+00 * x + 1.0D+00
  yp = 4.0D+00 * x + 3.0D+00
  ypp = 4.0D+00

  return
end
subroutine test101 ( )

!*****************************************************************************80
!
!! TEST101 tests R8POLY2_VAL2.
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

  integer ( kind = 4 ), parameter :: ndata = 5
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  real ( kind = 8 ) xdata(ndata)
  real ( kind = 8 ) xval
  real ( kind = 8 ) ydata(dim_num,ndata)
  real ( kind = 8 ) yval(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST101'
  write ( *, '(a)' ) '  R8POLY2_VAL2 evaluates parabolas through'
  write ( *, '(a)' ) '  3 points in a table'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our data tables will actually be parabolas:'
  write ( *, '(a)' ) '    A: 2*x*x + 3 * x + 1.'
  write ( *, '(a)' ) '    B: 4*x*x - 2 * x + 5.'
  write ( *, '(a)' ) ' '

  do i = 1, ndata
    xval = 2.0D+00 * real ( i, kind = 8 )
    xdata(i) = xval
    ydata(1,i) = 2.0D+00 * xval**2 + 3.0D+00 * xval + 1.0D+00
    ydata(2,i) = 4.0D+00 * xval**2 - 2.0D+00 * xval + 5.0D+00
    write ( *, '(2x,i8,3g14.6)' ) i, xdata(i), ydata(1,i), ydata(2,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolated data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LEFT, X, Y1, Y2'
  write ( *, '(a)' ) ' '

  do i = 0, 4
    xval = real ( 2 * i + 1, kind = 8 )
    left = max ( min ( i + 1, ndata - 2 ), 1 )
    call r8poly2_val2 ( dim_num, ndata, xdata, ydata, left, xval, yval )
    write ( *, '(2x,i8,3g14.6)' ) left, xval, yval(1), yval(2)
  end do

  return
end
subroutine test102 ( )

!*****************************************************************************80
!
!! TEST102 tests R8POLY2_ROOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ), dimension(test_num) :: a_test = (/ &
    2.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), dimension(test_num) :: b_test = (/ &
    -2.0D+00, -20.0D+00, -2.0D+00 /)
  real ( kind = 8 ) c
  real ( kind = 8 ), dimension(test_num) :: c_test = (/ &
    -24.0D+00, 100.0D+00, 10.0D+00 /)
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST102'
  write ( *, '(a)' ) '  R8POLY2_ROOT finds quadratic equation roots.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A         B         C     R1         R2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a = a_test(test)
    b = b_test(test)
    c = c_test(test)

    call r8poly2_root ( a, b, c, r1, r2 )

    write ( *, '(2x,3f8.1,4g14.6)' ) a, b, c, r1, r2

  end do

  return
end
subroutine test103 ( )

!*****************************************************************************80
!
!! TEST103 tests R8POLY3_ROOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) a
  real ( kind = 8 ), dimension(test_num) :: a_test = (/ &
    1.0D+00, 9.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), dimension(test_num) :: b_test = (/ &
    -6.0D+00, -36.0D+00, -5.0D+00, -8.0D+00  /)
  real ( kind = 8 ) c
  real ( kind = 8 ), dimension(test_num) :: c_test = (/ &
    11.0D+00, 54.0D+00, 8.0D+00, 25.0D+00  /)
  real ( kind = 8 ) d
  real ( kind = 8 ), dimension(test_num) :: r8_test = (/ &
    -6.0D+00, -27.0D+00, -4.0D+00, -26.0D+00  /)
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  integer ( kind = 4 ) test
!
!  1: Three distinct real roots, 1, 2, 3.
!  2: One repeated real root, 1.5, 1.5, 1.5.
!  3: Two real roots, one repeated, 1, 2, 2.
!  4: One real root, a complex conjugate pair, 2, 3+2I, 3-2I.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST103'
  write ( *, '(a)' ) '  R8POLY3_ROOT finds roots of cubic equations.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a = a_test(test)
    b = b_test(test)
    c = c_test(test)
    d = r8_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Polynomial coefficients A, B, C, D:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,4g14.6)' ) a, b, c, d

    call r8poly3_root ( a, b, c, d, r1, r2, r3 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Roots:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,2g14.6)' ) r1
    write ( *, '(2x,2g14.6)' ) r2
    write ( *, '(2x,2g14.6)' ) r3

  end do

  return
end
subroutine test104 ( )

!*****************************************************************************80
!
!! TEST104 tests R8POLY4_ROOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  real ( kind = 8 ) a
  real ( kind = 8 ), dimension(test_num) :: a_test = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), dimension(test_num) :: b_test = (/ &
    -10.0D+00, -5.0D+00, -22.0D+00, -16.0D+00, -20.0D+00, &
    2.0D+00, 0.0D+00 /)
  real ( kind = 8 ) c
  real ( kind = 8 ), dimension(test_num) :: c_test = (/ &
    35.0D+00, 1.0D+00, 141.0D+00, 72.0D+00, 150.0D+00, &
    1.0D+00, 13.0D+00 /)
  real ( kind = 8 ) d
  real ( kind = 8 ), dimension(test_num) :: r8_test = (/ &
    -50.0D+00, 21.0D+00, -220.0D+00, -128.0D+00, -500.0D+00, &
    8.0D+00, 0.0D+00 /)
  real ( kind = 8 ) e
  real ( kind = 8 ), dimension(test_num) :: e_test = (/ &
    24.0D+00, -18.0D+00, +100.0D+00, 80.0D+00, 625.0D+00, &
    -12.0D+00, 36.0D+00 /)
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  complex ( kind = 8 ) r4
  integer ( kind = 4 ) test
!
!  1: Four distinct real roots, 1, 2, 3, 4.
!  2: Three distinct real roots, 1, -2, 3, 3
!  3: Two distinct real roots, 1, 1, 10, 10.
!  4: Two distinct real roots, 2, 2, 2, 10
!  5: One real root, 5, 5, 5, 5
!  6: Two distinct real roots, one complex conjugate pair.
!  7: Two distinct complex conjugate pairs.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST104'
  write ( *, '(a)' ) '  R8POLY4_ROOT finds roots of quartic equations.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a = a_test(test)
    b = b_test(test)
    c = c_test(test)
    d = r8_test(test)
    e = e_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  A =', a
    write ( *, '(a,g14.6)' ) '  B =', b
    write ( *, '(a,g14.6)' ) '  C =', c
    write ( *, '(a,g14.6)' ) '  D =', d
    write ( *, '(a,g14.6)' ) '  E =', e

    call r8poly4_root ( a, b, c, d, e, r1, r2, r3, r4 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Roots:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,2g14.6)' ) r1
    write ( *, '(2x,2g14.6)' ) r2
    write ( *, '(2x,2g14.6)' ) r3
    write ( *, '(2x,2g14.6)' ) r4

  end do

  return
end
subroutine test105 ( )

!*****************************************************************************80
!
!! TEST105 tests R8ROW_MAX and R8ROW_MIN;
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) amax(m)
  real ( kind = 8 ) amin(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105'
  write ( *, '(a)' ) '  For a R8ROW (a matrix regarded as rows):'
  write ( *, '(a)' ) '  R8ROW_MAX computes maximums;'
  write ( *, '(a)' ) '  R8ROW_MIN computes minimums;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k )
    end do
  end do

  call r8mat_print ( m, n, a, '  The original matrix:' )

  call r8row_max ( m, n, a, amax )

  call r8row_min ( m, n, a, amin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Row maximum, minimum:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(i3,3x,4f10.4)' ) i, amax(i), amin(i)
  end do

  return
end
subroutine test106 ( )

!*****************************************************************************80
!
!! TEST106 tests R8ROW_MEAN and R8ROW_SUM;
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) mean(m)
  real ( kind = 8 ) rowsum(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST106'
  write ( *, '(a)' ) '  For a R8ROW (a matrix regarded as rows):'
  write ( *, '(a)' ) '  R8ROW_MEAN computes means;'
  write ( *, '(a)' ) '  R8ROW_SUM computes sums;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k )
    end do
  end do

  call r8mat_print ( m, n, a, '  The original matrix:' )

  call r8row_sum ( m, n, a, rowsum )

  call r8row_mean ( m, n, a, mean )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Row sum, mean:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(i3,3x,4f10.4)' ) i, rowsum(i), mean(i)
  end do

  return
end
subroutine test1064 ( )

!*****************************************************************************80
!
!! TEST1064 tests R8ROW_PART_QUICK_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ), dimension ( m, n ) :: a = reshape ( (/ &
     2.0D+00, 8.0D+00, 6.0D+00, 0.0D+00, 10.0D+00, &
    10.0D+00, 0.0D+00, 5.0D+00, &
     4.0D+00, 8.0D+00, 2.0D+00, 2.0D+00,  6.0D+00, &
     0.0D+00, 6.0D+00, 8.0D+00 /), (/ m, n /) )
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1064'
  write ( *, '(a)' ) '  For an R8ROW;'
  write ( *, '(a)' ) '  R8ROW_PART_QUICK_A partitions the matrix.'

  call r8mat_print ( m, n, a, '  The matrix:' )

  l = 2
  r = 4
  call r8row_part_quick_a ( m, n, a, l, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  L = ', l
  write ( *, '(a,i4)' ) '  R = ', r

  call r8mat_print ( m, n, a, '  The partitioned matrix:' )

  return
end
subroutine test1065 ( )

!*****************************************************************************80
!
!! TEST1065 tests R8ROW_SORT_HEAP_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
     2.0D+00,  4.0D+00, 1.0D+00,  3.0D+00, &
     6.0D+00,  8.0D+00, 5.0D+00,  7.0D+00, &
    10.0D+00, 12.0D+00, 9.0D+00, 11.0D+00 /), (/ m, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1065'
  write ( *, '(a)' ) '  For an R8ROW;'
  write ( *, '(a)' ) '  R8ROW_SORT_HEAP_A does an ascending heap sort'

  call r8mat_print ( m, n, a, '  The unsorted matrix:' )

  call r8row_sort_heap_a ( m, n, a )

  call r8mat_print ( m, n, a, '  The sorted matrix:' )

  return
end
subroutine test1066 ( )

!*****************************************************************************80
!
!! TEST1066 tests R8ROW_SORT_HEAP_INDEX_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 15
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension (m,n) :: a = reshape ( (/ &
    2.0D+00,  4.0D+00,  1.0D+00,  3.0D+00,  2.0D+00, &
    3.0D+00,  0.0D+00,  0.0D+00,  2.0D+00,  3.0D+00, &
    2.0D+00,  2.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
    6.0D+00,  8.0D+00,  5.0D+00,  7.0D+00,  6.0D+00, &
    4.0D+00,  0.0D+00,  6.0D+00,  6.0D+00,  7.0D+00, &
    0.0D+00,  6.0D+00,  5.0D+00,  5.0D+00,  5.1D+00, &
   10.0D+00, 12.0D+00,  9.0D+00, 11.0D+00,  0.0D+00, &
   18.0D+00,  0.0D+00, 10.0D+00, 10.0D+00, 11.0D+00, &
   10.0D+00, 10.0D+00,  9.0D+00,  9.1D+00,  9.0D+00 /), (/ m, n /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) indx(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1066'
  write ( *, '(a)' ) '  R8ROW_SORT_HEAP_INDEX_A computes an index vector which'
  write ( *, '(a)' ) '  ascending sorts an R8ROW.'

  call r8mat_transpose_print ( m, n, a, '  The unsorted R8ROW:' )

  call r8row_sort_heap_index_a ( m, n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The implicitly sorted R8ROW'
  write ( *, '(a)' ) ' '

  do i = 1, m
    i2 = indx(i)
    write ( *, '(2x,i4,a,2x,f10.1,2x,f10.1,2x,f10.1)' ) i2, ':', a(i2,1:n)
  end do

  return
end
subroutine test1067 ( )

!*****************************************************************************80
!
!! TEST1067 tests R8ROW_SORT_QUICK_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ), parameter :: b = 0.0D+00
  real ( kind = 8 ), parameter :: c = 10.0D+00
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1067'
  write ( *, '(a)' ) '  For an R8ROW;'
  write ( *, '(a)' ) '  R8ROW_SORT_QUICK_A does a quicksort.'

  seed = 123456789

  call r8mat_uniform_ab ( m, n, b, c, seed, a )

  call r8mat_print ( m, n, a, '  The unsorted matrix:' )

  call r8row_sort_quick_a ( m, n, a )

  call r8mat_print ( m, n, a, '  The sorted matrix:' )

  return
end
subroutine test107 ( )

!*****************************************************************************80
!
!! TEST107 tests R8ROW_SWAP;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST107'
  write ( *, '(a)' ) '  For a R8ROW (a matrix regarded as rows):'
  write ( *, '(a)' ) '  R8ROW_SWAP swaps two rows;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k )
    end do
  end do

  call r8mat_print ( m, n, a, '  The original matrix:' )

  row1 = 1
  row2 = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,i3,a,i3)' ) '  Swap rows ', row1, ' and ', row2

  call r8row_swap ( m, n, a, row1, row2 )

  call r8mat_print ( m, n, a, '  The modified matrix:' )

  return
end
subroutine test108 ( )

!*****************************************************************************80
!
!! TEST108 tests R8ROW_TO_R8VEC.
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST108'
  write ( *, '(a)' ) '  R8ROW_TO_R8VEC converts an array of rows into a vector.'

  do i = 1, m
    do j = 1, n
      a(i,j) = real ( 10 * i + j, kind = 8 )
    end do
  end do

  call r8mat_print ( m, n, a, '  The array of rows:' )

  call r8row_to_r8vec ( m, n, a, x )

  call r8vec_print ( m*n, x, '  The resulting vector of rows:' )

  return
end
subroutine test109 ( )

!*****************************************************************************80
!
!! TEST109 tests R8ROW_VARIANCE.
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

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) variance(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST109'
  write ( *, '(a)' ) '  For a R8ROW (a matrix regarded as rows):'
  write ( *, '(a)' ) '  R8ROW_VARIANCE computes variances;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k )
    end do
  end do

  call r8mat_print ( m, n, a, '  The original matrix:' )

  call r8row_variance ( m, n, a, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Row variances:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,i3,3x,f10.4)' ) i, variance(i)
  end do

  return
end
subroutine test110 ( )

!*****************************************************************************80
!
!! TEST110 tests R8SLMAT_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ), allocatable, dimension ( : ) :: a
  real ( kind = 8 ), dimension ( 21 ) :: a1 = (/ &
    21.0D+00, 31.0D+00, 41.0D+00, 51.0D+00, 61.0D+00, 71.0D+00, &
              32.0D+00, 42.0D+00, 52.0D+00, 62.0D+00, 72.0D+00, &
                        43.0D+00, 53.0D+00, 63.0D+00, 73.0D+00, &
                                  54.0D+00, 64.0D+00, 74.0D+00, &
                                            65.0D+00, 75.0D+00, &
                                                      76.0D+00 /)
  real ( kind = 8 ), dimension ( 15 ) :: a2 = (/ &
    21.0D+00, 31.0D+00, 41.0D+00, 51.0D+00, 61.0D+00, 71.0D+00, &
              32.0D+00, 42.0D+00, 52.0D+00, 62.0D+00, 72.0D+00, &
                        43.0D+00, 53.0D+00, 63.0D+00, 73.0D+00 /)
  real ( kind = 8 ), dimension ( 6 ) :: a3 = (/ &
    21.0D+00, 31.0D+00, 41.0D+00, &
              32.0D+00, 42.0D+00, &
                        43.0D+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension ( test_num ) :: m_test = (/ 7, 7, 4 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), dimension ( test_num ) :: n_test = (/ 7, 3, 7 /)
  integer ( kind = 4 ) size
  integer ( kind = 4 ), dimension ( test_num ) :: size_test = (/ 21, 15, 6 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST110'
  write ( *, '(a)' ) '  R8SLMAT_PRINT prints a strictly lower triangular matrix'
  write ( *, '(a)' ) '  stored compactly.  Only the (possibly) nonzero '
  write ( *, '(a)' ) '  elements are printed.'

  do test = 1, test_num

    m = m_test(test)
    n = n_test(test)
    size = size_test(test)

    allocate ( a(1:size) )

    if ( test == 1 ) then
      a(1:size) = a1(1:size)
    else if ( test == 2 ) then
      a(1:size) = a2(1:size)
    else if ( test == 3 ) then
      a(1:size) = a3(1:size)
    end if

    call r8slmat_print ( m, n, a, '  R8SLMAT' )

    deallocate ( a )

  end do

  return
end
subroutine test111 ( )

!*****************************************************************************80
!
!! TEST111 tests R8VEC_AMAX and R8VEC_AMIN;
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

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST111'
  write ( *, '(a)' ) '  For a R8VEC:'
  write ( *, '(a)' ) '  R8VEC_AMAX:      maximum magnitude entry;'
  write ( *, '(a)' ) '  R8VEC_AMIN:      minimum magnitude entry.'

  b = - real ( n, kind = 8 )
  c = real ( n, kind = 8 )

  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print ( n, a, '  Input vector:' )

  write ( *, '(a)' ) ' '

  call r8vec_amax ( n, a, aval )
  write ( *, '(a,g14.6)' ) '  Maximum absolute:         ', aval

  call r8vec_amin ( n, a, aval )
  write ( *, '(a,g14.6)' ) '  Minimum absolute:         ', aval

  return
end
subroutine test112 ( )

!*****************************************************************************80
!
!! TEST112 tests R8VEC_BRACKET.
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
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  integer ( kind = 4 ) test
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( test_num ) :: xtest = (/ &
    -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /)
  real ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST112'
  write ( *, '(a)' ) '  R8VEC_BRACKET finds a pair of entries in a'
  write ( *, '(a)' ) '  sorted R8VEC which bracket a value.'

  call r8vec_indicator ( n, x )
  x(6) = x(5)

  call r8vec_print ( n, x, '  Sorted array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    LEFT             RIGHT'
  write ( *, '(a)' ) '  X(LEFT)   XVAL   X(RIGHT)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    xval = xtest(test)

    call r8vec_bracket ( n, x, xval, left, right )

    write ( *, '(i14,14x,i14)' ) left, right

    if ( 1 <= left .and. 1 <= right ) then
      write ( *, '(2x,3g14.6)' ) x(left), xval, x(right)
    else if ( left < 1 .and. 1 <= right ) then
      write ( *, '(2x,14x,2g14.6)' )          xval, x(right)
    else if ( 1 <= left .and. right < 1 ) then
      write ( *, '(2x,2g14.6)' ) x(left), xval
    else if ( left < 1 .and. right < 1 ) then
      write ( *, '(2x,14x,g14.6)' )          xval
    end if

  end do

  return
end
subroutine test113 ( )

!*****************************************************************************80
!
!! TEST113 tests R8VEC_BRACKET2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  integer ( kind = 4 ) start
  integer ( kind = 4 ) test
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( test_num ) :: xtest = (/ &
    -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /)
  real ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST113'
  write ( *, '(a)' ) '  R8VEC_BRACKET2 finds a pair of entries in a'
  write ( *, '(a)' ) '  sorted R8VEC which bracket a value.'

  call r8vec_indicator ( n, x )
  x(6) = x(5)

  call r8vec_print ( n, x, '  Sorted array:' )

  left = 0

  do test = 1, test_num

    xval = xtest(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Search for XVAL = ', xval

    if ( 0 < left ) then
      start = left
    else
      start = ( n + 1 ) / 2
    end if

    write ( *, '(a,i8)' ) '  Start = ', start

    call r8vec_bracket2 ( n, x, xval, start, left, right )

    write ( *, '(a,i8)' ) '  Left = ', left
    write ( *, '(a,i8)' ) '  Right = ', right

    if ( 1 <= left ) then
      write ( *, '(a,g14.6)' ) '  X(LEFT)=', x(left)
    end if

    if ( 1 <= right ) then
      write ( *, '(a,g14.6)' ) '  X(RIGHT) = ', x(right)
    end if

  end do

  return
end
subroutine test114 ( )

!*****************************************************************************80
!
!! TEST114 tests R8VEC_BRACKET3.
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
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) left
  integer ( kind = 4 ) test
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( test_num ) :: x_test = (/ &
    -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /)
  real ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST114'
  write ( *, '(a)' ) '  R8VEC_BRACKET3 finds a pair of entries in a'
  write ( *, '(a)' ) '  sorted R8VEC which bracket a value.'

  call r8vec_indicator ( n, x )
  x(6) = x(5)

  call r8vec_print ( n, x, '  Sorted array:' )

  left = ( n + 1 ) / 2

  do test = 1, test_num

    xval = x_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Search for XVAL = ', xval

    write ( *, '(a,i8)' ) '  Starting guess for interval is = ', left

    call r8vec_bracket3 ( n, x, xval, left )

    write ( *, '(a)' ) '  Nearest interval:'
    write ( *, '(a,i8,a,g14.6)' ) '    X[', left,' ]= ', x(left)
    write ( *, '(a,i8,a,g14.6)' ) '    X[', left+1, ' ]= ', x(left+1)

  end do

  return
end
subroutine test1143 ( )

!*****************************************************************************80
!
!! TEST1143 tests R8VEC_BRACKET5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) left
  integer ( kind = 4 ) r8vec_bracket5
  integer ( kind = 4 ) right
  integer ( kind = 4 ) test
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( test_num ) :: xtest = (/ &
    -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /)
  real ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1143'
  write ( *, '(a)' ) '  R8VEC_BRACKET5 finds a pair of entries in a'
  write ( *, '(a)' ) '  sorted R8VEC which bracket a value.'

  call r8vec_indicator ( n, x )
  x(6) = x(5)

  call r8vec_print ( n, x, '  Sorted array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        LEFT                   RIGHT'
  write ( *, '(a)' ) '      X(LEFT)       XVAL     X(RIGHT)'
  write ( *, '(a)' ) ' ' 

  do test = 1, test_num

    xval = xtest(test)

    left = r8vec_bracket5 ( n, x, xval )

    if ( left == -1 ) then
      write ( *, '(2x,i10)' ) left
      write ( *, '(2x,10x,2x,f10.4,2x,a)' ) xval, '(Not bracketed!)'
    else
      right = left + 1
      write ( *, '(2x,i10,2x,10x,2x,i10)' ) left, right
      write ( *, '(2x,f10.4,2x,f10.4,2x,f10.4)' ) x(left), xval, x(right)
    end if

  end do

  return
end
subroutine test1145 ( )

!*****************************************************************************80
!
!! TEST1145 tests R8VEC_CHEBYSPACE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: r(:)
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1145'
  write ( *, '(a)' ) '  R8VEC_CHEBYSPACE computes N Chebyshev points in [R1,R2].'

  r1 = -1.0D+00
  r2 = +1.0D+00
  n = 5
  allocate ( r(1:n) )

  call r8vec_chebyspace ( n, r1, r2, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a,g14.6,a,g14.6)' ) '  N = ', n, '  R1 = ', r1, '  R2 = ', r2

  call r8vec_print ( n, r, '  Chebyshev points:' )

  deallocate ( r )

  r1 =   0.0D+00
  r2 = +10.0D+00
  n = 7

  allocate ( r(1:n) )

  call r8vec_chebyspace ( n, r1, r2, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a,g14.6,a,g14.6)' ) '  N = ', n, '  R1 = ', r1, '  R2 = ', r2

  call r8vec_print ( n, r, '  Chebyshev points:' )

  deallocate ( r )

  return
end
subroutine test1147 ( )

!*****************************************************************************80
!
!! TEST1147 tests R8VEC_CONVOLUTION
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension(m) :: x = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00 /)
  real ( kind = 8 ), dimension(n) :: y = (/ &
   -1.0D+00, 5.0D+00, 3.0D+00 /)
  real ( kind = 8 ) z(m+n-1)
  real ( kind = 8 ), dimension (m+n-1) :: z_correct = (/ &
    -1.0D+00, 3.0D+00, 10.0D+00, 17.0D+00, 29.0D+00, 12.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1147'
  write ( *, '(a)' ) '  R8VEC_CONVOLUTION computes the convolution'
  write ( *, '(a)' ) '  of two vectors.'

  call r8vec_print ( m, x, '  The factor X:' )
  call r8vec_print ( n, y, '  The factor Y:' )

  call r8vec_convolution ( m, x, n, y, z )

  call r8vec_print ( m + n - 1, z, '  The convolution z = x star y:' )

  call r8vec_print ( m + n - 1, z_correct, '  Correct answer:' )

  return
end
subroutine test115 ( )

!*****************************************************************************80
!
!! TEST115 tests R8VEC_CONVOLUTION_CIRC
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

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ), dimension(n) :: x = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00 /)
  real ( kind = 8 ), dimension(n) :: y = (/ &
    1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00 /)
  real ( kind = 8 ) z(n)
  real ( kind = 8 ), dimension ( n ) :: z_correct = (/ &
    37.0D+00, 44.0D+00, 43.0D+00, 26.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST115'
  write ( *, '(a)' ) '  R8VEC_CONVOLUTION_CIRC computes the circular convolution'
  write ( *, '(a)' ) '  of two vectors.'

  call r8vec_print ( n, x, '  The factor X:' )
  call r8vec_print ( n, y, '  The factor Y:' )

  call r8vec_convolution_circ ( n, x, y, z )

  call r8vec_print ( n, z, '  The circular convolution z = x CC y:' )

  call r8vec_print ( n, z_correct, '  Correct answer:' )

  return
end
subroutine test116 ( )

!*****************************************************************************80
!
!! TEST116 tests R8VEC_DIF.
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

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) cof(0:n)
  real ( kind = 8 ) fdif
  real ( kind = 8 ) :: h = 0.01D+00
  integer ( kind = 4 ) i
  real ( kind = 8 ) test116_f
  real ( kind = 8 ) :: x = 1.0D+00
  real ( kind = 8 ) xi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST116'
  write ( *, '(a)' ) '  R8VEC_DIF estimates derivatives.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Estimate the derivative of order N = ', n
  write ( *, '(a,g14.6)' ) '  Using H = ', h
  write ( *, '(a,g14.6)' ) '  at argument X = ', x
!
!  Get the coefficients.
!
  call r8vec_dif ( n, h, cof )

  call r8vec_print ( n+1, cof, '  The difference coefficients:' )

  fdif = 0.0D+00
  do i = 0, n
    xi = x + real ( 2 * i - n, kind = 8 ) * h
    fdif = fdif + cof(i) * test116_f ( xi )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimate is FDIF = ', fdif

  return
end
function test116_f ( x )

!*****************************************************************************80
!
!! TEST116_F evaluates the function used in TEST116.
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

  real ( kind = 8 ) test116_f
  real ( kind = 8 ) x

  test116_f = exp ( x )

  return
end
subroutine test1165 ( )

!*****************************************************************************80
!
!! TEST1165 tests R8VEC_DIRECT_PRODUCT.
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
  real ( kind = 8 ), allocatable, dimension ( : ) :: factor_value
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(factor_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1165'
  write ( *, '(a)' ) '  R8VEC_DIRECT_PRODUCT forms the entries of a'
  write ( *, '(a)' ) '  direct product of a given number of R8VEC factors.'

  x(1:factor_num,1:point_num) = 0.0D+00

  do factor_index = 1, factor_num

    if ( factor_index == 1 ) then
      factor_order = 4
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00 /)
    else if ( factor_index == 2 ) then
      factor_order = 3
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 50.0D+00, 60.0D+00, 70.0D+00 /)
    else if ( factor_index == 3 ) then
      factor_order = 2
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 800.0D+00, 900.0D+00 /)
    end if

    call r8vec_direct_product ( factor_index, factor_order, factor_value,  &
      factor_num, point_num, x )

    deallocate ( factor_value )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     J         X(1)      X(2)      X(3)'
  write ( *, '(a)' ) ' '

  do j = 1, point_num
    write ( *, '(2x,i4,4x,f8.1,2x,f8.1,2x,f8.1)' ) j, x(1:factor_num,j)
  end do

  return
end
subroutine test1166 ( )

!*****************************************************************************80
!
!! TEST1166 tests R8VEC_DIRECT_PRODUCT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2007
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
  real ( kind = 8 ), allocatable, dimension ( : ) :: factor_value
  integer ( kind = 4 ) j
  real ( kind = 8 ) w(point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1166'
  write ( *, '(a)' ) '  R8VEC_DIRECT_PRODUCT2 forms the entries of a'
  write ( *, '(a)' ) '  direct product of a given number of R8VEC factors.'

  w(1:point_num) = 1.0D+00

  do factor_index = 1, factor_num

    if ( factor_index == 1 ) then
      factor_order = 4
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 2.0D+00, 3.0D+00, 5.0D+00, 7.0D+00 /)
    else if ( factor_index == 2 ) then
      factor_order = 3
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 11.0D+00, 13.0D+00, 17.0D+00 /)
    else if ( factor_index == 3 ) then
      factor_order = 2
      allocate ( factor_value(1:factor_order) )
      factor_value = (/ 19.0D+00, 21.0D+00 /)
    end if

    call r8vec_direct_product2 ( factor_index, factor_order, factor_value,  &
      factor_num, point_num, w )

    deallocate ( factor_value )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     J         W(J)'
  write ( *, '(a)' ) ' '

  do j = 1, point_num
    write ( *, '(2x,i4,4x,f8.1)' ) j, w(j)
  end do

  return
end
subroutine test117 ( )

!*****************************************************************************80
!
!! TEST117 tests R8VEC_EVEN.
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

  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo

  xlo = 0.0D+00
  xhi = 99.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST117'
  write ( *, '(a)' ) '  R8VEC_EVEN computes N evenly spaced values'
  write ( *, '(a)' ) '  between XLO and XHI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  XLO = ', xlo
  write ( *, '(a,g14.6)' ) '  XHI = ', xhi
  write ( *, '(a,i8)' ) '  while N = ', n

  call r8vec_even ( n, xlo, xhi, x )

  call r8vec_print ( n, x, '  Resulting array:' )

  return
end
subroutine test118 ( )

!*****************************************************************************80
!
!! TEST118 tests R8VEC_EVEN2.
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

  integer ( kind = 4 ), parameter :: nold = 5
  integer ( kind = 4 ), parameter :: maxval = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) istar
  integer ( kind = 4 ) jstar
  integer ( kind = 4 ), dimension(nold-1) :: nfill = (/ 4, 3, 5, 0 /)
  integer ( kind = 4 ) nval
  real ( kind = 8 ), dimension(nold) :: xold = (/ &
    0.0D+00, 1.0D+00, 5.0D+00, 2.0D+00, 0.0D+00 /)
  real ( kind = 8 ) xval(maxval)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST118'
  write ( *, '(a)' ) '  R8VEC_EVEN2 interpolates a specified number of '
  write ( *, '(a)' ) '  points pairs of values in a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input data:'
  write ( *, '(a)' ) ' '
  do i = 1, nold
    write ( *, '(2x,g14.6)' ) xold(i)
    if ( i < nold ) then
      write ( *, '(2x,a,i10,a)' ) '(', nfill(i), ')'
    end if
  end do

  call r8vec_even2 ( maxval, nfill, nold, xold, nval, xval )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Resulting vector:'
  write ( *, '(a)' ) ' '

  istar = 1
  jstar = 1
  do i = 1, nval

    if ( i == istar ) then

      write ( *, '(2x,a1,g14.6)' ) '*', xval(i)

      if ( jstar < nold ) then
        istar = istar + nfill(jstar) + 1
        jstar = jstar + 1
      end if

    else

      write ( *, '(2x,g14.6)' ) xval(i)

    end if

  end do

  return
end
subroutine test119 ( )

!*****************************************************************************80
!
!! TEST119 tests R8VEC_EVEN3.
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

  integer ( kind = 4 ), parameter :: nold = 4
  integer ( kind = 4 ), parameter :: nval = 12

  real ( kind = 8 ), dimension(nold) :: xold = (/ &
    0.0D+00, 5.1D+00, 7.0D+00, 10.0D+00 /)
  real ( kind = 8 ) xval(nval)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST119'
  write ( *, '(a)' ) '  R8VEC_EVEN3 tries to evenly interpolate new data'
  write ( *, '(a)' ) '  between old values.'

  call r8vec_print ( nold, xold, '  Original vector:' )

  call r8vec_even3 ( nold, nval, xold, xval )

  call r8vec_print ( nval, xval, '  New vector:' )

  return
end
subroutine test120 ( )

!*****************************************************************************80
!
!! TEST120 tests R8VEC_EXPAND_LINEAR.
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

  integer ( kind = 4 ), parameter :: n = 6
  integer ( kind = 4 ), parameter :: fat = 3
  integer ( kind = 4 ), parameter :: nfat = ( n - 1 ) * ( fat + 1 ) + 1

  real ( kind = 8 ), dimension(n) :: x =  (/ &
    16.0D+00, 4.0D+00, 0.0D+00, 4.0D+00, 16.0D+00, 36.0D+00 /)
  real ( kind = 8 ) xfat(nfat)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST120'
  write ( *, '(a)' ) '  R8VEC_EXPAND_LINEAR linearly interpolates new data'
  write ( *, '(a)' ) '  between old values.'

  call r8vec_print ( n, x, '  Original vector:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Expansion factor is ', fat

  call r8vec_expand_linear ( n, x, fat, xfat )

  call r8vec_print ( nfat, xfat, '  Fattened vector:' )

  return
end
subroutine test121 ( )

!*****************************************************************************80
!
!! TEST121 tests R8VEC_FRAC.
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

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) afrac
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST121'
  write ( *, '(a)' ) '  R8VEC_FRAC: K-th smallest R8VEC entry;'

  seed = 123456789

  call r8vec_uniform_01 ( n, seed, a )

  call r8vec_print ( n, a, '  Array to search:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Fractile  Value '
  write ( *, '(a)' ) ' '

  do k = 1, n, n/2

    call r8vec_frac ( n, a, k, afrac )

    write ( *, '(2x,i8,2x,g14.6)' ) k, afrac

  end do

  return
end
subroutine test1215 ( )

!*****************************************************************************80
!
!! TEST1215 tests R8VEC_HEAP_D_EXTRACT, R8VEC_HEAP_D_INSERT and R8VEC_HEAP_D_MAX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  real ( kind = 8 ) a(n_max)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1215'
  write ( *, '(a)' ) '  For a heap descending sorted R8VEC,'
  write ( *, '(a)' ) '  R8VEC_HEAP_D_INSERT inserts a value into the heap.'
  write ( *, '(a)' ) '  R8VEC_HEAP_D_EXTRACT extracts the maximum value;'
  write ( *, '(a)' ) '  R8VEC_HEAP_D_MAX reports the maximum value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  These 3 operations are enough to model a priority queue.'

  n = 0

  seed = 123456789

  do i = 1, n_max

    b = 0.0D+00
    c = 10.0D+00

    value = r8_uniform_ab ( b, c, seed )

    call r8vec_heap_d_insert ( n, a, value )

    write ( *, '(a)' ) ' '
    write ( *, '(a,f10.4)' ) '  Inserting value          ', value

    call r8vec_heap_d_max ( n, a, value )

    write ( *, '(a,f10.4)' ) '  Current maximum value is ', value

  end do

  call r8vec_print ( n, a, '  Current heap as a vector:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now extract the maximum several times.'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call r8vec_heap_d_extract ( n, a, value )
    write ( *, '(a,f10.4)' ) '  Extracting maximum element = ', value
  end do

  call r8vec_print ( n, a, '  Current heap as a vector:' )

  return
end
subroutine test122 ( )

!*****************************************************************************80
!
!! TEST122 tests R8VEC_HISTOGRAM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: histo_num = 20
  integer ( kind = 4 ), parameter :: n = 1000

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) :: a_hi
  real ( kind = 8 ) :: a_lo
  real ( kind = 8 ) bin_hi
  real ( kind = 8 ) bin_lo
  integer ( kind = 4 ), dimension (0:histo_num+1) :: histo_gram
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST122'
  write ( *, '(a)' ) '  R8VEC_HISTOGRAM histograms a real vector.'

  do test = 1, test_num

    if ( test == 1 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Uniform data:'

      a_lo =  0.0D+00
      a_hi = +1.0D+00
      call r8vec_uniform_01 ( n, seed, a )

    else if ( test == 2 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Normal data:'
      a_lo = -3.0D+00
      a_hi = +3.0D+00
      call r8vec_normal_01 ( n, seed, a )

    end if

    call r8vec_histogram ( n, a, a_lo, a_hi, histo_num, histo_gram )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Histogram of data:'
    write ( *, '(a)' ) ' '

    do i = 0, histo_num+1

      if ( i == 0 ) then

        write ( *, '(2x,10x,2x,f10.4,2x,i8)' ) a_lo, histo_gram(i)

      else if ( i <= histo_num ) then

        bin_lo = ( real ( histo_num - i + 1, kind = 8 ) * a_lo   &
                 + real (             i - 1, kind = 8 ) * a_hi ) &
                 / real ( histo_num,         kind = 8 )

        bin_hi = ( real ( histo_num - i,     kind = 8 ) * a_lo   &
                 + real (             i,     kind = 8 ) * a_hi ) &
                 / real ( histo_num,         kind = 8 )

        write ( *, '(2x,f10.4,2x,f10.4,2x,i8)' ) bin_lo, bin_hi, histo_gram(i)

      else if ( i == histo_num+1 ) then

        write ( *, '(2x,f10.4,2x,10x,2x,i8)' ) a_hi, histo_gram(i)

      end if

    end do

  end do

  return
end
subroutine test123 ( )

!*****************************************************************************80
!
!! TEST123 tests R8VEC_INDEX_INSERT, R8VEC_INDEX_DELETE_ALL, R8VEC_INDEX_DELETE_DUPES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) xval

  n = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST123'
  write ( *, '(a)' ) '  R8VEC_INDEX_INSERT inserts values into an'
  write ( *, '(a)' ) '  index sorted array.'
  write ( *, '(a)' ) '  R8VEC_INDEX_DELETE_ALL deletes all copies of a'
  write ( *, '(a)' ) '  particular value.'
  write ( *, '(a)' ) '  R8VEC_INDEX_DELETE_ONE deletes one copies of a'
  write ( *, '(a)' ) '  particular value.'
  write ( *, '(a)' ) '  R8VEC_INDEX_DELETE_DUPES deletes duplicates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate some random values:'
  write ( *, '(a)' ) ' '

  xval = 8.0D+00
  call r8vec_index_insert ( n, x, indx, xval )

  xval = 7.0D+00
  call r8vec_index_insert ( n, x, indx, xval )

  seed = 123456789

  do i = 1, 20
    xval = r8_uniform_ab ( 0.0D+00, 20.0D+00, seed )
    xval = real ( nint ( xval ), kind = 8 )
    write ( *, '(4x,f6.2)' ) xval
    call r8vec_index_insert ( n, x, indx, xval )
  end do

  xval = 7.0D+00
  call r8vec_index_insert ( n, x, indx, xval )

  xval = 8.0D+00
  call r8vec_index_insert ( n, x, indx, xval )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,f6.2,9x,f6.2)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call R8VEC_INDEX_DELETE_ONE to delete one value of 8:'

  xval = 8.0D+00
  call r8vec_index_delete_one ( n, x, indx, xval, n, x, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call R8VEC_INDEX_DELETE_ALL to delete all values of 7:'

  xval = 7.0D+00
  call r8vec_index_delete_all ( n, x, indx, xval )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,f6.2,9x,f6.2)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call R8VEC_INDEX_DELETE_DUPES to delete duplicates:'

  call r8vec_index_delete_dupes ( n, x, indx, n, x, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of unique entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,f6.2)' ) i, indx(i), x(i)
  end do

  return
end
subroutine test124 ( )

!*****************************************************************************80
!
!! TEST124 tests R8VEC_INDEX_INSERT_UNIQUE and R8VEC_INDEX_ORDER.
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

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) xval

  n = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST124'
  write ( *, '(a)' ) '  R8VEC_INDEX_INSERT_UNIQUE inserts unique values into an'
  write ( *, '(a)' ) '  index sorted array.'
  write ( *, '(a)' ) '  R8VEC_INDEX_ORDER sorts an index sorted array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate some random values:'
  write ( *, '(a)' ) ' '
  seed = 123456789

  do i = 1, N_MAX
    xval = r8_uniform_ab ( 0.0D+00, 20.0D+00, seed )
    xval = real ( nint ( xval ), kind = 8 )
    write ( *, '(4x,f6.2)' ) xval
    call r8vec_index_insert_unique ( n, x, indx, xval )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of unique entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,f6.2,9x,f6.2)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call R8VEC_INDEX_ORDER to carry out the sorting:'

  call r8vec_index_order ( n, x, indx )

  call r8vec_print ( n, x, '  X:' )

  return
end
subroutine test125 ( )

!*****************************************************************************80
!
!! TEST125 tests R8VEC_INDEX_INSERT_UNIQUE and R8VEC_INDEX_SEARCH.
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

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) xval

  n = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST125'
  write ( *, '(a)' ) '  R8VEC_INDEX_INSERT_UNIQUE inserts unique values into an'
  write ( *, '(a)' ) '  index sorted array.'
  write ( *, '(a)' ) '  R8VEC_INDEX_SEARCH searches for an entry '
  write ( *, '(a)' ) '  with a given value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate some random values:'
  write ( *, '(a)' ) ' '

  b = 0.0D+00
  c = real ( n_max, kind = 8 )
  seed = 123456789

  do i = 1, n_max
    xval = r8_uniform_ab ( b, c, seed )
    xval = real ( nint ( xval ), kind = 8 )
    write ( *, '(4x,f6.2)' ) xval
    call r8vec_index_insert_unique ( n, x, indx, xval )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indexed list of entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,6x,i3,3x,f6.2,9x,f6.2)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Results of search for given XVAL:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XVAL  Less Equal More'
  write ( *, '(a)' ) ' '

  do i = 0, n_max
    xval = real ( i )
    call r8vec_index_search ( n, x, indx, xval, less, equal, more )
    write ( *, '(2x,f6.2,3x,i3,3x,i3,3x,i3)' ) xval, less, equal, more
  end do

  return
end
subroutine test1251 ( )

!*****************************************************************************80
!
!! TEST1251 tests R8VEC_INDEX_SORTED_RANGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  integer ( kind = 4 ) indx(n)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r_lo
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1251'
  write ( *, '(a)' ) '  R8VEC_INDEX_SORTED_RANGE seeks the range I_LO:I_HI'
  write ( *, '(a)' ) '  of entries of sorted indexed R so that'
  write ( *, '(a)' ) '  R_LO <= R(INDX(I)) <= R_HI for I_LO <= I <= I_HI.'

  seed = 123456789

  do test = 1, 5

    call r8vec_uniform_01 ( n, seed, r )

    call r8vec_print ( n, r, '  Array' )

    call r8vec_sort_heap_index_a ( n, r, indx )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I  INDX    R(INDX(I))'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, indx(i), r(indx(i))
    end do

    r_lo = r8_uniform_01 ( seed )
    r_hi = r8_uniform_01 ( seed )

    if ( r_hi < r_lo ) then
      t = r_lo
      r_lo = r_hi
      r_hi = t
    end if

    call r8vec_index_sorted_range ( n, r, indx, r_lo, r_hi, i_lo, i_hi )

    write ( *, '(a)' ) ' '
    if ( i_hi < i_lo ) then
      write ( *, '(2x,a4,2x,6x,g14.6)' ) 'R_LO', r_lo
      write ( *, '(2x,a4,2x,6x,g14.6)' ) 'R_HI', r_hi
      write ( *, '(a)' ) '  Empty range in R.'
    else

      write ( *, '(2x,a4,2x,6x,g14.6)' ) 'R_LO', r_lo
      do i = i_lo, i_hi
        write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, indx(i), r(indx(i))
      end do
      write ( *, '(2x,a4,2x,6x,g14.6)' ) 'R_HI', r_hi
    end if

  end do

  return
end
subroutine test1252 ( )

!*****************************************************************************80
!
!! TEST1252 tests R8VEC_INDEXED_HEAP_D;
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

  real ( kind = 8 ) :: a(m) = (/ &
    101.0D+00, 102.0D+00, 103.0D+00, 104.0D+00, 105.0D+00, &
    106.0D+00, 107.0D+00, 108.0D+00, 109.0D+00, 110.0D+00, &
    111.0D+00, 112.0D+00, 113.0D+00, 114.0D+00, 115.0D+00, &
    116.0D+00, 117.0D+00, 118.0D+00, 119.0D+00, 120.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: indx(n) = (/ &
    1, 11, 17, 5, 7, 13, 15, 3, 19, 9 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1252'
  write ( *, '(a)' ) '  R8VEC_INDEXED_HEAP_D creates a descending heap'
  write ( *, '(a)' ) '  from an indexed R8VEC.'
!
!  Print before.
!
  call r8vec_print ( m, a, '  The data vector:' )
  call i4vec_print ( n, indx, '  The index vector:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX):'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
  end do
!
!  Create the heap.
!
  call r8vec_indexed_heap_d ( n, a, indx )
!
!  Print afterwards.  Only INDX should change.
!
  call r8vec_print ( m, a, '  The data vector (should NOT change):' )
  call i4vec_print ( n, indx, '  The index vector (may change):' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX) is now a heap:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
  end do

  return
end
subroutine test1255 ( )

!*****************************************************************************80
!
!! TEST1255 tests R8VEC_INDEXED_HEAP_D_EXTRACT and related routines.
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

  real ( kind = 8 ) a(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) indx_extract
  integer ( kind = 4 ) indx_insert
  integer ( kind = 4 ) indx_max
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1255'
  write ( *, '(a)' ) '  For an indexed R8VEC,'
  write ( *, '(a)' ) '  R8VEC_INDEXED_HEAP_D_INSERT inserts a value into the heap.'
  write ( *, '(a)' ) '  R8VEC_INDEXED_HEAP_D_EXTRACT extracts the maximum value;'
  write ( *, '(a)' ) '  R8VEC_INDEXED_HEAP_D_MAX reports the maximum value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  These 3 operations are enough to model a priority queue.'
!
!  Set the data array.  To keep things easy, we will use the indicator vector.
!
  call r8vec_indicator ( m, a )
!
!  The index array will initially be a random subset of the numbers 1 to M,
!  in random order.
!
  n = 5
  indx(1:11) = (/ 9, 2, 8, 14, 5, 7, 15, 1, 19, 20, 3 /)

  call r8vec_print ( m, a, '  The data vector:' )
  call i4vec_print ( n, indx, '  The index vector:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX):'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
  end do
!
!  Create the descending heap.
!
  call r8vec_indexed_heap_d ( n, a, indx )

  call i4vec_print ( n, indx, '  The index vector after heaping:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX) after heaping:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
  end do
!
!  Insert five entries, and monitor the maximum.
!
  do i = 1, 5

    indx_insert = indx(n+1)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Inserting value ', a(indx_insert)

    call r8vec_indexed_heap_d_insert ( n, a, indx, indx_insert )

    call r8vec_indexed_heap_d_max ( n, a, indx, indx_max )

    write ( *, '(a,g14.6)' ) '  Current maximum is ', a(indx_max)

  end do
  call r8vec_print ( m, a, '  The data vector after insertions:' )
  call i4vec_print ( n, indx, '  The index vector after insertions:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX) after insertions:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
  end do
!
!  Extract the first 5 largest elements.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now extract the maximum several times.'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call r8vec_indexed_heap_d_extract ( n, a, indx, indx_extract )
    write ( *, '(a,i8,a,g14.6)' ) '  Extracting maximum element A(', &
      indx_extract,') = ', a(indx_extract)
  end do

  call r8vec_print ( m, a, '  The data vector after extractions:' )
  call i4vec_print ( n, indx, '  The index vector after extractions:' )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(INDX) after extractions:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
  end do

  return
end
subroutine test1256 ( )

!*****************************************************************************80
!
!! TEST1256 tests R8VEC_LEGENDRE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: r(:)
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1256'
  write ( *, '(a)' ) '  R8VEC_LEGENDRE computes N Legendre points in [R1,R2].'

  r1 = -1.0D+00
  r2 = +1.0D+00
  n = 5
  allocate ( r(1:n) )

  call r8vec_legendre ( n, r1, r2, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a,g14.6,a,g14.6)' ) '  N = ', n, '  R1 = ', r1, '  R2 = ', r2

  call r8vec_print ( n, r, '  Legendre points:' )

  deallocate ( r )

  r1 =   0.0D+00
  r2 = +10.0D+00
  n = 7

  allocate ( r(1:n) )

  call r8vec_legendre ( n, r1, r2, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a,g14.6,a,g14.6)' ) '  N = ', n, '  R1 = ', r1, '  R2 = ', r2

  call r8vec_print ( n, r, '  Legendre points:' )

  deallocate ( r )

  return
end
subroutine test1258 ( )

!*****************************************************************************80
!
!! TEST1258 tests R8VEC_LINSPACE, R8VEC_LINSPACE2 and R8VEC_MIDSPACE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1258'
  write ( *, '(a)' ) '  For a R8VEC:'
  write ( *, '(a)' ) '  R8VEC_LINSPACE: evenly spaced points between A and B;'
  write ( *, '(a)' ) '  R8VEC_LINSPACE2: R8VEC_LINSPACE, but no endpoints;'
  write ( *, '(a)' ) '  R8VEC_MIDSPACE: evenly spaced midpoints between A and B'

  a = 10.0D+00
  b = 20.0D+00

  call r8vec_linspace ( n, a, b, x )
  call r8vec_print ( n, x, '  r8vec_linspace ( 5, 10, 20 )' )

  call r8vec_linspace2 ( n, a, b, x )
  call r8vec_print ( n, x, '  r8vec_linspace2 ( 5, 10, 20 )' )

  call r8vec_midspace ( n, a, b, x )
  call r8vec_print ( n, x, '  r8vec_midspace ( 5, 10, 20 )' )


  return
end
subroutine test126 ( )

!*****************************************************************************80
!
!! TEST126 tests R8VEC_MAX and R8VEC_MIN.
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

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST126'
  write ( *, '(a)' ) '  For a R8VEC:'
  write ( *, '(a)' ) '  R8VEC_MAX:       maximum entry;'
  write ( *, '(a)' ) '  R8VEC_MIN:       minimum entry.'

  b = - real ( n, kind = 8 )
  c = real ( n, kind = 8 )

  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print ( n, a, '  Input vector:' )

  write ( *, '(a)' ) ' '

  call r8vec_max ( n, a, aval )
  write ( *, '(a,g14.6)' ) '  Maximum:                 ', aval

  call r8vec_min ( n, a, aval )
  write ( *, '(a,g14.6)' ) '  Minimum:                 ', aval

  return
end
subroutine test127 ( )

!*****************************************************************************80
!
!! TEST127 tests R8VEC_MAX_INDEX and R8VEC_MIN_INDEX.
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

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST127'
  write ( *, '(a)' ) '  For a R8VEC:'
  write ( *, '(a)' ) '  R8VEC_MAX_INDEX: index of maximum entry;'
  write ( *, '(a)' ) '  R8VEC_MIN_INDEX: index of minimum entry;'

  b = - real ( n, kind = 8 )
  c = real ( n, kind = 8 )

  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print ( n, a, '  Input vector:' )

  write ( *, '(a)' ) ' '

  call r8vec_max_index ( n, a, ival )
  write ( *, '(a,i8)' ) '  Maximum index:           ', ival

  call r8vec_min_index ( n, a, ival )
  write ( *, '(a,i8)' ) '  Minimum index:           ', ival

  return
end
subroutine test128 ( )

!*****************************************************************************80
!
!! TEST128 tests R8VEC_MEAN and R8VEC_MEDIAN;
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

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean
  real ( kind = 8 ) median
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST128'
  write ( *, '(a)' ) '  For a R8VEC:'
  write ( *, '(a)' ) '  R8VEC_MEAN:      mean value;'
  write ( *, '(a)' ) '  R8VEC_MEDIAN:    median value;'

  b = - real ( n, kind = 8 )
  c = real ( n, kind = 8 )

  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print ( n, a, '  Input vector:' )

  write ( *, '(a)' ) ' '

  call r8vec_mean ( n, a, mean )
  write ( *, '(a,g14.6)' ) '  Mean:                     ', mean
  call r8vec_median ( n, a, median )
  write ( *, '(a,g14.6)' ) '  Median:                   ', median

  return
end
subroutine test129 ( )

!*****************************************************************************80
!
!! TEST129 tests R8VEC_NORM_L1, R8VEC_NORM_L2, R8VEC_NORM_LI.
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

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r8vec_norm_l1
  real ( kind = 8 ) r8vec_norm_l2
  real ( kind = 8 ) r8vec_norm_li
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST129'
  write ( *, '(a)' ) '  For a R8VEC:'
  write ( *, '(a)' ) '  R8VEC_NORM_L1:   L1 norm.'
  write ( *, '(a)' ) '  R8VEC_NORM_L2:   L2 norm.'
  write ( *, '(a)' ) '  R8VEC_NORM_LI:   L-infinity norm.'

  b = - real ( n, kind = 8 )
  c = real ( n, kind = 8 )

  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print ( n, a, '  Input vector:' )

  write ( *, '(a)' ) ' '

  write ( *, '(a,g14.6)' ) '  L1 norm:                 ', r8vec_norm_l1 ( n, a )
  write ( *, '(a,g14.6)' ) '  L2 norm:                 ', r8vec_norm_l2 ( n, a )
  write ( *, '(a,g14.6)' ) '  L-Infinity norm:         ', r8vec_norm_li ( n, a )

  return
end
subroutine test130 ( )

!*****************************************************************************80
!
!! TEST130 tests R8VEC_NORMAL_01.
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

  integer ( kind = 4 ), parameter :: n_max = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_mean
  real ( kind = 8 ) x_min
  real ( kind = 8 ) x_var

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST130'
  write ( *, '(a)' ) '  R8VEC_NORMAL_01 computes a vector of normally'
  write ( *, '(a)' ) '  distributed random numbers.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed
!
!  Test 1:
!  Simply call 5 times for 1 value, and print.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #1: Call 5 times, 1 value each time.'
  write ( *, '(a)' ) ' '

  n = 1
  do i = 1, 5
    call r8vec_normal_01 ( n, seed, x )
    write ( *, '(2x,i8,g14.6)' ) i, x(1)
  end do
!
!  Test 2:
!  Restore the random number seed, and repeat.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #2: Restore the random number seed.'
  write ( *, '(a)' ) '  Call 5 times, 1 value each time.'
  write ( *, '(a)' ) '  The results should be identical.'
  write ( *, '(a)' ) ' '

  n = -1
  call r8vec_normal_01  ( n, seed, x )

  seed = 123456789

  n = 1
  do i = 1, 5
    call r8vec_normal_01 ( n, seed, x )
    write ( *, '(2x,i8,g14.6)' ) i, x(1)
  end do
!
!  Test 3:
!  Restore the random number seed, compute all 5 values at once.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #3: Restore the random number seed.'
  write ( *, '(a)' ) '  Call 1 time for 5 values.'
  write ( *, '(a)' ) '  The results should be identical.'
  write ( *, '(a)' ) ' '

  n = -1
  call r8vec_normal_01 ( n, seed, x )

  seed = 123456789

  n = 5
  call r8vec_normal_01 ( n, seed, x )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do
!
!  Test 4:
!  Restore the random number seed, compute all 5 values at once.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #4: Restore the random number seed.'
  write ( *, '(a)' ) '  Call for 2, 1, and 2 values.'
  write ( *, '(a)' ) '  The results should be identical.'
  write ( *, '(a)' ) ' '

  n = -1
  call r8vec_normal_01 ( n, seed, x )

  seed = 123456789

  n = 2
  call r8vec_normal_01 ( n, seed, x )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do

  n = 1
  call r8vec_normal_01 ( n, seed, x )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do

  n = 2
  call r8vec_normal_01 ( n, seed, x )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do
!
!  Test 5:
!  Determine the minimum, maximum, mean and variance.
!
  n = n_max
  call r8vec_normal_01 ( n, seed, x )
  x_min = minval ( x(1:n) )
  x_max = maxval ( x(1:n) )
  x_mean = sum ( x(1:n) ) / real ( n, kind = 8 )
  x_var = sum ( ( x(1:n) - x_mean )**2 ) / real ( n - 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #5:'
  write ( *, '(a,i12)' ) '  Number of samples was ', n
  write ( *, '(a,g14.6)' ) '  Minimum value was ', x_min
  write ( *, '(a,g14.6)' ) '  Maximum value was ', x_max
  write ( *, '(a,g14.6)' ) '  Average value was ', x_mean
  write ( *, '(a,g14.6)' ) '  Variance was      ', x_var
  write ( *, '(a,g14.6)' ) '  Expected average  ', 0.0D+00
  write ( *, '(a,g14.6)' ) '  Expected variance ', 1.0D+00

  return
end
subroutine test152 ( )

!*****************************************************************************80
!
!! TEST152 tests R8VEC_NORMALIZE_L1.
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

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST152'
  write ( *, '(a)' ) '  For a R8VEC:'
  write ( *, '(a)' ) '  R8VEC_NORMALIZE_L1:  make unit sum;'

  b = - real ( n, kind = 8 )
  c = real ( n, kind = 8 )

  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print ( n, a, '  Input vector:' )

  call r8vec_normalize_l1 ( n, a )

  call r8vec_print ( n, a, '  After calling R8VEC_NORMALIZE_L1:' )

  return
end
subroutine test131 ( )

!*****************************************************************************80
!
!! TEST131 tests R8VEC_ORDER_TYPE.
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

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) j
  integer ( kind = 4 ) order
  integer ( kind = 4 ) test
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( n, test_num ) :: x_test = reshape ( (/ &
    1.0D+00, 3.0D+00, 2.0D+00, 4.0D+00, &
    2.0D+00, 2.0D+00, 2.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00, 2.0D+00, 4.0D+00, &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, &
    4.0D+00, 4.0D+00, 3.0D+00, 1.0D+00, &
    9.0D+00, 7.0D+00, 3.0D+00, 0.0D+00 /), (/ n, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST131'
  write ( *, '(a)' ) '  R8VEC_ORDER_TYPE classifies a R8VEC as'
  write ( *, '(a)' ) '  -1: no order'
  write ( *, '(a)' ) '   0: all equal;'
  write ( *, '(a)' ) '   1: ascending;'
  write ( *, '(a)' ) '   2: strictly ascending;'
  write ( *, '(a)' ) '   3: descending;'
  write ( *, '(a)' ) '   4: strictly descending.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    x(1:n) = x_test(1:n,test)

    call r8vec_order_type ( n, x, order )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The following vector has order type ', order
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(i8,g14.6)' ) j, x(j)
    end do

  end do

  return
end
subroutine test132 ( )

!*****************************************************************************80
!
!! TEST132 tests R8VEC_PERMUTE.
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

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( n ) :: perm = (/ 2, 4, 5, 1, 3 /)
  real ( kind = 8 ), dimension (n) :: x = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST132'
  write ( *, '(a)' ) '  R8VEC_PERMUTE permutes a R8VEC in place.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, Perm(I), X(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,2i8,g14.6)' ) i, perm(i), x(i)
  end do

  call r8vec_permute ( n, perm, x )

  call r8vec_print ( n, x, '  Permuted array:' )

  return
end
subroutine test133 ( )

!*****************************************************************************80
!
!! TEST133 tests R8VEC_POLARIZE.
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

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension ( n ) :: a = (/ &
    1.0D+00, 2.0D+00,  3.0D+00 /)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a_normal(n)
  real ( kind = 8 ) a_parallel(n)
  real ( kind = 8 ) ap_norm
  real ( kind = 8 ), dimension ( n ) :: p = (/ &
    3.0D+00, 1.0D+00, -2.0D+00 /)
  real ( kind = 8 ) p_norm
  real ( kind = 8 ) pan
  real ( kind = 8 ) pap

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST133'
  write ( *, '(a)' ) '  R8VEC_POLARIZE decomposes a vector into'
  write ( *, '(a)' ) '  components parallel and normal to a direction.'

  call r8vec_print ( n, a, '  Original vector:' )

  call r8vec_print ( n, p, '  Direction vector:' )

  call r8vec_polarize ( n, a, p, a_normal, a_parallel )

  call r8vec_print ( n, a_normal, '  Normal component:' )

  call r8vec_print ( n, a_parallel, '  Parallel component:' )

  pan = dot_product ( p(1:n), a_normal(1:n) )

  p_norm = sqrt ( sum ( p(1:n)**2 ) )
  ap_norm = sqrt ( sum ( a_parallel(1:n)**2 ) )

  pap = dot_product ( p(1:n), a_parallel(1:n) ) / ( p_norm * ap_norm )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) &
    '  Dot product of P and A_normal (should be 0) ', pan
  write ( *, '(a,g14.6)' ) &
    '  Cosine of angle between P and A_parallel (should be 1 or -1) ', pap

  a2(1:n) = a_normal(1:n) + a_parallel(1:n)

  call r8vec_print ( n, a2, '  Sum of components (should equal A):' )

  return
end
subroutine test134 ( )

!*****************************************************************************80
!
!! TEST134 tests R8VEC_ROTATE.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ), dimension ( n ) :: a = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)
  integer ( kind = 4 ) m

  m = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST134'
  write ( *, '(a)' ) '  R8VEC_ROTATE rotates a R8VEC in place.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Rotate entries ', m, ' places to the right.'

  call r8vec_print ( n, a, '  Original array:' )

  call r8vec_rotate ( n, a, m )

  call r8vec_print ( n, a, '  Rotated array:' )

  return
end
subroutine test135 ( )

!*****************************************************************************80
!
!! TEST135 tests R8VEC_REVERSE.
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

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST135'
  write ( *, '(a)' ) '  R8VEC_REVERSE reverses a R8VEC.'

  call r8vec_indicator ( n, a )

  call r8vec_print ( n, a, '  Original array:' )

  call r8vec_reverse ( n, a )

  call r8vec_print ( n, a, '  Reversed array:' )

  a(1:n) = a(n:1:-1)

  call r8vec_print ( n, a, '  Re-reversed array using a(1:n) = a(n:1:-1)' )

  return
end
subroutine test136 ( )

!*****************************************************************************80
!
!! TEST136 tests R8VEC_SEARCH_BINARY_A;
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

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) index
  integer ( kind = 4 ) nc
  real ( kind = 8 ) search_val
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST136'
  write ( *, '(a)' ) '  For ascending order:'
  write ( *, '(a)' ) '  R8VEC_SEARCH_BINARY_A searches a sorted array;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r8vec_uniform_01 ( n, seed, a )

  search_val = a(1)

  call r8vec_sort_heap_a ( n, a )

  call r8vec_print ( n, a, '  Sorted vector A:' )
!
!  Now search the sorted array for a given value.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Search the array for the value ', search_val

  call r8vec_search_binary_a ( n, a, search_val, index )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEARCH RESULT:'
  write ( *, '(a)' ) ' '

  if ( 0 < index ) then
    write ( *, '(a,i8)' ) '    The value occurs in index ', index
  else
    write ( *, '(a)' ) '    The value does not occur in the array.'
  end if

  return
end
subroutine test137 ( )

!*****************************************************************************80
!
!! TEST137 tests R8VEC_SORT_BUBBLE_A.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST137'
  write ( *, '(a)' ) '  R8VEC_SORT_BUBBLE_A ascending sorts a R8VEC.'

  b = 0.0D+00
  c = 3.0D+00 * real ( n, kind = 8 )
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print_some ( n, a, 1, 10, '  Original array:' )

  call r8vec_sort_bubble_a ( n, a )

  call r8vec_print_some ( n, a, 1, 10, '  Ascending sorted array:' )

  return
end
subroutine test138 ( )

!*****************************************************************************80
!
!! TEST138 tests R8VEC_SORT_HEAP_A.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST138'
  write ( *, '(a)' ) '  R8VEC_SORT_HEAP_A ascending sorts a R8VEC.'

  b = 0.0D+00
  c = 3.0D+00 * real ( n, kind = 8 )
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print_some ( n, a, 1, 10, '  Original array:' )

  call r8vec_sort_heap_a ( n, a )

  call r8vec_print_some ( n, a, 1, 10, '  Ascending sorted array:' )

  return
end
subroutine test139 ( )

!*****************************************************************************80
!
!! TEST139 tests R8VEC_SORT_HEAP_D.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST139'
  write ( *, '(a)' ) '  R8VEC_SORT_HEAP_D descending sorts a R8VEC.'

  b = 0.0D+00
  c = 3.0D+00 * real ( n, kind = 8 )
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print_some ( n, a, 1, 10, '  Original array:' )

  call r8vec_sort_heap_d ( n, a )

  call r8vec_print_some ( n, a, 1, 10, '  Descending sorted array:' )

  return
end
subroutine test140 ( )

!*****************************************************************************80
!
!! TEST140 tests R8VEC_SORT_HEAP_INDEX_A and R8VEC_SORT_HEAP_INDEX_D.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST140'
  write ( *, '(a)' ) '  R8VEC_SORT_HEAP_INDEX_A creates an ascending'
  write ( *, '(a)' ) '  sort index for a R8VEC.'
  write ( *, '(a)' ) '  R8VEC_SORT_HEAP_INDEX_D creates a descending'
  write ( *, '(a)' ) '  sort index for a R8VEC.'

  b = 0.0D+00
  c = 3.0D+00 * real ( n, kind = 8 )
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print_some ( n, a, 1, 10, '  Unsorted array:' )

  call r8vec_sort_heap_index_a ( n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After indexed ascending sort:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), A(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i8,g14.6)' ) i, indx(i), a(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now use the index array to carry out the'
  write ( *, '(a)' ) '  permutation implicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INDX(I), A(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) indx(i), a(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call R8VEC_PERMUTE to carry out the permutation'
  write ( *, '(a)' ) '  explicitly.'

  call r8vec_permute ( n, indx, a )

  call r8vec_print ( n, a, '  I, A(I)' )

  call r8vec_sort_heap_index_d ( n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After indexed descending sort:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), A(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i8,g14.6)' ) i, indx(i), a(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now use the index array to carry out the'
  write ( *, '(a)' ) '  permutation implicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INDX(I), ARRAY(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) indx(i), a(indx(i))
  end do

  return
end
subroutine test141 ( )

!*****************************************************************************80
!
!! TEST141 tests R8VEC_SORT_HEAP_MASK_A.
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

  integer ( kind = 4 ), parameter :: mask_num = 10
  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(mask_num)
  integer ( kind = 4 ), dimension ( mask_num ) :: mask = (/ &
    2, 4, 7, 8, 9, 12, 13, 16, 18, 19 /)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST141'
  write ( *, '(a)' ) '  R8VEC_SORT_HEAP_MASK_A creates an ascending'
  write ( *, '(a)' ) '  sort index for a masked R8VEC.'

  b = 0.0D+00
  c = 3.0D+00 * real ( n, kind = 8 )
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print ( n, a, '  Unsorted array:' )

  call i4vec_print ( mask_num, mask, '  The mask array:' )

  call r8vec_mask_print ( n, a, mask_num, mask, '  The masked unsorted array:' )

  call r8vec_sort_heap_mask_a ( n, a, mask_num, mask, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After masked indexed ascending sort:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), MASK(INDX(I)), A(MASK(INDX(I)))'
  write ( *, '(a)' ) ' '
  do i = 1, mask_num
    write ( *, '(3i8,g14.6)' ) i, indx(i), mask(indx(i)), a(mask(indx(i)))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call I4VEC_PERMUTE to carry out the index permutation'
  write ( *, '(a)' ) '  explicitly on the MASK vector.'

  call i4vec_permute ( mask_num, indx, mask )
!
!  Essentially, INDX becomes the identity vector now.
!
  call i4vec_indicator ( mask_num, indx )

  call i4vec_print ( mask_num, mask, '  The reordered mask array:' )

  call r8vec_mask_print ( n, a, mask_num, mask, &
    '  The reordered masked sorted array:' )

  return
end
subroutine test142 ( )

!*****************************************************************************80
!
!! TEST142 tests R8VEC_SORT_INSERT_A.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST142'
  write ( *, '(a)' ) '  R8VEC_SORT_INSERT_A ascending sorts a R8VEC.'

  b = 0.0D+00
  c = 3.0D+00 * real ( n, kind = 8 )
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print_some ( n, a, 1, 10, '  Unsorted array:' )

  call r8vec_sort_insert_a ( n, a )

  call r8vec_print_some ( n, a, 1, 10, '  Sorted array:' )

  return
end
subroutine test143 ( )

!*****************************************************************************80
!
!! TEST143 tests R8VEC_SORT_INSERT_INDEX_A.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST143'
  write ( *, '(a)' ) '  R8VEC_SORT_INSERT_INDEX_A creates an ascending'
  write ( *, '(a)' ) '  sort index for a R8VEC.'

  b = 0.0D+00
  c = 3.0D+00 * real ( n, kind = 8 )
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print_some ( n, a, 1, 10, '  Unsorted array:' )

  call r8vec_sort_insert_index_a ( n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After indexed ascending sort:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), A(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i8,g14.6)' ) i, indx(i), a(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now use the index array to carry out the'
  write ( *, '(a)' ) '  permutation implicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), A(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8, i8,g14.6)' ) i, indx(i), a(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call R8VEC_PERMUTE to carry out the permutation'
  write ( *, '(a)' ) '  explicitly.'

  call r8vec_permute ( n, indx, a )

  call r8vec_print_some ( n, a, 1, 10, '  Permuted data' )

  return
end
subroutine test144 ( )

!*****************************************************************************80
!
!! TEST144 tests R8VEC_SORT_QUICK_A.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST144'
  write ( *, '(a)' ) '  R8VEC_SORT_QUICK_A sorts a R8VEC'
  write ( *, '(a)' ) '  using quick sort.'

  b = 0.0D+00
  c = 10.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  call r8vec_print ( n, a, '  Unsorted array:' )

  call r8vec_sort_quick_a ( n, a )

  call r8vec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test145 ( )

!*****************************************************************************80
!
!! TEST145 tests R8VEC_SORTED_MERGE_A;
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

  integer ( kind = 4 ), parameter :: na = 10
  integer ( kind = 4 ), parameter :: nb = 10

  real ( kind = 8 ) a(na)
  real ( kind = 8 ) b(nb)
  real ( kind = 8 ) c(na+nb)
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST145'
  write ( *, '(a)' ) '  For ascending order:'
  write ( *, '(a)' ) '  R8VEC_SORTED_MERGE_A merges two sorted R8VEC''s;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r8vec_uniform_01 ( na, seed, a )
  call r8vec_uniform_01 ( nb, seed, b )

  call r8vec_sort_heap_a ( na, a )

  call r8vec_sort_heap_a ( nb, b )

  call r8vec_print ( na, a, '  Sorted vector A:' )

  call r8vec_print ( nb, b, '  Sorted vector B:' )

  call r8vec_sorted_merge_a ( na, a, nb, b, nc, c )

  call r8vec_print ( nc, c, '  Merged vector C:' )

  return
end
subroutine test146 ( )

!*****************************************************************************80
!
!! TEST146 tests R8VEC_SORTED_NEAREST.
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

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_uniform_ab
  integer ( kind = 4 ) r8vec_sorted_nearest
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST146'
  write ( *, '(a)' ) '  R8VEC_SORTED_NEAREST finds the nearest entry'
  write ( *, '(a)' ) '  in a sorted R8VEC.'

  b = 0.0D+00
  c = 10.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, x )
  call r8vec_sort_heap_a ( n, x )

  call r8vec_print ( n, x, '  Sorted array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Test        Nearest'
  write ( *, '(a)' ) '     Value    Index   Value'
  write ( *, '(a)' ) ' '
  do i = 1, 10

    xval = r8_uniform_ab ( b, c, seed )

    j = r8vec_sorted_nearest ( n, x, xval )

    write ( *, '(2x,f8.4,4x,i8,2x,f8.4)' ) xval, j, x(j)

  end do

  return
end
subroutine test1465 ( )

!*****************************************************************************80
!
!! TEST1465 tests R8VEC_SORTED_RANGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) r_lo
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1465'
  write ( *, '(a)' ) '  R8VEC_SORTED_RANGE seeks the range of indices'
  write ( *, '(a)' ) '  in a sorted vector R so that'
  write ( *, '(a)' ) '  R_LO <= R(I_LO:I_HI) <= R_HI.'

  seed = 123456789

  do test = 1, 5

    call r8vec_uniform_01 ( n, seed, r )

    call r8vec_sort_heap_a ( n, r )

    call r8vec_print ( n, r, '  Sorted array R:' )

    r_lo = r8_uniform_01 ( seed )
    r_hi = r8_uniform_01 ( seed )

    if ( r_hi < r_lo ) then
      t = r_lo
      r_lo = r_hi
      r_hi = t
    end if

    call r8vec_sorted_range ( n, r, r_lo, r_hi, i_lo, i_hi )

    write ( *, '(a)' ) ' '
    if ( i_hi < i_lo ) then
      write ( *, '(2x,a4,2x,g14.6)' ) 'R_LO', r_lo
      write ( *, '(2x,a4,2x,g14.6)' ) 'R_HI', r_hi
      write ( *, '(2x,a)' ) '  Empty range in R.'
    else

      write ( *, '(2x,a4,2x,g14.6)' ) 'R_LO', r_lo
      do i = i_lo, i_hi
        write ( *, '(2x,i4,2x,g14.6)' ) i, r(i)
      end do
      write ( *, '(2x,a4,2x,g14.6)' ) 'R_HI', r_hi
    end if

  end do

  return
end
subroutine test147 ( )

!*****************************************************************************80
!
!! TEST147 tests R8VEC_SORTED_SPLIT and R8VEC_SPLIT.
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

  integer ( kind = 4 ), parameter :: n = 25

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i_gt
  integer ( kind = 4 ) i_lt
  integer ( kind = 4 ) isplit
  integer ( kind = 4 ) seed
  real ( kind = 8 ) split

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST147'
  write ( *, '(a)' ) '  R8VEC_SORTED_SPLIT splits a sorted vector into'
  write ( *, '(a)' ) '  entries less than and greater than a'
  write ( *, '(a)' ) '  splitting value.'
  write ( *, '(a)' ) '  R8VEC_SPLIT splits an unsorted vector'
  write ( *, '(a)' ) '  in the same way.'

  b = 0.0D+00
  c = 10.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  a(1:n) = real ( nint ( a(1:n) ), kind = 8 ) / 2.0D+00

  call r8vec_sort_heap_a ( n, a )

  split = 0.5D+00 * ( a(1) + a(n) )

  call r8vec_print ( n, a, '  The sorted array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Splitting value is ', split

  call r8vec_sorted_split ( n, a, split, i_lt, i_gt )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Lower index I_LT = ', i_lt
  write ( *, '(a,i8)' ) '  Upper index I_GT = ', i_gt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now repeat test with R8VEC_SPLIT.'

  call r8vec_permute_uniform ( n, a, seed )

  call r8vec_print ( n, a, '  The shuffled array:' )

  call r8vec_split ( n, a, split, isplit )

  call r8vec_print ( n, a, '  The split array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Array entries <= SPLIT up to index ', isplit

  return
end
subroutine test1475 ( )

!*****************************************************************************80
!
!! R8LIB_TEST1475 tests R8VEC_SORTED_UNDEX.
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
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) tol
  integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
  integer ( kind = 4 ) x_unique_num
  real ( kind = 8 ), dimension ( x_num ) :: x_val = (/ &
    11.0, 11.0, 11.0, 22.0, 22.0, 33.0, 33.0, 55.0, 55.0 /)
  integer ( kind = 4 ) xdnu(x_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: xu_val

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8LIB_TEST1475'
  write ( *, '(a)' ) '  R8VEC_SORTED_UNDEX produces index vectors which create a sorted'
  write ( *, '(a)' ) '  list of the unique elements of a sorted R8VEC,'
  write ( *, '(a)' ) '  and a map from the original vector to the (implicit)'
  write ( *, '(a)' ) '  vector of sorted unique elements.'

  call r8vec_print ( x_num, x_val, '  The vector X:' )

  tol = r8_epsilon ( )
  call r8vec_sorted_unique_count ( x_num, x_val, tol, x_unique_num )

  allocate ( undx(1:x_unique_num) )
  allocate ( xu_val(1:x_unique_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Tolerance for equality is ', tol
  write ( *, '(a,i8)' ) '  Number of unique entries in X is ', x_unique_num

  call r8vec_sorted_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to list the unique elements of X'
  write ( *, '(a)' ) '  in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX   X(UNDX)'
  write ( *, '(a)' ) ' '
  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,f8.1)' ) i, undx(i), x_val(undx(i))
  end do

  xu_val(1:x_unique_num) = x_val(undx(1:x_unique_num))

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to created XU, a copy of X'
  write ( *, '(a)' ) '  containing only the unique elements, in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX     XU(I)'
  write ( *, '(a)' ) ' '
  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,f8.1)' ) i, undx(i), xu_val(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' )'  XDNU can be used to match each element of X with one of the'
  write ( *, '(a)' )'  unique elements'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' )'     I  XDNU    X(I)       XU(XDNU(I))'
  write ( *, '(a)' ) ' '

  do i = 1, x_num
    write ( *, '(2x,i4,2x,i4,2x,f8.1,2x,f12.1)' ) &
      i, xdnu(i), x_val(i), xu_val(xdnu(i))
  end do

  deallocate ( undx )
  deallocate ( xu_val )

  return
end
subroutine test148 ( )

!*****************************************************************************80
!
!! TEST148 tests R8VEC_SORTED_UNIQUE;
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

  integer ( kind = 4 ), parameter :: n = 30

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ), parameter :: tol = 0.25D+00
  integer ( kind = 4 ) unique_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST148'
  write ( *, '(a)' ) '  R8VEC_SORTED_UNIQUE finds the unique entries in'
  write ( *, '(a)' ) '  a sorted R8VEC;'

  b = 0.0D+00
  c = real ( n, kind = 8 )
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  a(1:n) = real ( int ( a(1:n) ), kind = 8 )

  call r8vec_print_some ( n, a, 1, 10, '  Unsorted array:' )

  call r8vec_sort_heap_a ( n, a )

  call r8vec_sorted_unique ( n, a, tol, unique_num )

  call r8vec_print ( unique_num, a, '  Unique entries' )

  return
end
subroutine test149 ( )

!*****************************************************************************80
!
!! TEST149 tests R8VEC_SORTED_UNIQUE_COUNT;
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

  integer ( kind = 4 ), parameter :: n = 30

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ), parameter :: tol = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST149'
  write ( *, '(a)' ) '  R8VEC_SORTED_UNIQUE_COUNT counts the unique entries'
  write ( *, '(a)' ) '  of a sorted R8VEC;'

  b = 0.0D+00
  c = real ( n, kind = 8 )
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a )

  a(1:n) = real ( int ( a(1:n) ), kind = 8 )

  call r8vec_sorted_unique_count ( n, a, tol, unique_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Using a tolerance of ', tol
  write ( *, '(a,i8,a)' ) '  R8VEC_SORTED_UNIQUE_COUNT counts ', unique_num, &
    ' unique entries in A.'

  return
end
subroutine test150 ( )

!*****************************************************************************80
!
!! TEST150 tests R8VEC_SORTED_UNIQUE_HIST.
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

  integer ( kind = 4 ), parameter :: unique_max = 30
  integer ( kind = 4 ), parameter :: n = 30

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) acount(unique_max)
  real ( kind = 8 ) auniq(unique_max)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed
  real ( kind = 8 ), parameter :: tol = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST150'
  write ( *, '(a)' ) '  R8VEC_SORTED_UNIQUE_HIST makes a historgram of'
  write ( *, '(a)' ) '  the unique entries in a real vector.'

  b = 0.0D+00
  c = real ( n, kind = 8 )
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Using random number seed ', seed

  call r8vec_uniform_ab ( n, b, c, seed, a )

  a(1:n) = real ( int ( a(1:n) ), kind = 8 ) + 0.5D+00

  call r8vec_print ( n, a, '  Unsorted array:' )

  call r8vec_sort_bubble_a ( n, a )

  call r8vec_print ( n, a, '  Ascending sorted array:' )

  call r8vec_sorted_unique_hist ( n, a, tol, unique_max, unique_num, &
    auniq, acount )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  R8VEC_SORTED_UNIQUE_HIST counts ' , unique_num, &
    ' unique entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Value  Multiplicity'
  write ( *, '(a)' ) ' '
  do i = 1, unique_num
    write ( *, '(2x,i8,2x,g14.6,2x,i8)' ) i, auniq(i), acount(i)
  end do

  return
end
subroutine test1504 ( )

!*****************************************************************************80
!
!! TEST1504 tests R8VEC_TRANSPOSE_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1504'
  write ( *, '(a)' ) '  R8VEC_TRANSPOSE_PRINT prints an R8VEC "tranposed",'
  write ( *, '(a)' ) '  that is, placing multiple entries on a line.'

  call r8vec_uniform_01 ( n, seed, x )

  call r8vec_transpose_print ( n, x, '  The vector X:' )

  return
end
subroutine test1505 ( )

!*****************************************************************************80
!
!! TEST1505 tests R8VEC_UNDEX.
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
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) tol
  integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
  integer ( kind = 4 ) x_unique_num
  real ( kind = 8 ), dimension ( x_num ) :: x_val = (/ &
    33.0, 55.0, 11.0, 11.0, 55.0, 33.0, 22.0, 22.0, 11.0 /)
  integer ( kind = 4 ) xdnu(x_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: xu_val

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1505'
  write ( *, '(a)' ) '  R8VEC_UNDEX produces index vectors which create a sorted'
  write ( *, '(a)' ) '  list of the unique elements of an (unsorted) R8VEC,'
  write ( *, '(a)' ) '  and a map from the original vector to the (implicit)'
  write ( *, '(a)' ) '  vector of sorted unique elements.'

  call r8vec_print ( x_num, x_val, '  The vector X:' )

  tol = r8_epsilon ( )
  call r8vec_unique_count ( x_num, x_val, tol, x_unique_num )

  allocate ( undx(1:x_unique_num) )
  allocate ( xu_val(1:x_unique_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Tolerance for equality is ', tol
  write ( *, '(a,i8)' ) '  Number of unique entries in X is ', x_unique_num

  call r8vec_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to list the unique elements of X'
  write ( *, '(a)' ) '  in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX   X(UNDX)'
  write ( *, '(a)' ) ' '
  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,f8.1)' ) i, undx(i), x_val(undx(i))
  end do

  xu_val(1:x_unique_num) = x_val(undx(1:x_unique_num))

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  UNDX can be used to created XU, a copy of X'
  write ( *, '(a)' ) '  containing only the unique elements, in sorted order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  UNDX     XU(I)'
  write ( *, '(a)' ) ' '
  do i = 1, x_unique_num
    write ( *, '(2x,i4,2x,i4,2x,f8.1)' ) i, undx(i), xu_val(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' )'  XDNU can be used to match each element of X with one of the'
  write ( *, '(a)' )'  unique elements'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' )'     I  XDNU    X(I)       XU(XDNU(I))'
  write ( *, '(a)' ) ' '

  do i = 1, x_num
    write ( *, '(2x,i4,2x,i4,2x,f8.1,2x,f12.1)' ) i, xdnu(i), x_val(i), xu_val(xdnu(i))
  end do

  deallocate ( undx )
  deallocate ( xu_val )

  return
end
subroutine test151 ( )

!*****************************************************************************80
!
!! TEST151 tests R8VEC_UNIFORM_AB.
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

  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ), parameter :: b = 10.0D+00
  real ( kind = 8 ), parameter :: c = 20.0D+00
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST151'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_AB returns a random R8VEC '
  write ( *, '(a)' ) '  with entries in a given range [ B, C ]'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this problem:'
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a,g14.6)' ) '  C = ', c

  seed = 123456789

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Input SEED = ', seed

    call r8vec_uniform_ab ( n, b, c, seed, r )

    call r8vec_print_some ( n, r, 1, 10, '  Random vector:' )

  end do

  return
end
subroutine test1515 ( )

!*****************************************************************************80
!
!! TEST1515 tests R8VEC_UNIFORM_01.
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

  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1515'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_01 returns a random R8VEC '
  write ( *, '(a)' ) '  with entries in [0,1].'

  seed = 123456789

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Input SEED = ', seed

    call r8vec_uniform_01 ( n, seed, r )

    call r8vec_print_some ( n, r, 1, 10, '  Random vector:' )

  end do

  return
end
subroutine test153 ( )

!*****************************************************************************80
!
!! TEST153 tests R8VEC2_SORT_A and R8VEC2_SORT_D.
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

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST153'
  write ( *, '(a)' ) '  For a pair of R8VEC''s:'
  write ( *, '(a)' ) '  R8VEC2_SORT_A ascending sorts;'
  write ( *, '(a)' ) '  R8VEC2_SORT_D descending sorts;'

  b = 1.0D+00
  c = 3.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, b,  c, seed, a1 )

  b = 5.0D+00
  c = 10.0D+00

  call r8vec_uniform_ab ( n, b, c, seed, a2 )

  a1(3) = a1(1)
  a2(3) = a2(1)

  a1(6) = a1(2)
  a2(6) = a2(2)

  a1(9) = a1(1)
  a2(9) = a2(1)

  call r8vec2_print ( n, a1, a2, '  The pair of arrays:' )

  call r8vec2_sort_a ( n, a1, a2 )

  call r8vec2_print ( n, a1, a2, '  Arrays after ascending sort:' )

  call r8vec2_sort_d ( n, a1, a2 )

  call r8vec2_print ( n, a1, a2, '  Arrays after descending sort:' )

  return
end
subroutine test154 ( )

!*****************************************************************************80
!
!! TEST154 tests R8VEC2_SORT_HEAP_INDEX_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform_ab
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST154'
  write ( *, '(a)' ) '  R8VEC2_SORT_HEAP_INDEX_A creates a sort index'
  write ( *, '(a)' ) '  for an (X,Y) array.'

  do i = 1, n

    x(i) = real ( i4_uniform_ab ( 0, n, seed ), kind = 8 ) / real ( n, kind = 8 )
    y(i) = real ( i4_uniform_ab ( 0, n, seed ), kind = 8 ) / real ( n, kind = 8 )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The unsorted array:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, X(I), Y(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,6x,2g14.6)' ) i, x(i), y(i)
  end do

  call r8vec2_sort_heap_index_a ( n, x, y, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After sorting:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), X(I), Y(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i8,2g14.6)' ) i, indx(i), x(i), y(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now use the index array to carry out the'
  write ( *, '(a)' ) '  permutation implicitly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, INDX(I), X(INDX(I)), Y(INDX(I))'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i8,2g14.6)' ) i, indx(i), x(indx(i)), y(indx(i))
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8VEC_PERMUTE carries out the permutation.'

  call r8vec_permute ( n, indx, x )
  call r8vec_permute ( n, indx, y )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, X(I), Y(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,6x,2g14.6)' ) i, x(i), y(i)
  end do

  return
end
subroutine test155 ( )

!*****************************************************************************80
!
!! TEST155 tests R8VEC2_SORTED_UNIQUE.
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

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST155'
  write ( *, '(a)' ) '  For a pair of R8VEC''s:'
  write ( *, '(a)' ) '  R8VEC2_SORTED_UNIQUE counts unique entries.'

  b = 1.0D+00
  c = 3.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a1 )

  b = 5.0D+00
  c = 10.0D+00

  call r8vec_uniform_ab ( n, b, c, seed, a2 )

  a1(3) = a1(1)
  a2(3) = a2(1)

  a1(6) = a1(2)
  a2(6) = a2(2)

  a1(9) = a1(1)
  a2(9) = a2(1)

  call r8vec2_print ( n, a1, a2, '  The pair of arrays:' )

  call r8vec2_sort_a ( n, a1, a2 )

  call r8vec2_print ( n, a1, a2, '  Arrays after ascending sort:' )

  call r8vec2_sorted_unique ( n, a1, a2, unique_num )

  call r8vec2_print ( unique_num, a1, a2, '  UNIQed array:' )

  return
end
subroutine test156 ( )

!*****************************************************************************80
!
!! TEST156 tests R8VEC2_SORTED_UNIQUE_INDEX.
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

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST156'
  write ( *, '(a)' ) '  For a pair of R8VEC''s:'
  write ( *, '(a)' ) '  R8VEC2_SORTED_UNIQUE_INDEX indexes unique entries.'

  b = 1.0D+00
  c = 3.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a1 )

  b = 5.0D+00
  c = 10.0D+00

  call r8vec_uniform_ab ( n, b, c, seed, a2 )

  a1(3) = a1(1)
  a2(3) = a2(1)

  a1(6) = a1(2)
  a2(6) = a2(2)

  a1(9) = a1(1)
  a2(9) = a2(1)

  call r8vec2_sort_a ( n, a1, a2 )

  call r8vec2_print ( n, a1, a2, '  Sorted arrays:' )

  call r8vec2_sorted_unique_index ( n, a1, a2, unique_num, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of unique elements is ', unique_num

  call i4vec_print ( unique_num, indx, '  Index of Unique Elements:' )

  call r8vec_index_order ( unique_num, a1, indx )
  call r8vec_index_order ( unique_num, a2, indx )

  call r8vec2_print ( unique_num, a1, a2, '  After Indexed Nonunique Deletion.' )

  return
end
subroutine test157 ( )

!*****************************************************************************80
!
!! TEST157 tests R8VEC2_SUM_MAX_INDEX.
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

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST157'
  write ( *, '(a)' ) '  For a pair of R8VEC''s:'
  write ( *, '(a)' ) '  R8VEC2_SUM_MAX_INDEX: index of the sum vector'
  write ( *, '(a)' ) '  with maximum value.'

  b = 0.0D+00
  c = 10.0D+00
  seed = 123456789

  call r8vec_uniform_ab ( n, b, c, seed, a1 )

  b = 0.0D+00
  c = 5.0D+00

  call r8vec_uniform_ab ( n, b, c, seed, a2 )

  call r8vec2_print ( n, a1, a2, '  The pair of vectors:' )

  call r8vec2_sum_max_index ( n, a1, a2, ival )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Index of maximum in A+B: ', ival

  return
end
subroutine test158 ( )

!*****************************************************************************80
!
!! TEST158 tests R8VECS_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: na = 15

  real ( kind = 8 ) :: a(na) = (/ &
    11.0D+00, 12.0D+00, 13.0D+00, &
    21.0D+00, 22.0D+00, &
    31.0D+00, 32.0D+00, 33.0D+00, 34.0D+00, 35.0D+00, 36.0D+00, 37.0D+00, &
    41.0D+00, 42.0D+00, &
    51.0D+00 /)

  integer ( kind = 4 ) :: nvec(m+1) = (/ 1, 4, 6, 13, 15, 16 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST158'
  write ( *, '(a)' ) '  R8VECS_PRINT prints a packed R8VEC.'

  call r8vecs_print ( m, nvec, na, a, '  Packed R8VEC:' )

  return
end

