program main

!*****************************************************************************80
!
!! MAIN is the main program for SUBPAK_PRB.
!
!  Discussion:
!
!    SUBPAK_PRB tests the SUBPAK library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUBPAK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SUBPAK library.'

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
  call test225 ( )
  call test23 ( )
  call test24 ( )
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
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUBPAK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ANGLE_SHIFT.
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

  real ( kind = 8 ) alpha
  real ( kind = 8 ) angle_hi
  real ( kind = 8 ) angle_lo
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  ANGLE_SHIFT shifts an angle by multiples of'
  write ( *, '(a)' ) '  2 Pi until it lies between BETA and BETA+2Pi.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA     GAMMA   BETA+2Pi'
  write ( *, '(a)' ) ' '

  angle_lo = -4.0D+00 * pi
  angle_hi = +4.0D+00 * pi

  seed = 123456789

  do test = 1, test_num

    alpha = r8_uniform ( angle_lo, angle_hi, seed )

    beta = r8_uniform ( angle_lo, angle_hi, seed )

    call angle_shift ( alpha, beta, gamma )

    write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1)' ) &
      alpha, beta, gamma, beta + 2.0D+00 * pi

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests ANGLE_SHIFT_DEG.
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

  real ( kind = 8 ) alpha
  real ( kind = 8 ) :: angle_hi = +720.0D+00
  real ( kind = 8 ) :: angle_lo = -720.0D+00
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  ANGLE_SHIFT_DEG shifts an angle by multiples of'
  write ( *, '(a)' ) '  360 degrees until it lies between BETA and BETA+360.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA     GAMMA   BETA+360'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    alpha = r8_uniform ( angle_lo, angle_hi, seed )

    beta = r8_uniform ( angle_lo, angle_hi, seed )

    call angle_shift_deg ( alpha, beta, gamma )

    write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1)' ) &
      alpha, beta, gamma, beta + 360.0D+00

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests ANGLE_TO_RGB.
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

  real ( kind = 8 ) angle
  real ( kind = 8 ) :: angle_hi = 360.0D+00
  real ( kind = 8 ) :: angle_lo =   0.0D+00
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  ANGLE_TO_RGB converts an angle into an RGB color.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ANGLE        R         G         B'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    angle = r8_uniform ( angle_lo, angle_hi, seed )

    call angle_to_rgb ( angle, r, g, b )

    write ( *, '(2x,f8.1,2x,f8.3,2x,f8.3,2x,f8.3)' ) angle, r, g, b

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests AXIS_LIMITS.
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

  integer ( kind = 4 ) ndivs
  integer ( kind = 4 ) nticks
  real ( kind = 8 ) pxdiv
  real ( kind = 8 ) pxmax
  real ( kind = 8 ) pxmin
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  xmin = 67.3D+00
  xmax = 114.7D+00
  ndivs = 6

  call axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  AXIS_LIMITS adjusts plot limits'
  write ( *, '(a)' ) '    to "nicer" values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Input XMIN =    ', xmin
  write ( *, '(a,g14.6)' ) '  Input XMAX =    ', xmax
  write ( *, '(a,i8)' ) '  Input NDIVS =   ', ndivs
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Output PXMIN =  ', pxmin
  write ( *, '(a,g14.6)' ) '  Output PXMAX =  ', pxmax
  write ( *, '(a,g14.6)' ) '  Output PXDIV =  ', pxdiv
  write ( *, '(a,i8)' ) '  Output NTICKS = ', nticks

  xmin = - 26.0D+00
  xmax = + 26.0D+00
  ndivs = 10

  call axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Input XMIN =    ', xmin
  write ( *, '(a,g14.6)' ) '  Input XMAX =    ', xmax
  write ( *, '(a,i8)' ) '  Input NDIVS =   ', ndivs
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Output PXMIN =  ', pxmin
  write ( *, '(a,g14.6)' ) '  Output PXMAX =  ', pxmax
  write ( *, '(a,g14.6)' ) '  Output PXDIV =  ', pxdiv
  write ( *, '(a,i8)' ) '  Output NTICKS = ', nticks

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests AXIS_LIMITS.
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

  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) ndivs
  integer ( kind = 4 ) nticks
  real ( kind = 8 ) pxdiv
  real ( kind = 8 ) pxmax
  real ( kind = 8 ) pxmin
  integer ( kind = 4 ) test
  real ( kind = 8 ) xmax
  real ( kind = 8 ), dimension ( test_num ) :: xmax_test = (/ &
    9.0D+00, 4.125D+00, 193.75D+00, 2000.250D+00, 12.0D+00 /)
  real ( kind = 8 ) xmin
  real ( kind = 8 ), dimension ( test_num ) :: xmin_test = (/ &
    1.0D+00, 1.003D+00, 101.25D+00, 2000.125D+00, -7.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  AXIS_LIMITS computes "nice" limits for a graph'
  write ( *, '(a)' ) '    that must include a given range.'

  ndivs = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  All tests use NDIVS = ', ndivs
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      XMIN        XMAX       PXMIN       PXMAX       PXDIV    NTICKS'
  write ( *, '(a)') ' '

  do test = 1, test_num

    xmin = xmin_test(test)
    xmax = xmax_test(test)

    call axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )

    write ( *, '(2x,5g12.4,i8)' ) xmin, xmax, pxmin, pxmax, pxdiv, nticks

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests BAR_CHECK, BAR_CODE, BAR_DIGIT_CODE_LEFT, BAR_DIGIT_CODE_RIGHT.
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

  character ( len = 113 ) bar
  integer ( kind = 4 ) check
  integer ( kind = 4 ) digit(12)
  character ( len = 7 ) codel
  character ( len = 7 ) coder
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  BAR_CHECK checks digits for a barcode;'
  write ( *, '(a)' ) '  BAR_CODE computes the barcode for a string of'
  write ( *, '(a)' ) '    11 digits;'
  write ( *, '(a)' ) '  BAR_DIGIT_CODE_LEFT returns the left digit code.'
  write ( *, '(a)' ) '  BAR_DIGIT_CODE_RIGHT returns the right digit code.'

  do i = 1, 11
    digit(i) = mod ( i-1, 10 )
  end do

  call bar_check ( digit, check )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The check digit is ', check

  digit(12) = check

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The left and right digit codes:'
  write ( *, '(a)' ) ' '
  do i = 0, 9
    call bar_digit_code_left ( i, codel )
    call bar_digit_code_right ( i, coder )
    write ( *, '(2x,i2,2x,a7,2x,a7)' ) i, codel, coder
  end do

  call bar_code ( digit, bar )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bar code:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,a)' ) bar(1:9)
  write ( *, '(4x,a)' ) bar(10:12)
  write ( *, '(4x,a)' ) bar(13:19)
  write ( *, '(4x,a)' ) bar(20:26)
  write ( *, '(4x,a)' ) bar(27:33)
  write ( *, '(4x,a)' ) bar(34:40)
  write ( *, '(4x,a)' ) bar(41:47)
  write ( *, '(4x,a)' ) bar(48:54)
  write ( *, '(4x,a)' ) bar(55:59)
  write ( *, '(4x,a)' ) bar(60:66)
  write ( *, '(4x,a)' ) bar(67:73)
  write ( *, '(4x,a)' ) bar(74:80)
  write ( *, '(4x,a)' ) bar(81:87)
  write ( *, '(4x,a)' ) bar(88:94)
  write ( *, '(4x,a)' ) bar(95:101)
  write ( *, '(4x,a)' ) bar(102:104)
  write ( *, '(4x,a)' ) bar(105:113)

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests BMI_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bmi
  real ( kind = 8 ) bmi_english
  real ( kind = 8 ) c
  real ( kind = 8 ) r8_uniform
  real ( kind = 8 ) h
  real ( kind = 8 ) h_ft
  real ( kind = 8 ) h_in
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  real ( kind = 8 ) w

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  BMI_ENGLISH computes the Body Mass Index'
  write ( *, '(a)' ) '  given body measurements in English Units.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Weight         Height            BMI'
  write ( *, '(a)' ) '       (LB)     (FT           IN)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    b = 100.0D+00
    c = 250.0D+00

    w = r8_uniform ( b, c, seed )

    b = 4.0D+00
    c = 6.75D+00

    h = r8_uniform ( b, c, seed )

    h_ft = int ( h )
    h_in = real ( nint ( 12.0D+00 * ( h - h_ft ) ), kind = 8 )

    bmi = bmi_english ( w, h_ft, h_in )
    write ( *, '(2x,4f10.2)' ) w, h_ft, h_in, bmi

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests FAC_DIV, FAC_GCD, FAC_LCM, FAC_MUL, FAC_TO_I4, and I4_TO_FAC.
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

  integer ( kind = 4 ), parameter :: prime_num = 5

  integer ( kind = 4 ) bot
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) npower1(prime_num)
  integer ( kind = 4 ) npower2(prime_num)
  integer ( kind = 4 ) npower3(prime_num)
  integer ( kind = 4 ) top

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  For products of prime factors:'
  write ( *, '(a)' ) '  FAC_DIV computes a quotient;'
  write ( *, '(a)' ) '  FAC_MUL multiplies;'
  write ( *, '(a)' ) '  FAC_LCM computes the LCM;'
  write ( *, '(a)' ) '  FAC_GCD computes the GCD;'
  write ( *, '(a)' ) '  I4_TO_FAC converts an integer;'
  write ( *, '(a)' ) '  FAC_TO_I4 converts to an integer.'
  write ( *, '(a)' ) '  FAC_TO_RAT converts to a ratio.'

  i1 = 720
  i2 = 42

  call i4_to_fac ( i1, prime_num, npower1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Representation of I1 = ', i1
  write ( *, '(a)' ) ' '

  call fac_print ( prime_num, npower1 )

  call i4_to_fac ( i2, prime_num, npower2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Representation of I2 = ', i2
  write ( *, '(a)' ) ' '

  call fac_print ( prime_num, npower2 )

  call fac_lcm ( prime_num, npower1, npower2, npower3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LCM of I1, I2:'
  write ( *, '(a)' ) ' '

  call fac_print ( prime_num, npower3 )

  call fac_gcd ( prime_num, npower1, npower2, npower3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GCD of I1, I2:'
  write ( *, '(a)' ) ' '

  call fac_print ( prime_num, npower3 )

  call fac_mul ( prime_num, npower1, npower2, npower3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Product of I1, I2:'
  write ( *, '(a)' ) ' '

  call fac_print ( prime_num, npower3 )

  call fac_div ( prime_num, npower2, npower1, npower3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Quotient of I2 / I1:'
  write ( *, '(a)' ) ' '

  call fac_print ( prime_num, npower3 )

  call fac_to_rat ( prime_num, npower3, top, bot )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Quotient as a rational: ', top, ' / ', bot

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests GAUSS_SUM
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

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension ( n ) :: amplitude = (/ &
    10.0D+00, 5.0D+00, -3.0D+00 /)
  real ( kind = 8 ), dimension(dim_num,n) :: center = reshape ( (/ &
    2.0D+00, 3.0D+00, &
    5.0D+00, 8.0D+00, &
    7.0D+00, 5.0D+00 /), &
    (/ dim_num, n /) )
  real ( kind = 8 ) gauss_sum
  real ( kind = 8 ) gxy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension(n) :: width = (/ 1.0D+00, 2.0D+00, 4.0D+00 /)
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  GAUSS_SUM evaluates a function which is the sum of'
  write ( *, '(a)' ) '  Gaussian functions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of component Gaussians = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          Center    Amplitude  Width'
  write ( *, '(a)' ) '        X       Y'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(2x,i2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2)' ) &
      j, center(1:2,j), amplitude(j), width(j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X       Y        Gauss_Sum(X,Y)'
  write ( *, '(a)' ) ' '

  do i = 0, 10
    x(1) = real ( i, kind = 8 )
    do j = 0, 10
      x(2) = real ( j, kind = 8 )
      gxy = gauss_sum ( dim_num, n, amplitude, center, width, x )
      write ( *, '(2x,f6.2,2x,f6.2,2x,g14.6)' ) x(1), x(2), gxy
    end do
  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests GET_SEED.
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

  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_0
  integer ( kind = 4 ) seed_1
  integer ( kind = 4 ) seed_2
  integer ( kind = 4 ) seed_3
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  GET_SEED gets a seed for the random number'
  write ( *, '(a)' ) '  generator.  These values are computed from'
  write ( *, '(a)' ) '  the time and date.  Values computed nearby'
  write ( *, '(a)' ) '  in time will be near to each other, and'
  write ( *, '(a)' ) '  should be passed through a random number'
  write ( *, '(a)' ) '  generator a few times before use.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I      R(I)  R2(I)        R3(I)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call get_seed ( seed )
    seed_0 = seed
    x = r8_uniform_01 ( seed )
    seed_1 = seed
    x = r8_uniform_01 ( seed )
    seed_2 = seed
    x = r8_uniform_01 ( seed )
    seed_3 = seed
    write ( *, '(2x,4i12)' ) seed_0, seed_1, seed_2, seed_3
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests GRID1.
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

  integer ( kind = 4 ), parameter :: dim_num = 5
  integer ( kind = 4 ), parameter :: nstep = 11

  real ( kind = 8 ) x(dim_num,nstep)
  real ( kind = 8 ), dimension ( dim_num ) :: x1 = (/ &
    1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: x2 = (/ &
    1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  GRID1 computes a 1D grid between'
  write ( *, '(a)' ) '  two DIM_NUM dimensional points X1 and X2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Here, we will use ', nstep, ' steps'
  write ( *, '(a)' ) '  going from '
  write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
  write ( *, '(a)' ) '  to'
  write ( *, '(2x,5g12.4)' ) x2(1:dim_num)

  call grid1 ( dim_num, nstep, x1, x2, x )

  call r8mat_transpose_print ( dim_num, nstep, x, '  The grid matrix:' )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests GRID1N.
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

  integer ( kind = 4 ), parameter :: dim_num = 5
  integer ( kind = 4 ), parameter :: nstep = 11

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: x1 = (/ &
    1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x2 = (/ &
    1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  GRID1N computes a 1D grid between'
  write ( *, '(a)' ) '  two DIM_NUM dimensional points X1 and X2,'
  write ( *, '(a)' ) '  one point at a time.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Here, we will use ', nstep, ' steps'
  write ( *, '(a)' ) '  going from '
  write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
  write ( *, '(a)' ) '  to'
  write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
  write ( *, '(a)') ' '

  do i = 1, nstep
    call grid1n ( i, dim_num, nstep, x1, x2, x )
    write ( *, '(2x,i3,5g12.4)' ) i, x(1:dim_num)
  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests GRID2.
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

  integer ( kind = 4 ), parameter :: dim_num = 5
  integer ( kind = 4 ), parameter :: nstep = 20

  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) x(dim_num,nstep)
  real ( kind = 8 ), dimension(dim_num) :: x1 = (/ &
    1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x2 = (/ &
    1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /)

  j1 = 3
  j2 = 13

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  GRID2 computes a 1 D grid between'
  write ( *, '(a)' ) '  two DIM_NUM dimensional points X1 and X2,'
  write ( *, '(a)' ) '  computing X1 and X2 at user specified times.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Here, we will use ', nstep, ' steps,'
  write ( *, '(a,i8,a)' ) '  and on step ', j1, ' we will compute'
  write ( *,'(2x,5g12.4)' ) x1(1:dim_num)
  write ( *, '(a,i8,a)' ) '  and on step ', j2, ' we will compute'
  write ( *,'(2x,5g12.4)' ) x2(1:dim_num)

  call grid2 ( j1, j2, dim_num, nstep, x1, x2, x )

  call r8mat_print ( dim_num, nstep, x, '  The grid matrix:' )

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests GRID2N.
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

  integer ( kind = 4 ), parameter :: dim_num = 5

  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: x1 = (/ &
    1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x2 = (/ &
    1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /)

  j1 = 3
  j2 = 13

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  GRID2N computes points from a 1D grid'
  write ( *, '(a)' ) '  between two DIM_NUM dimensional points'
  write ( *, '(a)' ) '  X1 and X2, one at a time, with X1 and X2'
  write ( *, '(a)' ) '  having user specified J coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Here, on step ', j1, ' we would compute'
  write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
  write ( *, '(a,i8,a)' ) '  and on step ', j2, ' we would compute'
  write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
  write ( *, '(a)' ) ' '

  do j = 1, 20
    call grid2n ( j, j1, j2, dim_num, x1, x2, x )
    write ( *, '(2x,i3,5g12.4)' ) j, x(1:dim_num)
  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests GRID3.
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

  integer ( kind = 4 ), parameter :: dim_num = 5
  integer ( kind = 4 ), parameter :: nstep1 = 3
  integer ( kind = 4 ), parameter :: nstep2 = 6

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(dim_num,nstep1,nstep2)
  real ( kind = 8 ), dimension(dim_num) :: x1 = (/ &
    1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x2 = (/ &
    1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x3 = (/ &
    1.0D+00,  5.0D+00,  0.0D+00,  0.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  GRID3 computes a 2D grid in the plane'
  write ( *, '(a)' ) '  containing the DIM_NUM-dimensional'
  write ( *, '(a)' ) '  points X1, X2 and X3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Here, we will use ', nstep1, ' steps'
  write ( *, '(a)' ) '  going from '
  write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
  write ( *, '(a)' ) '  to'
  write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
  write ( *, '(a,i8,a)' ) '  and ', nstep2,' steps going to '
  write ( *, '(2x,5g12.4)' ) x3(1:dim_num)

  call grid3 ( dim_num, nstep1, nstep2, x1, x2, x3, x )

  do i = 1, nstep1
    write ( *, '(a)' ) ' '
    do j = 1, nstep2

      write ( *, '(2x,i3,i3,5g12.4)' ) i, j, x(1:dim_num,i,j)

    end do
  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests GRID3N.
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

  integer ( kind = 4 ), parameter :: dim_num = 5
  integer ( kind = 4 ), parameter :: nstep1 = 3
  integer ( kind = 4 ), parameter :: nstep2 = 6

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: x1 = (/ &
    1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x2 = (/ &
    1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x3 = (/ &
    1.0D+00,  5.0D+00,  0.0D+00,  0.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  GRID3N computes a point from a 2D'
  write ( *, '(a)' ) '  grid in the plane containing the '
  write ( *, '(a)' ) '  DIM_NUM-dimensional points X1, X2 and X3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  We use ', nstep1, ' steps from '
  write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
  write ( *, '(a)' ) '  to'
  write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
  write ( *, '(a,i8,a)' ) '  and ', nstep2, ' steps going to '
  write ( *, '(2x,5g12.4)' ) x3(1:dim_num)

  do j = 1, nstep1
    write ( *, '(a)' ) ' '
    do k = 1, nstep2

      call grid3n ( j, k, dim_num, nstep1, nstep2, x1, x2, x3, x )
      write ( *, '(2x,i3,i3,5g12.4)' ) j, k, x(1:dim_num)

    end do
  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests GRID4.
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

  integer ( kind = 4 ), parameter :: dim_num = 5
  integer ( kind = 4 ), parameter :: nstep1 = 6
  integer ( kind = 4 ), parameter :: nstep2 = 10

  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) x(dim_num,nstep1,nstep2)
  real ( kind = 8 ), dimension(dim_num) :: x1 = (/ &
    1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x2 = (/ &
    1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x3 = (/ &
    1.0D+00,  5.0D+00,  0.0D+00,  0.0D+00, 3.0D+00 /)

  j1 = 2
  j2 = 5
  k1 = 3
  k2 = 9

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  GRID4 computes a 2D planar grid'
  write ( *, '(a)' ) '  containing the DIM_NUM-dimensional'
  write ( *, '(a)' ) '  points X1, X2 and X3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We compute the points on the following steps:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2i8)' ) '  X1 on step ', j1, k1
  write ( *, '(a,2i8)' ) '  X2 on step ', j2, k1
  write ( *, '(a,2i8)' ) '  X3 on step ', j1, k2
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  We use ', nstep1, ' steps in the J direction'
  write ( *, '(a,i8,a)' ) '  and ', nstep2, ' steps in the K direction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The points X1, X2 and X3 are:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g12.4)' ) x3(1:dim_num)

  call grid4 ( j1, j2, k1, k2, dim_num, nstep1, nstep2, x1, x2, x3, x )

  do j = 1, nstep1
    write ( *, '(a)' ) ' '
    do k = 1, nstep2

      write ( *, '(2x,i3,i3,5g12.4)' ) j, k, x(1:dim_num,j,k)

    end do
  end do

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests GRID4N.
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

  integer ( kind = 4 ), parameter :: dim_num = 5
  integer ( kind = 4 ), parameter :: nstep1 = 6
  integer ( kind = 4 ), parameter :: nstep2 = 10

  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: x1 = (/ &
    1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x2 = (/ &
    1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: x3 = (/ &
    1.0D+00,  5.0D+00,  0.0D+00,  0.0D+00, 3.0D+00 /)

  j1 = 2
  j2 = 5
  k1 = 3
  k2 = 9

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  GRID4N computes, one at a time, points'
  write ( *, '(a)' ) '  on a 2D grid in the plane containing'
  write ( *, '(a)' ) '  the DIM_NUM-dimensional points X1, X2 and X3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We wish to compute the points on the following'
  write ( *, '(a)' ) '  steps:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2i8)' ) '  X1 on step ', j1, k1
  write ( *, '(a,2i8)' ) '  X2 on step ', j2, k1
  write ( *, '(a,2i8)' ) '  X3 on step ', j1, k2
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  We use ', nstep1, ' steps in the J direction'
  write ( *, '(a,i8,a)' ) '  and ', nstep2, ' steps in the K direction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The points X1, X2 and X3 are:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g12.4)' ) x3(1:dim_num)

  do j = 1, nstep1
    write ( *, '(a)' ) ' '
    do k = 1, nstep2

      call grid4n ( j, j1, j2, k, k1, k2, dim_num, nstep1, nstep2, &
        x1, x2, x3, x )

      write ( *, '(2x,i3,i3,5g12.4)' ) j, k, x(1:dim_num)

    end do
  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests INDEX1_COL, INDEX1_ROW, and related functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) in(n_max)
  integer ( kind = 4 ) in_max(n_max)
  integer ( kind = 4 ) in_min(n_max)
  integer ( kind = 4 ) index1_col
  integer ( kind = 4 ) index1_row
  integer ( kind = 4 ) index2_col
  integer ( kind = 4 ) index2_row
  integer ( kind = 4 ) index3_col
  integer ( kind = 4 ) index3_row
  integer ( kind = 4 ) index4_col
  integer ( kind = 4 ) index4_row
  integer ( kind = 4 ) indexn_col
  integer ( kind = 4 ) indexn_row
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) k_min
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  INDEX1_COL column indexes a 1D array,'
  write ( *, '(a)' ) '  INDEX1_ROW row indexes a 1D array,'
  write ( *, '(a)' ) '  and there are several more versions of these functions.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  By COLS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Imin     I  Imax  Xmin Index'
  write ( *, '(a)' ) ' '

  i_min = 1
  i = 3
  i_max = 5
  index_min = 0

  value = index1_col ( i_min, i, i_max, index_min )
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX1_COL', index_min, value

  n = 1
  in_min(1) = 1
  in(1) = 3
  in_max(1) = 5
  index_min = 0
  value = indexn_col ( n, in_min, in, in_max, index_min )
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_COL', index_min, value

  i_min = 1
  i = 3
  i_max = 5
  j_min = 1
  j = 2
  j_max = 4
  index_min = 0
  value = index2_col ( i_min, i, i_max, j_min, j, j_max, index_min )
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX2_COL', index_min, value

  n = 2
  in_min(1) = 1
  in(1) = 3
  in_max(1) = 5
  in_min(2) = 1
  in(2) = 2
  in_max(2) = 4
  index_min = 0
  value = indexn_col ( n, in_min, in, in_max, index_min )
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_COL', index_min, value

  i_min = 1
  i = 3
  i_max = 5
  j_min = 1
  j = 2
  j_max = 4
  k_min = 1
  k = 1
  k_max = 3
  index_min = 0
  value = index3_col ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, &
   index_min )
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) k_min, k, k_max
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX3_COL', index_min, value

  n = 3
  in_min(1) = 1
  in(1) = 3
  in_max(1) = 5
  in_min(2) = 1
  in(2) = 2
  in_max(2) = 4
  in_min(3) = 1
  in(3) = 1
  in_max(3) = 3
  index_min = 0
  value = indexn_col ( n, in_min, in, in_max, index_min )
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_COL', index_min, value

  i_min = 1
  i = 3
  i_max = 5
  j_min = 1
  j = 2
  j_max = 4
  k_min = 1
  k = 1
  k_max = 3
  l_min = 1
  l = 2
  l_max = 2
  index_min = 0
  value = index4_col ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, &
    l_min, l, l_max, index_min )
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) k_min, k, k_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) l_min, l, l_max
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX4_COL', index_min, value

  n = 4
  in_min(1) = 1
  in(1) = 3
  in_max(1) = 5
  in_min(2) = 1
  in(2) = 2
  in_max(2) = 4
  in_min(3) = 1
  in(3) = 1
  in_max(3) = 3
  in_min(4) = 1
  in(4) = 2
  in_max(4) = 2
  index_min = 0
  value = indexn_col ( n, in_min, in, in_max, index_min )
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_COL', index_min, value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  By ROWS:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Imin     I  Imax  Xmin Index'
  write ( *, '(a)' ) ' '

  i_min = 1
  i = 3
  i_max = 5
  index_min = 0
  value = index1_row ( i_min, i, i_max, index_min )
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX1_ROW', index_min, value

  n = 1
  in_min(1) = 1
  in(1) = 3
  in_max(1) = 5
  index_min = 0
  value = indexn_row ( n, in_min, in, in_max, index_min )
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_ROW', index_min, value

  i_min = 1
  i = 3
  i_max = 5
  j_min = 1
  j = 2
  j_max = 4
  index_min = 0
  value = index2_row ( i_min, i, i_max, j_min, j, j_max, index_min )
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX2_ROW', index_min, value

  n = 2
  in_min(1) = 1
  in(1) = 3
  in_max(1) = 5
  in_min(2) = 1
  in(2) = 2
  in_max(2) = 4
  index_min = 0
  value = indexn_row ( n, in_min, in, in_max, index_min )
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_ROW', index_min, value

  i_min = 1
  i = 3
  i_max = 5
  j_min = 1
  j = 2
  j_max = 4
  k_min = 1
  k = 1
  k_max = 3
  index_min = 0
  value = index3_row ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, &
    index_min )
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) k_min, k, k_max
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX3_ROW', index_min, value

  n = 3
  in_min(1) = 1
  in(1) = 3
  in_max(1) = 5
  in_min(2) = 1
  in(2) = 2
  in_max(2) = 4
  in_min(3) = 1
  in(3) = 1
  in_max(3) = 3
  index_min = 0
  value = indexn_row ( n, in_min, in, in_max, index_min )
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_ROW', index_min, value

  i_min = 1
  i = 3
  i_max = 5
  j_min = 1
  j = 2
  j_max = 4
  k_min = 1
  k = 1
  k_max = 3
  l_min = 1
  l = 2
  l_max = 2
  index_min = 0
  value = index4_row ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, &
    l_min, l, l_max, index_min )
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) k_min, k, k_max
  write ( *, '(2x,i4,2x,i4,2x,i4)' ) l_min, l, l_max
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX4_ROW', index_min, value

  n = 4
  in_min(1) = 1
  in(1) = 3
  in_max(1) = 5
  in_min(2) = 1
  in(2) = 2
  in_max(2) = 4
  in_min(3) = 1
  in(3) = 1
  in_max(3) = 3
  in_min(4) = 1
  in(4) = 2
  in_max(4) = 2
  index_min = 0
  value = indexn_row ( n, in_min, in, in_max, index_min )
  write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_ROW', index_min, value

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests ISBN_CHECK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: test_num = 8

  integer   ( kind = 4 ) check
  character ( len = 20 ) isbn
  character ( len = 20 ), dimension ( test_num ) :: isbn_test = (/ &
    '0-8493-9640-9', &
    '0-201-54275-7', &
    '0-521-35796-9', &
    '0-07-034025-0', &
    '0-7493-9640-9', &
    '0-201-54275-X', &
    '0-521-X5796-9', &
    '0-37-034025-0' /)
  integer   ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  ISBN_CHECK checks ISBN''s.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A correct ISBN has a checksum of 0.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    isbn = isbn_test(test)

    call isbn_check ( isbn, check )

    write ( *, '(2x,a20,5x,i2)' ) isbn, check

  end do

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests ISBN_FILL.
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

  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) check
  character ( len = 20 ) isbn
  character ( len = 20*test_num ) :: isbn_test = &
    '0-?493-9640-9' // &
    '0-201-5427?-7' // &
    '0-521-35796-?' // &
    '?-07-034025-0' // &
    '0-07-05?489-2'
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  ISBN_FILL can fill in a single missing digit'
  write ( *, '(a)' ) '  in an ISBN.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    isbn(1:20) = isbn_test(1+(test-1)*20:test*20)

    call isbn_fill ( isbn )

    call isbn_check ( isbn, check )

    write ( *, '(2x,a20,5x,a20,5x,i2)' ) &
      isbn_test(1+(test-1)*20:test*20), isbn, check

  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests LCM_12N.
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

  integer ( kind = 4 ) lcm_12n
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  LCM_12N computes the least common multiple of the'
  write ( *, '(a)' ) '  integers 1 through N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N       LCM_12N ( N )'
  write ( *, '(a)' ) ' '
  do n = 1, 12
    write ( *, '(2x,i3,2x,i8)' ) n, lcm_12n ( n )
  end do

  return
end
subroutine test225 ( )

!*****************************************************************************80
!
!! TEST225 tests LMAT_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 20
  integer ( kind = 4 ), parameter :: n = 50

  logical a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST225'
  write ( *, '(a)' ) '  LMAT_PRINT prints a logical matrix.'

  do i = 1, m
    do j = 1, n
      a(i,j) = ( mod ( i, j ) == 0 )
    end do
  end do

  call lmat_print ( m, n, a, '  A(I,J) = I is divisible by J' )

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests LUHN_CHECK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) check_sum
  integer ( kind = 4 ), dimension ( test_num ) :: check_sum_test = (/ &
     6, &
    20, &
    40, &
    80  &
  /)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: digit
  integer ( kind = 4 ) digit_num
  integer ( kind = 4 ), dimension ( test_num ) :: digit_num_test = (/ &
     4, &
     4, &
     9, &
    15  &
  /)
  integer ( kind = 4 ), dimension ( 32 ) :: digit_test = (/ &
    1, 1, 1, 1, &
    8, 7, 6, 3, &
    4, 4, 6, 6, 6, 7, 6, 5, 1, &
    3, 7, 7, 9, 5, 6, 5, 7, 0, 9, 4, 4, 7, 2, 6 &
  /)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  LUHN_CHECK computes the Luhn checksum'
  write ( *, '(a)' ) '  for a string of digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A correct string has a checksum divisible by 10.'

  ihi = 0

  do test = 1, test_num

    digit_num = digit_num_test(test)

    ilo = ihi + 1
    ihi = ihi + digit_num

    allocate ( digit(1:digit_num) )
    digit(1:digit_num) = digit_test(ilo:ihi)

    call luhn_check ( digit_num, digit, check_sum )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Test number ', test
    write ( *, '(a,i8)' ) '  Number of digits = ', digit_num
    write ( *, '(a,20i1)' ) '  Digits = ', digit(1:digit_num)
    write ( *, '(a,i8)' ) '  Computed check sum = ', check_sum
    write ( *, '(a,i8)' ) '  Correct check sum =  ', check_sum_test(test)

    deallocate ( digit )

  end do

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests PERM_INVERSE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ), dimension ( n ) :: p = (/ 4, 3, 5, 1, 7, 6, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  PERM_INVERSE inverts a permutation in place;'

  call perm_print ( n, p, '  The original permutation:' )

  call perm_inverse ( n, p )

  call perm_print ( n, p, '  The inverted permutation:' )

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests PRIME_GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) p
  integer ( kind = 4 ) prime_ge

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  PRIME_GE returns the smallest prime number greater'
  write ( *, '(a)' ) '    than or equal to N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     P'
  write ( *, '(a)' ) ' '

  do n = 1, 20

    p = prime_ge ( n )
    write ( *, '(2i8)' ) n, p

  end do

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests RANDOM_INITIALIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  RANDOM_INITIALIZE can make up a seed for the FORTRAN90'
  write ( *, '(a)' ) '  random number generator RANDOM_NUMBER, or use a'
  write ( *, '(a)' ) '  single SEED value from the user.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling RANDOM_INITIALIZE with a zero input value of SEED'
  write ( *, '(a)' ) '  tells the routine to make up a seed.  And, at least for'
  write ( *, '(a)' ) '  calls a few milliseconds apart, the output SEED should'
  write ( *, '(a)' ) '  be different.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In any case, if RANDOM_NUMBER is restarted by calling'
  write ( *, '(a)' ) '  RANDOM_INITIALIZE with a nonzero input SEED, then'
  write ( *, '(a)' ) '  the random number sequence should repeat.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call RANDOM_INITIALIZE 10 times.'
  write ( *, '(a)' ) '  Use a 0 SEED the first time, then GET_SEED after than.'
  write ( *, '(a)' ) '  Also, get the first three real random values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    SEED         Random 1, 2, 3'
  write ( *, '(a)' ) ' '

  seed = 0

  do i = 1, 10

    call random_initialize ( seed )
    call random_number ( harvest = r1 )
    call random_number ( harvest = r2 )
    call random_number ( harvest = r3 )
    write ( *, '(i12,2x,3f12.8)' ) seed, r1, r2, r3

    call get_seed ( seed )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call RANDOM_INITIALIZE with SEED = 5, 95, 5, 95.'
  write ( *, '(a)' ) '  We promise the random numbers will repeat the second time.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    SEED         Random 1, 2, 3'
  write ( *, '(a)' ) ' '
  seed = 5
  do i = 1, 4
    call random_initialize ( seed )
    call random_number ( harvest = r1 )
    call random_number ( harvest = r2 )
    call random_number ( harvest = r3 )
    write ( *, '(2x,i12,2x,3f12.8)' ) seed, r1, r2, r3
    seed = 100 - seed
  end do

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests RANDOM_NUMBER.
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

  integer ( kind = 4 ), parameter :: n_max = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed_size
  integer ( kind = 4 ), allocatable, dimension ( : ) :: seed_vec
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_mean
  real ( kind = 8 ) x_min
  real ( kind = 8 ) x_var

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  RANDOM_NUMBER is an intrinsic FORTRAN90 routine'
  write ( *, '(a)' ) '  to computer uniform random numbers.'
!
!  Determine the size of the random number seed vector.
!
  call random_seed ( size = seed_size )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Random number seed size = ', seed_size
!
!  Allocate the seed vector.
!
  allocate ( seed_vec(1:seed_size) )

  do i = 1, seed_size
    seed_vec(i) = 1234567 * i
  end do
!
!  Set the random number seed.
!
  call random_seed ( put = seed_vec )
!
!  Test 1:
!  Simply call 5 times for 1 value, and print.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #1: Call 5 times, 1 value each time.'
  write ( *, '(a)' ) ' '

  n = 1
  do i = 1, 5
    call random_number ( harvest = x(1) )
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

  do i = 1, seed_size
    seed_vec(i) = 1234567 * i
  end do

! call random_seed ( put = seed_vec )

  n = 1
  do i = 1, 5
    call random_number ( harvest = x(1) )
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

  do i = 1, seed_size
    seed_vec(i) = 1234567 * i
  end do

  call random_seed ( put = seed_vec )

  n = 5
  call random_number ( harvest = x(1:5) )

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

  do i = 1, seed_size
    seed_vec(i) = 1234567 * i
  end do

  call random_seed ( put = seed_vec )

  n = 2
  call random_number ( harvest = x(1:2) )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do

  n = 1
  call random_number ( harvest = x(1:1) )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do

  n = 2
  call random_number ( harvest = x(1:2) )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do
!
!  Test 5:
!  Determine the minimum, maximum, mean and variance.
!
  n = n_max
  call random_number (  harvest = x(1:n) )
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
  write ( *, '(a,g14.6)' ) '  Expected average  ', 0.5D+00
  write ( *, '(a,g14.6)' ) '  Expected variance ', 1.0D+00 / 12.0D+00

  deallocate ( seed_vec )

  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests RANDOM_NUMBER and RANDOM_SEED.
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

  integer ( kind = 4 ) i
  real ( kind = 8 ) mean
  integer ( kind = 4 ), parameter :: n = 1000
  integer ( kind = 4 ) seed_size
  integer ( kind = 4 ), allocatable, dimension ( : ) :: seed_vec
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  RANDOM_SEED is a FORTRAN90 routine which sets or gets'
  write ( *, '(a)' ) '  the random number set.'
  write ( *, '(a)' ) '  RANDOM_NUMBER returns a uniformly distributed random '
  write ( *, '(a)' ) '  value between 0 and 1.'

  call random_seed ( size = seed_size )

  allocate ( seed_vec(1:seed_size) )

  call random_seed ( get = seed_vec )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Current random seed vector:'
  write ( *, '(a)' ) ' '
  do i = 1, seed_size
    write ( *, '(2x,i4,2x,i20)' ) i, seed_vec(i)
  end do

  call random_number ( x(1:n) )

  call r8vec_mean ( n, x, mean )

  call r8vec_variance ( n, x, variance )

  write ( *, '(a,i8)' ) '  Number of values computed was N = ', n
  write ( *, '(a,g14.6)' ) '  Average value was ', mean
  write ( *, '(a,g14.6)' ) '  Variance was ', variance

  deallocate ( seed_vec )

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests RAT_FACTOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: factor_max = 10

  integer ( kind = 4 ) factor(factor_max)
  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mleft
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) power(factor_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  RAT_FACTOR factors a rational value.'

  m = 13 * 7 * 9 * 2
  n = 12

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Rational value is ', m, '/', n

  call rat_factor ( m, n, factor_max, factor_num, factor, power, mleft, nleft )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Prime representation:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, FACTOR(I), POWER(I)'
  write ( *, '(a)' ) ' '

  if ( mleft /= 1 .or. nleft /= 1 ) then
    write ( *, '(i8,i8,a,i8,a)' ) 0, mleft, ' / ', nleft, &
      ' (UNFACTORED PORTION)'
  end if

  do i = 1, factor_num
    write ( *, '(2x,3i8)' ) i, factor(i), power(i)
  end do

  return
end
subroutine test30 ( )

!*****************************************************************************80
!
!! TEST30 tests ROOTS_TO_R8POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  ROOTS_TO_R8POLY computes the coefficients of'
  write ( *, '(a)' ) '  a polynomial from its roots.'
  write ( *, '(a)' ) '  R8POLY_PRINT prints a polynomial.'

  call r8vec_indicator ( n, x )

  call r8vec_print ( n, x, '  Roots:' )

  call roots_to_r8poly ( n, x, c )

  call r8poly_print ( n, c, '  The polynomial' )

  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests SORT_HEAP_EXTERNAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ),  parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  SORT_HEAP_EXTERNAL sorts objects externally.'

  indx = 0
  i = 0
  j = 0
  isgn = 0

  b = 1
  c = n
  seed = 123456789

  call i4vec_uniform ( n, b, c, seed, a )

  call i4vec_print ( n, a, '  Unsorted array:' )

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( indx < 0 ) then

      isgn = 1
      if ( a(i) <= a(j) ) then
        isgn = -1
      end if

    else if ( 0 < indx ) then

      call i4_swap ( a(i), a(j) )

    else

      exit

    end if

  end do

  call i4vec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests TVEC_EVEN, TVEC_EVEN2 and TVEC_EVEN3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxt = 5

  integer ( kind = 4 ) nt
  real ( kind = 8 ) t(maxt)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  For evenly spaced angles between 0 and 2*PI:'
  write ( *, '(a)' ) '  TVEC_EVEN'
  write ( *, '(a)' ) '  TVEC_EVEN2'
  write ( *, '(a)' ) '  TVEC_EVEN3'

  nt = 4

  call tvec_even ( nt, t )

  call r8vec_print ( nt, t, '  TVEC_EVEN' )

  nt = 4

  call tvec_even2 ( nt, t )

  call r8vec_print ( nt, t, '  TVEC_EVEN2' )

  nt = 4

  call tvec_even3 ( nt, t )

  call r8vec_print ( nt, t, '  TVEC_EVEN3' )

  return
end
subroutine test33 ( )

!*****************************************************************************80
!
!! TEST33 tests TVEC_EVEN_BRACKET, TVEC_EVEN_BRACKET2, TVEC_EVEN_BRACKET3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxt = 5

  integer ( kind = 4 ) nt
  real ( kind = 8 ) t(maxt)
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  For evenly spaced angles between THETA1 and THETA2:'
  write ( *, '(a)' ) '  TVEC_EVEN_BRACKET'
  write ( *, '(a)' ) '  TVEC_EVEN_BRACKET2.'
  write ( *, '(a)' ) '  TVEC_EVEN_BRACKET3.'

  nt = 4
  theta1 = 30.0D+00
  theta2 = 90.0D+00

  call tvec_even_bracket ( nt, theta1, theta2, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '    NT = ', nt
  write ( *, '(a,g14.6)' ) '    THETA1 = ', theta1
  write ( *, '(a,g14.6)' ) '    THETA2 = ', theta2

  call r8vec_print ( nt, t, '  TVEC_EVEN_BRACKET' )

  nt = 5

  call tvec_even_bracket2 ( nt, theta1, theta2, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '    NT = ', nt
  write ( *, '(a,g14.6)' ) '    THETA1 = ', theta1
  write ( *, '(a,g14.6)' ) '    THETA2 = ', theta2

  call r8vec_print ( nt, t, '  TVEC_EVEN_BRACKET2' )

  nt = 3

  call tvec_even_bracket3 ( nt, theta1, theta2, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '    NT = ', nt
  write ( *, '(a,g14.6)' ) '    THETA1 = ', theta1
  write ( *, '(a,g14.6)' ) '    THETA2 = ', theta2

  call r8vec_print ( nt, t, '  TVEC_EVEN_BRACKET3' )

  return
end
subroutine test34 ( )

!*****************************************************************************80
!
!! TEST34 tests UPC_CHECK_DIGIT.
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

  integer ( kind = 4 ), parameter :: test_num = 2

  integer ( kind = 4 ) c
  integer ( kind = 4 ) l
  integer ( kind = 4 ), dimension ( test_num ) :: l_test = (/ 72890, 12345 /)
  integer ( kind = 4 ) p
  integer ( kind = 4 ), dimension ( test_num ) :: p_test = (/ 0, 0 /)
  integer ( kind = 4 ) r
  integer ( kind = 4 ), dimension ( test_num ) :: r_test = (/ 00011, 67890 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  UPC_CHECK_DIGIT determines the check digit for a UPC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  P-LLLLL-RRRRR-C'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p = p_test(test)
    l = l_test(test)
    r = r_test(test)

    call upc_check_digit ( p, l, r, c )

    write ( *, '(2x,i1,''-'',i5.5,''-'',i5.5,''-'',i1 )' ) p, l, r, c

  end do

  return
end
subroutine test35 ( )

!*****************************************************************************80
!
!! TEST35 calls VERSINE_PULSE.
!
!  Modified:
!
!    20 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) amp
  integer ( kind = 4 ) i
  real ( kind = 8 ) t
  real ( kind = 8 ) ta
  real ( kind = 8 ) tb
  real ( kind = 8 ) v
  real ( kind = 8 ) versine_pulse
  real ( kind = 8 ) v1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  VERSINE_PULSE adds a versine pulse to a constant signal.'
  write ( *, '(a)' ) ' '

  ta = 2.0D+00
  tb = 4.0D+00
  v1 = 1.0D+00
  amp = 3.0D+00

  do i = 0, 100
    t = real ( i, kind = 8 ) / 10.0D+00
    v = versine_pulse ( t, ta, tb, v1, amp )
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, t, v
  end do

  return
end
