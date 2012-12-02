program main

!*****************************************************************************80
!
!! MAIN is the main program for POLPAK_PRB.
!
!  Discussion:
!
!    POLPAK_PRB calls the POLPAK test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLPAK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the POLPAK library.'

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
  call test0115 ( )
  call test013 ( )
  call test012 ( )
  call test014 ( )
  call test015 ( )
  call test016 ( )
  call test017 ( )
  call test0175 ( )
  call test018 ( )
  call test0185 ( )
  call test019 ( )

  call test020 ( )
  call test021 ( )
  call test0215 ( )
  call test0216 ( )
  call test0217 ( )
  call test0218 ( )
  call test024 ( )
  call test02405 ( )
  call test0241 ( )
  call test0242 ( )
  call test0243 ( )
  call test01155 ( )
  call test025 ( )
  call test0255 ( )
  call test026 ( )
  call test0265 ( )
  call test028 ( )
  call test027 ( )
  call test029 ( )

  call test031 ( )
  call test032 ( )
  call test033 ( )
  call test034 ( )
  call test036 ( )
  call test037 ( )
  call test038 ( )

  call test040 ( )
  call test041 ( )
  call test042 ( )
  call test0425 ( )
  call test0427 ( )
  call test043 ( )
  call test044 ( )
  call test045 ( )
  call test046 ( )
  call test047 ( )
  call test048 ( )
  call test049 ( )

  call test050 ( )
  call test0505 ( )
  call test051 ( )
  call test052 ( )
  call test054 ( )
  call test055 ( )
  call test0552 ( )
  call test059 ( )
  call test0595 ( )
  call test060 ( )
  call test057 ( )
  call test058 ( )

  call test061 ( )
  call test0615 ( )
  call test062 ( )
  call test0623 ( )
  call test0625 ( )
  call test063 ( )
  call test0635 ( )
  call test064 ( )
  call test065 ( )
  call test066 ( )
  call test0665 ( )
  call test0667 ( )
  call test067 ( )
  call test0675 ( )
  call test0676 ( )
  call test06765 ( )
  call test022 ( )
  call test068 ( )
  call test0685 ( )
  call test06855 ( )
  call test06856 ( )
  call test069 ( )
  call test0695 ( )
  call test0696 ( )
  call test0697 ( )

  call test070 ( )
  call test071 ( )
  call test072 ( )
  call test073 ( )
  call test074 ( )
  call test075 ( )
  call test076 ( )
  call test077 ( )
  call test0773 ( )
  call test0775 ( )
  call test078 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLPAK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests AGM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) agm
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  AGM computes the arithmetic geometric mean.'
  write ( *, '(a)' ) '  AGM_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a,a)' ) '      A           B          ', &
    '   AGM                       AGM                   Diff'
  write ( *, '(a,a)' ) '                             ', &
    '  (Tabulated)                AGM(A,B)'
  write ( *, '(a)' ) ' '
     
  n_data = 0

  do

    call agm_values ( n_data, a, b, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = agm ( a, b )

    write ( *, '(2x,f10.6,2x,f10.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      a, b, fx, fx2, abs ( fx - fx2 )

  end do
     
  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests AGUD and GUD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) agud
  real ( kind = 8 ) gamma
  real ( kind = 8 ) gud
  integer ( kind = 4 ) i
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  AGUD computes the inverse Gudermannian;'
  write ( *, '(a)' ) '  GUD computes the Gudermannian.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              GUD(X)       AGUD(GUD(X))'
  write ( *, '(a)' ) ' '

  do i = 0, 10
    x = 1.0D+00 + real ( i, kind = 8 ) / 5.0D+00
    gamma = gud ( x )
    x2 = agud ( gamma )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, gamma, x2
  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests ALIGN_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m_max = 10
  integer ( kind = 4 ), parameter :: n_max = 10

  integer ( kind = 4 ) align_enum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) table(0:m_max,0:n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  ALIGN_ENUM counts the number of possible'
  write ( *, '(a)' ) '  alignments of two biological sequences.'

  do i = 0, m_max
    do j = 0, n_max
      table(i,j) = align_enum ( i, j )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Alignment enumeration table:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,5i5,6i8)' ) ( i, i = 0, n_max )
  write ( *, '(a)' ) ' '
  do i = 0, m_max
    write ( *, '(2x,i2,5i5,6i8)' ) i, table(i,0:n_max)
  end do

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests ARC_COSINE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arc_cosine
  integer ( kind = 4 ) i
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  ARC_COSINE computes the inverse cosine'
  write ( *, '(a)' ) '  of a given value, and chops out of bound arguments.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X     ARC_COSINE(X)     COS(ARC_COSINE(X))'
  write ( *, '(a)' ) ' '

  do i = -5, 5
    x = 1.0D+00 + real ( i, kind = 8 ) / 5.0D+00
    a = arc_cosine ( x )
    x2 = cos ( a )
    write ( *, '(2x,3g14.6)' ) x, a, x2
  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests ATAN4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  real ( kind = 8 ) atan4
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension ( test_num ) :: x_test = (/ &
     1.0D+00,  1.0D+00, 0.0D+00, -1.0D+00, &
    -1.0D+00, -1.0D+00, 0.0D+00,  1.0D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), dimension ( test_num ) :: y_test = (/ &
    0.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
    0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  ATAN4 computes the arc-tangent given Y and X;'
  write ( *, '(a)' ) '  ATAN2 is the system version of this routine.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X             Y          ATAN2(Y,X)   ATAN4(Y,X)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    x = x_test(test)
    y = y_test(test)
    write ( *, '(2x,4g14.6)' ) x, y, atan2 ( y, x ), atan4 ( y, x )
  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests BELL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2(0:10)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  BELL computes Bell numbers.'
  write ( *, '(a)' ) '  BELL_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N  exact C(I)  computed C(I)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call bell_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call bell ( n, c2 )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2(n)

  end do
 
  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests BENFORD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) benford
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  BENFORD(I) is the Benford probability of the'
  write ( *, '(a)' ) '  initial digit sequence I.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I,  BENFORD(I)'
  write ( *, '(a)' ) ' '

  do i = 1, 9
    write ( *, '(2x,i2,2x,g14.6)' )  i, benford(i)
  end do
 
  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests BERNOULLI_NUMBER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c0
  real ( kind = 8 ) c1(0:30)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  BERNOULLI_NUMBER computes Bernoulli numbers;'
  write ( *, '(a)' ) '  BERNOULLI_NUMBER_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I      Exact     Bernoulli'
  write ( *, '(a)' ) ' '
  
  n_data = 0

  do

    call bernoulli_number_values ( n_data, n, c0 )

    if ( n_data == 0 ) then
      exit
    end if

    call bernoulli_number ( n, c1 )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, c0, c1(n)

  end do
 
  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests BERNOULLI_NUMBER2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c0
  real ( kind = 8 ) c1(0:30)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  BERNOULLI_NUMBER2 computes Bernoulli numbers;'
  write ( *, '(a)' ) '  BERNOULLI_NUMBER_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I      Exact     Bernoulli2'
  write ( *, '(a)' ) ' '
  
  n_data = 0

  do

    call bernoulli_number_values ( n_data, n, c0 )

    if ( n_data == 0 ) then
      exit
    end if

    call bernoulli_number2 ( n, c1 )
 
    write ( *, '(2x,i4,2g14.6)' ) n, c0, c1(n)

  end do
 
  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests BERNOULLI_NUMBER3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  BERNOULLI_NUMBER3 computes Bernoulli numbers.'
  write ( *, '(a)' ) '  BERNOULLI_NUMBER_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I      Exact     BERNOULLI3'
  write ( *, '(a)' ) ' '
  
  n_data = 0

  do

    call bernoulli_number_values ( n_data, n, c0 )

    if ( n_data == 0 ) then
      exit
    end if

    call bernoulli_number3 ( n, c1 )

    write ( *, '(2x,i4,2g14.6)' ) n, c0, c1

  end do
 
  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests BERNOULLI_POLY;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bx
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 15
  real ( kind = 8 ) x

  x = 0.2D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  BERNOULLI_POLY evaluates Bernoulli polynomials;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I          BX'
  write ( *, '(a)' ) ' '
 
  do i = 1, n
    call bernoulli_poly ( i, x, bx )
    write ( *, '(2x,i2,2x,g16.8)' ) i, bx
  end do
 
  return
end
subroutine test0115 ( )

!*****************************************************************************80
!
!! TEST0115 tests BERNOULLI_POLY2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bx
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 15
  real ( kind = 8 ) x

  x = 0.2D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0115'
  write ( *, '(a)' ) '  BERNOULLI_POLY2 evaluates Bernoulli polynomials. '
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I          BX'
  write ( *, '(a)' ) ' '
 
  do i = 1, n
    call bernoulli_poly2 ( i, x, bx )
    write ( *, '(2x,i2,2x,2g16.8)' ) i, bx
  end do
 
  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests BERNSTEIN_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bvec(0:10)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013:'
  write ( *, '(a)' ) '  BERNSTEIN_POLY evaluates the Bernstein polynomials.'
  write ( *, '(a)' ) '  BERNSTEIN_POLY_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N   K   X       Exact        B(N,K)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernstein_poly_values ( n_data, n, k, x, b )

    if ( n_data == 0 ) then
      exit
    end if

    call bernstein_poly ( n, x, bvec )

    write ( *, '(2x,i4,i4,f7.4,2g14.6)' ) n, k, x, b, bvec(k)

  end do

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests BETA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) beta
  real ( kind = 8 ) fxy
  real ( kind = 8 ) fxy2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012:'
  write ( *, '(a)' ) '  BETA evaluates the Beta function.'
  write ( *, '(a)' ) '  BETA_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X      Y        Exact F       BETA(X,Y)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_values ( n_data, x, y, fxy )

    if ( n_data == 0 ) then
      exit
    end if

    fxy2 = beta ( x, y )

    write ( *, '(2x,2f8.4,2g14.6)' ) x, y, fxy, fxy2

  end do

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests BPAB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bern(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  BPAB evaluates Bernstein polynomials.'
  write ( *, '(a)' ) ' '

  x = 0.3D+00
  a = 0.0D+00
  b = 1.0D+00
  call bpab ( n, a, b, x, bern )
 
  write ( *, '(a,i4)' ) '  The Bernstein polynomials of degree ', n
  write ( *, '(a,g14.6)' ) '  based on the interval from ', a
  write ( *, '(a,g14.6)' ) '  to ', b
  write ( *, '(a,g14.6)' ) '  evaluated at X = ', x
  write ( *, '(a)' ) ' '
 
  do i = 0, n
    write ( *, '(2x,i4,2x,g14.6)' )  i, bern(i)
  end do
 
  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests CARDAN and CARDAN_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10

  real ( kind = 8 ) c(0:n_max)
  real ( kind = 8 ) cx1
  real ( kind = 8 ) cx2(0:n_max)
  integer ( kind = 4 ) n
  real ( kind = 8 ) s
  real ( kind = 8 ) x

  s = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  CARDAN_POLY_COEF returns the coefficients of a'
  write ( *, '(a)' ) '  Cardan polynomial.'
  write ( *, '(a)' ) '  CARDAN evaluates a Cardan polynomial directly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  We use the parameter S = ', s
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Table of polynomial coefficients:'
  write ( *, '(a)' ) ' '

  do n = 0, n_max
    call cardan_poly_coef ( n, s, c )
    write ( *, '(2x,i2,11f7.0)' ) n, c(0:n)
  end do

  s = 0.5D+00
  x = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare CARDAN_POLY_COEF + R8POLY_VAL_HORNER'
  write ( *, '(a)' ) '  versus CARDAN alone.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Evaluate polynomials at X = ', x
  write ( *, '(a,g14.6)' ) '  We use the parameter S = ', s
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Order, Horner, Direct'
  write ( *, '(a)' ) ' '

  call cardan ( n, x, s, cx2 )

  do n = 0, n_max

    call cardan_poly_coef ( n, s, c )
    call r8poly_val_horner ( n, c, x, cx1 )

    write ( *, '(2x,i2,2g14.6)' ) n, cx1, cx2(n)

  end do

  return
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests CATALAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2(0:10)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  CATALAN computes Catalan numbers.'
  write ( *, '(a)' ) '  CATALAN_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  N  exact C(I)  computed C(I)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call catalan_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call catalan ( n, c2 )

    write ( *, '(2x,i4,2i8)' ) n, c, c2(n)

  end do
 
  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests CATALAN_ROW_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) c(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First, compute row 7:'

  ido = 0
  i = 7
  call catalan_row_next ( ido, i, c )
  write ( *, '(2x,i2,2x,11i6)' ) i, c(0:i)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now compute rows one at a time:'
  write ( *, '(a)' ) ' '

  ido = 0
 
  do i = 0, n
    call catalan_row_next ( ido, i, c )
    ido = 1
    write ( *, '(2x,i2,2x,11i6)' ) i, c(0:i)
  end do
 
  return
end
subroutine test0175 ( )

!*****************************************************************************80
!
!! TEST0175 tests CHARLIER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) a
  real ( kind = 8 ), dimension ( test_num ) :: a_test = (/ &
    0.25D+00, 0.5D+00, 1.0D+00, 2.0D+00, 10.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ) value(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0175:'
  write ( *, '(a)' ) '  CHARLIER evaluates Charlier polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      A         X        P(N,A,X)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a = a_test(test)

    do j = 0, 5

      x = real ( j, kind = 8 ) / 2.0D+00

      call charlier ( n, a, x, value )

      write ( *, '(a)' ) ' '

      do i = 0, n

        write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,g14.6)' ) i, a, x, value(i)

      end do

    end do

  end do

  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests CHEBY_T_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2(0:n_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018:'
  write ( *, '(a)' ) '  CHEBY_T_POLY evaluates the Chebyshev T polynomial.'
  write ( *, '(a)' ) '  CHEBY_T_POLY_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N      X        Exact F       T(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_t_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call cheby_t_poly ( 1, n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine test0185 ( )

!*****************************************************************************80
!
!! TEST0185 tests CHEBY_T_POLY_ZERO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 4

  real ( kind = 8 ), allocatable :: fx(:,:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) z(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0185:'
  write ( *, '(a)' ) '  CHEBY_T_POLY_ZERO returns zeroes of T(N)(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      X        T(N)(X)'

  do n = 1, n_max

    call cheby_t_poly_zero ( n, z )

    allocate ( fx(1:n,0:n) )

    call cheby_t_poly ( n, n, z, fx )

    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i8,2x,f8.4,2x,g14.6)' ) n, z(i), fx(i,n)
    end do

    deallocate ( fx )

  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests CHEBY_T_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  CHEBY_T_POLY_COEF determines ' // &
    'the Chebyshev T polynomial coefficients.'

  call cheby_t_poly_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  T(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do
 
  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests CHEBY_U_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2(0:n_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020:'
  write ( *, '(a)' ) '  CHEBY_U_POLY evaluates the Chebyshev U polynomial.'
  write ( *, '(a)' ) '  CHEBY_U_POLY_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N      X        Exact F       U(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_u_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call cheby_u_poly ( n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests CHEBY_U_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  CHEBY_U_POLY_COEF determines ' // &
    'the Chebyshev U polynomial coefficients.'

  call cheby_u_poly_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  U(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do
 
  return
end
subroutine test0215 ( )

!*****************************************************************************80
!
!! TEST0215 tests CHEBY_U_POLY_ZERO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 4

  real ( kind = 8 ) fx(0:n_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) z(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0215:'
  write ( *, '(a)' ) '  CHEBY_U_POLY_ZERO returns zeroes of the U(N)(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      X        U(N)(X)'
  write ( *, '(a)' ) ' '

  do n = 1, n_max

    call cheby_u_poly_zero ( n, z )

    do i = 1, n

      call cheby_u_poly ( n, z(i), fx )

      write ( *, '(2x,i8,2x,f8.4,2x,g14.6)' ) n, z(i), fx(n)

    end do

    write ( *, '(a)' ) ' '

  end do

  return
end
subroutine test0216 ( )

!*****************************************************************************80
!
!! TEST0216 tests CHEBYSHEV_DISCRETE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) x
  real ( kind = 8 ) value(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0216:'
  write ( *, '(a)' ) &
    '  CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      M         X        T(N,M,X)'
  write ( *, '(a)' ) ' '

  m = 5

  do j = 0, 5

    x = real ( j, kind = 8 ) / 2.0D+00

    call chebyshev_discrete ( n, m, x, value )

    write ( *, '(a)' ) ' '

    do i = 0, n

      write ( *, '(2x,i8,2x,i8,2x,f8.4,2x,g14.6)' ) i, m, x, value(i)

    end do

  end do

  return
end
subroutine test0217 ( )

!*****************************************************************************80
!
!! TEST0217 tests COLLATZ_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) collatz_count
  integer ( kind = 4 ) count
  integer ( kind = 4 ) count2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0217:'
  write ( *, '(a)' ) '  COLLATZ_COUNT(N) counts the length of the'
  write ( *, '(a)' ) '  Collatz sequence beginning with N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       COUNT(N)     COUNT(N)'
  write ( *, '(a)' ) '              (computed)    (table)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call collatz_count_values ( n_data, n, count )

    if ( n_data == 0 ) then
      exit
    end if

    count2 = collatz_count ( n )

    write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, count, count2

  end do

  return
end
subroutine test0218 ( )

!*****************************************************************************80
!
!! TEST0218 tests COLLATZ_COUNT_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0218:'
  write ( *, '(a)' ) '  COLLATZ_COUNT_MAX(N) returns the length of the'
  write ( *, '(a)' ) '  longest Collatz sequence from 1 to N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N     I_MAX     J_MAX'
  write ( *, '(a)' ) ' '

  n = 10

  do while ( n <= 100000 )

    call collatz_count_max ( n, i_max, j_max )

    write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, i_max, j_max

    n = n * 10

  end do

  return
end
subroutine test024 ( )

!*****************************************************************************80
!
!! TEST024 tests COMB_ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) c(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024'
  write ( *, '(a)' ) '  COMB_ROW computes a row of Pascal''s triangle.'
  write ( *, '(a)' ) ' '
 
  ido = 0
 
  do i = 0, n
    call comb_row ( ido, i, c )
    ido = 1
    write ( *, '(2x,i2,2x,11i5)' ) i, c(0:i)
  end do
 
  return
end
subroutine test02405 ( )

!*****************************************************************************80
!
!! TEST02405 tests COMMUL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) factor(4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ncomb
  integer ( kind = 4 ) nfactor

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02405'
  write ( *, '(a)' ) '  COMMUL computes a multinomial coefficient.'
  write ( *, '(a)' ) ' '

  n = 8
  nfactor = 2
  factor(1) = 6
  factor(2) = 2

  call commul ( n, nfactor, factor, ncomb ) 

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i8)' ) '  Number of factors = ', nfactor
  do i = 1, nfactor
    write ( *, '(2x,i2,2x,i8)' ) i, factor(i)
  end do
  write ( *, '(a,i12)' ) '  Value of coefficient = ', ncomb

  n = 8
  nfactor = 3
  factor(1) = 2
  factor(2) = 2
  factor(3) = 4
  call commul ( n, nfactor, factor, ncomb )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i8)' ) '  Number of factors = ', nfactor
  do i = 1, nfactor
    write ( *, '(2x,i2,2x,i8)' ) i, factor(i)
  end do
  write ( *, '(a,i12)' ) '  Value of coefficient = ', ncomb

  n = 13
  nfactor = 4
  factor(1) = 5
  factor(2) = 3
  factor(3) = 3
  factor(4) = 2
  call commul ( n, nfactor, factor, ncomb )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i8)' ) '  Number of factors = ', nfactor
  do i = 1, nfactor
    write ( *, '(2x,i2,2x,i8)' ) i, factor(i)
  end do
  write ( *, '(a,i12)' ) '  Value of coefficient = ', ncomb

  return
end
subroutine test0241 ( )

!*****************************************************************************80
!
!! TEST0241 tests COS_DEG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) cos_deg
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0241'
  write ( *, '(a)' ) '  COS_DEG computes the cosine of an angle'
  write ( *, '(a)' ) '  given in degrees.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ANGLE    COS_DEG(ANGLE)'
  write ( *, '(a)' ) ' '
 
  do i = 0, 360, 10
    angle = real ( i, kind = 8 )
    write ( *, '(2x,f8.2,2x,g14.6)' )  angle, cos_deg ( angle )
  end do
 
  return
end
subroutine test0243 ( )

!*****************************************************************************80
!
!! TEST0243 tests COS_POWER_INT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cos_power_int
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0243:'
  write ( *, '(a)' ) '  COS_POWER_INT returns values of '
  write ( *, '(a)' ) '  the integral of COS(X)^N from A to B.'
  write ( *, '(a)' ) '  COS_POWER_INT_VALUES stores some selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      A         B          N      Exact           Computed'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cos_power_int_values ( n_data, a, b, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = cos_power_int ( a, b, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,i8,2x,g14.6,2x,g14.6)' ) a, b, n, fx, fx2

  end do

  return
end
subroutine test01155 ( )

!*****************************************************************************80
!
!! TEST01155 tests DELANNOY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 8

  integer ( kind = 4 ) a(0:m,0:n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01155'
  write ( *, '(a)' ) '  DELANNOY computes the Delannoy numbers A(0:M,0:N).'
  write ( *, '(a)' ) '  A(M,N) counts the paths from (0,0) to (M,N).'
  write ( *, '(a)' ) ' '

  call delannoy ( m, n, a )

  do i = 0, m
    write ( *, '(2x,i4,2x,5i4,3i8,i10)' )  i, a(i,0:n)
  end do
 
  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests ERROR_F.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error_f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025:'
  write ( *, '(a)' ) '  ERROR_F evaluates the error function.'
  write ( *, '(a)' ) '  ERF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X      Exact F       ERROR_F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = error_f ( x )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test0255 ( )

!*****************************************************************************80
!
!! TEST0255 tests ERROR_F_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error_f_inverse
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0255:'
  write ( *, '(a)' ) '  ERROR_F_INVERSE inverts the error function.'
  write ( *, '(a)' ) '  ERF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    FX            X    ERROR_F_INVERSE(FX)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x1, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = error_f_inverse ( fx )

    write ( *, '(2x,f8.4,2g14.6)' ) fx, x1, x2

  end do

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests EULER_NUMBER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c1
  integer ( kind = 4 ) c2(0:12)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  EULER_NUMBER computes Euler numbers.'
  write ( *, '(a)' ) '  EULER_NUMBER_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N       exact   EULER_NUMBER'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call euler_number_values ( n_data, n, c1 )

    if ( n_data == 0 ) then
      exit
    end if

    call euler_number ( n, c2 )

    write ( *, '(2x,i4,2i12,g14.6)' ) n, c1, c2(n)

  end do
 
  return
end
subroutine test0265 ( )

!*****************************************************************************80
!
!! TEST0265 tests EULER_NUMBER2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) euler_number2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0265'
  write ( *, '(a)' ) '  EULER_NUMBER2 computes Euler numbers.'
  write ( *, '(a)' ) '  EULER_NUMBER_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N       exact   EULER_NUMBER2'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call euler_number_values ( n_data, n, c1 )

    if ( n_data == 0 ) then
      exit
    end if

    c2 = euler_number2 ( n )

    write ( *, '(2x,i4,i12,g14.6)' ) n, c1, c2

  end do
 
  return
end
subroutine test028 ( )

!*****************************************************************************80
!
!! TEST028 tests EULER_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) euler_poly
  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 15
  real ( kind = 8 ) x

  x = 0.5D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028'
  write ( *, '(a)' ) '  EULER_POLY evaluates Euler polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N      X             F(X)'
  write ( *, '(a)' ) ' '
   
  do i = 0, n
    f = euler_poly ( i, x )
    write ( *, '(2x,i2,2x,2g14.6)' ) i, x, f
  end do
 
  return
end
subroutine test027 ( )

!*****************************************************************************80
!
!!  TEST027 tests EULERIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) e(n,n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  EULERIAN evaluates Eulerian numbers.'
  write ( *, '(a)' ) ' '
 
  call eulerian ( n, e )

  do i = 1, n
    write ( *, '(2x,10i6)' )  e(i,1:n)
  end do
 
  return
end
subroutine test029 ( )

!*****************************************************************************80
!
!! TEST029 tests F_HOFSTADTER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) f
  integer ( kind = 4 ) f_hofstadter
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029'
  write ( *, '(a)' ) '  F_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, '(a)' ) '  F function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N   F(N)'
  write ( *, '(a)' ) ' '

  do i = 0, 30
    f = f_hofstadter ( i )
    write ( *, '(2x,2i8)' ) i, f
  end do

  return
end
subroutine test031 ( )

!*****************************************************************************80
!
!! TEST031 tests FIBONACCI_DIRECT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST031'
  write ( *, '(a)' ) '  FIBONACCI_DIRECT evalutes a Fibonacci number directly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I        F(I)'
  write ( *, '(a)' ) ' '

  n = 20
 
  do i = 1, n
    call fibonacci_direct ( i, f )
    write ( *, '(2x,i8,i10)' ) i, f
  end do
 
  return
end
subroutine test032 ( )

!*****************************************************************************80
!
!! TEST032 tests FIBONACCI_FLOOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032'
  write ( *, '(a)' ) '  FIBONACCI_FLOOR computes the largest Fibonacci number'
  write ( *, '(a)' ) '  less than or equal to a given positive integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N  Fibonacci  Index'
  write ( *, '(a)' ) ' ' 

  do n = 1, 20
    call fibonacci_floor ( n, f, i )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, f, i
  end do
 
  return
end
subroutine test033 ( )

!*****************************************************************************80
!
!! TEST033 tests FIBONACCI_RECURSIVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) f(n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST033'
  write ( *, '(a)' ) '  FIBONACCI_RECURSIVE computes the Fibonacci sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       F(N)'
  write ( *, '(a)' ) ' '
 
  call fibonacci_recursive ( n, f )
 
  do i = 1, n
    write ( *, '(2x,i8,i10)' ) i, f(i)
  end do
 
  return
end
subroutine test034 ( )

!*****************************************************************************80
!
!! TEST034 tests G_HOFSTADTER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) g
  integer ( kind = 4 ) g_hofstadter
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST034'
  write ( *, '(a)' ) '  G_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, '(a)' ) '  G function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N   G(N)'
  write ( *, '(a)' ) ' '

  do i = 0, 30
    g = g_hofstadter ( i )
    write ( *, '(2x,2i8)' ) i, g
  end do

  return
end
subroutine test036 ( )

!*****************************************************************************80
!
!! TEST036 tests GAMMA_LOG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036:'
  write ( *, '(a)' ) '  R8_GAMMA_LOG evaluates the logarithm of the '
  write ( *, '(a)' ) '  Gamma function.'
  write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       GAMMA_LOG(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_gamma_log ( x )

    write ( *, '(2x,f8.4,2g18.10)' ) x, fx, fx2

  end do

  return
end
subroutine test037 ( )

!*****************************************************************************80
!
!! TEST037 tests GEGENBAUER_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), allocatable, dimension ( : ) :: c
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037:'
  write ( *, '(a)' ) '  GEGENBAUER_POLY computes values of '
  write ( *, '(a)' ) '  the Gegenbauer polynomials.'
  write ( *, '(a)' ) '  GEGENBAUER_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Gegenbauer polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N        A           X       GPV      GEGENBAUER'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gegenbauer_poly_values ( n_data, n, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( c(0:n) )

    call gegenbauer_poly ( n, a, x, c )
    fx2 = c(n)

    write ( *, '(2x,i8,2x,f10.4,2x,f10.4,2g14.6)' ) n, a, x, fx, fx2

    deallocate ( c )

  end do

  return
end
subroutine test052 ( )

!*****************************************************************************80
!
!! TEST052 tests GEN_LAGUERRE_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test = (/ &
    0.0D+00, 0.0D+00, 0.1D+00, 0.1D+00, 0.5D+00, 1.0D+00 /)
  real ( kind = 8 ) c(0:n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.5D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST052'
  write ( *, '(a)' ) '  GEN_LAGUERRE_POLY evaluates the generalized Laguerre '
  write ( *, '(a)' ) '  polynomials.'

  do test = 1, test_num

    x = x_test(test)
    alpha = alpha_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Table of L(N,ALPHA)(X) for'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '    N(max) = ', n
    write ( *, '(a,g14.6)' ) '    ALPHA =  ', alpha
    write ( *, '(a,g14.6)' ) '    X =      ', x
    write ( *, '(a)' ) ' '
  
    call gen_laguerre_poly ( n, alpha, x, c )
 
    do j = 0, n
      write ( *, '(2x,i8,2x,g14.6)' ) j, c(j)
    end do

  end do
 
  return
end
subroutine test038 ( )

!*****************************************************************************80
!
!! TEST038 tests GUD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  real ( kind = 8 ) gud
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST038:'
  write ( *, '(a)' ) '  GUD evaluates the Gudermannian function.'
  write ( *, '(a)' ) '  GUD_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X      Exact F       GUD(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gud_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = gud ( x )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test040 ( )

!*****************************************************************************80
!
!! TEST040 tests H_HOFSTADTER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) h
  integer ( kind = 4 ) h_hofstadter
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST040'
  write ( *, '(a)' ) '  H_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, '(a)' ) '  H function.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N   H(N)'
  write ( *, '(a)' ) ' '

  do i = 0, 30
    h = h_hofstadter ( i )
    write ( *, '(2x,i8,2x,i8)' ) i, h
  end do

  return
end
subroutine test041 ( )

!*****************************************************************************80
!
!! TEST041 tests HERMITE_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2(0:n_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST041:'
  write ( *, '(a)' ) '  HERMITE_POLY evaluates the Hermite polynomial.'
  write ( *, '(a)' ) '  HERMITE_POLY_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N    X      Exact F       H(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hermite_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call hermite_poly ( n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine test042 ( )

!*****************************************************************************80
!
!! TEST042 tests HERMITE_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST042'
  write ( *, '(a)' ) &
    '  HERMITE_POLY_COEF determines the Hermite polynomial coefficients.'

  call hermite_poly_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  H(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do
 
  return
end
subroutine test0427 ( )

!*****************************************************************************80
!
!! TEST0427 tests I4_CHOOSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cnk
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0427'
  write ( *, '(a)' ) '  I4_CHOOSE evaluates C(N,K).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     K    CNK'
  write ( *, '(a)' ) ' '
 
  do n = 0, 4
    do k = 0, n
      cnk = i4_choose ( n, k )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, k, cnk
    end do
  end do
 
  return
end
subroutine test043 ( )

!*****************************************************************************80
!
!! TEST043 tests I4_FACTORIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) fn2
  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST043:'
  write ( *, '(a)' ) '  I4_FACTORIAL evaluates the factorial function.'
  write ( *, '(a)' ) '  I4_FACTORIAL_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       I4_FACTORIAL(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i4_factorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    fn2 = i4_factorial ( n )

    write ( *, '(2x,i4,2i12)' ) n, fn, fn2

  end do

  return
end
subroutine test044 ( )

!*****************************************************************************80
!
!! TEST044 tests I4_FACTORIAL2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) fn2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) i4_factorial2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST044:'
  write ( *, '(a)' ) '  I4_FACTORIAL2 evaluates the double factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   Exact  I4_FACTORIAL2(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i4_factorial2_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    fn2 = i4_factorial2 ( n )

    write ( *, '(2x,i4,2i8)' ) n, fn, fn2

  end do

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests PARTITION_COUNT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045:'
  write ( *, '(a)' ) '  For the number of partitions of an integer,'
  write ( *, '(a)' ) '  PARTITION_COUNT_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           N       Exact F'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call partition_count_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,2i10)' ) n, c

  end do

  return
end
subroutine test046 ( )

!*****************************************************************************80
!
!! TEST046 tests I4_PARTITION_DISTINCT_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), parameter :: n_max = 20

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST046:'
  write ( *, '(a)' ) '  For the number of partitions of an integer'
  write ( *, '(a)' ) '  into distinct parts,'
  write ( *, '(a)' ) '  I4_PARTITION_DISTINCT_COUNT'
  write ( *, '(a)' ) '  computes any value.'
  write ( *, '(a)' ) '  PARTITION_DISTINCT_COUNT_VALUES '
  write ( *, '(a)' ) '  returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           N       Exact F    Q(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call partition_distinct_count_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    if ( n_max < n ) then
      cycle
    end if

    call i4_partition_distinct_count ( n, c2 )

    write ( *, '(2x,3i10)' ) n, c, c2

  end do

  return
end
subroutine test047 ( )

!*****************************************************************************80
!
!! TEST047 tests I4_POCHHAMMER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_pochhammer
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST047:'
  write ( *, '(a)' ) '  I4_POCHHAMMER evaluates the integer Pochhammer function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I   J   I4_Pochhammer(I,J)'
  write ( *, '(a)' ) ' '

  i = 3
  j = 3
  k = i4_pochhammer ( i, j, k )

  write ( *, '(2x,i4,i4,i4)' ) i, j,  k

  i = 3
  j = 4
  k = i4_pochhammer ( i, j, k )
  write ( *, '(2x,i4,i4,i4)' ) i, j,  k

  i = 3
  j = 5
  k = i4_pochhammer ( i, j, k )
  write ( *, '(2x,i4,i4,i4)' ) i, j,  k

  return
end
subroutine test048 ( )

!*****************************************************************************80
!
!! TEST048 tests I4_IS_TRIANGULAR, I4_TO_TRIANGLE and TRIANGLE_TO_I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  logical i4_is_triangular
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical l

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST048'
  write ( *, '(a)' ) '  I4_TO_TRIANGLE converts a linear index to a'
  write ( *, '(a)' ) '  triangular one.'
  write ( *, '(a)' ) '  TRIANGLE_TO_I4 converts a triangular index to a'
  write ( *, '(a)' ) '  linear one.'
  write ( *, '(a)' ) '  I4_IS_TRIANGULAR returns T or F depending on'
  write ( *, '(a)' ) '  whether I is triangular.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  =>   J   K  =>   I   T/F'
  write ( *, '(a)' ) ' '

  do i = 0, 20

    call i4_to_triangle ( i, j, k )

    call triangle_to_i4 ( j, k, i2 )

    l = i4_is_triangular ( i )

    write ( *, '(2x,i4,4x,i4,i4,4x,i4,4x,l1)' )  i, j, k, i2, l

  end do
 
  return
end
subroutine test049 ( )

!*****************************************************************************80
!
!! TEST049 tests JACOBI_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: c(:)
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST049:'
  write ( *, '(a)' ) '  JACOBI_POLY computes values of '
  write ( *, '(a)' ) '  the Jacobi polynomial.'
  write ( *, '(a)' ) '  JACOBI_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Jacobi polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       A       B      X       JPV      JACOBI'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_poly_values ( n_data, n, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( c(n+1) )

    call jacobi_poly ( n, a, b, x, c )
    fx2 = c(n+1)

    write ( *, '(2x,i8,2x,f8.4,2x,f8.4,f10.4,2g14.6)' ) n, a, b, x, fx, fx2

    deallocate ( c )

  end do

  return
end
subroutine test050 ( )

!*****************************************************************************80
!
!! TEST050 tests JACOBI_SYMBOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) l
  integer ( kind = 4 ) p
  integer ( kind = 4 ), dimension ( test_num ) :: p_test = (/ 3, 9, 10, 12 /)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST050'
  write ( *, '(a)' ) '  JACOBI_SYMBOL computes the Jacobi symbol'
  write ( *, '(a)' ) '  (Q/P), which records if Q is a quadratic '
  write ( *, '(a)' ) '  residue modulo the number P.'

  do test = 1, test_num
    p = p_test(test)
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Jacobi Symbols for P = ', p
    write ( *, '(a)' ) ' '
    do q = 0, p
      call jacobi_symbol ( q, p, l )
      write ( *, '(2x,3i8)' ) p, q, l
    end do
  end do

  return
end
subroutine test0505 ( )

!*****************************************************************************80
!
!! TEST0505 tests KRAWTCHOUK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: test_num = 2

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) p
  real ( kind = 8 ), dimension ( test_num ) :: p_test = (/ &
    0.25D+00, 0.5D+00 /)
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ) value(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0505:'
  write ( *, '(a)' ) '  KRAWTCHOUK evaluates Krawtchouk polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N      P         X           M      K(N,P,X,M)'
  write ( *, '(a)' ) ' '

  m = 5

  do test = 1, test_num

    p = p_test(test)

    do j = 0, 5

      x = real ( j, kind = 8 ) / 2.0D+00

      call krawtchouk ( n, p, x, m, value )

      write ( *, '(a)' ) ' '

      do i = 0, n

        write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,i8,2x,g14.6)' ) &
          i, p, x, m, value(i)

      end do

    end do

  end do

  return
end
subroutine test051 ( )

!*****************************************************************************80
!
!! TEST051 tests LAGUERRE_ASSOCIATED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6
  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) c(0:n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension ( test_num ) :: m_test = (/ 0, 0, 1, 2, 3, 1 /)
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.5D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST051'
  write ( *, '(a)' ) '  LAGUERRE_ASSOCIATED evaluates the associated Laguerre'
  write ( *, '(a)' ) '  polynomials.'

  do test = 1, test_num

    m = m_test(test)
    x = x_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Table of L(N,M)(X) for'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  N(max) = ', n
    write ( *, '(a,i4)' ) '  M      = ', m
    write ( *, '(a,g14.6)' ) '  X =      ', x
    write ( *, '(a)' ) ' '
 
    call laguerre_associated ( n, m, x, c )
 
    do j = 0, n
      write ( *, '(2x,i8,g14.6)' ) j, c(j)
    end do
 
  end do

  return
end
subroutine test054 ( )

!*****************************************************************************80
!
!! TEST054 tests LAGUERRE_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2(0:n_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST054:'
  write ( *, '(a)' ) '  LAGUERRE_POLY evaluates the Laguerre polynomial.'
  write ( *, '(a)' ) '  LAGUERRE_POLYNOMIAL_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N    X      Exact F       L(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call laguerre_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call laguerre_poly ( n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine test055 ( )

!*****************************************************************************80
!
!! TEST055 tests LAGUERRE_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) c(0:n,0:n)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST055'
  write ( *, '(a)' ) &
    '  LAGUERRE_POLY_COEF determines the Laguerre polynomial coefficients.'

  call laguerre_poly_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  L(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do
 
  fact = 1.0D+00

  do i = 0, n

    if ( 0 < i ) then
      fact = fact * real ( i, kind = 8 )
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  Factorially scaled L(', i, ')'
    write ( *, '(a)' ) ' '

    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) fact * c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) fact * c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) fact * c(i,j), ' * x**', j
      end if
    end do
    
  end do

  return
end
subroutine test0552 ( )

!*****************************************************************************80
!
!! TEST0552 tests LAMBERT_W, LAMBERT_W_CRUDE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  real ( kind = 8 ) fx3
  real ( kind = 8 ) lambert_w
  real ( kind = 8 ) lambert_w_crude
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0552:'
  write ( *, '(a)' ) '  LAMBERT_W estimates the Lambert W function.'
  write ( *, '(a)' ) '  LAMBERT_W_CRUDE makes a crude estimate of the'
  write ( *, '(a)' ) '  Lambert W function.'
  write ( *, '(a)' ) '  LAMBERT_W_VALUES returns some tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           X           W(X)        W(X)         W(X)'
  write ( *, '(a)' ) '                   Tabulated       Crude     Estimate'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lambert_w_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = lambert_w_crude ( x )

    fx3 = lambert_w ( x )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) x, fx, fx2, fx3

  end do

  return
end
subroutine test059 ( )

!*****************************************************************************80
!
!! TEST059 tests LEGENDRE_ASSOCIATED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx2(0:n_max)
  real ( kind = 8 ) fx
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST059:'
  write ( *, '(a)' ) &
    '  LEGENDRE_ASSOCIATED evaluates associated Legendre functions.'
  write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       M    X     Exact F       PNM(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call legendre_associated ( n, m, x, fx2 )

    write ( *, '(2x,i8,2x,i8,f8.4,2g14.6)' ) n, m, x, fx, fx2(n)

  end do

  return
end
subroutine test0595 ( )

!*****************************************************************************80
!
!! TEST0595 tests LEGENDRE_ASSOCIATED_NORMALIZED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx2(0:n_max)
  real ( kind = 8 ) fx
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST059:'
  write ( *, '(a)' ) &
    '  LEGENDRE_ASSOCIATED_NORMALIZED evaluates associated Legendre functions.'
  write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       M    X     Exact F       PNM(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_normalized_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call legendre_associated_normalized ( n, m, x, fx2 )

    write ( *, '(2x,i8,2x,i8,f8.4,2g14.6)' ) n, m, x, fx, fx2(n)

  end do

  return
end
subroutine test060 ( )

!*****************************************************************************80
!
!! TEST060 tests LEGENDRE_FUNCTION_Q.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2(0:n_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST060:'
  write ( *, '(a)' ) '  LEGENDRE_FUNCTION_Q evaluates the Legendre Q function.'
  write ( *, '(a)' ) '  LEGENDRE_FUNCTION_Q_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N    X      Exact F       Q(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_function_q_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call legendre_function_q ( n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine test057 ( )

!*****************************************************************************80
!
!! TEST057 tests LEGENDRE_POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 12

  real ( kind = 8 ) fx
  real ( kind = 8 ) fp2(0:n_max)
  real ( kind = 8 ) fx2(0:n_max)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST057:'
  write ( *, '(a)' ) '  LEGENDRE_POLY evaluates the Legendre PN function.'
  write ( *, '(a)' ) '  LEGENDRE_POLY_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N    X      Exact F       P(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call legendre_poly ( n, x, fx2, fp2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine test058 ( )

!*****************************************************************************80
!
!! TEST058 tests LEGENDRE_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST058'
  write ( *, '(a)' ) &
    '  LEGENDRE_POLY_COEF returns Legendre polynomial coefficients.'

  call legendre_poly_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  P(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do

  return
end
subroutine test061 ( )

!*****************************************************************************80
!
!! TEST061 tests LEGENDRE_SYMBOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) l
  integer ( kind = 4 ) p
  integer ( kind = 4 ), dimension ( test_num ) :: p_test = (/ 7, 11, 13, 17 /)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST061'
  write ( *, '(a)' ) '  LEGENDRE_SYMBOL computes the Legendre'
  write ( *, '(a)' ) '  symbol (Q/P) which records whether Q is '
  write ( *, '(a)' ) '  a quadratic residue modulo the prime P.'

  do test = 1, test_num
    p = p_test(test)
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Legendre Symbols for P = ', p
    write ( *, '(a)' ) ' '
    do q = 0, p
      call legendre_symbol ( q, p, l )
      write ( *, '(2x,3i8)' ) p, q, l
    end do
  end do

  return
end
subroutine test0615 ( )

!*****************************************************************************80
!
!! TEST0615 tests LERCH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  real ( kind = 8 ) lerch
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) s
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0615'
  write ( *, '(a)' ) '  LERCH computes the Lerch function.'
  write ( *, '(a)' ) '  LERCH_VALUES returns some tabulated values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       Z       S       A         Lerch           Lerch'
  write ( *, '(a)' ) '                             Tabulated        Computed'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call lerch_values ( n_data, z, s, a, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = lerch ( z, s, a )

    write ( *, '(2x,f8.4,2x,i4,2x,f8.4,2x,g14.6,2x,g14.6)' ) z, s, a, fx, fx2

  end do
 
  return
end
subroutine test062 ( )

!*****************************************************************************80
!
!! TEST062 tests LOCK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(0:n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST062'
  write ( *, '(a)' ) '  LOCK counts the combinations on a button lock.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I      LOCK(I)'
  write ( *, '(a)' ) ' '

  call lock ( n, a )

  do i = 0, n
    write ( *, '(2x,i8,2x,i10)' )  i, a(i)
  end do
 
  return
end
subroutine test0623 ( )

!*****************************************************************************80
!
!! TEST0623 tests MEIXNER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) beta
  real ( kind = 8 ) :: beta_test(test_num) = (/ &
    0.5D+00, 1.0D+00, 2.0D+00 /)
  real ( kind = 8 ) c
  real ( kind = 8 ) :: c_test(test_num) = (/ &
    0.125D+00, 0.25D+00, 0.5D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) test
  real ( kind = 8 ) v(0:n)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLPAK_TEST0623:'
  write ( *, '(a)' ) '  MEIXNER evaluates Meixner polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      BETA         C         X        M(N,BETA,C,X)'

  do test = 1, test_num

    beta = beta_test(test)
    c = c_test(test)

    do j = 0, 5

      x = real ( j, kind = 8 ) / 2.0D+00

      call meixner ( n, beta, c, x, v )

      write ( *, '(a)' ) ' '

      do i = 0, n

        write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) &
          i, beta, c, x, v(i)

      end do

    end do

  end do

  return
end
subroutine test0625 ( )

!*****************************************************************************80
!
!! TEST0625 tests MERTENS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) mertens
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0625'
  write ( *, '(a)' ) '  MERTENS computes the Mertens function.'
  write ( *, '(a)' ) '  MERTENS_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '         N     Exact   MERTENS(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call mertens_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    c2 = mertens ( n )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine test063 ( )

!*****************************************************************************80
!
!! TEST063 tests MOEBIUS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST063'
  write ( *, '(a)' ) '  MOEBIUS computes the Moebius function.'
  write ( *, '(a)' ) '  MOEBIUS_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '         N     Exact   MOEBIUS(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call moebius_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call moebius ( n, c2 )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine test0635 ( )

!*****************************************************************************80
!
!! TEST0635 tests MOTZKIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(0:n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0635'
  write ( *, '(a)' ) '  MOTZKIN computes the Motzkin numbers A(0:N).'
  write ( *, '(a)' ) '  A(N) counts the paths from (0,0) to (N,0).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         A(I)'
  write ( *, '(a)' ) ' '

  call motzkin ( n, a )

  do i = 0, n
    write ( *, '(2x,i8,2x,i10)' )  i, a(i)
  end do
 
  return
end
subroutine test064 ( )

!*****************************************************************************80
!
!! TEST064 tests OMEGA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST064'
  write ( *, '(a)' ) '  OMEGA counts the distinct prime divisors of an integer N.'
  write ( *, '(a)' ) '  OMEGA_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '             N      Exact   OMEGA(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call omega_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call omega ( n, c2 )

    write ( *, '(2x,i12,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests PENTAGON_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  PENTAGON_NUM computes the pentagonal numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I      Pent(I)'
  write ( *, '(a)' ) ' '

  do n = 1, 10
    call pentagon_num ( n, p )
    write ( *, '(2x,i8,2x,i8)' ) n, p
  end do
 
  return
end
subroutine test066 ( )

!*****************************************************************************80
!
!! TEST066 tests PHI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST066'
  write ( *, '(a)' ) '  PHI computes the PHI function.'
  write ( *, '(a)' ) '  PHI_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N     Exact     PHI(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call phi_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call phi ( n, c2 )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine test0665 ( )

!*****************************************************************************80
!
!! TEST0665 tests POLY_BERNOULLI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) b
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0665'
  write ( *, '(a)' ) '  POLY_BERNOULLI computes the poly-Bernoulli numbers'
  write ( *, '(a)' ) '  of negative index, B_n^(-k)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     K    B_N^(-K)'
  write ( *, '(a)' ) ' '

  do k = 0, 6
    write ( *, '(a)' ) ' '
    do n = 0, 6

      call poly_bernoulli ( n, k, b )

      write ( *, '(2x,i4,2x,i4,2x,i12)' ) n, k, b

    end do
  end do

  return
end
subroutine test0667 ( )

!*****************************************************************************80
!
!! TEST0667 tests POLY_COEF_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) degree
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) poly_coef_count

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0667'
  write ( *, '(a)' ) '  POLY_COEF_COUNT counts the number of coefficients'
  write ( *, '(a)' ) '  in a polynomial of degree DEGREE and dimension DIM'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Dimension    Degree     Count'

  do dim = 1, 10, 3
    write ( *, '(a)' ) ' '
    do degree = 0, 5
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
       dim, degree, poly_coef_count ( dim, degree )
    end do
  end do
 
  return
end
subroutine test067 ( )

!*****************************************************************************80
!
!! TEST067 tests PYRAMID_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) pyramid_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST067'
  write ( *, '(a)' ) '  PYRAMID_NUM computes the pyramidal numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I    PYR(I)'
  write ( *, '(a)' ) ' '

  do n = 1, 10
    write ( *, '(2x,i8,2x,i8)' ) n, pyramid_num ( n )
  end do
 
  return
end
subroutine test0675 ( )

!*****************************************************************************80
!
!! TEST0675 tests R8_ACOSH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_acosh
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0675'
  write ( *, '(a)' ) '  R8_ACOSH computes the inverse hyperbolic cosine'
  write ( *, '(a)' ) '  of a given value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X        R8_ACOSH(X)     COSH(R8_ACOSH(X))'
  write ( *, '(a)' ) ' '

  do i = 0, 10
    x = 1.0D+00 + real ( i, kind = 8 ) / 5.0D+00
    a = r8_acosh ( x )
    x2 = cosh ( a )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, a, x2
  end do

  return
end
subroutine test0676 ( )

!*****************************************************************************80
!
!! TEST0676 tests R8_ASINH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_asinh
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0676'
  write ( *, '(a)' ) '  R8_ASINH computes the inverse hyperbolic sine'
  write ( *, '(a)' ) '  of a given value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       R8_ASINH(X)     SINH(R8_ASINH(X))'
  write ( *, '(a)' ) ' '

  do i = 0, 10
    x = 1.0D+00 + real ( i, kind = 8 ) / 5.0D+00
    a = r8_asinh ( x )
    x2 = sinh ( a )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, a, x2
  end do

  return
end
subroutine test06765 ( )

!*****************************************************************************80
!
!! TEST06765 tests R8_ATANH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_atanh
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06765'
  write ( *, '(a)' ) '  R8_ATANH computes the inverse hyperbolic tangent'
  write ( *, '(a)' ) '  of a given value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X       R8_ATANH(X)     TANH(R8_ATANH(X))'
  write ( *, '(a)' ) ' '

  do i = -2, 9
    x = real ( i, kind = 8 ) / 10.0D+00
    a = r8_atanh ( x )
    x2 = tanh ( a )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, a, x2
  end do

  return
end
subroutine test0242 ( )

!*****************************************************************************80
!
!! TEST0242 tests R8_CAS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle_rad
  real ( kind = 8 ) r8_cas
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0242'
  write ( *, '(a)' ) '  R8_CAS computes the "casine" of an angle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ANGLE    R8_CAS(ANGLE)'
  write ( *, '(a)' ) ' '
 
  do i = 0, 360, 10
    angle = real ( i, kind = 8 )
    angle_rad = angle * pi / 180.0D+00
    write ( *, '(2x,f8.2,2x,g14.6)' )  angle, r8_cas ( angle_rad )
  end do
 
  return
end
subroutine test022 ( )

!*****************************************************************************80
!
!! TEST022 tests R8_CHOOSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cnk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_choose

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022'
  write ( *, '(a)' ) '  R8_CHOOSE evaluates C(N,K).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N       K      CNK'
  write ( *, '(a)' ) ' '
 
  do n = 0, 4
    do k = 0, n
      cnk = r8_choose ( n, k )
      write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) n, k, cnk
    end do
  end do
 
  return
end
subroutine test068 ( )

!*****************************************************************************80
!
!! TEST068 tests R8_FACTORIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fn
  real ( kind = 8 ) fn2
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST068:'
  write ( *, '(a)' ) '  R8_FACTORIAL evaluates the factorial function.'
  write ( *, '(a)' ) '  R8_FACTORIAL_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N       Exact F       R8_FACTORIAL(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_factorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    fn2 = r8_factorial ( n )

    write ( *, '(2x,i4,2g14.6)' ) n, fn, fn2

  end do

  return
end
subroutine test0685 ( )

!*****************************************************************************80
!
!! TEST0685 tests R8_FACTORIAL_LOG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fn
  real ( kind = 8 ) fn2
  real ( kind = 8 ) r8_factorial_log
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0685:'
  write ( *, '(a)' ) '  R8_FACTORIAL_LOG evaluates the logarithm of the '
  write ( *, '(a)' ) '  factorial function.'
  write ( *, '(a)' ) '  R8_FACTORIAL_LOG_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N	   Exact F	 R8_FACTORIAL_LOG(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_factorial_log_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    fn2 = r8_factorial_log ( n )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, fn, fn2

  end do

  return
end
subroutine test06855 ( )

!*****************************************************************************80
!
!! TEST06855 tests R8_GAMMA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06855:'
  write ( *, '(a)' ) '  R8_GAMMA evaluates the Gamma function.'
  write ( *, '(a)' ) '  GAMMA_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X         Gamma(X)                   Gamma(X)' &
  // '               DIFF'
  write ( * , '(a)' ) '               (Tabulated)                (R8_GAMMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_gamma ( x )

    write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test0425 ( )

!*****************************************************************************80
!
!! TEST0425 tests R8_HYPER_2F1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0425:'
  write ( *, '(a)' ) '  R8_HYPER_2F1 evaluates the hypergeometric 2F1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      A       B       C       X      ', &
  ' 2F1                       2F1                     DIFF'
  write ( *, '(a,a)' ) '                                     ', &
  '(tabulated)               (computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hyper_2f1_values ( n_data, a, b, c, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call r8_hyper_2f1 ( a, b, c, x, fx2 )

    write ( *, &
    '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, b, c, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test06856 ( )

!*****************************************************************************80
!
!! TEST06856 tests R8_PSI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06856:'
  write ( *, '(a)' ) '  R8_PSI evaluates the Psi function.'
  write ( *, '(a)' ) '  PSI_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X         Psi(X)                     Psi(X)  ' &
  // '               DIFF'
  write ( * , '(a)' ) '               (Tabulated)                (R8_PSI)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_psi ( x )

    write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test069 ( )

!*****************************************************************************80
!
!! TEST069 tests SIGMA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST069'
  write ( *, '(a)' ) '  SIGMA computes the SIGMA function.'
  write ( *, '(a)' ) '  SIGMA_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N     Exact   SIGMA(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call sigma_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call sigma ( n, c2 )

    write ( *, '(2x,i4,2i10)' ) n, c, c2

  end do
 
  return
end
subroutine test0695 ( )

!*****************************************************************************80
!
!! TEST0695 tests SIN_POWER_INT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sin_power_int

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0695:'
  write ( *, '(a)' ) '  SIN_POWER_INT returns values of '
  write ( *, '(a)' ) '  the integral of SIN(X)^N from A to B.'
  write ( *, '(a)' ) '  SIN_POWER_INT_VALUES stores some selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      A         B          N      Exact           Computed'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sin_power_int_values ( n_data, a, b, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = sin_power_int ( a, b, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,i8,2x,g14.6,2x,g14.6)' ) a, b, n, fx, fx2

  end do

  return
end
subroutine test0696 ( )

!*****************************************************************************80
!
!! TEST0696 tests SLICE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_max = 5
  integer ( kind = 4 ), parameter :: slice_max = 8

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) p(dim_max,slice_max)
  integer ( kind = 4 ) piece_num
  integer ( kind = 4 ) slice_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0696:'
  write ( *, '(a)' ) '  SLICE determines the maximum number of pieces created'
  write ( *, '(a)' ) '  by SLICE_NUM slices in a DIM_NUM space.'

  do dim_num = 1, dim_max
    do slice_num = 1, slice_max
      call slice ( dim_num, slice_num, piece_num )
      p(dim_num,slice_num) = piece_num
    end do
  end do

  call i4mat_print ( dim_max, slice_max, p, '  Slice Array:' )

  return
end
subroutine test0697 ( )

!*****************************************************************************80
!
!! TEST0697 tests SPHERICAL_HARMONIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) c(0:n_max)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) phi
  real ( kind = 8 ) s(0:n_max)
  real ( kind = 8 ) theta
  real ( kind = 8 ) yi
  real ( kind = 8 ) yi2
  real ( kind = 8 ) yr
  real ( kind = 8 ) yr2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0697:'
  write ( *, '(a)' ) '  SPHERICAL_HARMONIC evaluates spherical harmonic'
  write ( *, '(a)' ) '  functions.'
  write ( *, '(a)' ) '  SPHERICAL_HARMONIC_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       L       M   THETA    PHI     C              S'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call spherical_harmonic_values ( n_data, l, m, theta, phi, yr, yi )

    if ( n_data == 0 ) then
      exit
    end if

    call spherical_harmonic ( l, m, theta, phi, c, s )

    yr2 = c(l)
    yi2 = s(l)

    write ( *, '(2x,i8,2x,i6,2f8.4,2g14.6)' ) l, m, theta, phi, yr,  yi
    write ( *, '(2x,8x,2x,6x,16x,  2g14.6)' )                   yr2, yi2

  end do

  return
end
subroutine test070 ( )

!*****************************************************************************80
!
!! TEST070 tests STIRLING1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) s1(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST070'
  write ( *, '(a)' ) '  STIRLING1: Stirling numbers of first kind.'
  write ( *, '(a,i8)' ) '  Get rows 1 through ', m
  write ( *, '(a)' ) ' '
 
  call stirling1 ( m, n, s1 )
 
  do i = 1, m
    write ( *, '(2x,i8,8i8)' ) i, s1(i,1:n)
  end do
 
  return
end
subroutine test071 ( )

!*****************************************************************************80
!
!! TEST071 tests STIRLING2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) s2(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST071'
  write ( *, '(a)' ) '  STIRLING2: Stirling numbers of second kind.'
  write ( *, '(a,i4)' ) '  Get rows 1 through ', m
  write ( *, '(a)' ) ' '
 
  call stirling2 ( m, n, s2 )
 
  do i = 1, m
    write ( *, '(2x,i8,8i8)' ) i, s2(i,1:n)
  end do
 
  return
end
subroutine test072 ( )

!*****************************************************************************80
!
!! TEST072 tests TAU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST072'
  write ( *, '(a)' ) '  TAU computes the Tau function.'
  write ( *, '(a)' ) '  TAU_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '         N  exact C(I)  computed C(I)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call tau_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call tau ( n, c2 )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine test073 ( )

!*****************************************************************************80
!
!! TEST073 tests TETRAHEDRON_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) tetrahedron_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST073'
  write ( *, '(a)' ) '  TETRAHEDRON_NUM computes the tetrahedron numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I    TETR(I)'
  write ( *, '(a)' ) ' '

  do n = 1, 10
    write ( *, '(2x,i8,2x,i8)' ) n, tetrahedron_num ( n )
  end do
 
  return
end
subroutine test074 ( )

!*****************************************************************************80
!
!! TEST074 tests TRIANGLE_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) triangle_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST074'
  write ( *, '(a)' ) '  TRIANGLE_NUM computes the triangular numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I    TRI(I)'
  write ( *, '(a)' ) ' '
 
  do n = 1, 10
    write ( *, '(2x,i8,2x,i8)' ) n, triangle_num ( n )
  end do
 
  return
end
subroutine test075 ( )

!*****************************************************************************80
!
!! TEST075 tests V_HOFSTADTER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) v
  integer ( kind = 4 ) v_hofstadter

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST075'
  write ( *, '(a)' ) '  V_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, '(a)' ) '  V function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N   V(N)'
  write ( *, '(a)' ) ' '

  do i = 0, 30
    v = v_hofstadter ( i )
    write ( *, '(2x,i8,2x,i8)' ) i, v
  end do

  return
end
subroutine test076 ( )

!*****************************************************************************80
!
!! TEST076 tests VIBONACCI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20
  integer ( kind = 4 ), parameter :: n_time = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) v(n,n_time)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST076'
  write ( *, '(a)' ) '  VIBONACCI computes a Vibonacci sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of times we compute the series: ', n_time
  write ( *, '(a)' ) ' '

  seed = 123456789

  do j = 1, n_time
    call vibonacci ( n, seed, v(1,j) ) 
  end do

  do i = 1, n
    write ( *, '(2x,i8,2x,3i8)' ) i, v(i,1:n_time)
  end do
 
  return
end
subroutine test0773 ( )

!*****************************************************************************80
!
!! TEST0773 tests ZERNIKE_POLY and ZERNIKE_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 5

  real ( kind = 8 ) c(0:n_max)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) rho
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0773'
  write ( *, '(a)' ) '  ZERNIKE_POLY_COEF returns the coefficients of a'
  write ( *, '(a)' ) '  Zernike polynomial.'
  write ( *, '(a)' ) '  ZERNIKE_POLY evaluates a Zernike polynomial directly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Table of polynomial coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   M'
  write ( *, '(a)' ) ' '

  do n = 0, 5

    write ( *, '(a)' ) ' '

    do m = 0, n
      call zernike_poly_coef ( m, n, c )
      write ( *, '(2x,i2,2x,i2,2x,11f7.0)' ) n, m, c(0:n)
    end do

  end do

  rho = 0.987654321D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Z1: Compute polynomial coefficients,'
  write ( *, '(a)' ) '  then evaluate by Horner''s method;'
  write ( *, '(a)' ) '  Z2: Evaluate directly by recursion.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   M       Z1              Z2'
  write ( *, '(a)' ) ' '

  do n = 0, 5

    write ( *, '(a)' ) ' '

    do m = 0, n

      call zernike_poly_coef ( m, n, c )
      call r8poly_val_horner ( n, c, rho, z1 )

      call zernike_poly ( m, n, rho, z2 )

      write ( *, '(2x,i2,2x,i2,2x,g16.8,2x,g16.8)' ) n, m, z1, z2

    end do

  end do

  return
end
subroutine test077 ( )

!*****************************************************************************80
!
!! TEST077 tests ZECKENDORF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m_max = 20

  integer ( kind = 4 ) i_list(m_max)
  integer ( kind = 4 ) f_list(m_max)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST077'
  write ( *, '(a)' ) '  ZECKENDORF computes the Zeckendorf decomposition of'
  write ( *, '(a)' ) '  an integer N into nonconsecutive Fibonacci numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N Sum M Parts'
  write ( *, '(a)' ) ' '

  do n = 1, 100

    call zeckendorf ( n, m_max, m, i_list, f_list )

    write ( *, '(2x,i8,2x,15i4)' ) n, f_list(1:m)

  end do

  return
end
subroutine test0775 ( )

!*****************************************************************************80
!
!! TEST0775 tests ZERNIKE_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) c(0:n)
  integer ( kind = 4 ) m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0775'
  write ( *, '(a)' ) '  ZERNIKE_POLY_COEF determines the Zernike'
  write ( *, '(a)' ) '  polynomial coefficients.'

  do m = 0, n

    call zernike_poly_coef ( m, n, c )
 
    call r8poly_print ( n, c, '  Zernike polynomial' )

  end do

  return
end
subroutine test078 ( )

!*****************************************************************************80
!
!! TEST078 tests ZETA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 June 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) n_real
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2
  real ( kind = 8 ) zeta

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST078'
  write ( *, '(a)' ) '  ZETA computes the Zeta function.'
  write ( *, '(a)' ) '  ZETA_VALUES returns some exact values.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N    exact Zeta    computed Zeta'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call zeta_values ( n_data, n, z1 )

    if ( n_data == 0 ) then
      exit
    end if

    n_real = real ( n, kind = 8 )

    z2 = zeta ( n_real )

    write ( *, '(2x,i8,2x,g20.12,2x,g20.12)' ) n, z1, z2

  end do

  return
end
