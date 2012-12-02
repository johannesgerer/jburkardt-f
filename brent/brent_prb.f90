program main

!*****************************************************************************80
!
!! MAIN is the main program for BRENT_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BRENT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the BRENT library.'

  call test_zero_all ( )
  call test_zero_rc_all ( )
  call test_local_min_all ( )
  call test_local_min_rc_all ( )
  call test_glomin_all ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BRENT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test_zero_all ( )

!*****************************************************************************80
!
!! TEST_ZERO_ALL tests Brent's zero finding routine on all test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f_01
  real ( kind = 8 ), external :: f_02
  real ( kind = 8 ), external :: f_03
  real ( kind = 8 ), external :: f_04
  real ( kind = 8 ), external :: f_05
  real ( kind = 8 ) machep
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) t
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ZERO_ALL'
  write ( *, '(a)' ) '  Test the Brent ZERO routine, which seeks'
  write ( *, '(a)' ) '  a root of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B].'

  machep = epsilon ( machep )
  t = machep

  a = 1.0D+00
  b = 2.0D+00

  call test_zero_one ( a, b, machep, t, f_01, &
    'f_01(x) = sin ( x ) - x / 2' )

  a = 0.0D+00
  b = 1.0D+00

  call test_zero_one ( a, b, machep, t, f_02, &
    'f_02(x) = 2 * x - exp ( - x )' )

  a = -1.0D+00
  b =  0.5D+00

  call test_zero_one ( a, b, machep, t, f_03, &
    'f_03(x) = x * exp ( - x )' )

  a =  0.0001D+00
  b =  20.0D+00

  call test_zero_one ( a, b, machep, t, f_04, &
    'f_04(x) = exp ( x ) - 1 / ( 100 * x * x )' )

  a = -5.0D+00
  b =  2.0D+00

  call test_zero_one ( a, b, machep, t, f_05, &
    'f_05(x) = (x+3) * (x-1) * (x-1)' )

  return
end
subroutine test_zero_rc_all ( )

!*****************************************************************************80
!
!! TEST_ZERO_RC_ALL tests Brent's zero finding routine on all test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f_01
  real ( kind = 8 ), external :: f_02
  real ( kind = 8 ), external :: f_03
  real ( kind = 8 ), external :: f_04
  real ( kind = 8 ), external :: f_05
  real ( kind = 8 ) machep
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) t
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ZERO_RC_ALL'
  write ( *, '(a)' ) '  Test the reverse communication version of the Brent ZERO'
  write ( *, '(a)' ) '  routine, which seeks a root of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B].'

  machep = epsilon ( machep )
  t = machep

  a = 1.0D+00
  b = 2.0D+00

  call test_zero_rc_one ( a, b, machep, t, f_01, &
    'f_01(x) = sin ( x ) - x / 2' )

  a = 0.0D+00
  b = 1.0D+00

  call test_zero_rc_one ( a, b, machep, t, f_02, &
    'f_02(x) = 2 * x - exp ( - x )' )

  a = -1.0D+00
  b =  0.5D+00

  call test_zero_rc_one ( a, b, machep, t, f_03, &
    'f_03(x) = x * exp ( - x )' )

  a =  0.0001D+00
  b =  20.0D+00

  call test_zero_rc_one ( a, b, machep, t, f_04, &
    'f_04(x) = exp ( x ) - 1 / ( 100 * x * x )' )

  a = -5.0D+00
  b =  2.0D+00

  call test_zero_rc_one ( a, b, machep, t, f_05, &
    'f_05(x) = (x+3) * (x-1) * (x-1)' )

  return
end
subroutine test_local_min_all ( )

!*****************************************************************************80
!
!! TEST_LOCAL_MIN_ALL tests Brent's local minimizer on all test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: g_01
  real ( kind = 8 ), external :: g_02
  real ( kind = 8 ), external :: g_03
  real ( kind = 8 ), external :: g_04
  real ( kind = 8 ), external :: g_05
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) t
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LOCAL_MIN_ALL'
  write ( *, '(a)' ) '  Test the Brent LOCAL_MIN routine, which seeks'
  write ( *, '(a)' ) '  a local minimizer of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B].'

  eps = 10.0D+00 * sqrt ( epsilon ( eps ) )
  t = eps

  a = 0.0D+00
  b = 3.141592653589793D+00

  call test_local_min_one ( a, b, eps, t, g_01, &
    'g_01(x) = ( x - 2 ) * ( x - 2 ) + 1' )

  a = 0.0D+00
  b = 1.0D+00

  call test_local_min_one ( a, b, eps, t, g_02, &
    'g_02(x) = x * x + exp ( - x )' )

  a = -2.0D+00
  b =  2.0D+00

  call test_local_min_one ( a, b, eps, t, g_03, &
    'g_03(x) = x^4 + 2x^2 + x + 3' )

  a =  0.0001D+00
  b =  1.0D+00

  call test_local_min_one ( a, b, eps, t, g_04, &
    'g_04(x) = exp ( x ) + 1 / ( 100 x )' )

  a =  0.0002D+00
  b = 2.0D+00

  call test_local_min_one ( a, b, eps, t, g_05, &
    'g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)' )

  return
end
subroutine test_local_min_rc_all ( )

!*****************************************************************************80
!
!! TEST_LOCAL_MIN_RC_ALL tests LOCAL_MIN_RC on all test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: g_01
  real ( kind = 8 ), external :: g_02
  real ( kind = 8 ), external :: g_03
  real ( kind = 8 ), external :: g_04
  real ( kind = 8 ), external :: g_05
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) t
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LOCAL_MIN_RC_ALL'
  write ( *, '(a)' ) '  Test the reverse communication version of '
  write ( *, '(a)' ) '  the Brent LOCALM routine, which seeks'
  write ( *, '(a)' ) '  a local minimizer of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B].'

  t = 10.0D+00 * sqrt ( epsilon ( t ) )

  a = 0.0D+00
  b = 3.141592653589793D+00

  call test_local_min_rc_one ( a, b, t, g_01, &
    'g_01(x) = ( x - 2 ) * ( x - 2 ) + 1' )

  a = 0.0D+00
  b = 1.0D+00

  call test_local_min_rc_one ( a, b, t, g_02, &
    'g_02(x) = x * x + exp ( - x )' )

  a = -2.0D+00
  b =  2.0D+00

  call test_local_min_rc_one ( a, b, t, g_03, &
    'g_03(x) = x^4 + 2x^2 + x + 3' )

  a =  0.0001D+00
  b =  1.0D+00

  call test_local_min_rc_one ( a, b, t, g_04, &
    'g_04(x) = exp ( x ) + 1 / ( 100 x )' )

  a =  0.0002D+00
  b = 2.0D+00

  call test_local_min_rc_one ( a, b, t, g_05, &
    'g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)' )

  return
end
subroutine test_glomin_all ( )

!*****************************************************************************80
!
!! TEST_GLOMIN_ALL tests the Brent global minimizer on all test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) e
  real ( kind = 8 ), external :: h_01
  real ( kind = 8 ), external :: h_02
  real ( kind = 8 ), external :: h_03
  real ( kind = 8 ), external :: h_04
  real ( kind = 8 ), external :: h_05
  real ( kind = 8 ) m
  real ( kind = 8 ) machep
  real ( kind = 8 ) t
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_GLOMIN_ALL'
  write ( *, '(a)' ) '  Test the Brent GLOMIN routine, which seeks'
  write ( *, '(a)' ) '  a global minimizer of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A,B],'
  write ( *, '(a)' ) '  given some upper bound M for F".'

  machep = epsilon ( machep )
  e = sqrt ( machep )
  t = sqrt ( machep )

  a = 7.0D+00
  b = 9.0D+00
  c = ( a + b ) / 2.0D+00
  m = 0.0D+00

  call test_glomin_one ( a, b, c, m, machep, e, t, h_01, &
    'h_01(x) = 2 - x' )

  a = 7.0D+00
  b = 9.0D+00
  c = ( a + b ) / 2.0D+00
  m = 100.0D+00

  call test_glomin_one ( a, b, c, m, machep, e, t, h_01, &
    'h_01(x) = 2 - x' )

  a = -1.0D+00
  b = +2.0D+00
  c = ( a + b ) / 2.0D+00
  m = 2.0D+00

  call test_glomin_one ( a, b, c, m, machep, e, t, h_02, &
    'h_02(x) = x * x' )

  a = -1.0D+00
  b = +2.0D+00
  c = ( a + b ) / 2.0D+00
  m = 2.1D+00

  call test_glomin_one ( a, b, c, m, machep, e, t, h_02, &
    'h_02(x) = x * x' )

  a = -0.5D+00
  b =  +2.0D+00
  c = ( a + b ) / 2.0D+00
  m = 14.0D+00

  call test_glomin_one ( a, b, c, m, machep, e, t, h_03, &
    'h_03(x) = x^3 + x^2' )

  a = -0.5D+00
  b =  +2.0D+00
  c = ( a + b ) / 2.0D+00
  m = 28.0D+00

  call test_glomin_one ( a, b, c, m, machep, e, t, h_03, &
    'h_03(x) = x^3 + x^2' )

  a = -10.0D+00
  b = +10.0D+00
  c = ( a + b ) / 2.0D+00
  m = 72.0D+00

  call test_glomin_one ( a, b, c, m, machep, e, t, h_04, &
    'h_04(x) = ( x + sin(x) ) * exp(-x*x)' )

  a = -10.0D+00
  b = +10.0D+00
  c = ( a + b ) / 2.0D+00
  m = 72.0D+00

  call test_glomin_one ( a, b, c, m, machep, e, t, h_05, &
    'h_05(x) = ( x - sin(x) ) * exp(-x*x)' )

  return
end
subroutine test_zero_one ( a, b, machep, t, f, title )

!*****************************************************************************80
!
!! TEST_ZERO_ONE tests Brent's zero finding routine on one test function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the two endpoints of the change of sign
!    interval.
!
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Input, character ( LEN = * ) TITLE, a title for the problem.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fz
  real ( kind = 8 ) machep
  real ( kind = 8 ) t
  character ( len = *  ) title
  real ( kind = 8 ) z
  real ( kind = 8 ) zero
  
  z = zero ( a, b, machep, t, f )
  fz = f ( z )
  fa = f ( a )
  fb = f ( b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A                 Z             B'
  write ( *, '(a)' ) '    F(A)              F(Z)          F(B)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,f14.8,2x,f14.8,2x,f14.8)' ) a,  z,  b
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) fa, fz, fb

  return
end
subroutine test_zero_rc_one ( a, b, machep, t, f, title )

!*****************************************************************************80
!
!! TEST_ZERO_RC_ONE tests Brent's zero finding routine on one test function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the two endpoints of the change of sign
!    interval.
!
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Input, character ( LEN = * ) TITLE, a title for the problem.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) machep
  integer ( kind = 4 ) status
  real ( kind = 8 ) t
  character ( len = *  ) title
  real ( kind = 8 ) value
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STATUS      X               F(X)'
  write ( *, '(a)' ) ' '

  status = 0

  do 

    call zero_rc ( a, b, t, arg, status, value )

    if ( status < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ZERO_RC returned an error flag!'
      exit
    end if

    value = f ( arg )

    write ( *, '(2x,i8,2x,g14.8,2x,g14.8)' ) status, arg, value

    if ( status == 0 ) then
      exit 
    end if

  end do

  return
end
subroutine test_local_min_one ( a, b, eps, t, f, title )

!*****************************************************************************80
!
!! TEST_LOCAL_MIN_ONE tests Brent's local minimizer on one test function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, real ( kind = 8 ) EPS, a positive relative error tolerance.
!
!    Input, real ( kind = 8 ) T, a positive absolute error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    Input, character ( LEN = * ) TITLE, a title for the problem.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fx
  real ( kind = 8 ) local_min
  real ( kind = 8 ) t
  character ( len = * ) title
  real ( kind = 8 ) x

  fx = local_min ( a, b, eps, t, f, x )
  fa = f ( a )
  fb = f ( b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A                 X             B'
  write ( *, '(a)' ) '    F(A)              F(X)          F(B)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,f14.8,2x,f14.8,2x,f14.8)' ) a,  x,  b
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) fa, fx, fb

  return
end
subroutine test_local_min_rc_one ( a, b, t, f, title )

!*****************************************************************************80
!
!! TEST_LOCAL_MIN_RC_ONE tests LOCAL_MIN_RC on one test function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, real ( kind = 8 ) T, a positive absolute error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    Input, character ( LEN = * ) TITLE, a title for the problem.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fx
  integer ( kind = 4 ) status
  integer ( kind = 4 ) step
  real ( kind = 8 ) t
  character ( len = * ) title
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step      X                          F(X)'
  write ( *, '(a)' ) ' '
  step = 0

  arg = a
  value = f ( arg )
  write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) step, arg, value

  arg = b
  value = f ( arg )
  write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) step, arg, value

  a2 = a
  b2 = b
  status = 0

  do

    call local_min_rc ( a2, b2, arg, status, value )
 
    if ( status < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_LOCAL_MIN_RC_ONE - Fatal error!'
      write ( *, '(a)' ) '  LOCAL_MIN_RC returned negative status.'
      exit
    end if

    value = f ( arg )

    step = step + 1
    write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) step, arg, value

    if ( status == 0 ) then
      exit
    end if 

  end do

  return
end
subroutine test_glomin_one ( a, b, c, m, machep, e, t, f, title )

!*****************************************************************************80
!
!! TEST_GLOMIN_ONE tests the Brent global minimizer on one test function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, real ( kind = 8 ) C, an initial guess for the global
!    minimizer.  If no good guess is known, C = A or B is acceptable.
!
!    Input, real ( kind = 8 ) M, the bound on the second derivative.
!
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real ( kind = 8 ) E, a positive tolerance, a bound for the
!    absolute error in the evaluation of F(X) for any X in [A,B].
!
!    Input, real ( kind = 8 ) T, a positive absolute error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose global minimum is being sought.
!
!    Input, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) e
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fx
  real ( kind = 8 ) glomin
  real ( kind = 8 ) m
  real ( kind = 8 ) machep
  real ( kind = 8 ) t
  character ( len = * ) title
  real ( kind = 8 ) x

  fx = glomin ( a, b, c, m, machep, e, t, f, x )
  fa = f ( a )
  fb = f ( b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A                 X             B'
  write ( *, '(a)' ) '    F(A)              F(X)          F(B)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,f14.8,2x,f14.8,2x,f14.8)' ) a,  x,  b
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) fa, fx, fb

  return
end
function f_01 ( x )

!*****************************************************************************80
!
!! F_01 evaluates sin ( x ) - x / 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_01, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_01
  real ( kind = 8 ) x

  f_01 = sin ( x ) - 0.5D+00 * x

  return
end
function f_02 ( x )

!*****************************************************************************80
!
!! F_02 evaluates 2*x-exp(-x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_02, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_02
  real ( kind = 8 ) x

  f_02 = 2.0D+00 * x - exp ( - x )

  return
end
function f_03 ( x )

!*****************************************************************************80
!
!! F_03 evaluates x*exp(-x).
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_03, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_03
  real ( kind = 8 ) x

  f_03 = x * exp ( - x )

  return
end
function f_04 ( x )

!*****************************************************************************80
!
!! F_04 evaluates exp(x) - 1 / (100*x*x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_04, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_04
  real ( kind = 8 ) x

  f_04 = exp ( x ) - 1.0D+00 / 100.0D+00 / x / x

  return
end
function f_05 ( x )

!*****************************************************************************80
!
!! F_05 evaluates (x+3)*(x-1)*(x-1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_05, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_05
  real ( kind = 8 ) x

  f_05 = ( x + 3.0D+00 ) * ( x - 1.0D+00 ) * ( x - 1.0D+00 )

  return
end
function g_01 ( x )

!*****************************************************************************80
!
!! G_01 evaluates (x-2)^2 + 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) G_01, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) g_01
  real ( kind = 8 ) x

  g_01 = ( x - 2.0D+00 ) * ( x - 2.0D+00 ) + 1.0D+00

  return
end
function g_02 ( x )

!*****************************************************************************80
!
!! G_02 evaluates x^2 + exp ( - x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) G_02, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) g_02
  real ( kind = 8 ) x

  g_02 = x * x + exp ( - x )

  return
end
function g_03 ( x )

!*****************************************************************************80
!
!! G_03 evaluates x^4+2x^2+x+3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) G_03, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) g_03
  real ( kind = 8 ) x

  g_03 = ( ( x * x + 2.0D+00 ) * x + 1.0D+00 ) * x + 3.0D+00

  return
end
function g_04 ( x )

!*****************************************************************************80
!
!! G_04 evaluates exp(x)+1/(100X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) G_04, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) g_04
  real ( kind = 8 ) x

  g_04 = exp ( x ) + 0.01D+00 / x

  return
end
function g_05 ( x )

!*****************************************************************************80
!
!! G_05 evaluates exp(x) - 2x + 1/(100x) - 1/(1000000x^2)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) G_05, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) g_05
  real ( kind = 8 ) x

  g_05 = exp ( x ) - 2.0D+00 * x + 0.01D+00 / x - 0.000001D+00 / x / x

  return
end
function h_01 ( x )

!*****************************************************************************80
!
!! H_01 evaluates 2 - x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) H_01, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) h_01
  real ( kind = 8 ) x

  h_01 = 2.0D+00 - x

  return
end
function h_02 ( x )

!*****************************************************************************80
!
!! H_02 evaluates x^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) H_02, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) h_02
  real ( kind = 8 ) x

  h_02 = x * x

  return
end
function h_03 ( x )

!*****************************************************************************80
!
!! H_03 evaluates x^3+x^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) H_03, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) h_03
  real ( kind = 8 ) x

  h_03 = x * x * ( x + 1.0D+00 )

  return
end
function h_04 ( x )

!*****************************************************************************80
!
!! H_04 evaluates ( x + sin ( x ) ) * exp ( - x * x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) H_04, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) h_04
  real ( kind = 8 ) x

  h_04 = ( x + sin ( x ) ) * exp ( - x * x )

  return
end
function h_05 ( x )

!*****************************************************************************80
!
!! H_05 evaluates ( x - sin ( x ) ) * exp ( - x * x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) H_05, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) h_05
  real ( kind = 8 ) x

  h_05 = ( x - sin ( x ) ) * exp ( - x * x )

  return
end
