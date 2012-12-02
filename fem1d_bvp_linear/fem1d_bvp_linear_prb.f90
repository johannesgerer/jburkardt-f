program main

!*****************************************************************************80
!
!! FEM1D_BVP_LINEAR_PRB tests the routines in FEM1D_BVP_LINEAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FEM1D_BVP_LINEAR library.'

  call fem1d_bvp_linear_test01 ( )
  call fem1d_bvp_linear_test02 ( )
  call fem1d_bvp_linear_test03 ( )
  call fem1d_bvp_linear_test04 ( )
  call fem1d_bvp_linear_test05 ( )
  call fem1d_bvp_linear_test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine fem1d_bvp_linear_test01 ( )

!*****************************************************************************80
!
!! FEM1D_BVP_LINEAR_TEST01 carries out test case #1.
!
!  Discussion:
!
!    Use A1, C1, F1, EXACT1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ), external :: a1
  real ( kind = 8 ), external :: c1
  real ( kind = 8 ), external :: exact1
  real ( kind = 8 ), external :: exact_ux1
  real ( kind = 8 ), external :: f1
  integer ( kind = 4 ) i
  real ( kind = 8 ) l2_error
  real ( kind = 8 ) seminorm_error
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) uexact
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_first
  real ( kind = 8 ) x_last

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST01'
  write ( *, '(a)' ) '  A1(X)  = 1.0'
  write ( *, '(a)' ) '  C1(X)  = 0.0'
  write ( *, '(a)' ) '  F1(X)  = X * ( X + 3 ) * exp ( X )'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
!
!  Geometry definitions.
!
  x_first = 0.0D+00
  x_last = 1.0D+00
  call r8vec_even ( n, x_first, x_last, x )

  call fem1d_bvp_linear ( n, a1, c1, f1, x, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    uexact = exact1 ( x(i) )
    write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, x(i), u(i), uexact, abs ( u(i) - uexact )
  end do

  call compute_l2_error ( n, x, u, exact1, l2_error )
  call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
  write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

  return
end
subroutine fem1d_bvp_linear_test02 ( )

!*****************************************************************************80
!
!! FEM1D_BVP_LINEAR_TEST02 carries out test case #2.
!
!  Discussion:
!
!    Use A1, C2, F2, EXACT1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ), external :: a1
  real ( kind = 8 ), external :: c2
  real ( kind = 8 ), external :: exact1
  real ( kind = 8 ), external :: exact_ux1
  real ( kind = 8 ), external :: f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) l2_error
  real ( kind = 8 ) seminorm_error
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) uexact
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_first
  real ( kind = 8 ) x_last

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST02'
  write ( *, '(a)' ) '  A1(X)  = 1.0'
  write ( *, '(a)' ) '  C2(X)  = 2.0'
  write ( *, '(a)' ) '  F2(X)  = X * ( 5 - X ) * exp ( X )'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
!
!  Geometry definitions.
!
  x_first = 0.0D+00
  x_last = 1.0D+00
  call r8vec_even ( n, x_first, x_last, x )

  call fem1d_bvp_linear ( n, a1, c2, f2, x, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    uexact = exact1 ( x(i) )
    write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, x(i), u(i), uexact, abs ( u(i) - uexact )
  end do

  call compute_l2_error ( n, x, u, exact1, l2_error )
  call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
  write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

  return
end
subroutine fem1d_bvp_linear_test03 ( )

!*****************************************************************************80
!
!! FEM1D_BVP_LINEAR_TEST03 carries out test case #3.
!
!  Discussion:
!
!    Use A1, C3, F3, EXACT1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ), external :: a1
  real ( kind = 8 ), external :: c3
  real ( kind = 8 ), external :: exact1
  real ( kind = 8 ), external :: exact_ux1
  real ( kind = 8 ), external :: f3
  integer ( kind = 4 ) i
  real ( kind = 8 ) l2_error
  real ( kind = 8 ) seminorm_error
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) uexact
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_first
  real ( kind = 8 ) x_last

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST03'
  write ( *, '(a)' ) '  A1(X)  = 1.0'
  write ( *, '(a)' ) '  C3(X)  = 2.0 * X'
  write ( *, '(a)' ) '  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
!
!  Geometry definitions.
!
  x_first = 0.0D+00
  x_last = 1.0D+00
  call r8vec_even ( n, x_first, x_last, x )

  call fem1d_bvp_linear ( n, a1, c3, f3, x, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    uexact = exact1 ( x(i) )
    write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, x(i), u(i), uexact, abs ( u(i) - uexact )
  end do

  call compute_l2_error ( n, x, u, exact1, l2_error )
  call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
  write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

  return
end
subroutine fem1d_bvp_linear_test04 ( )

!*****************************************************************************80
!
!! FEM1D_BVP_LINEAR_TEST04 carries out test case #4.
!
!  Discussion:
!
!    Use A2, C1, F4, EXACT1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ), external :: a2
  real ( kind = 8 ), external :: c1
  real ( kind = 8 ), external :: exact1
  real ( kind = 8 ), external :: exact_ux1
  real ( kind = 8 ), external :: f4
  integer ( kind = 4 ) i
  real ( kind = 8 ) l2_error
  real ( kind = 8 ) seminorm_error
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) uexact
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_first
  real ( kind = 8 ) x_last

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST04'
  write ( *, '(a)' ) '  A2(X)  = 1.0 + X * X'
  write ( *, '(a)' ) '  C1(X)  = 0.0'
  write ( *, '(a)' ) '  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
!
!  Geometry definitions.
!
  x_first = 0.0D+00
  x_last = 1.0D+00
  call r8vec_even ( n, x_first, x_last, x )

  call fem1d_bvp_linear ( n, a2, c1, f4, x, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    uexact = exact1 ( x(i) )
    write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, x(i), u(i), uexact, abs ( u(i) - uexact )
  end do

  call compute_l2_error ( n, x, u, exact1, l2_error )
  call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
  write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

  return
end
subroutine fem1d_bvp_linear_test05 ( )

!*****************************************************************************80
!
!! FEM1D_BVP_LINEAR_TEST05 carries out test case #5.
!
!  Discussion:
!
!    Use A3, C1, F5, EXACT1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ), external :: a3
  real ( kind = 8 ), external :: c1
  real ( kind = 8 ), external :: exact1
  real ( kind = 8 ), external :: exact_ux1
  real ( kind = 8 ), external :: f5
  integer ( kind = 4 ) i
  real ( kind = 8 ) l2_error
  real ( kind = 8 ) seminorm_error
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) uexact
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_first
  real ( kind = 8 ) x_last

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST05'
  write ( *, '(a)' ) '  A3(X)  = 1.0 + X * X for X <= 1/3'
  write ( *, '(a)' ) '         = 7/9 + X     for      1/3 < X'
  write ( *, '(a)' ) '  C1(X)  = 0.0'
  write ( *, '(a)' ) '  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )'
  write ( *, '(a)' ) '                       for X <= 1/3'
  write ( *, '(a)' ) '         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) .* exp ( X )'
  write ( *, '(a)' ) '                       for      1/3 <= X'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
!
!  Geometry definitions.
!
  x_first = 0.0D+00
  x_last = 1.0D+00
  call r8vec_even ( n, x_first, x_last, x )

  call fem1d_bvp_linear ( n, a3, c1, f5, x, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    uexact = exact1 ( x(i) )
    write ( *, '(2x,i4,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, x(i), u(i), uexact, abs ( u(i) - uexact )
  end do

  call compute_l2_error ( n, x, u, exact1, l2_error )
  call compute_seminorm_error ( n, x, u, exact_ux1, seminorm_error )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', l2_error
  write ( *, '(a,g14.6)' ) '  Seminorm of error = ', seminorm_error

  return
end
subroutine fem1d_bvp_linear_test06 ( )

!*****************************************************************************80
!
!! FEM1D_BVP_LINEAR_TEST06 does an error analysis.
!
!  Discussion:
!
!    Use A1, C1, F6, EXACT4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), external :: a1
  real ( kind = 8 ), external :: c1
  real ( kind = 8 ), external :: exact4
  real ( kind = 8 ), external :: exact_ux4
  real ( kind = 8 ), external :: f6
  integer ( kind = 4 ) i
  real ( kind = 8 ) l2_error
  integer ( kind = 4 ) n
  real ( kind = 8 ) seminorm_error
  real ( kind = 8 ), allocatable :: u(:)
  real ( kind = 8 ) uexact
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x_first
  real ( kind = 8 ) x_last

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_BVP_LINEAR_TEST06'
  write ( *, '(a)' ) '  A1(X)  = 1.0 '
  write ( *, '(a)' ) '  C1(X)  = 0.0'
  write ( *, '(a)' ) '  F6(X)  = pi*pi*sin(pi*X)'
  write ( *, '(a)' ) '  U4(X)  = sin(pi*x)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Compute L2 norm and seminorm of error for various N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N        L2 error      Seminorm error'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '

  n = 11

  do i = 0, 4
!
!  Geometry definitions.
!
    x_first = 0.0D+00
    x_last = 1.0D+00
    allocate ( x(1:n) )
    call r8vec_even ( n, x_first, x_last, x )

    allocate ( u(1:n) )
    call fem1d_bvp_linear ( n, a1, c1, f6, x, u )

    call compute_l2_error ( n, x, u, exact4, l2_error )
    call compute_seminorm_error ( n, x, u, exact_ux4, seminorm_error )

    write ( *, '(2x,i4,2x,f14.6,2x,f14.6)' ) n, l2_error, seminorm_error

    deallocate ( u )
    deallocate ( x )

    n = 2 * ( n - 1 ) + 1

  end do

  return
end
function a1 ( x )

!*****************************************************************************80
!
!! A1 evaluates A function #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) A1, the value of A(X).
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) x

  a1 = 1.0D+00

  return
end
function a2 ( x )

!*****************************************************************************80
!
!! A2 evaluates A function #2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) A2, the value of A(X).
!
  implicit none

  real ( kind = 8 ) a2
  real ( kind = 8 ) x

  a2 = 1.0D+00 + x * x

  return
end
function a3 ( x )

!*****************************************************************************80
!
!! A3 evaluates A function #3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) A3, the value of A(X).
!
  implicit none

  real ( kind = 8 ) a3
  real ( kind = 8 ) x

  if ( x <= 1.0D+00 / 3.0D+00 ) then
    a3 = 1.0D+00 + x * x
  else
    a3 = x + 7.0D+00 / 9.0D+00
  end if

  return
end
function c1 ( x )

!*****************************************************************************80
!
!! C1 evaluates C function #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) C1, the value of C(X).
!
  implicit none

  real ( kind = 8 ) c1
  real ( kind = 8 ) x

  c1 = 0.0D+00

  return
end
function c2 ( x )

!*****************************************************************************80
!
!! C2 evaluates C function #2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) C2, the value of C(X).
!
  implicit none

  real ( kind = 8 ) c2
  real ( kind = 8 ) x

  c2 = 2.0D+00

  return
end
function c3 ( x )

!*****************************************************************************80
!
!! C3 evaluates C function #3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) C3, the value of C(X).
!
  implicit none

  real ( kind = 8 ) c3
  real ( kind = 8 ) x

  c3 = 2.0D+00 * x

  return
end
function exact1 ( x )

!*****************************************************************************80
!
!! EXACT1 evaluates exact solution #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) EXACT1, the value of U(X).
!
  implicit none

  real ( kind = 8 ) exact1
  real ( kind = 8 ) x

  exact1 = x * ( 1.0D+00 - x ) * exp ( x )

  return
end
function exact_ux1 ( x )

!*****************************************************************************80
!
!! EXACT_UX1 evaluates the derivative of exact solution #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) EXACT_UX1, the value of dUdX(X).
!
  implicit none

  real ( kind = 8 ) exact_ux1
  real ( kind = 8 ) x

  exact_ux1 = ( 1.0D+00 - x - x * x ) * exp ( x )

  return
end
function exact2 ( x )

!*****************************************************************************80
!
!! EXACT2 returns exact solution #2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) EXACT2, the value of U(X).
!
  implicit none

  real ( kind = 8 ) exact2
  real ( kind = 8 ) x

  if ( x <= 2.0D+00 / 3.0D+00 ) then
    exact2 =  x * ( 1.0D+00 - x ) * exp ( x )
  else
    exact2 = x * ( 1.0D+00 - x )  * exp ( 2.0D+00 / 3.0D+00 )
  end if

  return
end
function exact3 ( x )

!*****************************************************************************80
!
!! EXACT3 returns exact solution #3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) EXACT3, the value of U(X).
!
  implicit none

  real ( kind = 8 ) exact3
  real ( kind = 8 ) x

  if ( x <= 2.0D+00 / 3.0D+00 ) then
    exact3 = x * ( 1.0D+00 - x ) * exp ( x )
  else
    exact3 = x * ( 1.0D+00 - x )
  end if

  return
end
function exact4 ( x )

!*****************************************************************************80
!
!! EXACT4 evaluates exact solution #4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) EXACT4, the value of U(X).
!
  implicit none

  real ( kind = 8 ) exact4
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  exact4 = sin ( pi * x )

  return
end
function exact_ux4 ( x )

!*****************************************************************************80
!
!! EXACT_UX4 evaluates the derivative of exact solution #4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) EXACT_UX4, the value of dUdX(X).
!
  implicit none

  real ( kind = 8 ) exact_ux4
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  exact_ux4 = pi * cos ( pi * x )

  return
end
function f1 ( x )

!*****************************************************************************80
!
!! F1 evaluates right hand side function #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) F1, the value of F(X).
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) x

  f1 = x * ( x + 3.0D+00 ) * exp ( x )

  return
end
function f2 ( x )

!*****************************************************************************80
!
!! F2 evaluates right hand side function #2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) F2, the value of F(X).
!
  implicit none

  real ( kind = 8 ) f2
  real ( kind = 8 ) x

  f2 = x * ( 5.0D+00 - x ) * exp ( x )

  return
end
function f3 ( x )

!*****************************************************************************80
!
!! F3 evaluates right hand side function #3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) F3, the value of F(X).
!
  implicit none

  real ( kind = 8 ) f3
  real ( kind = 8 ) x

  f3 = - x * ( 2.0D+00 * x * x - 3.0D+00 * x - 3.0D+00 ) * exp ( x )

  return
end
function f4 ( x )

!*****************************************************************************80
!
!! F4 evaluates right hand side function #4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) F4, the value of F(X).
!
  implicit none

  real ( kind = 8 ) f4
  real ( kind = 8 ) x

  f4 = ( x + 3.0D+00 * x * x + 5.0D+00 * x * x * x + x * x * x * x ) * exp ( x )

  return
end
function f5 ( x )

!*****************************************************************************80
!
!! F5 evaluates right hand side function #5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) F5, the value of F(X).
!
  implicit none

  real ( kind = 8 ) f5
  real ( kind = 8 ) x

  if ( x <= 1.0D+00 / 3.0D+00 ) then
    f5 = ( x + 3.0D+00 * x * x + 5.0D+00 * x * x * x + x * x * x * x ) &
      * exp ( x )
  else
    f5 = ( - 1.0D+00 + ( 10.0D+00 / 3.0D+00 ) * x &
      + ( 43.0D+00 / 9.0D+00 ) * x * x + x * x * x ) * exp ( x )
  end if

  return
end
function f6 ( x )

!*****************************************************************************80
!
!! F6 evaluates right hand side function #6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) F6, the value of F(X).
!
  implicit none

  real ( kind = 8 ) f6
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  f6 = pi * pi * sin ( pi * x )

  return
end
