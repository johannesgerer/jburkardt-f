program main

!*****************************************************************************80
!
!! MAIN is the main program for FD1D_BVP_PRB.
!
!  Discussion:
!
!    FD1D_BVP_PRB tests the routines in FD1D_BVP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BVP_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FD1D_BVP library.'

  call fd1d_bvp_test01 ( )
  call fd1d_bvp_test02 ( )
  call fd1d_bvp_test03 ( )
  call fd1d_bvp_test04 ( )
  call fd1d_bvp_test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BVP_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine fd1d_bvp_test01 ( )

!*****************************************************************************80
!
!! FD1D_BVP_TEST01 carries out test case #1.
!
!  Discussion:
!
!    Use A1, C1, F1, EXACT1.
!
!    Repeat computation using variable mesh spacing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21

  real ( kind = 8 ), external :: a1
  real ( kind = 8 ), external :: a1prime
  real ( kind = 8 ), external :: c1
  real ( kind = 8 ), external :: f1
  character ( len = 255 ) filename
  integer ( kind = 4 ) i
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) u2(2,n)
  real ( kind = 8 ) uexact(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x1 = 0.0D+00
  real ( kind = 8 ) :: x2 = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BVP_TEST01'
  write ( *, '(a)' ) '  A1(X)  = 1.0'
  write ( *, '(a)' ) '  A1''(X) = 0.0'
  write ( *, '(a)' ) '  C1(X)  = 0.0'
  write ( *, '(a)' ) '  F1(X)  = X * ( X + 3 ) * exp ( X )'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
  write ( *, '(a,g14.6)' ) '  X1 = ', x1
  write ( *, '(a,g14.6)' ) '  X2 = ', x2

  call r8vec_even ( n, x1, x2, x )

  call fd1d_bvp ( n, a1, a1prime, c1, f1, x, u )

  call exact1 ( n, x, uexact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) &
      i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat, using a nonuniform mesh.'

  call r8vec_even ( n, x1, x2, x )

  x(1:n) = sqrt ( x(1:n) )

  call fd1d_bvp ( n, a1, a1prime, c1, f1, x, u )

  call exact1 ( n, x, uexact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) &
      i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
  end do
!
!  Write the data to files.
!
  filename = 'fd1d_bvp_test01_nodes.txt'
  call r8mat_write ( filename, 1, n, x );

  u2(1,1:n) = u(1:n)
  u2(2,1:n) = uexact(1:n)

  filename = 'fd1d_bvp_test01_values.txt'
  call r8mat_write ( filename, 2, n, u2 )

  return
end
subroutine fd1d_bvp_test02 ( )

!*****************************************************************************80
!
!! FD1D_BVP_TEST02 carries out test case #2.
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
!    15 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ), external :: a1
  real ( kind = 8 ), external :: a1prime
  real ( kind = 8 ), external :: c2
  real ( kind = 8 ), external :: f2
  character ( len = 255 ) filename
  integer ( kind = 4 ) i
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) u2(2,n)
  real ( kind = 8 ) uexact(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x1 = 0.0D+00
  real ( kind = 8 ) :: x2 = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BVP_TEST02'
  write ( *, '(a)' ) '  A1(X)  = 1.0'
  write ( *, '(a)' ) '  A1''(X) = 0.0'
  write ( *, '(a)' ) '  C2(X)  = 2.0'
  write ( *, '(a)' ) '  F2(X)  = X * ( 5 - X ) * exp ( X )'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
  write ( *, '(a,g14.6)' ) '  X1 = ', x1
  write ( *, '(a,g14.6)' ) '  X2 = ', x2

  call r8vec_even ( n, x1, x2, x )

  call fd1d_bvp ( n, a1, a1prime, c2, f2, x, u )

  call exact1 ( n, x, uexact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) &
      i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
  end do
!
!  Write the data to files.
!
  filename = 'fd1d_bvp_test02_nodes.txt'
  call r8mat_write ( filename, 1, n, x );

  u2(1,1:n) = u(1:n)
  u2(2,1:n) = uexact(1:n)

  filename = 'fd1d_bvp_test02_values.txt'
  call r8mat_write ( filename, 2, n, u2 )

  return
end
subroutine fd1d_bvp_test03 ( )

!*****************************************************************************80
!
!! FD1D_BVP_TEST03 carries out test case #3.
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
!    15 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ), external :: a1
  real ( kind = 8 ), external :: a1prime
  real ( kind = 8 ), external :: c3
  real ( kind = 8 ), external :: f3
  character ( len = 255 ) filename
  integer ( kind = 4 ) i
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) u2(2,n)
  real ( kind = 8 ) uexact(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x1 = 0.0D+00
  real ( kind = 8 ) :: x2 = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BVP_TEST03'
  write ( *, '(a)' ) '  A1(X)  = 1.0'
  write ( *, '(a)' ) '  A1''(X) = 0.0'
  write ( *, '(a)' ) '  C3(X)  = 2.0 * X'
  write ( *, '(a)' ) '  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
  write ( *, '(a,g14.6)' ) '  X1 = ', x1
  write ( *, '(a,g14.6)' ) '  X2 = ', x2

  call r8vec_even ( n, x1, x2, x )

  call fd1d_bvp ( n, a1, a1prime, c3, f3, x, u )

  call exact1 ( n, x, uexact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) &
      i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
  end do
!
!  Write the data to files.
!
  filename = 'fd1d_bvp_test03_nodes.txt'
  call r8mat_write ( filename, 1, n, x );

  u2(1,1:n) = u(1:n)
  u2(2,1:n) = uexact(1:n)

  filename = 'fd1d_bvp_test03_values.txt'
  call r8mat_write ( filename, 2, n, u2 )

  return
end
subroutine fd1d_bvp_test04 ( )

!*****************************************************************************80
!
!! FD1D_BVP_TEST04 carries out test case #4.
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
!    15 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ), external :: a2
  real ( kind = 8 ), external :: a2prime
  real ( kind = 8 ), external :: c1
  real ( kind = 8 ), external :: f4
  character ( len = 255 ) filename
  integer ( kind = 4 ) i
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) u2(2,n)
  real ( kind = 8 ) uexact(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x1 = 0.0D+00
  real ( kind = 8 ) :: x2 = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BVP_TEST04'
  write ( *, '(a)' ) '  A2(X)  = 1.0 + X * X'
  write ( *, '(a)' ) '  A2''(X) = 2.0 * X'
  write ( *, '(a)' ) '  C1(X)  = 0.0'
  write ( *, '(a)' ) '  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
  write ( *, '(a,g14.6)' ) '  X1 = ', x1
  write ( *, '(a,g14.6)' ) '  X2 = ', x2

  call r8vec_even ( n, x1, x2, x )

  call fd1d_bvp ( n, a2, a2prime, c1, f4, x, u )

  call exact1 ( n, x, uexact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) &
      i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
  end do
!
!  Write the data to files.
!
  filename = 'fd1d_bvp_test04_nodes.txt'
  call r8mat_write ( filename, 1, n, x );

  u2(1,1:n) = u(1:n)
  u2(2,1:n) = uexact(1:n)

  filename = 'fd1d_bvp_test04_values.txt'
  call r8mat_write ( filename, 2, n, u2 )

  return
end
subroutine fd1d_bvp_test05 ( )

!*****************************************************************************80
!
!! FD1D_BVP_TEST05 carries out test case #5.
!
!  Discussion:
!
!    Use A3, C1, F5, EXACT1.
!
!    Also, run with a variable mesh spacing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real ( kind = 8 ), external :: a3
  real ( kind = 8 ), external :: a3prime
  real ( kind = 8 ), external :: c1
  real ( kind = 8 ), external :: f5
  character ( len = 255 ) filename
  integer ( kind = 4 ) i
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) u2(2,n)
  real ( kind = 8 ) uexact(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) :: x1 = 0.0D+00
  real ( kind = 8 ) :: x2 = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BVP_TEST05'
  write ( *, '(a)' ) '  A3(X)  = 1.0 + X * X for X <= 1/3'
  write ( *, '(a)' ) '         = 7/9 + X     for      1/3 < X'
  write ( *, '(a)' ) '  A3''(X) = 2.0 * X     for X <= 1/3'
  write ( *, '(a)' ) '           1           for      1/3 < X'
  write ( *, '(a)' ) '  C1(X)  = 0.0'
  write ( *, '(a)' ) '  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )'
  write ( *, '(a)' ) '                       for X <= 1/3'
  write ( *, '(a)' ) '         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) * exp ( X )'
  write ( *, '(a)' ) '                       for      1/3 <= X'
  write ( *, '(a)' ) '  U1(X)  = X * ( 1 - X ) * exp ( X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
  write ( *, '(a,g14.6)' ) '  X1 = ', x1
  write ( *, '(a,g14.6)' ) '  X2 = ', x2

  call r8vec_even ( n, x1, x2, x )

  call fd1d_bvp ( n, a3, a3prime, c1, f5, x, u )

  call exact1 ( n, x, uexact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    X         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) &
      i, x(i), u(i), uexact(i), abs ( u(i) - uexact(i) )
  end do
!
!  Write the data to files.
!
  filename = 'fd1d_bvp_test05_nodes.txt'
  call r8mat_write ( filename, 1, n, x );

  u2(1,1:n) = u(1:n)
  u2(2,1:n) = uexact(1:n)

  filename = 'fd1d_bvp_test05_values.txt'
  call r8mat_write ( filename, 2, n, u2 )

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
!    16 May 2009
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
function a1prime ( x )

!*****************************************************************************80
!
!! A1PRIME evaluates A' function #1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) A1PRIME, the value of A'(X).
!
  implicit none

  real ( kind = 8 ) a1prime
  real ( kind = 8 ) x

  a1prime = 0.0D+00

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
!    16 May 2009
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
function a2prime ( x )

!*****************************************************************************80
!
!! A2PRIME evaluates A' function #2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) A2PRIME, the value of A'(X).
!
  implicit none

  real ( kind = 8 ) a2prime
  real ( kind = 8 ) x

  a2prime = 2.0D+00 * x

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
!    16 May 2009
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
function a3prime ( x )

!*****************************************************************************80
!
!! A3PRIME evaluates A' function #3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) A3PRIME, the value of A'(X).
!
  implicit none

  real ( kind = 8 ) a3prime
  real ( kind = 8 ) x

  if ( x <= 1.0D+00 / 3.0D+00 ) then
    a3prime = 2.0D+00 * x
  else
    a3prime = 1.0D+00
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
!    16 May 2009
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
!    16 May 2009
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
!    16 May 2009
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
!    16 May 2009
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
!    16 May 2009
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
!    16 May 2009
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
!    16 May 2009
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

  f4 = ( x + 3.0D+00 * x**2 + 5.0D+00 * x**3 + x**4 ) * exp ( x )

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
!    16 May 2009
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
    f5 = ( x + 3.0D+00 * x**2 + 5.0D+00 * x**3 + x**4 ) * exp ( x )
  else
    f5 = ( - 1.0D+00 + ( 10.0D+00 / 3.0D+00 ) * x &
      + ( 43.0D+00 / 9.0D+00 ) * x**2 + x**3 ) * exp ( x )
  end if

  return
end
subroutine exact1 ( n, x, uexact )

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
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) UEXACT(N), the values of U(X(1:N)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) uexact(n)
  real ( kind = 8 ) x(n)

  do i = 1, n
    uexact(i) = x(i) * ( 1.0D+00 - x(i) ) * exp ( x(i) )
  end do

  return
end
subroutine exact2 ( n, x, uexact )

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
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) UEXACT(N), the values of U(X(1:N)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) uexact(n)
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) <= 2.0D+00 / 3.0D+00 ) then
      uexact(i) = x(i) * ( 1.0D+00 - x(i) ) * exp ( x(i) )
    else
      uexact(i) = x(i) * ( 1.0D+00 - x(i) ) * exp ( 2.0D+00 / 3.0D+00 )
    end if
  end do

  return
end
subroutine exact3 ( n, x, uexact )

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
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) UEXACT(N), the values of U(X(1:N)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) uexact(n)
  real ( kind = 8 ) x(n)

  do i = 1, n
    if ( x(i) <= 2.0D+00 / 3.0D+00 ) then
      uexact(i) = x(i) * ( 1.0D+00 - x(i) ) * exp ( x(i) )
    else
      uexact(i) = x(i) * ( 1.0D+00 - x(i) )
    end if
  end do

  return
end
