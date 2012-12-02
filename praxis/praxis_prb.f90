program main

!*****************************************************************************80
!
!! MAIN is the main program for PRAXIS_PRB.
!
!  Discussion:
!
!    PRAXIS_PRB tests PRAXIS.
!
!  Modified:
!
!    11 February 2007
!
!  Author:
!
!    John Burkardt
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRAXIS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PRAXIS library.'

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
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRAXIS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 calls PRAXIS for the Beale function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ), external :: f_01
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  The Beale function.'

  t0 = 0.00001D+00
  h0 = 0.25D+00
  prin = 0

  x(1:n) = (/ 0.1D+00, 0.1D+00 /)

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_01 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_01 )

  call r8vec_print ( n, x, '  Computed minimizer:' )
  
  write ( *, '(a,g14.6)' ) '  Function value = ', f_01 ( x, n )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 calls PRAXIS for the Box function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), external :: f_02
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  The Box function.'

  t0 = 0.00001D+00
  h0 = 20.0D+00
  prin = 0

  x(1:n) = (/ 0.0D+00, 10.0D+00, 20.0D+00 /)

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_02 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_02 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_02 ( x, n )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 calls PRAXIS for the Chebyquad function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 8

  real ( kind = 8 ), external :: f_03
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  The Chebyquad function.'

  t0 = 0.00001D+00
  h0 = 0.1D+00
  prin = 0

  do i = 1, n
    x(i) = real ( i ) / real ( n + 1 )
  end do

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_03 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_03 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_03 ( x, n )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 calls PRAXIS for the Cube function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ), external :: f_04
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  The Cube function.'

  t0 = 0.00001D+00
  h0 = 1.0D+00
  prin = 0

  x(1:n) = (/ -1.2D+00, -1.0D+00 /)

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_04 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_04 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_04 ( x, n )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 calls PRAXIS for the Fletcher-Powell Helix function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), external :: f_05
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  The Fletcher-Powell Helix function.'

  t0 = 0.00001D+00
  h0 = 1.0D+00
  prin = 0

  x(1:n) = (/ -1.0D+00, 0.0D+00, 0.0D+00 /)

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_05 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_05 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_05 ( x, n )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 calls PRAXIS for the Hilbert function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ), external :: f_06
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  The Hilbert function.'

  t0 = 0.00001D+00
  h0 = 10.0D+00
  prin = 0

  x(1:n) = 1.0D+00

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_06 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_06 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_06 ( x, n )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 calls PRAXIS for the Powell 3D function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), external :: f_07
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  The Powell 3D function.'

  t0 = 0.00001D+00
  h0 = 1.0D+00
  prin = 0

  x(1:n) = (/ 0.0D+00, 1.0D+00, 2.0D+00 /)

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_07 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_07 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_07 ( x, n )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 calls PRAXIS for the Rosenbrock function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ), external :: f_08
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  The Rosenbrock function.'

  t0 = 0.00001D+00
  h0 = 1.0D+00
  prin = 0

  x(1:n) = ( -1.2D+00, 1.0D+00 )

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_08 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_08 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_08 ( x, n )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 calls PRAXIS for the Powell Singular function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ), external :: f_09
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  The Powell Singular function.'

  t0 = 0.00001D+00
  h0 = 1.0D+00
  prin = 0

  x(1:n) = (/ 3.0D+00, -1.0D+00, 0.0D+00, 1.0D+00 /)

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_09 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_09 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_09 ( x, n )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 calls PRAXIS for the Tridiagonal function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ), external :: f_10
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  The Tridiagonal function.'

  t0 = 0.00001D+00
  h0 = 8.0D+00
  prin = 0

  x(1:n) = 0.0D+00

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_10 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_10 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_10 ( x, n )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 calls PRAXIS for the Watson function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ), external :: f_11
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  The Watson function.'

  t0 = 0.00001D+00
  h0 = 1.0D+00
  prin = 0

  x(1:n) = 0.0D+00

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_11 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_11 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_11 ( x, n )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 calls PRAXIS for the Wood function.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ), external :: f_12
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  real ( kind = 8 ) pr
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) t0
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  The Wood function.'

  t0 = 0.00001D+00
  h0 = 10.0D+00
  prin = 0

  x(1:n) = (/ -3.0D+00, -1.0D+00, -3.0D+00, -1.0D+00 /)

  call r8vec_print ( n, x, '  Initial point:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_12 ( x, n )

  pr = praxis ( t0, h0, n, prin, x, f_12 )

  call r8vec_print ( n, x, '  Computed minimizer:' )

  write ( *, '(a,g14.6)' ) '  Function value = ', f_12 ( x, n )

  return
end
function f_01 ( x, n )

!*****************************************************************************80
!
!! F_01 evaluates the Beale function.
!
!  Discussion:
!
!    The function is the sum of the squares of three functions.
!
!    This function has a valley approaching the line X(2) = 1.
!
!  Reference:
!
!    E Beale,
!    On an Iterative Method for Finding a Local Minimum of a Function
!    of More than One Variable,
!    Technical Report 25, Statistical Techniques Research Group,
!    Princeton University, 1958.
!
!    Richard Brent,
!    Algorithms for Finding Zeros and Extrema of Functions Without
!    Calculating Derivatives,
!    Stanford University Technical Report STAN-CS-71-198.
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_01, the value of the objective function.
!
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: c1 = 1.5D+00
  real ( kind = 8 ), parameter :: c2 = 2.25D+00
  real ( kind = 8 ), parameter :: c3 = 2.625D+00
  real ( kind = 8 ) f_01
  real ( kind = 8 ) x(n)

  f_01 =    ( c1 - x(1) * ( 1.0D+00 - x(2)    ) )**2 &
          + ( c2 - x(1) * ( 1.0D+00 - x(2)**2 ) )**2 &
          + ( c3 - x(1) * ( 1.0D+00 - x(2)**3 ) )**2

  return
end
function f_02 ( x, n )

!*****************************************************************************80
!
!! F_02 evaluates the Box function.
!
!  Discussion:
!
!    The function is formed by the sum of squares of 10 separate terms.
!
!  Modified:
!
!    04 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_02, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) f_02
  real ( kind = 8 ) fi
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  f_02 = 0.0D+00

  do i = 1, 10

    c = - real ( i ) / 10.0D+00

    fi = exp ( c * x(1) ) - exp ( c * x(2) ) - x(3) * &
      ( exp ( c ) - exp ( 10.0D+00 * c ) )
   
    f_02 = f_02 + fi**2

  end do

  return
end
function f_03 ( x, n )

!*****************************************************************************80
!
!! F_03 evaluates the Chebyquad function.
!
!  Discussion:
!
!    The function is formed by the sum of squares of N separate terms.
!
!  Modified:
!
!    26 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_03, the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_03
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) th
  real ( kind = 8 ) x(n)

  fvec(1:n) = 0.0D+00

  do j = 1, n
    t1 = 1.0D+00
    t2 = 2.0D+00 * x(j) - 1.0D+00
    t = 2.0D+00 * t2
    do i = 1, n
      fvec(i) = fvec(i) + t2
      th = t * t2 - t1
      t1 = t2
      t2 = th
    end do
  end do

  do i = 1, n
    fvec(i) = fvec(i) / real ( n )
    if ( mod ( i, 2 ) == 0 ) then
      fvec(i) = fvec(i) + 1.0D+00 / ( real ( i )**2 - 1.0D+00 )
    end if
  end do
!
!  Compute F.
!
  f_03 = sum ( fvec(1:n)**2 )

  return
end
function f_04 ( x, n )

!*****************************************************************************80
!
!! F_04 evaluates the Cube function.
!
!  Discussion:
!
!    The function is the sum of the squares of two functions.
!
!  Modified:
!
!    27 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_04, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_04
  real ( kind = 8 ) x(n)

  f_04 = ( 10.0D+00 * ( x(2) - x(1)**3 ) )**2 + ( 1.0D+00 - x(1) )**2

  return
end
function f_05 ( x, n )

!*****************************************************************************80
!
!! F_05 evaluates the Helix function.
!
!  Discussion:
!
!    The function is the sum of the squares of three functions.
!
!  Modified:
!
!    27 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_05, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_05
  real ( kind = 8 ), parameter :: pi = 3.14159265D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)

  r = sqrt ( x(1)**2 + x(2)**2 )

  if ( x(1) >= 0.0D+00 ) then
    theta = 0.5D+00 * atan2 ( x(2), x(1) ) / pi
  else if ( x(1) < 0.0D+00 ) then
    theta = 0.5D+00 * ( atan2 ( x(2), x(1) ) + pi ) / pi
  end if

  f_05 = &
      ( 10.0D+00 * ( x(3) - 10.0D+00 * theta ) )**2 &
    + ( 10.0D+00 * ( r - 1.0D+00 ) )**2 &
    + x(3)**2

  return
end
function f_06 ( x, n )

!*****************************************************************************80
!
!! F_06 evaluates the Hilbert function.
!
!  Discussion:
!
!    The function is a positive definite quadratic function of 
!    the form
!
!      f(x) = x' A x
!
!    where A is the Hilbert matrix.
!
!  Modified:
!
!    27 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_06, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_06
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  f_06 = 0.0D+00

  do i = 1, n
    do j = 1, n
      f_06 = f_06 + x(i) * x(j) / real ( i + j - 1 )
    end do
  end do

  return
end
function f_07 ( x, n )

!*****************************************************************************80
!
!! F_07 evaluates the Powell 3D function.
!
!  Reference:
!
!    M J D Powell,
!    An Efficient Method for Finding the Minimum of a Function of
!      Several Variables Without Calculating Derivatives,
!    Computer Journal, 
!    Volume 7, Number 2, pages 155-162, 1964.
!    
!  Modified:
!
!    03 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_07, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_07
  real ( kind = 8 ), parameter :: pi = 3.14159265D+00
  real ( kind = 8 ) x(n)

  f_07 = 3.0D+00 - 1.0D+00 / ( 1.0D+00 + ( x(1) - x(2) )**2 ) &
    - sin ( 0.5D+00 * pi * x(2) * x(3) ) &
    - exp ( - ( ( x(1) - 2.0D+00 * x(2) + x(3) ) / x(2) )**2 )

  return
end
function f_08 ( x, n )

!*****************************************************************************80
!
!! F_08 evaluates the Rosenbrock function.
!
!  Modified:
!
!    01 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_08, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_08
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  f_08 = 0.0D+00
  do j = 1, n
    if ( mod ( j, 2 ) == 1 ) then
      f_08 = f_08 + ( 1.0D+00 - x(j) )**2
    else
      f_08 = f_08 + 100.0D+00 * ( x(j) - x(j-1)**2 )**2
    end if
  end do

  return
end
function f_09 ( x, n )

!*****************************************************************************80
!
!! F_09 evaluates the Powell Singular function.
!
!  Modified:
!
!    01 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_09, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_09
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xjp1
  real ( kind = 8 ) xjp2
  real ( kind = 8 ) xjp3

  f_09 = 0.0D+00

  do j = 1, n, 4

    if ( j + 1 <= n ) then
      xjp1 = x(j+1)
    else
      xjp1 = 0.0D+00
    end if

    if ( j + 2 <= n ) then
      xjp2 = x(j+2)
    else
      xjp2 = 0.0D+00
    end if

    if ( j + 3 <= n ) then
      xjp3 = x(j+3)
    else
      xjp3 = 0.0D+00
    end if
 
    f1 = x(j) + 10.0D+00 * xjp1
    if ( j + 1 <= n ) then
      f2 = xjp2 - xjp3
    else
      f2 = 0.0D+00
    end if
    if ( j + 2 <= n ) then
      f3 = xjp1 - 2.0D+00 * xjp2
    else
      f3 = 0.0D+00
    end if
    if ( j + 3 <= n ) then
      f4 = x(j) - xjp3
    else
      f4 = 0.0D+00
    end if

    f_09 = f_09 + f1**2 + 5.0D+00 * f2**2 + f3**4 + 10.0D+00 * f4**4

  end do

  return
end
function f_10 ( x, n )

!*****************************************************************************80
!
!! F_10 evaluates the tridiagonal function.
!
!  Modified:
!
!    01 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_10, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_10
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  f_10 = x(1)**2 + 2.0D+00 * sum ( x(2:n)**2 )

  do i = 1, n-1
    f_10 = f_10 - 2.0D+00 * x(i) * x(i+1)
  end do

  f_10 = f_10 - 2.0D+00 * x(1)

  return
end
function f_11 ( x, n )

!*****************************************************************************80
!
!! F_11 evaluates the Watson function.
!
!  Modified:
!
!    01 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_11, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  real ( kind = 8 ) f_11
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) x(n)

  f_11 = 0.0D+00
  do i = 1, 29

    s1 = 0.0D+00
    d = 1.0D+00
    do j = 2, n
      s1 = s1 + real ( j - 1 ) * d * x(j)
      d = d * real ( i ) / 29.0D+00
    end do

    s2 = 0.0D+00
    d = 1.0D+00
    do j = 1, n
      s2 = s2 + d * x(j)
      d = d * real ( i ) / 29.0D+00
    end do

    f_11 = f_11 + ( s1 - s2**2 - 1.0D+00 )**2

  end do

  f_11 = f_11 + x(1)**2 + ( x(2) - x(1)**2 - 1.0D+00 )**2

  return
end
function f_12 ( x, n )

!*****************************************************************************80
!
!! F_12 evaluates the Wood function.
!
!  Modified:
!
!    01 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(N), the argument of the objection function.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Output, real ( kind = 8 ) F_12, the value of the objective function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f_12
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  real ( kind = 8 ) f5
  real ( kind = 8 ) f6
  real ( kind = 8 ) x(n)

  f1 = x(2) - x(1)**2
  f2 = 1.0D+00 - x(1)
  f3 = x(4) - x(3)**2
  f4 = 1.0D+00 - x(3)
  f5 = x(2) + x(4) - 2.0D+00
  f6 = x(2) - x(4)

  f_12 = 100.0D+00 * f1**2 + f2**2 + 90.0D+00 * f3**2 + f4**2 &
    + 10.0D+00 * f5**2 + 0.1D+00 * f6**2

  return
end
