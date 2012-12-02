program main

!*****************************************************************************80
!
!! SINE_TRANSFORM_TEST tests SINE_TRANSFORM.
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

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINE_TRANSFORM_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the SINE_TRANSFORM library.'

  call sine_transform_test01 ( )
  call sine_transform_test02 ( )
  call sine_transform_test03 ( )
  call sine_transform_test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINE_TRANSFORM_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine sine_transform_test01 ( )

!*****************************************************************************80
!
!! SINE_TRANSFORM_TEST01 demonstrates that the transform is its own inverse.
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

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINE_TRANSFORM_TEST01:'
  write ( *, '(a)' ) '  SINE_TRANSFORM_DATA does a sine transform of data'
  write ( *, '(a)' ) '  defined by a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate that the transform is its own inverse.'
  write ( *, '(a)' ) '  Let R be a random N vector.'
  write ( *, '(a)' ) '  Let S be the transform of D.'
  write ( *, '(a)' ) '  Let T be the transform of E.'
  write ( *, '(a)' ) '  Then R and T will be equal.'

  call r8vec_uniform_01 ( n, seed, r )
  call sine_transform_data ( n, r, s )
  call sine_transform_data ( n, s, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I      R(I)        S(I)        T(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g10.4,2x,g10.4,2x,g10.4)' ) i, r(i), s(i), t(i)
  end do

  return
end
subroutine sine_transform_test02 ( )

!*****************************************************************************80
!
!! SINE_TRANSFORM_TEST02 uses the functional form of the routine.
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

  integer ( kind = 4 ), parameter :: n = 9

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) f1(n)
  real ( kind = 8 ) f2(n)
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  integer ( kind = 4 ) i
  real ( kind = 8 ), external :: poly5
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) x(n)

  a = 1.0D+00
  b = 3.0D+00
!
!  Evenly spaced points between A and B, but omitting
!  A and B themselves.
!
  do i = 1, n
    x(i) = ( real ( n - i + 1, kind = 8 ) * a   &
           + real (     i,     kind = 8 ) * b ) &
           / real ( n     + 1, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINE_TRANSFORM_TEST02:'
  write ( *, '(a)' ) '  SINE_TRANSFORM_FUNCTION does a sine transform of data'
  write ( *, '(a)' ) '  defined by a function F(X) evaluated at equally spaced'
  write ( *, '(a)' ) '  points in an interval [A,B].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate that the transform is its own inverse.'
  write ( *, '(a)' ) '  Let X(0:N+1) be N+2 equally spaced points in [A,B].'
  write ( *, '(a)' ) '  Let S be the transform of F(X(1:N)).'
  write ( *, '(a)' ) '  Let F1 be the linear interpolant of (A,F(A)), (B,F(B)).'
  write ( *, '(a)' ) '  Let F2 be the transform of S.'
  write ( *, '(a)' ) '  Then F(X(1:N)) = F1(X(1:N)) + F2(1:N).'

  call sine_transform_function ( n, a, b, poly5, s )

  fa = poly5 ( a )
  fb = poly5 ( b )
  f1(1:n) = ( ( b - x(1:n)     ) * fa   &
            + (     x(1:n) - a ) * fb ) &
            / ( b          - a )

  call sine_transform_data ( n, s, f2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I      X(I)      F(X(I))       ' // &
    'S           F1          F2          F1+F2'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
    0, a, poly5 ( a ), 0.0D+00, fa, 0.0D+00, fa

  do i = 1, n
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
        i, x(i), poly5 ( x(i) ), s(i), f1(i), f2(i), f1(i) + f2(i)
  end do

  write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
    n + 1, b, poly5 ( b ), 0.0D+00, fb, 0.0D+00, fb

  return
end
subroutine sine_transform_test03 ( )

!*****************************************************************************80
!
!! SINE_TRANSFORM_TEST03 evaluates the sine transform interpolant.
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

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: n2 = 1 + 2 * ( n + 1 )

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) f2(n2)
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  integer ( kind = 4 ) i
  real ( kind = 8 ), external :: poly5
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x2(n2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINE_TRANSFORM_TEST03:'
  write ( *, '(a)' ) '  SINE_TRANSFORM_FUNCTION does a sine transform of data'
  write ( *, '(a)' ) '  defined by a function F(X) evaluated at N equally spaced'
  write ( *, '(a)' ) '  points in an interval [A,B].'
  write ( *, '(a)' ) '  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The interpolant will be 0 at the 0th and (N+1)-th points.'
  write ( *, '(a)' ) '  It equals the function at points 1 through N.'
  write ( *, '(a)' ) '  In between, it can approximate smooth functions,'
  write ( *, '(a)' ) '  and the approximation improves with N.'
!
!  N determines the number of data points, indexed by 1 to N.  
!  However, we essentially have N+2 data points, indexed 0 to N+1,
!  with the data value being 0 at the first and last auxilliary points.
!
  a = 1.0D+00
  b = 4.0D+00
!
!  Evenly spaced points between A and B, but omitting
!  A and B themselves.
!
  do i = 1, n
    x(i) = ( real ( n - i + 1, kind = 8 ) * a   &
           + real (     i,     kind = 8 ) * b ) &
           / real ( n     + 1, kind = 8 )
  end do
!
!  Determine the interpolant coefficients.
!
  call sine_transform_function ( n, a, b, poly5, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I      X(I)      F(X(I))        S(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, x(i), poly5 ( x(i) ), s(i)
  end do
!
!  Evaluate the interpolant.
!
  fa = poly5 ( a )
  fb = poly5 ( b )
!
!  Evenly spaced points between A and B, including A and B,
!  and twice the density of the previous set of points.
!
  do i = 1, n2
    x2(i) = ( real ( n2 - i,     kind = 8 ) * a   &
            + real (      i - 1, kind = 8 ) * b ) &
            / real ( n2     - 1, kind = 8 )
  end do

  call sine_transform_interpolant ( n, a, b, fa, fb, s, n2, x2, f2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I      X            F(X)        FHAT(X)'
  write ( *, '(a)' ) ' '

  do i = 1, n2
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, x2(i), poly5 ( x2(i) ), f2(i)
  end do

  return
end
subroutine sine_transform_test04 ( )

!*****************************************************************************80
!
!! SINE_TRANSFORM_TEST04 evaluates the sine transform interpolant.
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

  integer ( kind = 4 ), parameter :: n = 15
  integer ( kind = 4 ), parameter :: n2 = 1 + 5 * ( n + 1 )

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: cosine_sum
  real ( kind = 8 ) f2(n2)
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  integer ( kind = 4 ) i
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x2(n2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINE_TRANSFORM_TEST04:'
  write ( *, '(a)' ) '  SINE_TRANSFORM_FUNCTION does a sine transform of data'
  write ( *, '(a)' ) '  defined by a function F(X) evaluated at N equally spaced'
  write ( *, '(a)' ) '  points in an interval [A,B].'
  write ( *, '(a)' ) '  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The interpolant will be 0 at the 0th and (N+1)-th points.'
  write ( *, '(a)' ) '  It equals the function at points 1 through N.'
  write ( *, '(a)' ) '  In between, it can approximate smooth functions,'
  write ( *, '(a)' ) '  and the approximation improves with N.'
!
!  N determines the number of data points, indexed by 1 to N.  
!  However, we essentially have N+2 data points, indexed 0 to N+1,
!  with the data value being 0 at the first and last auxilliary points.
!
  a = 0.0D+00
  b = 7.0D+00
!
!  Evenly spaced points between A and B, but omitting
!  A and B themselves.
!
  do i = 1, n
    x(i) = ( real ( n - i + 1, kind = 8 ) * a   &
           + real (     i,     kind = 8 ) * b ) &
           / real ( n     + 1, kind = 8 )
  end do
!
!  Determine the interpolant coefficients.
!
  call sine_transform_function ( n, a, b, cosine_sum, s )
!
!  Evaluate the interpolant.
!
  fa = cosine_sum ( a )
  fb = cosine_sum ( b )
!
!  Evenly spaced points between A and B, including A and B,
!  and twice the density of the previous set of points.
!
  do i = 1, n2
    x2(i) = ( real ( n2 - i,     kind = 8 ) * a   &
            + real (      i - 1, kind = 8 ) * b ) &
            / real ( n2     - 1, kind = 8 )
  end do

  call sine_transform_interpolant ( n, a, b, fa, fb, s, n2, x2, f2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expect exact agreement every 5th sample.'
  write ( *, '(a)' ) ' '

  do i = 1, n2
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, x2(i), cosine_sum ( x2(i) ), f2(i)
  end do

  return
end
function cosine_sum ( x )

!*****************************************************************************80
!
!! COSINE_SUM evaluates a function which is a cosine sum.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) COSINE_SUM, the value.
!
  implicit none

  real ( kind = 8 ) cosine_sum
  real ( kind = 8 ) x

  cosine_sum =  cos (           x ) &
    + 5.0D+00 * cos ( 1.6D+00 * x ) &
    - 2.0D+00 * cos ( 2.0D+00 * x ) &
    + 5.0D+00 * cos ( 4.5D+00 * x ) &
    + 7.0D+00 * cos ( 9.0D+00 * x )

  return
end
function poly5 ( x )

!*****************************************************************************80
!
!! POLY5 evaluates a particular fifth-degree polynomial.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) POLY5, the value of the polynomial at X.
!
  implicit none

  real ( kind = 8 ) poly5
  real ( kind = 8 ) x

  poly5 = ( x - 0.1D+00 ) * &
          ( x - 0.2D+00 ) * &
          ( x - 0.4D+00 ) * &
          ( x - 2.1D+00 ) * &
          ( x - 3.0D+00 )

  return
end
