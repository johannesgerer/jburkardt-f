program main

!*****************************************************************************80
!
!! MAIN is the main program for KRONROD_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KRONROD_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the KRONROD library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KRONROD_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests the code for the odd case N = 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real    ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  real    ( kind = 8 ) s
  real    ( kind = 8 ) w1(n+1)
  real    ( kind = 8 ) w2(n+1)
  real    ( kind = 8 ) :: wg(n) = (/ &
    0.555555555555555555556D+00, &
    0.888888888888888888889D+00, &
    0.555555555555555555556D+00 /)
  real    ( kind = 8 ) :: wk(2*n+1) = (/ &
    0.104656226026467265194D+00, &
    0.268488089868333440729D+00, &
    0.401397414775962222905D+00, &
    0.450916538658474142345D+00, &
    0.401397414775962222905D+00, &
    0.268488089868333440729D+00, &
    0.104656226026467265194D+00 /)
  real    ( kind = 8 ) x(n+1)
  real    ( kind = 8 ) :: xg(n) = (/ &
   -0.77459666924148337704D+00, &
    0.0D+00, &
    0.77459666924148337704D+00 /)
  real    ( kind = 8 ) :: xk(2*n+1) = (/ &
   -0.96049126870802028342D+00, &
   -0.77459666924148337704D+00, &
   -0.43424374934680255800D+00, &
    0.0D+00, &
    0.43424374934680255800D+00, &
    0.77459666924148337704D+00, &
    0.96049126870802028342D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Request KRONROD to compute the Gauss rule'
  write ( *, '(a)' ) '  of order 3, and the Kronrod extension of'
  write ( *, '(a)' ) '  order 3+4=7.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare to exact data.'

  eps = 0.000001D+00

  call kronrod ( n, eps, x, w1, w2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i2)' ) '  KRONROD returns 3 vectors of length ', n + 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I      X               WK              WG'
  write ( *, '(a)' ) ' '
  do i = 1, n + 1
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w1(i), w2(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               Gauss Abscissas'
  write ( *, '(a)' ) '            Exact           Computed'
  write ( *, '(a)' ) ' '
  do i = 1, n
    if ( 2 * i <= n + 1 ) then
      i2 = 2 * i
      s = -1.0D+00
    else
      i2 = 2 * ( n + 1 ) - 2 * i
      s = +1.0D+00
    end if
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, xg(i), s * x(i2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               Gauss Weights'
  write ( *, '(a)' ) '            Exact           Computed'
  write ( *, '(a)' ) ' '
  do i = 1, n
    if ( 2 * i <= n + 1 ) then
      i2 = 2 * i
    else
      i2 = 2 * ( n + 1 ) - 2 * i
    end if
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, wg(i), w2(i2)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             Gauss Kronrod Abscissas'
  write ( *, '(a)' ) '            Exact           Computed'
  write ( *, '(a)' ) ' '
  do i = 1, 2 * n + 1
    if ( i <= n + 1 ) then
      i2 = i
      s = -1.0D+00
    else
      i2 = 2 * ( n + 1 ) - i
      s = +1.0D+00
    end if
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, xk(i), s * x(i2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             Gauss Kronrod Weights'
  write ( *, '(a)' ) '            Exact           Computed'
  write ( *, '(a)' ) ' '
  do i = 1, 2 * n + 1
    if ( i <= n + 1 ) then
      i2 = i
    else
      i2 = 2 * ( n + 1 ) - i
    end if
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, wk(i), w1(i2)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests the code for the even case N = 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real    ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  real    ( kind = 8 ) s
  real    ( kind = 8 ) w1(n+1)
  real    ( kind = 8 ) w2(n+1)
  real    ( kind = 8 ) x(n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Request KRONROD to compute the Gauss rule'
  write ( *, '(a)' ) '  of order 4, and the Kronrod extension of'
  write ( *, '(a)' ) '  order 4+5=9.'

  eps = 0.000001D+00

  call kronrod ( n, eps, x, w1, w2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i2)' ) '  KRONROD returns 3 vectors of length ', n + 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I      X               WK              WG'
  write ( *, '(a)' ) ' '
  do i = 1, n + 1
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w1(i), w2(i)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!  Purpose:
!
!    TEST03 uses the program to estimate an integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) eps
  real ( kind = 8 ) exact
  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  real ( kind = 8 ) i1
  real ( kind = 8 ) i2
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: w1(:)
  real ( kind = 8 ), allocatable :: w2(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Call Kronrod to estimate the integral of a function.'
  write ( *, '(a)' ) '  Keep trying until the error is small.'

  exact = 1.5643964440690497731D+00
!
!  EPS just tells KRONROD how carefully it must compute X, W1 and W2.
!  It is NOT a statement about the accuracy of your integral estimate!
!
  eps = 0.000001D+00
!
!  Start the process with a 1 point rule.
!
  n = 1

  do

    allocate ( x(n+1) )
    allocate ( w1(n+1) )
    allocate ( w2(n+1) )

    call kronrod ( n, eps, x, w1, w2 )
!
!  Compute the estimates.
!  There are two complications here:
!
!  1) Both rules use all the points.  However, the lower order rule uses
!     a zero weight for the points it doesn't need.
!
!  2) The points X are all positive, and are listed in descending order.
!     this means that 0 is always in the list, and always occurs as the
!     last member.  Therefore, the integral estimates should use the
!     function value at 0 once, and the function values at the other
!     X values "twice", that is, once at X and once at -X.
!
    i1 = w1(n+1) * f ( x(n+1) )
    i2 = w2(n+1) * f ( x(n+1) )

    do i = 1, n
      i1 = i1 + w1(i) * ( f ( - x(i) ) + f ( x(i) ) )
      i2 = i2 + w2(i) * ( f ( - x(i) ) + f ( x(i) ) )
    end do

    if ( abs ( i1 - i2 ) < 0.0001D+00 ) then 
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Error tolerance satisfied with N = ', n
      write ( *, '(a,g14.6)' ) '  Coarse integral estimate = ', i1
      write ( *, '(a,g14.6)' ) '  Fine   integral estimate = ', i2
      write ( *, '(a,g14.6)' ) '  Error estimate = ', abs ( i2 - i1 )
      write ( *, '(a,g14.6)' ) '  Actual error = ', abs ( exact - i2 )
      exit
    end if

    if ( 25 < n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Error tolerance failed even for n = ', n
      write ( *, '(a)' ) '  Canceling iteration, and accepting bad estimates!'
      write ( *, '(a,g14.6)' ) '  Coarse integral estimate = ', i1
      write ( *, '(a,g14.6)' ) '  Fine   integral estimate = ', i2
      write ( *, '(a,g14.6)' ) '  Error estimate = ', abs ( i2 - i1 )
      write ( *, '(a,g14.6)' ) '  Actual error = ', abs ( exact - i2 )
      exit
    end if

    deallocate ( x )
    deallocate ( w1 )
    deallocate ( w2 )

    n = 2 * n + 1

  end do

  return
end
function f ( x )

!*****************************************************************************80
!
!! F is a function whose integral from -1 to +1 is to be estimated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) F, the value.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) x

  f = 1.0D+00 / ( x * x + 1.005D+00 )

  return
end
