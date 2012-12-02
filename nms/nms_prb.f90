program main

!*****************************************************************************80
!
!! MAIN is the main program for NMS_PRB.
!
!  Discussion:
!
!    NMS_PRB calls the NMS tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NMS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NMS library.'

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
  call test024 ( )
  call test025 ( )
  call test026 ( )
  call test0265 ( )
  call test027 ( )
! call test028 ( )
  call test029 ( )

  call test030 ( )
  call test031 ( )
  call test032 ( )
  call test033 ( )
  call test034 ( )
  call test035 ( )
  call test0355 ( )
  call test036 ( )
  call test037 ( )
  call test038 ( )
  call test039 ( )

  call test040 ( )
  call test041 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NMS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests ALNGAM and GAMMA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alngam
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001:'
  write ( *, '(a)' ) '  ALNGAM evaluates the log of the Gamma function.'
  write ( *, '(a)' ) '  GAMMA_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       ALNGAM(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    if ( x <= 0.0D+00 ) then
      cycle
    end if

    fx = log ( fx )

    fx2 = alngam ( x )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests BESI0 and BESSEL_I0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) besi0
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002:'
  write ( *, '(a)' ) '  BESI0 evaluates the Bessel I0 function.'
  write ( *, '(a)' ) '  BESSEL_I0_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X	   Exact F	 BESI0(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call bessel_i0_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0D+00 ) then
      cycle
    end if

    fx2 = besi0 ( x )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests BESJ and BESSEL_J0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: alpha = 0.0D+00
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2(1)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003:'
  write ( *, '(a)' ) '  BESJ evaluates the Bessel function.'
  write ( *, '(a)' ) '  BESSEL_J0_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       BESJ(0)(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call bessel_j0_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0D+00 ) then
      cycle
    end if

    call besj ( x, alpha, 1, fx2(1), nz )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2(1)

  end do

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests BESJ and BESSEL_J1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: alpha = 1.0D+00
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2(1)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004:'
  write ( *, '(a)' ) '  BESJ evaluates the Bessel function.'
  write ( *, '(a)' ) '  BESSEL_J1_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       BESJ(1)(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call bessel_j1_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0D+00 ) then
      cycle
    end if

    call besj ( x, alpha, 1, fx2(1), nz )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2(1)

  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests BESJ and BESSEL_JN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2(1)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nu
  integer ( kind = 4 ) nz
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005:'
  write ( *, '(a)' ) '  BESJ evaluates the Bessel function.'
  write ( *, '(a)' ) '  BESSEL_JN_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NU     X       Exact F       BESJ(NU)(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call bessel_jn_values ( n, nu, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0D+00 ) then
      cycle
    end if

    alpha = real ( nu, kind = 8 )
    call besj ( x, alpha, 1, fx2(1), nz )

    write ( *, '(2x,i4,f8.4,2g14.6)' ) nu, x, fx, fx2(1)

  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests BP01 and BERNSTEIN_POLY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
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
  write ( *, '(a)' ) 'TEST006:'
  write ( *, '(a)' ) '  BP01 evaluates the Bernstein polynomials.'
  write ( *, '(a)' ) '  BERNSTEIN_POLY_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   K   X   Exact   B(N,K)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernstein_poly_values ( n_data, n, k, x, b )

    if ( n_data == 0 ) then
      exit
    end if

    call bp01 ( n, x, bvec )

    write ( *, '(2x,i4,i4,f7.4,2g14.6)' ) n, k, x, b, bvec(k)

  end do

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests CHKDER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: m = n
  integer ( kind = 4 ), parameter :: ldfjac = n

  real ( kind = 8 ) err(m)
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(m)
  real ( kind = 8 ) fvecp(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) mode
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xp(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  CHKDER compares a user supplied jacobian'
  write ( *, '(a)' ) '  and a finite difference approximation to it'
  write ( *, '(a)' ) '  and judges whether the jacobian is correct.'

  do ido = 1, 2

    if ( ido == 1 ) then

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) '  On the first test, use a correct jacobian.'

    else if ( ido == 2 ) then

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) '  Repeat the test, but use a "bad" jacobian'
       write ( *, '(a)' ) '  and see if the routine notices!'
       write ( *, '(a)' ) ' '

     end if
!
!  Set the point at which the test is to be made:
!
    x(1:n) = 0.5D+00

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Evaluation point X:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,g14.6)' ) x(i)
    end do

    mode = 1
    call chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )

    iflag = 1

    call f33 ( n, x, fvec, fjac, ldfjac, iflag )
    call f33 ( n, xp, fvecp, fjac, ldfjac, iflag )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Sampled function values F(X) and F(XP)'
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(2x,i3,2g14.6)' ) i, fvec(i), fvecp(i)
    end do

    iflag = 2
    call f33 ( n, x, fvec, fjac, ldfjac, iflag )
!
!  Here's where we put a mistake into the jacobian, on purpose.
!
    if ( ido == 2 ) then
      fjac(1,1) = 1.01D+00 * fjac(1,1)
      fjac(2,3) = - fjac(2,3)
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Computed jacobian'
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(2x,5g14.6)' ) fjac(i,1:n)
    end do

    mode = 2
    call chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  CHKDER error estimates:'
    write ( *, '(a)' ) '     > 0.5, gradient component is probably correct.'
    write ( *, '(a)' ) '     < 0.5, gradient component is probably incorrect.'
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(2x,i6,g14.6)' ) i, err(i)
    end do

  end do

  return
end
subroutine f33 ( n, x, fvec, fjac, ldfjac, iflag )

!*****************************************************************************80
!
!! F33 is a function/jacobian routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the variable values.
!
!    Output, real ( kind = 8 ) FVEC(N), the function values at X, if IFLAG = 1.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), the N by N jacobian at X,
!    if IFLAG = 2.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC, which must
!    be at least N.
!
!    Input, integer ( kind = 4 ) IFLAG:
!    1, please compute F(I) (X).
!    2, please compute FJAC(I,J) (X).
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) prod
  real ( kind = 8 ) x(n)
!
!  If IFLAG is 1, we are supposed to evaluate F(X).
!
  if ( iflag == 1 ) then

    do i = 1, n-1
      fvec(i) = x(i) - real ( n + 1, kind = 8 ) + sum ( x(1:n) )
    end do

    fvec(n) = product ( x(1:n) ) - 1.0D+00
!
!  If IFLAG is 2, we are supposed to evaluate FJAC(I,J) = d F(I)/d X(J)
!
  else if ( iflag == 2 ) then

    fjac(1:n-1,1:n) = 1.0D+00

    do i = 1, n-1
      fjac(i,i) = 2.0D+00
    end do

    prod = product ( x(1:n) )

    fjac(n,1:n) = prod / x(1:n)

  end if

  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests COST1I, COST1F and COST1B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096

  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) wsave(3*n+15)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  For real fast cosine transforms, 1D,'
  write ( *, '(a)' ) '  COSTI initializes the transforms,'
  write ( *, '(a)' ) '  COST does a forward or backward transform'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  seed = 1973

  call r8vec_uniform_01 ( n, seed, r )

  call r8vec_print_some ( n, r, 1, 10, '  First 10 data values:' )
!
!  Allocate and initialize the WSAVE array.
!
  call costi ( n, wsave )
!
!  Compute the FFT coefficients.
!
  call cost ( n, r, wsave )

  call r8vec_print_some ( n, r, 1, 10, '  First 10 FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cost ( n, r, wsave )
!
!  Normalize the data
!
  r(1:n) = r(1:n) / ( 2.0D+00 * real ( n - 1, kind = 8 ) )

  call r8vec_print_some ( n, r, 1, 10, '  First 10 retrieved data values:' )

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests DNOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nbins = 32

  real ( kind = 8 ), parameter :: a = -3.0D+00
  real ( kind = 8 ), parameter :: b = 3.0D+00
  real ( kind = 8 ) dnor
  real ( kind = 8 ) dstart
  integer ( kind = 4 ) h(nbins)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inbin
  integer ( kind = 4 ) iseed
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: nr = 10000
  real ( kind = 8 ) r
  real ( kind = 8 ) rseed
  real ( kind = 8 ) width

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  DNOR, normal random number generator.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of normal values to compute is ', nr
  write ( *, '(a,i6)' ) '  Number of bins is ', nbins

  iseed = 305
  rseed = dstart ( iseed )
  width = ( b - a ) / real  ( nbins - 2, kind = 8 )

  h(1:nbins) = 0

  do i = 1, nr
    r = dnor ( )
    j = inbin ( r, nbins, a, b, width )
    h(j) = h(j) + 1
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Histogram for DNOR: number in bin 1,...,32'
  write ( *, '(a)' ) '  (-infinity,-3],(-3,-2.8],...,(2.8,3],(3,infinity)'
  write ( *, '(a)' ) '  (values are slightly computer dependent)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,9i8)' ) h(1:nbins)

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests DNOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) dnor
  real ( kind = 8 ) dstart
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iseed
  real ( kind = 8 ) r
  real ( kind = 8 ) rseed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  DNOR generates random normal numbers.'
  write ( *, '(a)' ) ' '
!
!  Set the initial seed.
!
  iseed = 305
  rseed = dstart ( iseed )
!
!  DSTART returns a floating echo of ISEED.
!
  write ( *, '(a,i20)' ) '  ISEED = ', iseed
  write ( *, '(a,g14.6)' ) '  RSEED = ', rseed
  write ( *, '(a)' ) ' '

  do i = 1, 5
    r = dnor ( )
    write ( *, '(2x,g14.6)' ) r
  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests DDRIV1.
!
!  Discussion:
!
!    An example of the use of the ODE solver DDRIV1.
!
!    Here we solve the simple system
!
!      Y1' = Y2
!      Y2' = -Y1
!
!    with initial conditions
!
!      Y1(0) = 0
!      Y2(0) = 1
!
!    with exact solution
!
!      Y1(T) = SIN(T)
!      Y2(T) = COS(T)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  integer ( kind = 4 ), parameter :: lenw = n * n + 11 * n + 225

  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mstate
  integer ( kind = 4 ) nstep
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) tout
  real ( kind = 8 ) work(lenw)
  real ( kind = 8 ) y(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  DDRIV1 is a simple interface to the ODE solver.'
  write ( *, '(a)' ) ' '
!
!  Set the error tolerance.
!
  eps = 0.00001D+00
!
!  Set the initial time.
!
  t = 0.0D+00
  tout = t
!
!  Set the initial conditions
!
  y(1) = 0.0D+00
  y(2) = 1.0D+00
!
!  Set the number of steps we will take in the DO loop.
!
  nstep = 12
!
!  Tell DDRIV1 that this is the first call for this problem.
!
  mstate = 1
!
!  Print a header for the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Results'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   T        Y(1)      Y(2)'
  write ( *, '(a)' ) '            SIN(T)    COS(T)'
  write ( *, '(a)' ) '            Error     Error'
!
!  Call DDRIV1 NSTEP+1 times.
!
  do i = 0, nstep

    tout = real ( 2 * i, kind = 8 ) * pi / real ( nstep, kind = 8 )

    call ddriv1 ( n, t, y, tout, mstate, eps, work, lenw )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,3f11.5)' ) t, y(1), y(2)
    write ( *, '(13x,2f11.5)' ) sin(t), cos(t)
    write ( *, '(13x,2f11.5)' ) y(1)-sin(t), y(2)-cos(t)
!
!  Cancel the computation if we get any output code but 1 or 2.
!
    if ( mstate /= 1 .and. mstate /= 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16 - Fatal error!'
      write ( *, '(a,i4)' ) '  DDRIV1 returned MSTATE = ', mstate
      write ( *, '(a)' ) '  The computation is being cancelled.'
      return
    end if

  end do

  return
end
subroutine f ( n, t, y, ydot )

!*****************************************************************************80
!
!! F evaluates the right hand sides of the ODE's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) t
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ydot(n)

  ydot(1) = y(2)
  ydot(2) = -y(1)

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests DDRIV2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: nroot = 1

  integer ( kind = 4 ), parameter :: lw = n*n+10*n+2*nroot+204
  integer ( kind = 4 ), parameter :: liw = 23

  real ( kind = 8 ) eps
  real ( kind = 8 ) ewt
  external fsub
  real ( kind = 8 ), external :: gfun
  real ( kind = 8 ), parameter :: h = 10.0D+00
  integer ( kind = 4 ) iw(liw)
  real ( kind = 8 ), parameter :: mass = 0.125D+00
  integer ( kind = 4 ), parameter :: mint = 2
  integer ( kind = 4 ) ms
  real ( kind = 8 ) t
  real ( kind = 8 ) tout
  real ( kind = 8 ) w(lw)
  real ( kind = 8 ) y(n+1)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  DDRIV2 is an ODE solver.'
  write ( *, '(a)' ) ' '

  eps = 1.0D-05
!
!  Set initial point
!
  t = 0.0D+00
  tout = t
!
!  Set for pure relative error
!
  ewt = 0.0D+00
!
!  Set the initial conditions.
!
  y(1) = h
  y(2) = 0.0D+00
!
!  Set the parameter value.
!
  y(3) = mass
!
!  Set MS to 1, signaling the beginning of the run.
!
  ms = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DDRIV2 results'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   t,         y(1),      y(2),     ms '

  do

    call ddriv2 ( n, t, y, fsub, tout, ms, nroot, eps, ewt, mint, w, &
      lw, iw, liw, gfun )

    tout = tout + 0.1D+00

    if ( ms == 5 ) then
      write ( *, '(2x,3f11.5,i4,a,f11.5)' ) &
        t, y(1), y(2), ms, ' <-- y=0 at t= ', t
      exit
    else
      write ( *, '(2x,3f11.5,i4)' ) t,y(1),y(2),ms
!
!  Stop if any output code but 1 or 2.
!
      if ( 2 < ms ) then
        exit
      end if

    end if

  end do

  return
end
subroutine fsub ( n, t, y, ydot )

!*****************************************************************************80
!
!! FSUB returns the right hand side of an ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: g = 32.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) y(n+1)
  real ( kind = 8 ) ydot(n)

  ydot(1) = y(2)
  ydot(2) = - g - y(2) / y(3)

  return
end
function gfun ( n, t, y, iroot )

!*****************************************************************************80
!
!! GFUN ??
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) gfun
  integer ( kind = 4 ) iroot
  real ( kind = 8 ) t
  real ( kind = 8 ) y(n)

  gfun = y(1)

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests DNSQE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: lw = 19

  external f18
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iopt
  external j18
  integer ( kind = 4 ) nprint
  real ( kind = 8 ) tol
  real ( kind = 8 ) w(lw)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  DNSQE, nonlinear equation system solver.'
!
!  Set the parameters for DNSQE.
!
  tol = 1.0E-05

  x(1:2) = (/ 2.0D+00, 3.0D+00 /)

  call f18 ( n, x, fvec, iflag )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial solution estimate X0:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,2g14.6)' ) x(1:2)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Function value F(X0):'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,2g14.6)' ) fvec(1:2)

  iopt = 2
  nprint = 0
!
!  Solve the nonlinear equations.
!
  call dnsqe ( f18, j18, iopt, n, x, fvec, tol, nprint, info, w, lw )
!
!  Print the results.
!
  if ( info /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DNSQE INFO flag = ', info
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DNSQE solution estimate X:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,2g14.6)' ) x(1:2)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Function value F(X):'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,2g14.6)' ) fvec(1:2)

  return
end
subroutine f18 ( n, x, fvec, iflag )

!*****************************************************************************80
!
!! F18 evaluates a set of nonlinear equations whose zero is sought.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)

  fvec(1) = x(1) * x(2) - x(2)**3 - 1.0D+00
  fvec(2) = x(1)**2 * x(2) + x(2) - 5.0D+00

  return
end
subroutine j18 ( n, x, fvec, fjac, ldfjac, iflag )

!*****************************************************************************80
!
!! J18 is a dummy routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests DQRLS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mm = 5
  integer ( kind = 4 ), parameter :: nn = 3

  real ( kind = 8 ) a(mm,nn)
  real ( kind = 8 ) b(mm)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) itask
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jpvt(nn)
  integer ( kind = 4 ) kr
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) qraux(nn)
  real ( kind = 8 ) tol
  real ( kind = 8 ) work(nn)
  real ( kind = 8 ) x(nn)
!
!  Set up least-squares problem
!  quadratic model, equally-spaced points
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  DQRLS solves linear systems in the least squares sense.'
  write ( *, '(a)' ) ' '

  m = 5
  n = 3
  do i = 1, m
    a(i,1) = 1.0D+00
    do j = 2, n
      a(i,j) = a(i,j-1) * i
    end do
  end do

  b(1:5) = (/ 1.0D+00, 2.3D+00, 4.6D+00, 3.1D+00, 1.2D+00 /)

  tol = 1.0D-06

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coefficient matrix'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,3f12.6)' ) a(i,1:n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Right-hand side'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5f12.6)' ) b(1:m)
!
!  Solve least-squares problem
!
  itask = 1
  call dqrls ( a, mm, m, n, tol, kr, b, x, b, jpvt, qraux, work, itask, ind )
!
!  Print results
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Error code =', ind
  write ( *, '(a,i4)' ) '  Estimated matrix rank =', kr
  write ( *, '(a)' ) '  Parameters'
  write ( *, '(2x,3f12.6)' ) x(1:n)
  write ( *, '(a)' ) '  Residuals'
  write ( *, '(2x,5f12.6)' ) b(1:m)

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests DSVDC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ldx = 8
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: p = 3
  integer ( kind = 4 ), parameter :: ldu = n
  integer ( kind = 4 ), parameter :: ldv = p
  integer ( kind = 4 ), parameter :: job = 11

  real ( kind = 8 ) c(p)
  real ( kind = 8 ) e(p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ), dimension ( n ) :: pop = (/ 75.994575D+00, 91.972266D+00, &
    105.710620D+00, 122.775046D+00, 131.669275D+00, 150.697361D+00, &
    179.323175D+00, 203.235298D+00 /)
  real ( kind = 8 ) pop80
  real ( kind = 8 ) r
  real ( kind = 8 ) relerr
  real ( kind = 8 ) s(p)
  real ( kind = 8 ) sum2
  real ( kind = 8 ) tol
  real ( kind = 8 ) u(ldu,ldu)
  real ( kind = 8 ) v(ldv,p)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(ldx,p)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) year
!
!  C contains coefficients of the polynomial
!
!    c(1)*1+c(2)*t+c(3)*t*t
!
!  where t=year (1900 etc.).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  DSVDC computes the singular value decomposition.'

  do i = 1, 8
    y(i) = 1900.0D+00 + real ( ( i - 1 ) * 10, kind = 8 )
  end do

  x(1:8,1) = 1.0D+00
  x(1:8,2) = y(1:8)
  x(1:8,3) = y(1:8)**2

  call dsvdc ( x, ldx, n, p, s, e, u, ldu, v, ldv, w, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed singular values: '
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g12.4)' ) s(1:p)

  c(1:p) = 0.0D+00
!
!  RELERR reflects number of accurate digits in data
!  e.g. 6 digits ==> relerr=1.0E-06, and so on.
!  Making relerr larger increases residuals.
!
  relerr = 1.0D-06
  tol = relerr * s(1)
!
!  Multiply U' * pop, and solve for coefficients c(i)
!
  do j = 1, p

    if ( tol < s(j) ) then

      sum2 = dot_product ( pop(1:n), u(1:n,j) ) / s(j)

      c(1:p) = c(1:p) + sum2 * v(1:p,j)

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed polynomial coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g12.4)' ) c(1:p)
!
!  Evaluate the model using Horner's rule, and residuals at
!  year =1900,...,1980
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           Model       True'
  write ( *, '(a)' ) '  Year  Population  Population  Error'
  write ( *, '(a)') ' '

  r = 0.0D+00

  do i = 1, 9

    year = 1900.0D+00 + real ( i - 1, kind = 8 ) * 10.0D+00

    pop80 = 0.0D+00
    do j = p, 1, -1
      pop80 = year * pop80 + c(j)
    end do

    if ( i < 9 ) then
      r = r + ( pop(i) - pop80 )**2
      write ( *, '(2x,i4,3f10.2)' ) int(year), pop80, pop(i), pop(i) - pop80
    else
      write ( *, '(2x,i4,f10.3)' ) int(year), pop80
    end if

  end do

  r = sqrt ( r )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  RMS error is ', r

  return
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests DGEFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: lda = 10

  real ( kind = 8 ) a(lda,lda)
  real ( kind = 8 ) b(lda)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) itask
  integer ( kind = 4 ) iwork(lda)
  integer ( kind = 4 ) n
  real ( kind = 8 ) rcond
  real ( kind = 8 ) work(lda)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  DGEFS solves a system of linear equations.'
  write ( *, '(a)' ) ' '
!
!  Set the number of equations.
!
  n = 3

  itask = 1
!
!  Set the coefficient matrix A.
!
  a(1,1) = 10.0D+00
  a(2,1) = -3.0D+00
  a(3,1) =  5.0D+00
  a(1,2) = -7.0D+00
  a(2,2) =  2.0D+00
  a(3,2) = -1.0D+00
  a(1,3) =  0.0D+00
  a(2,3) =  6.0D+00
  a(3,3) =  5.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coefficient matrix A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,3f12.6)' ) a(i,1:n)
  end do
!
!  Set the right hand side vector B.
!
  b(1:3) = (/ 7.0D+00, 4.0D+00, 6.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Right-hand side B:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,3f12.6)') b(1:n)
!
!  Solve the linear system A*x=b.
!
  call dgefs ( a, lda, n, b, itask, ind, work, iwork, rcond )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DGEFS results:'
  write ( *, '(a)' ) ' '

  if ( ind == -10 ) then
    write ( *, '(a,i4)' ) '  Error code =',ind
  else if ( ind < 0 ) then
    write ( *, '(a,i4)' ) '  Error code =',ind
    return
  else
    write ( *, '(a,i4)' ) '  Estimated number of accurate digits =', ind
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,3f12.6)' ) b(1:n)

  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests ERROR_F and ERF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error_f
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017:'
  write ( *, '(a)' ) '  ERROR_F evaluates the Error function.'
  write ( *, '(a)' ) '  ERF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       ERROR_F(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call erf_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    fx2 = error_f ( x )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests ERROR_FC and ERF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error_fc
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018:'
  write ( *, '(a)' ) '  ERROR_FC evaluates the Complementary Error function.'
  write ( *, '(a)' ) '  ERF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       ERROR_FC(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call erf_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    fx = 1.0D+00 - fx
    fx2 = error_fc ( x )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests EZFFTI, EZFFTF and EZFFTB.
!
!  Discussion:
!
!    Find the autocorrelation to El Nino data using real FFT methods.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 168

  real ( kind = 8 ) a(2*n)
  real ( kind = 8 ) acovr(0:2*n-1)
  real ( kind = 8 ) azero
  real ( kind = 8 ) b(2*n)
  real ( kind = 8 ) el(0:2*n-1)
  real ( kind = 8 ) el_sum
  logical ex
  integer ( kind = 4 ) i
  real ( kind = 8 ) wsave(4*(2*n)+15)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  For the FFT of a real data sequence:'
  write ( *, '(a)' ) '  EZFFTI initializes,'
  write ( *, '(a)' ) '  EZFFTF does forward transforms,'
  write ( *, '(a)' ) '  EZFFTB does backward transforms.'
!
!  Check to see if the required data file exists.
!
  inquire ( file = 'elnino.dat', exist = ex )

  if ( .not. ex ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Cannot find the data file: elnino.dat '
    return
  end if

  open ( unit = 8, file = 'elnino.dat', status = 'old' )
!
!  Read the data, find the mean of the data.
!
  do i = 0, n-1
    read ( 8, * ) el(i)
  end do

  close ( unit = 8 )

  el_sum = sum ( el(0:n-1) )
!
!  Subtract the mean, and append N zeroes.
!
  el(0:n-1) = el(0:n-1) - el_sum / real ( n, kind = 8 )

  el(n:2*n-1) = 0.0D+00
!
!  fft approach (real).
!
!  Compute FFT of data of length 2n.
!  EZFFTF produces correctly scaled A's and B's so no extra scaling
!  is needed to get transform.
!
  call ezffti ( 2*n, wsave )

  call ezfftf ( 2*n, el, azero, a, b, wsave )
!
!  Compute array of square of each frequency component and place
!  in cosine array (a's) to be back transformed. set b's to 0.
!  There are n a's, and n b's.
!  Note that care must be taken to compute magnitude correctly,
!    0.5*(a(i)**2+b(i)**2) for i < n, twice that for i=n.
!
  azero = azero**2

  a(1:n-1) = ( a(1:n-1)**2 + b(1:n-1)**2 ) / 2.0D+00
  a(n) = a(n)**2 + b(n)**2

  b(1:n) = 0.0D+00
!
!  Compute the back transform, throwing away its second half.
!
  call ezfftb ( 2*n, acovr, azero, a, b, wsave )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Autocorrelation (real fft) output reduced.'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5e14.6)' ) acovr(0:19) / acovr(0)

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests EZFFTI, EZFFTF and EZFFTB.
!
!  Discussion:
!
!    Using the real discrete fourier transform, find the approximate
!    Fourier coefficients to Runge's function on [-1,1] with n=16 and
!    n=17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mcoef = 17

  real ( kind = 8 ) a(mcoef/2)
  real ( kind = 8 ) azero
  real ( kind = 8 ) b(mcoef/2)
  real ( kind = 8 ) c(mcoef/2)
  logical, parameter :: debug = .true.
  real ( kind = 8 ) del
  real ( kind = 8 ) dfta(mcoef/2)
  real ( kind = 8 ) dftb(mcoef/2)
  real ( kind = 8 ) error
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r(mcoef)
  real ( kind = 8 ), external :: runge
  real ( kind = 8 ) s(mcoef/2)
  real ( kind = 8 ) tn
  real ( kind = 8 ) wsave(3*mcoef+15)
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) xj

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020'
  write ( *, '(a)' ) '  The "EZ" FFT package:'
  write ( *, '(a)' ) '  EZFFTI initializes,'
  write ( *, '(a)' ) '  EZFFTF does forward transforms,'
  write ( *, '(a)' ) '  EZFFTB does backward transforms.'

  x0 = -1.0D+00

  do n = mcoef-1, mcoef

    call ezffti ( n, wsave )
!
!  Function assumed to be periodic on [-1,1], of length 2.
!
    del = 2.0D+00 / real ( n, kind = 8 )
    f = 2.0D+00 * pi / ( real ( n, kind = 8 ) * del )
!
!  The first sample point is at -1, the last at 1-del.
!
    do j = 1, n

      xj = (-1.0D+00) + (j-1) * del
      r(j) = runge ( xj )
!
!  Compute sines and cosines to adjust output of EZFFTF to give
!  approximate Fourier coefficients.
!
      if ( j <= n/2 ) then
        c(j) = cos ( j * f * x0 )
        s(j) = sin ( j * f * x0 )
      end if

    end do

    call ezfftf ( n, r, azero, a, b, wsave )
!
!  As a convenience this loop can go to N/2.  If N is even, the last B is zero.
!
    do j = 1, n/2
      dfta(j) = a(j) * c(j) - b(j) * s(j)
      dftb(j) = a(j) * s(j) + b(j) * c(j)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  EZFFTF results'
    write ( *, '(a,i4)' ) '  N = ' , n
    write ( *, '(a,g14.6)' )  '  AZERO = ', azero
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  j             dfta(j)              dftb(j) '
    write ( *, '(a)' ) ' '

    do j = 1, n/2
      write ( *, '(2x,i6,2g14.6)' ) j, dfta(j), dftb(j)
    end do
!
!  Evaluate interpolant at points on [-1,1]
!
    if ( debug ) then

      m = 21

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) ' Trigonometric polynomial order n= ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X    Interpolant    Runge     Error'
      write ( *, '(a)' ) ' '

      do k = 1, m

        x = - 1.0D+00 + 2.0D+00 * ( k - 1.0D+00 ) / real ( m - 1, kind = 8 )

        tn = azero
        do j = 1, n/2
          tn = tn + dfta(j) * cos(j*f*x) + dftb(j) * sin(j*f*x)
        end do

        error = tn - runge ( x )
        write ( *, '(2x,4g14.6)' ) x, tn, runge(x), error

      end do

    end if

  end do

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests FMIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fmin
  real ( kind = 8 ) fp05
  real ( kind = 8 ), external :: fx05
  real ( kind = 8 ) tol
  real ( kind = 8 ) xstar

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  FMIN, function minimizer.'
  write ( *, '(a)' ) '  Find a minimizer of F(X) = X^3 - 2 * X - 5.'

  a = 0.1D+00
  b = 0.9D+00
  tol = 1.0E-07

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The initial interval is [A,B]:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A =      ', a
  write ( *, '(a,g14.6)' ) '  F(A) =   ', fx05 ( a)
  write ( *, '(a,g14.6)' ) '  F''(A) = ', fp05 ( a )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  B =      ', b
  write ( *, '(a,g14.6)' ) '  F(B) =   ', fx05 ( b)
  write ( *, '(a,g14.6)' ) '  F''(B) = ', fp05 ( b )

  xstar = fmin ( a, b, fx05, tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The final interval [A,B] and minimizer X*:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A =      ', a
  write ( *, '(a,g14.6)' ) '  F(A) =   ', fx05 ( a)
  write ( *, '(a,g14.6)' ) '  F''(A) = ', fp05 ( a )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  B =      ', b
  write ( *, '(a,g14.6)' ) '  F(B) =   ', fx05 ( b)
  write ( *, '(a,g14.6)' ) '  F''(B) = ', fp05 ( b )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X* =      ', xstar
  write ( *, '(a,g14.6)' ) '  F(X*) =   ', fx05 ( xstar)
  write ( *, '(a,g14.6)' ) '  F''(X*) = ', fp05 ( xstar )

  return
end
subroutine test022 ( )

!*****************************************************************************80
!
!! TEST022 tests FMIN_RC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) fp05
  real ( kind = 8 ) fx05
  integer ( kind = 4 ) status
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022'
  write ( *, '(a)' ) '  FMIN_RC, function minimizer with reverse communication.'
  write ( *, '(a)' ) '  Find a minimizer of F(X) = X^3 - 2 * X - 5.'

  a = 0.1D+00
  b = 0.9D+00
  status = 0
  arg = 0.0D+00
  value = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The initial interval is [A,B]:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A =      ', a
  write ( *, '(a,g14.6)' ) '  F(A) =   ', fx05 ( a)
  write ( *, '(a,g14.6)' ) '  F''(A) = ', fp05 ( a )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  B =      ', b
  write ( *, '(a,g14.6)' ) '  F(B) =   ', fx05 ( b)
  write ( *, '(a,g14.6)' ) '  F''(B) = ', fp05 ( b )

  do

    call fmin_rc ( a, b, arg, status, value )

    if ( status == 0 ) then
      exit
    end if

    value = fx05 ( arg )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The final interval [A,B] and minimizer X*:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A =      ', a
  write ( *, '(a,g14.6)' ) '  F(A) =   ', fx05 ( a)
  write ( *, '(a,g14.6)' ) '  F''(A) = ', fp05 ( a )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  B =      ', b
  write ( *, '(a,g14.6)' ) '  F(B) =   ', fx05 ( b)
  write ( *, '(a,g14.6)' ) '  F''(B) = ', fp05 ( b )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X* =      ', arg
  write ( *, '(a,g14.6)' ) '  F(X*) =   ', fx05 ( arg )
  write ( *, '(a,g14.6)' ) '  F''(X*) = ', fp05 ( arg )

  return
end
function fx05 ( x )

!*****************************************************************************80
!
!! FX05 is a function to be minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx05
  real ( kind = 8 ) x

  fx05 = x * ( x * x - 2.0D+00 ) - 5.0D+00

  return
end
function fp05 ( x )

!*****************************************************************************80
!
!! FP05 is the derivative of a function to be minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fp05
  real ( kind = 8 ) x

  fp05 = 3.0D+00 * x * x - 2.0D+00

  return
end
subroutine test023 ( )

!*****************************************************************************80
!
!! TEST023 tests FZERO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ae
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), external :: fx06
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) re

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023'
  write ( *, '(a)' ) '  FZERO, single nonlinear equation solver.'
  write ( *, '(a)' ) '  F(X) = X^3 - 2 * X - 5'
  write ( *, '(a)' ) ' '

  b = 2.0D+00
  c = 3.0D+00
  ae = 1.0E-06
  re = 1.0E-06

  write ( *, '(a)' ) '  Initial interval: '
  write ( *, '(2x,2g14.6)' ) b, c
  write ( *, '(a,g14.6)' ) '  Absolute error tolerance=', ae
  write ( *, '(a,g14.6)' ) '  Relative error tolerance=', re

  call fzero ( fx06, b, c, c, re, ae, iflag )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FZERO results'
  write ( *, '(a)' ) ' '

  if ( iflag /= 1 ) then
    write ( *, '(a,i4)' ) '  FZERO returned error code =', iflag
  end if

  write ( *, '(a,g14.6)' ) '  Estimate of zero = ', b
  write ( *, '(a,g14.6)' ) '  Function value=    ', fx06(b)

  return
end
function fx06 ( x )

!*****************************************************************************80
!
!! FX06 is a function whose zero is desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx06
  real ( kind = 8 ) x

  fx06 = x * ( x * x - 2.0D+00 ) - 5.0D+00

  return
end
subroutine test024 ( )

!*****************************************************************************80
!
!! TEST024 tests GAMMA and GAMMA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024:'
  write ( *, '(a)' ) '  GAMMA evaluates the Gamma function.'
  write ( *, '(a)' ) '  GAMMA_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       GAMMA(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call gamma_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    fx2 = gamma ( x )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests PCHEZ and PCHEV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21
  integer ( kind = 4 ), parameter :: nwk = 42
  integer ( kind = 4 ), parameter :: ne = 101

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) diff
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fd(ne)
  real ( kind = 8 ) fe(ne)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  real ( kind = 8 ), external :: runge
  real ( kind = 8 ), external :: rungep
  logical spline
  real ( kind = 8 ) wk(nwk)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xe(ne)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  PCHEZ carries out piecewise cubic '
  write ( *, '(a)' ) '    spline or Hermite interpolation.'
  write ( *, '(a)' ) '  PCHEV evaluates the interpolant.'
  write ( *, '(a)' ) ' '
!
!  Compute Runge's function at N points in [-1,1].
!
  do i = 1, n
    x(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / 10.0D+00
    f(i) = runge ( x(i) )
  end do
!
!  Setting SPLINE = FALSE means we are requesting a piecewise
!  cubic Hermite interpolant.
!
  spline = .false.
!
!  PCHEZ takes the data in X and F, and constructs a table in D
!  that defines the interpolant.
!
  call pchez ( n, x, f, d, spline, wk, nwk, ierr )

  if ( ierr < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a,i4)' ) '  PCHEZ returned error code IERR = ', ierr
    return
  end if
!
!  Evaluate the interpolant and derivative at NE points from -1 to 0.
!
  do i = 1, ne
    xe(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / real ( ne - 1, kind = 8 )
  end do

  call pchev ( n, x, f, d, ne, xe, fe, fd, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a,i4)' ) '  PCHEV returned error code IERR = ', ierr
    return
  end if
!
!  Print the table of X, F(exact) and F(interpolated)
!
  do i = 1, ne
    diff = fe(i) - runge ( xe(i) )
    write ( *, '(2x,f8.4,2x,f10.6,2x,f10.6,2x,g14.6)' ) &
      xe(i), runge ( xe(i) ), fe(i), diff
  end do

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests PCHEZ and PCHQA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21
  integer ( kind = 4 ), parameter :: nwk = 42
  integer ( kind = 4 ), parameter :: ne = 101

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) pchqa
  real ( kind = 8 ) q
  real ( kind = 8 ), external :: runge
  real ( kind = 8 ), external :: rungep
  logical spline
  real ( kind = 8 ) wk(nwk)
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  PCHEZ carries out piecewise cubic '
  write ( *, '(a)' ) '    spline or Hermite interpolation.'
  write ( *, '(a)' ) '  PCHQA integrates the interpolant.'
  write ( *, '(a)' ) ' '
!
!  Compute Runge's function at N points in [-1,1].
!
  do i = 1, n
    x(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / 10.0D+00
    f(i) = runge ( x(i) )
  end do
!
!  Setting SPLINE = FALSE means we are requesting a piecewise
!  cubic Hermite interpolant.
!
  spline = .false.
!
!  PCHEZ takes the data in X and F, and constructs a table in D
!  that defines the interpolant.
!
  call pchez ( n, x, f, d, spline, wk, nwk, ierr )

  if ( ierr < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a,i4)' ) '  PCHEZ returned error code IERR = ', ierr
    return
  end if
!
!  Compute the integral over the interval [0,1].
!
  a = 0.0D+00
  b = 1.0D+00
  q = pchqa ( n, x, f, d, a, b, ierr )

  if ( ierr < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a,i4)' ) '  PCHQA returned error code IERR = ', ierr
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PCHQA estimates the integral from A to B.'
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a,g14.6)' ) '  Integral estimate = ', q
  write ( *, '(a,i4)' ) '  Return code IERR = ', ierr

  return
end
subroutine test0265 ( )

!*****************************************************************************80
!
!! TEST0265 tests PCHIM and PCHFE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21
  integer ( kind = 4 ), parameter :: ne = 101

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) diff
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fd(ne)
  real ( kind = 8 ) fe(ne)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), parameter :: incfd = 1
  real ( kind = 8 ), external :: runge
  real ( kind = 8 ), external :: rungep
  logical skip
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xe(ne)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0265'
  write ( *, '(a)' ) '  PCHIM carries out piecewise cubic '
  write ( *, '(a)' ) '    Hermite interpolation.'
  write ( *, '(a)' ) '  PCHFE evaluates the interpolant.'
  write ( *, '(a)' ) ' '
!
!  Compute Runge's function at N points in [-1,1].
!
  do i = 1, n
    x(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / 10.0D+00
    f(i) = runge ( x(i) )
  end do
!
!  PCHIM takes the data in X and F, and constructs a table in D
!  that defines the interpolant.
!
  call pchim ( n, x, f, d, incfd, ierr )

  if ( ierr < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a,i4)' ) '  PCHIM returned error code IERR = ', ierr
    return
  end if
!
!  Evaluate the interpolant and derivative at NE points from -1 to 0.
!
  do i = 1, ne
    xe(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / real ( ne - 1, kind = 8 )
  end do

  skip = .false.

  call pchfe ( n, x, f, d, incfd, skip, ne, xe, fe, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a,i4)' ) '  PCHEV returned error code IERR = ', ierr
    return
  end if
!
!  Print the table of X, F(exact) and F(interpolated)
!
  do i = 1, ne
    diff = fe(i) - runge ( xe(i) )
    write ( *, '(2x,f8.4,2x,f10.6,2x,f10.6,2x,g14.6)' ) &
      xe(i), runge ( xe(i) ), fe(i), diff
  end do

  return
end
subroutine test027 ( )

!*****************************************************************************80
!
!! TEST027 tests Q1DA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) e
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: f10
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) kf
  real ( kind = 8 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  Q1DA, quadrature routine.'
  write ( *, '(a)' ) ' '

  a = 0.0D+00
  b = 1.0D+00

  eps = 0.001D+00

  call q1da ( f10, a, b, eps, r, e, kf, iflag )

  if ( iflag < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a,i4)' ) '  Q1DA returned error code IFLAG = ', iflag
    return
  end if

  write ( *, '(a)' ) '  Q1DA results: a, b, eps, r, e, kf, iflag'
  write ( *, '(2x,3f7.4,2e16.8,2i4)' ) a,b,eps,r,e,kf,iflag

  return
end
function f10 ( x )

!*****************************************************************************80
!
!! F10 is a function to be integrated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f10
  real ( kind = 8 ) x

  f10 = sin ( 2.0D+00 * x ) - sqrt ( x )

  return
end
subroutine test028 ( )

!*****************************************************************************80
!
!! TEST028 tests Q1DA.
!
!  Discussion:
!
!    Compute a double integral by two one dimensional subroutines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a_x
  real ( kind = 8 ) b_x
  real ( kind = 8 ) eps_x
  real ( kind = 8 ) err
  real ( kind = 8 ), external :: g13
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) kf
  real ( kind = 8 ) result
  real ( kind = 8 ) x_fixed

  common /comm13/ x_fixed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028'
  write ( *, '(a)' ) '  Demonstration of two-dimensional quadrature'
  write ( *, '(a)' ) '  using one-dimensional techniques.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Q1DA is called to integrate'
  write ( *, '(a)' ) '  G(X) = Integral ( F(X,Y) dY ).'

  a_x = 0.0D+00
  b_x = 1.0D+00
  eps_x = 1.0E-04

  call q1da ( g13, a_x, b_x, eps_x, result, err, kf, iflag )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  EPS_X =  ', eps_x
  write ( *, '(a,g14.6)' ) '  Result = ', result


  return
end
function g13 ( x )

!*****************************************************************************80
!
!! G13 returns G13(X) = Integral ( F(X,Y) dY )
!
!  Discussion:
!
!    The routine is used as part of an effort to integrate F(X,Y)
!    over a two dimensional region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a_y
  real ( kind = 8 ) b_y
  real ( kind = 8 ) eps
  real ( kind = 8 ) err
  real ( kind = 8 ), external :: f13
  real ( kind = 8 ) g13
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) kf
  real ( kind = 8 ) result
  real ( kind = 8 ) x
  real ( kind = 8 ) x_fixed

  common /comm13/ x_fixed

  x_fixed = x

  a_y = 0.0D+00
  b_y = 2.0D+00
  eps = 1.0E-04

  call q1da ( f13, a_y, b_y, eps, result, err, kf, iflag )

  g13 = result

  return
end
function f13 ( y )

!*****************************************************************************80
!
!! F13 is a function of two variables to be integrated.
!
!  Discussion:
!
!    The function is being integrated over a Y range, with a fixed
!    value of X.  The fixed value of X is passed through a common block.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f13
  real ( kind = 8 ) x_fixed
  real ( kind = 8 ) y

  common /comm13/ x_fixed

  f13 = exp ( - x_fixed**2 * y**2 )

  return
end
subroutine test029 ( )

!*****************************************************************************80
!
!! TEST029 tests Q1DAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 50

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) e
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: f105
  real ( kind = 8 ) fmax
  real ( kind = 8 ) fmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) kf
  real ( kind = 8 ) r
  logical rst
  real ( kind = 8 ) w(nmax,6)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029'
  write ( *, '(a)' ) '  Q1DAX estimates the integral of a function over a'
  write ( *, '(a)' ) '  a finite interval, allowing more flexibility than Q1DA.'

  a = 0.0D+00
  b = 1.0D+00
!
!  Set up an initial partition of 2 intervals, with an internal
!  partition point at 0.3.
!
  w(1,1) = a
  w(2,1) = 0.3D+00
  w(3,1) = b

  int_num = 2

  do i = 1, 2

    if ( i == 1 ) then
      eps = 0.001D+00
      rst = .false.
    else
      eps = 0.0001D+00
      rst = .true.
    end if

    call q1dax ( f105, a, b, eps, r, e, int_num, rst, w, nmax, fmin, fmax, &
      kf, iflag )

    if ( 3 <= iflag ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Q1DAX error flag = ', iflag
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Error tolerance      ', eps
    write ( *, '(a,g14.6)' ) '  Integral estimate    ', r
    write ( *, '(a,g14.6)' ) '  Error estimate       ', e

  end do

  return
end
function f105 ( x )

!*****************************************************************************80
!
!! F105 is a function to be integrated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f105
  real ( kind = 8 ) x

  if ( x < 0.3D+00 ) then
    f105 = x**( 0.2D+00 ) * log ( x )
  else
    f105 = sin ( x )
  end if

  return
end
subroutine test030 ( )

!*****************************************************************************80
!
!! TEST030 tests QAGI.
!
!  Discussion:
!
!    Compute integral of exp(-x)*cos(x*x)**2 on [0,infinity)
!    Correct result is 0.70260...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 100
  integer ( kind = 4 ), parameter :: lenw = limit * 4

  real ( kind = 8 ) abserr
  real ( kind = 8 ) bound
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ), external :: f11
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inf
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) work(lenw)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST030'
  write ( *, '(a)' ) '  QAGI estimates an integral on a semi-infinite interval.'

  bound = 0.0D+00
  inf = 1
  epsabs = 0.0D+00
  epsrel = 1.0D-05

  call qagi ( f11, bound, inf, epsabs, epsrel, result, abserr, neval, &
    ier, limit, lenw, last, iwork, work )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimated integral =   ', result
  write ( *, '(a,g14.6)' ) '  (Correct value)    =   ', 0.70260D+00
  write ( *, '(a,g14.6)' ) '  Estimated error =      ', abserr
  write ( *, '(a,i4)' ) '  Function evaluations = ', neval
  write ( *, '(a,i4)' ) '  Return code IER =      ', ier

  return
end
function f11 ( x )

!*****************************************************************************80
!
!! F11 is a function to be integrated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f11
  real ( kind = 8 ) x

  f11 = exp ( -x ) * cos ( x**2 )**2

  return
end
subroutine test031 ( )

!*****************************************************************************80
!
!! TEST031 tests QK15.
!
!  Discussion:
!
!    Compute erf(1), i.e. integral of 2/sqrt(pi) * exp(-x*x) from 0
!    to 1.0D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f12
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST031'
  write ( *, '(a)' ) '  QK15 estimates an integral using '
  write ( *, '(a)' ) '  Gauss-Kronrod integration.'

  a = 0.0D+00
  b = 1.0D+00

  call qk15 ( f12, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  QK15 estimate of ERF(1) '
  write ( *, '(a)' ) '  2 / sqrt ( PI ) * result,      abserr'
  write ( *, '(2x,2g14.6)' ) 2.0D+00 / sqrt ( pi ) * result, abserr

  return
end
function f12 ( x )

!*****************************************************************************80
!
!! F12 is a function to be integrated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f12
  real ( kind = 8 ) x

  f12 = exp ( - x**2 )

  return
end
subroutine test032 ( )

!*****************************************************************************80
!
!! TEST032 tests SINT and SINTI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096

  real ( kind = 8 ), parameter :: ahi = 5.0D+00
  real ( kind = 8 ), parameter :: alo = 0.0D+00
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) wsave((5*n+30)/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032'
  write ( *, '(a)' ) '  For sine analysis of real data,'
  write ( *, '(a)' ) '  SINTI initializes the FFT routines.'
  write ( *, '(a)' ) '  SINT does a forward or backward FFT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  seed = 123456789
  call r8vec_uniform_01 ( n, seed, c )

  c(1:n) = alo + c(1:n) * ( ahi - alo )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first 10 data values:'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,g14.6)' ) c(i)
  end do
!
!  Initialize the WSAVE array.
!
  call sinti ( n, wsave )
!
!  Compute the coefficients.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute the sine coefficients from data.'

  call sint ( n, c, wsave )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first 10 sine coefficients:'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,g14.6)' ) c(i)
  end do
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Retrieve data from coeficients.'

  call sint ( n, c, wsave )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first 10 data values:'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,g14.6)' ) c(i) / real ( 2 * ( n + 1 ), kind = 8 )
  end do

  return
end
subroutine test033 ( )

!*****************************************************************************80
!
!! TEST033 tests SINT1I, SINT1F and SINT1B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096

  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) wsave((5*n+30)/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST033'
  write ( *, '(a)' ) '  For real fast sine transforms, 1D,'
  write ( *, '(a)' ) '  SINTI initializes the transforms,'
  write ( *, '(a)' ) '  SINT does a forward or backward transform'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  seed = 1973

  call r8vec_uniform_01 ( n, seed, r )

  call r8vec_print_some ( n, r, 1, 10, '  First 10 data values:' )
!
!  Allocate and initialize the WSAVE array.
!
  call sinti ( n, wsave )
!
!  Compute the FFT coefficients.
!
  call sint ( n, r, wsave )

  call r8vec_print_some ( n, r, 1, 10, '  Fist 10 FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call sint ( n, r, wsave )
!
!  Normalize the data
!
  r(1:n) = r(1:n) / ( 2.0D+00 * real ( n + 1, kind = 8 ) )

  call r8vec_print_some ( n, r, 1, 10, '  First 10 retrieved data values:' )

  return
end
subroutine test034 ( )

!*****************************************************************************80
!
!! TEST034 tests UNCMIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: lwork = n*(n+10)

  real ( kind = 8 ) f
  external f8
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) work(lwork)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST034'
  write ( *, '(a)' ) '  UNCMIN, unconstrained minimization code.'
  write ( *, '(a)' ) ' '
!
!  Specify an initial estimate of the solution.
!
  x0(1) = 1.0D+00
  x0(2) = 1.0D+00
!
!  Minimize function
!
  call uncmin ( n, x0, f8, x, f, ierror, work, lwork )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  UNCMIN return code =', ierror
  write ( *, '(a,g14.6)' ) '  f(x*) =', f
  write ( *, '(a)' ) '  x* ='
  write ( *, '(2x,4g14.6)' ) x(1:n)

  return
end
subroutine f8 ( n, x, f )

!*****************************************************************************80
!
!! F8 is a function to be minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter, dimension ( m ) :: b = &
    (/ 20.0D+00, 9.0D+00, 3.0D+00, 1.0D+00 /)
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter, dimension ( m ) :: t = &
    (/ 0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ) x(n)

  f = 0.0D+00
  do j = 1, m
    f = f + ( b(j) - x(1) * exp ( x(2) * t(j) ) )**2
  end do

  return
end
subroutine test035 ( )

!*****************************************************************************80
!
!! TEST035 tests UNCMIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: lwork = n*(n+10)

  real ( kind = 8 ) f
  external f21
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) work(lwork)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  UNCMIN carries out unconstrained minimization'
  write ( *, '(a)' ) '  of a scalar function of several variables.'

  do i = 1, n
    x0(i) = real ( i, kind = 8 ) / real ( n + 1, kind = 8 )
  end do
!
!  Minimize the function.
!
  call uncmin ( n, x0, f21, x, f, ierror, work, lwork )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Return code =', ierror
  write ( *, '(a,g14.6)' ) '  f(x*) =', f
  write ( *, '(a)' ) '  x* ='
  write ( *, '(2x,5f12.6)' ) x(1:n)

  return
end
subroutine f21 ( n, x, f )

!*****************************************************************************80
!
!! F21 is a function to be minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the X vector.
!
!    Input, real ( kind = 8 ) X(N), the value of the X vector.
!
!    Output, real ( kind = 8 ) F, the value of the function F(X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) x(n)

  t1 = 0.0D+00
  do i = 2, n
    t1 = t1 + ( x(i) - x(i-1)**2 )**2
  end do

  t2 = sum ( ( 1.0D+00 - x(1:n-1) )**2 )

  f = 1.0D+00 + 100.0D+00 * t1 + t2

  return
end
subroutine test0355 ( )

!*****************************************************************************80
!
!! TEST0355 tests UNCMIN.
!
!  Discussion:
!
!    In this problem, we start with a distance chart, and try to assign
!    XY positions to cities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: city_num = 4
  integer ( kind = 4 ), parameter :: n = 2 * city_num
  integer ( kind = 4 ), parameter :: lwork = n*(n+10)

  real    ( kind = 8 ) d(city_num,city_num)
  real    ( kind = 8 ) f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  external             map
  real    ( kind = 8 ) work(lwork)
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) x0(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0355'
  write ( *, '(a)' ) '  UNCMIN carries out unconstrained minimization'
  write ( *, '(a)' ) '  of a scalar function of several variables.'

  do i = 1, n
    x0(i) = real ( i, kind = 8 ) / real ( n + 1, kind = 8 )
  end do

! x0 = (/ 0.0, 0.0, 3.0, 0.0, 4.0, 3.0, 0.0, 4.0 /)
!
!  Minimize the function.
!
  call uncmin ( n, x0, map, x, f, ierror, work, lwork )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Return code =', ierror
  write ( *, '(a,g14.6)' ) '  f(x*) =', f
  write ( *, '(a)' ) '  x* ='

  do i = 1, n / 2
    write ( *, '(2x,i8,2x,f12.6,2x,f12.6)' ) i, x(2*i-1), x(2*i)
  end do

  do i = 1, n/2
    do j = 1, n/2
      d(i,j) = sqrt ( ( x(i*2-1) - x(j*2-1) )**2 &
                    + ( x(i*2  ) - x(j*2  ) )**2 )
    end do
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) d(i,1:n/2)
  end do

  return
end
subroutine map ( n, x, f )

!*****************************************************************************80
!
!! MAP is a function to be minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the X vector.
!
!    Input, real ( kind = 8 ) X(N), the value of the X vector.
!
!    Output, real ( kind = 8 ) F, the value of the function F(X).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) city_num
  integer ( kind = 4 ) city1
  integer ( kind = 4 ) city2
  real    ( kind = 8 ) d
  real    ( kind = 8 ) :: distance(4,4) = reshape ( (/ &
   0.0D+00,   3.0D+00,       5.0D+00,        4.0D+00, &
   3.0D+00,   0.0D+00,       3.162277D+00,   5.0D+00, &
   5.0D+00,   3.162277D+00,  0.0D+00,        4.1231055D+00, &
   4.0D+00,   5.0D+00,       4.1231055D+00,  0.0D+00 /), (/ 4, 4 /) )
  real    ( kind = 8 ) f
  integer ( kind = 4 ) i
  real    ( kind = 8 ) r
  real    ( kind = 8 ) x(n)

  f = 0.0D+00
!
!  Force the first city to be at (0,0) and the second city to have
!  Y coordinate 0.
!
  f = f + x(1)**2 + x(2)**2 + x(4)**2

  city_num = n / 2

  do city1 = 1, city_num
    do city2 = city1 + 1, city_num

      d = sqrt ( ( x(city1*2-1) - x(city2*2-1) )**2 &
               + ( x(city1*2  ) - x(city2*2  ) )**2 )

      r = distance(city1,city2) - d

      f = f + r * r

    end do
  end do

  return
end
subroutine test036 ( )

!*****************************************************************************80
!
!! TEST036 tests UNI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iseed
  real ( kind = 8 ) u
  real ( kind = 8 ) uni
  real ( kind = 8 ) ustart
  real ( kind = 8 ) useed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036'
  write ( *, '(a)' ) '  UNI is a uniform random number generator.'
!
!  Set the initial seed.
!
  iseed = 305
  useed = ustart ( iseed )
!
!  USTART returns floating echo of iseed.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i20)' ) '  The seed value ISEED is ', iseed
  write ( *, '(a,g14.6)' ) '  The starting value is ', useed

  do i = 1, 1000
    u = uni ( )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The 1000-th random number generated is ', u

  return
end
subroutine test037 ( )

!*****************************************************************************80
!
!! TEST037 computes an autocorrelation using a direct method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 168

  real ( kind = 8 ) acov(0:n-1)
  real ( kind = 8 ) el(0:2*n-1)
  real ( kind = 8 ) el_sum
  logical ex
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037'
  write ( *, '(a)' ) '  Compute the autocorrelation of El Nino data'
  write ( *, '(a)' ) '  using a direct method.'
!
!  Check to see if the required data file exists.
!
  inquire ( file = 'elnino.dat', exist = ex )

  if ( .not. ex ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST004 - Fatal error!'
    write ( *, '(a)' ) '  Cannot find the data file: elnino.dat '
    return
  end if

  open ( unit = 8, file = 'elnino.dat', status = 'old' )
!
!  Read the data, find the mean of the data.
!
  do i = 0, n-1
    read ( 8, * ) el(i)
  end do

  close ( unit = 8 )

  el_sum = sum ( el(0:n-1) )
!
!  Subtract the mean, and append N zeroes.
!
  el(0:n-1) = el(0:n-1) - el_sum / real ( n, kind = 8 )

  el(n:2*n-1) = 0.0D+00
!
!  Direct calculation.
!  Only sum as far as there is data.
!  Simple, but slow.
!
  do j = 0, n-1
    acov(j) = 0.0D+00
    do m = 0, n-1-j
      acov(j) = acov(j) + el(m) * el(m+j)
    end do
  end do
!
!  Write the scaled autocorrelation.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Autocorrelation by the direct method.'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5e14.6)' ) acov(0:19) / acov(0)

  return
end
subroutine test038 ( )

!*****************************************************************************80
!
!! TEST038 tests ZFFTI, ZFFTF and ZFFTB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) wsave(4*n+15)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST038'
  write ( *, '(a)' ) '  For complex fast Fourier transforms, 1D,'
  write ( *, '(a)' ) '  ZFFTI initializes the transforms,'
  write ( *, '(a)' ) '  ZFFTF does a forward transforms;'
  write ( *, '(a)' ) '  ZFFTB does a backward transforms.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  seed = 1973

  call c8vec_uniform_01 ( n, seed, c )

  call c8vec_print_some ( n, c, 1, 10, '  The original data:' )
!
!  Initialize the WSAVE array.
!
  call zffti ( n, wsave )
!
!  Compute the FFT coefficients.
!
  call zfftf ( n, c, wsave )

  call c8vec_print_some ( n, c, 1, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.
!
  call zfftb ( n, c, wsave )
!
!  Normalize the data.
!
  c(1:n) = c(1:n) / real ( n, kind = 8 )
!
!  Now the data should be equal to its input value.
!
  call c8vec_print_some ( n, c, 1, 10, '  The retrieved data:' )

  return
end
subroutine test039 ( )

!*****************************************************************************80
!
!! TEST039 tests ZFFTI, ZFFTF and ZFFTB.
!
!  Discussion:
!
!    Find the autocorrelation to El Nino data using FFT methods.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 168

  real ( kind = 8 ) acov(0:n-1)
  complex ( kind = 8 ) cel(0:2*n-1)
  complex ( kind = 8 ) corr(0:2*n-1)
  real ( kind = 8 ) el(0:2*n-1)
  real ( kind = 8 ) el_sum
  logical ex
  integer ( kind = 4 ) i
  real ( kind = 8 ) wsave(4*(2*n)+15)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST039'
  write ( *, '(a)' ) '  For Fourier transforms of complex data,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ZFFTI initializes,'
  write ( *, '(a)' ) '  ZFFTF forward transforms data,'
  write ( *, '(a)' ) '  ZFFTB backward transforms coefficient.'
!
!  Check to see if the required data file exists.
!
  inquire ( file = 'elnino.dat', exist = ex )

  if ( .not. ex ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST006 - Fatal error!'
    write ( *, '(a)' ) '  Cannot find the data file: elnino.dat '
    return
  end if

  open ( unit = 8, file = 'elnino.dat', status = 'old' )
!
!  Read the data, find the mean of the data.
!
  do i = 0, n-1
    read ( 8, * ) el(i)
  end do

  close ( unit = 8 )

  el_sum = sum ( el(0:n-1) )
!
!  Subtract the mean, and append N zeroes.
!
  el(0:n-1) = el(0:n-1) - el_sum / real ( n, kind = 8 )

  el(n:2*n-1) = 0.0D+00
!
!  Make a complex copy of EL.
!
  cel(0:2*n-1) = cmplx ( el(0:2*n-1), 0.0D+00, kind = 8 )
!
!  Compute FFT of data of length 2*N.
!
!  Compute square of magnitude of transform components and place
!  in complex array as real parts.
!
!  Compute inverse transform, throwing away second half and
!  imaginary parts (which are zero), and multiply by length of
!  sequence, 2*N.
!
  call zffti ( 2*n, wsave )

  call zfftf ( 2*n, cel, wsave )
!
!  ZFFTF returns unscaled transforms.
!  The actual transforms are output divided by 2*N.
!
  corr(0:2*n-1) = abs ( cel(0:2*n-1) / real ( 2 * n, kind = 8 ) ) **2
!
!  Since we compute transform times its conjugate, we must divide by
!  (2n) for each, i.e., (2n)**2.
!
  call zfftb ( 2*n, corr, wsave )

  acov(0:n-1) = real ( corr(0:n-1), kind = 8 ) * real ( 2 * n, kind = 8 )
!
!  Autocovariance is the inverse transform times the sequence length, 2*N.
!
!  Normally, all the scaling would be done once
!  by dividing by 2*N.  We've broken it up for exposition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Autocorrelation  by the complex FFT method.'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5e14.6)' ) acov(0:19) / acov(0)

  return
end
subroutine test040 ( )

!*****************************************************************************80
!
!! TEST040 tests ZFFTF_2D and ZFFTB_2D.
!
!  Discussion:
!
!    Plot the image and transform of an 8 by 8 unit source
!    in a 64 by 64 array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 64
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) dat
  real ( kind = 8 ) ermax
  real ( kind = 8 ) err
  integer ( kind = 4 ) i
  complex ( kind = 8 ) image(lda,n)
  complex ( kind = 8 ) image2(lda,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) wsave(4*n+15)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST040'
  write ( *, '(a)' ) '  For two dimensional complex data:'
  write ( *, '(a)' ) '  ZFFTF_2D computes the forward FFT transform;'
  write ( *, '(a)' ) '  ZFFTB_2D computes the backward FFT transform.'
  write ( *, '(a)' ) ' '
!
!  Initialize WSAVE.
!
  call zffti ( n, wsave )
!
!  Set up the data.
!
!  IMAGE is the original data.
!
!  IMAGE2 is IMAGE scaled by (-1)**(I+J), to place the Fourier coefficients
!  in the correct place for viewing.
!
  do i = 1, n
    do j = 1, n

      if ( i >= (n/2-4) .and. i <= (n/2+4) .and. &
           j >= (n/2-4) .and. j <= (n/2+4) ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if

      image(i,j) = cmplx ( a(i,j), 0.0D+00, kind = 8 )
      image2(i,j) = image(i,j) * (-1)**(i-1+j-1)

    end do
  end do
!
!  Compute the forward Fourier transform of IMAGE and IMAGE2.
!
  call zfftf_2d ( lda, n, image, wsave )

  call zfftf_2d ( lda, n, image2, wsave )
!
!  Compute the magnitude of the components of transforms.
!  The actual transforms are unscaled and need to be divided by N*N
!  to be correct.
!
  a(1:n,1:n) = abs ( image(1:n,1:n) )
!
!  Compute the inverse Fourier transform of IMAGE.
!
  call zfftb_2d ( lda, n, image, wsave )
!
!  The transforms need to be divided by N*N to be correct.
!
  image(1:n,1:n) = image(1:n,1:n) / real ( n**2, kind = 8 )
!
!  See if the result agrees with the original data.
!
  ermax = 0.0D+00

  do i = 1, n
    do j = 1, n

      if ( (n/2-4) <= i .and. i <= (n/2+4) .and. &
           (n/2-4) <= j .and. j <= (n/2+4)) then
        dat = 1.0D+00
      else
        dat = 0.0D+00
      end if

      err = abs ( dat - abs ( image(i,j) ) )

      ermax = max ( ermax, err )

    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Maximum error in ZFFT2D calculation:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6)' ) ermax

  return
end
subroutine test041 ( )

!*****************************************************************************80
!
!! TEST041 tests ZFFTF and ZFFTI.
!
!  Discussion:
!
!    Using complex discrete Fourier transform, find the approximate
!    Fourier coefficients to Runge's function on [-1,1] with N=16 and
!    N=17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 8 ) coeff(0:16)
  real ( kind = 8 ) del
  real ( kind = 8 ) f
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ), external :: runge
  complex ( kind = 8 ) sqtm1
  real ( kind = 8 ) wsave(150)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xj

  sqtm1 = cmplx ( 0.0D+00, -1.0D+00, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST041'
  write ( *, '(a)' ) '  ZFFTI initializes the complex FFT routines.'
  write ( *, '(a)' ) '  ZFFTF does a forward Fourier transform on'
  write ( *, '(a)' ) '    complex data.'
  write ( *, '(a)' ) ' '

  x0 = -1.0D+00

  do n = 16, 17

    call zffti ( n, wsave )
!
!  Function assumed to be periodic on [-1,1], an interval of length 2.
!
    del = 2.0D+00 / real ( n, kind = 8 )
    f = 2.0D+00 * pi / ( real ( n, kind = 8 ) * del )
!
!  First sample point at -1, last at 1-del
!
    do j = 0, n-1
      xj = -1.0D+00 + real ( j, kind = 8 ) * del
      coeff(j) = cmplx ( runge(xj), 0.0D+00, kind = 8 )
    end do

    call zfftf ( n, coeff, wsave )
!
!  Returned coefficients must be divided by N for correct
!  normaliziation.  Also, note repetition after N/2 in original
!  coefficients.  Scaling because X0 not at origin destroys this
!  to some extent.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Results for N = ' , n
    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) ' czero=', coeff(0)/n*2
    write ( *, '(a)' ) &
      ' j          output from ZFFTF,             scaled coeffs'

    do j = 1, n-1
       write ( *, '(2x,i5,2e15.6,5x,2e15.6)' ) &
         j, coeff(j), exp(-sqtm1*j*f*x0) * coeff(j)/n *2
    end do

  end do

  return
end
