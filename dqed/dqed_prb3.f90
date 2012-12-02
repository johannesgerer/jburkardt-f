program main

!*****************************************************************************80
!
!! MAIN is the main program for DQED_PRB3.
!
!  Discussion:
!
!    DQED_PRB3 tests DQED.
!
!    This program illustrates the use of the Hanson-Krogh nonlinear 
!    least squares solver DQED for fitting two exponentials to data.
! 
!    The problem is to find values for the four variables x(1),...,x(4),
!    which specify the model function
! 
!      h(t) = x(1)*exp(x(2)*t) + x(3)*exp(x(4)*t)
!
!    which mininize the sum of the squares of the error between h(t) and
!    recorded values of h at five values of t:
!
!      t = 0.05, 0.1, 0.4, 0.5, and 1.0.
!
!    We also have problem constraints that 
!
!          0 < =  x(1)
!      -25.0 < =  x(2) <= 0
!          0 < =  x(3)
!      -25.0 < =  x(4) <= 0. 
!
!    and a minimal separation requirement between x(2) and x(4) that
!    can be expressed as
!
!       0.05 < =  x(2) - x(4)
!
!
!  Output:
! 
!    model is h(t) = x(1)*exp(-t*x(2)) + x(3)*exp(t*x(4))
!    x(1),x(2),x(3),x(4)  = 
!        1.999475    -.999801     .500057   -9.953988
!     residual after the fit =   4.2408d-04
!     number of evaluations of parameter model =    14
!     output flag from solver =                      4
! 
  implicit none

  integer ( kind = 4 ), parameter :: liwork = 84
  integer ( kind = 4 ), parameter :: lwork = 640
  integer ( kind = 4 ), parameter :: mcon = 1
  integer ( kind = 4 ), parameter :: mequa = 5
  integer ( kind = 4 ), parameter :: nvars = 4

  integer ( kind = 4 ), parameter :: ldfj = mcon + mequa

  real ( kind = 8 ) bl(nvars+mcon)
  real ( kind = 8 ) bu(nvars+mcon)
  real ( kind = 8 ), dimension ( mequa ) :: f = (/ &
    2.206D+00, 1.994D+00, 1.350D+00, 1.216D+00, 0.7358D+00 /)
  real ( kind = 8 ) fj(ldfj,nvars+1)
  real ( kind = 8 ) fnorm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) igo
  integer ( kind = 4 ) ind(nvars+mcon)
  integer ( kind = 4 ) iopt(24)
  integer ( kind = 4 ) iwork(liwork)
  integer ( kind = 4 ) niters
  real ( kind = 8 ) ropt(1)
  external dqedev
  real ( kind = 8 ), dimension ( mequa ) :: t = (/ &
    0.05D+00, 0.10D+00, 0.40D+00, 0.50D+00, 1.00D+00 /)
  real ( kind = 8 ) work(lwork)
  real ( kind = 8 ) x(nvars)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQED_PRB3'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A test for DQED'
!
!  Define the bounding constraints on the variables X(1) through X(4),
!  and then the linear constraints on the separation condition.
!
!  1: 0.0 < =  X(1)
!
  ind(1) = 1
  bl(1) = 0.0D+00
  bu(1) = 0.0D+00
!
!  2: -25.0 < =  X(2) <= 0.0
!
  ind(2) = 3
  bl(2) = -25.0D+00
  bu(2) = 0.0D+00
!
!  3: 0.0 < =  X(3)
!
  ind(3) = 1
  bl(3) = 0.0D+00
  bu(3) = 0.0D+00
!
!  4: -25.0 < =  X(4) <= 0.0
!
  ind(4) = 3
  bl(4) = -25.0D+00
  bu(4) = 0.0D+00
!
!  5: 0.05 < =  X(2) - X(4)
!
  ind(5) = 1
  bl(5) = 0.05D+00
  bu(5) = 0.0D+00
!
!  Define the initial values of the variables.
!  We don't know anything at all, so all variables are set zero.
!
  x(1:nvars) = 0.0D+00
!
!  Tell how much storage we gave the solver.
!
  iwork(1) = lwork
  iwork(2) = liwork
!
!  Initialize the call counter.
!
  niters = 0
!
!  Use reverse commumication to evaluate the derivatives.
!
  iopt(1) = 16
  iopt(2) = 1
!
!  No more options.
!
  iopt(3) = 99

  do

    call dqed ( dqedev, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, fnorm, &
      igo, iopt, ropt, iwork, work )

    if ( 1 < igo ) then
      exit
    end if
!
!  Count function evaluations.
!
    niters = niters + 1
! 
!  Set the value of the constraint, and the residual functions.
!
    fj(1,5) = x(2) - x(4)

    do i = 1, mequa

      fj(mcon+i,5) = x(1)*exp(x(2)*t(i)) + x(3)*exp(x(4)*t(i)) - f(i)

    end do
!
!  If IGO is nonzero, compute the derivatives.
!
    if ( igo /= 0 ) then
 
      fj(1,1) = 0.0D+00
      fj(1,2) = 1.0D+00
      fj(1,3) = 0.0D+00
      fj(1,4) = -1.0D+00

      do i = 1,mequa

        fj(mcon+i,1) = exp(x(2)*t(i))
        fj(mcon+i,2) = x(1)*t(i)*exp(x(2)*t(i))
        fj(mcon+i,3) = exp(x(4)*t(i))
        fj(mcon+i,4) = x(3)*t(i)*exp(x(4)*t(i))

      end do

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Model:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    h(t) = x(1)*exp(-t*x(2)) + x(3)*exp(t*x(4))'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Constrained minimizing X:'
  write ( *, '(a)' ) ' '
  write ( *, '(g14.6)' ) x(1:nvars)
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) ' Residual of the fitted data, FNORM = ',fnorm
  write ( *, '(a,i6)' ) '  Number of model evaluations NITERS = ',niters
  write ( *, '(a,i6)' ) '  DQED output flag IGO = ',igo
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQED_PRB3'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
