program main

!*****************************************************************************80
!
!! DQED_PRB4 tests DQED.
!
!  Discussion:
!
!    This program illustrates the use of the Hanson-Krogh nonlinear least
!    squares solver DQED.
! 
!    The variables X(1:6) represent the dietary proportion of 6 prey
!    animals in the diet of an alligator.
!
!    The percentages of isotopic Carbon 13 and Nitrogen 15 are known
!    for each of these prey animals, as well as for the alligator.
!
!    The problem is to find values for the variables x(1:6)
!    so that the carbon and nitrogen levels in the alligator can be
!    explained purely in terms of diet.  In other words, we know
!    C1 through C6 and CA, N1 through N6 and NA, and we seek X1 through
!    X6 so that, as closely as possible, it is true that:
!
!                              (X1)
!                              (X2)
!      ( C1 C2 C3 C4 C5 C6 ) * (X3)  = (CA NA)
!      ( N1 N2 N3 N4 N5 N6 )   (X4)
!                              (X5)
!                              (X6)
!
!    We have constraints on our variables, namely, that each X is nonnegative,
!    and that the X's sum to 1.
!
!    The impetus for this problem, and the field data, came from
!
!    Matt Aresco,
!    Biology Department,
!    Florida State University,
!    aresco@bio.fsu.edu
!
!  Modified:
!
!    21 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Laurel Saito, Brett Johnson, John Bartholow, Blair Hanna,
!    Assessing Ecosystem Effects of Reservoir Operations Using Food
!      Web-Energy Transfer and Water Quality Models,
!    Ecosystems, 
!    Volume 4, Number 2, pages 105-125, 2001.
!
  implicit none

  integer ( kind = 4 ), parameter :: mcon = 1
  integer ( kind = 4 ), parameter :: mequa = 2
  integer ( kind = 4 ), parameter :: npmax = 5
  integer ( kind = 4 ), parameter :: nvars = 6

  integer ( kind = 4 ), parameter :: nall = mcon + 2 * nvars + npmax + 1

  integer ( kind = 4 ), parameter :: nplus = 3 * nall + 2

  integer ( kind = 4 ), parameter :: liwork = 3 * mcon + 9 * nvars + 4 * npmax + nall + 11
  integer ( kind = 4 ), parameter :: ldfj = mcon + mequa
  integer ( kind = 4 ), parameter :: lwork = nall * nall + 4 * nall + nvars * npmax &
    + 33 * nvars + mequa * npmax + nvars * npmax + 13 * npmax &
    + 9 * mcon + 26 + nplus

  real ( kind = 8 ) bl(nvars+mcon)
  real ( kind = 8 ) bu(nvars+mcon)
  real ( kind = 8 ), dimension ( nvars ) :: c = (/ &
    -27.60D+00, -25.40D+00, -25.63D+00, -27.17D+00, -24.81D+00, -26.22D+00 /)
  real ( kind = 8 ) :: c_total = -25.34
  real ( kind = 8 ) fj(ldfj,nvars+1)
  real ( kind = 8 ) fnorm
  integer ( kind = 4 ) igo
  integer ( kind = 4 ) ind(nvars+mcon)
  integer ( kind = 4 ) iopt(24)
  integer ( kind = 4 ) iwork(liwork)
  real ( kind = 8 ), dimension ( nvars ) :: n = (/ &
    2.79D+00, 6.93D+00, 8.21D+00, 8.34D+00, 1.97D+00, 4.36D+00 /)
  real ( kind = 8 ) :: n_total = 5.10D+00
  integer ( kind = 4 ) niters
  real ( kind = 8 ) ropt(1)
  external dqedev
  real ( kind = 8 ) work(lwork)
  real ( kind = 8 ) x(nvars)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQED_PRB4'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A test for DQED'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This example uses reverse communication.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The stable isotope problem for the alligator.'
  write ( *, '(a)' ) '  We assume an alligator has 6 kinds of prey, and '
  write ( *, '(a)' ) '  that the concentrations of isotopic C13 and N15 '
  write ( *, '(a)' ) '  are known in the prey and in the alligator.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We want to determine dietary proportions X1 through X6'
  write ( *, '(a)' ) '  so that the isotope values in the alligator can '
  write ( *, '(a)' ) '  be explained'
  write ( *, '(a)' ) '  in terms of how much of what prey it eats.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our model:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    (C_All) = (C1 C2 C3 C4 C5 C6 ) * (X1 X2 X3 X4 X5 X6)'''
  write ( *, '(a)' ) '    (N_All) = (N1 N2 N3 N4 N5 N6 )'
!
!  Define the bounding constraints on the variables.
!
!  1: 0.0 <= X(1)
!
  ind(1) = 1
  bl(1) = 0.0D+00
  bu(1) = 0.0D+00
!
!  2: 0.0 <= X(2)
!
  ind(2) = 1
  bl(2) = 0.0D+00
  bu(2) = 0.0D+00
!
!  3: 0.0 <= X(3)
!
  ind(3) = 1
  bl(3) = 0.0D+00
  bu(3) = 0.0D+00
!
!  4: 0.0 <= X(4)
!
  ind(4) = 1
  bl(4) = 0.0D+00
  bu(4) = 0.0D+00
!
!  5: 0.0 <= X(5)
!
  ind(5) = 1
  bl(5) = 0.0D+00
  bu(5) = 0.0D+00
!
!  6: 0.0 <= X(6)
!
  ind(6) = 1
  bl(6) = 0.0D+00
  bu(6) = 0.0D+00
!
!  7: the linear constraint 1 = X(1) + X(2) + X(3) + X(4) + X(5) + X(6)
!
  ind(7) = 3
  bl(7) = 1.0D+00
  bu(7) = 1.0D+00
!
!  Define the initial values of the variables.
!  We don't know anything at all, so all variables are set equal.
!
  x(1:nvars) = 1.0D+00 / real ( nvars, kind = 8 )
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
!  Set the value of the constraint.
!
    fj(1,nvars+1) = sum ( x(1:nvars) )
! 
!  Set the value of the residual functions.
!
    fj(mcon+1,nvars+1) = dot_product ( x(1:nvars), c(1:nvars) ) - c_total
    fj(mcon+2,nvars+1) = dot_product ( x(1:nvars), n(1:nvars) ) - n_total
!
!  If IGO is nonzero, compute the derivatives.
!
    if ( igo == 0 ) then
      cycle
    end if
!
!  Partial derivatives of the constraint equation.
!
    fj(1,1:nvars) = 1.0D+00
!
!  Partial derivatives of the least squares equations.
!
    fj(mcon+1,1:nvars) = c(1:nvars)
    fj(mcon+2,1:nvars) = n(1:nvars)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our computed answer X:'
  write ( *, '(a)' ) ' '
  write ( *, '(g14.6)' ) x(1:nvars)
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The residual of the fitted data, FNORM = ', fnorm
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Sum of variables = ', sum ( x(1:nvars) )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of model evaluations NITERS = ', niters
  write ( *, '(a,i6)' ) '  DQED output flag IGO = ', igo
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQED_PRB4'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
