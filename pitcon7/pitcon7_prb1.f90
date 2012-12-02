program main

!*****************************************************************************80
!
!! MAIN is the main program for PITCON7_PRB1.
!
!  Discussion:
!
!    PITCON66_PRB1 treats a system based on the Freudenstein-Roth function.
!
!    The function F(X) is of the form
!
!      FX(1) = X1 - X2**3 + 5*X2**2 -  2*X2 - 13 + 34*(X3-1)
!      FX(2) = X1 + X2**3 +   X2**2 - 14*X2 - 29 + 10*(X3-1)
!
!    Starting from the point (15,-2,0), the program is required to produce
!    solution points along the curve until it reaches a solution point
!    (*,*,1).  It also may be requested to look for limit points in the
!    first or third components.
!
!    The correct value of the solution at X3=1 is (5,4,1).
!
!    Limit points in the first variable occur at:
!
!      (14.28309, -1.741377,  0.2585779)
!      (61.66936,  1.983801, -0.6638797)
!
!    Limit points for the third variable occur at:
!
!      (20.48586, -0.8968053, 0.5875873)
!      (61.02031,  2.230139, -0.6863528)
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    F Freudenstein, B Roth,
!    Numerical Solutions of Nonlinear Equations,
!    Journal of the Association for Computing Machinery,
!    Volume 10, 1963, Pages 550-556.
!
  implicit none

  integer, parameter :: nvar = 3
  integer, parameter :: liw = nvar + 29
  integer, parameter :: lrw = 29 + ( 6 + nvar ) * nvar

  external dfroth
  external dge_slv
  double precision fpar(1)
  external fxroth
  integer i
  integer ierror
  integer ipar(1)
  integer iwork(liw)
  integer j
  character ( len = 12 ) name
  double precision rwork(lrw)
  double precision xr(nvar)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PITCON7_PRB1:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PITCON test problem'
  write ( *, '(a)' ) '  Freudenstein-Roth function'
!
!  Set work arrays to zero:
!
  iwork(1:liw) = 0
  rwork(1:lrw) = 0.0E+00
!
!  Set some entries of work arrays.
!
!  IWORK(1)=0 ; This is a startup
!  IWORK(2)=2 ; Use X(2) for initial parameter
!  IWORK(3)=0 ; Program may choose parameter index
!  IWORK(4)=0 ; Update jacobian every newton step
!  IWORK(5)=3 ; Seek target values for X(3)
!  IWORK(6)=1 ; Seek limit points in X(1)
!  IWORK(7)=1 ; Control amount of output.
!  IWORK(9)=0 ; Jacobian choice.
!
  iwork(1) = 0
  iwork(2) = 2
  iwork(3) = 0
  iwork(4) = 0
  iwork(5) = 3
  iwork(6) = 1
  iwork(7) = 3
  iwork(9) = 0
!
!  RWORK(1)=0.00001; Absolute error tolerance
!  RWORK(2)=0.00001; Relative error tolerance
!  RWORK(3)=0.01   ; Minimum stepsize
!  RWORK(4)=10.0   ; Maximum stepsize
!  RWORK(5)=0.3    ; Starting stepsize
!  RWORK(6)=1.0    ; Starting direction
!  RWORK(7)=1.0    ; Target value (Seek solution with X(3)=1)
!
  rwork(1) = 0.00001E+00
  rwork(2) = 0.00001E+00
  rwork(3) = 0.01E+00
  rwork(4) = 10.0E+00
  rwork(5) = 0.3E+00
  rwork(6) = 1.0E+00
  rwork(7) = 1.0E+00
!
!  Set the starting point.
!
  xr(1:3) =  (/ 15.0E+00, -2.0E+00, 0.0E+00 /)

  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of equations is ', nvar - 1
  write ( *, '(a,i8)' ) '  Number of variables is ', nvar
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Step  Type of point     X(1)         X(2)         X(3)'
  write ( *, '(a)'    ) ' '

  i = 0
  name = 'Start point  '
  write ( *, '(i3,2x,a12,2x,3g14.6)' ) i, name, xr(1:nvar)

  do i = 1, 30

    call pitcon ( dfroth, fpar, fxroth, ierror, ipar, iwork, liw, &
      nvar, rwork, lrw, xr, dge_slv )

    if ( iwork(1) == 1 ) then
      name = 'Corrected    '
    else if ( iwork(1) == 2 ) then
      name = 'Continuation '
    else if ( iwork(1) == 3 ) then
      name = 'Target point '
    else if ( iwork(1) == 4 ) then
      name = 'Limit point  '
    else if ( iwork(1) < 0 ) then
      name = 'Jacobian   '
    end if

    write ( *, '(i3,2x,a12,2x,3g14.6)' ) i, name, xr(1:nvar)

    if ( iwork(1) == 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PITCON reached the target point.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'The computation succeeded.'
      exit
    end if

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PITCON returned an error code:'
      write ( *, '(a,i6)' ) 'IERROR = ', ierror
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'The computation failed.'
      exit
    end if

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PITCON7_PRB1:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine fxroth ( nvar, fpar, ipar, x, f )

!*****************************************************************************80
!
!! FXROTH evaluates the function F(X) at X.
!
!  Function:
!
!    ( X1 - ((X2-5.0)*X2+2.0)*X2 - 13.0 + 34.0*(X3-1.0)  )
!    ( X1 + ((X2+1.0)*X2-14.0)*X2 - 29.0 + 10.0*(X3-1.0) )
!
  implicit none

  integer nvar

  double precision f(*)
  double precision fpar(*)
  integer ipar(*)
  double precision x(nvar)

  f(1) = x(1) - ( ( x(2) - 5.0E+00 ) * x(2) + 2.0E+00 ) * x(2) - 13.0E+00 &
         + 34.0E+00 * ( x(3) - 1.0E+00 )

  f(2) = x(1) + ( ( x(2) + 1.0E+00 ) * x(2) - 14.0E+00 ) * x(2) - 29.0E+00 &
         + 10.0E+00 * ( x(3) - 1.0E+00 )

  return
end
subroutine dfroth ( nvar, fpar, ipar, x, fjac )

!*****************************************************************************80
!
!! DFROTH evaluates the Jacobian J(X) at X.
!
!  Jacobian:
!
!    ( 1.0   (-3.0*X(2)+10.0)*X(2)- 2.0   34.0  )
!    ( 1.0   ( 3.0*X(2)+ 2.0)*X(2)-14.0   10.0  )
!
  implicit none

  integer nvar

  double precision fjac(nvar,nvar)
  double precision fpar(*)
  integer ipar(*)
  double precision x(nvar)

  fjac(1,1) = 1.0E+00
  fjac(1,2) = ( - 3.0E+00 * x(2) + 10.0E+00 ) * x(2) - 2.0E+00
  fjac(1,3) = 34.0E+00

  fjac(2,1) = 1.0E+00
  fjac(2,2) = ( 3.0E+00 * x(2) + 2.0E+00 ) * x(2) - 14.0E+00
  fjac(2,3) = 10.0E+00

  return
end
