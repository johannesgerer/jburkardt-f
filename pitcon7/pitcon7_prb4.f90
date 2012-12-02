program main

!*****************************************************************************80
!
!! PCPRB4 solves a problem involving the Freudenstein-Roth function.
!
!  Modified:
!
!    10 November 1999
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
!  Discussion:
!
!    This version of the Freudenstein-Roth test problem is used to
!    demonstrate the use of the fixed parameterization option.
!
!    Six runs are made:
!
!    IWORK(2)  IWORK(3)  LIM
!
!       2         0       0  Vary index, no limit points.
!       2         1       0  Index=2, no limit points.
!       2         0       1  Vary index, find limit points in index 1.
!       2         1       1  Index=2, find limit points in index 1.
!       2         0       3  Vary index, find limit points in index 3.
!       2         1       3  Index=2, find limit points in index 3.
!
!  The function:
!
!    FX(1) = X1 - X2**3 + 5*X2**2 -  2*X2 - 13 + 34*(X3-1)
!    FX(2) = X1 + X2**3 +   X2**2 - 14*X2 - 29 + 10*(X3-1)
!
!  Starting from the point (15,-2,0), the program is required to produce
!  solution points along the curve until it reaches a solution point
!  (*,*,1).  It also may be requested to look for limit points in the
!  first or third components.
!
!  The correct value of the solution at X3=1 is (5,4,1).
!
!  Limit points in the first variable occur at:
!
!    (14.28309, -1.741377,  0.2585779)
!    (61.66936,  1.983801, -0.6638797)
!
!  Limit points for the third variable occur at:
!
!    (20.48586, -0.8968053, 0.5875873)
!    (61.02031,  2.230139, -0.6863528)
!
  integer, parameter :: nvar=3
  integer, parameter :: liw=nvar+29
  integer, parameter :: lrw=29+(6+nvar)*nvar
!
  external fxroth
  external dfroth
  external dge_slv
  external pitcon
!
  double precision fpar(1)
  integer i
  integer ierror
  integer ipar(1)
  integer itry
  integer iwork(liw)
  integer j
  character ( len = 12 ) name
  double precision rwork(lrw)
  double precision xr(nvar)
!
  call timestamp ( )
  write ( *, * ) ' '
  write ( *, * ) 'PITCON7_PRB4:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  PITCON test problem'
  write ( *, * ) '  Freudenstein-Roth function'
  write ( *, * ) ' '
  write ( *, '(a,i8)' ) '  Number of equations is ', nvar - 1
  write ( *, '(a,i8)' ) '  Number of variables is ', nvar

  itry = 0

10    continue

  itry=itry+1
  write ( *, * ) ' '
  write ( *, * ) 'This is run number ',itry
!
!  Set work arrays to zero:
!
  iwork(1:liw) = 0
  rwork(1:lrw) = 0.0
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
!  IWORK(9)=2 ; Jacobian choice.
!
  iwork(1)=0
  iwork(2)=2

  if(mod(itry,2)==1)then
    iwork(3)=0
    write ( *, * ) 'PITCON is free to choose parameterization.'
  else
    iwork(3)=1
    write ( *, * ) 'The user fixes the parameterization.'
  end if

  iwork(4)=0
  iwork(5)=3
  if(itry==1.or.itry==2)then
    iwork(6)=0
    write ( *, * ) 'No limit points are sought.'
  elseif(itry==3.or.itry==4)then
    iwork(6)=1
    write ( *, * ) 'Seek limit points in index ',iwork(6)
  elseif(itry==5.or.itry==6)then
    iwork(6)=3
    write ( *, * ) 'Seek limit points in index ',iwork(6)
  end if

  iwork(7)=1
  iwork(9)=0
!
!  RWORK(1)=0.00001; Absolute error tolerance
!  RWORK(2)=0.00001; Relative error tolerance
!  RWORK(3)=0.01   ; Minimum stepsize
!  RWORK(4)=20.0   ; Maximum stepsize
!  RWORK(5)=0.3    ; Starting stepsize
!  RWORK(6)=1.0    ; Starting direction
!  RWORK(7)=1.0    ; Target value (Seek solution with X(3)=1)
!
  rwork(1)=0.00001
  rwork(2)=0.00001
  rwork(3)=0.01
  rwork(4)=10.0
  rwork(5)=0.3
  rwork(6)=1.0
  rwork(7)=1.0
!
!  Set starting point.
!
  xr(1)=15.0
  xr(2)=-2.0
  xr(3)=0.0

  write ( *, * ) ' '
  write ( *, * ) 'Step  Type of point     X(1)         X(2)         X(3)'
  write ( *, * ) ' '
  i=0
  name='Start point  '
  write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,(xr(j),j=1,nvar)

  do i=1,30

    call pitcon(dfroth,fpar,fxroth,ierror,ipar,iwork,liw, &
      nvar,rwork,lrw,xr,dge_slv)

    if(iwork(1)==1)then
      name='Corrected    '
    elseif(iwork(1)==2)then
      name='Continuation '
    elseif(iwork(1)==3)then
      name='Target point '
    elseif(iwork(1)==4)then
      name='Limit point  '
    elseif(iwork(1)<0)then
      name='Jacobian   '
    end if

    write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,(xr(j),j=1,nvar)

    if(iwork(1)==3)then
      write ( *, * ) ' '
      write ( *, * ) 'We have reached the point we wanted.'
      write ( *, * ) 'The code may stop now.'
      go to 60
    end if

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PITCON returned an error code:'
      write ( *, * ) 'IERROR = ', ierror
      write ( *, * ) ' '
      write ( *, * ) 'The computation failed.'
      go to 50
    end if

  end do

50    continue

  write ( *, * ) ' '
  write ( *, * ) 'PITCON did not reach the point of interest.'

60    continue

  if(itry<6)go to 10
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PITCON7_PRB4:'
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
!  The function has the form:
!
!    ( X1 - ((X2-5.0)*X2+2.0)*X2 - 13.0 + 34.0*(X3-1.0)  )
!    ( X1 + ((X2+1.0)*X2-14.0)*X2 - 29.0 + 10.0*(X3-1.0) )
!
  integer nvar

  double precision f(*)
  double precision fpar(*)
  integer ipar(*)
  double precision x(nvar)

  f(1) = x(1) - ( ( x(2) - 5.0 ) * x(2) + 2.0 ) * x(2) - 13.0 &
         + 34.0 * ( x(3) - 1.0 )

  f(2) = x(1) + ( ( x(2) + 1.0 ) * x(2) - 14.0 ) * x(2) - 29.0 &
         + 10.0 * ( x(3) - 1.0 )

  return
end
subroutine dfroth ( nvar, fpar, ipar, x, fjac )

!*****************************************************************************80
!
!! DFROTH evaluates the Jacobian J(X) at X.
!
!  The jacobian has the form:
!
!    ( 1.0   (-3.0*X(2)+10.0)*X(2)- 2.0   34.0  )
!    ( 1.0   ( 3.0*X(2)+ 2.0)*X(2)-14.0   10.0  )
!
  integer nvar

  double precision fjac(nvar,nvar)
  double precision fpar(*)
  integer ipar(*)
  double precision x(nvar)

  fjac(1,1) = 1.0
  fjac(1,2) = ( - 3.0 * x(2) + 10.0 ) * x(2) - 2.0
  fjac(1,3) = 34.0

  fjac(2,1) = 1.0
  fjac(2,2) = ( 3.0 * x(2) + 2.0 ) * x(2) - 14.0
  fjac(2,3) = 10.0

  return
end
