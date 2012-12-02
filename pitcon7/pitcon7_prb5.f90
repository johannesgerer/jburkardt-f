program main

!*****************************************************************************80
!
!! PCPRB5 solves a discretized two point boundary value problem with parameter.
!
!  Modified:
!
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!  The function:
!
!    Y'' + LAMBDA * EXP ( Y ) = 0
!    Y(0) = 0
!    Y'(1) = 0
!
!  Discussion:
!
!    This test compares the use of:
!
!    * the full Newton method, which updates the Jacobian on every
!      Newton iteration,
!
!    * the modified Newton method, which updates the Jacobian only at
!      specified steps in the Newton iteration,
!
!    * the "cheap" Newton method, which updates the Jacobian only
!      when convergence fails.
!
  integer, parameter :: nvar = 22
  integer, parameter :: lrw = 29 + ( nvar + 6 ) * nvar
  integer, parameter :: liw = nvar + 29

  external dge_slv
  external fpfull
  external fxbend
  external pitcon

  double precision fpar(1)
  integer i
  integer ierror
  integer ipar(2)
  integer itry
  integer iwork(liw)
  integer j
  integer modcon
  character ( len = 12 ) name
  integer nit
  double precision rwork(lrw)
  double precision xr(nvar)

  call timestamp ( )
  write ( *, * ) ' '
  write ( *, * ) 'PITCON7_PRB5'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  PITCON test problem'
  write ( *, * ) '  Two point boundary value problem.'
  write ( *, * ) ' '
  write ( *, '(a,i8)' ) '  Number of equations is ', nvar - 1
  write ( *, '(a,i8)' ) '  Number of variables is ', nvar
  write ( *, * ) ' '
  write ( *, * ) '  This program will be run three times.'

  itry = 0

10    continue

  itry = itry + 1
  write(*,*)' '
  write(*,*)'This is run number ',itry
  nit=0

  write(*,*)' '

  iwork(1:liw) = 0
  rwork(1:lrw) = 0.0

  ierror=0
!
!  Set input quantities
!
  if(itry==1)modcon=0
  if(itry==2)modcon=1
  if(itry==3)modcon=2
!
!  IWORK(1)=0      ; This is a startup
!  IWORK(2)=NVAR   ; Use X(NVAR) for initial parameter
!  IWORK(3)=0      ; Program may choose parameter index
!  IWORK(4)=MODCON ; Control frequency of jacobian update in newton iteration
!  IWORK(5)=NVAR   ; Seek target values for X(NVAR)
!  IWORK(6)=NVAR   ; Seek limit points in X(NVAR)
!  IWORK(7)=1      ; Control amount of output.
!  IWORK(9)=0      ; Jacobian choice.
!
  iwork(1)=0
  iwork(2)=nvar
  iwork(3)=0
  iwork(4)=modcon
  iwork(5)=nvar
  iwork(6)=nvar
  iwork(7)=1
  iwork(9)=0
!
!  Jacobian upper/lower bandwidths:
!
  ipar(1)=1
  ipar(2)=1
!
!  RWORK(1)=0.0001; Absolute error tolerance
!  RWORK(2)=0.0001; Relative error tolerance
!  RWORK(3)=0.01   ; Minimum stepsize
!  RWORK(4)=0.25   ; Maximum stepsize
!  RWORK(5)=0.05   ; Starting stepsize
!  RWORK(6)=1.0    ; Starting direction
!  RWORK(7)=0.80   ; Target value (Seek solution with X(NVAR)=0.80
!
  rwork(1)=0.0001
  rwork(2)=0.0001
  rwork(3)=0.01
  rwork(4)=0.25
  rwork(5)=0.05
  rwork(6)=1.0
  rwork(7)=0.80
!
!  Set starting point
!
  do i=1,nvar
    xr(i)=0.0
  end do

  if(itry==1)then
    write(*,*)' '
    write(*,*)'Using the full Newton method.'
    write(*,*)'The jacobian is updated on every iteration.'
  elseif(itry==2)then
    write(*,*)' '
    write(*,*)'Using the modified Newton method.'
    write(*,*)'The jacobian is held fixed while correcting a '
    write(*,*)'point.'
  else
    write(*,*)' '
    write(*,*)'Using the "cheap" Newton method.'
    write(*,*)'The jacobian is held fixed as long as possible,'
    write(*,*)'perhaps over multiple points, until'
    write(*,*)'convergence fails.'
  end if

  write(*,*)' '
  write(*,*)'Step  Type of point     Lambda'
  write(*,*)' '

  i=0
  name='Start point  '
  write(*,'(1x,i3,2x,a12,2x,g14.6)')i,name,xr(nvar)

  do i=1,50

    call pitcon(fpfull,fpar,fxbend,ierror,ipar,iwork,liw, &
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

    write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,xr(nvar)

    if(iwork(1)==3)then
      nit=nit+1
      write(*,*)' '
      write(*,*)'Complete value for target point ',nit
      write(*,*)' '
      write(*,'(1x,5g14.6)')(xr(j),j=1,nvar)
      if(nit>=2)go to 60
    end if

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PITCON returned an error code:'
      write ( *, * ) 'IERROR = ', ierror
      write ( *, * ) ' '
      write ( *, * ) 'The computation failed.'
      stop
    end if

  end do

60    continue

  write(*,*)' '
  write(*,*)'Jacobians      ',iwork(19)
  write(*,*)'Factorizations ',iwork(20)
  write(*,*)'Solves         ',iwork(21)
  write(*,*)'Functions      ',iwork(22)

  if(itry<3)go to 10
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PITCON7_PRB5:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine fxbend(nvar,fpar,ipar,x,fx )

!*****************************************************************************80
!
!! FXBEND evaluates the nonlinear function for two point boundary value problem.
!
  integer nvar

  double precision fpar(1)
  double precision fx(nvar)
  double precision hsq
  integer i
  integer ipar(2)
  double precision x(nvar)

  hsq=1.0/(nvar-2)**2

  fx(1)=x(1)

  do i=2,nvar-2
    fx(i)=x(i-1)-2.0*x(i)+x(i+1)+x(nvar)*exp(x(i))*hsq
  end do

  fx(nvar-1)=x(nvar-1)-x(nvar-2)

  return
end
subroutine fpfull ( nvar, fpar, ipar, x, fprime )

!*****************************************************************************80
!
!! FPFULL evaluates the full storage jacobian.
!
  integer nvar

  double precision fpar(*)
  double precision fprime(nvar,nvar)
  double precision hsq
  integer i
  integer ipar(2)
  double precision x(nvar)

  hsq = 1.0 / dble ( nvar - 2 )**2

  fprime(1,1) = 1.0

  do i=2,nvar-2
    fprime(i,i-1) = 1.0
    fprime(i,i) = - 2.0 + x(nvar) * exp ( x(i) ) * hsq
    fprime(i,i+1) = 1.0
    fprime(i,nvar) = exp ( x(i) ) * hsq
  end do

  fprime(nvar-1,nvar-2) = - 1.0
  fprime(nvar-1,nvar-1) = 1.0

  return
end
