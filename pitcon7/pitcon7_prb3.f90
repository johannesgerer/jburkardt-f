program main

!*****************************************************************************80
!
!! PCPRB3 solves a second order boundary value problem with a parameter.
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
!  We use a finite difference approximation, with 21 equally spaced
!  nodes.  The value of Y at node (I) is X(I), for I=1 to 21.
!  X(22) is the value of LAMBDA.
!
!  We expect a limit point in LAMBDA at roughly LAMBDA=0.878.
!
!  NVAR-1 is the number of grid points between 0 and 1.  This problem may
!  be set up with arbitrarily many grid points, simply by increasing the
!  value of NVAR.  No other change is required.
!
!
!  This problem is solved six times:
!
!    With full storage user jacobian,
!    With full storage forward difference jacobian,
!    With full storage central difference jacobian.
!    With band storage user jacobian,
!    With band storage forward difference jacobian,
!    With band storage central difference jacobian.
!
  integer, parameter :: nvar = 22
  integer, parameter :: lrw = 29 + ( nvar + 6 ) * nvar
  integer, parameter :: liw = nvar + 29

  external dgb_slv
  external dge_slv
  external fpband
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
  integer jac
  character ( len = 12 ) name
  integer nit
  double precision rwork(lrw)
  double precision xr(nvar)

  call timestamp ( )
  write ( *, * ) ' '
  write ( *, * ) 'PITCON7_PRB3'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  PITCON sample program.'
  write ( *, * ) '  Two point boundary value problem.'
  write ( *, * ) ' '
  write ( *, '(a,i8)' ) '  Number of equations is ', nvar - 1
  write ( *, '(a,i8)' ) '  Number of variables is ', nvar
  write ( *, * ) ' '
  write ( *, * ) '  Seek limit points in lambda.'
  write ( *, * ) ' '
  write ( *, * ) '  This problem will be run six times.'
  write ( *, * ) ' '

  itry=0

10    continue

  itry=itry+1

  write ( *, * ) ' '
  write ( *, * ) 'This is run number ', itry
  write ( *, * ) ' '

  iwork(1:liw) = 0
  rwork(1:lrw) = 0.0

  ierror=0

  nit=0
!
!  For runs 1 and 4, supply Jacobian.
!  For runs 2 and 5, use forward difference approximation,
!  For runs 3 and 6, use central difference approximation.
!
  if(mod(itry,3)==1)jac=0
  if(mod(itry,3)==2)jac=1
  if(mod(itry,3)==0)jac=2
!
!  IWORK(1)=0    ; This is a startup
!  IWORK(2)=NVAR ; Use X(NVAR) for initial parameter
!  IWORK(3)=0    ; Program may choose parameter index
!  IWORK(4)=1    ; Cut down evaluations of jacobian on newton steps.
!  IWORK(5)=NVAR ; Seek target values for X(NVAR)
!  IWORK(6)=NVAR ; Seek limit points in X(NVAR)
!  IWORK(7)=1    ; Control amount of output.
!  IWORK(9)=*    ; Jacobian choice.  0=user, 1=forward, 2=central.
!
  iwork(1)=0
  iwork(2)=nvar
  iwork(3)=0
  iwork(4)=1
  iwork(5)=nvar
  iwork(6)=nvar
  iwork(7)=1
  iwork(9)=jac
!
!  Pass the lower and upper bandwidths of the Jacobian
!  in the integer parameter array IPAR.
!
  ipar(1)=1
  ipar(2)=1
!
!  RWORK(1)=0.0001 ; Absolute error tolerance
!  RWORK(2)=0.0001 ; Relative error tolerance
!  RWORK(3)=0.01   ; Minimum stepsize
!  RWORK(4)=0.25   ; Maximum stepsize
!  RWORK(5)=0.05   ; Starting stepsize
!  RWORK(6)=1.0    ; Starting direction
!  RWORK(7)=0.80   ; Target value (Seek solution with X(NVAR)=0.80)
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
  do i = 1, nvar
    xr(i) = 0.0
  end do

  if(itry<=3)then
    write ( *, * ) 'This run uses the full linear solver DENSLV.'
  else
    write ( *, * ) 'This run uses the banded linear solver BANSLV.'
  end if

  if(jac==0)then
    write ( *, * ) 'This run assumes that the user supplies the'
    write ( *, * ) 'jacobian matrix via a subroutine.'
  else if (jac==1)then
    write ( *, * ) 'This run assumes that PITCON will approximate'
    write ( *, * ) 'the jacobian using forward differences.'
  else if ( jac == 2 ) then
    write ( *, * ) 'This run assumes that PITCON will approximate'
    write ( *, * ) 'the jacobian using central differences.'
  end if

  write ( *, * ) ' '
  write ( *, * ) 'Step  Type of point     Lambda'
  write ( *, * ) ' '
  i=0
  name='Start point  '
  write(*,'(1x,i3,2x,a12,2x,g14.6)')i,name,xr(nvar)
!
!  Take steps along the curve
!
  do i=1,50

    if(itry<=3)then
      call pitcon(fpfull,fpar,fxbend,ierror,ipar,iwork,liw, &
        nvar,rwork,lrw,xr,dge_slv)
    else
      call pitcon(fpband,fpar,fxbend,ierror,ipar,iwork,liw, &
        nvar,rwork,lrw,xr,dgb_slv)
    end if

    if(iwork(1)==1)then
      name='Corrected    '
    else if (iwork(1)==2)then
      name='Continuation '
    else if (iwork(1)==3)then
      name='Target point '
    else if (iwork(1)==4)then
      name='Limit point  '
    else if (iwork(1)<0)then
      name='Jacobian   '
    end if

    write(*,'(1x,i3,2x,a12,2x,g14.6)')i,name,xr(nvar)

    if(iwork(1)==3)then
      write ( *, * ) ' '
      nit=nit+1
      write ( *, * ) 'Value of target point ',nit
      write(*,'(1x,5g14.6)')(xr(j),j=1,nvar)
      write ( *, * ) ' '
      if(nit>=2)go to 60
    end if

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PITCON returned an error code:'
      write ( *, * ) 'IERROR = ', ierror
      write ( *, * ) ' '
      write ( *, * ) 'The computation failed.'
      go to 60
    end if

  end do

60    continue

  write ( *, * ) ' '
  write ( *, * ) 'Jacobians:      ',iwork(19)
  write ( *, * ) 'Factorizations: ',iwork(20)
  write ( *, * ) 'Solves:         ',iwork(21)
  write ( *, * ) 'Functions:      ',iwork(22)

  if ( itry < 6 ) then
    go to 10
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PITCON7_PRB3:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine fxbend(nvar,fpar,ipar,x,fx )

!*****************************************************************************80
!
!! FXBEND evaluates the nonlinear function associated with the two point BVP.
!
!  Equation 1:
!
!    Y(0) = 0
!
!  Equations 2 through NVAR-2:
!
!    Y'' + LAMBDA * EXP ( Y ) = 0
!
!  Equation NVAR-1:
!
!    Y'(1) = 0
!
  integer nvar

  double precision fpar(*)
  double precision fx(nvar)
  double precision hsq
  integer i
  integer ipar(2)
  double precision x(nvar)

  hsq = 1.0 / ( nvar - 2 )**2
!
!  Y(0) = 0
!
  fx(1) = x(1)
!
!  Y'' + LAMBDA * EXP ( Y ) = 0
!
  do i = 2, nvar-2
    fx(i) = x(i-1) - 2.0 * x(i) + x(i+1) + x(nvar) * exp ( x(i) ) * hsq
  end do
!
!  Y'(1) = 0
!
  fx(nvar-1) = x(nvar-1) - x(nvar-2)

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
subroutine fpband ( nvar, fpar, ipar, x, fprime )

!*****************************************************************************80
!
!! FPBAND evaluates the band storage jacobian.
!
  integer nvar

  double precision fpar(*)
  double precision fprime(*)
  double precision hsq
  integer i
  integer indx
  integer ipar(2)
  integer j
  integer ml
  integer mu
  integer nband
  double precision x(nvar)

  hsq=1.0/(nvar-2)**2

  ml=ipar(1)
  mu=ipar(2)
  nband=2*ml+mu+1

  i = 1
  j = 1
  indx=1+j*(nband-1)-ml
  fprime(indx)=1.0

  do i = 2, nvar - 2

    j = i - 1
    indx=i+j*(nband-1)-ml
    fprime(indx)=1.0

    j = i
    indx=i+j*(nband-1)-ml
    fprime(indx)=-2.0+x(nvar)*exp(x(i))*hsq

    j = i + 1
    indx=i+j*(nband-1)-ml
    fprime(indx)=1.0

    j = nvar
    indx=(nvar-1)*nband+i
    fprime(indx)=exp(x(i))*hsq

  end do

  i = nvar - 1
  j = nvar - 2
  indx=i+j*(nband-1)-ml
  fprime(indx)=-1.0

  j = nvar - 1
  indx = i + j*(nband-1)-ml
  fprime(indx) = 1.0

  return
end
