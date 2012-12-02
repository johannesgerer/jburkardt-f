program main

!*****************************************************************************80
!
!! PCPRB7 solves the materially nonlinear rod problem.
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
!    Ivo Babuska, Werner C Rheinboldt,
!    Reliable Error Estimations and Mesh adaptation for the Finite Element Method,
!    In: Computation Methods in Nonlinear Mechanics,
!    Edited by J T Oden,
!    North Holland Publishing Company,
!    Amsterdam 1980, Pages 67-108.
!
!  The general problem:
!
!    -D/DX PHI(X,U,DU/DX,T)  +  PSI(X,U,DU/DX,T) = 0,
!
!    U(0) = U0
!    U(1) = U1
!
!  The particular problem:
!
!    U'' + T*SIN(U+U**2+U**3) = 0
!
!    U(0) = 0.0
!    U(1) = 0.0
!
!  U is approximated by dividing the interval into NINT subintervals,
!  and using the NPOLYS Legendre functions of orders 0 through NPOLYS-1
!  on each subinterval.  The NINT*NPOLYS coefficients of these basis
!  functions form the representation of U.
!
!  Moreover, there are NINT*NPOLYS finite element equations set up,
!  of the form
!
!  INTEGRAL(Subinterval I) -PHI(X,U,U',T)*PL'(J) + PSI(X,U,U',T)*PL(J) = 0
!
!  There is a complication however, in that Lagrange multipliers are
!  used to enforce continuity of the solution U at the interface
!  nodes between neighboring subintervals.  These conditions create
!  new variables, and new equations, and even modify the above standard
!  finite element equations.
!
!  Options:
!
!  ICHOOZ=1  Piecewise linear, 1 continuity condition
!  ICHOOZ=2  Piecewise cubic, 1 continuity condition
!  ICHOOZ=3  Piecewise cubic, 2 continuity conditions
!  ICHOOZ=4  Piecewise quintic, 1 continuity condition
!  ICHOOZ=5  Piecewise quintic, 2 continuity conditions
!  ICHOOZ=6  Piecewise quintic, 3 continuity conditions
!
!  All options use 8 intervals.
!
  integer, parameter :: ngauss = 8
  integer, parameter :: maxpol = 10

  integer, parameter :: mint = 8
  integer, parameter :: mvar = mint * 6 + 2 + ( mint - 1 ) * 3
  integer, parameter :: liw = mvar + 29
  integer, parameter :: lrw = 29 + ( 6 + mvar ) * mvar

  external bval
  external dge_slv
  external fp0011
  external fx0011
  external init
  external pitcon
  external uval

  intrinsic mod

  double precision bcone
  double precision bczero
  double precision dbcodt
  double precision dbczdt
  double precision fpar(1)
  double precision gcoef
  double precision gpoint
  integer i
  integer ichooz
  integer ierror
  integer ihold
  integer ii
  integer ipar(1)
  integer iskip
  integer iwork(liw)
  character ( len = 12 ) name
  integer nbcz
  integer nbco
  integer ndrv
  integer nint
  integer npolys
  integer nvar
  integer nvary
  integer nvarz
  double precision pl
  double precision pld
  double precision rwork(lrw)
  double precision theta
  double precision u
  double precision uprym
  double precision xg
  double precision xl
  double precision xr(mvar)
  double precision xrr
!
  save /intmem/
  save /relmem/
!
  common /intmem/ npolys,ndrv,nvary,nvarz,nbcz,nbco,nint

  common /relmem/ bczero(8),dbczdt(8),bcone(8),dbcodt(8), &
    pl(8),pld(8),gcoef(8),gpoint(8),theta(10,10)
!
!  Zero out work arrays
!
  iwork(1:liw) = 0
  rwork(1:lrw) = 0.0
!
!  Initialize COMMON arrays
!
  call init(gcoef,gpoint,maxpol,ngauss,theta)
!
!  Set number of variables
!
  ichooz=1
  nint=8
  if(ichooz==1)npolys=2
  if(ichooz==2.or.ichooz==3)npolys=4
  if(ichooz==4.or.ichooz==5.or.ichooz==6)npolys=6
  nvary=nint*npolys
  nbcz=1
  if(ichooz==1.or.ichooz==2.or.ichooz==4)ndrv=1
  if(ichooz==3.or.ichooz==5)ndrv=2
  if(ichooz==6)ndrv=3
  nbco=1
  nvarz=nbcz+(nint-1)*ndrv+nbco
  nvar=nvary+nvarz+1
!
!  Starting point
!
  do i=1,nvar
    xr(i)=0.0
  end do
!
!  We perturb the zero solution, by setting the 5-th coefficient
!  to some nonzero value.  We will then also require that the program
!  find an exact solution by correcting this starting point, but
!  without altering the 5-th coefficient, thus guaranteeing that the
!  solution is not all zero.
!
  ihold=5
  xr(ihold)=0.001
!
!  IWORK(1)=0     ; This is a startup
!  IWORK(2)=IHOLD ; Use index IHOLD for first parameter
!  IWORK(3)=0     ; Program may choose index
!  IWORK(4)=0     ; Update jacobian every newton step
!  IWORK(5)=0     ; Seek no target values.
!  IWORK(6)=NVAR  ; Seek limit points in index NVAR
!  IWORK(7)=1     ; Moderate amount of output
!  IWORK(9)=0     ; Use user's jacobian routine
!
!  IWORK(17)=20   ; (Unusual setting)  Increase the allowed number of Newton
!                   iterations to 20.
!
  iwork(1)=0
  iwork(2)=ihold
  iwork(3)=0
  iwork(4)=0
  iwork(5)=0
  iwork(6)=nvar
  iwork(7)=2
  iwork(9)=0
  iwork(17)=20
!
!  RWORK(1)=0.0001 ; Absolute error tolerance
!  RWORK(2)=0.0001 ; Relative error tolerance
!  RWORK(3)=0.5    ; Minimum stepsize
!  RWORK(4)=0.5    ; Maximum stepsize
!  RWORK(5)=0.25   ; Starting stepsize
!  RWORK(6)=1.0    ; Starting direction
!  RWORK(7)=0.0    ; Target value
!
  rwork(1)=0.0001
  rwork(2)=0.0001
  rwork(3)=0.001
  rwork(4)=0.5
  rwork(5)=0.25
  rwork(6)=1.0
  rwork(7)=0.0

  call timestamp ( )
  write ( *, * ) ' '
  write ( *, * ) 'PITCON7_PRB7:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  PITCON sample program'
  write ( *, * ) ' '
  write ( *, '(a,i8)' ) '  Number of equations is ', nvar - 1
  write ( *, '(a,i8)' ) '  Number of variables is ', nvar
  write ( *, * ) ' '
  write ( *, * ) '  The materially nonlinear problem.'
  write ( *, * ) ' '

  if(ichooz==1)then
    write ( *, * ) 'Piecewise linears'
  elseif(ichooz==2.or.ichooz==3)then
    write ( *, * ) 'Piecewise cubics'
  else
    write ( *, * ) 'Piecewise quintics.'
  end if

  if(ichooz==6)then
    write ( *, * ) '3 continuity conditions at breakpoints'
  elseif(ichooz==5.or.ichooz==3)then
    write ( *, * ) '2 continuity conditions at breakpoints'
  else
    write ( *, * ) '1 continuity condition at breakpoints.'
  end if

  write ( *, * ) ' '
  write ( *, * ) 'Step  Type of point     Lambda'
  write ( *, * ) ' '
  i=0
  name='Start point  '
  write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,xr(nvar)
!
!  Take another step.
!
  do i=1,30

    call pitcon(fp0011,fpar,fx0011,ierror,ipar,iwork,liw, &
      nvar,rwork,lrw,xr,dge_slv)

    if(iwork(1)==1)then
      name='corrected    '
    elseif(iwork(1)==2)then
      name='continuation '
    elseif(iwork(1)==3)then
      name='target point '
    elseif(iwork(1)==4)then
      name='limit point  '
    elseif(iwork(1)<0)then
      name='jacobian   '
    end if

    write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,xr(nvar)
!
!  Print out shape of rod every 5 steps
!
    if(mod(i,5)==0)then
      write ( *, * ) ' '
      write ( *, * ) 'Current rod coordinates:'
      write ( *, * ) ' '
      write ( *, * ) 'X,  U(X),  U''(X)'
      write ( *, * ) ' '

      do ii=1,nint

        xl=(ii-1)/real(nint)
        xrr=(ii)/real(nint)
        iskip=(ii-1)*npolys
        xg=0.5*(xl+xrr)
        call bval(npolys,pl,pld,xg,xl,xrr)
        call uval(npolys,pl,pld,xg,u,uprym,xr,iskip)
        write(*,'(1x,3g14.6)')xg,u,uprym
      end do

      write ( *, * ) ' '

    end if

    if(iwork(1)==3)then
      exit
    end if

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PITCON returned an error code:'
      write ( *, * ) 'IERROR = ', ierror
      write ( *, * ) ' '
      write ( *, * ) 'The computation failed.'
      exit
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'Jacobians      ',iwork(19)
  write ( *, * ) 'Factorizations ',iwork(20)
  write ( *, * ) 'Solves         ',iwork(21)
  write ( *, * ) 'Functions      ',iwork(22)
!
! Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PITCON7_PRB7:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine init(gcoef,gpoint,maxpol,ngauss,theta)
!
!*****************************************************************************80
!
!! INIT initializes the values of the COMMON block variables
!  GCOEF, GPOINT, and THETA, which are the 8 point Gauss integration
!  coefficients and abscissas, and the value and derivatives of the Legendre
!  polynomials at 1.
!
  integer maxpol
  integer ngauss
!
  double precision gcoef(ngauss)
  double precision gpoint(ngauss)
  integer i
  integer j
  double precision theta(maxpol,maxpol)
!
  gcoef(1)=.1012285363
  gcoef(2)=.2223810345
  gcoef(3)=.3137066459
  gcoef(4)=.3626837834
  gcoef(5)=.3626837834
  gcoef(6)=.3137066459
  gcoef(7)=.2223810345
  gcoef(8)=.1012285363
!
  gpoint(1)= .9602898565
  gpoint(2)= .7966664774
  gpoint(3)= .5255324099
  gpoint(4)= .1834346425
  gpoint(5)=-.1834346425
  gpoint(6)=-.5255324099
  gpoint(7)=-.7966664774
  gpoint(8)=-.9602898565
!
!  THETA(I,J) is derivative (I-1) of Legendre polynomial (J-1) at X=1.
!
  do j=1,maxpol

    theta(1,j)=1.0
    do i=2,maxpol
      theta(i,j)=(j+i-2)*(j-i+1)*theta(i-1,j)/(2*i-2)
    end do

  end do

  return
end
subroutine fx0011(nvar,fpar,ipar,x,fx )
!
!*****************************************************************************80
!
!! FX0011 evaluates the nonlinear function which constrains the shape
!  of the rod.
!
  external bounds
  external bval
  external fval
  external uval
!
  integer nvar
!
  double precision bcone
  double precision bczero
  double precision dbcodt
  double precision dbczdt
  double precision dtdx
  double precision dtdxl
  double precision dtdxr
  double precision fpar(*)
  double precision fx(nvar)
  double precision gcoef
  double precision gpoint
  double precision h2i
  double precision h2il
  double precision h2ir
  integer i
  integer ieqn
  integer indx
  integer ipar(*)
  integer iskip
  integer ivar
  integer j
  integer k
  integer khi
  integer l
  integer lhil
  integer lhir
  integer nbcz
  integer nbco
  integer ncl
  integer ncr
  integer ndrv
  integer ndsum
  integer neqn
  integer nint
  integer nodes
  integer npolys
  integer npsum
  integer nvary
  integer nvarz
  double precision phi
  double precision phipt
  double precision phipu
  double precision phipup
  double precision pl
  double precision pld
  double precision psi
  double precision psipt
  double precision psipu
  double precision psipup
  double precision s
  double precision tcon
  double precision term
  double precision terml
  double precision termr
  double precision theta
  double precision u
  double precision uprym
  double precision x(nvar)
  double precision xc
  double precision xg
  double precision xl
  double precision xr
!
  save /intmem/
  save /relmem/
!
  common /intmem/ npolys,ndrv,nvary,nvarz,nbcz,nbco,nint

  common /relmem/ bczero(8),dbczdt(8),bcone(8),dbcodt(8), &
    pl(8),pld(8),gcoef(8),gpoint(8),theta(10,10)
!
!  Zero out FX
!
  tcon=x(nvar)
  call bounds
  neqn=nvar-1

  fx(1:neqn)=0.0
!
!  1. Set up the terms (A*Y) involving the bivariate form.
!
!  For each interval I:
!
  do 40 i=1,nint

    xl=(i-1)/real(nint)
    xr=(i)/real(nint)
    iskip=(i-1)*npolys
!
!  For each Gauss point J, evaluate the integrand.
!
    do 30 j=1,8

      xg=xl+(gpoint(j)+1.0)*(xr-xl)/2.0
      call bval(npolys,pl,pld,xg,xl,xr)
      call uval(npolys,pl,pld,xg,u,uprym,x,iskip)
      call fval(xg,u,uprym,tcon,phi,phipu,phipup,phipt,psi,psipu,psipup,psipt)
!
!  L   PROJECT ONTO EACH TEST FUNCTION PL(L) AND PLD(L)
!
      ieqn=iskip

      do l=1,npolys
        ieqn=ieqn+1
        term=gcoef(j)*(xr-xl)*(psi*pl(l)+phi*pld(l))/2.0
        fx(ieqn)=fx(ieqn)+term
      end do

   30     continue
   40   continue
!
!  2. ADD THE TERMS INVOLVING THE CONTINUITY OF THE TEST FUNCTIONS
!  WHICH ARE THE TERMS B*Z IN  F=A*Y + B*Z
!
  ncl=nvary
  ncr=nvary+nbcz
!
!  I  FOR EACH INTERVAL
!
  do 110 i=1,nint

    xl=(i-1)/real(nint)
    xr=(i)/real(nint)
    dtdx=2.0/(xr-xl)
!
!  J   FOR THE POLYNOMIALS USED IN APPROXIMATING EACH U,
!  COUNT CONDITIONS AT LEFT ENDPOINT, LHIL, AND AT RIGHT, LHIR
!  IF WE ARE IN THE FIRST OR LAST INTERVAL, ONE OF
!  THESE WILL BE BOUNDARY CONDITIONS
!
    lhil=ndrv
    lhir=ndrv
    if (i==1) lhil=nbcz
    if (i==nint) lhir=nbco
    if (lhil<=0.and.lhir<=0) go to 100
!
!  K  FOR EACH TEST FUNCTION PL(K)
!
    do 90 k=1,npolys
      terml=0.0
      termr=0.0
      s=(-1.0)**(k+1)
      ieqn=(i-1)*npolys+k
!
!  L  FOR DERIVATIVES 0 THRU NBCZ SPECIFIED AT 0.0 TO BE ZERO
!  OR DERIVATIVES 0 THRU NDRV SPECIFIED TO BE CONTINUOUS
!  AT INTERIOR NODES
!  OR DERIVATIVES 0 THRU NBCO SPECIFIED AT 1.0 TO BE ZERO
!
      h2i=1.0
      do l=1,lhil
        s=-s
        indx=ncl+l
        terml=terml+s*h2i*theta(l,k)*x(indx)
        h2i=h2i*dtdx
      end do

      h2i=1.0
      do l=1,lhir
        indx=ncr+l
        termr=termr+h2i*theta(l,k)*x(indx)
        h2i=h2i*dtdx
      end do

      fx(ieqn)=fx(ieqn)+terml+termr
   90     continue
    if (lhil>0) ncl=ncl+lhil
    if (lhir>0) ncr=ncr+lhir
  100   continue
  110   continue
!
!  3. CREATE THE TERMS FOR THE U FUNCTIONS AND THEIR DERIVATIVES
!  THE MATRIX TERMS  ( C*Y )
!  ONE EQUATION IS GENERATED FOR COMPONENT AND CONDITION
!
  npsum=0
  dtdxr=0.0
  dtdxl=0.0
!
!  I  FOR EACH NODE
!
  ndsum=nvary
  nodes=nint+1
  do 210 i=1,nodes

    if (i>1) xl=(i-2)/real(nint)
    xc=(i-1)/real(nint)
    if (i<nodes) xr=(i)/real(nint)
    if (xc/=xl) dtdxl=2.0/(xc-xl)
    if (xr/=xc) dtdxr=2.0/(xr-xc)
    h2il=1.0
    h2ir=1.0
!
!  K   FOR DERIVATIVES 0 THRU NBCZ SPECIFIED AT 0.0
!  OR DERIVATIVES 0 THRU NDRV SPECIFIED TO BE CONTINUOUS
!  AT INTERIOR NODES
!  OR DERIVATIVES 0 THRU NBCO SPECIFIED AT 1.0
!
    khi=ndrv
    if (i==1) khi=nbcz
    if (i==nodes) khi=nbco
    if (khi<=0) go to 190
    do 180 k=1,khi
      s=(-1.0)**(k+1)
!
!  L   SET UP THE TERM FROM THE LEFT HAND INTERVAL
!
      terml=0.0
      if (i>1) go to 120
      terml=bczero(k)
      go to 140
  120 continue
      do l=1,npolys
        ivar=npsum+l-npolys
        terml=terml+x(ivar)*h2il*theta(k,l)
      end do
!
!  L   SET UP THE TERM FROM THE RIGHT HAND INTERVAL
!
  140     if (i<nodes) go to 150
      termr=-bcone(k)
      go to 170
  150     termr=0.0
      do l=1,npolys
        ivar=npsum+l
        s=-s
        termr=termr+s*x(ivar)*h2ir*theta(k,l)
      end do
  170     ieqn=ndsum+k
      fx(ieqn)=terml+termr
      h2il=h2il*dtdxl
      h2ir=h2ir*dtdxr
  180     continue
    ndsum=ndsum+khi
  190   npsum=npsum+npolys
  200   continue
  210   continue
  return
end
subroutine fp0011(nvar,fpar,ipar,x,fprime )
!
!*****************************************************************************80
!
!! FP0011 computes the jacobian of the materially nonlinear rod problem.
!
  external bounds
  external bval
  external fval
  external uval
!
  integer nvar
!
  double precision bcone
  double precision bczero
  double precision dbcodt
  double precision dbczdt
  double precision dtdx
  double precision fpar(*)
  double precision fprime(nvar,nvar)
  double precision gcoef
  double precision gpoint
  double precision h2i
  integer i
  integer ieqn
  integer ipar(*)
  integer iskip
  integer ivar
  integer j
  integer k
  integer khil
  integer khir
  integer l
  integer lhil
  integer lhir
  integer n
  integer nbcz
  integer nbco
  integer ncl
  integer ncr
  integer ndrv
  integer nint
  integer npolys
  integer npsum
  integer nvary
  integer nvarz
  double precision phi
  double precision phipt
  double precision phipu
  double precision phipup
  double precision pl
  double precision pld
  double precision psi
  double precision psipt
  double precision psipu
  double precision psipup
  double precision s
  double precision sum
  double precision tcon
  double precision term
  double precision term1
  double precision term2
  double precision term3
  double precision term4
  double precision theta
  double precision u
  double precision uprym
  double precision x(nvar)
  double precision xg
  double precision xl
  double precision xr
!
  save /intmem/
  save /relmem/
!
  common /intmem/ npolys,ndrv,nvary,nvarz,nbcz,nbco,nint

  common /relmem/ bczero(8),dbczdt(8),bcone(8),dbcodt(8), &
    pl(8),pld(8),gcoef(8),gpoint(8),theta(10,10)
!
!  ZERO OUT THE MATRIX
!
  tcon=x(nvar)
  call bounds

  fprime(1:nvar,1:nvar)=0.0
!
!  1.  SET UP THE TERMS FROM THE BIVARIATE FORM (A*Y )
!
!      I  FOR EACH INTERVAL (XL,XR)
!
  do i=1,nint

    xl=(i-1)/real(nint)
    xr=(i)/real(nint)
    iskip=(i-1)*npolys
!
!  J  FOR EACH GAUSS POINT XG IN THE INTERVAL
!  NOTE THAT THE LEGENDRE POLYNOMIALS AND DERIVATIVES ARE SET UP
!  BY THE CALL TO VALUE
!
    do j=1,8
      xg=xl+0.5*(xr-xl)*(gpoint(j)+1.0)
      call bval(npolys,pl,pld,xg,xl,xr)
      call uval(npolys,pl,pld,xg,u,uprym,x,iskip)
      call fval(xg,u,uprym,tcon,phi,phipu,phipup,phipt,psi,psipu,psipup,psipt)
!
!  L  FOR EACH LEGENDRE POLYNOMIAL COEFFICIENT
!
      ieqn=iskip
      do l=1,npolys
        ieqn=ieqn+1
        term=gcoef(j)*0.5*(xr-xl)*(psipt*pl(l)+phipt*pld(l))
        fprime(ieqn,nvar)=fprime(ieqn,nvar)+term
!
!  N   FOR EACH Y-COEFFICIENT OF A U
!
        do n=1,npolys
          ivar=npolys*(i-1)+n
          term1=psipu*pl(n)*pl(l)
          term2=psipup*pld(n)*pl(l)
          term3=phipu*pl(n)*pld(l)
          term4=phipup*pld(n)*pld(l)
          sum=term1+term2+term3+term4
          fprime(ieqn,ivar)=fprime(ieqn,ivar)+gcoef(j)*0.5*(xr-xl)*sum
        end do
      end do
    end do
  end do
!
!  2. ADD THE TERMS INVOLVING THE CONTINUITY OF THE TEST FUNCTIONS
!  WHICH ARE THE TERMS B*Z IN  F=A*Y+B*Z
!
!
!     I  FOR EACH INTERVAL,
!
  do i=1,nint
    ncl=nvary
    if (i>1) ncl=nvary+nbcz+(i-2)*ndrv
    ncr=nvary+nbcz+(i-1)*ndrv

    xl = real ( i - 1 ) / real ( nint )
    xr = real ( i ) / real ( nint )
    dtdx=2.0/(xr-xl)
!
!  J   FOR THE POLYNOMIALS USED IN APPROXIMATING EACH U,
!      COUNT CONDITIONS AT LEFT ENDPOINT, LHIL, AND AT RIGHT, LHIR
!      IF WE ARE IN THE FIRST OR LAST INTERVAL, ONE OF
!      THESE WILL BE BOUNDARY CONDITIONS
!
      lhil=ndrv
      lhir=ndrv
      if (i==1) lhil=nbcz
      if (i==nint) lhir=nbco
      if (lhil<=0.and.lhir<=0) go to 130
!
!  K  FOR EACH TEST FUNCTION PL(K)
!
      do k=1,npolys
        s=(-1.0)**(k+1)
        ieqn=(i-1)*npolys+k
!
!  L  FOR DERIVATIVES 0 THRU NBCZ SPECIFIED AT 0.0 TO BE ZERO
!     OR DERIVATIVES 0 THRU NDRV SPECIFIED TO BE CONTINUOUS
!     AT INTERIOR NODES
!     OR DERIVATIVES 0 THRU NBCO SPECIFIED AT 1.0 TO BE ZERO
!
!     EVALUATE CONTRIBUTION FROM LEFT ENDPOINT
!
        h2i=1.0
        do l=1,lhil
          s=-s
          ivar=ncl+l
          fprime(ieqn,ivar)=s*h2i*theta(l,k)
          h2i=h2i*dtdx
        end do
!
!  EVALUATE CONTRIBUTION FROM RIGHT ENDPOINT
!
        h2i=1.0
        do l=1,lhir
          ivar=ncr+l
          fprime(ieqn,ivar)=h2i*theta(l,k)
          h2i=h2i*dtdx
        end do
      end do
      if (lhil>0) ncl=ncl+lhil
      if (lhir>0) ncr=ncr+lhir
  130     continue
  end do
!
!  3. CREATE THE TERMS FOR THE U FUNCTIONS AND THEIR DERIVATIVES
!  THE MATRIX TERMS  ( C*Y )
!  ONE EQUATION IS GENERATED FOR COMPONENT AND CONDITION
!
!
!  I  FOR EACH INTERVAL
!
  do i=1,nint
    ncl=nvary
    if (i>1) ncl=nvary+nbcz+(i-2)*ndrv
    ncr=nvary+nbcz+(i-1)*ndrv
    npsum=(i-1)*npolys

    xl=(i-1)/real(nint)
    xr=(i)/real(nint)
    dtdx=2.0/(xr-xl)
    h2i=1.0
!
!  K   FOR DERIVATIVES 0 THRU NBCZ SPECIFIED AT 0.0
!      OR DERIVATIVES 0 THRU NDRV SPECIFIED TO BE CONTINUOUS
!      AT INTERIOR NODES
!      OR DERIVATIVES 0 THRU NBCO SPECIFIED AT 1.0
!
    khil=ndrv
    if (i==1) khil=nbcz
    do k=1,khil
      ieqn=ncl+k
!
!  L   SET UP THE TERM FROM THE LEFT HAND ENDPOINT
!
      if (i==1) fprime(ieqn,nvar)=dbczdt(k)
      s=(-1.0)**(k+1)
      do l=1,npolys
        ivar=npsum+l
        s=-s
        fprime(ieqn,ivar)=s*h2i*theta(k,l)
      end do
      h2i=h2i*dtdx
    end do
    ncl=ncl+khil
    h2i=1.0
    khir=ndrv
    if (i==nint) khir=nbco
    do k=1,khir
      ieqn=ncr+k
      fprime(ieqn,nvar)=0.0
      if (i==nint) fprime(ieqn,nvar)=-dbcodt(k)
      do l=1,npolys
        ivar=npsum+l
        fprime(ieqn,ivar)=h2i*theta(k,l)
      end do
      h2i=h2i*dtdx
    end do
    ncr=ncr+khir
    npsum=npsum+npolys
  end do

  return
end
subroutine bval(npolys,pl,pld,xg,xl,xr)
!
!*****************************************************************************80
!
!  EVALUATE LEGENDRE POLYNOMIALS AND DERIVATIVES
!
  double precision a
  integer i
  integer npolys
  double precision pl(8)
  double precision pld(8)
  double precision t
  double precision xg
  double precision xl
  double precision xr
!
  t=-1.0+2.0*(xg-xl)/(xr-xl)
  pl(1)=1.0
  pld(1)=0.0
  pl(2)=t
  pld(2)=1.0
  a=0.0
  do i=3,8
    a=a+1.0
    pl(i)=((a+a+1.0)*t*pl(i-1)-a*pl(i-2))/(a+1.0)
    pld(i)=((a+a+1.0)*(t*pld(i-1)+pl(i-1))-a*pld(i-2))/(a+1.0)
  end do

  do i=1,8
    pld(i)=2.0*pld(i)/(xr-xl)
  end do
  return
end
subroutine uval(npolys,pl,pld,xg,u,uprym,x,iskip)
!
!*****************************************************************************80
!
!  UVAL evaluates the solution U and its spatial derivative UPRYM
!  by summing the coefficients times the polynomial values.
!
  integer npolys
!
  integer indx
  integer iskip
  integer j
  double precision pl(npolys)
  double precision pld(npolys)
  double precision u
  double precision uprym
  double precision x(*)
  double precision xg
!
  u=0.0
  uprym=0.0
  do j=1,npolys
    indx=iskip+j
    u=u+x(indx)*pl(j)
    uprym=uprym+x(indx)*pld(j)
  end do

  return
end
subroutine fval(xg,u,uprym,tcon,phi,phipu,phipup,phipt,psi,psipu,psipup,psipt)

!*****************************************************************************80
!
!! FVAL: AUXILLIARY ROUTINE FOR MATERIALLY NONLINEAR PROBLEM.
!
  double precision phi
  double precision phipu
  double precision phipup
  double precision phipt
  double precision psi
  double precision psipu
  double precision psipup
  double precision psipt
  double precision tcon
  double precision u
  double precision uprym
  double precision xg

  phi=-uprym
  phipu=0.0
  phipup=-1.0
  phipt=0.0

  psi=tcon*sin(u*(1.0+u*(1.0+u)))
  psipu=tcon*(1.0+u*(2.0+3.0*u))*cos(u*(1.0+u*(1.0+u)))
  psipup=0.0
  psipt=sin(u*(1.0+u*(1.0+u)))

  return
end
subroutine bounds

!*****************************************************************************80
!
!! BOUNDS sets boundary values for the materially nonlinear problem.
!
  double precision bcone
  double precision bczero
  double precision dbcodt
  double precision dbczdt
  double precision gcoef
  double precision gpoint
  double precision pl
  double precision pld
  double precision theta

  save /relmem/

  common /relmem/ bczero(8),dbczdt(8),bcone(8),dbcodt(8), &
    pl(8),pld(8),gcoef(8),gpoint(8),theta(10,10)

  bcone(1)=0.0
  bczero(1)=0.0
  dbcodt(1)=0.0
  dbczdt(1)=0.0

  return
end
