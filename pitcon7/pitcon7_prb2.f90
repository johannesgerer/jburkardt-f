program main

!*****************************************************************************80
!
!! PCPRB2 runs a problem involving the aircraft stability problem.
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
!    Raman Mehra, William Kessel, James Carroll,
!    Global stability and contral analysis of aircraft at high angles of attack,
!    Technical Report CR-215-248-1, -2, -3,
!    Office of Naval Research, June 1977.
!
!    Rami Melhem, Werner Rheinboldt,
!    A Comparison of Methods for Determining Turning Points of Nonlinear Equations,
!    Computing,
!    Volume 29, Number 3, September 1982, pages 201-226.
!
!    Albert Schy, Margery Hannah,
!    Prediction of Jump Phenomena in Roll-coupled Maneuvers of Airplanes,
!    Journal of Aircraft,
!    Volume 14, Number 4, 1977,  pages 375-382.
!
!    John Young, Albert Schy, Katherine Johnson,,
!    Prediction of Jump Phenomena in Aircraft Maneuvers, Including
!    Nonlinear Aerodynamic Effects,
!    Journal of Guidance and Control,
!    Volume 1, Number 1, 1978, pages 26-31.
!
!  The variables X(I) are control parameters for an aircraft under
!  the guidance of a pilot:
!
!    X(1) = Roll rate
!    X(2) = Pitch rate
!    X(3) = Yaw rate
!    X(4) = Incremental angle of attack
!    X(5) = Sideslip angle
!    X(6) = Elevator angle
!    X(7) = Aileron angle
!    X(8) = Rudder angle
!
!  The function:
!
!    For I=1 to 5,
!
!      F(I) = Sum ( 1<=J<=8 ) B(I,J) * X(J)
!           + Sum ( 1<=J<=8 ) Sum(1<=K<=8) PHI(I,J,K) * X(J) * X(K)
!
!    F(6) = X(IFIX1) - VAL1
!    F(7) = X(IFIX2) - VAL2
!
!  Options:
!
!    ICHOOZ  IFIX1  VAL1   IFIX2  VAL2  LIM  BARRAY
!
!       1      6  -0.0500    8    0.0    7   Rheinboldt
!       2      6  -0.0080    8    0.0    7   Rheinboldt
!       3      6   0.0000    8    0.0    7   Rheinboldt
!       4      6   0.0500    8    0.0    7   Rheinboldt
!       5      6   0.1000    8    0.0    7   Rheinboldt
!      11      6  -0.0500    8    0.0    7   Melhem
!      12      6  -0.0080    8    0.0    7   Melhem
!      13      6   0.0000    8    0.0    7   Melhem
!      14      6   0.0500    8    0.0    7   Melhem
!      15      6   0.1000    8    0.0    7   Melhem
!
!
!  Limit points:
!
!
!  Melhem lists the following limit points in X7.  Note that Melhem
!  has BARRAY(4,1)=1.0, BARRAY(4,2)=0.0:
!
!
!     X1       X2        X3        X4        X5      X6     X7      X8
!
!  -2.9691  0.83074  -0.072748  0.41029  -0.26880  -.05  0.50919    0.0
!  -2.8158 -0.17481  -0.089469  0.026319  0.070951 -.008 0.20442    0.0
!  -3.7571 -0.64911  -0.39350   0.091810  0.19685  -.008 -.0038238  0.0
!  -4.1637  0.092284 -0.092610  0.022402 -0.017106 -.008 0.37823    0.0
!  -2.5839 -0.22128  -0.054079  0.013524  0.090871 0.0   0.18608    0.0
!  -3.9007 -1.1421   -0.57863   0.13284   0.32685  0.0  -0.50703    0.0
!  -2.3610 -0.72360   0.032739 -0.039108  0.29347  0.05  0.29272    0.0
!  -2.2982  1.4033    0.063244 -0.079383  0.58336  0.10  0.58336    0.0
!
!
!  Rheinboldt lists the following limit points, using
!  BARRAY(4,1)=0.0, BARRAY(4,2)=1.0:
!
!
!     X1       X2        X3        X4        X5      X6     X7      X8
!
!   2.9648  0.82556   0.07366   0.041309  0.26734 -0.050 -0.050481  0.0
!   2.8173 -0.17628   0.08992   0.026429 -0.07147 -0.008 -0.204973  0.0
!   3.7579 -0.65541   0.38658   0.092520 -0.19867 -0.008  0.006200  0.0
!   4.1638  0.08913   0.09480   0.022888  0.16232 -0.008 -0.377660  0.0
!   2.5873 -0.22354   0.05468   0.013676 -0.09168  0.000 -0.186908  0.0
!   3.9005 -1.14815   0.58156   0.133516 -0.32858  0.000  0.510158  0.0
!   2.3639 -0.72974  -0.31604  -0.038785 -0.29583  0.050 -0.295772  0.0
!   2.2992 -1.41023  -0.06184  -0.079009 -0.58629  0.100 -0.689717  0.0
!
!  Bifurcation points:
!
!    Rheinboldt lists:
!
!    X1      X2        X3        X4        X5       X6         X7      X8
!   4.482   0.1632    0.02373   0.006205  0.03527 -0.0006177 -0.3986  0.0
!   3.319  -0.1869    0.1605    0.04379  -0.06888 -0.01250   -0.2374  0.0
!   4.466   0.1467    0.04045   0.009777  0.03089 -0.006129  -0.3995  0.0
!  -3.325   0.1880   -0.1614    0.04395   0.06911 -0.01247    0.2367  0.0
!
  integer, parameter :: nvar = 8
  integer, parameter :: liw = 29 + nvar
  integer, parameter :: lrw = 29 + nvar * ( nvar + 6 )
!
  external dge_slv
  external fpair
  external fxair
  external pitcon
!
  double precision barray(5,8)
  double precision hmax
  integer i
  integer ichooz
  integer ierror
  integer ifix1
  integer ifix2
  integer ipar(2)
  integer iwork(liw)
  integer j
  integer k
  integer lim
  character ( len = 12 ) name
  double precision rpar(42)
  double precision rwork(lrw)
  double precision val1
  double precision val2
  double precision xr(nvar)
!
!  Set options.
!
  ichooz = 1
  ifix1 = 6
  ifix2 = 8
!
!  Initialize work arrays.
!
  iwork(1:liw) = 0
  rwork(1:lrw) = 0.0D+00
!
!  Set VAL1, VAL2 based on ICHOOZ
!
  if ( ichooz == 1 .or. ichooz == 11 ) then
    val1 = -0.050D+00
  else if ( ichooz == 2 .or. ichooz == 12 ) then
    val1 = -0.008D+00
  else if ( ichooz == 3 .or. ichooz == 13 ) then
    val1 =  0.000D+00
  else if ( ichooz == 4 .or. ichooz == 14 ) then
    val1 =  0.050D+00
  else if ( ichooz == 5 .or. ichooz == 15 ) then
    val1 =  0.100D+00
  else if ( ichooz == 6 .or. ichooz == 16 ) then
    val1 = -0.01250D+00
  end if

  val2 = 0.0D+00
!
!  Set some parameters
!
  ipar(1) = ifix1
  ipar(2) = ifix2

  if (ichooz/=6.and.ichooz/=16)then
    hmax=0.4D+00
  else
    hmax=0.2D+00
  end if

  if(ichooz/=6.and.ichooz/=16)then
    lim=7
  else
    lim=0
  end if
!
!  Set starting point
!
  do i=1,nvar
    xr(i)=0.0D+00
  end do

  xr(ifix1)=val1
  xr(ifix2)=val2
!
!  Set BARRAY
!
  barray(1,1)=-3.933D+00
  barray(2,1)=0.0D+00
  barray(3,1)=0.002D+00
  if(ichooz<10)then
    barray(4,1)=0.0D+00
  else
    barray(4,1)=1.0D+00
  end if
  barray(5,1)=0.0D+00
  barray(1,2)=0.107D+00
  barray(2,2)=-0.987D+00
  barray(3,2)=0.0D+00
  if(ichooz<10)then
    barray(4,2)=1.0D+00
  else
    barray(4,2)=0.0D+00
  end if
  barray(5,2)=0.0D+00
  barray(1,3)=0.126D+00
  barray(2,3)=0.0D+00
  barray(3,3)=-0.235D+00
  barray(4,3)=0.0D+00
  barray(5,3)=-1.0D+00
  barray(1,4)=0.0D+00
  barray(2,4)=-22.95D+00
  barray(3,4)=0.0D+00
  barray(4,4)=-1.0D+00
  barray(5,4)=0.0D+00
  barray(1,5)=-9.99D+00
  barray(2,5)=0.0D+00
  barray(3,5)=5.67D+00
  barray(4,5)=0.0D+00
  barray(5,5)=-0.196D+00
  barray(1,6)=0.0D+00
  barray(2,6)=-28.37D+00
  barray(3,6)=0.0D+00
  barray(4,6)=-0.168D+00
  barray(5,6)=0.0D+00
  barray(1,7)=-45.83D+00
  barray(2,7)=0.0D+00
  barray(3,7)=-0.921D+00
  barray(4,7)=0.0D+00
  barray(5,7)=-0.0071D+00
  barray(1,8)=-7.64D+00
  barray(2,8)=0.0D+00
  barray(3,8)=-6.51D+00
  barray(4,8)=0.0D+00
  barray(5,8)=0.0D+00
!
!  Copy BARRAY, VAL1, VAL2 into RPAR.
!
  k=0
  do i=1,5
    do j=1,nvar
      k=k+1
      rpar(k)=barray(i,j)
    end do
  end do

  k=k+1
  rpar(k)=val1
  k=k+1
  rpar(k)=val2
!
!  Set work array input
!
!  IWORK(1)=0 ; This is a startup
!  IWORK(2)=7 ; Use X(7) for initial parameter
!  IWORK(3)=0 ; Program may choose parameter index
!  IWORK(4)=0 ; Update jacobian every newton step
!  IWORK(5)=0 ; Seek no target values.
!  IWORK(6)=LIM ; Seek limit points in X(LIM)
!  IWORK(7)=1 ; Control amount of output.
!  IWORK(9)=0 ; Jacobian choice.
!
  iwork(1)=0
  iwork(2)=7
  iwork(3)=0
  iwork(4)=0
  iwork(5)=0
  iwork(6)=lim
  iwork(7)=1
  iwork(9)=0
!
!  RWORK(1)=0.0001 ; Absolute error tolerance
!  RWORK(2)=0.0001 ; Relative error tolerance
!  RWORK(3)=0.0001 ; Minimum stepsize
!  RWORK(4)=HMAX   ; Maximum stepsize
!  RWORK(5)=0.1    ; Starting stepsize
!  RWORK(6)=-1.0   ; Starting direction
!  RWORK(7)=0.0    ; Target value
!
  rwork(1)=0.0001D+00
  rwork(2)=0.0001D+00
  rwork(3)=0.0001D+00
  rwork(4)=hmax
  rwork(5)=0.1D+00
  rwork(6)=-1.0D+00
  rwork(7)=0.0D+00

  call timestamp ( )
  write ( *, * ) ' '
  write ( *, * ) 'PITCON7_PRB2:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  PITCON sample program.'
  write ( *, * ) '  The aircraft stability problem'
  write ( *, * ) ' '
  write ( *, '(a,i8)' ) '  Number of equations is ', nvar - 1
  write ( *, '(a,i8)' ) '  Number of variables is ', nvar
  write ( *, * ) ' '
  write ( *, * ) 'using option ',ichooz
  write ( *, * ) 'fix variable ',ifix1,' at ',val1
  write ( *, * ) 'fix variable ',ifix2,' at ',val2
  if(ichooz<10)then
    write ( *, * ) 'rheinboldt version of barray'
  else
    write ( *, * ) 'melhem version of barray'
  end if
  write ( *, * ) ' '
  write ( *, * ) 'step  type of point     '
  write ( *, * ) ' '
  i=0
  name='start point  '
  write(*,'(1x,i3,2x,a12,2x,8f7.3)')i,name,(xr(j),j=1,nvar)

  do i = 1, 50

    call pitcon(fpair,rpar,fxair,ierror,ipar,iwork,liw, &
      nvar,rwork,lrw,xr,dge_slv)

    if (iwork(1)==1) then
      name='corrected    '
    else if (iwork(1)==2) then
      name='continuation '
    else if (iwork(1)==3) then
      name='target point '
    else if (iwork(1)==4) then
      name='limit point  '
    else if (iwork(1)<0) then
      name='jacobian   '
    end if

    write(*,'(1x,i3,2x,a12,2x,8f7.3)')i,name,(xr(j),j=1,nvar)

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PITCON returned an error code:'
      write ( *, * ) 'IERROR = ', ierror
      write ( *, * ) 'The computation is terminated.'
      exit
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'Jacobians:     ', iwork(19)
  write ( *, * ) 'Factorizations:', iwork(20)
  write ( *, * ) 'Solves:        ', iwork(21)
  write ( *, * ) 'Functions:     ', iwork(22)
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PITCON7_PRB2:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine fxair ( nvar, rpar, ipar, x, fx )

!*****************************************************************************80
!
!! FXAIR evaluates the function for the aircraft stability problem.
!
  integer nvar

  double precision fx(nvar)
  integer i
  integer ifix1
  integer ifix2
  integer ipar(2)
  integer j
  integer k
  double precision phi
  double precision rpar(42)
  double precision val1
  double precision val2
  double precision x(nvar)
!
!  Compute linear terms
!
  k=0
  do i=1,5
    fx(i)=0.0D+00
    do j=1,8
      k=k+1
      fx(i)=fx(i)+rpar(k)*x(j)
    end do
  end do
!
!  Compute nonlinear terms
!
  phi=-0.727D+00*x(2)*x(3)+8.39D+00*x(3)*x(4)-684.4D+00*x(4)*x(5)+63.5D+00*x(4)*x(7)
  fx(1)=fx(1)+phi

  phi=0.949D+00*x(1)*x(3)+0.173D+00*x(1)*x(5)
  fx(2)=fx(2)+phi

  phi=-0.716D+00*x(1)*x(2)-1.578D+00*x(1)*x(4)+1.132D+00*x(4)*x(7)
  fx(3)=fx(3)+phi

  phi=-x(1)*x(5)
  fx(4)=fx(4)+phi

  phi=x(1)*x(4)
  fx(5)=fx(5)+phi
!
!  Two function values restrict two variables:
!
  ifix1=ipar(1)
  val1=rpar(41)
  fx(6)=x(ifix1)-val1

  ifix2=ipar(2)
  val2=rpar(42)
  fx(7)=x(ifix2)-val2

  return
end
subroutine fpair(nvar,rpar,ipar,x,fprime)

!*****************************************************************************80
!
!! FPAIR evaluates the jacobian of the aircraft stability function.
!
  integer nvar

  double precision fprime(nvar,nvar)
  integer i
  integer ifix1
  integer ifix2
  integer ipar(2)
  integer j
  integer k
  double precision rpar(42)
  double precision x(nvar)
!
!  The linear part.
!
  k=0
  do i=1,5
    do j=1,8
      k=k+1
      fprime(i,j)=rpar(k)
    end do
  end do

  do i=6,7
    do j=1,8
      fprime(i,j)=0.0D+00
    end do
  end do

  ifix1=ipar(1)
  ifix2=ipar(2)
  fprime(6,ifix1)=1.0D+00
  fprime(7,ifix2)=1.0D+00
!
!  The nonlinear part.
!
  fprime(1,2)=fprime(1,2)-0.727D+00*x(3)
  fprime(1,3)=fprime(1,3)-0.727D+00*x(2)+8.39D+00*x(4)
  fprime(1,4)=fprime(1,4)+8.39D+00*x(3)-684.4D+00*x(5)+63.5D+00*x(7)
  fprime(1,5)=fprime(1,5)-684.4D+00*x(4)
  fprime(1,7)=fprime(1,7)+63.5D+00*x(4)

  fprime(2,1)=fprime(2,1)+0.949D+00*x(3)+0.173D+00*x(5)
  fprime(2,3)=fprime(2,3)+0.949D+00*x(1)
  fprime(2,5)=fprime(2,5)+0.173D+00*x(1)

  fprime(3,1)=fprime(3,1)-0.716D+00*x(2)-1.578D+00*x(4)
  fprime(3,2)=fprime(3,2)-0.716D+00*x(1)
  fprime(3,4)=fprime(3,4)-1.578D+00*x(1)+1.132D+00*x(7)
  fprime(3,7)=fprime(3,7)+1.132D+00*x(4)

  fprime(4,1)=fprime(4,1)-x(5)
  fprime(4,5)=fprime(4,5)-x(1)

  fprime(5,1)=fprime(5,1)+x(4)
  fprime(5,4)=fprime(5,4)+x(1)

  return
end
