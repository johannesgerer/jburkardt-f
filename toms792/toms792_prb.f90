program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS792_PRB.
!
!                          CS2TST
!                         11/20/98
!
!   This program computes interpolation errors using the
! scattered data package CSHEP2D for each of ten test
! functions and a 33 by 33 uniform grid of interpolation
! points in the unit square.
!
!   This program uses Subroutines TESTDT and TSTFN1 from
! ACM Algorithm SURVEY to generate a node set and and the
! test function values.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) NFUN, the number of test functions.
!
!    Local, integer ( kind = 4 ) NSET, the number of node sets.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 100
  integer ( kind = 4 ), parameter :: nrmax = 10
  integer ( kind = 4 ), parameter :: ni = 33

  real ( kind = 8 ) a(9,nmax)
  real ( kind = 8 ) ft(ni,ni)
  integer ( kind = 4 ) lcell(nrmax,nrmax)
  integer ( kind = 4 ) lnext(nmax)
  integer ( kind = 4 ), parameter :: nfun = 10
  integer ( kind = 4 ), parameter :: nset = 5
  real ( kind = 8 ) p(ni)
  real ( kind = 8 ) rw(nmax)
  real ( kind = 8 ) w(nmax)
  real ( kind = 8 ) x(nmax)
  real ( kind = 8 ) y(nmax)

      double precision dum, dx, dy, ermax, ermean, pw, &
                       rmax, ssa, sse, ssm, sum, xmin, ymin
      double precision cs2val
      integer          i, ier, j, k, kf, kff, kfl, ks, &
                       n, nc, nfun, np, nr, nset, nw, nwmax
!
! Input format:
!
  100 FORMAT (I2)
!
! Get a user-specified node set number KS.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS792_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS792 library.'
  write ( *, '(a)' ) ' '

    1 WRITE (*,110) NSET
  110 FORMAT (///13X,'CS2TST:  CSHEP2D Test Program'// &
              5X,'Specify a node set number in the range 1', &
                 ' to ',I2,':'/)
      READ (*,100,END=999,ERR=1) KS
      IF (KS .LT. 1  .OR.  KS .GT. NSET) GO TO 1
!
! Copy N and the nodal coordinates for node set KS.
!
      CALL TESTDT ( KS, N, X, Y )

      IF (N .LT. 10  .OR.  N .GT. NMAX) then
        WRITE (*,500) N, NMAX
        STOP
      end if
!
! Allow the user to specify a range of function numbers.
!
      WRITE (*,120) NFUN
  120 FORMAT (//5X,'Specify the first test function ', &
                   '(1 to ',I2,'):'/)
      READ (*,100,ERR=1) KFF
      IF (KFF .LT. 1  .OR.  KFF .GT. NFUN) GO TO 1
      write ( *, * ) '  KFF = ', kff

      WRITE (*,130) KFF, NFUN
  130 FORMAT (//5X,'Specify the last test function (', &
                   I2,' to ',I2,'):'/)
      READ (*,100,ERR=1) KFL
      IF (KFL .LT. KFF  .OR.  KFL .GT. NFUN) GO TO 1
      write ( *, * ) '  KFL = ', kfl

      NFUN = KFL-KFF+1
!
! Input NC, NW, and NR from the console.
!
      NWMAX = MIN(40,N-1)
    2 WRITE (*,140) N
  140 FORMAT (//5X,'N =',I4//5X, &
              'Specify the number of nodes NC for the ', &
              'least squares fits.'/5X,'NC = 17 is ', &
              'recommended.  NC GE 9.'/)
      READ (*,100,ERR=2) NC
      IF (NC .LT. 9  .OR.  NC .GT. NWMAX) GO TO 2
      write ( *, * ) '  NC = ', nc

    3 WRITE (*,150)
  150 FORMAT (///5X,'Specify the number of nodes NW for ', &
              'the weights.  NW = 30 is'/5X,'recommended. ', &
              ' 1  LE  NW  LE  MIN(40,N-1).')
      READ (*,100,ERR=2) NW
      IF (1 .GT. NW  .OR.  NW .GT. NWMAX) GO TO 3
      write ( *, * ) '  NW = ', nw
!
    4 WRITE (*,160) NRMAX
  160 FORMAT (///5X,'Specify the number of rows and column', &
              's NR in the uniform grid'/5X,'of cells used', &
              ' to locate nearest neighbors.  NR = Sqrt(N/', &
              '3) is'/5X,'recommended.  1 LE NR LE ',I2)
      READ (*,100,ERR=3) NR
      IF (NR .LT. 1  .OR.  NR .GT. NRMAX) GO TO 4
      write ( *, * ) '  NR = ', nr
!
! Set up uniform grid points.
!
      DO I = 1,NI
        P(I) = DBLE ( I - 1 ) / DBLE ( NI - 1 )
      end do
!
! Initialize the average SSE/SSM value to zero.
!
      SSA = 0.0D+00
!
! Print a heading and loop on test functions.
!
      WRITE (*,200) KS, N, NI, NC, NW, NR

      DO KF = KFF, KFL
!
!   Compute true function values at the nodes.
!
        DO K = 1, N
          CALL TSTFN1 (KF,X(K),Y(K),0, W(K),DUM,DUM)
        end do
!
!   Compute true function values FT on the uniform grid, and
!     accumulate the sum of values SUM and sum of squared
!     values SSM.
!
        sum = 0.0D+00
        ssm = 0.0D+00
        do i = 1,ni
          do j = 1,ni
            call tstfn1 (kf,p(i),p(j),0, ft(i,j),dum,dum)
            sum = sum + ft(i,j)
            ssm = ssm + ft(i,j)**2
          end do
        end do
!
!   Compute the sum of squared deviations from the mean SSM.
!
        SSM = SSM - SUM * SUM / DBLE ( NI * NI )
!
!   Compute parameters A and RW defining the interpolant.
!
        CALL CSHEP2 (N,X,Y,W,NC,NW,NR, LCELL,LNEXT,XMIN, &
                     YMIN,DX,DY,RMAX,RW,A,IER)

        IF ( IER .NE. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TOMS792_PRB - Fatal error!'
          write ( *, '(a,i4)' ) '  CSHEP2 returns IER = ', ier
          IF ( IER == 2 ) then
            write ( *, '(a)' ) '  Duplicate nodes encountered.'
          end if
          if ( IER == 3 ) then
            write ( *, '(a)' ) '  All nodes are collinear.'
          end if
          stop
        end if
!
!   Compute interpolation errors.
!
        ermean = 0.0D+00
        ermax = 0.0D+00
        sse = 0.0D+00

        do i = 1,ni
          do j = 1,ni
            pw = cs2val (p(i),p(j),n,x,y,w,nr,lcell,lnext, &
                         xmin,ymin,dx,dy,rmax,rw,a) - &
                 ft(i,j)
            ermean = ermean + abs(pw)
            ermax = max(ermax,abs(pw))
            sse = sse + pw*pw
          end do
        end do

        np = ni*ni
        ermean = ermean/dble(np)
        sse = sse/ssm
        ssa = ssa + sse
        write (*,210) kf, ermax, ermean, sse

      end do
!
! Print the average SSE/SSM value (averaged over the test
!   functions).
!
      WRITE (*,220) SSA/DBLE(NFUN)
      go to 1
!
!  Terminate.
!
999   continue
      write ( *, * ) ' '
      write ( *, '(a)' ) 'TOMS792_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      stop
!
! Print formats:
!
  200 FORMAT (30X,'CS2TST Output:'// &
              1X,16X,'CSHEP2D Interpolation Errors for ', &
              'N nodes and an'/ &
              1X,6X,'NI by NI Uniform Grid of Interpolation', &
              'ion Points in the Unit Square'//1X, &
              6X,'Node set ',I2,4X,'N =',I4,4X,'NI = ',I2, &
              4X,'NC = ',I2,4X,'NW = ',I2,4X,'NR = ',I2/// &
              1X,16X,'Function',4X,'Max Error',4X, &
              'Mean Error',4X,'SSE/SSM'/)
  210 FORMAT (1X,19X,I2,9X,F7.4,6X,F8.5,2X,F9.6)
  220 FORMAT (//1X,11X,'Average SSE/SSM over the test ', &
              'functions = ',F9.6)
!
! Error message formats:
!
  500 FORMAT (///1X,10X,'*** Error in data -- N = ',I4, &
              ', Maximum value =',I4,' ***')
  510 FORMAT (///1X,14X,'*** Error in CSHEP2 -- duplicate ', &
              'nodes encountered ***')
  520 FORMAT (///1X,14X,'*** Error in CSHEP2 -- all nodes ', &
              'are collinear ***')
      END
