      program main

!*******************************************************************************
!
!! SNSIMP is a simple program to call ARPACK for a symmetric eigenproblem.
!
!     This example program is intended to illustrate the
!     simplest case of using ARPACK in considerable detail.
!     This code may be used to understand basic usage of ARPACK
!     and as a template for creating an interface to ARPACK.
!
!     This code shows how to use ARPACK to find a few eigenvalues
!     (lambda) and corresponding eigenvectors (x) for the standard
!     eigenvalue problem:
!
!                        A*x = lambda*x
!
!     where A is a n by n real nonsymmetric matrix.
!
!     The main points illustrated here are
!
!        1) How to declare sufficient memory to find NEV
!           eigenvalues of largest magnitude.  Other options
!           are available.
!
!        2) Illustration of the reverse communication interface
!           needed to utilize the top level ARPACK routine SNAUPD
!           that computes the quantities needed to construct
!           the desired eigenvalues and eigenvectors(if requested).
!
!        3) How to extract the desired eigenvalues and eigenvectors
!           using the ARPACK routine SNEUPD.
!
!     The only thing that must be supplied in order to use this
!     routine on your problem is to change the array dimensions
!     appropriately, to specify WHICH eigenvalues you want to compute
!     and to supply a matrix-vector product
!
!                         w <-  Av
!
!     in place of the call to AV( )  below.
!
!     Once usage of this routine is understood, you may wish to explore
!     the other available options to improve convergence, to solve generalized
!     problems, etc.  Look at the file ex-nonsym.doc in DOCUMENTS directory.
!     This codes implements
!
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!         where A is obtained from the standard central difference
!         discretization of the convection-diffusion operator
!                 (Laplacian u) + rho*(du / dx)
!         on the unit square, with zero Dirichlet boundary condition.
!
!     ... OP = A  and  B = I.
!     ... Assume "call av (nx,x,y)" computes y = A*x
!     ... Use mode 1 of SNAUPD.
!
!\BeginLib
!
!\Routines called:
!     snaupd  ARPACK reverse communication interface routine.
!     sneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     av      Matrix vector multiplication routine that computes A*x.
!     tv      Matrix vector multiplication routine that computes T*x,
!             where T is a tridiagonal matrix.  It is used in routine
!             av.
!
!\Author
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: nsimp.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!---------------------------------------------------------------------------
!
!     %------------------------------------------------------%
!     | Storage Declarations:                                |
!     |                                                      |
!     | The maximum dimensions for all arrays are            |
!     | set here to accommodate a problem size of            |
!     | N <= MAXN                                          |
!     |                                                      |
!     | NEV is the number of eigenvalues requested.          |
!     |     See specifications for ARPACK usage below.       |
!     |                                                      |
!     | NCV is the largest number of basis vectors that will |
!     |     be used in the Implicitly Restarted Arnoldi      |
!     |     Process.  Work per major iteration is            |
!     |     proportional to N*NCV*NCV.                       |
!     |                                                      |
!     | You must set:                                        |
!     |                                                      |
!     | MAXN:   Maximum dimension of the A allowed.          |
!     | MAXNEV: Maximum NEV allowed.                         |
!     | MAXNCV: Maximum NCV allowed.                         |
!     %------------------------------------------------------%
!
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=12, maxncv=30, ldv=maxn)
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Real &
                        ax(maxn), d(maxncv,3), resid(maxn), &
                        v(ldv,maxncv), workd(3*maxn), &
                        workev(3*maxncv), &
                        workl(3*maxncv*maxncv+6*maxncv)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, ierr, &
                        j, ishfts, maxitr, mode1, nconv
      Real &
                        tol, sigmar, sigmai
      logical           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                        zero
      parameter         (zero = 0.0E+0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Real &
                        slapy2, snrm2
      external          slapy2, snrm2, saxpy
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic         abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The following include statement and assignments |
!     | initiate trace output from the internal         |
!     | actions of ARPACK.  See debug.doc in the        |
!     | DOCUMENTS directory for usage.  Initially, the  |
!     | most useful information will be a breakdown of  |
!     | time spent in the various stages of computation |
!     | given by setting mnaupd = 1.                    |
!     %-------------------------------------------------%
!
      include 'debug.h'

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SNSIMP:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A simple ARPACK calling program.'
  write ( *, '(a)' ) '  This program defines an eigenproblem for a'
  write ( *, '(a)' ) '  nonsymmetric matrix.'

      ndigit = -3
      logfil = 6
      mnaitr = 0
      mnapps = 0
      mnaupd = 1
      mnaup2 = 0
      mneigh = 0
      mneupd = 0
!
!     %-------------------------------------------------%
!     | The following sets dimensions for this problem. |
!     %-------------------------------------------------%
!
      nx    = 10
      n     = nx*nx
!
!     %-----------------------------------------------%
!     |                                               |
!     | Specifications for ARPACK usage are set       |
!     | below:                                        |
!     |                                               |
!     |    1) NEV = 4  asks for 4 eigenvalues to be   |
!     |       computed.                               |
!     |                                               |
!     |    2) NCV = 20 sets the length of the Arnoldi |
!     |       factorization.                          |
!     |                                               |
!     |    3) This is a standard problem.             |
!     |         (indicated by bmat  = 'I')            |
!     |                                               |
!     |    4) Ask for the NEV eigenvalues of          |
!     |       largest magnitude.                      |
!     |         (indicated by which = 'LM')           |
!     |       See documentation in SNAUPD for the     |
!     |       other options SM, LR, SR, LI, SI.       |
!     |                                               |
!     | Note: NEV and NCV must satisfy the following  |
!     | conditions:                                   |
!     |              NEV <= MAXNEV                    |
!     |          NEV + 2 <= NCV <= MAXNCV             |
!     |                                               |
!     %-----------------------------------------------%
!
      nev   = 4
      ncv   = 20
      bmat  = 'I'
      which = 'LM'
!
      if ( n > maxn ) then
         print *, ' ERROR with _NSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev > maxnev ) then
         print *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv > maxncv ) then
         print *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
!
!     %-----------------------------------------------------%
!     |                                                     |
!     | Specification of stopping rules and initial         |
!     | conditions before calling SNAUPD                    |
!     |                                                     |
!     | TOL  determines the stopping criterion.             |
!     |                                                     |
!     |      Expect                                         |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!     |               computed   true                       |
!     |                                                     |
!     |      If TOL <= 0,  then TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION parameter         |
!     |      used to specify actions to be taken on return  |
!     |      from SNAUPD. (see usage below)                 |
!     |                                                     |
!     |      It MUST initially be set to 0 before the first |
!     |      call to SNAUPD.                                |
!     |                                                     |
!     | INFO on entry specifies starting vector information |
!     |      and on return indicates error codes            |
!     |                                                     |
!     |      Initially, setting INFO=0 indicates that a     |
!     |      random starting vector is requested to         |
!     |      start the ARNOLDI iteration.  Setting INFO to  |
!     |      a nonzero value on the initial call is used    |
!     |      if you want to specify your own starting       |
!     |      vector (This vector must be placed in RESID).  |
!     |                                                     |
!     | The work array WORKL is used in SNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.                                  |
!     |                                                     |
!     %-----------------------------------------------------%
!
      lworkl  = 3*ncv**2+6*ncv
      tol    = zero
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | Specification of Algorithm Mode:                  |
!     |                                                   |
!     | This program uses the exact shift strategy        |
!     | (indicated by setting IPARAM(1) = 1).             |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of SNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | SNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode1 = 1
!
      iparam(1) = ishfts
!
      iparam(3) = maxitr
!
      iparam(7) = mode1
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine SNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call snaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
                       v, ldv, iparam, ipntr, workd, workl, lworkl, &
                       info )
!
         if (ido == -1 .or. ido == 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- Op*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               |
!           %-------------------------------------------%
!
            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call SNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         endif
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info < 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in SNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ',info
         print *, ' Check the documentation of _naupd'
         print *, ' '
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using SNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        |                                           |
!        | The routine SNEUPD now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1,)                                   |
!        |                                           |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call sneupd ( rvec, 'A', select, d, d(1,2), v, ldv, &
              sigmar, sigmai, workev, bmat, n, which, nev, tol, &
              resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
              lworkl, ierr )
!
!        %------------------------------------------------%
!        | The real parts of the eigenvalues are returned |
!        | in the first column of the two dimensional     |
!        | array D, and the IMAGINARY part are returned   |
!        | in the second column of D.  The corresponding  |
!        | eigenvectors are returned in the first         |
!        | NCONV (= IPARAM(5)) columns of the two         |
!        | dimensional array V if requested.  Otherwise,  |
!        | an orthogonal basis for the invariant subspace |
!        | corresponding to the eigenvalues in D is       |
!        | returned in V.                                 |
!        %------------------------------------------------%
!
         if ( ierr /= 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of SNEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _neupd, info = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '
!
         else
!
            first = .true.
            nconv =  iparam(5)
            do 20 j=1, nconv
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |   ||  A*x - lambda*x ||   |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (IPARAM(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
               if (d(j,2) == zero)  then
!
!  Ritz value is real.
!
                  call av(nx, v(1,j), ax)
                  call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,3) = snrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
!
               else if (first) then
!
!                 %------------------------%
!                 | Ritz value is complex. |
!                 | Residual of one Ritz   |
!                 | value of the conjugate |
!                 | pair is computed.      |
!                 %------------------------%
!
                  call av(nx, v(1,j), ax)
                  call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  call saxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                  d(j,3) = snrm2(n, ax, 1)
                  call av(nx, v(1,j+1), ax)
                  call saxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                  call saxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                  d(j,3) = slapy2( d(j,3), snrm2(n, ax, 1) )
                  d(j,3) = d(j,3) / slapy2(d(j,1),d(j,2))
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if

 20         continue
!
!  Display computed residuals.
!
            call smout(6, nconv, 3, d, maxncv, -6, &
                 'Ritz values (Real, Imag) and residual residuals')
         end if
!
!  Print additional convergence information.
!
         if ( info == 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info == 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit', &
                      ' Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
!
         print *, ' '
         print *, ' _NSIMP '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
                    nconv
         print *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '

      end if

 9000 continue

      end
!
!==========================================================================
!
!     matrix vector subroutine
!
!     The matrix used is the 2 dimensional convection-diffusion
!     operator discretized using central difference.
!
      subroutine av (nx, v, w)

!*******************************************************************************
!
      integer           nx, j, lo
      Real &
                        v(nx*nx), w(nx*nx), one, h2
      parameter         (one = 1.0E+0)
      external          saxpy, tv
!
!     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
!     tridiagonal matrix
!
!                  | T -I          |
!                  |-I  T -I       |
!             OP = |   -I  T       |
!                  |        ...  -I|
!                  |           -I T|
!
!     derived from the standard central difference discretization
!     of the 2 dimensional convection-diffusion operator
!     (Laplacian u) + rho*(du/dx) on a unit square with zero boundary
!     condition.
!
!     When rho*h/2 <= 1, the discrete convection-diffusion operator
!     has real eigenvalues.  When rho*h/2 > 1, it has complex
!     eigenvalues.
!
!     The subroutine TV is called to computed y<---T*x.
!
!
      h2 = one / real((nx+1)*(nx+1))

      call tv(nx,v(1),w(1))
      call saxpy(nx, -one/h2, v(nx+1), 1, w(1), 1)

      do j = 2, nx-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call saxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
         call saxpy(nx, -one/h2, v(lo+nx+1), 1, w(lo+1), 1)
      end do

      lo = (nx-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call saxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)

      return
      end
      subroutine tv (nx, x, y)

!*******************************************************************************
!
      integer           nx, j
      Real &
                        x(nx), y(nx), h, dd, dl, du, h2
!
      Real &
                        one, rho
      parameter         (one = 1.0E+0, rho = 1.0E+2)
!
!     Compute the matrix vector multiplication y<---T*x
!     where T is a nx by nx tridiagonal matrix with DD on the
!     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
!
!     When rho*h/2 <= 1, the discrete convection-diffusion operator
!     has real eigenvalues.  When rho*h/2 > 1, it has complex
!     eigenvalues.
!
      h   = one / real(nx+1)
      h2  = h*h
      dd  = 4.0E+0 / h2
      dl  = -one/h2 - 5.0E-1*rho/h
      du  = -one/h2 + 5.0E-1*rho/h

      y(1) =  dd*x(1) + du*x(2)
      do j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1)
      end do
      y(nx) =  dl*x(nx-1) + dd*x(nx)
      return
      end
