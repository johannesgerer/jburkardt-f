program main

!*******************************************************************************
!
!! SSSIMP is a simple program to call ARPACK for a symmetric eigenproblem.
!
!  Discussion:
!
!    This example program is intended to illustrate the
!    simplest case of using ARPACK in considerable detail.
!    This code may be used to understand basic usage of ARPACK
!    and as a template for creating an interface to ARPACK.
!
!    This code shows how to use ARPACK to find a few eigenvalues
!    LAMBDA and corresponding eigenvectors X for the standard
!    eigenvalue problem:
!
!      A * X = LAMBDA * X
!
!    where A is an N by N real symmetric matrix.
!
!    The main points illustrated here are:
!
!    1) How to declare sufficient memory to find NEV
!       eigenvalues of largest magnitude.  Other options
!       are available.
!
!    2) Illustration of the reverse communication interface
!       needed to utilize the top level ARPACK routine SSAUPD
!       that computes the quantities needed to construct
!       the desired eigenvalues and eigenvectors (if requested).
!
!    3) How to extract the desired eigenvalues and eigenvectors
!       using the ARPACK routine SSEUPD.
!
!    The only things that must be supplied in order to use this
!    routine on your problem is:
!
!    * to change the array dimensions appropriately, 
!    * to specify WHICH eigenvalues you want to compute
!    * to supply a matrix-vector product
!      w <- A * v
!      in place of the call to AV( ) below.
!
!    Once usage of this routine is understood, you may wish to explore
!    the other available options to improve convergence, to solve generalized
!    problems, etc.  Look at the file ex-sym.doc in DOCUMENTS directory.
!
!
!    For this problem, we want to solve A*x = lambda*x in regular mode.
!
!    A is derived from the central difference discretization
!    of the 2-dimensional Laplacian on the unit square with
!    zero Dirichlet boundary condition.
!
!    OP = A  and  B = I.
!    Assume "call av (n,x,y)" computes y = A*x
!    Use mode 1 of SSAUPD.
!
!  Author:
!
!    Richard Lehoucq, Danny Sorensen, Chao Yang,
!    Department of Computational and Applied Mathematics,
!    Rice University,
!    Houston, Texas.
!
!  Storage:
! 
!    The maximum dimensions for all arrays are set here to accommodate 
!    a problem size of N <= MAXN
!
!    NEV is the number of eigenvalues requested.
!    See specifications for ARPACK usage below.
!
!    NCV is the largest number of basis vectors that will be used in 
!    the Implicitly Restarted Arnoldi Process.  Work per major iteration is
!    proportional to N*NCV*NCV.
!
!    You must set: 
! 
!    MAXN:   Maximum dimension of the A allowed. 
!    MAXNEV: Maximum NEV allowed. 
!    MAXNCV: Maximum NCV allowed. 
!
  implicit none

  integer, parameter :: maxn = 256
  integer, parameter :: maxnev = 10
  integer, parameter :: maxncv = 25

  integer, parameter :: ldv = maxn

  intrinsic abs
  real ax(maxn)
  character bmat  
  real d(maxncv,2)
  integer ido
  integer ierr
  integer info
  integer iparam(11)
  integer ipntr(11)
  integer ishfts
  integer j
  integer lworkl
  integer maxitr
  integer mode1
  integer n
  integer nconv
  integer ncv
  integer nev
  integer nx
  real resid(maxn)
  logical rvec
  external saxpy
  logical select(maxncv)
  real sigma
  real, external :: snrm2
  real tol
  real v(ldv,maxncv)
  character ( len = 2 ) which
  real workl(maxncv*(maxncv+8))
  real workd(3*maxn)
  real, parameter :: zero = 0.0E+00
!
!  The following include statement and assignments control trace output 
!  from the internal actions of ARPACK.  See debug.doc in the
!  DOCUMENTS directory for usage.  
!
!  Initially, the most useful information will be a breakdown of
!  time spent in the various stages of computation given by setting 
!  msaupd = 1.
!
  include 'debug.h'

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SSSIMP:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A simple ARPACK calling program.'
  write ( *, '(a)' ) '  This program defines an eigenproblem for a'
  write ( *, '(a)' ) '  symmetric matrix.'

  ndigit = -3
  logfil = 6
  msgets = 0
  msaitr = 0
  msapps = 0
  msaupd = 1
  msaup2 = 0
  mseigt = 0
  mseupd = 0
!
!  Set dimensions for this problem.
!
  nx = 10
  n = nx * nx
!
!  Specifications for ARPACK usage are set below:
!
!  1) NEV = 4 asks for 4 eigenvalues to be computed.                            !
!  2) NCV = 20 sets the length of the Arnoldi factorization.
!
!  3) This is a standard problem(indicated by bmat  = 'I')
!
!  4) Ask for the NEV eigenvalues of largest magnitude
!     (indicated by which = 'LM')
!
!  See documentation in SSAUPD for the other options SM, LA, SA, LI, SI.
!
!  NEV and NCV must satisfy the following conditions:
!
!    NEV <= MAXNEV
!    NEV + 1 <= NCV <= MAXNCV
!
  nev = 4
  ncv = 20
  bmat = 'I'
  which = 'LM'

  if ( maxn < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SSSIMP - Fatal error!'
    write ( *, '(a)' ) '  N is greater than MAXN '
    stop
  else if ( maxnev < nev ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SSSIMP - Fatal error!'
    write ( *, '(a)' ) '  NEV is greater than MAXNEV '
    stop
  else if ( maxncv < ncv ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SSSIMP - Fatal error!'
    write ( *, '(a)' ) '  NCV is greater than MAXNCV '
    stop
  end if
!
!  Specification of stopping rules and initial
!  conditions before calling SSAUPD
!
!  TOL determines the stopping criterion.  Expect
!    abs(lambdaC - lambdaT) < TOL*abs(lambdaC)
!  computed   true
!  If TOL <= 0, then TOL <- macheps (machine precision) is used.
!
!  IDO is the REVERSE COMMUNICATION parameter
!  used to specify actions to be taken on return
!  from SSAUPD. (See usage below.)
!  It MUST initially be set to 0 before the first
!  call to SSAUPD.
!
!  INFO on entry specifies starting vector information
!  and on return indicates error codes
!  Initially, setting INFO=0 indicates that a
!  random starting vector is requested to 
!  start the ARNOLDI iteration.  Setting INFO to
!  a nonzero value on the initial call is used 
!  if you want to specify your own starting 
!  vector. (This vector must be placed in RESID.)
!
!  The work array WORKL is used in SSAUPD as workspace.  Its dimension
!  LWORKL is set as illustrated below. 
!
  lworkl = ncv * ( ncv + 8 )
  tol = zero
  info = 0
  ido = 0
!
!  Specification of Algorithm Mode:
!
!  This program uses the exact shift strategy
!  (indicated by setting PARAM(1) = 1).
!
!  IPARAM(3) specifies the maximum number of Arnoldi iterations allowed.  
!
!  Mode 1 of SSAUPD is used (IPARAM(7) = 1). 
!
!  All these options can be changed by the user.  For details see the
!  documentation in SSAUPD.
!
  ishfts = 1
  maxitr = 300
  mode1 = 1

  iparam(1) = ishfts

  iparam(3) = maxitr

  iparam(7) = mode1
!
!  MAIN LOOP (Reverse communication loop)
!
!  Repeatedly call SSAUPD and take actions indicated by parameter 
!  IDO until convergence is indicated or MAXITR is exceeded.
!
  do

    call ssaupd ( ido, bmat, n, which, nev, tol, resid, &
      ncv, v, ldv, iparam, ipntr, workd, workl, &
      lworkl, info )

    if ( ido /= -1 .and. ido /= 1 ) then
      exit
    end if
!
!  Perform matrix vector multiplication
!
!    y <--- OP*x
!
!  The user supplies a matrix-vector multiplication routine that takes
!  workd(ipntr(1)) as the input, and return the result to workd(ipntr(2)).
!
    call av ( nx, workd(ipntr(1)), workd(ipntr(2)) )

   end do
!
!  Either we have convergence or there is an error.
!
  if ( info < 0 ) then
!
!  Error message. Check the documentation in SSAUPD.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SSSIMP - Fatal error!'
    write ( *, '(a,i6)' ) '  Error with SSAUPD, INFO = ', info
    write ( *, '(a)' ) '  Check documentation in SSAUPD.'

  else
!
!  No fatal errors occurred.
!  Post-Process using SSEUPD.
!
!  Computed eigenvalues may be extracted.
!
!  Eigenvectors may be also computed now if
!  desired.  (indicated by rvec = .true.)
!
!  The routine SSEUPD now called to do this
!  post processing (Other modes may require
!  more complicated post processing than mode1.)
!
    rvec = .true.

    call sseupd ( rvec, 'All', select, d, v, ldv, sigma, &
      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
      iparam, ipntr, workd, workl, lworkl, ierr )
!
!  Eigenvalues are returned in the first column of the two dimensional 
!  array D and the corresponding eigenvectors are returned in the first 
!  NCONV (=IPARAM(5)) columns of the two dimensional array V if requested.
!
!  Otherwise, an orthogonal basis for the invariant subspace corresponding 
!  to the eigenvalues in D is returned in V.
!
    if ( ierr /= 0 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SSSIMP - Fatal error!'
      write ( *, '(a,i6)' ) '  Error with SSEUPD, IERR = ', ierr
      write ( *, '(a)' ) '  Check the documentation of SSEUPD.'
!
!  Compute the residual norm
!
!    ||  A*x - lambda*x ||
! 
!  for the NCONV accurately computed eigenvalues and 
!  eigenvectors.  (iparam(5) indicates how many are 
!  accurate to the requested tolerance)
!
    else

      nconv =  iparam(5)

      do j = 1, nconv
        call av ( nx, v(1,j), ax )
        call saxpy ( n, -d(j,1), v(1,j), 1, ax, 1 )
        d(j,2) = snrm2 ( n, ax, 1)
        d(j,2) = d(j,2) / abs ( d(j,1) )
      end do
!
!  Display computed residuals.
!
      call smout ( 6, nconv, 2, d, maxncv, -6, &
        'Ritz values and relative residuals' )

    end if
!
!  Print additional convergence information.
!
    if ( info == 1) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Maximum number of iterations reached.'
    else if ( info == 3) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No shifts could be applied during implicit' &
        // ' Arnoldi update, try increasing NCV.'
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SSSIMP:'
    write ( *, '(a)' ) '====== '
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Size of the matrix is ', n
    write ( *, '(a,i6)' ) '  The number of Ritz values requested is ', nev
    write ( *, '(a,i6)' ) &
      '  The number of Arnoldi vectors generated (NCV) is ', ncv
    write ( *, '(a)' ) '  What portion of the spectrum: ' // which
    write ( *, '(a,i6)' ) &
      '  The number of converged Ritz values is ', nconv
    write ( *, '(a,i6)' ) &
      '  The number of Implicit Arnoldi update iterations taken is ', iparam(3)
    write ( *, '(a,i6)' ) '  The number of OP*x is ', iparam(9)
    write ( *, '(a,g14.6)' ) '  The convergence criterion is ', tol

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SSSIMP:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine av ( nx, v, w )

!*******************************************************************************
!
!! AV computes w <- A * V where A is a discretized Laplacian.
!
!  Discussion:
!
!    The matrix used is the 2 dimensional discrete Laplacian on unit
!    square with zero Dirichlet boundary condition.
!
!    Compute w <--- OP * v, where OP is the nx*nx by nx*nx block
!    tridiagonal matrix
!
!      | T -I           |
!      |-I  T -I        |
!      |   -I  T        |
!      |        ...  -I |
!      |           -I T |
!
!  Parameters:
!
!    Input, integer NX, the length of the vectors.
!
!    Input, real V(NX), the vector to be operated on by A.
!
!    Output, real W(NX), the result of A*V.
!
  implicit none

  integer nx

  real h2
  integer j
  integer lo
  integer n2
  real, parameter :: one = 1.0E+00
  real v(nx*nx)
  real w(nx*nx)

  call tv ( nx, v(1), w(1) )

  call saxpy ( nx, -one, v(nx+1), 1, w(1), 1 )

  do j = 2, nx-1
    lo = (j-1) * nx
    call tv ( nx, v(lo+1), w(lo+1) )
    call saxpy ( nx, -one, v(lo-nx+1), 1, w(lo+1), 1 )
    call saxpy ( nx, -one, v(lo+nx+1), 1, w(lo+1), 1 )
  end do

  lo = (nx-1) * nx
  call tv ( nx, v(lo+1), w(lo+1) )
  call saxpy ( nx, -one, v(lo-nx+1), 1, w(lo+1), 1 )
!
!  Scale the vector W by (1/H^2), where H is the mesh size.
!
  n2 = nx * nx
  h2 = one / real ( ( nx + 1 ) * ( nx + 1 ) )
  call sscal ( n2, one/h2, w, 1 )

  return
end
subroutine tv ( nx, x, y )

!*******************************************************************************
!
!! TV computes y <-- T*x.
!
!  Discussion:
!
!    The full matrix A is the 2 dimensional discrete Laplacian on unit
!    square with zero Dirichlet boundary condition.
!
!    A can be represented as the nx*nx by nx*nx block tridiagonal matrix
!
!      | T -I          |
!      |-I  T -I       |
!      |   -I  T       |
!      |        ...  -I|
!      |           -I T|
!
!    It is the job of this routine to compute products of the form T*x.
!
!  Parameters:
!
!    Input, integer NX, the length of the vectors.
!
!    Input, real X(NX), the vector to be multiplied by T.
!
!    Output, real Y(NX), the product T*X.
!
  implicit none

  integer nx

  integer j
  real dd
  real dl
  real du
  real x(nx)
  real y(nx)
!
!  Compute the matrix vector multiplication y<---T*x
!  where T is a nx by nx tridiagonal matrix with DD on the
!  diagonal, DL on the subdiagonal, and DU on the superdiagonal.
!
  dd  = 4.0E+00
  dl  = -1.0E+00
  du  = -1.0E+00

  y(1) =  dd * x(1) + du * x(2)

  do j = 2, nx-1
     y(j) = dl * x(j-1) + dd * x(j) + du * x(j+1)
  end do

  y(nx) = dl * x(nx-1) + dd * x(nx)

  return
end
